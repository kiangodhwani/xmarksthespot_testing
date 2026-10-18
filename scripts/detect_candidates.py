"""Detect independently; then evaluate spatial repeatability using held-out truth."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import json
import hashlib
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt
from landmarks.candidates import detect_candidates
from landmarks.masking import histology_hysteresis
from landmarks.transforms import apply_points
from landmarks.repeatability import repeatability


def inside(points,mask):
    p=np.rint(points).astype(int)
    valid=(p[:,0]>=0)&(p[:,0]<mask.shape[1])&(p[:,1]>=0)&(p[:,1]<mask.shape[0])
    result=np.zeros(len(p),bool); q=p[valid]; result[valid]=mask[q[:,1],q[:,0]]
    return result


def random_points(mask,count,rng,min_spacing=8):
    rows,cols=np.where(mask); selected=[]
    for _ in range(100000):
        if len(selected)==count: break
        i=rng.integers(len(rows)); p=np.array([cols[i],rows[i]],float)
        if not selected or np.min(np.sum((np.asarray(selected)-p)**2,axis=1))>=min_spacing**2: selected.append(p)
    if len(selected)!=count: raise ValueError('Could not sample matched-density random baseline')
    return np.asarray(selected).reshape(-1,2)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/candidates/harris_y8')
    args=p.parse_args()
    if args.output.exists(): raise FileExistsError('Choose a new --output directory')
    paths=[ROOT/'results/transforms/coronal_y8_demo/inputs.npz',ROOT/'results/transforms/coronal_y8_demo/truth/original_histology.npz',
           ROOT/'results/masks/coronal_y8_otsu/masks.npz',ROOT/'results/masks/coronal_y8_hysteresis/masks.npz']
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    with np.load(paths[0]) as d:
        images={'mri':d['mri'],'histology_transformed':d['histology']}
        valid={'mri':d['mri_valid'],'histology_transformed':d['histology_valid']}
    # Original histology is a separate diagnostic image, never supplied to the
    # transformed-image detector; no affine or planted transform enters detection.
    with np.load(paths[1]) as d: images['histology_original']=d['image']; valid['histology_original']=d['valid']
    with np.load(paths[2]) as d: tissue={'mri':d['mri_cleaned']}
    with np.load(paths[3]) as d: tissue['histology_transformed']=d['weak_0.01']
    m,_=histology_hysteresis(images['histology_original'],valid['histology_original'],weak_fraction=.01)
    tissue['histology_original']=m['mask']
    points={}; stats={}
    args.output.mkdir(parents=True)
    for name,image in images.items():
        points[name],stats[name]=detect_candidates(image,valid[name],tissue[name])
        with (args.output/(name+'.csv')).open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['col_px','row_px','relative_response','sigma_px','tile_id']);writer.writerows(points[name])
        print(name,stats[name]['count'],'candidates',flush=True)
    sensitivity={}
    for name,image in images.items():
        variants=[]
        for label,config in [('lower_threshold',{'threshold_rel':.005}),('higher_threshold',{'threshold_rel':.02}),('fine_scales',{'scales':(2.,4.)}),('coarse_scales',{'scales':(4.,8.)})]:
            q,detail=detect_candidates(image,valid[name],tissue[name],**config)
            overlap=repeatability(points[name][:,:2],q[:,:2],(2.,))[0]
            variants.append(dict(variant=label,count=len(q),fraction_near_default=overlap['fraction_of_smaller_set'],tolerance_working_pixels=2.))
        sensitivity[name]=variants
    (args.output/'sensitivity.json').write_text(json.dumps(sensitivity,indent=2))
    # Freeze detections before opening planted transform metadata.
    (args.output/'detector_report.json').write_text(json.dumps(dict(parameters=stats,input_sha256=hashes,
        note='No supplied labels, affine or planted transform used for candidate generation'),indent=2))
    geometry=json.loads((ROOT/'results/transforms/coronal_y8_demo/truth/geometry.json').read_text())
    t=np.asarray(geometry['original_to_transformed']); inverse=np.linalg.inv(t)
    canonical={name:(inverse if name=='histology_transformed' else np.eye(3)) for name in images}
    fig,axs=plt.subplots(1,3,figsize=(17,7),layout='constrained')
    for ax,name in zip(axs,['mri','histology_original','histology_transformed']):
        image=images[name]; lo,hi=np.percentile(image[valid[name]],[1,99])
        ax.imshow(np.ma.array(image,mask=~valid[name]),cmap='gray',vmin=lo,vmax=hi)
        q=points[name]
        ax.scatter(q[:,0],q[:,1],s=12,facecolors='none',edgecolors='#00d5ff',linewidths=.6)
        ax.scatter(q[:,0],q[:,1],s=2,color='#ff0055')
        # Display only crop from independent tissue extent; coordinates unchanged.
        rr,cc=np.where(tissue[name]); ax.set_xlim(cc.min()-25,cc.max()+25);ax.set_ylim(rr.max()+25,rr.min()-25)
        ax.set_title(name.replace('_',' ').title()+f'\n{len(q)} independent candidates')
        ax.set_xlabel('Column (px)'); ax.set_ylabel('Row (px)');ax.set_facecolor('#ddd')
    fig.suptitle('Multiscale Harris baseline — pink dots mark exact centres; cyan rings are display markers, not tissue regions')
    fig.get_layout_engine().set(rect=(0,.05,1,.95))
    fig.text(.5,.02,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(args.output/'candidates.png',dpi=160);plt.close(fig)
    rng=np.random.default_rng(20261018); evaluation={}
    for a,b in [('histology_original','histology_transformed'),('mri','histology_original'),('mri','histology_transformed')]:
        ca,cb=canonical[a],canonical[b]
        pa=points[a][:,:2];pb=points[b][:,:2]
        eligible_a=inside(apply_points(pa,np.linalg.inv(cb)@ca),valid[b])
        eligible_b=inside(apply_points(pb,np.linalg.inv(ca)@cb),valid[a])
        actual=repeatability(apply_points(pa[eligible_a],ca)*.2,apply_points(pb[eligible_b],cb)*.2)
        baseline=[]
        for _ in range(50):
            random={}
            for name,count in [(a,len(pa)),(b,len(pb))]:
                safe=distance_transform_edt(np.pad(valid[name],1))[1:-1,1:-1]>15
                random[name]=random_points(tissue[name]&safe,count,rng)
            ra,rb=random[a],random[b]
            ea=inside(apply_points(ra,np.linalg.inv(cb)@ca),valid[b])
            eb=inside(apply_points(rb,np.linalg.inv(ca)@cb),valid[a])
            baseline.append(repeatability(apply_points(ra[ea],ca)*.2,apply_points(rb[eb],cb)*.2))
        for i,item in enumerate(actual):
            chance=np.array([trial[i]['fraction_of_smaller_set'] for trial in baseline])
            item['random_mean_fraction']=float(chance.mean());item['random_5_95_percentile']=np.percentile(chance,[5,95]).tolist()
        evaluation[a+'__'+b]=dict(source_count=int(eligible_a.sum()),target_count=int(eligible_b.sum()),
            results=actual,meaning='One-to-one spatial proximity under reference geometry, not descriptor matches or anatomical precision')
    report=dict(evaluation=evaluation,random_trials=50,random_seed=20261018,
        baseline='Independent uniform tissue sampling with same per-image point counts and 8px minimum spacing; >15px from invalid support. Does not reproduce detector spatial tile occupancy.',
        tolerances_mm=[1,2,5],truth_usage='Evaluation only; transformed detections mapped into original 0.2mm grid',
        caveats=['Single slice and single transformation; no held-out validation',
            'Reference MRI–histology alignment is approximate',
            '5mm proximity can occur often by chance; does not establish correct correspondence',
            'Multiscale Harris not intrinsically scale invariant; candidate score is not confidence'])
    (args.output/'repeatability.json').write_text(json.dumps(report,indent=2))
    for path in paths:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==hashes[str(path.relative_to(ROOT))]
    print(json.dumps(report,indent=2)); print('Preview:',args.output/'candidates.png')

if __name__=='__main__': main()
