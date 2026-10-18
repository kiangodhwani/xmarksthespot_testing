"""Evaluation-only comparison and contextual MIND diagnostics; never selects matches."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from scipy.ndimage import gaussian_filter
from landmarks import mind
from landmarks.transforms import apply_points
from run_matching_baseline import error_mm


def scoring_patch(image,valid,point,radius,angle):
    lo,hi=np.percentile(image[valid],[1,99])
    normal=np.clip((image.astype(float)-lo)/max(hi-lo,1e-12),0,1)
    normal[~valid]=0
    smooth=gaussian_filter(normal,max(.7,radius/20*.7),truncate=4.,mode='constant')
    return mind.sample_patch(smooth,point,radius,angle)


def box(ax,point,radius,angle,label,color):
    theta=np.deg2rad(angle);c,s=np.cos(theta),np.sin(theta)
    rotation=np.array([[c,-s],[s,c]])
    corners=np.array([[-1,-1],[1,-1],[1,1],[-1,1]])*radius@rotation.T+point
    ax.add_patch(Polygon(corners,fill=False,edgecolor=color,lw=1.8))
    ax.plot(*point,'o',color=color,ms=3)
    ax.annotate(label,corners[0],xytext=(2,6),textcoords='offset points',color=color,
        fontsize=10,bbox=dict(facecolor='white',edgecolor=color,alpha=.9,pad=2))


def context(pair_dir,out):
    # Saved inference is read first. All truth access below is evaluation-only.
    with np.load(out/'predictions.npz') as d: pred={k:d[k] for k in d.files}
    with np.load(pair_dir/'inputs.npz') as d: images={k:d[k] for k in d.files}
    legal=np.flatnonzero((pred['similarity']>=0).any(axis=1))
    query=int(legal[len(legal)//2])  # fixed median usable index, never best correctness
    pa=pred['mri_candidates'][query];j=int(np.argmax(pred['similarity'][query]))
    pb=pred['histology_candidates'][j];chosen=pred['best_radius_angle_indices'][query,j]
    reference=json.loads((pair_dir/'truth/geometry.json').read_text())
    expected=apply_points(pa[None],np.asarray(reference['original_to_transformed']))[0]
    a,av=mind.describe(images['mri'],images['mri_valid'],pa[None])
    b,bv=mind.describe(images['histology'],images['histology_valid'],expected[None])
    reference_score,reference_choice=mind.similarity(a,av,b,bv,return_choices=True)
    options=[('Reference diagnostic',expected,reference_choice[0,0],float(reference_score[0,0]),'#0072b2','A1','B'),
             ('Algorithm top-scoring candidate',pb,chosen,float(pred['similarity'][query,j]),'#d55e00','A2','C')]
    fig,axs=plt.subplots(3,2,figsize=(14,17),gridspec_kw=dict(height_ratios=[1.25,1,1]))
    for ax,name,title in zip(axs[0],['mri','histology'],['MRI: query context','Transformed histology: reference B / prediction C']):
        image=images[name];valid=images[name+'_valid'];lo,hi=np.percentile(image[valid],[1,99])
        ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
        ax.set_facecolor('#ddd');ax.set_title(title);ax.set_xlabel('Column (working px)');ax.set_ylabel('Row (working px)')
    records=[]
    for row,(label,target,choice,score,color,left_label,right_label) in enumerate(options,1):
        if score<0:
            for ax in axs[row]:ax.text(.5,.5,'No fully supported reference patch',ha='center');ax.axis('off')
            continue
        ka,kb,angle_index=map(int,choice)
        ra,rb=mind.RADII[ka],mind.RADII[kb];angle=mind.ANGLES[angle_index]
        box(axs[0,0],pa,ra,0,left_label,color)
        box(axs[0,1],target,rb,angle,right_label,color)
        patches=[scoring_patch(images['mri'],images['mri_valid'],pa,ra,0),
                 scoring_patch(images['histology'],images['histology_valid'],target,rb,angle)]
        check=float(mind.patch_descriptor(patches[0])@mind.patch_descriptor(patches[1]))
        if not np.isclose(check,score,atol=1e-6):raise ValueError('Displayed patches do not reproduce score')
        for ax,patch,radius,rotation,tag,modality in zip(axs[row],patches,[ra,rb],[0,angle],[left_label,right_label],['MRI','histology']):
            ax.imshow(patch,cmap='gray',vmin=0,vmax=1,extent=[-20.5*radius/20,20.5*radius/20,20.5*radius/20,-20.5*radius/20])
            # Scoring-context centres, with one actual local SSD patch shown at centre.
            centres=np.array([-12,-6,0,6,12])*radius/20
            xx,yy=np.meshgrid(centres,centres);ax.plot(xx.ravel(),yy.ravel(),'.',color=color,ms=3)
            half=2*radius/20
            ax.add_patch(Polygon([[-half,-half],[half,-half],[half,half],[-half,half]],fill=False,edgecolor=color,lw=1.5))
            ax.set_title(f'{tag}: {modality}; radius {radius:g} px; sampled at {rotation}° CW\n41×41 samples; SSD 5×5 spans {4*radius/20:g} px between outer centres',fontsize=10)
            ax.set_xlabel('Canonical horizontal offset (source-image px)');ax.set_ylabel('Canonical vertical offset (source-image px)')
        error=float(error_mm(pa[None],target[None],reference)[0])
        axs[row,1].text(.5,-.21,f'{label}: cosine {score:.4f}; reference error {error:.2f} mm',transform=axs[row,1].transAxes,ha='center',color=color,fontsize=10)
        records.append(dict(label=label,source_radius_px=ra,target_radius_px=rb,target_sampling_clockwise_degrees=angle,
            score=score,reference_error_mm=error,source_centre=pa.tolist(),target_centre=target.tolist()))
    accepted=bool(np.any(np.all(pred['distinctive_pairs']==[query,j],axis=1)))
    retained=bool(np.any(np.all(pred['retained_pairs']==[query,j],axis=1)))
    fig.suptitle(f'{pair_dir.name}: fixed median usable MRI candidate #{query}; top histology candidate #{j}\nBlue = reference diagnostic only; orange = algorithm preference (filter accepted: {accepted}; retained: {retained})',fontsize=12)
    fig.subplots_adjust(left=.07,right=.98,top=.90,bottom=.14,hspace=.65,wspace=.24)
    fig.text(.5,.015,'Scores search 3×3 radii and 12 rotations, identically for reference and candidate. Cosine is not a probability.\nDots: 25 context centres; central box: one 5×5 SSD patch (also compared at eight offsets). Gaussian prefilter adds support.\nReference B is an approximate-alignment expected location, not verified anatomy or an inference input.\nBigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(out/'evaluation/descriptor_context.png',dpi=120);plt.close(fig)
    (out/'evaluation/descriptor_context.json').write_text(json.dumps(dict(query_index=query,target_index=j,examples=records),indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mind',type=Path,default=ROOT/'results/matching/mind_v1')
    parser.add_argument('--polar',type=Path,default=ROOT/'results/matching/polar_backend_check')
    args=parser.parse_args(); rows=[];reports={}
    for method,base in [('polar',args.polar),('mind',args.mind)]:
        reports[method]=json.loads((base/'summary.json').read_text())
        for report in reports[method]:
            row=dict(method=method,pair=report['pair_id'],inference_seconds=report['inference_seconds'])
            for stage,data in report['stages'].items():
                for key in ['count','count_within_2mm','count_within_5mm','median_error_mm']:
                    row[f'{stage}_{key}']=data[key]
            rows.append(row)
    for a,b in zip(reports['polar'],reports['mind']):
        assert a['pair_id']==b['pair_id']
        for key in ['mri_candidates','histology_candidates']:
            with np.load(args.polar/a['pair_id']/'predictions.npz') as x,np.load(args.mind/b['pair_id']/'predictions.npz') as y:
                np.testing.assert_array_equal(x[key],y[key])
    with (args.mind/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    controls={k:json.loads((base/'histology_control.json').read_text()) for k,base in [('polar',args.polar),('mind',args.mind)]}
    (args.mind/'comparison.json').write_text(json.dumps(dict(rows=rows,controls=controls,candidates_identical=True,
        thresholds_unchanged=True,note='One slice, six augmentations; runtimes descriptive, runs overlapped; no held-out validation'),indent=2))
    fig,axs=plt.subplots(1,3,figsize=(15,5))
    for ax,stage,title in zip(axs,['mutual','distinctive','retained'],['Mutual proposals','After distinctiveness','After RANSAC']):
        x=np.arange(6)
        for method,offset,color in [('polar',-.18,'#0072b2'),('mind',.18,'#d55e00')]:
            values=[r['stages'][stage]['count_within_2mm'] for r in reports[method]]
            bars=ax.bar(x+offset,values,.36,label=method,color=color);ax.bar_label(bars,padding=2)
        ax.set_xticks(x,['Fixed','001','002','003','004','005']);ax.set_ylabel('Pairs within 2 mm of reference');ax.set_title(title);ax.set_ylim(0,14);ax.legend()
    fig.suptitle('Fixed comparison: correct-pair counts at each stage (all pairs included)')
    fig.subplots_adjust(left=.05,right=.99,top=.83,bottom=.26,wspace=.28)
    fig.text(.5,.02,'Histology control: polar 117/117 retained within 2 mm; MIND-style 65/65.\nApproximate anatomical reference; six transformations of ONE slice. Same candidates and thresholds.\nBigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(args.mind/'comparison.png',dpi=140);plt.close(fig)
    for report in reports['mind']:
        name=report['pair_id'];pair_dir=ROOT/'results/transforms'/name if name=='coronal_y8_demo' else ROOT/'results/random_pairs/y8_seed42'/name
        context(pair_dir,args.mind/name)
    print('Comparison and six context figures saved to',args.mind)

if __name__=='__main__':main()
