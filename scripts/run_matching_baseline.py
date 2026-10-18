"""One end-to-end baseline: image-only inference, then separate reference evaluation."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import argparse
import csv
import hashlib
import json
import os
import tempfile
import time
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np
from scipy.spatial import ConvexHull,QhullError
from landmarks.masking import otsu_mask,histology_hysteresis
from landmarks.candidates import detect_candidates
from landmarks.descriptors import describe,similarity
from landmarks.matching import propose,geometric_filter
from landmarks.planes import pixel_to_world
from landmarks.transforms import apply_points


def infer(images,descriptor="polar"):
    """Images/support only. No truth directory, seed, affine or transformation."""
    from landmarks import mind
    describe_fn = mind.describe if descriptor == "mind" else describe
    masks={};points={};banks={};usable={};detector={}
    masks['mri']=otsu_mask(images['mri'],images['mri_valid'],True)[0]['cleaned']
    masks['histology']=histology_hysteresis(images['histology'],images['histology_valid'],weak_fraction=.01)[0]['mask']
    for name in ['mri','histology']:
        points[name],detector[name]=detect_candidates(images[name],images[name+'_valid'],masks[name])
        banks[name],usable[name]=describe_fn(images[name],images[name+'_valid'],points[name][:,:2])
    choices=None
    if descriptor == 'mind':
        scores,choices=mind.similarity(banks['mri'],usable['mri'],banks['histology'],usable['histology'],return_choices=True)
    else:
        scores=similarity(banks['mri'],usable['mri'],banks['histology'],usable['histology'])
    mutual,accepted,ratios=propose(scores,max_ratio=.9,min_similarity=.75)
    pa=points['mri'][:,:2];pb=points['histology'][:,:2]
    model,keep,geometry=geometric_filter(pa[accepted[:,0]],pb[accepted[:,1]])
    return dict(points=points,masks=masks,scores=scores,mutual=mutual,accepted=accepted,ratios=ratios,
        model=model,keep=keep,geometry=geometry,detector=detector,choices=choices,
        usable_counts={k:int(v.any(axis=tuple(range(1,v.ndim))).sum()) for k,v in usable.items()})


def error_mm(pa,pb,reference):
    a=pixel_to_world(pa,np.asarray(reference['original_pixel_to_world_ras_mm']))
    b=pixel_to_world(pb,np.asarray(reference['transformed_pixel_to_world_ras_mm']))
    return np.linalg.norm(a-b,axis=1)


def spread(points):
    if len(points)<3: return dict(hull_area_mm2=0.,minor_major_spread_ratio=0.)
    singular=np.linalg.svd(points-points.mean(axis=0),compute_uv=False)
    ratio=float(singular[1]/singular[0]) if singular[0]>0 else 0.
    try: area=float(ConvexHull(points).volume)
    except QhullError: area=0.
    return dict(hull_area_mm2=area,minor_major_spread_ratio=ratio)


def evaluate_pairs(pairs,pa,pb,reference,mask_area,rng):
    errors=error_mm(pa[pairs[:,0]],pb[pairs[:,1]],reference)
    result=dict(count=len(pairs),count_within_2mm=int((errors<=2).sum()),count_within_5mm=int((errors<=5).sum()),
        fraction_within_2mm=float((errors<=2).mean()) if len(errors) else None,
        fraction_within_5mm=float((errors<=5).mean()) if len(errors) else None,
        median_error_mm=float(np.median(errors)) if len(errors) else None)
    correct=pa[pairs[errors<=2,0]]*.2
    result['correct_within_2mm_spread']=spread(correct)
    result['correct_hull_fraction_of_mri_mask_area']=result['correct_within_2mm_spread']['hull_area_mm2']/mask_area
    if len(pairs):
        chance=[]
        for _ in range(200):
            a=rng.choice(len(pa),size=len(pairs),replace=False);b=rng.choice(len(pb),size=len(pairs),replace=False)
            e=error_mm(pa[a],pb[b],reference);chance.append([(e<=2).mean(),(e<=5).mean()])
        chance=np.asarray(chance)
        result['random_candidate_pairing']=dict(trials=200,mean_within_2mm=float(chance[:,0].mean()),mean_within_5mm=float(chance[:,1].mean()),
            note='Same count, independent unique candidate pairing; not a null distribution conditioned on RANSAC selection')
    else: result['random_candidate_pairing']=None
    return result,errors


def show_pair(images,masks,pa,pb,stages,stage_errors,report,path):
    fig=plt.figure(figsize=(14,15),dpi=140)
    positions=[(.04,.67,.43,.24),(.53,.67,.43,.24),(.04,.37,.43,.24),(.53,.37,.43,.24),(.04,.07,.43,.24),(.53,.07,.43,.24)]
    names=['Mutual nearest-neighbour proposals','After distinctiveness filter','After RANSAC (final retained)']
    for row,((stage,pairs),errors,title) in enumerate(zip(stages.items(),stage_errors.values(),names)):
        axes=[fig.add_axes(positions[2*row]),fig.add_axes(positions[2*row+1])]
        for ax,name in zip(axes,['mri','histology']):
            im=images[name];lo,hi=np.percentile(im[images[name+'_valid']],[1,99])
            ax.imshow(np.ma.array(im,mask=~images[name+'_valid']),cmap='gray',vmin=lo,vmax=hi)
            rr,cc=np.where(masks[name]);ax.set_xlim(cc.min()-20,cc.max()+20);ax.set_ylim(rr.max()+20,rr.min()-20)
            ax.set_xticks([]);ax.set_yticks([]);ax.set_facecolor('#ddd')
            ax.set_title(('MRI' if name=='mri' else 'Transformed histology')+f' — {title}',fontsize=10)
        # Fixed display subset: highest descriptor scores, not selected by correctness.
        display=np.arange(min(12,len(pairs)))
        if len(pairs)>12:
            scores=[report['scores_for_plot'][stage][k] for k in range(len(pairs))]
            display=np.argsort(-np.asarray(scores))[:12]
        for k in display:
            i,j=pairs[k];color='#168a35' if errors[k]<=2 else '#d28b00' if errors[k]<=5 else '#d12c3f'
            for ax,point in zip(axes,[pa[i],pb[j]]):
                ax.plot(*point,'o',ms=4,mfc=color,mec='white',mew=.5)
                ax.annotate(str(k+1),point,xytext=(4,4),textcoords='offset points',fontsize=7,color=color,
                    bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=.4))
            fig.add_artist(ConnectionPatch(xyA=pa[i],coordsA=axes[0].transData,xyB=pb[j],coordsB=axes[1].transData,
                color=color,lw=.8,alpha=.6,clip_on=False))
        summary=report['stages'][stage]
        fig.text(.5,positions[2*row][1]-.025,f'{len(pairs)} total pairs; {summary["count_within_2mm"]} within 2 mm; {summary["count_within_5mm"]} within 5 mm. Showing {len(display)} highest-score proposals.',ha='center',fontsize=10)
        if not len(pairs): axes[1].text(.5,.5,'No pairs retained\n'+report['geometric_filter']['status'].replace('_',' '),transform=axes[1].transAxes,ha='center',va='center',bbox=dict(facecolor='white',alpha=.95),fontsize=12)
    fig.suptitle(report.get('descriptor','polar')+' descriptor: proposed pairs → distinctiveness → geometric filtering\nLine colours are EVALUATION ONLY: green ≤2 mm; amber 2–5 mm; red >5 mm.',y=.98,fontsize=13)
    fig.text(.5,.015,'Reference alignment is approximate. A geometric inlier is not automatically anatomically correct.\nBigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(path);plt.close(fig)


def run(pair_dir,out,descriptor="polar"):
    from landmarks import mind
    path=pair_dir/'inputs.npz';digest=hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path) as loaded: images={k:loaded[k] for k in loaded.files}
    started=time.perf_counter()
    predicted=infer(images,descriptor)
    elapsed=time.perf_counter()-started
    out.mkdir(parents=True)
    pa=predicted['points']['mri'][:,:2];pb=predicted['points']['histology'][:,:2]
    stages={'mutual':predicted['mutual'],'distinctive':predicted['accepted'],'retained':predicted['accepted'][predicted['keep']]}
    # Freeze outputs BEFORE opening reference metadata. No oracle-ranked patches.
    np.savez_compressed(out/'predictions.npz',mri_candidates=pa,histology_candidates=pb,
        mutual_pairs=stages['mutual'],distinctive_pairs=stages['distinctive'],retained_pairs=stages['retained'],
        similarity=predicted['scores'],
        **({'best_radius_angle_indices':predicted['choices']} if predicted['choices'] is not None else {}))
    prediction=dict(descriptor=descriptor,descriptor_config=mind.CONFIG if descriptor=='mind' else None,
        inference_seconds=elapsed,usable_candidate_counts=predicted['usable_counts'],detector=predicted['detector'],geometric_filter=predicted['geometry'],
        estimated_mri_to_histology=predicted['model'].tolist() if predicted['model'] is not None else None,
        parameters=dict(ratio_max=.9,min_similarity=.75,ransac_residual_histology_px=6.,minimum_consensus=4),
        input_sha256=digest,truth_access='None during inference')
    (out/'prediction.json').write_text(json.dumps(prediction,indent=2))
    reference=json.loads((pair_dir/'truth/geometry.json').read_text())
    evaluation_dir=out/'evaluation';evaluation_dir.mkdir()
    report=dict(descriptor=descriptor,inference_seconds=elapsed,usable_candidate_counts=predicted['usable_counts'],pair_id=pair_dir.name,candidate_counts=dict(mri=len(pa),histology=len(pb)),
        geometric_filter=predicted['geometry'],stages={},scores_for_plot={},
        scope='Actual independent candidate matches; reference used only after predictions saved',
        caveats=['One specimen and one plane with multiple augmentations, not independent held-out examples',
            'Reference alignment approximate; errors are reference-coordinate errors, not verified anatomical truth',
            'Fixed thresholds not tuned on evaluation; baseline may produce no final pairs'])
    errors={};rng=np.random.default_rng(7)
    area=float(predicted['masks']['mri'].sum()*.2**2)
    for name,pairs in stages.items():
        report['stages'][name],errors[name]=evaluate_pairs(pairs,pa,pb,reference,area,rng)
        report['scores_for_plot'][name]=predicted['scores'][pairs[:,0],pairs[:,1]].tolist()
    recovery=None
    if predicted['model'] is not None:
        estimated=predicted['model']; planted=np.asarray(reference['original_to_transformed'])
        theta=lambda m:float(np.rad2deg(np.arctan2(m[1,0],m[0,0])))
        angle_error=abs((theta(estimated)-theta(planted)+180)%360-180)
        scale=lambda m:float(np.sqrt(np.linalg.det(m[:2,:2])))
        eworld=pixel_to_world(apply_points(pa,estimated),np.asarray(reference['transformed_pixel_to_world_ras_mm']))
        aworld=pixel_to_world(pa,np.asarray(reference['original_pixel_to_world_ras_mm']))
        recovery=dict(rotation_absolute_error_degrees=angle_error,scale_absolute_error=abs(scale(estimated)-scale(planted)),
            effective_matrix_translation_error_px=float(np.linalg.norm(estimated[:2,2]-planted[:2,2])),
            median_whole_candidate_transform_error_mm=float(np.median(np.linalg.norm(eworld-aworld,axis=1))),
            note='Canvas-aware planted matrix; not raw shift alone')
    report['transform_recovery']=recovery
    # Final CSV implements the README's proposed external coordinate format.
    with (out/'landmarks.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['id','name','mri_x_mm','mri_y_mm','mri_z_mm','hist_col_px','hist_row_px','score'])
        for k,(i,j) in enumerate(stages['retained'],1):
            world=pixel_to_world(pa[i:i+1],np.asarray(reference['original_pixel_to_world_ras_mm']))[0]
            writer.writerow([k,'candidate_'+str(i),*world,*pb[j],float(predicted['scores'][i,j])])
    with (evaluation_dir/'all_proposals.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['stage','id','mri_col','mri_row','hist_col','hist_row','descriptor_similarity','reference_error_mm'])
        for name,pairs in stages.items():
            for k,((i,j),error) in enumerate(zip(pairs,errors[name]),1):writer.writerow([name,k,*pa[i],*pb[j],predicted['scores'][i,j],error])
    show_pair(images,predicted['masks'],pa,pb,stages,errors,report,out/'matches.png')
    report.pop('scores_for_plot')
    (evaluation_dir/'report.json').write_text(json.dumps(report,indent=2))
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    print(pair_dir.name, 'mutual',len(stages['mutual']),'distinctive',len(stages['distinctive']),
          'retained',len(stages['retained']),'status',predicted['geometry']['status'],flush=True)
    return report


def histology_control(output,descriptor="polar"):
    """Positive control: same descriptor/matcher/RANSAC on two histology images."""
    from landmarks import mind
    describe_fn=mind.describe if descriptor=="mind" else describe
    similarity_fn=mind.similarity if descriptor=="mind" else similarity
    fixed=ROOT/'results/transforms/coronal_y8_demo'
    with np.load(fixed/'truth/original_histology.npz') as d: original=d['image'];ov=d['valid']
    with np.load(fixed/'inputs.npz') as d: transformed=d['histology'];tv=d['histology_valid']
    points=[];banks=[];flags=[]
    for image,valid in [(original,ov),(transformed,tv)]:
        mask=histology_hysteresis(image,valid,weak_fraction=.01)[0]['mask']
        q,_=detect_candidates(image,valid,mask);bank,usable=describe_fn(image,valid,q[:,:2])
        points.append(q[:,:2]);banks.append(bank);flags.append(usable)
    scores=similarity_fn(banks[0],flags[0],banks[1],flags[1])
    mutual,pairs,_=propose(scores);model,keep,status=geometric_filter(points[0][pairs[:,0]],points[1][pairs[:,1]])
    np.savez_compressed(output/'histology_control_predictions.npz',source_candidates=points[0],target_candidates=points[1],
        distinctive_pairs=pairs,retained_pairs=pairs[keep],similarity=scores,
        estimated_model=model if model is not None else np.empty((0,0)))
    reference=json.loads((fixed/'truth/geometry.json').read_text())
    errors=error_mm(points[0][pairs[:,0]],points[1][pairs[:,1]],reference)
    result=dict(scope='Same-modality histology control; independent detections, not oracle locations',
        mutual=len(mutual),distinctive=len(pairs),retained=int(keep.sum()),
        retained_within_2mm=int((errors[keep]<=2).sum()),geometric_filter=status)
    (output/'histology_control.json').write_text(json.dumps(result,indent=2))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--descriptor',choices=['polar','mind'],default='polar')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output is None: args.output=ROOT/'results/matching'/('mind_v1' if args.descriptor=='mind' else 'polar_baseline')
    if args.output.exists():raise FileExistsError('Choose a new output directory')
    pairs=[ROOT/'results/transforms/coronal_y8_demo']+sorted((ROOT/'results/random_pairs/y8_seed42').glob('pair_*'))
    if len(pairs)!=6 or any(not (p/'inputs.npz').is_file() for p in pairs):
        raise FileNotFoundError('Expected fixed example and exactly five random pairs')
    args.output.mkdir(parents=True)
    if args.descriptor=='mind':
        from landmarks.mind import CONFIG
        (args.output/'experiment_config.json').write_text(json.dumps(CONFIG,indent=2))
    reports=[]
    for pair in pairs: reports.append(run(pair,args.output/pair.name,args.descriptor))
    with (args.output/'summary.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['pair','mutual','distinctive','retained','retained_within_2mm','retained_within_5mm','status'])
        for r in reports:
            stage=r['stages'];writer.writerow([r['pair_id'],stage['mutual']['count'],stage['distinctive']['count'],stage['retained']['count'],
                stage['retained']['count_within_2mm'],stage['retained']['count_within_5mm'],r['geometric_filter']['status']])
    (args.output/'summary.json').write_text(json.dumps(reports,indent=2))
    histology_control(args.output,args.descriptor)
    print('Saved:',args.output)

if __name__=='__main__':main()
