"""Freeze three proposal-filter/RANSAC predictions, then evaluate against reference."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
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
from landmarks.matching import geometric_filter
from landmarks.masking import otsu_mask,histology_hysteresis
from landmarks.proposal_ablation import VARIANTS,CONFIG,select_variants
from landmarks.planes import pixel_to_world
from landmarks.transforms import apply_points
from run_matching_baseline import evaluate_pairs,error_mm

LABELS=dict(current='Current: score + ratio',similarity_only='Score only; no ratio',all_mutual='All mutual; no score or ratio')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze(source,out,control=False):
    """Reads saved image-derived candidates/scores only. No pair reference access."""
    with np.load(source) as d: saved={k:d[k] for k in d.files}
    pa=saved['source_candidates' if control else 'mri_candidates']
    pb=saved['target_candidates' if control else 'histology_candidates']
    scores=saved['similarity']
    mutual,variants,audit=select_variants(scores)
    np.testing.assert_array_equal(variants['current'],saved['distinctive_pairs'])
    if not control:np.testing.assert_array_equal(mutual,saved['mutual_pairs'])
    out.mkdir(parents=True)
    results={}
    for variant,pairs in variants.items():
        started=time.perf_counter()
        model,keep,status=geometric_filter(pa[pairs[:,0]],pb[pairs[:,1]])
        elapsed=time.perf_counter()-started
        if variant=='current':np.testing.assert_array_equal(pairs[keep],saved['retained_pairs'])
        path=out/variant;path.mkdir()
        np.savez_compressed(path/'predictions.npz',source_candidates=pa,target_candidates=pb,
            input_pairs=pairs,retained_pairs=pairs[keep],similarity=scores,
            estimated_model=model if model is not None else np.empty((0,0)))
        result=dict(variant=variant,description=CONFIG['variants'][variant],geometric_filter=status,
            ransac_seconds=elapsed,source_prediction_sha256=sha(source),truth_access='None before predictions saved')
        (path/'prediction.json').write_text(json.dumps(result,indent=2))
        results[variant]=(pairs,pairs[keep],model,status,elapsed)
    (out/'proposal_audit.json').write_text(json.dumps(audit,indent=2))
    return dict(pa=pa,pb=pb,scores=scores,mutual=mutual,audit=audit,variants=results)


def recovery(model,pa,reference):
    if model is None:return None
    planted=np.asarray(reference['original_to_transformed'])
    theta=lambda m:float(np.rad2deg(np.arctan2(m[1,0],m[0,0])))
    scale=lambda m:float(np.sqrt(np.linalg.det(m[:2,:2])))
    eworld=pixel_to_world(apply_points(pa,model),np.asarray(reference['transformed_pixel_to_world_ras_mm']))
    aworld=pixel_to_world(pa,np.asarray(reference['original_pixel_to_world_ras_mm']))
    return dict(rotation_absolute_error_degrees=abs((theta(model)-theta(planted)+180)%360-180),
        scale_absolute_error=abs(scale(model)-scale(planted)),
        effective_matrix_translation_error_px=float(np.linalg.norm(model[:2,2]-planted[:2,2])),
        median_whole_candidate_transform_error_mm=float(np.median(np.linalg.norm(eworld-aworld,axis=1))))


def draw(images,data,reports,path,method,case):
    fig,axs=plt.subplots(3,2,figsize=(14,15))
    fig.subplots_adjust(left=.035,right=.98,top=.91,bottom=.09,hspace=.32,wspace=.22)
    ids={tuple(pair):i+1 for i,pair in enumerate(data['mutual'])}
    for row,variant in enumerate(VARIANTS):
        pairs=data['variants'][variant][1];report=reports[variant]
        errors=np.asarray(report['retained_reference_errors_mm'])
        order=np.argsort(-data['scores'][pairs[:,0],pairs[:,1]],kind='stable')[:12]
        for ax,name in zip(axs[row],['mri','histology']):
            im=images[name];valid=images[name+'_valid'];lo,hi=np.percentile(im[valid],[1,99])
            ax.imshow(np.ma.array(im,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
            ax.set_facecolor('#ddd');ax.set_xticks([]);ax.set_yticks([])
            ax.set_title(f'{"MRI" if name=="mri" else "Transformed histology"} — {LABELS[variant]}',fontsize=10)
        for k in order:
            i,j=pairs[k];color='#168a35' if errors[k]<=2 else '#d28b00' if errors[k]<=5 else '#d12c3f'
            for ax,point in zip(axs[row],[data['pa'][i],data['pb'][j]]):
                ax.plot(*point,'o',ms=4,color=color)
                ax.annotate(str(ids[(i,j)]),point,xytext=(3,4),textcoords='offset points',fontsize=7,color=color,
                    bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=.3))
            fig.add_artist(ConnectionPatch(xyA=data['pa'][i],coordsA=axs[row,0].transData,
                xyB=data['pb'][j],coordsB=axs[row,1].transData,color=color,lw=.8,alpha=.65,clip_on=False))
        final=report['stages']['retained'];inp=report['stages']['input']
        axs[row,0].text(.5,-.09,f'RANSAC input: {inp["count"]} pairs; {inp["count_within_2mm"]} within 2 mm',
            transform=axs[row,0].transAxes,ha='center',fontsize=9)
        axs[row,1].text(.5,-.09,f'Final: {final["count"]} pairs; {final["count_within_2mm"]} ≤2 mm; {final["count_within_5mm"]} ≤5 mm. Showing {len(order)}.',
            transform=axs[row,1].transAxes,ha='center',fontsize=9)
        if not len(pairs):axs[row,1].text(.5,.5,'No model: insufficient consensus',transform=axs[row,1].transAxes,
            ha='center',bbox=dict(facecolor='white',alpha=.9))
    fig.suptitle(f'{method}: {case} — final RANSAC pairs under three fixed proposal filters\nEvaluation-only colours: green ≤2 mm; amber 2–5 mm; red >5 mm.',fontsize=12)
    fig.text(.5,.02,'Stable labels identify original mutual proposals across rows. If >12 retained, display highest cosine scores; metrics include ALL pairs.\nSame descriptor/candidates and RANSAC settings; reference geometry only used after prediction. Cosine is not confidence.\nApproximate reference alignment; six cases share ONE slice. BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(path,dpi=120);plt.close(fig)


def evaluate(data,reference,mask_area,out,method,case):
    pa,pb=data['pa'],data['pb'];reports={}
    mutual_errors=error_mm(pa[data['mutual'][:,0]],pb[data['mutual'][:,1]],reference)
    audit=[dict(item,reference_error_mm=float(error)) for item,error in zip(data['audit'],mutual_errors)]
    with (out/'proposal_audit_evaluation.csv').open('w',newline='') as f:
        fields=['source_index','target_index','similarity','distance_ratio','ratio_defined','passes_ratio','passes_similarity','reference_error_mm']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(audit)
    for variant,(pairs,retained,model,status,elapsed) in data['variants'].items():
        report=dict(method=method,case=case,variant=variant,geometric_filter=status,ransac_seconds=elapsed,stages={})
        rng=np.random.default_rng(7)
        for label,p in [('input',pairs),('retained',retained)]:
            report['stages'][label],errors=evaluate_pairs(p,pa,pb,reference,mask_area,rng)
            if label=='retained':report['retained_reference_errors_mm']=errors.tolist()
        report['transform_recovery']=recovery(model,pa,reference)
        # Explicit diagnostic flag, not a new model-acceptance rule.
        final=report['stages']['retained']
        report['model_with_majority_pairs_beyond_5mm']=bool(model is not None and final['fraction_within_5mm']<.5)
        report['rejected_mutual_within_2mm']=dict(
            ratio_only=sum(r['reference_error_mm']<=2 and not r['passes_ratio'] and r['passes_similarity'] for r in audit),
            score_only=sum(r['reference_error_mm']<=2 and r['passes_ratio'] and not r['passes_similarity'] for r in audit),
            both=sum(r['reference_error_mm']<=2 and not r['passes_ratio'] and not r['passes_similarity'] for r in audit))
        (out/variant/'evaluation.json').write_text(json.dumps(report,indent=2));reports[variant]=report
    return reports


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mind',type=Path,default=ROOT/'results/matching/mind_v1')
    p.add_argument('--polar',type=Path,default=ROOT/'results/matching/polar_backend_check')
    p.add_argument('--output',type=Path,default=ROOT/'results/matching/filter_ablation_v1')
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output directory')
    sources={'mind':args.mind,'polar':args.polar}
    cases=['coronal_y8_demo']+[f'pair_{i:03d}' for i in range(1,6)]
    required=[base/case/'predictions.npz' for base in sources.values() for case in cases]
    required += [base/'histology_control_predictions.npz' for base in sources.values()]
    if any(not path.exists() for path in required):raise FileNotFoundError('Run both descriptor baselines with control predictions first')
    hashes={str(path):sha(path) for path in required}
    args.output.mkdir(parents=True)
    (args.output/'experiment_config.json').write_text(json.dumps(CONFIG,indent=2))
    frozen={}
    # Complete and freeze EVERY model before opening ANY evaluation reference.
    for method,base in sources.items():
        for case in cases+['histology_control']:
            control=case=='histology_control'
            source=base/'histology_control_predictions.npz' if control else base/case/'predictions.npz'
            frozen[(method,case)]=freeze(source,args.output/method/case,control)
    print('All predictions frozen; current variants exactly reproduce saved retained pairs. Starting evaluation.',flush=True)
    rows=[];full=[]
    for (method,case),data in frozen.items():
        control=case=='histology_control'
        pair_dir=ROOT/'results/transforms/coronal_y8_demo' if control or case=='coronal_y8_demo' else ROOT/'results/random_pairs/y8_seed42'/case
        reference=json.loads((pair_dir/'truth/geometry.json').read_text())
        with np.load(pair_dir/'inputs.npz') as d:images={k:d[k] for k in d.files}
        mask=otsu_mask(images['mri'],images['mri_valid'],True)[0]['cleaned']
        if control:
            with np.load(pair_dir/'truth/original_histology.npz') as d:original=d['image'];valid=d['valid']
            mask=histology_hysteresis(original,valid,weak_fraction=.01)[0]['mask']
        reports=evaluate(data,reference,float(mask.sum()*.2**2),args.output/method/case,method,case)
        if not control:draw(images,data,reports,args.output/method/case/'final_comparison.png',method,case)
        for variant,r in reports.items():
            row=dict(method=method,case=case,variant=variant)
            for stage in ['input','retained']:
                for key in ['count','count_within_2mm','count_within_5mm','fraction_within_2mm','median_error_mm']:
                    row[f'{stage}_{key}']=r['stages'][stage][key]
            row['median_transform_error_mm']=r['transform_recovery']['median_whole_candidate_transform_error_mm'] if r['transform_recovery'] else None
            row['model_with_majority_pairs_beyond_5mm']=r['model_with_majority_pairs_beyond_5mm']
            rows.append(row);full.append(r)
            print(method,case,variant,'retained',row['retained_count'],'<=2mm',row['retained_count_within_2mm'],flush=True)
    with (args.output/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (args.output/'reports.json').write_text(json.dumps(full,indent=2))
    for path,digest in hashes.items():assert sha(Path(path))==digest
    (args.output/'checks.json').write_text(json.dumps(dict(all_saved_predictions_unchanged=True,
        current_retained_pairs_reproduced=True,all_variants_frozen_before_any_reference_read=True,
        input_sha256=hashes),indent=2))
    print('Saved',args.output)

if __name__=='__main__':main()
