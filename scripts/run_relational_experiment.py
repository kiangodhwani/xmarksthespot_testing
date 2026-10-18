"""Fixed MIND top-3/relational experiment: predictions first, reference scoring later."""
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
from landmarks.relational import CONFIG,topk_proposals,relational_filter,exhaustive_triples
from landmarks.proposal_ablation import select_variants
from landmarks.matching import geometric_filter
from landmarks.masking import otsu_mask,histology_hysteresis
from run_matching_baseline import evaluate_pairs,error_mm
from run_filter_ablation import recovery

VARIANTS=('mutual_ransac','mutual_relational','top3_relational')
LABELS=dict(mutual_ransac='Mutual + seeded RANSAC',mutual_relational='Mutual + relational geometry',
    top3_relational='Top 3 + relational geometry',mutual_exhaustive='Mutual + exhaustive three-point fits')


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze(source,out,control=False,diagnostic=False):
    with np.load(source) as d:saved={k:d[k] for k in d.files}
    pa=saved['source_candidates' if control else 'mri_candidates']
    pb=saved['target_candidates' if control else 'histology_candidates'];scores=saved['similarity']
    mutual=select_variants(scores)[1]['similarity_only'];top3=topk_proposals(scores,k=3,min_similarity=.75)
    methods=list(VARIANTS)+(['mutual_exhaustive'] if diagnostic else [])
    out.mkdir(parents=True);results={}
    for variant in methods:
        pairs=top3 if variant=='top3_relational' else mutual
        start=time.perf_counter()
        if variant=='mutual_ransac':
            model,keep,status=geometric_filter(pa[pairs[:,0]],pb[pairs[:,1]])
            retained=pairs[keep]
        elif variant=='mutual_exhaustive':model,retained,status=exhaustive_triples(pa,pb,pairs)
        else:model,retained,status=relational_filter(pa,pb,pairs,scores)
        elapsed=time.perf_counter()-start
        path=out/variant;path.mkdir()
        np.savez_compressed(path/'predictions.npz',source_candidates=pa,target_candidates=pb,
            input_pairs=pairs,retained_pairs=retained,estimated_model=model if model is not None else np.empty((0,0)))
        prediction=dict(variant=variant,description=LABELS[variant],geometry=status,inference_seconds=elapsed,
            source_prediction_sha256=digest(source),reference_access='None before all experiment predictions frozen')
        (path/'prediction.json').write_text(json.dumps(prediction,indent=2))
        results[variant]=dict(pairs=pairs,retained=retained,model=model,status=status,inference_seconds=elapsed)
        print('Frozen',out.name,variant,len(pairs),'proposals ->',len(retained),'pairs',f'{elapsed:.2f}s',flush=True)
    return dict(pa=pa,pb=pb,scores=scores,results=results)


def ranks(data,reference,out):
    pa,pb,scores=data['pa'],data['pb'],data['scores'];n,m=scores.shape
    ii,jj=np.indices((n,m));errors=error_mm(pa[ii.ravel()],pb[jj.ravel()],reference).reshape(n,m)
    rows=[];summary=dict(source_candidates=n,target_candidates=m)
    for tol in [2,5]:
        available=(errors<=tol).any(axis=1)
        usable=((errors<=tol)&(scores>=.75)).any(axis=1)
        summary[f'sources_with_any_target_within_{tol}mm']=int(available.sum())
        summary[f'sources_with_score_eligible_target_within_{tol}mm']=int(usable.sum())
        for k in [1,3,5]:
            pairs=topk_proposals(scores,k,.75)
            found=np.unique(pairs[errors[pairs[:,0],pairs[:,1]]<=tol,0])
            summary[f'top{k}_sources_with_partner_within_{tol}mm']=len(found)
    for i in range(n):
        order=np.argsort(-scores[i],kind='stable');legal=scores[i,order]>=.75
        close=(errors[i,order]<=2)&legal;indices=np.flatnonzero(close)
        rows.append(dict(source_index=i,any_target_within_2mm=bool((errors[i]<=2).any()),
            best_score_eligible_reference_close_rank=int(indices[0]+1) if len(indices) else None))
    with (out/'reference_partner_ranks.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'retrieval_diagnostic.json').write_text(json.dumps(summary,indent=2))
    return summary


def plot(data,images,reports,out,case):
    names=list(data['results']);fig,axs=plt.subplots(len(names),2,figsize=(14,len(names)*4.5))
    fig.subplots_adjust(left=.04,right=.98,top=.92,bottom=.08,hspace=.38,wspace=.20)
    for row,name in enumerate(names):
        item=data['results'][name];pairs=item['retained'];report=reports[name]
        show=pairs if len(pairs) else item['pairs']
        mode='FINAL' if len(pairs) else 'INPUT: no model'
        if len(show):
            order=np.argsort(-data['scores'][show[:,0],show[:,1]],kind='stable')[:12];shown=show[order]
        else:shown=show
        # Errors are evaluation-only; do not use them for display selection.
        from_report=report['retained_errors_mm'] if len(pairs) else report['input_errors_mm']
        errors=np.asarray(from_report)[order] if len(show) else []
        for ax,modality,title in zip(axs[row],['mri','histology'],['MRI','Transformed histology']):
            im=images[modality];valid=images[modality+'_valid'];lo,hi=np.percentile(im[valid],[1,99])
            ax.imshow(np.ma.array(im,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
            ax.set_xticks([]);ax.set_yticks([]);ax.set_facecolor('#ddd')
            ax.set_title(f'{title} — {LABELS[name]}',fontsize=10)
        for k,((i,j),error) in enumerate(zip(shown,errors),1):
            color='#168a35' if error<=2 else '#d28b00' if error<=5 else '#d12c3f'
            for ax,point in zip(axs[row],[data['pa'][i],data['pb'][j]]):
                ax.plot(*point,'o',color=color,ms=4)
                ax.annotate(str(k),point,xytext=(3,4),textcoords='offset points',fontsize=7,color=color,
                    bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=.3))
            fig.add_artist(ConnectionPatch(xyA=data['pa'][i],coordsA=axs[row,0].transData,
                xyB=data['pb'][j],coordsB=axs[row,1].transData,color=color,alpha=.55,lw=.8,clip_on=False))
        final=report['retained'];inp=report['input']
        axs[row,0].text(.5,-.10,f'Input: {inp["count"]} pairs; {inp["count_within_2mm"]} ≤2 mm',transform=axs[row,0].transAxes,ha='center',fontsize=9)
        axs[row,1].text(.5,-.10,f'{mode}; final {final["count_within_2mm"]}/{final["count"]} ≤2 mm; {final["count_within_5mm"]} ≤5 mm. Showing {len(shown)}.',transform=axs[row,1].transAxes,ha='center',fontsize=9)
    fig.suptitle(f'{case}: MIND proposals + geometry\nEvaluation-only colours: green ≤2 mm; amber 2–5 mm; red >5 mm.',fontsize=13)
    fig.text(.5,.018,'When no model returns, lines show INPUT proposals, not accepted matches. At most 12 highest-cosine pairs shown; metrics include ALL pairs.\nLabels are local to each row. Cosine is not confidence. Descriptor search unchanged: 3×3 radii; 12 rotations; no reference geometry in inference.\nApproximate reference, one slice/specimen. BigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=9)
    fig.savefig(out/'matches.png',dpi=120);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=ROOT/'results/matching/mind_v1')
    p.add_argument('--output',type=Path,default=ROOT/'results/matching/relational_v1')
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('Choose a new output directory')
    cases=['coronal_y8_demo']+[f'pair_{i:03d}' for i in range(1,6)]+['histology_control']
    sources={case:args.input/('histology_control_predictions.npz' if case=='histology_control' else f'{case}/predictions.npz') for case in cases}
    if any(not s.is_file() for s in sources.values()):raise FileNotFoundError('Complete MIND saved predictions required')
    hashes={case:digest(path) for case,path in sources.items()}
    args.output.mkdir(parents=True)
    (args.output/'experiment_config.json').write_text(json.dumps(CONFIG,indent=2))
    frozen={case:freeze(path,args.output/case,case=='histology_control',case=='coronal_y8_demo') for case,path in sources.items()}
    # Compare seeded control against previous ablation before reference access.
    for case,data in frozen.items():
        previous=ROOT/'results/matching/filter_ablation_v1/mind'/case/'similarity_only/predictions.npz'
        if previous.is_file():
            with np.load(previous) as d:
                np.testing.assert_array_equal(data['results']['mutual_ransac']['retained'],d['retained_pairs'])
                actual=data['results']['mutual_ransac']['model']
                np.testing.assert_array_equal(actual if actual is not None else np.empty((0,0)),d['estimated_model'])
    print('All predictions frozen. Starting reference-only evaluation.',flush=True)
    rows=[];all_reports=[];retrieval=[]
    for case,data in frozen.items():
        out=args.output/case;control=case=='histology_control'
        pair_dir=ROOT/'results/transforms/coronal_y8_demo' if case=='coronal_y8_demo' or control else ROOT/'results/random_pairs/y8_seed42'/case
        reference=json.loads((pair_dir/'truth/geometry.json').read_text())
        with np.load(pair_dir/'inputs.npz') as d:images={k:d[k] for k in d.files}
        mask=otsu_mask(images['mri'],images['mri_valid'],True)[0]['cleaned']
        if control:
            with np.load(pair_dir/'truth/original_histology.npz') as d:im=d['image'];valid=d['valid']
            mask=histology_hysteresis(im,valid,weak_fraction=.01)[0]['mask']
        retrieval.append(dict(case=case,**ranks(data,reference,out)))
        reports={}
        for name,item in data['results'].items():
            pairs=item['pairs'];retained=item['retained'];pa,pb=data['pa'],data['pb']
            errors=error_mm(pa[pairs[:,0]],pb[pairs[:,1]],reference)
            # Input proposals are many-to-many: cannot use the evaluator's unique
            # random-pair comparator with count > number of source candidates.
            inp=dict(count=len(pairs),count_within_2mm=int((errors<=2).sum()),count_within_5mm=int((errors<=5).sum()),
                unique_sources_within_2mm=len(np.unique(pairs[errors<=2,0])),chance_note='No unique-pair random comparator for many-to-many proposal counts')
            final,final_errors=evaluate_pairs(retained,pa,pb,reference,float(mask.sum()*.2**2),np.random.default_rng(7))
            rec=recovery(item['model'],pa,reference)
            report=dict(case=case,variant=name,input=inp,retained=final,input_errors_mm=errors.tolist(),
                retained_errors_mm=final_errors.tolist(),transform_recovery=rec,geometry=item['status'],inference_seconds=item['inference_seconds'])
            (out/name/'evaluation.json').write_text(json.dumps(report,indent=2));reports[name]=report;all_reports.append(report)
            rows.append(dict(case=case,variant=name,input_count=len(pairs),input_within_2mm=inp['count_within_2mm'],
                input_unique_sources_within_2mm=inp['unique_sources_within_2mm'],retained=final['count'],within_2mm=final['count_within_2mm'],
                within_5mm=final['count_within_5mm'],median_pair_error_mm=final['median_error_mm'],
                median_transform_error_mm=rec['median_whole_candidate_transform_error_mm'] if rec else None,
                inference_seconds=item['inference_seconds']))
            print(case,name,f'{final["count_within_2mm"]}/{final["count"]} <=2mm',flush=True)
        if not control:plot(data,images,reports,out,case)
    for name,records in [('summary',rows),('retrieval',retrieval)]:
        with (args.output/f'{name}.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (args.output/'reports.json').write_text(json.dumps(all_reports,indent=2))
    assert all(digest(path)==hashes[case] for case,path in sources.items())
    (args.output/'checks.json').write_text(json.dumps(dict(predictions_frozen_before_any_reference_access=True,
        unchanged_source_hashes=hashes,seeded_mutual_control_reproduced=True),indent=2))
    print('Saved',args.output)

if __name__=='__main__':main()
