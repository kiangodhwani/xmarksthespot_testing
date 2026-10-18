"""Frozen MIND top-3 validation on two visually accepted additional slices.

Inference receives image arrays and validity only. Every prediction, including
controls, is frozen before any reference geometry is opened for evaluation.
"""
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
import numpy as np
from landmarks import mind,relational
from landmarks.candidates import detect_candidates
from landmarks.masking import otsu_mask,histology_hysteresis
from landmarks.planes import pixel_to_world
from run_matching_baseline import evaluate_pairs,error_mm
from run_filter_ablation import recovery


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def infer_frozen(images,source_kind='mri'):
    """Pure image/support input boundary; no paths, y values, transforms or seeds."""
    masks={};points={};banks={};usable={};detector={}
    for name,kind in [('mri',source_kind),('histology','histology')]:
        image,valid=images[name],images[name+'_valid']
        if kind=='mri':mask=otsu_mask(image,valid,True)[0]['cleaned']
        else:mask=histology_hysteresis(image,valid,weak_fraction=.01)[0]['mask']
        masks[name]=mask
        points[name],detector[name]=detect_candidates(image,valid,mask)
        banks[name],usable[name]=mind.describe(image,valid,points[name][:,:2])
    scores,choices=mind.similarity(banks['mri'],usable['mri'],banks['histology'],usable['histology'],return_choices=True)
    pa,pb=points['mri'][:,:2],points['histology'][:,:2]
    proposals=relational.topk_proposals(scores,k=3,min_similarity=.75)
    model,retained,status=relational.relational_filter(pa,pb,proposals,scores,residual_px=6.)
    return dict(pa=pa,pb=pb,scores=scores,choices=choices,proposals=proposals,model=model,
        retained=retained,status=status,masks=masks,detector=detector,
        usable_counts={name:int(v.any(axis=tuple(range(1,v.ndim))).sum()) for name,v in usable.items()})


def freeze_job(study,slice_id,case,out):
    folder=study/slice_id;control=case=='histology_control'
    pair_dir=folder/'fixed' if control or case=='fixed' else folder/'random_seed42'/case
    input_path=pair_dir/'inputs.npz';input_hash=digest(input_path)
    with np.load(input_path) as d:images={k:d[k] for k in d.files}
    source_hash=None
    if control:
        original=folder/'pair/images.npz';source_hash=digest(original)
        with np.load(original) as d:images['mri']=d['histology'];images['mri_valid']=d['histology_valid']
    start=time.perf_counter();pred=infer_frozen(images,'histology' if control else 'mri');elapsed=time.perf_counter()-start
    out.mkdir(parents=True)
    np.savez_compressed(out/'predictions.npz',source_candidates=pred['pa'],target_candidates=pred['pb'],
        input_pairs=pred['proposals'],retained_pairs=pred['retained'],similarity=pred['scores'],
        best_radius_angle_indices=pred['choices'],
        estimated_model=pred['model'] if pred['model'] is not None else np.empty((0,0)))
    np.savez_compressed(out/'search_support.npz',source=pred['masks']['mri'],target=pred['masks']['histology'])
    metadata=dict(slice_id=slice_id,case=case,source_modality='original_histology' if control else 'mri',
        candidate_counts=dict(source=len(pred['pa']),target=len(pred['pb'])),usable_counts=pred['usable_counts'],
        detector=pred['detector'],geometry=pred['status'],inference_seconds=elapsed,
        algorithm='Frozen MIND-style + source top3 + one-to-one relational similarity consensus',
        mind_config=mind.CONFIG,relational_config=relational.CONFIG,
        input_sha256=input_hash,control_source_sha256=source_hash,
        truth_access='None; all 14 predictions saved before any reference geometry is read')
    (out/'prediction.json').write_text(json.dumps(metadata,indent=2))
    assert digest(input_path)==input_hash
    if control:assert digest(folder/'pair/images.npz')==source_hash
    print(slice_id,case,'frozen:',len(pred['proposals']),'proposals ->',len(pred['retained']),'retained',f'{elapsed:.2f}s',flush=True)
    return dict(slice_id=slice_id,case=case,pair_dir=pair_dir,input_path=input_path,
        input_sha256=input_hash,control_source_sha256=source_hash,out=out,pred=pred,images=images,inference_seconds=elapsed)


def evaluate_job(job):
    pred=job['pred'];pa,pb=pred['pa'],pred['pb'];out=job['out']
    reference=json.loads((job['pair_dir']/'truth/geometry.json').read_text())
    final,errors=evaluate_pairs(pred['retained'],pa,pb,reference,
        float(pred['masks']['mri'].sum()*.2**2),np.random.default_rng(7))
    pairs=pred['proposals'];input_errors=error_mm(pa[pairs[:,0]],pb[pairs[:,1]],reference)
    inp=dict(count=len(pairs),count_within_2mm=int((input_errors<=2).sum()),count_within_5mm=int((input_errors<=5).sum()),
        unique_sources_within_2mm=len(np.unique(pairs[input_errors<=2,0])),
        chance_note='Many-to-many proposals; unique-pair chance comparator only applies to final one-to-one matches')
    report=dict(slice_id=job['slice_id'],case=job['case'],input=inp,retained=final,
        input_errors_mm=input_errors.tolist(),retained_errors_mm=errors.tolist(),
        transform_recovery=recovery(pred['model'],pa,reference),geometry=pred['status'],
        inference_seconds=job['inference_seconds'],
        caveats=['Two new slices of same specimen; transform settings reused; approximate released alignment',
            'Final random pairing is not conditioned on consensus selection and is not a p-value'])
    (out/'evaluation.json').write_text(json.dumps(report,indent=2))
    with (out/'all_proposals_evaluation.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['source_index','target_index','cosine_similarity','reference_error_mm','retained'])
        retained={tuple(p) for p in pred['retained']}
        for pair,error in zip(pairs,input_errors):
            i,j=pair;w.writerow([i,j,pred['scores'][i,j],error,tuple(pair) in retained])
    control=job['case']=='histology_control'
    with (out/'landmarks.csv').open('w',newline='') as f:
        w=csv.writer(f)
        prefix='histology_original' if control else 'mri'
        w.writerow(['id','name',prefix+'_x_mm',prefix+'_y_mm',prefix+'_z_mm','hist_col_px','hist_row_px','score'])
        for k,(i,j) in enumerate(pred['retained'],1):
            world=pixel_to_world(pa[i:i+1],np.asarray(reference['original_pixel_to_world_ras_mm']))[0]
            w.writerow([k,'candidate_'+str(i),*world,*pb[j],float(pred['scores'][i,j])])
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',type=Path,default=ROOT/'results/validation/two_slices_v1')
    p.add_argument('--output',type=Path)
    args=p.parse_args();study=args.study
    out=args.output or study/'matching'
    if out.exists():raise FileExistsError('Choose a new --output directory')
    plan=json.loads((study/'study_plan.json').read_text())
    for path,expected in plan['frozen_source_sha256'].items():
        if digest(ROOT/path)!=expected:raise ValueError(f'Frozen algorithm changed: {path}')
    # Authorisation is recorded separately from the immutable predeclared plan.
    approval=study/'preparation_user_review.json'
    if not approval.is_file() or json.loads(approval.read_text()).get('accepted') is not True:
        raise RuntimeError('Prepared images require explicit user visual review before matching')
    if json.loads(json.dumps(mind.CONFIG))!=plan['mind_config'] or json.loads(json.dumps(relational.CONFIG))!=plan['relational_config']:
        raise ValueError('Configuration differs from predeclared study plan')
    out.mkdir(parents=True)
    (out/'run_provenance.json').write_text(json.dumps(dict(study_plan_sha256=digest(study/'study_plan.json'),
        preparation_review_sha256=digest(approval),frozen_source_sha256=plan['frozen_source_sha256'],
        runner_sha256=digest(Path(__file__)),thresholds_mm=[2.,5.],configuration_unchanged=True),indent=2))
    jobs=[]
    for slice_id in ['y_m20','y35']:
        for case in ['fixed']+[f'pair_{i:03d}' for i in range(1,6)]+['histology_control']:
            jobs.append(freeze_job(study,slice_id,case,out/slice_id/case))
    print('ALL 14 predictions frozen. Beginning reference-only evaluation.',flush=True)
    reports=[evaluate_job(job) for job in jobs]
    (out/'reports.json').write_text(json.dumps(reports,indent=2))
    rows=[]
    for r,job in zip(reports,jobs):
        final=r['retained'];rec=r['transform_recovery']
        rows.append(dict(slice_id=r['slice_id'],case=r['case'],source_candidates=len(job['pred']['pa']),target_candidates=len(job['pred']['pb']),
            proposals=r['input']['count'],proposal_sources_within_2mm=r['input']['unique_sources_within_2mm'],
            retained=final['count'],within_2mm=final['count_within_2mm'],within_5mm=final['count_within_5mm'],
            fraction_within_2mm=final['fraction_within_2mm'],median_pair_error_mm=final['median_error_mm'],
            correct_hull_fraction=final['correct_hull_fraction_of_mri_mask_area'],
            median_transform_error_mm=rec['median_whole_candidate_transform_error_mm'] if rec else None,
            inference_seconds=r['inference_seconds']))
        assert digest(job['input_path'])==job['input_sha256']
        if job['control_source_sha256'] is not None:
            assert digest(study/job['slice_id']/'pair/images.npz')==job['control_source_sha256']
        print(r['slice_id'],r['case'],f'{final["count_within_2mm"]}/{final["count"]} within2mm',flush=True)
    with (out/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for path,expected in plan['frozen_source_sha256'].items():assert digest(ROOT/path)==expected
    (out/'checks.json').write_text(json.dumps(dict(all_14_predictions_before_reference_evaluation=True,
        all_inputs_and_frozen_algorithm_hashes_unchanged=True,configuration_and_thresholds_unchanged=True),indent=2))
    print('Saved',out)

if __name__=='__main__':main()
