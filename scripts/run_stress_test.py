"""Frozen image-only inference on an accepted stress case, then reference scoring."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import json
import time
import numpy as np
from landmarks import mind,relational
from landmarks.matching import noncollinear
from run_two_slice_validation import infer_frozen,evaluate_job,digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,default=ROOT/'results/validation/rotation_mismatch_stress_v1')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();study=args.study;out=args.output or study/'matching'
    if out.exists():raise FileExistsError('Choose a new --output directory')
    plan=json.loads((study/'study_plan.json').read_text())
    approval=json.loads((study/'preparation_user_review.json').read_text())
    if approval.get('accepted') is not True or approval.get('preview_sha256')!=digest(study/'pair/preview.png'):
        raise RuntimeError('Explicit user review of this exact preparation is required')
    for path,expected in plan['frozen_source_sha256'].items():
        if digest(ROOT/path)!=expected:raise ValueError(f'Frozen algorithm changed: {path}')
    if json.loads(json.dumps(mind.CONFIG))!=plan['mind_config'] or json.loads(json.dumps(relational.CONFIG))!=plan['relational_config']:
        raise ValueError('Frozen configuration mismatch')
    source=ROOT/plan['source_pair']/'images.npz';target=study/'pair/inputs.npz'
    hashes={str(p):digest(p) for p in [source,target]}
    out.mkdir(parents=True)
    (out/'run_provenance.json').write_text(json.dumps(dict(study_plan_sha256=digest(study/'study_plan.json'),
        preparation_review_sha256=digest(study/'preparation_user_review.json'),frozen_source_sha256=plan['frozen_source_sha256'],
        runner_sha256=digest(Path(__file__)),input_sha256=hashes),indent=2))
    jobs=[]
    for case,source_kind in [('stress','mri'),('histology_control','histology')]:
        with np.load(target) as d:images={k:d[k] for k in d.files}
        if source_kind=='histology':
            with np.load(source) as d:images['mri']=d['histology'];images['mri_valid']=d['histology_valid']
        started=time.perf_counter();pred=infer_frozen(images,source_kind);elapsed=time.perf_counter()-started
        path=out/case;path.mkdir()
        np.savez_compressed(path/'predictions.npz',source_candidates=pred['pa'],target_candidates=pred['pb'],
            input_pairs=pred['proposals'],retained_pairs=pred['retained'],similarity=pred['scores'],
            best_radius_angle_indices=pred['choices'],
            estimated_model=pred['model'] if pred['model'] is not None else np.empty((0,0)))
        np.savez_compressed(path/'search_support.npz',source=pred['masks']['mri'],target=pred['masks']['histology'])
        (path/'prediction.json').write_text(json.dumps(dict(case=case,candidate_counts=dict(source=len(pred['pa']),target=len(pred['pb'])),
            geometry=pred['status'],inference_seconds=elapsed,source_kind=source_kind,mind_config=mind.CONFIG,
            relational_config=relational.CONFIG,truth_access='None; both predictions saved before reference scoring'),indent=2))
        jobs.append(dict(slice_id='y8',case=case,pair_dir=study/'pair',input_path=target,input_sha256=hashes[str(target)],
            out=path,pred=pred,images=images,inference_seconds=elapsed))
        print(case,'frozen:',len(pred['proposals']),'proposals ->',len(pred['retained']),'retained',flush=True)
    print('BOTH predictions frozen. Starting reference-only evaluation.',flush=True)
    reports=[];rows=[];audits=[]
    for job in jobs:
        report=evaluate_job(job)
        report['caveats']=['One predeclared stress transformation of previously used y=8 slice; not a worst-case guarantee',
            'Released alignment approximate; random pairing is not conditioned on consensus selection',
            'All inference settings frozen; no retuning or selection of a replacement case']
        (job['out']/'evaluation.json').write_text(json.dumps(report,indent=2));reports.append(report)
        pred=job['pred'];pairs=pred['retained'];model=pred['model']
        if len(pairs):
            assert len(np.unique(pairs[:,0]))==len(pairs)==len(np.unique(pairs[:,1]))
            assert len(pairs)>=4 and noncollinear(pred['pa'][pairs[:,0]],pred['pb'][pairs[:,1]])
            assert np.isfinite(model).all() and np.linalg.det(model[:2,:2])>0
            residual=np.linalg.norm(pred['pa'][pairs[:,0]]@model[:2,:2].T+model[:2,2]-pred['pb'][pairs[:,1]],axis=1)
            assert (residual<6).all() and (pred['scores'][pairs[:,0],pairs[:,1]]>=.75).all()
        final=report['retained'];rec=report['transform_recovery']
        rows.append(dict(case=job['case'],source_candidates=len(pred['pa']),target_candidates=len(pred['pb']),
            proposals=report['input']['count'],proposal_sources_within_2mm=report['input']['unique_sources_within_2mm'],
            retained=final['count'],within_2mm=final['count_within_2mm'],within_5mm=final['count_within_5mm'],
            median_pair_error_mm=final['median_error_mm'],correct_hull_ratio=final['correct_hull_fraction_of_mri_mask_area'],
            median_transform_error_mm=rec['median_whole_candidate_transform_error_mm'] if rec else None,
            inference_seconds=job['inference_seconds']))
        audits.append(dict(case=job['case'],one_to_one_and_residual_rules_verified=True,
            max_model_residual_target_px=float(residual.max()) if len(pairs) else None))
        print(job['case'],f'{final["count_within_2mm"]}/{final["count"]} within 2 mm;',final['count_within_5mm'],'within 5 mm',flush=True)
    (out/'reports.json').write_text(json.dumps(reports,indent=2))
    with (out/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for path,expected in hashes.items():assert digest(Path(path))==expected
    for path,expected in plan['frozen_source_sha256'].items():assert digest(ROOT/path)==expected
    (out/'checks.json').write_text(json.dumps(dict(both_predictions_before_reference_geometry=True,
        input_and_algorithm_hashes_unchanged=True,configuration_and_thresholds_unchanged=True,assignment_audits=audits),indent=2))
    print('Saved',out)

if __name__=='__main__':main()
