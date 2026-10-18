"""Export saved matches and success metrics without rerunning or changing inference."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import csv
import hashlib
import json
import shutil
import tempfile
import zipfile
import numpy as np
from landmarks.planes import pixel_to_world
from landmarks.matching import noncollinear
from run_matching_baseline import error_mm,spread
from run_filter_ablation import recovery

HEADER=['id','name','mri_x_mm','mri_y_mm','mri_z_mm','hist_col_px','hist_row_px','score']
CONTROL_HEADER=['id','name','histology_original_x_mm','histology_original_y_mm','histology_original_z_mm','hist_col_px','hist_row_px','score']


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def identify(path,results):
    parts=path.relative_to(results).parts
    control='histology_control' in parts or path.name=='histology_control_predictions.npz'
    record=dict(prediction=path,control=control,primary=False,slice_id='y8',descriptor='mind',
        score_source=path,case='histology_control' if control else None)
    if parts[0]=='matching':
        run=parts[1];record['experiment']=run
        if run in ['polar_baseline','polar_backend_check','mind_v1']:
            case=record['case'] or parts[2];record.update(case=case,variant='current_ratio',descriptor='mind' if run=='mind_v1' else 'polar')
            record['evaluation']=path.parent/'histology_control.json' if control else path.parent/'evaluation/report.json'
        elif run=='filter_ablation_v1':
            record.update(descriptor=parts[2],case=parts[3],variant=parts[4],evaluation=path.parent/'evaluation.json')
        elif run=='relational_v1':
            record.update(case=parts[2],variant=parts[3],evaluation=path.parent/'evaluation.json')
            record['primary']=parts[3]=='top3_relational' and not control
            record['score_source']=results/'matching/mind_v1'/('histology_control_predictions.npz' if control else parts[2]+'/predictions.npz')
        else:raise ValueError(f'Unrecognised matching study: {path}')
        case=record['case']
        record['pair_dir']=results/'transforms/coronal_y8_demo' if control or case=='coronal_y8_demo' else results/'random_pairs/y8_seed42'/case
        record['case_id']='fixed' if case=='coronal_y8_demo' else case
    elif parts[0]=='validation' and parts[1]=='two_slices_v1' and parts[2]=='matching':
        record.update(experiment='two_slices_v1',slice_id=parts[3],case=parts[4],case_id=parts[4],variant='top3_relational',primary=not control,evaluation=path.parent/'evaluation.json')
        base=results/'validation/two_slices_v1'/parts[3]
        record['pair_dir']=base/'fixed' if control or record['case']=='fixed' else base/'random_seed42'/record['case']
    elif parts[0]=='validation' and parts[1]=='rotation_mismatch_stress_v1' and parts[2]=='matching':
        record.update(experiment='rotation_mismatch_stress_v1',case=parts[3],case_id='histology_control' if control else 'rotation_mismatch_stress',
            variant='top3_relational',primary=not control,evaluation=path.parent/'evaluation.json',pair_dir=results/'validation/rotation_mismatch_stress_v1/pair')
    else:raise ValueError(f'Unrecognised prediction artifact, do not silently omit: {path}')
    record['pair_id']=record['slice_id']+'/'+record['case_id']
    if record['primary']:record['output_relative']=Path('final_pipeline')/record['pair_id']
    else:
        suffix=path.parent.relative_to(results)
        if path.name!='predictions.npz':suffix=suffix/'histology_control'
        record['output_relative']=Path('controls' if control else 'historical')/suffix
    return record


def normalize_prediction(record):
    with np.load(record['prediction']) as d:saved={k:d[k] for k in d.files}
    pa=saved['source_candidates' if 'source_candidates' in saved else 'mri_candidates']
    pb=saved['target_candidates' if 'target_candidates' in saved else 'histology_candidates']
    pairs=np.asarray(saved['retained_pairs'],int).reshape(-1,2)
    source=record['score_source']
    if 'similarity' in saved:scores=saved['similarity']
    else:
        with np.load(source) as d:
            scores=d['similarity']
            np.testing.assert_array_equal(pa,d['source_candidates' if record['control'] else 'mri_candidates'])
            np.testing.assert_array_equal(pb,d['target_candidates' if record['control'] else 'histology_candidates'])
    if 'estimated_model' in saved:model=saved['estimated_model'] if saved['estimated_model'].size else None
    else:
        metadata=json.loads((record['prediction'].parent/'prediction.json').read_text())
        value=metadata['estimated_mri_to_histology'];model=np.asarray(value) if value is not None else None
    proposed=saved.get('input_pairs',saved.get('distinctive_pairs',np.empty((0,2),int)))
    assert scores.shape==(len(pa),len(pb))
    assert len(np.unique(pairs[:,0]))==len(pairs)==len(np.unique(pairs[:,1]))
    if len(pairs):
        assert (pairs>=0).all() and (pairs[:,0]<len(pa)).all() and (pairs[:,1]<len(pb)).all()
        assert np.isfinite(scores[pairs[:,0],pairs[:,1]]).all()
    return pa,pb,pairs,scores,model,proposed


def landmark_rows(pa,pb,pairs,scores,matrix):
    """Export predictions unchanged; reference matrix only converts source units."""
    world=pixel_to_world(pa[pairs[:,0]],np.asarray(matrix))
    return [[k,'candidate_'+str(int(i)),*map(float,w),*map(float,pb[j]),float(scores[i,j])]
        for k,((i,j),w) in enumerate(zip(pairs,world),1)]


def compare_existing(path,header,rows):
    if not path.is_file():return False
    with path.open(newline='') as f:actual=list(csv.reader(f))
    if actual[0]!=header or len(actual)!=len(rows)+1:raise ValueError(f'Existing CSV schema/count mismatch: {path}')
    for a,b in zip(actual[1:],rows):
        if int(a[0])!=b[0] or a[1]!=b[1]:raise ValueError(f'Existing IDs differ: {path}')
        np.testing.assert_allclose(np.asarray(a[2:],float),np.asarray(b[2:],float),rtol=0,atol=1e-10)
    return True


def score_metrics(pa,pb,pairs,reference,model,saved):
    errors=error_mm(pa[pairs[:,0]],pb[pairs[:,1]],reference)
    n=len(pairs);correct=pairs[errors<=2]
    world=pixel_to_world(pa[correct[:,0]],np.asarray(reference['original_pixel_to_world_ras_mm']))
    geometry=spread(world[:,[0,2]])
    old=saved.get('stages',{}).get('retained')
    if old is None and isinstance(saved.get('retained'),dict):old=saved['retained']
    result=dict(count=n,count_within_2mm=int((errors<=2).sum()),count_within_5mm=int((errors<=5).sum()),
        fraction_within_2mm=float((errors<=2).mean()) if n else None,
        fraction_within_5mm=float((errors<=5).mean()) if n else None,
        median_error_mm=float(np.median(errors)) if n else None,
        mean_error_mm=float(errors.mean()) if n else None,max_error_mm=float(errors.max()) if n else None,
        correct_within_2mm_spread=geometry,
        correct_pairs_ge4_noncollinear=bool(len(correct)>=4 and noncollinear(pa[correct[:,0]],pb[correct[:,1]])))
    if old:
        for key in ['count','count_within_2mm','count_within_5mm']:
            if result[key]!=old[key]:raise ValueError(f'Saved count mismatch: {key}')
        for key in ['fraction_within_2mm','fraction_within_5mm','median_error_mm']:
            if result[key] is None:assert old[key] is None
            else:np.testing.assert_allclose(result[key],old[key],rtol=0,atol=1e-9)
        for key,value in geometry.items():np.testing.assert_allclose(value,old['correct_within_2mm_spread'][key],rtol=0,atol=1e-8)
        result['correct_hull_area_divided_by_source_mask_area']=old.get('correct_hull_fraction_of_mri_mask_area')
    elif isinstance(saved.get('retained'),int):
        assert n==saved['retained'] and result['count_within_2mm']==saved['retained_within_2mm']
    chance=old.get('random_candidate_pairing') if old else None
    if n and chance is None:
        samples=[];rng=np.random.default_rng(7)
        for _ in range(200):
            a=rng.choice(len(pa),n,replace=False);b=rng.choice(len(pb),n,replace=False)
            e=error_mm(pa[a],pb[b],reference);samples.append([(e<=2).mean(),(e<=5).mean()])
        means=np.mean(samples,axis=0)
        chance=dict(trials=200,mean_within_2mm=float(means[0]),mean_within_5mm=float(means[1]),
            seed=7,note='Same-count independent unique candidate pairing; not conditioned on consensus selection')
        chance_origin='computed_for_missing_legacy_control_metric'
    else:chance_origin='preserved_saved_evaluation' if chance else 'not_applicable_empty_output'
    result['random_candidate_pairing']=chance;result['chance_metric_origin']=chance_origin
    result['observed_exceeds_random_mean_2mm']=result['fraction_within_2mm']>chance['mean_within_2mm'] if chance else None
    result['observed_exceeds_random_mean_5mm']=result['fraction_within_5mm']>chance['mean_within_5mm'] if chance else None
    rec=recovery(model,pa,reference)
    if saved.get('transform_recovery') is not None:
        for key,value in rec.items():np.testing.assert_allclose(value,saved['transform_recovery'][key],rtol=0,atol=1e-9)
    theta=lambda m:float(np.rad2deg(np.arctan2(m[1,0],m[0,0])))
    scale=lambda m:float(np.sqrt(np.linalg.det(m[:2,:2])))
    planted=np.asarray(reference['original_to_transformed'])
    result['transform_recovery']=rec
    result['transform_parameters']=dict(planted_clockwise_degrees=theta(planted),planted_scale=scale(planted),
        estimated_clockwise_degrees=theta(model) if model is not None else None,
        estimated_scale=scale(model) if model is not None else None,
        estimated_source_to_target_matrix=model.tolist() if model is not None else None)
    return result,errors


def fingerprint(record,pa,pb,pairs,scores,model,proposed):
    h=hashlib.sha256((record['pair_id']+str(record['control'])).encode())
    for array in [pa,pb,pairs,scores,np.empty((0,0)) if model is None else model,proposed]:
        v=np.asarray(array,dtype='<f8');h.update(str(v.shape).encode());h.update(v.tobytes())
    return h.hexdigest()


def write_csv(path,header,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/exports/landmark_release_v1')
    args=p.parse_args();requested=args.output;results=ROOT/'results'
    if requested.exists() or requested.with_suffix('.zip').exists():raise FileExistsError('Choose a fresh output directory; historical artifacts are immutable')
    requested.parent.mkdir(parents=True,exist_ok=True)
    out=Path(tempfile.mkdtemp(prefix=requested.name+'.building-',dir=requested.parent))
    paths=sorted(results.rglob('*predictions.npz'))
    records=[identify(path,results) for path in paths]
    if len({str(r['output_relative']) for r in records})!=len(records):raise ValueError('Output collision')
    summary=[];inventory=[];hashes={};seen={};reused=0
    for record in records:
        path=record['prediction'];reference_path=record['pair_dir']/'truth/geometry.json'
        for source in [path,record['score_source'],record['evaluation'],reference_path]:hashes[str(source)]=digest(source)
        reference=json.loads(reference_path.read_text());saved=json.loads(record['evaluation'].read_text())
        pa,pb,pairs,scores,model,proposed=normalize_prediction(record)
        rows=landmark_rows(pa,pb,pairs,scores,reference['original_pixel_to_world_ras_mm'])
        header=CONTROL_HEADER if record['control'] else HEADER
        destination=out/record['output_relative'];destination.mkdir(parents=True)
        existing=path.parent/'landmarks.csv'
        if compare_existing(existing,header,rows):
            hashes[str(existing)]=digest(existing);shutil.copyfile(existing,destination/'landmarks.csv');reused+=1;origin='copied_existing_verified_csv'
        else:write_csv(destination/'landmarks.csv',header,rows);origin='exported_from_saved_prediction'
        metrics,errors=score_metrics(pa,pb,pairs,reference,model,saved)
        evaluation=destination/'evaluation';evaluation.mkdir()
        (evaluation/'metrics.json').write_text(json.dumps(metrics,indent=2))
        shutil.copyfile(reference_path,evaluation/'reference_geometry.json')
        shutil.copyfile(record['evaluation'],evaluation/'original_report.json')
        write_csv(evaluation/'landmark_errors.csv',['id','source_candidate_index','target_candidate_index','reference_error_mm','within_2mm','within_5mm'],
            [[k,int(i),int(j),float(e),bool(e<=2),bool(e<=5)] for k,((i,j),e) in enumerate(zip(pairs,errors),1)])
        semantic=fingerprint(record,pa,pb,pairs,scores,model,proposed)
        duplicate=seen.get(semantic);seen.setdefault(semantic,str(record['output_relative']))
        chance=metrics['random_candidate_pairing'];rec=metrics['transform_recovery'];tp=metrics['transform_parameters']
        row=dict(experiment=record['experiment'],descriptor=record['descriptor'],variant=record['variant'],
            slice_id=record['slice_id'],pair_id=record['pair_id'],source_modality='original_histology' if record['control'] else 'mri',
            final_pipeline=record['primary'],control=record['control'],
            status='model_returned' if model is not None else 'no_consensus',count=len(rows),
            count_within_2mm=metrics['count_within_2mm'],count_within_5mm=metrics['count_within_5mm'],
            fraction_within_2mm=metrics['fraction_within_2mm'],fraction_within_5mm=metrics['fraction_within_5mm'],
            median_error_mm=metrics['median_error_mm'],max_error_mm=metrics['max_error_mm'],
            random_fraction_within_2mm=chance['mean_within_2mm'] if chance else None,
            random_fraction_within_5mm=chance['mean_within_5mm'] if chance else None,
            random_trials=chance['trials'] if chance else None,chance_metric_origin=metrics['chance_metric_origin'],
            observed_exceeds_random_mean_2mm=metrics['observed_exceeds_random_mean_2mm'],
            observed_exceeds_random_mean_5mm=metrics['observed_exceeds_random_mean_5mm'],
            correct_hull_area_mm2=metrics['correct_within_2mm_spread']['hull_area_mm2'],
            correct_minor_major_spread_ratio=metrics['correct_within_2mm_spread']['minor_major_spread_ratio'],
            correct_pairs_ge4_noncollinear=metrics['correct_pairs_ge4_noncollinear'],
            planted_clockwise_degrees=tp['planted_clockwise_degrees'],estimated_clockwise_degrees=tp['estimated_clockwise_degrees'],
            planted_scale=tp['planted_scale'],estimated_scale=tp['estimated_scale'],
            rotation_absolute_error_degrees=rec['rotation_absolute_error_degrees'] if rec else None,
            scale_absolute_error=rec['scale_absolute_error'] if rec else None,
            median_candidate_transform_error_mm=rec['median_whole_candidate_transform_error_mm'] if rec else None,
            prediction_duplicate_of=duplicate,csv=str(record['output_relative']/'landmarks.csv'))
        summary.append(row)
        provenance=dict(source_prediction=str(path.relative_to(ROOT)),source_prediction_sha256=hashes[str(path)],
            score_source=str(record['score_source'].relative_to(ROOT)),reference_geometry=str(reference_path.relative_to(ROOT)),
            evaluated_image_pair=str((record['pair_dir']/'inputs.npz').relative_to(ROOT)),
            source_evaluation=str(record['evaluation'].relative_to(ROOT)),csv_origin=origin,
            csv_sha256=digest(destination/'landmarks.csv'),source_modality=row['source_modality'],
            original_pairs=[[int(i),int(j)] for i,j in pairs],prediction_duplicate_of=duplicate,
            note='Rows are algorithm-retained predictions, including reference-incorrect pairs; no truth-conditioned filtering')
        (destination/'provenance.json').write_text(json.dumps(provenance,indent=2));inventory.append({**row,**provenance})
    header=list(summary[0])
    for name,rows in [('all_tests_summary',summary),('final_pipeline_summary',[r for r in summary if r['final_pipeline']]),
        ('historical_summary',[r for r in summary if not r['control'] and not r['final_pipeline']]),('controls_summary',[r for r in summary if r['control']])]:
        write_csv(out/(name+'.csv'),header,[[r[k] for k in header] for r in rows])
    # Original polar control has metrics but no frozen prediction artifact. A
    # verified exact replay provides its coordinates; label as alias, never invent.
    original=results/'matching/polar_baseline/histology_control.json'
    replay=results/'matching/polar_backend_check/histology_control.json'
    aliases=[]
    if original.exists() and replay.exists():
        assert json.loads(original.read_text())==json.loads(replay.read_text())
        hashes[str(original)]=digest(original);hashes[str(replay)]=digest(replay)
        aliases.append(dict(source=str(original.relative_to(ROOT)),status='metrics_only_original_with_exact_replay',
            csv='controls/matching/polar_backend_check/histology_control/landmarks.csv',
            reason='Original saved control report exactly equals backend replay; coordinates from replay explicitly, no invented original predictions'))
    for path,sha in hashes.items():assert digest(Path(path))==sha
    primary=[r for r in summary if r['final_pipeline']]
    aggregate=dict(image_pairs=len(primary),retained=sum(r['count'] for r in primary),
        within_2mm=sum(r['count_within_2mm'] for r in primary),within_5mm=sum(r['count_within_5mm'] for r in primary),
        scope='Three anatomical planes, one specimen; reused transformations; not independent subjects or all possible transformations')
    total=aggregate['retained']
    for threshold in [2,5]:
        aggregate[f'fraction_within_{threshold}mm']=aggregate[f'within_{threshold}mm']/total if total else None
        aggregate[f'same_count_random_expected_fraction_within_{threshold}mm']=sum(r['count']*(r[f'random_fraction_within_{threshold}mm'] or 0) for r in primary)/total if total else None
    manifest=dict(records=inventory,metrics_only_aliases=aliases,final_pipeline=aggregate,
        prepared_but_not_matched=['results/validation/boundary_stress_v1/pair'],
        exclusions='Candidate detection/repeatability and oracle-location descriptor diagnostics are not accepted match lists and are not exported as registration landmarks',
        source_hashes_unchanged=hashes)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    (out/'checks.json').write_text(json.dumps(dict(saved_prediction_artifacts=len(paths),csv_files=len(summary),
        final_pipeline_csvs=len(primary),historical_csvs=sum(not r['control'] and not r['final_pipeline'] for r in summary),
        control_csvs=sum(r['control'] for r in summary),existing_csvs_verified_and_reused=reused,
        all_saved_metrics_verified=True,all_input_hashes_unchanged=True,no_matching_or_refitting_run=True,
        primary_aggregate=aggregate),indent=2))
    readme=f'''# Landmark outputs and success metrics

Start with **final_pipeline_summary.csv** and **final_pipeline/**: {len(primary)}
image pairs tested with the frozen MIND top-3 + relational pipeline. Each pair
has one landmarks.csv in the exact requested eight-column MRI/histology schema.
MRI coordinates are RAS world mm; histology coordinates are (column,row) pixel
CENTRES in the saved transformed section, including its canvas offset. Names
are candidate IDs, not verified anatomical labels. Score is max-bank cosine,
not a probability. All retained predictions are exported; reference errors
never select CSV rows.

Final pipeline: {aggregate['within_2mm']}/{total} retained pairs within 2 mm;
{aggregate['within_5mm']}/{total} within 5 mm. This pools repeated transformations
of three slices of ONE specimen. See per-case errors, count and spread; pooled
precision alone does not establish registration quality/generalisation.

- all_tests_summary.csv: all {len(summary)} saved prediction artifacts, including
  baselines, ablations, repeated runs, final tests, failures and controls.
- historical_summary.csv / historical/: earlier MRI/histology variants.
- controls_summary.csv / controls/: histology-only diagnostics, with source
  columns named histology_original_x_mm/y_mm/z_mm. They are intentionally NOT
  falsely labelled as MRI coordinates or mixed into registration outputs.
- prediction_duplicate_of flags exact repeated prediction semantics (including
  candidate/score/proposal/model arrays). Do not count repeated variants as new
  anatomical tests. An empty no-consensus run still has a header-only CSV;
  accuracy fractions and transform errors are blank, not 100% or zero error.
- Each pair's evaluation/ folder contains metrics.json, original_report.json,
  evaluation-only reference_geometry.json and per-row landmark_errors.csv.
- provenance.json maps CSV row IDs to original candidate indices and files.
- manifest.json inventories every artifact, aliases, exclusions and source hashes.

Metrics: retained count, fractions <=2/5 mm, median/max error, correct-pair convex
hull area (mm²), minor/major spread and >=4 noncollinear correct-pair flag;
planted/estimated rotation and scale and their absolute errors; candidate-field
median transform recovery error. Saved chance metrics are preserved; missing
legacy-control chance metrics are computed from frozen candidate pools using
200 same-count unique random pairings and seed7. Chance is NOT conditioned on
consensus selection or a p-value. observed_exceeds_random_mean is a descriptive
comparison, not a significance test. Convex-hull/mask area ratios in original
reports are spread ratios and can exceed1; they are not tissue coverage.

Original polar-baseline control has no saved coordinate artifact. Its report
exactly matches a saved backend replay; manifest records that replay as an
explicit alias, rather than claiming recovered original coordinates. Preliminary
candidate/descriptor oracle diagnostics are not final match predictions. The
-30° boundary preparation was never matched and therefore has no invented CSV.

No matching, fitting, threshold changes or source modifications were performed.
World coordinates and success metrics were independently checked against saved
predictions/reference mappings and existing CSVs/reports. Same source specimen
and approximate released alignment remain limitations. Derived BigBrain data:
Amunts et al., Science (2013), CC BY-NC-SA 4.0.

Recreate into a new output directory:
`.venv/bin/python scripts/export_all_landmarks.py --output results/exports/landmark_release_v2`
'''
    (out/'README.md').write_text(readme)
    out.rename(requested);out=requested
    archive=out.with_suffix('.zip')
    if archive.exists():raise FileExistsError('Archive exists; choose a fresh output')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(out.rglob('*')):
            if path.is_file():z.write(path,arcname=Path(out.name)/path.relative_to(out))
    print(json.dumps(json.loads((out/'checks.json').read_text()),indent=2))
    print('Bundle:',archive)

if __name__=='__main__':main()
