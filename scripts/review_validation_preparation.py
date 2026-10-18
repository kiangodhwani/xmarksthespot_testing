"""Inspect new-slice preparation before user review; no detection or matching."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import hashlib
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from landmarks.masking import otsu_mask,histology_hysteresis
from landmarks.planes import pixel_to_world,world_to_pixel


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,default=ROOT/'results/validation/two_slices_v1')
    args=parser.parse_args();root=args.study
    plan=json.loads((root/'study_plan.json').read_text())
    for path,expected in plan['frozen_source_sha256'].items():
        if digest(ROOT/path)!=expected:raise ValueError(f'Frozen pipeline changed: {path}')
    fig,axs=plt.subplots(2,3,figsize=(16,11))
    fig.subplots_adjust(left=.055,right=.99,top=.89,bottom=.13,hspace=.24,wspace=.14)
    report=dict(status='numerically_checked_user_visual_review_pending',matching_run=False,
        frozen_pipeline_hashes_verified=True,study_plan_sha256=digest(root/'study_plan.json'),slices=[])
    expected_sources=dict(mri='098edc7b3a787890c39f12fb8fe31e6eaca703b341fd60ef6c653e0e269cf691',
        histology='69f19eb044d23a09c0c316e8ce985fcd6822914d19ac0ec36e72aa8e03026bef')
    for row,(name,y) in enumerate([('y_m20',-20.),('y35',35.)]):
        path=root/name
        with np.load(path/'pair/images.npz') as d:arrays={k:d[k] for k in d.files}
        meta=json.loads((path/'pair/geometry.json').read_text())
        assert meta['y_mm']==y and meta['spacing_mm']==.2
        for modality,sha in expected_sources.items():assert meta['sources'][modality]['sha256']==sha
        shape=arrays['mri'].shape
        assert shape==arrays['histology'].shape==(631,700)
        rr,cc=np.indices(shape);points=np.stack([cc,rr],axis=-1)
        matrix=np.asarray(meta['pixel_centres_to_world_ras_mm'])
        world=pixel_to_world(points,matrix)
        assert np.all(world[...,1]==y)
        error=float(np.max(np.abs(world_to_pixel(world,matrix)-points)))
        assert error<1e-9
        with np.load(path/'fixed/inputs.npz') as d:fixed={k:d[k] for k in d.files}
        checks=[]
        variants=[path/'fixed']+sorted((path/'random_seed42').glob('pair_*'))
        assert len(variants)==6
        for variant in variants:
            geometry=json.loads((variant/'truth/geometry.json').read_text())
            assert geometry['checks']['no_footprint_clipping']
            assert geometry['checks']['full_grid_roundtrip_max_px']<1e-9
            assert geometry['checks']['full_grid_world_roundtrip_max_mm']<1e-9
            with np.load(variant/'inputs.npz') as d:
                np.testing.assert_array_equal(d['mri'],arrays['mri'])
                np.testing.assert_array_equal(d['mri_valid'],arrays['mri_valid'])
                assert d['histology'].shape==d['histology_valid'].shape
                assert np.isfinite(d['histology']).all()
            checks.append(dict(case=variant.name,**geometry['checks']))
        manifest=json.loads((path/'random_seed42/truth/manifest.json').read_text())
        assert manifest['status']=='complete' and manifest['seed']==42
        assert [p['parameters'] for p in manifest['pairs']]==plan['random_parameters']
        mri_mask,mri_report=otsu_mask(arrays['mri'],arrays['mri_valid'],True)
        hist_mask,hist_report=histology_hysteresis(arrays['histology'],arrays['histology_valid'],weak_fraction=.01)
        trans_mask,trans_report=histology_hysteresis(fixed['histology'],fixed['histology_valid'],weak_fraction=.01)
        displayed=[(arrays['mri'],arrays['mri_valid'],mri_mask['cleaned'],'Original MRI\ncyan: Otsu search support'),
            (arrays['histology'],arrays['histology_valid'],hist_mask['mask'],'Original histology\ncyan: 1% hysteresis search support'),
            (fixed['histology'],fixed['histology_valid'],trans_mask['mask'],'Fixed transformed histology\ncyan: independently computed 1% support')]
        for ax,(image,valid,mask,title) in zip(axs[row],displayed):
            assert mask.shape==valid.shape and not np.any(mask&~valid)
            lo,hi=(0,250) if 'MRI' in title else (0,65535)
            ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=lo,vmax=hi)
            ax.contour(mask.astype(float),levels=[.5],colors=['#00a6bf'],linewidths=.55)
            ax.set_title(title,fontsize=10);ax.set_facecolor('#ddd');ax.set_xticks([]);ax.set_yticks([])
        axs[row,0].set_ylabel(f'Coronal y={y:g} mm',fontsize=12)
        np.savez_compressed(path/'preparation_support_masks.npz',mri=mri_mask['cleaned'],
            histology_original=hist_mask['mask'],histology_fixed=trans_mask['mask'])
        report['slices'].append(dict(id=name,y_mm=y,shape_rows_cols=list(shape),spacing_mm=.2,
            full_grid_pixel_world_roundtrip_max_px=error,source_hashes_verified=True,
            arrays_sha256=digest(path/'pair/images.npz'),six_transform_checks=checks,
            mask_reports=dict(mri=mri_report,histology_original=hist_report,histology_fixed=trans_report),
            caution='Mask outlines are provisional search support, not tissue-only segmentation or matching predictions'))
    fig.suptitle('Two additional slices — PREPARATION REVIEW ONLY; no matching has been run\nFirst two columns share the original grid; third column uses the fixed +20° /1.15×/(+40,+25) px transform',fontsize=13)
    fig.text(.5,.028,'Cyan outlines: independently computed provisional search support. Images retain original intensities and internal structures.\nGrey = invalid sampling support. Panels fit their own canvases; do not compare apparent scale with the third column.\nReference alignment remains approximate. These are additional slices of the SAME specimen, not independent subjects.\nBigBrain — Amunts et al., Science (2013), CC BY-NC-SA 4.0.',ha='center',fontsize=10)
    fig.savefig(root/'preparation_review.png',dpi=140);plt.close(fig)
    (root/'preparation_checks.json').write_text(json.dumps(report,indent=2))
    print('Verified frozen pipeline; both full grids; all 12 transformations; source hashes, unchanged MRI, masks and seed parameters.')
    print('Saved',root/'preparation_review.png')
    print('Matching NOT run. User visual review pending.')

if __name__=='__main__':main()
