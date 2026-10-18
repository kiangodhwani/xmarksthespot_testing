"""Validation-only comparison with BigBrain's supplied brain mask."""
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
import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates
from landmarks.planes import pixel_to_world

URL='https://bigbrain-ftp.loris.ca/bigbrain-ftp/BigBrainRelease.2015/3D_Classified_Volumes/Histological_Space/nii/full_cls_200um_mask.nii.gz'

def metrics(pred,reference,valid):
    p=pred & valid; r=reference & valid
    tp=int((p&r).sum()); fp=int((p&~r).sum()); fn=int((~p&r).sum())
    return dict(dice=2*tp/(2*tp+fp+fn),iou=tp/(tp+fp+fn),
        reference_foreground_recall=tp/(tp+fn),reference_foreground_precision=tp/(tp+fp),
        agreement_foreground_pixels=tp,prediction_only_pixels=fp,reference_only_pixels=fn,
        caveat='Agreement with supplied brain mask, not tissue-only accuracy; CSF/envelope semantics differ')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/masks/coronal_y8_reference_check')
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Choose a new output directory')
    path=ROOT/'data/bigbrain/full_cls_200um_mask.nii.gz'
    img=nib.load(path)
    extraction=json.loads((ROOT/'results/pairs/coronal_y8_repeat/geometry.json').read_text())
    geometry=json.loads((ROOT/'results/transforms/coronal_y8_demo/truth/geometry.json').read_text())
    hist_affine=np.asarray(extraction['sources']['histology']['affine'])
    assert np.allclose(img.affine,hist_affine,atol=1e-7)
    assert list(img.shape)==extraction['sources']['histology']['shape']
    assert img.header.get_xyzt_units()[0]=='mm'
    a=img.affine
    assert np.allclose(a[:3,:3],np.diag(np.diag(a[:3,:3])),atol=1e-7)
    y=extraction['y_mm']; j=int(round((y-a[1,3])/a[1,1]))
    # Categorical label resampling must be nearest neighbour, never trilinear.
    section=np.asarray(img.dataobj[:,j,:])
    labels,counts=np.unique(section,return_counts=True)
    if len(labels)!=2 or labels[0]!=0: raise ValueError('Unexpected reference labels; inspect before binarising')
    section=(section>0).astype(np.uint8)
    np.testing.assert_array_equal(np.unique(section),[0,1])
    with np.load(ROOT/'results/transforms/coronal_y8_demo/inputs.npz') as d:
        image=d['histology']; support=d['histology_valid']
    with np.load(ROOT/'results/masks/coronal_y8_otsu/masks.npz') as d:
        pred=d['histology_raw']; cleaned=d['histology_cleaned']
    rows,cols=np.indices(image.shape,dtype=float)
    m=np.asarray(geometry['transformed_pixel_to_world_ras_mm'])
    world=pixel_to_world(np.stack([cols,rows],axis=-1),m)
    inverse=np.linalg.inv(a); ijk=world@inverse[:3,:3].T+inverse[:3,3]
    assert np.allclose(ijk[...,1],(y-a[1,3])/a[1,1],atol=1e-9)
    ref_valid=(ijk[...,0]>=0)&(ijk[...,0]<=img.shape[0]-1)&(ijk[...,2]>=0)&(ijk[...,2]<=img.shape[2]-1)
    reference=map_coordinates(section,np.stack([ijk[...,0],ijk[...,2]]),order=0,mode='constant',cval=0,prefilter=False).astype(bool)
    valid=support&ref_valid
    report=dict(source_url=URL,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        source_shape=list(img.shape),source_affine=a.tolist(),source_spacing_mm=list(map(float,img.header.get_zooms()[:3])),
        source_slope=float(img.dataobj.slope),source_intercept=float(img.dataobj.inter),
        scaled_labels_on_sampled_slice=labels.tolist(),label_counts_on_sampled_slice=counts.tolist(),
        foreground_rule='scaled mask value > 0; binary foreground/background only, no tissue-class distinction',
        documentation='BigBrain Updates.txt Feb 27 2018 describes tissue/CSF mask separating CSF from background; do not assume tissue-only',
        extraction='Direct nearest-neighbour source-label sampling at inverse-transformed target pixel centres; uses saved full canvas mapping',
        requested_plane_y_mm=y,nearest_categorical_source_y_mm=float(a[1,3]+j*a[1,1]),
        raw_otsu=metrics(pred,reference,valid),cleaned_otsu=metrics(cleaned,reference,valid),
        evaluation_only=True,attribution=geometry['attribution'])
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'reference_comparison.npz',reference=reference,comparison_valid=valid,
                        prediction_only=pred&~reference&valid,reference_only=reference&~pred&valid)
    (args.output/'report.json').write_text(json.dumps(report,indent=2))
    fig,axes=plt.subplots(1,4,figsize=(17,6),layout='constrained')
    axes[0].imshow(np.ma.array(image,mask=~support),cmap='gray',vmin=0,vmax=65535)
    axes[0].set_title('Transformed histology')
    axes[1].imshow(pred&valid,cmap='gray',vmin=0,vmax=1); axes[1].set_title('Independent raw Otsu')
    axes[2].imshow(reference&valid,cmap='gray',vmin=0,vmax=1); axes[2].set_title('Official brain-mask foreground\nMay include CSF/internal spaces')
    axes[3].imshow(np.ma.array(image,mask=~support),cmap='gray',vmin=0,vmax=65535)
    color=np.zeros((*image.shape,4))
    color[reference&~pred&valid]=[1,.1,.1,.75]
    color[pred&~reference&valid]=[0,.55,1,.9]
    axes[3].imshow(color); axes[3].set_title('Red: reference only\nBlue: Otsu only')
    # Same crop on all panels, only for display; saved arrays keep original indices.
    rr,cc=np.where(reference&valid)
    for ax in axes:
        ax.set_xlim(cc.min()-30,cc.max()+30); ax.set_ylim(rr.max()+30,rr.min()-30)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor('#dddddd')
    fig.suptitle('Validation reference only — disagreement includes both missed pale tissue and mask-definition differences')
    fig.get_layout_engine().set(rect=(0,.05,1,.95))
    fig.text(.5,.015,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=8)
    fig.savefig(args.output/'comparison.png',dpi=150); plt.close(fig)
    print(json.dumps(report,indent=2)); print('Saved:',args.output/'comparison.png')

if __name__=='__main__': main()
