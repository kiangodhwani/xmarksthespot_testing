"""Histology hysteresis experiment; official mask read only AFTER prediction."""
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
from landmarks.masking import histology_hysteresis
from compare_reference_mask import metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/masks/coronal_y8_hysteresis')
    args=p.parse_args()
    if args.output.exists(): raise FileExistsError('Choose a new --output directory')
    path=ROOT/'results/transforms/coronal_y8_demo/inputs.npz'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(path) as d: image=d['histology']; valid=d['histology_valid']
    # No supplied labels, truth transforms, or MRI enter these computations.
    masks,stats=histology_hysteresis(image,valid)
    variants={}
    for weak in [.0025,.005,.01]:
        v,_=histology_hysteresis(image,valid,weak_fraction=weak,strong_fraction=.05)
        variants[str(weak)]=v['mask']
    with np.load(ROOT/'results/masks/coronal_y8_otsu/masks.npz') as d: otsu=d['histology_raw']
    # Evaluation only, after predictions and fixed sensitivity settings are generated.
    with np.load(ROOT/'results/masks/coronal_y8_reference_check/reference_comparison.npz') as d:
        reference=d['reference']; comparison_valid=d['comparison_valid']
    assert all(m.shape==image.shape and m.dtype==bool and not np.any(m&~valid) for m in [*masks.values(),*variants.values()])
    report=dict(status='Experimental; user review pending; not calibrated against tissue labels',
        parameters=stats,input_sha256=digest,reference_usage='Evaluation and display only; no parameter optimisation',
        otsu_envelope_agreement=metrics(otsu,reference,comparison_valid),
        hysteresis_envelope_agreement=metrics(masks['mask'],reference,comparison_valid),
        added_vs_otsu=int((masks['mask']&~otsu).sum()),lost_vs_otsu=int((otsu&~masks['mask']).sum()),
        weak_margin_sensitivity=[dict(weak_fraction=float(k),tissue_pixels=int(v.sum()),
            changed_vs_default=int((v!=masks['mask']).sum())) for k,v in variants.items()],
        caveat='Envelope includes CSF/spaces; overlap improvement is not tissue accuracy',
        coordinates='Unchanged; masks share transformed histology indices',
        attribution='BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.')
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'masks.npz',strong=masks['strong'],possible=masks['possible'],
                        hysteresis=masks['mask'],**{'weak_'+k:v for k,v in variants.items()})
    # Display crop is derived from image support only; no reference-based crop.
    rr,cc=np.where(valid); bounds=(cc.min()-5,cc.max()+5,rr.max()+5,rr.min()-5)
    def setup(ax,title):
        ax.set_title(title); ax.set_xlim(*bounds[:2]); ax.set_ylim(*bounds[2:])
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor('#dedede')
    def original(ax):
        ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=0,vmax=65535,interpolation='nearest')
    fig,axs=plt.subplots(2,3,figsize=(15,11),layout='constrained')
    original(axs[0,0]); setup(axs[0,0],'Original transformed histology')
    for ax,im,title in [(axs[0,1],otsu,'Otsu baseline'),(axs[0,2],masks['mask'],'Background-based hysteresis'),
                        (axs[1,0],reference&comparison_valid,'BigBrain brain envelope\nIncludes internal spaces')]:
        ax.imshow(im,cmap='gray',vmin=0,vmax=1,interpolation='nearest'); setup(ax,title)
    original(axs[1,1]); axs[1,1].contour(masks['mask'],levels=[.5],colors=['#e60070'],linewidths=.6)
    setup(axs[1,1],'Hysteresis boundary overlay (pink)')
    original(axs[1,2]); color=np.zeros((*image.shape,4))
    color[masks['mask']&~otsu]=[0,.7,.25,.6]; color[otsu&~masks['mask']]=[1,0,.2,.7]
    axs[1,2].imshow(color); setup(axs[1,2],'Green: added vs Otsu\nRed: removed vs Otsu')
    fig.suptitle('Histology masking experiment — fixed background thresholds; no hole filling or closing')
    fig.get_layout_engine().set(rect=(0,.04,1,.96)); fig.text(.5,.01,report['attribution'],ha='center',fontsize=9)
    fig.savefig(args.output/'comparison_panel.png',dpi=150); plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(15,6),layout='constrained')
    for ax,(k,mask) in zip(axs,variants.items()):
        original(ax); ax.contour(mask,levels=[.5],colors=['#e60070'],linewidths=.6)
        setup(ax,f'Weak margin: {float(k)*100:g}% of intensity range\n{int(mask.sum()):,} included pixels')
    fig.suptitle('Sensitivity check — smaller margin permits paler connected tissue')
    fig.get_layout_engine().set(rect=(0,.05,1,.95)); fig.text(.5,.01,report['attribution'],ha='center',fontsize=9)
    fig.savefig(args.output/'sensitivity.png',dpi=150); plt.close(fig)
    with np.load(args.output/'masks.npz') as saved:
        assert np.array_equal(saved['hysteresis'],masks['mask'])
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
    (args.output/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2)); print('Saved:',args.output/'comparison_panel.png')

if __name__=='__main__': main()
