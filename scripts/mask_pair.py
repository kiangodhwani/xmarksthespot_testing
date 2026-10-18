"""Try Otsu independently on MRI and transformed histology; save review plots."""
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
import skimage
from landmarks.masking import otsu_mask


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'results/transforms/coronal_y8_demo/inputs.npz')
    parser.add_argument('--output',type=Path,default=ROOT/'results/masks/coronal_y8_otsu')
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError('Output exists; use a new --output directory')
    source_hash=hashlib.sha256(args.input.read_bytes()).hexdigest()
    with np.load(args.input) as loaded: data={k:loaded[k] for k in loaded.files}
    masks={}; reports={}
    for name in ['mri','histology']:
        m,report=otsu_mask(data[name],data[name+'_valid'],tissue_is_bright=name=='mri')
        for key,value in m.items():
            assert value.shape==data[name].shape and value.dtype==bool
            assert not np.any(value & ~data[name+'_valid'])
            masks[name+'_'+key]=value
        reports[name]=report
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'masks.npz',**masks)
    fig,axes=plt.subplots(2,4,figsize=(16,9),layout='constrained')
    sensitivity,saxes=plt.subplots(2,3,figsize=(13,9),layout='constrained')
    histograms,haxes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for row,name in enumerate(['mri','histology']):
        image=data[name]; valid=data[name+'_valid']; values=image[valid]
        raw=masks[name+'_raw']; clean=masks[name+'_cleaned']
        lo,hi=float(values.min()),float(values.max())
        display=np.ma.array(image,mask=~valid)
        for col in [0,3]:
            axes[row,col].imshow(display,cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
        axes[row,0].set_title(f'{name.capitalize()} original')
        axes[row,1].imshow(raw,cmap='gray',vmin=0,vmax=1,interpolation='nearest'); axes[row,1].set_title('Raw Otsu tissue mask')
        axes[row,2].imshow(clean,cmap='gray',vmin=0,vmax=1,interpolation='nearest'); axes[row,2].set_title('After removing ≤25-pixel specks')
        axes[row,3].contour(raw,levels=[.5],colors=['#e60060'],linewidths=.6)
        axes[row,3].set_title('Raw mask boundary (pink)')
        for col,key in enumerate(['lower_threshold','raw','upper_threshold']):
            saxes[row,col].imshow(display,cmap='gray',vmin=lo,vmax=hi,interpolation='nearest')
            saxes[row,col].contour(masks[name+'_'+key],levels=[.5],colors=['#e60060'],linewidths=.6)
            saxes[row,col].set_title(f'{name}: '+['threshold − 2% range','Otsu','threshold + 2% range'][col])
        haxes[row].hist(values,bins=256,color='#555555')
        haxes[row].axvline(reports[name]['threshold'],color='#e60060',label='Otsu')
        haxes[row].set_yscale('log'); haxes[row].set_title(name.capitalize()+' valid-pixel histogram')
        haxes[row].set_xlabel('Scaled intensity'); haxes[row].set_ylabel('Pixel count (log)'); haxes[row].legend()
    for ax in axes.flat:
        ax.set_facecolor('#dddddd'); ax.set_xticks([]); ax.set_yticks([])
    for ax in saxes.flat:
        ax.set_facecolor('#dddddd'); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle('Otsu baseline — inspect pale tissue omissions; masks are not yet accepted')
    sensitivity.suptitle('Threshold sensitivity: pink contours; no anatomical truth used')
    for figure,file in [(fig,'review.png'),(sensitivity,'sensitivity.png'),(histograms,'histograms.png')]:
        figure.get_layout_engine().set(rect=(0,.04,1,.96))
        figure.text(.5,.01,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0.',ha='center',fontsize=8)
        figure.savefig(args.output/file,dpi=140); plt.close(figure)
    report=dict(status='Baseline for visual review, not validated tissue segmentation',
        input=str(args.input.resolve()),input_sha256=source_hash,
        input_access='Only images and support masks; no truth geometry or supplied tissue labels',
        versions=dict(numpy=np.__version__,skimage=skimage.__version__),
        cleanup='Remove 8-connected components of at most 25 working-image pixels; no physical-scale claim; retain raw masks',
        coordinates='Unchanged: masks match corresponding image array shape and index positions',results=reports)
    if hashlib.sha256(args.input.read_bytes()).hexdigest()!=source_hash: raise ValueError('Input modified')
    with np.load(args.output/'masks.npz') as saved:
        assert all(np.array_equal(saved[k],v) for k,v in masks.items())
    (args.output/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2)); print('Saved review:',args.output/'review.png')

if __name__=='__main__': main()
