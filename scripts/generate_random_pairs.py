"""Day 1: seeded random similarity test pairs from one accepted source plane."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import argparse
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from transform_pair import generate_pair
from landmarks.transforms import sample_parameters


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--count',type=int,default=5)
    p.add_argument('--input',type=Path,default=ROOT/'results/pairs/coronal_y8_repeat')
    p.add_argument('--output',type=Path,default=ROOT/'results/random_pairs/y8_seed42')
    args=p.parse_args()
    params=sample_parameters(args.seed,args.count)
    if args.output.exists(): raise FileExistsError('Batch exists; choose a new --output directory')
    if not (args.input/'images.npz').is_file(): raise FileNotFoundError('Extract the shared plane first')
    print('Randomisation plan (no fitting): seed',args.seed,flush=True)
    for i,values in enumerate(params,1): print(f'pair_{i:03d}: {values}',flush=True)
    args.output.mkdir(parents=True); truth=args.output/'truth';truth.mkdir()
    source_geometry=json.loads((args.input/'geometry.json').read_text())
    manifest=dict(seed=args.seed,count=args.count,generator='numpy.default_rng / PCG64',
        numpy_version=np.__version__,sampling=dict(angle_degrees_uniform=[-30,30],
        scale_uniform=[.85,1.2],translation_pixels_inclusive=[-50,50]),
        source_plane=dict(y_mm=source_geometry['y_mm'],note='One fixed plane; pairs are augmentations, not independent specimens'),
        status='in_progress',pairs=[])
    (truth/'manifest.json').write_text(json.dumps(manifest,indent=2))
    for i,values in enumerate(params,1):
        directory=args.output/f'pair_{i:03d}'
        generate_pair(args.input,directory,**values,randomisation=dict(seed=args.seed,index=i,
            angle_range=[-30,30],scale_range=[.85,1.2],shift_inclusive_range=[-50,50]))
        with np.load(directory/'inputs.npz') as d:
            assert d['histology'].shape==d['histology_valid'].shape
        report=json.loads((directory/'truth/geometry.json').read_text())
        manifest['pairs'].append(dict(id=directory.name,parameters=values,
            histology_shape_rows_cols=report['histology_shape_rows_cols'],checks=report['checks'],
            antialias_sigma_source_px=report['antialias_sigma_source_px']))
        (truth/'manifest.json').write_text(json.dumps(manifest,indent=2))
    manifest['status']='complete';(truth/'manifest.json').write_text(json.dumps(manifest,indent=2))
    with np.load(args.input/'images.npz') as d:
        original=d['histology'];original_valid=d['histology_valid'];mri=d['mri'];mri_valid=d['mri_valid']
    cols=3;rows=int(np.ceil((args.count+2)/cols))
    fig,axs=plt.subplots(rows,cols,figsize=(15,rows*4.6),layout='constrained')
    axs=np.asarray(axs).ravel()
    axs[0].imshow(np.ma.array(mri,mask=~mri_valid),cmap='gray',vmin=0,vmax=250);axs[0].set_title('Fixed MRI — identical for every pair')
    axs[1].imshow(np.ma.array(original,mask=~original_valid),cmap='gray',vmin=0,vmax=65535);axs[1].set_title('Original histology — starting image')
    for ax,item in zip(axs[2:],manifest['pairs']):
        with np.load(args.output/item['id']/'inputs.npz') as d:
            image=d['histology'];valid=d['histology_valid']
        ax.imshow(np.ma.array(image,mask=~valid),cmap='gray',vmin=0,vmax=65535)
        q=item['parameters'];dx,dy=q['shift']
        ax.set_title(f'{item["id"]}: {q["angle"]:+.1f}°; scale {q["scale"]:.3f}×\nshift {dx:+d} col / {dy:+d} row px')
    for i,ax in enumerate(axs):
        ax.set_facecolor('#ddd');ax.set_xticks([]);ax.set_yticks([])
        if i>=args.count+2: ax.set_visible(False)
    fig.suptitle(f'Day 1 random test pairs — seed {args.seed}; + angle clockwise, +col right, +row down\nPanels fit their own canvases: compare exact shifts using each pair’s full preview, not panel positions')
    fig.get_layout_engine().set(rect=(0,.04,1,.96))
    fig.text(.5,.01,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0. Grey = invalid support; no flips/crops/tilt.',ha='center',fontsize=9)
    fig.savefig(args.output/'overview.png',dpi=140);plt.close(fig)
    print('Batch complete. Overview:',args.output/'overview.png')
    print('Exact parameters and check results:',truth/'manifest.json')

if __name__=='__main__': main()
