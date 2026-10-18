"""Save one reproducible shared-plane pair, without synthetic transformations."""
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
import scipy
from landmarks.planes import make_grid,extract,pixel_to_world,world_to_pixel


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--y',type=float,default=8.)
    parser.add_argument('--output',type=Path,default=ROOT/'results'/'pairs'/'coronal_y8')
    args=parser.parse_args()
    sources={name:ROOT/'data'/'bigbrain'/file for name,file in [
        ('mri','mri_to_blk_aligned.nii.gz'),('histology','full8_200um_optbal.nii.gz')]}
    volumes={name:nib.load(path) for name,path in sources.items()}
    for img in volumes.values():
        if img.header.get_xyzt_units()[0]!='mm' or nib.aff2axcodes(img.affine)!=('R','A','S'):
            raise ValueError('Expected RAS volumes in millimetres')
    spacing=.2  # No spatial downsampling of either prescribed source volume.
    for img in volumes.values():
        if np.any(np.asarray(img.header.get_zooms()[:3]) < spacing-1e-7):
            raise ValueError('Finer source detected: implement antialiasing before downsampling')
    shape,matrix=make_grid([im.affine for im in volumes.values()],
                           [im.shape for im in volumes.values()],args.y,spacing)
    points=np.array([[0,0],[shape[1]-1,shape[0]-1],[53.25,71.75]])
    roundtrip=float(np.max(np.abs(world_to_pixel(pixel_to_world(points,matrix),matrix)-points)))
    if roundtrip>1e-9: raise ValueError('Coordinate round trip failed')
    # Refuse overwrites: each saved experiment is a fixed artifact.
    if args.output.exists(): raise FileExistsError(f'{args.output} exists; choose a new --output directory')
    arrays={}; provenance={}; ranges={}
    for name,img in volumes.items():
        print(f'Extracting {name} at y={args.y:g} mm...',flush=True)
        data=img.get_fdata(dtype=np.float32)
        if not np.isfinite(data).all(): raise ValueError('Nonfinite source intensities')
        ranges[name]=[float(data.min()),float(data.max())]
        image,valid=extract(data,img.affine,shape,matrix,background=65535 if name=='histology' else 0)
        arrays[name]=image; arrays[name+'_valid']=valid
        provenance[name]=dict(file=sources[name].name,sha256=sha256(sources[name]),
            shape=list(img.shape),affine=img.affine.tolist(),spacing_mm=list(map(float,img.header.get_zooms()[:3])),
            slope=float(img.dataobj.slope),intercept=float(img.dataobj.inter),
            scaled_intensity_range=ranges[name],valid_fraction=float(valid.mean()))
        del data
    arrays['joint_valid']=arrays['mri_valid'] & arrays['histology_valid']
    if not arrays['joint_valid'].any(): raise ValueError('No shared valid support')
    args.output.mkdir(parents=True)
    np.savez_compressed(args.output/'images.npz',**arrays)
    metadata=dict(schema_version=1,plane='coronal',y_mm=args.y,spacing_mm=spacing,
        shape_rows_cols=list(shape),pixel_point_order=['col','row'],array_index_order=['row','col'],
        pixel_centres_to_world_ras_mm=matrix.tolist(),
        inverse_formula='col=(x-matrix[0,2])/matrix[0,0]; row=(z-matrix[2,2])/matrix[2,1]; requires y=y_mm',
        field_of_view='Union of source voxel-centre bounds, expanded to 0.2 mm grid',
        interpolation='Trilinear; no downsampling; upsampled MRI contains no new anatomical information',
        validity='Inside source interpolation domain; not a tissue mask',sources=provenance,
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,nibabel=nib.__version__),
        checks=dict(max_pixel_roundtrip_error=roundtrip,shared_plane=True,shared_grid=True),
        reference_alignment='Preliminary review accepted by user; residual anatomical uncertainty remains; no global 2 mm guarantee',
        transform='None',attribution='BigBrain: Amunts et al., Science 2013; CC BY-NC-SA 4.0')
    (args.output/'geometry.json').write_text(json.dumps(metadata,indent=2))
    fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
    extent=[matrix[0,2]-.1,matrix[0,2]+(shape[1]-.5)*spacing,
            matrix[2,2]-(shape[0]-.5)*spacing,matrix[2,2]+.1]
    for ax,name in zip(axes,['mri','histology']):
        lo,hi=ranges[name]
        display=np.ma.array(arrays[name],mask=~arrays[name+'_valid'])
        ax.imshow(display,cmap='gray',origin='upper',extent=extent,vmin=lo,vmax=hi,interpolation='nearest')
        ax.set_title(f'{name.capitalize()}: y={args.y:g} mm')
        ax.set_xlabel('x (mm), increasing toward anatomical right'); ax.set_ylabel('z (mm)')
    fig.suptitle(f'Saved experiment pair: shared {spacing:g} mm grid; no synthetic transformation')
    fig.get_layout_engine().set(rect=(0,.04,1,.96))
    fig.text(.5,.01,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0. Masked margins are outside source support.',ha='center',fontsize=8)
    fig.savefig(args.output/'pair.png',dpi=150); plt.close(fig)
    # Check persisted data rather than only in-memory state.
    with np.load(args.output/'images.npz') as saved:
        for key in arrays:
            if not np.array_equal(saved[key],arrays[key]): raise ValueError('Saved array mismatch')
    print(f'Saved {shape[0]} rows × {shape[1]} columns to {args.output}')
    print(f'Coordinate round-trip error: {roundtrip:.3g} pixels; saved-array reload check passed.')

if __name__=='__main__': main()
