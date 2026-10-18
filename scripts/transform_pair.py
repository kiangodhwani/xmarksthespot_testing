"""Apply the agreed synthetic histology transform, retaining evaluation geometry."""
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
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from landmarks.transforms import expanded_transform,apply_points,warp,footprint,antialias_sigma
from landmarks.planes import pixel_to_world


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate_pair(input_dir,output_dir,angle=20.,scale=1.15,shift=(40.,25.),randomisation=None):
    if output_dir.exists(): raise FileExistsError('Output exists; choose another --output')
    source=input_dir/'images.npz'; geometry=input_dir/'geometry.json'
    hashes={p.name:digest(p) for p in [source,geometry]}
    meta=json.loads(geometry.read_text())
    with np.load(source) as saved: data={k:saved[k] for k in saved.files}
    hist=data['histology']; shape=hist.shape
    if data['mri'].shape!=shape or not np.isfinite(hist).all(): raise ValueError('Invalid source pair')
    t,outshape,offset,requested=expanded_transform(shape,angle,scale,shift)
    inverse=np.linalg.inv(t)
    transformed,valid=warp(hist,data['histology_valid'],t,outshape)
    # Test round trips on EVERY original pixel centre, plus image-footprint corners.
    rows,cols=np.indices(shape,dtype=float)
    original=np.stack([cols,rows],axis=-1)
    moved=apply_points(original,t)
    back=apply_points(moved,inverse)
    error=float(np.max(np.abs(back-original)))
    if error>1e-9: raise ValueError('Point roundtrip check failed')
    world_matrix=np.asarray(meta['pixel_centres_to_world_ras_mm'])
    moved_to_world=world_matrix@inverse
    original_world=pixel_to_world(original,world_matrix)
    moved_world=pixel_to_world(moved,moved_to_world)
    world_error=float(np.max(np.linalg.norm(moved_world-original_world,axis=-1)))
    if world_error>1e-9: raise ValueError('World mapping check failed')
    output_dir.mkdir(parents=True); truth=output_dir/'truth'; truth.mkdir()
    # Future matcher receives images and support only. Evaluation metadata is separate.
    np.savez_compressed(output_dir/'inputs.npz',mri=data['mri'],mri_valid=data['mri_valid'],
                        histology=transformed,histology_valid=valid)
    np.savez_compressed(truth/'original_histology.npz',image=hist,valid=data['histology_valid'])
    report=dict(schema_version=1,parameters=dict(clockwise_degrees=float(angle),scale=float(scale),
        translation_col_row_px=list(map(float,shift)),translation_in_original_grid_mm=(np.asarray(shift)*meta['spacing_mm']).tolist(),
        translation_mm_direction='Positive column shift is anatomical right; positive row shift is inferior; before canvas offset',
        centre_col_row_px=[(shape[1]-1)/2,(shape[0]-1)/2],padding_px=10),
        operations='Scale and clockwise rotation about original centre, translation, canvas offset; one bilinear image warp',
        point_convention='[col,row,1]; row increases downward',
        requested_original_to_transformed_without_canvas=requested.tolist(),canvas_offset_col_row_px=offset.tolist(),
        original_to_transformed=t.tolist(),transformed_to_original=inverse.tolist(),
        original_pixel_to_world_ras_mm=world_matrix.tolist(),transformed_pixel_to_world_ras_mm=moved_to_world.tolist(),
        mri_shape_rows_cols=list(shape),histology_shape_rows_cols=list(outshape),
        transformed_source_footprint=apply_points(footprint(shape),t).tolist(),
        canvas_policy='Original and zero-translation rotated/scaled footprints plus symmetric translation margins; translation not cropped away',
        interpolation='Gaussian antialias prefilter for reduction, then bilinear; invalid filter/interpolation support excluded; background 65535',
        antialias_sigma_source_px=antialias_sigma(t),randomisation=randomisation,
        scale_effect=f'Underlying sampling {meta["spacing_mm"]/scale:.8g} mm along rotated axes; use saved world mapping',
        ground_truth_scope='Synthetic coordinate truth only; original MRI/histology anatomical discrepancies remain',
        inputs=dict(directory=str(input_dir.resolve()),sha256=hashes),
        versions=dict(numpy=np.__version__,scipy=scipy.__version__),
        checks=dict(full_grid_roundtrip_max_px=error,full_grid_world_roundtrip_max_mm=world_error,no_footprint_clipping=True),
        attribution=meta['attribution'])
    (truth/'geometry.json').write_text(json.dumps(report,indent=2))
    # These are synthetic coordinate probes, not anatomical or detector landmarks.
    probes=np.array([[200.,180.],[480.,200.],[350.,420.]])
    destinations=apply_points(probes,t)
    with (truth/'coordinate_probes.csv').open('w',newline='') as f:
        writer=csv.writer(f); writer.writerow(['id','original_col','original_row','transformed_col','transformed_row','world_x_mm','world_y_mm','world_z_mm'])
        for i,(p,q,w) in enumerate(zip(probes,destinations,pixel_to_world(probes,world_matrix)),1):
            writer.writerow([i,*p,*q,*w])
    # Common display canvas preserves the visible translation relative to original.
    fig,axs=plt.subplots(1,3,figsize=(16,7),layout='constrained')
    for ax in axs:
        ax.set_xlim(-.5,outshape[1]-.5); ax.set_ylim(outshape[0]-.5,-.5)
        ax.set_aspect('equal'); ax.set_xlabel('Canvas column (px)'); ax.set_ylabel('Canvas row (px)')
        ax.set_facecolor('#e6e6e6')
    ext=[offset[0]-.5,offset[0]+shape[1]-.5,offset[1]+shape[0]-.5,offset[1]-.5]
    axs[0].imshow(np.ma.array(data['mri'],mask=~data['mri_valid']),cmap='gray',origin='upper',extent=ext,vmin=0,vmax=meta['sources']['mri']['scaled_intensity_range'][1])
    axs[1].imshow(np.ma.array(hist,mask=~data['histology_valid']),cmap='gray',origin='upper',extent=ext,vmin=0,vmax=65535)
    axs[2].imshow(np.ma.array(transformed,mask=~valid),cmap='gray',origin='upper',vmin=0,vmax=65535)
    for ax,title in zip(axs,['MRI — unchanged','Histology — original','Histology — transformed']): ax.set_title(title)
    for i,(p,q) in enumerate(zip(probes+offset,destinations),1):
        for ax,point in [(axs[1],p),(axs[2],q)]:
            ax.plot(*point,'o',color=['#e58a00','#0096d6','#df458c'][i-1],markersize=6)
            ax.annotate(str(i),point,xytext=(5,5),textcoords='offset points',color='black',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
    fig.suptitle(f'{angle:+.2f}° (positive = clockwise) · {scale:.3f}× scale · shift (col,row) = ({shift[0]:+g}, {shift[1]:+g}) px\nDots are synthetic coordinate probes, not anatomical matches')
    fig.get_layout_engine().set(rect=(0,.06,1,.96))
    fig.text(.5,.02,'BigBrain — Amunts et al., Science (2013). CC BY-NC-SA 4.0. Grey = outside valid support.',ha='center',fontsize=9)
    fig.savefig(output_dir/'preview.png',dpi=140); plt.close(fig)
    with np.load(output_dir/'inputs.npz') as saved:
        if not np.array_equal(saved['mri'],data['mri']): raise ValueError('MRI changed')
        if not np.array_equal(saved['histology'],transformed): raise ValueError('Saved histology differs')
    if {p.name:digest(p) for p in [source,geometry]}!=hashes: raise ValueError('Source modified')
    print(f'Saved: {output_dir}')
    print(f'Histology canvas: {outshape}; roundtrip: {error:.3g} px / {world_error:.3g} mm')
    print('PASS: full-grid mappings, footprint containment, MRI preservation, source hashes and saved-array reload.')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'results/pairs/coronal_y8_repeat')
    parser.add_argument('--output',type=Path,default=ROOT/'results/transforms/coronal_y8_demo')
    parser.add_argument('--angle',type=float,default=20.)
    parser.add_argument('--scale',type=float,default=1.15)
    parser.add_argument('--shift',type=float,nargs=2,default=(40.,25.),metavar=('COL','ROW'))
    args=parser.parse_args()
    generate_pair(args.input,args.output,args.angle,args.scale,args.shift)

if __name__=='__main__': main()
