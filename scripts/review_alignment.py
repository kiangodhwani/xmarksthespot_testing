"""Alignment review only: exact-plane previews and optional manual annotation.
Run with --annotate in PyCharm to select independently defined anatomical pairs.
"""
from pathlib import Path
import argparse
import csv
from datetime import datetime
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'xmarks-mpl'))
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'alignment_review'
SPACING = .4
X = np.arange(-70, 69.00001, SPACING)
Z = np.arange(67, -58.00001, -SPACING)

def sample(data, affine, world):
    ijk = nib.affines.apply_affine(np.linalg.inv(affine), world)
    valid = np.all((ijk >= 0) & (ijk <= np.array(data.shape)-1), axis=-1)
    values = map_coordinates(data, ijk.reshape(-1, 3).T, order=1,
                             mode='constant', cval=0, prefilter=False)
    return values.reshape(world.shape[:-1]), valid

def check_sampling():
    # A linear world-space signal is reproduced exactly by trilinear sampling.
    a = np.array([[.7,0,0,-3],[0,.8,0,2],[0,0,.9,-1],[0,0,0,1.]])
    indices = np.moveaxis(np.indices((12,13,14), dtype=float),0,-1)
    w = nib.affines.apply_affine(a, indices)
    data = 2*w[...,0] + 3*w[...,1] - w[...,2] + 7
    points = nib.affines.apply_affine(a, np.array([[2.2,3.4,4.6],[7.1,8.2,9.3]]))
    got, valid = sample(data,a,points)
    expected = 2*points[:,0]+3*points[:,1]-points[:,2]+7
    assert valid.all() and np.allclose(got,expected,atol=1e-10)
    assert np.allclose(nib.affines.apply_affine(np.linalg.inv(a),points),
                       [[2.2,3.4,4.6],[7.1,8.2,9.3]],atol=1e-10)
    return float(np.max(np.abs(got-expected)))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotate', action='store_true')
    parser.add_argument('--backend', choices=['TkAgg', 'MacOSX'], help='Interactive window backend')
    args = parser.parse_args()
    import matplotlib
    if not args.annotate:
        matplotlib.use('Agg')
    elif args.backend:
        matplotlib.use(args.backend, force=True)
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    print(f'Matplotlib {matplotlib.__version__}; backend: {matplotlib.get_backend()}', flush=True)
    error = check_sampling()
    xx, zz = np.meshgrid(X,Z)
    ys = [-20.,8.,35.]
    images = {}
    for name,file in [('MRI','mri_to_blk_aligned.nii.gz'),('Histology','full8_200um_optbal.nii.gz')]:
        img = nib.load(ROOT/'data'/'bigbrain'/file)
        assert img.header.get_xyzt_units()[0]=='mm'
        data = img.get_fdata(dtype=np.float32)
        for y in ys:
            world = np.stack([xx,np.full_like(xx,y),zz],axis=-1)
            section, valid = sample(data,img.affine,world)
            np.savez_compressed(OUT/f'{name.lower()}_y_{y:g}.npz',
                image=section,valid=valid,x_mm=X,z_mm=Z,y_mm=y)
            low, high = float(data.min()), float(data.max())
            display = np.clip((section-low)/(high-low),0,1)
            if name=='Histology':
                display=1-display  # Display-only inversion; arrays untouched.
            display[~valid]=0
            images[name,y]=display
        del data
    extent=[X[0]-SPACING/2,X[-1]+SPACING/2,Z[-1]-SPACING/2,Z[0]+SPACING/2]
    fig,axs=plt.subplots(3,3,figsize=(14,13),layout='constrained')
    for row,y in enumerate(ys):
        m,h=images['MRI',y],images['Histology',y]
        rgb=np.stack([m,h,h],axis=-1)
        for col,(im,title) in enumerate([(m,'MRI'),(h,'Histology (display inverted)'),(rgb,'Overlay: MRI red / histology cyan')]):
            ax=axs[row,col]
            ax.imshow(im,cmap='gray',vmin=0,vmax=1,origin='upper',extent=extent)
            ax.set_title(f'{title}\ny = {y:g} mm')
            ax.set_xlabel('x (mm)'); ax.set_ylabel('z (mm)')
    fig.suptitle('Alignment review — same physical plane; colour mismatch is not a distance measurement')
    fig.savefig(OUT/'review.png',dpi=130)
    plt.close(fig)
    report=dict(stage='alignment review, not production pair generation',
        spacing_mm=SPACING,planes_y_mm=ys,interpolation='linear',
        analytic_sampling_max_error=error,
        anatomical_error_mm=None,
        status='Awaiting independent anatomical annotations; no 2 mm accuracy claim',
        attribution='BigBrain, Amunts et al. Science 2013; CC BY-NC-SA 4.0')
    if args.annotate:
        records=[]
        session = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        csv_path = OUT / f'manual_landmarks_{session}.csv'
        with csv_path.open('x', newline='') as f:
            csv.writer(f).writerow(['name','y_mm','mri_x_mm','mri_z_mm','hist_x_mm','hist_z_mm','in_plane_error_mm'])
        print(f'Each completed pair will be saved immediately to: {csv_path}', flush=True)
        print('Define each landmark anatomically BEFORE clicking. Select MRI then histology. Close/Enter to finish each plane.')
        for y in ys:
            while True:
                label=input(f'y={y:g}: anatomical landmark name (blank ends plane): ').strip()
                if not label: break
                points=[]
                for name in ['MRI','Histology']:
                    fig,ax=plt.subplots(figsize=(8,8))
                    ax.imshow(images[name,y],cmap='gray',origin='upper',extent=extent,vmin=0,vmax=1)
                    ax.set_title(f'{name}: {label}; click once to select')
                    try:
                        chosen=fig.ginput(1,timeout=0)
                    finally:
                        plt.close(fig)
                    if not chosen: break
                    points.append(chosen[0])
                if len(points)==2:
                    p,q=np.array(points)
                    record = [label,y,*p,*q,float(np.linalg.norm(p-q))]
                    with csv_path.open('a', newline='') as f:
                        csv.writer(f).writerow(record)
                        f.flush()
                        os.fsync(f.fileno())
                    records.append(record)
                    print(f'Saved: {label}, separation {record[-1]:.3f} mm', flush=True)
        report['landmark_csv'] = str(csv_path)
        if records:
            e=np.array([r[-1] for r in records])
            report['manual_in_plane_errors']=dict(count=len(e),mean=float(e.mean()),median=float(np.median(e)),maximum=float(e.max()),fraction_within_2mm=float((e<=2).mean()),fraction_within_5mm=float((e<=5).mean()))
            report['status']='Manual 2D sample only; repeat annotations and 3D review required'
    report_path = OUT / (f'manual_report_{session}.json' if args.annotate else 'report.json')
    report_path.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    print('Saved preview:',OUT/'review.png')

if __name__=='__main__': main()
