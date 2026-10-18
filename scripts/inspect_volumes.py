"""Inspect source headers and native coronal slices; no resampling or registration."""
from pathlib import Path
import json
import hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'inspection'
OUT.mkdir(parents=True, exist_ok=True)
FILES = {'MRI': 'mri_to_blk_aligned.nii.gz',
         'Histology': 'full8_200um_optbal.nii.gz'}
planes = [-20., 8., 35.]
fig, axes = plt.subplots(3, 2, figsize=(12, 14), constrained_layout=True)
report = {}
for col, (name, filename) in enumerate(FILES.items()):
    path = ROOT / 'data' / 'bigbrain' / filename
    img = nib.load(path)
    a = img.affine
    assert np.isfinite(a).all() and abs(np.linalg.det(a[:3, :3])) > 1e-10
    assert img.header.get_xyzt_units()[0] == 'mm'
    assert nib.aff2axcodes(a) == ('R', 'A', 'S')
    assert np.allclose(a[:3, :3], np.diag(np.diag(a[:3, :3])), atol=1e-5)
    data = img.get_fdata(dtype=np.float32)
    assert np.isfinite(data).all()
    report[name] = dict(file=filename, shape=list(img.shape),
        spacing_mm=list(map(float, img.header.get_zooms()[:3])),
        orientation=list(nib.aff2axcodes(a)), affine=a.tolist(),
        units=img.header.get_xyzt_units()[0],
        scaled_range=[float(data.min()), float(data.max())],
        slope=float(img.dataobj.slope), intercept=float(img.dataobj.inter),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(), slices=[])
    for row, y in enumerate(planes):
        j = int(round((y-a[1, 3])/a[1, 1]))
        assert 0 <= j < img.shape[1]
        actual_y = float(a[1, 3]+j*a[1, 1])
        # Transpose x,z to rows/cols, then flip so superior is at top.
        section = data[:, j, :].T[::-1, :]
        np.save(OUT / f'{name.lower()}_y_{y:g}.npy', section)
        extent = [a[0,3]-.5*a[0,0], a[0,3]+(img.shape[0]-.5)*a[0,0],
                  a[2,3]-.5*a[2,2], a[2,3]+(img.shape[2]-.5)*a[2,2]]
        ax = axes[row, col]
        ax.imshow(section, cmap='gray', origin='upper', extent=extent,
                  vmin=float(data.min()), vmax=float(data.max()), interpolation='nearest')
        ax.set_xlim(-75, 75)
        ax.set_ylim(-65, 85)
        ax.set_title(f'{name}: requested y={y:g}; native y={actual_y:.3f} mm\nindex j={j}; spacing {a[0,0]:.3f} mm')
        ax.set_xlabel('x (mm): left → right anatomy')
        ax.set_ylabel('z (mm): inferior → superior')
        report[name]['slices'].append(dict(requested_y_mm=y, actual_y_mm=actual_y, index=j))
    del data
fig.suptitle('Native coronal slices — inspection only; no geometric transformations', fontsize=15)
fig.text(.5, .002, 'BigBrain: Amunts et al., Science (2013). CC BY-NC-SA 4.0. Native planes differ slightly; not evaluation pairs.', ha='center', fontsize=8)
fig.get_layout_engine().set(rect=(0, 0.025, 1, 0.975))
fig.savefig(OUT / 'native_coronal_preview.png', dpi=140)
plt.close(fig)
(OUT / 'headers.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
