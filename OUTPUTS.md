# Registration-ready landmark outputs

The complete release is `results/exports/landmark_release_v1/`, with a downloadable
archive at `results/exports/landmark_release_v1.zip` (about 0.75 MB).

Start with:

- `final_pipeline/`: one CSV for each of the 19 final-pipeline MRI/histology pairs.
- `final_pipeline_summary.csv`: success metrics for those 19 cases.
- `all_tests_summary.csv`: all 100 saved prediction artifacts and their metrics.
- `historical/` and `historical_summary.csv`: 67 earlier baseline/ablation variants.
- `controls/` and `controls_summary.csv`: 14 histology-only controls.
- `manifest.json`: full provenance, copied/generated CSV status, exact duplicate
  predictions, exclusions and an explicitly labelled legacy-control replay alias.

MRI/histology CSVs use exactly:

```csv
id,name,mri_x_mm,mri_y_mm,mri_z_mm,hist_col_px,hist_row_px,score
```

World coordinates are RAS mm, computed from the saved source pixel-centre
mapping. Histology coordinates are column/row pixel centres in the transformed
section, including its expanded canvas offset. Names such as candidate_12 are
algorithm IDs, not anatomical labels. Score is max-bank cosine, not probability.
Rows are final algorithm predictions; reference error never filters export rows.

Histology-only controls have correctly named `histology_original_x_mm/y_mm/z_mm`
source columns in a separate folder; these are not mislabelled as MRI.
All failed attempts are retained as header-only CSVs, with null precision/error
fields. The 47 exact repeated prediction records are flagged in summary/provenance,
not represented as independent anatomical evidence. The old polar-baseline
control report has no saved prediction artifact; an exact-report backend replay
is recorded as an alias with its real coordinate provenance.

Every case includes separate `evaluation/metrics.json`, `landmark_errors.csv`,
original_report.json and evaluation-only reference_geometry.json, plus
provenance.json. Summary fields cover:

- Counts, fractions within 2/5 mm and median/max landmark errors.
- Same-count random candidate-pair fractions and descriptive comparisons.
- Correct-pair hull area, minor/major spread and >=4 noncollinear correct-pair flag.
- Planted/estimated rotation and scale, their absolute errors, and median
  source-candidate transform-recovery error.

All 19 final cases have >=4 noncollinear <=2mm pairs. Across them, 204/207 retained
pairs are within 2 mm (98.55%) and 207/207 within 5 mm. Corresponding count-weighted
random means are 0.35% and 1.24%. These are repeated transformations of three slices
of ONE specimen, not independent subject statistics or a RANSAC-conditioned null.
The approximate source alignment remains distinct from exact synthetic geometry.

The package contains 100 landmark CSVs (1829 rows across all variants/controls,
with repeats), including 41 empty failed attempts. 34 existing CSVs were verified
and reused; missing files were reconstructed from frozen predictions. Saved
metrics were preserved and independently checked. Only missing legacy-control
chance metrics were computed, using 200 same-count unique pairings and seed 7.
No matching, fitting, thresholds or source artifacts changed.

Rebuild into a fresh directory:

```bash
.venv/bin/python scripts/export_all_landmarks.py --output results/exports/landmark_release_v2
.venv/bin/python -m unittest discover -s tests
```

34 tests pass, including export row/column/world mapping, prediction preservation
and header-only failure checks. Every exported row's world coordinates, target
pixels, score and reference error were independently verified against saved arrays
using direct matrix multiplication. Archive CRC, file inventory and every file's
bytes were checked. No benchmark thresholds changed. Preliminary candidate-only
and oracle-patch descriptor diagnostics are explicitly excluded from accepted
registration matches; the unexecuted -30 degree preparation has no invented CSV.
