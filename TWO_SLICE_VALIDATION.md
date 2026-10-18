# Two additional slices: frozen-pipeline validation

User requested testing the working MIND top-3 + relational pipeline on two
additional slices. Selected coronal y=-20 mm and y=35 mm before matching, because
these were already used for broad alignment inspection. They are two additional
planes of the SAME specimen, not independent subjects.

## Current stage

Preparation was numerically checked and visually accepted by the user (“yes they
look good!”). Frozen matching and reference evaluation are now COMPLETE.
All twelve new MRI/histology cases returned models; 133/133 retained pairs are
within 2 mm. Results await user review. No parameters were tuned.

The preparation review image/checks and initial study-plan status are historical
snapshots from before matching. Explicit acceptance is recorded separately in
`preparation_user_review.json`; matching status and results are under `matching/`.
Do not rerun preparation review to overwrite the acceptance record.

The frozen study plan is `results/validation/two_slices_v1/study_plan.json`.
It records the unchanged MIND/relational configuration, source-code SHA-256
values, slice selection, 0.2 mm sampling, fixed transform and five seed-42
parameter sets. Save/verify this plan before running matching. Do not tune to
new-slice results. The same transformations as y=8 control for transform setting
while testing different anatomical content. They are not new transform seeds.

Preparation output root: `results/validation/two_slices_v1/`.
Each of `y_m20/` and `y35/` contains:

- `pair/`: images.npz, geometry.json and pair.png (untransformed MRI/histology).
- `fixed/`: approved +20 degree CW /1.15x/(+40,+25) px example, preview.png,
  input images/validity, and separate truth/ metadata.
- `random_seed42/`: pair_001 through pair_005 with inputs, validity, separate
  truth/, overview.png and the batch manifest.
- `preparation_support_masks.npz`: unchanged independent MRI Otsu and histology
  1% hysteresis support for original and fixed-transformed images, for review.

Combined `preparation_review.png` shows both slices with full image context:
original MRI, original histology and fixed transformed histology. Cyan lines
are provisional search-support boundaries, NOT precise tissue segmentation or
predictions. First two columns share the original grid; third fits a different
canvas. Source intensities are preserved. Grey means invalid source support.
The -20 mm slice includes cerebellar anatomy with very different cross-modality
texture. The +35 mm slice shows small detached right/inferior components; they
have not been manually removed or used to change masks. No new alignment claim.

`preparation_checks.json` verifies both full pixel/world grids, exact y planes,
all twelve transforms, source hash metadata, MRI preservation, finite arrays,
validity shapes, no clipping and seed replay parameter equality. Both grids are
631x700, spacing 0.2 mm. Full-grid pixel/world roundtrip error 1.14e-13 pixels;
transform roundtrip errors below 1e-9 px/mm. All 32 tests pass. Frozen algorithm
source hashes match. These are numerical checks, not anatomical validation.

## Preparation reproduction

Use a fresh output study directory and write its frozen study plan first (copy
settings and hashes from the declared plan; preserve the previous artifacts).
Existing extraction/transform generators refuse to overwrite outputs.

```bash
.venv/bin/python scripts/extract_pair.py --y -20 --output results/validation/two_slices_v1/y_m20/pair
.venv/bin/python scripts/extract_pair.py --y 35 --output results/validation/two_slices_v1/y35/pair
.venv/bin/python scripts/transform_pair.py --input results/validation/two_slices_v1/y_m20/pair --output results/validation/two_slices_v1/y_m20/fixed
.venv/bin/python scripts/generate_random_pairs.py --input results/validation/two_slices_v1/y_m20/pair --output results/validation/two_slices_v1/y_m20/random_seed42 --seed 42 --count 5
.venv/bin/python scripts/transform_pair.py --input results/validation/two_slices_v1/y35/pair --output results/validation/two_slices_v1/y35/fixed
.venv/bin/python scripts/generate_random_pairs.py --input results/validation/two_slices_v1/y35/pair --output results/validation/two_slices_v1/y35/random_seed42 --seed 42 --count 5
.venv/bin/python scripts/review_validation_preparation.py
.venv/bin/python -m unittest discover -s tests
```

The preparation-review script can regenerate its own figure/masks/report inside
the supplied --study folder; it does not invoke detection or matching.

After explicit user review: run the unchanged image-only MIND descriptors,
source top 3 score>=.75 (no ratio), and relational one-to-one similarity solver
with residual<6 target px/min4 on all 12 cases. Freeze predictions before loading
references; then use the same 2/5 mm evaluator, coverage and transform recovery.
Include all failures, reference error lines, same-modality controls, and keep
outcomes separate by source slice. Review results before any modification.


## Completed frozen validation results

Final pairs within 2 mm / total retained:

| Case | y=-20 mm | y=35 mm |
|---|---:|---:|
| Fixed | 8/8 | 9/9 |
| 001 | 9/9 | 15/15 |
| 002 | 8/8 | 12/12 |
| 003 | 7/7 | 6/6 |
| 004 | 17/17 | 18/18 |
| 005 | 11/11 | 13/13 |
| Total MRI/histology | 60/60 | 73/73 |
| Histology control (separate) | 89/89 | 128/128 |

All twelve cases returned models; every retained MRI/histology pair is within
2 mm of the approximate reference, maximum 1.98889 mm. These are reference
coordinate metrics, not independently confirmed anatomical identity. Retained
pairs range 6-18 per case, not dense registration. Six transformations per plane
are repeated observations of the same anatomy, not twelve independent slices.

Whole-MRI-candidate median transform recovery errors:

| Case | y=-20 mm | y=35 mm |
|---|---:|---:|
| Fixed | 0.867 mm | 0.653 mm |
| 001 | 1.106 mm | 0.351 mm |
| 002 | 0.919 mm | 0.593 mm |
| 003 | 1.226 mm | 0.463 mm |
| 004 | 0.789 mm | 0.301 mm |
| 005 | 0.929 mm | 0.330 mm |

These measure transform recovery over detected source candidates, not the whole
pixel field, anatomical registration error or a worst-case bound. Retained-pair
median errors range 1.055-1.324 mm at -20 and 0.519-0.792 mm at +35.
Correct-pair convex hull area / MRI search-mask area ranges .489-.884 and
.379-.849 respectively. This hull-area ratio is a spread diagnostic, not a
bounded tissue-coverage fraction: convex hulls bridge holes/nonconvex boundaries;
in controls it can exceed 1. It is preserved from the existing evaluator.

All fourteen predictions, including controls, were saved before ANY reference
geometry was opened for evaluation. Input and frozen-algorithm hashes were
verified before/after, and MIND/relational configs match the predeclared plan.
The new runner was first checked against y=8 fixed: candidates, scores, selected
bank indices, retained pairs and model reproduce exactly. Final outputs pass
one-to-one source/target, min4 noncollinearity, finite similarity model,
cosine>=.75 and residual<6 target-pixel checks. All 32 existing tests pass.
No inference, geometry, descriptor, detector, masking or evaluator thresholds
were changed. No new dependency/GPU. SuperGlue was not needed for this stage.

### Matching reproduction

After review of the exact preparation being used:

```bash
.venv/bin/python scripts/run_two_slice_validation.py --output results/validation/two_slices_v1/matching_repeat
.venv/bin/python scripts/review_two_slice_validation.py --input results/validation/two_slices_v1/matching_repeat
.venv/bin/python -m unittest discover -s tests
```

The runner reads the existing plan and verifies frozen hashes. It requires the
explicit user-review record; do not manufacture acceptance for different images.
Keep the immutable preparation/plan and use a fresh --output for repeated runs.
The review script only renders frozen predictions and evaluation reports; it
does not refit or select models.

Artifacts in `results/validation/two_slices_v1/matching/`:

- `run_provenance.json`, `checks.json`, `final_assignment_audit.json`.
- `summary.csv`, `reports.json`, `slice_summary.json`, `comparison.png`.
- Each `{y_m20,y35}/{case}/`: image-only `predictions.npz`, `prediction.json`,
  search_support.npz, separate `evaluation.json`,
  `all_proposals_evaluation.csv`, `landmarks.csv`, and `matches.png`.
- Control CSVs explicitly label original-histology world coordinates rather than
  falsely calling them MRI coordinates. Histology controls use untransformed
  original histology as a diagnostic source input, without reading truth/ images.

Each full-image figure shows INPUT proposals (including incorrect red lines)
above FINAL one-to-one matches. Reference colours are evaluation-only:
green<=2 mm, amber2-5, red>5. Up to twelve highest-score pairs with distinct
endpoints are displayed, skipping duplicate endpoints for readability. This
is a DISPLAY choice, not a proposal filter; all numerical metrics use all pairs.
Numbers index pairs within each row. Existing MIND context documentation records
the unchanged scoring-patch sizes and 3x3 radius/12-rotation search.

Assistant inspected comparison, both fixed-example figures and -20/pair003;
labels and summaries fit. Current evidence supports success on two additional
slices of the same specimen with reused transformation settings. It does not
establish new-subject generalisation, arbitrary section damage, nonrigid recovery,
unknown out-of-plane matching or universally <=2 mm anatomical alignment.
Review this result before advancing to further validation or model changes.
