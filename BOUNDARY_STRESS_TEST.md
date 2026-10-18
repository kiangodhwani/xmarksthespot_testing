# Predeclared boundary transformation stress case

User requested one very tough histology transformation within the existing
parameter ranges. Prepared one case on the accepted y=8 mm source pair:

- Rotation -30 degrees (counterclockwise).
- Scale 0.85 (15% reduction), with the existing Gaussian antialias prefilter.
- Translation -50 columns /+50 rows (10 mm left /10 mm inferior in original grid).

These values are within the declared [-30,+30] degree, [.85,1.2] scale and
inclusive [-50,+50] pixel shift ranges. This is a combined boundary condition,
not a proven worst case. In particular -30 degrees lies exactly on one of the
MIND rotation samples, and translation is less challenging with preserved full
support. No crop, flip, tear, noise, nonrigid deformation or out-of-plane tilt
has been added. No matching outcome was used to choose this case.

Preparation: `results/validation/boundary_stress_v1/pair/`.
Preview: `preview.png`; image-only matcher input: `inputs.npz`; reference
matrices, original histology and synthetic probes: separate `truth/` directory.
The three preview dots are synthetic coordinate probes, NOT detected or
anatomically established landmarks. Grey is invalid support. Expanded canvas
883x904 contains the transformed footprint without clipping. MRI is unchanged.

`study_plan.json` was saved before preparation/matching with exact parameters,
unchanged inference configurations, algorithm hashes and 2/5 mm thresholds.
`preparation_checks.json` verifies requested parameters, antialiasing, finite
arrays, unchanged MRI and algorithm hashes. Full-grid inverse error 2.27e-13 px;
world mapping error 7.11e-14 mm. All 32 existing tests pass. Assistant inspected
preview: labels fit, full anatomy visible and no clipping.

**Current stage: new transformed image awaits user visual review. No matching
has been run.** AGENTS.md requires review of prepared images before matching.
The previously accepted original source does not replace review of this new warp.

Preparation reproduction (choose another --output to avoid overwriting):

```bash
.venv/bin/python scripts/transform_pair.py --input results/pairs/coronal_y8_repeat --output results/validation/boundary_stress_v1/pair --angle -30 --scale 0.85 --shift -50 50
.venv/bin/python -m unittest discover -s tests
```

After explicit review, run unchanged `infer_frozen` from
scripts/run_two_slice_validation.py on images/support only, save predictions
before reading reference geometry, then score with the existing evaluator.
Use saved transformed-pixel-to-world mapping for mm, not 0.2 mm times target
pixel error. Report failures/incorrect lines, spatial spread and transform
recovery honestly; no retuning or replacement case selected by its outcome.

## Requested rotation-mismatch variant (current preparation)

User chose -15 degrees (halfway between the 30-degree descriptor samples),
0.85 scale and the same (-50,+50) pixel shift. Saved separately under
`results/validation/rotation_mismatch_stress_v1/`; the -30 degree case is preserved.
The new plan was recorded before image generation and no matching has run.
The new preview awaits user visual review. Canvas793x834; full-grid roundtrip
3.41e-13px /7.94e-14mm; all32 tests pass; MRI/pipeline hashes unchanged.

```bash
.venv/bin/python scripts/transform_pair.py --input results/pairs/coronal_y8_repeat --output results/validation/rotation_mismatch_stress_v1/pair --angle -15 --scale 0.85 --shift -50 50
```

## Rotation-mismatch stress result (completed)

User visually accepted the -15 degree /0.85 scale /(-50,+50) pixel case and
explicitly authorised execution (“looks good execute”). Exact preview hash and
acceptance are saved in `rotation_mismatch_stress_v1/preparation_user_review.json`.
Initial study-plan and preparation-status fields are preserved as historical
snapshots; they do not describe the completed matching stage below.

The unchanged frozen MIND source-top3/relational pipeline returned 11 matches,
all 11 within 2 mm (maximum 1.5895 mm). Median retained reference error .9190 mm.
Median transform recovery error over all 86 MRI candidates .5700 mm; absolute
rotation error .4308 degrees and absolute scale error .009066. Correct-pair
convex-hull area 5191.48 mm2; hull-area/search-mask-area ratio .6315, with
minor/major spread .6301. This is spread evidence, not dense tissue coverage.

There were 258 input proposals, 12 within 2 mm; geometric selection retained 11
of those reference-close pairs. Histology control: 117 retained, all 117 within 2mm,
median pair error .1751 mm; median whole-candidate recovery error .06338 mm.
Control counts are separate from MRI/histology results.

Both predictions were saved BEFORE either reference transform was opened.
Algorithm/input hashes, configurations and benchmark 2/5 mm thresholds unchanged.
Final assignments pass one-to-one uniqueness, >=4 noncollinear support,
cosine>=.75 and residual<6 target-pixel checks. All 32 existing tests pass.
No retuning or replacement-case selection. The earlier -30 degree preparation
remains untested, since this case replaced it at the user's request.

The stress condition passed; it does NOT prove every allowed transform succeeds,
that it was the hardest case, or that the approximate anatomical alignment is
universally within 2mm. Large translations preserve full support; shrinkage uses
antialias filtering. No crop, flip, tear, noise, nonrigid change or tilt was tested.

Reproduce from the exact accepted preparation, choosing a new output directory:

```bash
.venv/bin/python scripts/run_stress_test.py --output results/validation/rotation_mismatch_stress_v1/matching_repeat
.venv/bin/python scripts/review_stress_test.py --input results/validation/rotation_mismatch_stress_v1/matching_repeat
.venv/bin/python -m unittest discover -s tests
```

Default --study is the rotation-mismatch study. The runner checks its explicit
user-review record and exact preview hash. Do not reuse acceptance for a different
warp. Original histology for the control comes from accepted source images,
not from an inference reference transform. The same evaluator is reused.

Artifacts in `results/validation/rotation_mismatch_stress_v1/matching/`:
`summary.csv`, `reports.json`, `run_provenance.json`, `checks.json`, plus
`stress/` and `histology_control/` with frozen predictions, detector search
support, separate evaluation/CSV exports and `matches.png`. Each figure shows
incorrect INPUT proposal lines above FINAL matches. At most 12 score-selected
pairs with distinct endpoints are shown; ALL pairs contribute to metrics.
Reference colours: green<=2, amber2-5, red>5 mm. Actual descriptor settings remain
radii 16/24/36 pixels, 41x41 context samples, local 5x5 SSD, 12 rotation samples.
Assistant inspected both figures for readability/clipping; no overlap requiring
correction found. Results await user review before advancing further.
