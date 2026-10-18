# MRI–histology foundations

## Scope and gates
1. Use the released aligned BigBrain MRI and Histological Space histology.
   We do not estimate a new registration. Header checks establish coordinate
   compatibility, not biological accuracy; visually review shared anatomy.
2. Resample both volumes onto exactly the same coronal world plane and grid.
3. Transform histology, retaining a complete invertible point mapping.
4. Review outputs before implementing detection or matching.

## Coordinate contract
- Source voxels are (i,j,k); NIfTI affine maps them to RAS world millimetres.
- Output images are arrays [row,col]. Point vectors are (col,row,1).
- col increases toward anatomical right; row increases downward (inferior).
- Grid pixel centres: x=x_min+spacing*col, y=chosen_y,
  z=z_max-spacing*row. Display with origin='upper'. This is anatomical right
  on screen right, not the conventional radiological display.
- Use the full inverse affine to sample source volumes. Never use equal slice
  indices, assume array axes, or round to different native planes.
- Histology transformation T maps original grid points to transformed points.
  Image resampling uses inverse(T). Save both matrices and the canvas offset.
- Ground-truth transformed histology -> MRI grid is inverse(T); then use the
  saved grid-to-world matrix for millimetres. Never measure error in pixels.

## Mandatory checks
- Finite, invertible affines; spatial units mm; orientations, spacings, extents,
  source hashes and NIfTI intensity scaling recorded.
- Chosen plane lies inside both volumes. Sampling validity masks saved; invalid
  regions excluded from evaluation. Valid support does not imply tissue.
- Both modalities use the same plane, field of view, spacing and dimensions.
- Identity sampling and analytic world-coordinate phantom checks pass.
- Identity, translation, rotation, scale and combined point round trips pass.
- Image-warp direction checked with an asymmetric impulse phantom.
- Expanded canvas contains every transformed source pixel footprint; no silent
  clipping. Record effective transform after canvas expansion. Translation may
  change margins, so do not treat canvas centring as anatomical information.
- Display normalisation is for PNGs only; numerical arrays retain scaled values.
- Inspect aligned side-by-side views and edge overlays at multiple planes.
  Header agreement and synthetic checks do not prove anatomical alignment.

## Matching and evaluation pitfalls (later stages)
- Fine histology features may have no MRI counterpart; compare physical scales.
- Contrast is modality dependent; inversion does not solve tissue differences.
- Open/closed sulci and processing deformation can move apparent landmarks.
- Bilateral symmetry and repetitive folds cause plausible false matches.
- RANSAC consistency is not anatomical truth. Permit unmatched candidates.
- Report precision within 2/5 mm, counts and spatial coverage, not only inliers.
- Use tissue-restricted, comparable random baselines and held-out slices/seeds.
- Keep affine, planted transform and truth points out of matcher inputs.
- BigBrain is one specimen; synthetic transforms do not establish generalisation
  to new subjects, damaged sections or unknown 3D planes.
- Released alignment is a reference, not perfect microscopic ground truth.
- Keep source data immutable; record interpolation, versions and parameters.
- Derived BigBrain imagery retains the data attribution/licence in README.

## Running and viewing: current inspection stage
From this directory:
```
.venv/bin/python scripts/inspect_volumes.py
```
Open results/inspection/native_coronal_preview.png and headers.json.
This script shows nearest native slices for orientation and appearance only.
They differ slightly in world y and MUST NOT be used as exact evaluation pairs.
The native inspection script performs no registration or common-grid extraction.
For the subsequently implemented experiment extraction, see the current status below.

## Inspection completed 2026-10-18
Both downloaded volumes load with the README's expected dimensions:
MRI (315,362,284), spacing 0.444444 mm; histology (696,770,605),
spacing 0.2 mm. Both are axis-aligned RAS in mm with invertible affines.
Histology scaling is slope 257, intercept 0, observed scaled range 0–65535.
MRI range is 0–249.880844. Sources and SHA-256 hashes are in headers.json.
Preview includes native planes near y=-20,8,35 mm. Broad anatomy is visually
consistent, but residual local differences remain; this is not a quantitative
validation of the released registration. At requested y=8, native MRI y is
7.985556 and native histology y is 8.089998 (difference 0.104441 mm).
No common-plane resampling or synthetic transformation has been performed.

## Alignment review update
Read PROJECT_PATH.md for current stage gates and evidence. Review-only common
plane sampling now exists in scripts/review_alignment.py, at 0.4 mm spacing.
The production shared-plane extraction is implemented as described below;
synthetic transformation remains pending. The official
BigBrain update log explicitly warns the prescribed MRI does not fit perfectly.
Anatomical residual error remains unmeasured until independent annotations.

## Shared-plane extraction status (current)
Alignment review closed as preliminary at user's request; broader anatomical
uncertainty remains. Run `python -m unittest discover -s tests` for geometry
checks. scripts/extract_pair.py saves a y=8 mm pair at 0.2 mm spacing; see
PROJECT_PATH.md for outputs and reproducibility. Numerical arrays preserve
scaled intensities, display uses original grayscale, invalid margins are masked.
Grid is 631x700; [col,row,1] maps to RAS world via the saved 4x3 matrix.
No downsampling or antialias filtering needed for these two prescribed sources;
finer inputs are rejected rather than silently downsampled. No synthetic
transformation or matching exists yet. User and assistant visually reviewed the saved pair; extraction accepted.
Next is controlled histology transformation, not yet implemented.

## Transformation status (current; supersedes earlier pending notes)
landmarks/transforms.py and scripts/transform_pair.py implement the approved
20-degree-clockwise / 1.15-scale / (40,25)-pixel translation example.
See PROJECT_PATH.md for output paths and checks. Effective T includes a canvas
offset; use truth/geometry.json, never only the requested rotation and shift.
Bilinear sampling validity requires all contributing samples valid. Grey preview
margins are invalid support; their shape is not anatomy. Future matchers must
not use artificial rectangular support/canvas edges as anatomical landmarks.
No downscaling is used here; anti-aliasing must be added before supporting
substantial scale reductions in image generation. Synthetic point mappings are
invertible; a forward/back image warp need not reconstruct intensities exactly.
Seven geometry tests pass and real full-grid coordinate checks pass. Assistant
inspected preview; user review is pending. Algorithm work remains deferred.

## Masking contract and current finding
Image-derived Boolean tissue masks are separate from sampling-validity masks.
Same image shape/index coordinates; no crop, move or resize. Always intersect
with valid support; histogram thresholds exclude invalid padding. Compute each
modality independently, without truth geometry or supplied labels. Retain raw
masks and log all cleanup. Never fill holes or bridge sulci indiscriminately.
Do not equate the largest connected component with all legitimate tissue.
BigBrain supplied labels are a potential validation reference, with semantics
and quality to inspect first; do not use them to force MRI/histology agreement.

Otsu baseline implemented and tested (see PROJECT_PATH.md). MRI looks plausible
on this slice; histology loses pale tissue and remains unaccepted. Cleanup
cannot recover threshold-excluded tissue. Pixel-count threshold sensitivity
is stability evidence, not accuracy. Next step is user review, not detection.

## Official reference-mask semantics established
full_cls_200um_mask.nii.gz is a brain-envelope mask on our reviewed slice: it
fills ventricles and many sulci. Do not score tissue-only masks as though all
reference foreground were tissue. Otsu Dice 0.864 against it is envelope
agreement only. Label values are scaled to 0 / ~0.00392; >0.5 would wrongly
produce an empty mask. Use nearest-neighbour for categorical reference labels,
respect full source affine and saved transformed-pixel mapping, and exclude
invalid support. These reference labels remain evaluation-only.

## Hysteresis experiment update
Background-relative hysteresis now implemented and reviewed by assistant; user
review pending. See PROJECT_PATH.md. Preserve distinction between low-gradient
white background, pale tissue and genuine non-tissue spaces. Near-white
intensities can be ambiguous. Brightest-decile background assumption is suitable
only if background is actually represented there; test other slices before
claiming robustness. Unknown synthetic scale must not be supplied to masking.
Fixed relative intensity margins are initial experimental settings, not learned
or proven universal thresholds. No BigBrain labels enter mask generation.
Do not attribute all improvement over Otsu to connectivity: most improvement in
this example is threshold choice; almost all possible pixels connect to seeds.

## Accepted provisional support (2026-10-18)
Use histology key weak_0.01 (1%) for provisional search support, not the 0.5%
default hysteresis mask. MRI Otsu is provisional. Grayscale images retain internal
structure; masks need not encode every visible line. Do not use the approximate
mask alone to claim reliable internal groove landmarks. Candidate detection is
next, with validation separate from matching and truth excluded from detection.

## Candidate baseline status (current)
Multiscale Harris implemented; see PROJECT_PATH.md. Tissue mask only filters
centres; original grayscale is not multiplied by mask. Discard detector support
near artificial invalid edges. Detection uses no truth or counterpart data.
Eval one-to-one proximity is repeatability/candidate coverage, NOT landmark
matching accuracy; report chance and counts. Baseline on one slice is promising
but changing smoothing scales changes many points. User review pending.

## Peak selection order
Find local maxima in the true grayscale response BEFORE applying tissue/support
eligibility. Masking a response map can create artificial maxima at mask edges.
Regression test now covers this. Current candidates are in harris_y8_peakfix;
previous counts/results are superseded. Inspect exact centre dots, not ring
footprints; provisional mask membership is not proof of anatomical tissue.

## Descriptor baseline status
Polar unsigned-gradient descriptor implemented; see latest PROJECT_PATH.md.
Known-location descriptor validation must remain separate from detected-candidate
matching. Truth may choose test patch centres but must not align their orientation
or provide scale to the descriptor. Reference-positive patches remain approximate
across modalities. Nearby negatives excluded in initial test; do not present top1
retrieval as final landmark precision. Global inversion tolerance does not imply
modality independence. Strong histology-only test and weak MRI/histology result
require review before selecting the next representation/matching approach.

## Random-pair generation status
README Day 1 random similarity examples now generated; see PROJECT_PATH.md.
Shrinking images applies Gaussian antialias prefilter and conservative support
exclusion. This supersedes prior notes saying reduction filtering is pending.
Forward/inverse coordinate mappings include canvas offsets. Store seed AND
sampled parameters/versions; same-seed replay checked on all five saved images.
Random images from the same source slice are augmentations, not independent
anatomical examples. Do not count them as independent specimens or randomly
split them across future training/test sets. Source-plane splits come first.

## Matching baseline status (current)
Complete first descriptor-match/RANSAC/reference-evaluation pipeline implemented.
See PROJECT_PATH.md for exact thresholds, outputs and weak results. Freeze
predictions before reading truth. Max descriptor score is not confidence; RANSAC
consensus is not anatomical truth. Return no model if too little evidence; avoid
misleading 100% precision on empty outputs. Label oracle-patch tests separately
from actual detected-candidate matching. Reference-conditioned point geometry
is used ONLY to score and colour lines, not to select matches or fit the model.
Report failures across all random cases, not just the one successful example.

## MIND-style comparison status
One fixed sparse 2D self-similarity adaptation is implemented; see
MIND_EXPERIMENT.md for exact method and limitations. It is not original dense
3D MIND/MIND-SSC. Keep rotation sampling angles distinct from recovered transform
angles. All 3x3 radius pairs and 12 rotations are searched without planted scale,
orientation or affine. Gaussian prefilter and interpolation support are validated.
More reference-correct initial proposals did not improve final <=2mm matching.
Benchmark thresholds and evaluator remain unchanged. Review before retuning;
score distributions across descriptors are not calibrated equivalently.

## Proposal-filter ablation update
User-authorised ablation shows ratio filtering discards useful MIND evidence:
removing only ratio increases model returns 1/6->4/6 on the same reviewed batch.
See FILTER_ABLATION.md. Score>=.75 removes no mutual proposals here, so no-ratio
and all-mutual variants coincide. Keep original default intact pending review.
All variants freeze predictions before reference scoring. Additional consensus
is assessed via pair errors, coverage AND whole-field transform recovery; e.g.
pair001 retained pairs look good but median candidate-field error remains 2.126 mm.
This is exploratory reuse of one slice, not held-out performance validation.

## Expanded proposals and relational selection
Fixed source top-3 MIND alternatives plus one-to-one similarity consensus now
returns models in all 6/6 augmentations; see RELATIONAL_EXPERIMENT.md. This is
new hypothesis/assignment logic with the same residual/minimum-support thresholds,
not unchanged RANSAC. Keep source/target uniqueness: alternatives are hypotheses,
not independent landmarks. Rank metrics count each source once. Input many-to-many
proposal counts cannot use the baseline unique-pair random sampler; final
one-to-one pairs retain that comparator. Pruning is heuristic, not global proof.
All predictions precede reference evaluation; diagnostic expected locations never
enter candidate galleries. Case004 recovery degrades despite extra matches;
review full-field error as well as pair counts. Single-slice validation remains.

## Two additional slices tested with frozen settings
After user accepted preparation, y=-20 and y=35 mm were evaluated with the
unchanged top-3 MIND/relational pipeline. All 12 cases return models; all 133
retained cross-modality pairs are within 2 mm of reference. See
TWO_SLICE_VALIDATION.md for controls, recovery/coverage, provenance and scope.
This extends evidence to two new anatomical planes, still one specimen and the
same transform settings. Do not equate retained-pair accuracy with dense image
accuracy or verified anatomy. Convex-hull area divided by mask area is a spread
ratio that may exceed 1, not a bounded tissue-coverage fraction. All predictions
precede reference evaluation; no thresholds/settings tuned to these slices.
