# Project path and working agreement

## Goal
Build a reproducible two-day proof of concept that finds corresponding anatomical
landmarks between known corresponding MRI and histology sections. Current work
is ONLY the foundations: reference alignment, shared-plane extraction, then
controlled histology transformation with retained geometry. Matching comes later.

## Working agreement
Work one stage at a time. Explain the purpose, show the output and checks, and
review it with the user before advancing. Read this file and FOUNDATIONS.md at
the start of project work. Keep the status honest: proposed, implemented,
validated numerically, and reviewed anatomically are different states.

## Stages and acceptance gates
1. Reference alignment assessment — CLOSED AS PRELIMINARY (user accepted).
   Use the released aligned MRI; no new registration is currently planned.
   Confirm headers and provenance, compare anatomy at multiple locations,
   and measure independently selected corresponding landmarks in world mm.
   Report regional errors and annotation uncertainty. Do not claim <=2 mm
   alignment from compatible headers, attractive overlays, or voxel size.
2. Exact shared-plane extraction — COMPLETE; numerically checked and visually reviewed by user and assistant.
   Save a common physical plane and grid with explicit image-to-world mapping.
   Validate against analytic synthetic data and inspect anatomical outputs.
   A review image may use common-plane sampling to assess stage 1, but that
   does not mean the production extraction stage is complete.
3. Controlled histology transformation — COMPLETE; user and assistant visually reviewed.
   Rotation, scale, translation, explicit forward/inverse matrices, expanded
   canvas and validity masks. Check known points, round trips and warp direction.
   Show before/after with known correspondence. Flips/crops/tilt are deferred.
4. Candidate detection — BASELINE IMPLEMENTED, USER REVIEW PENDING; matching DEFERRED.
   Start only after stages 1–3 are reviewed. Build and validate incrementally.
5. Evaluation — DEFERRED.
   Report errors in mm, fractions <=2/5 mm, count, coverage and fair chance
   baseline on held-out slices/settings. Separate numerical transformation
   truth from uncertainty in the released anatomical alignment.

## Current evidence (2026-10-18)
- Downloaded both prescribed volumes; source hashes and headers recorded.
- MRI shape (315,362,284), spacing 0.444444 mm; histology shape
  (696,770,605), spacing 0.2 mm. Both axis-aligned RAS in mm.
- Native-slice preview generated near y=-20,8,35 mm.
- Broad anatomy appears consistent. Residual anatomical alignment error is
  NOT measured. No verified 2 mm accuracy claim exists.
- 2 mm and 5 mm are README evaluation thresholds, not measured alignment bounds.
- Synthetic transformation is now implemented (see current transformation handoff); matching is not implemented.

## Optional stronger alignment assessment protocol
Use at least 10 clear anatomical point pairs distributed over both hemispheres,
internal structures and cortex, across at least three planes. Select each image's
point independently; do not generate its partner from shared coordinates. Use
precise anatomical definitions and repeat annotation (preferably a second rater)
to estimate selection variability. Record ambiguous/unidentifiable structures.
Compute Euclidean world-mm distances; report all values, median, mean, maximum,
and fractions <=2/5 mm. These describe sampled landmarks, not a global bound.
Cross-sectional ambiguity and out-of-plane residuals require 3D inspection;
2D landmark errors alone can underestimate full 3D error. If uncertain, retain
reference-alignment error as unresolved rather than manufacturing a number.

## Entry points
- scripts/inspect_volumes.py: native slices and header inspection.
- results/inspection/: initial preview, headers and provenance.
- FOUNDATIONS.md: coordinate conventions and mandatory pitfall checks.

## Alignment evidence and review tools (2026-10-18)
Official release notes, Aug 3 2017, explicitly state that the MRI in 3D_MRI
"does not fit perfectly" and recommend it over Raw_Data MRI:
https://ftp.bigbrainproject.org/bigbrain-ftp/BigBrainRelease.2015/Updates.txt
This confirms residual error but provides no numerical bound.
A primary study also calls BigBrain's MRI–histology alignment approximate:
https://pmc.ncbi.nlm.nih.gov/articles/PMC6742511/ (section 3.1.3).
Do not substitute error numbers from registration to MNI/PD25 or other specimens.

Run `.venv/bin/python scripts/review_alignment.py` for exact-plane inspection
at y=-20,8,35 mm with a 0.4 mm review grid and linear interpolation. Original
scaled arrays and validity masks are retained; histology is inverted ONLY in
the display. Analytic sampling test maximum intensity error was 3.55e-15;
this tests coordinate sampling, NOT anatomical accuracy.

Run `.venv/bin/python scripts/review_alignment.py --annotate` locally with an
interactive Matplotlib backend. Define a landmark, select its position in MRI,
then independently in histology; repeat across planes. Enter blank to advance.
Writes results/alignment_review/manual_landmarks.csv and report.json.
Manual GUI interaction has not yet been exercised. No anatomical annotations
have been supplied and no numerical anatomical error has been established.

## Annotation window troubleshooting and persistence
User successfully recorded points, then reported `SystemError: null argument to
internal routine` using Matplotlib 3.11.2 / macosx. No traceback was available;
backend involvement is a hypothesis, not a confirmed root cause.
Try `python scripts/review_alignment.py --annotate --backend TkAgg` from the
repository root. Tkinter 8.6 imports successfully in the project environment;
GUI operation with Tk still requires a user-side check.
Each run now writes its own timestamped manual_landmarks_*.csv, flushing each
completed pair immediately. Prior files are not overwritten. Summaries go to
manual_report_*.json. This supersedes the earlier fixed manual CSV/report paths.
An incomplete pair cannot be saved. Keep descriptive anatomy labels and exclude
random practice points. Current user-reported three points at y=-20 have
separations 1.23993, 0.80950, 0.22452 mm; labels are generic, anatomical identities
are unverified and additional planes/repeats are pending.

## Current handoff: shared-plane experiment pair
User chose to close alignment review as preliminary and proceed. The reported
ambiguous y=8 point (~5.02 mm separation) is NOT a verified anatomical pair;
record it as counterpart not clearly visible, excluded from confirmed errors.
Further manual annotation is not a requirement for moving forward.

Implemented scripts/extract_pair.py and landmarks/planes.py. Run:
```
python scripts/extract_pair.py --y 8
```
The y=8 output already exists; the script refuses overwrites. To reproduce into
another directory use `--output results/pairs/coronal_y8_repeat`.

Saved results/pairs/coronal_y8/ contains:
- images.npz: MRI, histology, separate validity masks and joint validity.
- geometry.json: shared grid, 4x3 pixel-centre-to-world mapping, source affines,
  hashes, intensity scaling, package versions and numerical check results.
- pair.png: untransformed pair for visual review.

Both arrays are 631 rows x 700 columns on an exact y=8 mm coronal plane,
spacing 0.2 mm. The review-only preview used 0.4 mm; the experiment uses 0.2 mm
so neither source is spatially downsampled. MRI upsampling adds no information.
Images preserve modality intensities; no inversion or matching preprocessing.
Union field of view includes invalid margins, explicitly masked. Valid support
is not a tissue segmentation. Original volumes remain unchanged.

Three tests pass: identity native-slice extraction/orientation, analytic linear
world signal through an oblique source affine plus inverse point mapping, and
out-of-volume validity/plane rejection. Saved real arrays reload exactly;
pixel-coordinate round-trip maximum error 7.11e-15 pixels. These are numerical
checks, not a biological alignment claim. Next: user reviews pair.png, then
controlled histology transformation with saved ground truth. No transform or
matching algorithm has been implemented.

## Extraction acceptance check (2026-10-18)
User inspected coronal_y8_repeat/pair.png and accepted its appearance.
Assistant independently inspected the repeat preview and verified that all saved
arrays and geometry metadata exactly match coronal_y8. All 441700 grid pixels
are at y=8 mm; x increases right and z decreases down by 0.2 mm per pixel.
Full-grid pixel/world round-trip maximum error is 1.14e-13 pixels. Separate
validity masks and their intersection are consistent; finite arrays and all
three synthetic geometry tests pass. Valid-source support fractions are
MRI 0.993982 and histology 0.950374; these do not measure tissue coverage.
Extraction stage complete. Next stage is controlled histology transformation;
not yet started. Released anatomical alignment uncertainty remains unchanged.

## Current handoff: synthetic transformation (2026-10-18)
Applied the user-approved fixed example to the accepted coronal_y8_repeat pair:
20 degrees clockwise around original image centre, scale 1.15, then translate
40 columns right / 25 rows down. This is 8 mm right / 5 mm inferior in the
original grid. One bilinear image resampling step. Full source footprint fits
inside a 1029-row x 1106-column output canvas; original data remain immutable.

Entry point: `python scripts/transform_pair.py` (already executed).
To reproduce without overwriting use
`python scripts/transform_pair.py --output results/transforms/coronal_y8_demo_repeat`.
Default input is results/pairs/coronal_y8_repeat; override with --input.

Output results/transforms/coronal_y8_demo/:
- inputs.npz: unchanged MRI, transformed histology, separate validity masks.
- preview.png: common display canvas showing MRI/original/transformed histology.
- truth/geometry.json: requested and effective forward/inverse mappings, canvas
  offset, centre, source hashes, original and transformed pixel-to-world mapping.
- truth/original_histology.npz: original histology and its support mask.
- truth/coordinate_probes.csv: three synthetic probe mappings; NOT anatomical
  correspondences or algorithm predictions. Kept out of future matcher inputs.

Canvas includes the original and zero-translation rotated/scaled footprints,
plus symmetric translation allowance and 10 px padding. It is not a tight crop:
translation remains visible relative to the padded original image. Preview
positions MRI on this display canvas only; saved MRI pixels/shape are unchanged.
Canvas dimensions/offset depend on parameters and are evaluation metadata.

Seven tests pass: three plane tests and four transform tests, covering identity,
asymmetric image translation, clockwise image rotation, scale and combined point
roundtrips, negative shifts, footprint containment and conservative validity.
Full real-grid roundtrip maximum errors: 2.27e-13 px and 7.14e-14 mm.
MRI equality, saved transformed-array equality, and source hashes checked.
Visual preview inspected by assistant; user visual review remains pending.
Image resampling is not lossless even though coordinate mappings are invertible.
Scale enlargement changes underlying physical sampling to 0.2/1.15 mm along
rotated axes; use saved transformed-pixel-to-world mapping for distances.
No landmark detector or matcher has been implemented.

## Current stage: independent tissue masking baseline (2026-10-18)
User accepted transformation preview. Prep stages 1–3 complete with stated limits.
User chose image-derived masks for BOTH modalities; BigBrain histology labels
are for validation only, not threshold selection or matcher input. Supplied mask
has not yet been downloaded or its exact labels inspected.

Implemented landmarks/masking.py and scripts/mask_pair.py. Run:
`python scripts/mask_pair.py --output results/masks/coronal_y8_otsu_repeat`.
Default results/masks/coronal_y8_otsu already exists; refuse overwrites.
Reads ONLY images/support in inputs.npz; no truth matrices, counterpart image
information, or supplied labels enter either mask calculation. Outputs:
masks.npz (raw, cleaned, threshold-sensitivity Boolean masks per modality),
review.png, sensitivity.png, histograms.png and report.json.
Shapes/coordinates unchanged; source intensities and input file hash unchanged.
Thresholds use 256-bin Otsu on valid samples; MRI >84.617752, histology <57168.7734.
Cleanup removes only <=25-working-pixel 8-connected components; no hole filling,
closing, erosion or smoothing. Pixel-size cleanup is a diagnostic setting,
not physically calibrated or assumed invariant to scale. Raw masks retained.
MRI raw/clean tissue pixels: 205536; cleanup removes 0.
Histology raw 224539 / clean 224120; cleanup removes 419 pixels.

Assistant visual finding: MRI plausible as a starting mask on this slice;
not quantitatively validated. Histology Otsu excludes substantial visible pale
interior tissue; not accepted for landmark extraction. +/-2% of valid intensity
range changes 886/939 MRI pixels and 19859/16953 histology pixels. Sensitivity
is not an accuracy score. No reference Dice/error claim has been made.
Nine tests pass, including masking checks for padding exclusion, threshold
independence from padding, polarity, cavity retention, speck handling and source
preservation. Mask shape/dtype/subset checks and saved reload checks passed.
Next: review masks with user, then decide a histology alternative; do not proceed
to landmark detection on the failed baseline. Potential next experiment is
background-based thresholding or hysteresis, retaining pale connected tissue;
not yet implemented or selected as successful.

## Official mask comparison (2026-10-18)
User requested checking Otsu against full_cls_200um_mask.nii.gz. Downloaded from
https://bigbrain-ftp.loris.ca/bigbrain-ftp/BigBrainRelease.2015/3D_Classified_Volumes/Histological_Space/nii/
Source hash a5f55a76a5a01be256c4283a9f89b01e08afda8c3bfabc61e066b5118b9b24a7.
Shape, 0.2 mm spacing and affine match histology. Scaled labels on sampled slice
are 0 and 0.003921568859; foreground is >0, not >0.5! Stored scaling must be
honoured. This mask is a brain ENVELOPE: visual inspection shows filled
ventricles/sulci/internal spaces. It is not a tissue-only reference.

scripts/compare_reference_mask.py samples categorical reference labels with
nearest-neighbour directly at inverse-transformed target pixel centres using
the saved geometry. Native categorical y=8.089998 mm is nearest to requested
y=8 mm; this unavoidable label discretisation is documented. Invalid support
excluded. Mask kept evaluation-only, independent mask algorithm unchanged.
Outputs results/masks/coronal_y8_reference_check/: comparison.png, report.json,
reference_comparison.npz. Script refuses overwrites; use --output for repeats.

Raw Otsu vs envelope: Dice 0.864316, IoU 0.761054; 224539 foreground-agreement
pixels, 70498 envelope-only, 0 Otsu-only. This means Otsu covers 76.1% of the
reference envelope, NOT that it misses 23.9% of true tissue. Disagreement
combines pale tissue omissions, correctly excluded internal spaces, and boundary
convention/interpolation differences. Cleaned-mask Dice 0.863400. Envelope
comparison cannot establish tissue-only accuracy. The 9-class volume may give
a more suitable reference after inspecting documented label meanings; it has
not been downloaded or used. Next: review with user and decide target mask
semantics before modifying histology segmentation. No detection/matching yet.

## Current experiment: background-based histology hysteresis (2026-10-18)
User authorised independent hysteresis and a panel showing original, Otsu,
hysteresis and BigBrain reference. Implemented histology_hysteresis in
landmarks/masking.py and scripts/try_hysteresis.py. No landmark detection yet.
Run from repository root:
`python scripts/try_hysteresis.py --output results/masks/coronal_y8_hysteresis_repeat`.
Default results/masks/coronal_y8_hysteresis already exists (refuse overwrite).

Prediction uses transformed histology + valid support ONLY. Estimate white
background from median of brightest valid decile (65535). Estimate bright-tail
MAD noise (0 in this processed slice); this is not a calibrated noise model.
Strong seed darkness >5% of valid intensity range (threshold I<62338.3533).
Possible tissue darkness >0.5% range (I<65215.3353), with 3*MAD noise floor.
Propagate strong seeds within possible tissue using 8-connectivity. No hole
filling, closing or cleanup. Fixed initial parameters NOT optimised to BigBrain.
MRI Otsu result and all source image data remain unchanged.

Saved masks.npz, comparison_panel.png, sensitivity.png and report.json.
Panel includes original transformed histology, Otsu, hysteresis, official brain
envelope, hysteresis boundary overlay and added/removed-pixel overlay. Only
plot axes are cropped, identically across panels; saved mask coordinates and
array shape are unchanged. Reference is read only after prediction generation.

Result: 261876 tissue pixels vs 224539 Otsu (+37337, no losses). Strong seeds
already include 253445 pixels, with 8431 pale connected pixels added by
propagation. Only 26 possible pixels rejected by connectivity: the primary
improvement is the background-relative thresholds, not hysteresis alone.
Weak-margin 0.25% gives 263213 pixels; 1% gives 260223 (1337/1653 pixels differ
from default). Assistant visual inspection: large pale interior area recovered,
major cavities retained; no obvious broad white-background leak. Narrow sulci
and boundary interpolation remain possible weaknesses. User review pending;
no tissue-only reference validation or cross-slice generalisation established.
Envelope Dice improves from 0.8643 to 0.9405, but this is NOT tissue accuracy:
the reference intentionally includes internal spaces. Do not tune to maximise it.

Eleven tests pass. New tests check pale connected tissue inclusion, cavity
preservation, isolated weak patch exclusion, input preservation, and invariance
to invalid padding. Runtime mask/support/shape checks, source SHA preservation
and saved-mask reload check pass. Next: user reviews this panel; only then
consider accepting mask baseline or revising parameters. No matching stage yet.

## User acceptance and next stage (2026-10-18)
User accepted the 1% hysteresis variant as PROVISIONAL HISTOLOGY SEARCH SUPPORT,
not precise or quantitatively validated tissue segmentation. Selected saved
array: results/masks/coronal_y8_hysteresis/masks.npz key weak_0.01 (260223 pixels).
Default hysteresis key uses 0.5%; do not accidentally substitute it for 1%.
MRI baseline remains results/masks/coronal_y8_otsu/masks.npz key mri_cleaned.
Internal grooves may be overincluded. Preserve grayscale intensities and detect
internal features from image content, not only this approximate mask contour.
Full-size images and side-by-side original/1% preview are in the hysteresis
output folder. No geometry or image data changed during mask visualisation.

Next stage: propose a grayscale candidate detector, then implement only that
stage after agreeing the method. Generate candidate overlays and validate first
on original vs transformed histology to isolate geometric repeatability, then
MRI vs histology for cross-modality repeatability. Transform geometry is used
by evaluator only. Candidate repeatability via nearby points is not anatomical
matching accuracy; compare against density/chance and report coverage/counts.
No descriptors, candidate detector, matching, RANSAC or final evaluation exist
yet. Later stages: describe candidates, propose matches allowing unmatched
points, reject inconsistent pairs, evaluate <=2/5 mm reference errors, coverage,
transform recovery and chance baseline, then test held-out slices/transforms.
DINOv2 is an optional descriptor comparison, not currently installed or used.
Continue showing and reviewing one stage at a time with the user.

## Current handoff: multiscale Harris candidate baseline (2026-10-18)
User authorised detector experiment. Implemented landmarks/candidates.py,
landmarks/repeatability.py and scripts/detect_candidates.py. Current complete
output is results/candidates/harris_y8_checked/ (harris_y8 is first equivalent
run, without sensitivity.json). Reproduce with:
`python scripts/detect_candidates.py --output results/candidates/harris_y8_repeat`.
No descriptor matching, RANSAC or anatomical labelling is implemented.

Detection: independent valid-pixel 1st/99th-percentile normalisation; Harris k=.04
at sigma 2/4/8 working pixels, with Gaussian presmoothing sigma/2. Candidate
centres filtered by MRI Otsu / histology 1% provisional support; original histology
mask independently recomputed using identical hysteresis method. Grayscale is
never multiplied by tissue mask. Exclude footprint near invalid samples/image
edges, remove duplicates within 8 px across scales, cap 180, quota 10 per tile
on independent 6x6 tissue-bounding-box grid. Scores normalised per scale/image;
not calibrated confidence. Unknown transform is not supplied to detector.

Outputs: mri.csv, histology_original.csv, histology_transformed.csv (col,row,
relative_response,sigma,tile); candidates.png; detector_report.json;
repeatability.json; sensitivity.json. Source files hash-preserved. All images
produce 180 candidates, spanning 34/35/30 of 36 bbox tiles respectively.
Assistant inspected preview: internal and contour candidates, no obvious frame
corners; visual anatomical quality remains to be reviewed with user.

Evaluation-only known geometry maps points to original 0.2 mm grid. One-to-one
maximum-cardinality proximity assignments at 1/2/5 mm; BOTH sets restricted to
counterpart valid support. Not descriptor matches or proof of anatomical identity.
Histology original/transformed: 152/180 within 1 mm (84.4%, random 7.1%);
156/180 within 2 mm (86.7%, random 23.1%). MRI/transformed histology: 36/180
within 1 mm (20.0%, random 6.1%); 90/180 within 2 mm (50.0%, random 22.0%).
MRI/original histology within 2 mm: 95/180 (52.8%, random 23.0%).
At 5 mm chance is high (~66%); nearby candidate does not establish a correct
match. These are on ONE slice/transform with approximate anatomical reference.
Baseline uses 50 seeded independent uniform tissue samples with same counts
and 8px spacing / valid border exclusion. Does not reproduce Harris tile
occupancy; interpret as an initial density-controlled comparator, not a perfect
null model or significance proof.

Sensitivity: halving/doubling response threshold (.005/.02 vs .01) leaves capped
selected sets within 2px unchanged here. Changing scales changes many points:
fine-only retains 44–60% near default; coarse-only 77–88%. Scale sensitivity
is unresolved; do not claim detector universally stable. Full scale vs transform
sweep and held-out slices are deferred.

14 tests pass, including synthetic corner/coordinate extraction, absence of
mask-edge features on constant image, invalid-padding invariance, and one-to-one
repeatability/tolerance sanity checks. Next: user reviews candidate overlays;
decide whether to keep baseline or compare contour detector, then descriptors
and matching. Do not mistake evaluator proximity associations for solved task.

## Candidate peak-selection correction (2026-10-18; current output)
User noticed circles over white areas. Audit of old candidates: no centres
outside supplied provisional support/tissue masks, but 154/180 MRI and 99/180
transformed-histology centres were within 2 pixels of mask edges; markers extend
into non-tissue and masks remain approximate. This does not prove anatomical
validity. Found actual algorithm issue: masking corner response before finding
peaks could manufacture local maxima on mask borders. Fixed to detect genuine
unmasked image-response maxima FIRST, then apply location eligibility filters.

Current output: results/candidates/harris_y8_peakfix/. Original run retained for
comparison; old counts/repeatability superseded. Added tiny pink centre dots
inside smaller cyan display rings. New counts: MRI 86, original histology 180,
transformed histology 180. Detector does not force an equal count.
15 tests pass, including new regression: restricting square image to a mask
excluding its true corners must not create substitute mask-edge peaks.
New histology repeatability within 1/2 mm =149/180 (82.8%)/154/180 (85.6%).
MRI/transformed histology proximity: 40/86 MRI candidates within 2 mm (46.5%),
which is 40/180 histology candidates (22.2%); not matching accuracy. Updated
random baselines and sensitivity are in repeatability.json/sensitivity.json.
Boundary points are not automatically false: anatomical cavity corners are
useful. User review of exact centres and possible overinclusive mask regions
still needed before accepting detector. Do not erode all boundaries blindly.

## Current handoff: first descriptor implemented (2026-10-18)
User accepted corrected detector provisionally despite some clustering and
authorised descriptor implementation. Built landmarks/descriptors.py and
scripts/test_descriptors.py; no new dependencies/GPU required. First baseline
is unsigned polar-gradient histogram: 4 radial rings x 12 angular sectors x
8 gradient orientations relative to radial direction =384 numbers per radius.
Square patches sampled at radii 16/24/36 working pixels, normalised to 41x41;
Gaussian prefilter, unsigned directions (global contrast reversal tolerance),
angular smoothing, square-root/L2 normalisation. Compare max cosine over all
3x3 radius combinations and 12 circular sector shifts (30-degree resolution).
Parameters do not use planted rotation/scale. Neither continuous scale nor
arbitrary contrast/deformation invariance is guaranteed.

Run `python scripts/test_descriptors.py --output results/descriptors/polar_gradient_y8_repeat`.
Default results/descriptors/polar_gradient_y8 already exists. Saves independent
candidate descriptors/points/usable-patch flags for MRI/original/transformed
histology BEFORE reading truth geometry. Original histology is a diagnostic
image; not a counterpart input to transformed-histology descriptor. Validity
masks exclude unsupported patch samples. Grayscale image intensities untouched.

Validation is DESCRIPTOR-ONLY, not final candidate matching: source candidates
are projected by evaluation geometry to target locations to construct reference
positive patches, described independently without orientation/scale from truth.
Gallery consists of these known-location target patches. Negatives <=5 mm from
query reference position are excluded; nearby ambiguity therefore untested.
Comparisons search radii/rotation equally for positives/negatives. True partner
ranks are diagnostic retrieval, not accepted anatomical matches or end-to-end
precision. No RANSAC, ratio-test matching or final registration output exists.

Results on one slice/transform:
- Original -> transformed histology: 170/180 reference partners rank first
  (94.4%, uniform random top1 chance ~0.56%, median rank1).
- MRI -> original histology: 21/86 rank first (24.4%, chance ~1.18%, median9.5).
- MRI -> transformed histology: 9/86 rank first (10.5%, chance ~1.18%, median15.5).
Strong within-modality geometric performance but poor cross-modality retrieval;
this descriptor is a baseline, NOT accepted as a reliable matching solution.
Approximate biological alignment also affects reference-positive patch identity.

Outputs: *_descriptors.npz; *_evaluation.npz (evaluation-only known locations and
score matrices); report.json; patch_examples.png; similarities.png. Example
query is fixed median CSV index, not selected best example. Patch plot shows
radius24 only, while scores search full bank; label explains this. No accepted
candidate pair list is exported. Source inputs SHA unchanged.
17 tests pass. Descriptor tests cover identity/global contrast reversal/90-degree
rotation behaviour, unsupported patch rejection, zero-information patches, and
input preservation. Next: user reviews descriptor results; discuss improving
cross-modality representation (e.g. self-similarity/MIND-style or learned feature
comparison) before proceeding to final matching. No automatic next-stage work.
Latest descriptor changes are in Mac workspace; ensure workstation transfer
includes these additions if ZIP was created earlier.

## Interpretable descriptor inspection (2026-10-18)
User could not interpret patch-only figure and requires future diagnostics to
include context. Added scripts/explain_descriptor.py and recorded requirement
in AGENTS.md. Run `python scripts/explain_descriptor.py` to regenerate views
without changing detection, descriptor or metrics. Current context outputs:
results/descriptors/polar_gradient_y8/context_histology_original.png and
context_histology_transformed.png, with context_examples.json.

Top panels show full MRI and histology with coloured boxes. Blue A1/B identify
MRI/reference crops; orange A2/C show the same MRI centre and the different
histology location preferred by descriptor. A1/A2 sizes can differ because
similarity searches radii. The actual best radius combinations are now shown,
rather than misleading fixed-radius examples with full-bank scores. Explicitly
label B as alignment-based expected location, not proven anatomical truth, and
C as different location, not accepted match. Native image orientations retained;
angular-search scores not probabilities. Best scores checked against saved
score matrices to 1e-6. Visualisation-only update; 17 prior tests still apply.
User review of baseline failures remains current stage; no new matcher added.

## Day 1 milestone completed: reproducible random pairs (2026-10-18)
User explicitly prioritised finishing README Day 1 before any other work.
No GitHub creation/push occurred; that task remains deferred. No new matching.

Generalised scripts/transform_pair.py with --angle, --scale, --shift and reusable
generate_pair(); fixed-example defaults remain 20 degrees /1.15/(40,25).
Added sample_parameters(seed,count) and scripts/generate_random_pairs.py.
Generated 5 examples at seed42 from accepted y=8 source plane. Random draws:
rotation uniform[-30,+30] degrees (positive clockwise); scale uniform[.85,1.2];
independent integer col/row shifts inclusive[-50,50] pixels. Starting experiment
ranges, not established realistic limits. Seed and realised values saved with
NumPy version. Image reduction now uses Gaussian antialias prefilter; filtering
support touching invalid data is excluded conservatively. Enlargements retain
previous one-bilinear-warp behaviour. Coordinates retain exact invertible mapping
regardless of intensity resampling; numerical invertibility is not lossless image
reconstruction. MRI unchanged; canvas containment checks prevent clipping.

Entry point from repository root:
`python scripts/generate_random_pairs.py --seed 42 --count 5`
Default results/random_pairs/y8_seed42 already exists; script refuses overwrite.
To reproduce: add `--output results/random_pairs/y8_seed42_repeat`.
To change random draws: change --seed and choose a new output directory.
Results include overview.png and pair_001 through pair_005, each with inputs.npz,
preview.png and truth/ (forward/inverse geometry, source hashes, original histology,
synthetic coordinate probes). Batch truth/manifest.json saves sampled values,
checks, seed and completion status; no parameters supplied to future matcher.
Overview panels fit independent canvas sizes, so do not interpret panel positions
as equal physical display scale. Each pair's preview uses a common canvas for
MRI/original/transformed histology and makes translation visible.

Realised pairs: angle/scale/(col,row):
001 +16.4374 /1.003607/(-7,36)
002 +11.8421 /0.882962/(3,48)
003 +15.6684 /1.125123/(1,-38)
004 -2.97684/0.979779/(-32,43)
005 +8.63191/1.137967/(5,-6).
19 tests pass, including same-seed reproducibility, different-seed difference,
random range checks, random geometry containment/inversion and reduction-filter
sanity/input preservation. Every real pair passes full-grid point roundtrip,
world mapping, footprint containment, source SHA preservation, MRI equality and
saved-array reload checks. Maximum full-grid pixel roundtrip across batch
3.41e-13 px; world error <8.7e-14 mm. Independent seed replay regenerated all
5 images and support masks exactly (truth/replay_check.json).

README Day 1: volumes/shared plane DONE; random unregistered pairs with retained
geometry DONE; independent initial detectors DONE on original/fixed transformed
example. Detectors have NOT yet been evaluated across all 5 random pairs. Five
augmentations are not five independent brains or a held-out test set. Flips,
crops and out-of-plane tilt remain later extensions. Next: user reviews batch
preview and confirms Day 1 status; then Day 2 baseline matching/RANSAC/metrics.

## Current handoff: end-to-end matching baseline (2026-10-18)
User authorised connecting current components to matching/RANSAC/evaluation before
feature improvements. Implemented landmarks/matching.py and
scripts/run_matching_baseline.py. No new learned model or GPU dependency.
Run from repo root: `python scripts/run_matching_baseline.py --output results/matching/polar_baseline_repeat`.
Default results/matching/polar_baseline already exists; refuse overwrite.

Actual inference uses images + validity masks only. Independently recomputes
MRI Otsu/histology 1% hysteresis, corrected Harris candidates and polar-gradient
descriptors. Compare actual independent candidate points, not known-location
oracle patches. Mutual nearest neighbours, source distance ratio<.9 and score>=.75.
Distance =sqrt(2-2*max-bank-cosine). Similarity RANSAC in MRI->histology pixel
coordinates, 3 noncollinear samples, residual<6 target working pixels, 5000
trials, rng0. Require >=4 noncollinear consensus pairs to return a model.
Thresholds fixed before evaluation, not selected using planted geometry.
No expectation of true rotation/scale/translation supplied to inference.

Predictions saved BEFORE reading pair truth: predictions.npz and prediction.json.
Then reference-only evaluation writes evaluation/report.json, all_proposals.csv,
matched-pair visualisation matches.png, and README-format landmarks.csv for final
retained pairs (header only if none). CSV score is similarity, not probability.
Batch summary.csv/summary.json cover fixed demo and all five random pairs.
Display rows: mutual proposals, distinctiveness filter, final RANSAC pairs.
Colours use reference error only: <=2 mm green, 2–5 amber, >5 red. At most 12
highest descriptor scores shown per row; all metrics include ALL pairs.

Results (mutual / distinctive / final retained / final <=2 mm):
fixed demo:29/11/0/0; pair001:27/5/0/0; pair002:35/10/0/0;
pair003:30/9/0/0; pair004:35/7/4/4; pair005:33/10/0/0.
Five of six cases lack >=4 geometrically consistent proposals. One random case
(pair004, mild ~-3 degree rotation) retains four pairs, all <=2 mm. Median pair
reference error .9684 mm; rotation error .6972 degrees; scale error .004605;
whole-MRI-candidate median transform error .7320 mm. Not a robust success claim.
No ratios are reported as 100% when zero pairs; precision/error are null.

Evaluation includes <=2/5 mm counts/fractions, median error, convex-hull area and
minor/major spread of <=2mm pairs, transform recovery, 200 same-count random
candidate pairings. This chance comparator does NOT condition on RANSAC selection;
not a p-value or geometric-fit null model. All are reference-coordinate metrics
with approximate anatomical alignment; six augmentations share ONE slice.

Positive control uses same pipeline descriptors/matching/RANSAC on independently
detected original vs transformed histology: 137 mutual, 121 distinctive, 117
retained, all 117 <=2 mm. Saves histology_control.json; generated by batch script.
This and synthetic tests support functioning matching geometry; MRI/histology
representation/proposal quality remains the practical bottleneck.
22 tests pass (new tests: mutual/ratio/invalid cases, known similarity recovery
with outliers, rejection of small/degenerate sets). Batch input hashes unchanged.
Assistant inspected fixed failure and pair004 success panels. Review with user
before changing representation or thresholds. End-to-end baseline now exists,
but project goal of reliable MRI–histology landmark correspondence is NOT met.
Next: interpret failures, compare one stronger representation/context strategy
using the same baseline evaluator; do not silently tune thresholds on this batch.
