# Fixed 2D MIND-style experiment

Implemented and evaluated on the existing y=8 mm benchmark. User reviewed the
regenerated shared plane and transformations, then the reproduced polar baseline,
and authorised this experiment. Results below await user review.

This is a sparse 2D adaptation of the self-similarity idea in Heinrich et al.
(2012), https://doi.org/10.1016/j.media.2012.05.008. It is not a reproduction of
the original dense 3D MIND or MIND-SSC. A negative result applies to this
configuration and pipeline, not every use of MIND.

## Fixed method

`landmarks/mind.py` computes independent descriptors from an image, validity mask
and candidate coordinates. Masks and Harris detection remain unchanged. The
configuration was written before real-data inference; there was no parameter
sweep or post-evaluation tuning.

- Independent valid-image 1st/99th percentile normalisation, clipped to [0,1].
- Radii 16, 24, 36 working pixels; each sampled to a 41x41 context patch. Outer
  sample-centre spans are 32, 48, 72 working pixels. At the original MRI grid
  these are 6.4, 9.6, 14.4 mm; inference never uses transformed reference spacing.
- Gaussian prefilter sigma=max(0.7, radius/20*0.7), truncate=4. Validity requires
  all prefilter and bilinear sampling support, not just eligible patch centres.
- 25 context centres on a 5x5 grid at canonical offsets -12,-6,0,6,12.
- Eight comparison offsets: cardinal and diagonal displacements with components
  +/-4 canonical samples. Local SSD uses a 5x5 Gaussian kernel, sigma=1,
  truncate=2. Its outer-centre span is 3.2, 4.8 or 7.2 working pixels depending
  on context radius, in addition to the prefilter support.
- Local self-similarity exp(-SSD/V), V=mean four cardinal SSDs, with numerical
  floor max(1e-10, 0.001*mean V). Divide each location's channels by their maximum;
  concatenate 25x8 values and L2 normalise. Reject constant/no-information patches.
- Max cosine over all 3x3 source/target radius pairs and 12 target patch sampling
  rotations (0 through 330 degrees in 30 degree steps). Source orientation is
  fixed to zero. All comparisons use the same search, without reference angles
  or scales. Target sampling angle is not a recovered whole-image transform.

Mutual nearest neighbours, ratio<0.9, cosine>=0.75, distance=sqrt(2-2*cosine),
RANSAC residual<6 target pixels, minimum consensus four, seeds and 2/5 mm
reference thresholds are unchanged. Similarity distributions differ between
representations; unchanged score thresholds are not calibrated equivalently.
Full rotation/scale search can also increase similarity to unrelated patches.

## Results

Counts within 2 mm / total pairs, using ALL pairs at each stage:

| Case | Polar mutual | MIND mutual | Polar distinctive | MIND distinctive | Polar final | MIND final |
|---|---:|---:|---:|---:|---:|---:|
| Fixed | 2/29 | 4/29 | 1/11 | 1/9 | 0/0 | 0/0 |
| 001 | 2/27 | 6/30 | 0/5 | 1/6 | 0/0 | 0/0 |
| 002 | 3/35 | 5/36 | 1/10 | 1/4 | 0/0 | 0/0 |
| 003 | 1/30 | 3/30 | 0/9 | 1/7 | 0/0 | 0/0 |
| 004 | 7/35 | 11/32 | 4/7 | 5/11 | 4/4 | 4/5 |
| 005 | 5/33 | 5/29 | 3/10 | 0/4 | 0/0 | 0/0 |

Both methods return no model in five cases. MIND improves initial correct-pair
counts in five cases, but this does not improve final <=2 mm counts. Pair004's
fifth retained MIND pair is between 2 and 5 mm. Pair004 median retained error:
polar 0.9684 mm, MIND 1.2728 mm. Median whole-candidate transform recovery error:
polar 0.7320 mm, MIND 1.6524 mm. Rotation absolute error: 0.6972 vs 0.7051 deg;
scale absolute error: 0.004605 vs 0.034312. The <=2 mm pair hull area increases
from 1025.76 to 1647.42 mm2 (12.48% to 20.04% of MRI mask area), with spread
minor/major ratio 0.197 vs 0.450; better coverage does not remove accuracy limits.

Histology control: polar retains 117, MIND 65; all retained control pairs are
within 2 mm. The implementation retains useful same-modality geometric behavior
but loses matches relative to polar. MIND inference took 3.39-4.10 s/case;
polar backend-check inference 0.62-1.05 s/case on this machine. Runs overlapped,
so these are descriptive measurements, not a controlled speed benchmark.

The fixed median usable MRI query diagnostic prefers an unrelated fixed-example
candidate (33.47 mm reference error, cosine 0.9407) over the reference-location
patch (cosine 0.9086). This is an illustration, not a population score analysis.
Its preferred candidate is rejected by filtering. Reference-location patches
are evaluation-only diagnostics and never inserted into the matching gallery.

This configuration is not an overall improvement. Do not change thresholds or
start another representation experiment until these results are reviewed.
Six augmentations share one slice/specimen. Alignment remains approximate.
Random-pairing comparison is not conditioned on RANSAC selection or a p-value.

## Reproduction and checks

From the repository root, after accepted preparation exists:

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/run_matching_baseline.py --descriptor mind --output results/matching/mind_v1_repeat
.venv/bin/python scripts/run_matching_baseline.py --output results/matching/polar_backend_check_repeat
.venv/bin/python scripts/review_mind_experiment.py --mind results/matching/mind_v1_repeat --polar results/matching/polar_backend_check_repeat
```

Inference refuses existing output directories. Review diagnostics regenerate
inside the supplied MIND output directory. No GPU or new dependency is needed.
Exact installed versions are saved in
`results/matching/mind_v1/requirements-resolved.txt`.

26 tests pass, including identity/contrast reversal/90-degree rotation,
synthetic 1.5x scale, invalid/prefilter-support exclusion, invalid-padding
independence, constant/empty inputs, and score-choice reconstruction.
All original polar prediction arrays were checked for exact equality after the
backend addition. All six MRI/histology candidate arrays are identical between
polar and MIND. Usable candidate counts are also identical (MRI86; histology180,
except pair002 with 177 usable).

Predictions and selected radius/angle indices are saved before reference geometry
is loaded. The histology-control predictions are now also saved before its
reference geometry is loaded. Its original histology image is explicitly the
control input, not a source for cross-modality inference. Evaluation functions
and benchmark thresholds are unchanged. Runtime input hash checks pass.

Outputs in `results/matching/mind_v1/` include configuration, per-case predictions,
original evaluator reports (2/5 mm accuracy, coverage, recovery and chance),
comparison.csv/json/png, histology_control.json and its predictions, and six
`evaluation/descriptor_context.png/json` diagnostics. Context figures use fixed
median usable query index 43, explicitly label reference vs algorithm preference,
and display the actual smoothed scoring patches and sampling angles/radii.
Every displayed score was reconstructed from its two displayed patches to 1e-6.
The comparison, all six context figures, and pair004 match panel were visually
inspected; overlapping titles were corrected before presentation.
