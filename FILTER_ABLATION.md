# Proposal-filter ablation

User reviewed the MIND result and observed that the distinctiveness filter was
removing useful matches before RANSAC. This fixed ablation tests that hypothesis;
it does not replace the default baseline. Results await user review.

## Method and isolation

Reuse the saved, image-only candidates and similarity matrices from `mind_v1`
and `polar_backend_check`, plus both saved histology controls. No descriptors,
images, masks, candidates or reference-derived locations are changed/recomputed
for match inference. Fit three variants:

1. Current: mutual nearest neighbours, cosine >=0.75, source distance ratio <0.9.
2. Similarity only: mutual nearest neighbours and cosine >=0.75, no ratio test.
3. All mutual: mutual nearest neighbours only, no score or ratio cutoff.

All variants use the SAME existing `geometric_filter`: similarity model,
three noncollinear samples, residual <6 target pixels, minimum four supporters,
5000 trials, seed 0, stop probability .999. There is no RANSAC tuning/seed search.
The original `propose` and `geometric_filter` implementations remain unchanged.

All 42 variant predictions (two representations x seven inputs x three variants)
were frozen before opening ANY reference geometry for evaluation. The evaluator
uses the same existing `evaluate_pairs`, 2/5 mm thresholds, coverage measures,
and 200 random candidate pairings. Transform recovery follows the baseline's
formulas. Random pairing is not a RANSAC-conditioned null or p-value.
The control's source image is original histology, not MRI; its mask area is
computed independently from that original histology for coverage normalization.

No-ratio variants have identical proposal and prediction arrays in all 14 input
sets because every mutual proposal already scores >=0.75. This isolates the ratio
condition: the absolute cosine cutoff removes no mutual proposals in this run.

## Main results

Entries are final pairs within 2 mm / total retained. Both no-ratio variants
produce the same results, so they share a column.

| Case | MIND current | MIND no ratio | Polar current | Polar no ratio |
|---|---:|---:|---:|---:|
| Fixed | 0/0 | 0/0 | 0/0 | 0/0 |
| 001 | 0/0 | 5/6 | 0/0 | 0/0 |
| 002 | 0/0 | 5/5 | 0/0 | 0/0 |
| 003 | 0/0 | 0/0 | 0/0 | 0/0 |
| 004 | 4/5 | 11/11 | 4/4 | 6/6 |
| 005 | 0/0 | 3/4 | 0/0 | 4/4 |
| Histology control | 65/65 | 96/96 | 117/117 | 127/127 |

MIND consensus increases from one to four of six cases. Polar increases from
one to two. All final retained pairs, including the controls, are within 5 mm
of the reference in this run. There is no observed >5 mm retained-pair failure;
this is not evidence that false consensus cannot occur on other inputs/seeds.
The report's `model_with_majority_pairs_beyond_5mm` flag is just a coarse
EVALUATION diagnostic (more than half final pairs beyond 5 mm), not an inference
acceptance criterion. It is false throughout this run. Use all errors and
transform recovery to assess models rather than relying on this flag alone.

For MIND pair004, the ratio test removes exactly SIX of ELEVEN <=2 mm mutual
proposals. None of those six fail the similarity cutoff. Feeding all mutual
proposals lets unchanged RANSAC retain all eleven <=2 mm pairs, while excluding
the previous 2-5 mm pair. Whole-MRI-candidate median transform error improves
from 1.6524 to 0.8564 mm. Matched-point coverage and spread are in each report.

MIND no-ratio whole-candidate median transform errors: pair001 2.1260 mm,
pair002 0.8569 mm, pair004 0.8564 mm, pair005 1.3426 mm. In pair001, five of six
retained pairs lie within 2 mm, but this does not establish equally accurate
recovery throughout the candidate field. Pair001 and pair005 retain one
2-5 mm pair each. Fixed and pair003 still return no model.

Interpretation: the current ratio heuristic is a material bottleneck for this
MIND-style representation. The previous finding that MIND did not improve final
matching applies to the ORIGINAL filter; this ablation shows the representation
can yield useful additional geometric evidence when ratio filtering is removed.
This is an exploratory result on the same reviewed single-slice batch, not a
held-out improvement claim. Neither the default filter nor benchmark thresholds
have been silently changed. Review before choosing a new default or extending
validation. Six augmentations do not constitute six independent brains/slices.

## Artifacts and reproduction

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/run_filter_ablation.py --output results/matching/filter_ablation_v1_repeat
.venv/bin/python scripts/review_filter_ablation.py --input results/matching/filter_ablation_v1_repeat
```

Prerequisites are both saved baseline runs, including their
`histology_control_predictions.npz` files; see MIND_EXPERIMENT.md. Inputs can be
overridden with --mind and --polar. The runner refuses an existing output path;
the review script regenerates only figures inside the supplied output folder.
No new dependency/GPU needed. No data preparation or geometry changes.

`results/matching/filter_ablation_v1/` contains:

- `experiment_config.json`: fixed named variants and unchanged settings.
- `summary.csv`, `reports.json`: all cases, controls, 2/5 mm accuracy,
  coverage, transform recovery, chance comparison and RANSAC time.
- `checks.json`: source hashes and inference/evaluation separation assertions.
- `{mind,polar}/{case}/{variant}/predictions.npz`, `prediction.json`,
  `evaluation.json`: frozen predictions and separate reference evaluation.
- Per-case `proposal_audit.json`: image-only ratio/score decisions; separate
  `proposal_audit_evaluation.csv` adds reference errors after prediction.
- `comparison.png`: all six cases, both representations, all three variants.
- `pair004_filter_audit.png`: all 32 MIND mutual proposals, score versus ratio,
  coloured by evaluation error. Cosine is not a probability.
- Per MRI/histology case `final_comparison.png`: full-image context for final
  RANSAC pairs under all three variants. Stable mutual-proposal IDs across rows,
  explicit evaluation-only colours, display cap 12 highest scores, ALL metrics
  include all pairs. Existing MIND context figures document actual scoring
  patches and the fixed 3x3 radii/12 rotation search; no patch settings changed.

28 tests pass: prior 26 plus independent condition-removal and tie/invalid/empty
cases. The current variant reproduces all 14 saved models AND retained-pair
arrays exactly. No-ratio variants are array-identical in all 14 sets. All input
prediction hashes are unchanged. Comparison, pair004 audit and case001/004
full-image figures were inspected; crowded audit labels were removed.
