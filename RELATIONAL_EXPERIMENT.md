# MIND proposals and relational geometry

Implemented after the user reviewed ratio-filter ablation and authorised keeping
multiple partners plus relational selection. No descriptor or broader-context
change is included: this isolates proposal retention and geometric selection.
SuperGlue remains a possible later comparison, not installed or run here.
Results await user review.

## Fixed setup

Reuse `results/matching/mind_v1` image-only candidates/similarities without
recomputing descriptors, changing prepared images, or adding reference locations.
Keep cosine >=0.75, residual <6 target working pixels, minimum four supporters,
and the same noncollinearity check. Evaluation thresholds remain <=2/5 mm.
There is no ratio test in any of these comparisons.

- **mutual_ransac:** the previous no-ratio pipeline, seed 0 and 5000 trials.
- **mutual_relational:** same mutual proposals, new deterministic geometry solver.
- **top3_relational:** up to three source-ranked partners per MRI candidate,
  score >=0.75, no mutual-nearest restriction; same new geometry solver.
- **mutual_exhaustive:** FIXED CASE ONLY diagnostic, enumerate all 3654 mutual
  three-pair subsets using the original noncollinearity/consensus/refit checks;
  select count then squared residual. No reference geometry selects the model.

The solver control separates some effects of geometry from expanding proposals.
Top3 also removes the mutual restriction, so this does not separately attribute
its benefit to k=3 versus allowing non-mutual alternatives. Top1/3/5 source
retrieval is measured separately for diagnosis; k=5 is NOT fitted or tuned.

All 22 predictions (three methods x seven inputs, plus fixed-case diagnostic)
were saved before ANY reference geometry was read. The reference is used only
for evaluation/visualization after inference. The seeded mutual control matches
all seven previous no-ratio models and retained arrays exactly. Original MIND
prediction files pass unchanged SHA-256 checks.

## Relational solver

Two proposed correspondences with distinct source and target IDs define an
orientation-preserving similarity hypothesis. Enumerate every such pair; exact tested and valid counts are recorded per case. For each hypothesis,
find proposals with residual <6 target pixels. Solve maximum-cardinality
one-to-one assignment, then minimize squared residual. Dummy assignments permit
unmatched sources. Alternatives from one source cannot inflate consensus;
one target can support only one source.

Hypotheses are ranked by initial unique-source/target support upper bounds.
Skip bounds below the best final count, then perform at most three assignment /
least-squares refit iterations. Final hypotheses require >=4 noncollinear pairs
and every returned pair is within the same residual bound. Select by final
count, then squared residual, then summed cosine. No scale/angle range prior,
canvas-centering cue, anatomical mask agreement or reference transform enters.

**Solver limitation:** support bounds hold before refitting; pruning can skip a
hypothesis that would grow after refitting. This is a deterministic consensus
heuristic, not a globally optimal graph-matching proof. Two-point hypotheses also
behave differently from three-point RANSAC: the mutual solver control loses one
pair in case001. RANSAC itself and the baseline default are unchanged; the new
solver is an explicit experimental alternative.

## Results

Final pairs within 2 mm / total retained:

| Case | Mutual + RANSAC | Mutual + relational | Top3 + relational |
|---|---:|---:|---:|
| Fixed | 0/0 | 0/0 | 8/8 |
| 001 | 5/6 | 4/5 | 10/11 |
| 002 | 5/5 | 5/5 | 11/11 |
| 003 | 0/0 | 0/0 | 9/9 |
| 004 | 11/11 | 11/11 | 12/14 |
| 005 | 3/4 | 3/4 | 10/10 |
| Histology control | 96/96 | 96/96 | 137/137 |

Top3 relational returns models in all six cases: 60/63 retained within 2 mm,
all 63 within 5 mm. These summed counts are descriptive repeated-slice results,
not independent-case statistical evidence. The final three 2-5 mm pairs are one
in 001 and two in 004. No retained pair exceeds 5 mm in this run, but additional
proposals can create false consensus on other inputs. All precision, errors,
coverage, hull fraction/spread and chance comparator are saved per method/case.

Whole-MRI-candidate median transform errors for top3:
fixed .9186 mm; 001 .6276 mm; 002 .7929 mm; 003 .8339 mm;
004 1.2010 mm; 005 .7580 mm. Case004 worsens from .8564 mm with mutual RANSAC,
despite gaining one <=2 mm pair; increased consensus does not guarantee better
transform recovery. Median retained pair errors range .8414-1.2736 mm.
Control retains 137 correct pairs vs 96 in no-ratio MIND and 117 in original polar;
this is a same-modality synthetic control, not cross-modality generalisation.

### What happened in the fixed case?

Exhaustive mutual triples test 3654 subsets, of which 3411 are noncollinear.
A four-pair model exists that seed 0 RANSAC missed: three pairs <=2 mm, all four
<=5 mm, median whole-candidate recovery error 2.1067 mm. So random sampling did
contribute to the earlier no-model return; this is a weak four-pair solution.
With expanded top3 proposals, the relational solver instead finds 8/8 <=2 mm and
.9186 mm recovery error. Sampling was not the only bottleneck.

### Reference-partner rank diagnostic

Each case has 86 MRI candidates. Count a source once if at least one target in
its first k score-eligible choices is within 2 mm (source retrieval, not mutual
matching). No expected-position patches are inserted into the candidate gallery.

| Case | Top1 | Top3 | Top5 diagnostic only | Any detected target within 2 mm |
|---|---:|---:|---:|---:|
| Fixed | 7 | 12 | 16 | 41 |
| 001 | 8 | 15 | 16 | 44 |
| 002 | 7 | 12 | 15 | 36 |
| 003 | 5 | 11 | 14 | 42 |
| 004 | 14 | 17 | 19 | 40 |
| 005 | 9 | 14 | 19 | 40 |

This supports retaining alternatives rather than treating local top1 ambiguity
as rejection. Candidate coverage and descriptor ranking still limit how many
landmarks can be recovered. A nearby reference candidate is not proven anatomy.

## Reproduction, validation and figures

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/run_relational_experiment.py --output results/matching/relational_v1_repeat
.venv/bin/python scripts/review_relational_experiment.py --input results/matching/relational_v1_repeat
```

The experiment requires saved MIND predictions and its histology-control
predictions (override with --input). Use a new output directory; existing runs
are preserved. No GPU/new dependencies. Installed scikit-image 0.26 supports the
new nondeprecated `SimilarityTransform.from_estimate` call.

32 tests pass, including stable top-k/invalid handling, synthetic recovery from
incorrect first-choice alternatives, one-to-one maximum-cardinality assignment
with unmatched sources, exhaustive triple recovery and degenerate rejection.
No geometry/preparation changes. Input source hashes and seeded control exact
reproduction are checked at runtime. Scores are unchanged and already encode
the MIND-style 3x3 radius / 12-rotation search documented in MIND_EXPERIMENT.md.

Outputs under `results/matching/relational_v1`:

- `experiment_config.json`, `checks.json`, `summary.csv`, `reports.json`.
- `retrieval.csv`, per-case `reference_partner_ranks.csv` and
  `retrieval_diagnostic.json`: evaluation-only top1/3/5 diagnostics.
- Per case/method `predictions.npz`, `prediction.json`, separate `evaluation.json`.
- `comparison.png`: all six cases and the solver control.
- `retrieval.png`: rank counts with explicit candidate-availability denominator.
- Per MRI/histology case `matches.png`: full-image context. **If a method fails,
  show its input proposals with incorrect lines; if it succeeds, show final
  retained pairs.** Labels are local to each row; colours are evaluation-only
  green<=2, amber2-5, red>5 mm. Show up to 12 highest-cosine pairs per row; metrics
  include every pair. Fixed-case plot adds the exhaustive diagnostic fourth row.

Assistant inspected comparison, retrieval and fixed/003 full-context figures.
Figures label input vs final and reference-only colours. No new default,
held-out claims, semantic anatomical context or SuperGlue stage has been added.
Review this stage before deciding whether to freeze and validate on new slices.
