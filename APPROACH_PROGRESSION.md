# From the first baseline to MIND with relational geometry

This document follows the evidence that led to our current MRI–histology landmark matcher. It links existing figures; no new images or experiments were created for this account.

The working solution combines **independent multiscale Harris candidates, multiscale MIND-style structural descriptors, multiple possible partners, and one-to-one geometric consensus**. Its success is a controlled proof of concept for known corresponding sections, rather than a solution to locating an unknown histology plane in a 3D MRI volume.

## 1. Start with a trustworthy experiment

We first established the coordinate conventions, extracted MRI and histology onto exactly the same coronal plane, and applied known rotation, scale and translation to histology. We retained the complete transform, including expanded-canvas offsets, for later evaluation. Prepared images were reviewed before matching.

We used image-derived search support independently in each modality. Early preparation showed that a simple histology Otsu threshold excluded pale interior tissue; background-relative hysteresis provided provisional search support. These masks restrict candidate locations; they are not anatomical labels or validated tissue segmentations. Grayscale structure remains available to the detector.

**Existing preparation figures:**

- [Original shared-plane pair](results/pairs/coronal_y8_repeat/pair.png): compare the underlying anatomy and contrast before synthetic distortion.
- [Fixed transformed histology preview](results/transforms/coronal_y8_demo/preview.png): see the planted rotation, enlargement and translation. Preview probes are known coordinate examples, not detected anatomical matches.
- [Five random transformations](results/random_pairs/y8_seed42/overview.png): the initial benchmark includes these five cases plus the fixed example, all from the same y=8 mm slice.

We kept reference geometry separate from inference. Predictions were saved before reference transforms were loaded for scoring. Reference transforms determine evaluation errors and line colours; they do not tell the matcher where to look. See [FOUNDATIONS.md](FOUNDATIONS.md) for the coordinate and evaluation contract.

## 2. First baseline: does ordinary local matching transfer across modalities?

We detected candidates separately using multiscale Harris responses, then described local structure with a polar unsigned-gradient descriptor. The first matcher used:

1. Mutual nearest-neighbour descriptor proposals.
2. Cosine similarity at least 0.75 and a best/second-best distance ratio below 0.9.
3. Similarity-model RANSAC, allowing unmatched candidates.

The histology-to-transformed-histology control retained **117 matches, all within 2 mm**. MRI-to-histology was much weaker: **five of six cases had no consensus**, while pair004 retained **four matches, all within 2 mm**.

That contrast suggested that the cross-modality representation and proposals were major bottlenecks. The control supported useful same-modality operation of the pipeline; it did not prove every detector or geometric choice was optimal.

**Compare these existing figures:**

- [Polar baseline: successful pair004](results/matching/polar_baseline/pair_004/matches.png).
- [Polar baseline: failed fixed example](results/matching/polar_baseline/coronal_y8_demo/matches.png).

Each panel follows proposals → distinctiveness filtering → final RANSAC pairs. Green means reference error ≤2 mm, amber means >2 to ≤5 mm, and red means >5 mm. These are evaluation-only colours. An empty final row means insufficient consensus, not perfect accuracy.

## 3. Why choose MIND?

MRI and histology do not have a reliable one-to-one intensity relationship. Some boundaries and folds are shared, but their brightness, sharpness and apparent shape differ. We therefore wanted a representation based more on **internal structural relationships** than on matching absolute intensities.

MIND describes how a small neighbourhood resembles neighbouring patches within the *same image*. The cross-image comparison then asks whether those patterns of self-similarity agree. That made it a plausible next experiment after the weak gradient baseline.

It also fit the existing framework: a classical, inspectable method requiring no training data, GPU or additional learned-model dependency. **We did not establish that MIND was the best of every candidate method.** We selected one motivated structural representation to test cleanly. DINOv2, learned matchers and other descriptors were not exhaustively benchmarked; SuperGlue remained a possible fallback and was not run.

Our implementation is a **sparse 2D MIND-style adaptation**, not a reproduction of dense 3D MIND or MIND-SSC. Its fixed configuration uses radii 16, 24 and 36 working pixels, sampled as 41×41 context patches, with local 5×5 SSD comparisons. Matching searches all 3×3 radius combinations and 12 target rotations at 30° intervals. The maximum cosine score is a similarity measure, not a probability. No planted rotation or scale enters that search.

**Evidence and implementation:**

- [MIND versus polar comparison](results/matching/mind_v1/comparison.png).
- [MIND pair004: proposals through final matches](results/matching/mind_v1/pair_004/matches.png).
- [Actual scoring patches and full-image context](results/matching/mind_v1/pair_004/evaluation/descriptor_context.png): distinguishes reference-location diagnostics from the algorithm's preferred candidate and documents the sampled patches.
- [Exact experiment specification](MIND_EXPERIMENT.md), [descriptor implementation](landmarks/mind.py). The specification cites the original MIND paper.

MIND improved the number of reference-correct mutual proposals in five of six cases. For pair004, that count rose from **7 to 11**. However, under the original filter, MIND still returned a model in only **one of six cases**. Better initial evidence was not yet becoming better final performance.

## 4. The ratio test was discarding useful evidence

The user noticed that many green lines disappeared between the first and second rows of the MIND figure. We tested this directly rather than changing several components together.

The ratio test compares the best descriptor distance with the second-best distance. A ratio near one means two possible partners look similarly good. The original rule required a ratio below 0.9, so it discarded ambiguous candidates even when their best partner was reference-correct.

Repeated folds and bilateral structure can create such ambiguity. A locally ambiguous match can still become useful when other correspondences establish the global geometry.

We compared the original filter, score-only filtering, and all mutual proposals, with **the same saved descriptors, candidates and RANSAC settings**. For pair004:

- There were 11 mutual proposals within 2 mm.
- Five passed the ratio rule; six failed it.
- None of those six failed the cosine cutoff.
- Removing the ratio test let unchanged RANSAC retain **11/11 within 2 mm**, compared with **4/5** under the original filter.

Across all six MIND cases, model returns increased from **1/6 to 4/6**. Fixed and pair003 still failed.

**The clearest existing figures:**

- [Why the filter rejected correct proposals](results/matching/filter_ablation_v1/pair004_filter_audit.png): all 32 mutual proposals, with ratio on the horizontal axis and cosine similarity on the vertical axis. Correct green points to the right of the 0.9 boundary fail the ratio rule.
- [Pair004 final matches under all three filters](results/matching/filter_ablation_v1/mind/pair_004/final_comparison.png): see the jump from four to eleven ≤2 mm matches.
- [All-case filter comparison](results/matching/filter_ablation_v1/comparison.png).
- [Fixed case still failing after removing the ratio](results/matching/filter_ablation_v1/mind/coronal_y8_demo/final_comparison.png).

**“No ratio” means removing that ambiguity test only.** It preserves mutual nearest neighbours, the cosine cutoff and geometric filtering. Score-only and all-mutual variants happened to coincide because every mutual proposal already scored at least 0.75. See [FILTER_ABLATION.md](FILTER_ABLATION.md).

## 5. Why add alternative partners and relational geometry?

Removing the ratio could preserve an existing proposal, but it could not recover a correct partner that mutual nearest-neighbour selection never proposed. Some correct partners were among a candidate's next-best choices, rather than its first mutual choice.

We therefore retained **up to three score-eligible histology partners per MRI candidate**, removing the mutual-nearest restriction, and asked which alternatives fit together geometrically.

The new solver forms similarity-transform hypotheses from two proposed correspondences, then finds a **one-to-one assignment** consistent with each hypothesis. It permits unmatched candidates, refits, and selects consensus primarily by match count and then squared residual. A source or target cannot be counted repeatedly to inflate support.

We preserved the residual bound of **<6 target working pixels**, minimum **four noncollinear supporters**, and evaluation thresholds of **2 mm and 5 mm**. This is an explicit alternative to the original RANSAC solver, not simply RANSAC with another name. Its pruning/refitting procedure is a deterministic heuristic, not a proof of global optimality.

The experiment included a mutual-proposals/new-solver control. That control did not solve the remaining two failures: **expanded alternatives mattered**, not just swapping the solver. Conversely, top-three also removed mutuality, so this experiment does not isolate the benefit of k=3 from the benefit of allowing non-mutual partners. Multiscale context was already part of the MIND representation; this stage did not add a new descriptor.

| y=8 mm case | Original polar | MIND, original filter | MIND, no ratio + RANSAC | MIND, top-three + relational |
|---|---:|---:|---:|---:|
| Fixed | 0/0 | 0/0 | 0/0 | 8/8 |
| 001 | 0/0 | 0/0 | 5/6 | 10/11 |
| 002 | 0/0 | 0/0 | 5/5 | 11/11 |
| 003 | 0/0 | 0/0 | 0/0 | 9/9 |
| 004 | 4/4 | 4/5 | 11/11 | 12/14 |
| 005 | 0/0 | 0/0 | 3/4 | 10/10 |

Entries are **pairs within 2 mm / total retained**, not recall. The final method returned models in all six cases: **60/63 within 2 mm; 63/63 within 5 mm**.

**Existing figures:**

- [Reference-partner rank diagnostic](results/matching/relational_v1/retrieval.png): shows why retaining alternatives can help; reference locations are evaluation-only and never inserted into the candidate gallery.
- [All-case geometry comparison](results/matching/relational_v1/comparison.png).
- [Fixed case: failure becomes eight correct matches](results/matching/relational_v1/coronal_y8_demo/matches.png).
- [Pair003: the other remaining failure becomes nine correct matches](results/matching/relational_v1/pair_003/matches.png).
- [Pair004: additional matches with remaining error](results/matching/relational_v1/pair_004/matches.png).

An exhaustive diagnostic also found a weak four-pair mutual model in the fixed case that seeded RANSAC missed. Thus sampling contributed to that failure, but was not the only issue. Expanded proposals produced a stronger eight-pair result. More matches also do not guarantee better transform recovery: pair004's median candidate-field recovery error worsened from about 0.856 to 1.201 mm despite gaining matches. See [RELATIONAL_EXPERIMENT.md](RELATIONAL_EXPERIMENT.md) and [landmarks/relational.py](landmarks/relational.py).

## 6. Freeze the method and test additional anatomy

We then froze the method and tested two additional coronal slices, y=−20 and y=35 mm, using the same six transformation settings after preparation review. All **12 cases returned models**, with **133/133 retained cross-modality matches within 2 mm**.

- [Two-slice results overview](results/validation/two_slices_v1/matching/comparison.png).
- [y=−20 mm, fixed example](results/validation/two_slices_v1/matching/y_m20/fixed/matches.png).
- [y=35 mm, fixed example](results/validation/two_slices_v1/matching/y35/fixed/matches.png).
- [Validation protocol and metrics](TWO_SLICE_VALIDATION.md).

We also tested a predeclared targeted challenge: **−15° rotation, scale 0.85, and translation (−50,+50) pixels**. The angle lies halfway between descriptor rotation samples, and shrinkage is at the permitted minimum. The unchanged pipeline retained **11/11 within 2 mm**, with median candidate-field transform recovery error about **0.570 mm**.

- [Stress preparation](results/validation/rotation_mismatch_stress_v1/pair/preview.png).
- [Stress proposals and final matches](results/validation/rotation_mismatch_stress_v1/matching/stress/matches.png): incorrect input lines remain visible above the final result.
- [Stress protocol and interpretation](BOUNDARY_STRESS_TEST.md).

This passed a targeted challenge; it was not proven to be the hardest possible transform. The earlier −30° boundary preparation was not matched and is excluded from results.

## 7. Current outcome and its limits

Across the **19 final-pipeline cases**, we retained **207 pairs: 204 within 2 mm (98.55%) and all 207 within 5 mm**. Every case had at least four noncollinear reference-correct pairs within 2 mm. Count-weighted random-pairing means were approximately **0.35% within 2 mm** and **1.24% within 5 mm**, using the same candidate pools and output counts.

That clearly beats random pairing on these tests. The random comparator does not pass random proposals through equivalent geometric selection, so it is not a geometry-conditioned statistical significance test.

The evidence covers three slices of one specimen and synthetic similarity transforms. It does not establish performance on new subjects, damaged or cropped sections, flips, nonrigid distortion, or unknown out-of-plane sections. Released MRI–histology alignment is approximate; reference-close pairs are not independently verified anatomical identities. Retained-pair precision is not dense registration accuracy or detection recall.

**Reviewable outputs:** [final-case metrics](results/exports/landmark_release_v1/final_pipeline_summary.csv), [complete historical inventory](results/exports/landmark_release_v1/all_tests_summary.csv), and [CSV format/provenance guide](OUTPUTS.md). Standard per-pair CSVs contain MRI world-mm coordinates, transformed-histology pixel coordinates and descriptor similarity scores.

The progression was evidence-led: a strong same-modality control exposed the cross-modality gap; MIND supplied better structural proposals; the filter ablation showed useful evidence was being discarded; expanded alternatives and one-to-one geometry resolved the remaining benchmark failures; frozen additional-slice and stress tests checked whether that improvement extended beyond the development cases.

*Existing derived imagery uses BigBrain data: Amunts et al., Science (2013), CC BY-NC-SA 4.0; see the repository README for attribution and licensing.*
