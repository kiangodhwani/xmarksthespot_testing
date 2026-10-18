# Landmark outputs and success metrics

Start with **final_pipeline_summary.csv** and **final_pipeline/**: 19
image pairs tested with the frozen MIND top-3 + relational pipeline. Each pair
has one landmarks.csv in the exact requested eight-column MRI/histology schema.
MRI coordinates are RAS world mm; histology coordinates are (column,row) pixel
CENTRES in the saved transformed section, including its canvas offset. Names
are candidate IDs, not verified anatomical labels. Score is max-bank cosine,
not a probability. All retained predictions are exported; reference errors
never select CSV rows.

Final pipeline: 204/207 retained pairs within 2 mm;
207/207 within 5 mm. This pools repeated transformations
of three slices of ONE specimen. See per-case errors, count and spread; pooled
precision alone does not establish registration quality/generalisation.

- all_tests_summary.csv: all 100 saved prediction artifacts, including
  baselines, ablations, repeated runs, final tests, failures and controls.
- historical_summary.csv / historical/: earlier MRI/histology variants.
- controls_summary.csv / controls/: histology-only diagnostics, with source
  columns named histology_original_x_mm/y_mm/z_mm. They are intentionally NOT
  falsely labelled as MRI coordinates or mixed into registration outputs.
- prediction_duplicate_of flags exact repeated prediction semantics (including
  candidate/score/proposal/model arrays). Do not count repeated variants as new
  anatomical tests. An empty no-consensus run still has a header-only CSV;
  accuracy fractions and transform errors are blank, not 100% or zero error.
- Each pair's evaluation/ folder contains metrics.json, original_report.json,
  evaluation-only reference_geometry.json and per-row landmark_errors.csv.
- provenance.json maps CSV row IDs to original candidate indices and files.
- manifest.json inventories every artifact, aliases, exclusions and source hashes.

Metrics: retained count, fractions <=2/5 mm, median/max error, correct-pair convex
hull area (mm²), minor/major spread and >=4 noncollinear correct-pair flag;
planted/estimated rotation and scale and their absolute errors; candidate-field
median transform recovery error. Saved chance metrics are preserved; missing
legacy-control chance metrics are computed from frozen candidate pools using
200 same-count unique random pairings and seed7. Chance is NOT conditioned on
consensus selection or a p-value. observed_exceeds_random_mean is a descriptive
comparison, not a significance test. Convex-hull/mask area ratios in original
reports are spread ratios and can exceed1; they are not tissue coverage.

Original polar-baseline control has no saved coordinate artifact. Its report
exactly matches a saved backend replay; manifest records that replay as an
explicit alias, rather than claiming recovered original coordinates. Preliminary
candidate/descriptor oracle diagnostics are not final match predictions. The
-30° boundary preparation was never matched and therefore has no invented CSV.

No matching, fitting, threshold changes or source modifications were performed.
World coordinates and success metrics were independently checked against saved
predictions/reference mappings and existing CSVs/reports. Same source specimen
and approximate released alignment remain limitations. Derived BigBrain data:
Amunts et al., Science (2013), CC BY-NC-SA 4.0.

Recreate into a new output directory:
`.venv/bin/python scripts/export_all_landmarks.py --output results/exports/landmark_release_v2`
