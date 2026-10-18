# Project foundations
Before changing data preparation, matching, evaluation, or visualisation, read
PROJECT_PATH.md and FOUNDATIONS.md. Preserve the coordinate conventions and run the foundation
checks after any geometry change. Do not start landmark matching until the
prepared images have been visually reviewed with the user.

# Interpretability
User requires outputs to be interpretable. For local-feature diagnostics show
full-image context with labelled crop locations, consistent colours, explicit
reference-vs-prediction labels, and what numeric scores mean. Show the actual
patch sizes producing scores, explain any search over rotations/scales, and
separate evaluation reference geometry from algorithm predictions. Inspect
rendered figures for overlap/clipping before presenting. Do not advance stages
without explaining and reviewing the current result with the user.
