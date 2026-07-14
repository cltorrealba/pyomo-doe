# Legacy fermentation-model artifacts

Legacy material is preserved for traceability and comparison. It is not the
recommended starting point for calibration or MBDoE.

## `development_2026/`

This directory contains the progressive laboratory/model-development workflow
that preceded the current integrated pilot notebook:

- `notebooks/initial_calibration/`: original calibration and effective
  reformulation notebooks.
- `notebooks/model_development/`: data loading, glycerol, secondary-state and
  fit-capacity notebooks.
- `notebooks/preliminary_design/`: early extended and aroma DOE notebooks.
- `scripts/`: runners and utilities that generated those iterations.

No source file was deleted during archival. Historical notebooks may contain
absolute paths or references to their original root-level location; use them as
an audit record. The authoritative current paths are documented in
`../README.md`.

Before reviving an archived runner, verify its input paths and write outputs to
a new results directory. Do not overwrite current priors.
