# Superseded development results

These artifacts were previously mixed with active outputs under
`fermentation_model/results/`. They are preserved here with the archived 2026
development runners and notebooks that produced or consumed them.

## Initial calibration

`initial_calibration/` contains the former flat CSV, JSON and PNG outputs for
the initial calibration, PSO/multistart trials, FIM perturbations and profile
likelihood. They are retained for traceability and are not the current parameter
endpoint.

## Archived result families

- `curve_validation/`
- `aroma_campaign_doe/`
- `aroma_joint_campaign_doe/`
- `aroma_surrogate_regression_smoke/`
- `deep_model_selection/`
- `doe_experiment_design/`
- `extended_campaign_doe/`
- `extended_campaign_oed_refinement/`
- `fit_strategy_analysis/`
- `medium_transfer_diagnostics/`
- `operational_campaign_schedule/`
- `secondary_fit_capacity/`
- `literature_text/`

Do not import these artifacts into an active calibration by an implicit path.
If one must be revived, register it explicitly in `campaigns/workflows.csv` and
write its new output to an owned active result directory.
