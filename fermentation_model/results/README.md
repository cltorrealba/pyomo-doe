# Results directory status

This directory contains both active priors and retained evidence. It is not a
single disposable build directory.

## Active model inputs and priors

These directories are read by current model or validation workflows:

- `new_must_data_loading/`
- `identifiability_reduction/`
- `new_must_glycerol_estimability_doe/`
- `new_must_glycerol_overnight_validation/`
- `secondary_metabolite_data_review/`
- `secondary_joint_campaign_doe/`
- `secondary_v2_model_evaluation/`

The CSV and PNG files directly under `results/` are the initial calibration,
PSO, profile-likelihood and q-sensitivity baseline. They are retained because
they are tracked evidence and some contain user modifications.

The PSO and broad multistart files are classified as transient development
evidence. The preferred method family is local IPOPT estimation followed by
weighted FIM/eigen-analysis and profile likelihood.

## Current operational and Lot 1 outputs

- `final_operational_doe_v2/`: active baseline for the constrained design.
- `final_operational_doe_volume_constrained/`: selected operational campaign.
- `design_execution_bundle_2026-06-09/`: frozen execution handoff.
- `lot1_data_preview/`
- `lot1_actual_mbdoe_reassessment/`
- `lot1_campaign_reassessment/`
- `lot1_pulse_timing_mbdoe/`
- `lot1_express_optimal_sampling/`
- `estimability_old_vs_lot1/`
- `vc_doe_nb_plots/`

## Retained development evidence

The following outputs are superseded as standalone workflows but remain here
because current reports or bundles reference them:

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

Do not move or delete a result family until `rg` confirms that no current
runner, notebook generator or evidence-bundle script references it.

Administrative `AXX` bundles are not scientific result families. Their retained
files are documented under `../legacy/rendicion/` and
`../pilot_2025/bundles/`.
