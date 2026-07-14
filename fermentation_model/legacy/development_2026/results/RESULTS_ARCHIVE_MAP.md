# Results directory status (pre-reorganization snapshot)

This document preserves the navigation note that accompanied the former mixed
`fermentation_model/results/` directory. Its paths describe the repository
before the ownership reorganization and are retained only as historical
provenance. Use `fermentation_model/REPOSITORY_MAP.md` for current navigation.

## Active model inputs and priors at the time

These directories were read by model or validation workflows:

- `new_must_data_loading/`
- `identifiability_reduction/`
- `new_must_glycerol_estimability_doe/`
- `new_must_glycerol_overnight_validation/`
- `secondary_metabolite_data_review/`
- `secondary_joint_campaign_doe/`
- `secondary_v2_model_evaluation/`

The CSV and PNG files directly under the old `results/` directory were the
initial calibration, PSO, profile-likelihood and q-sensitivity baseline. They
were retained because they are tracked evidence and some contain user
modifications.

The PSO and broad multistart files were classified as transient development
evidence. The preferred method family was local IPOPT estimation followed by
weighted FIM/eigen-analysis and profile likelihood.

## Current operational and Lot 1 outputs at the time

- `final_operational_doe_v2/`
- `final_operational_doe_volume_constrained/`
- `design_execution_bundle_2026-06-09/`
- `lot1_data_preview/`
- `lot1_actual_mbdoe_reassessment/`
- `lot1_campaign_reassessment/`
- `lot1_pulse_timing_mbdoe/`
- `lot1_express_optimal_sampling/`
- `estimability_old_vs_lot1/`
- `vc_doe_nb_plots/`

## Retained development evidence at the time

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

Administrative `AXX` bundles were not scientific result families. Their files
remain under `fermentation_model/legacy/rendicion/` and
`fermentation_model/pilot_2025/bundles/`.
