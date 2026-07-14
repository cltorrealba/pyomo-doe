# Version de codigo y trazabilidad

Generado: 2026-06-22T14:57:10.841188-04:00
Fecha de cierre documental usada para la rendicion: 2026-06-22

## Git

- Rama: `ctorrealba_fermentation`
- Commit HEAD: `361ecdcc8a6ae9744e8ae10c4b5b127e57001e44`
- Commit corto: `361ecdc`

## Estado de trabajo

```text
M fermentation_model/results/fim_weighted_relative_eigendirections_nominal_pm20.csv
 M fermentation_model/results/fim_weighted_relative_eigendirections_wsse_fit_pm20.csv
 M fermentation_model/results/fim_weighted_relative_eigenvalues_nominal_pm20.csv
 M fermentation_model/results/fim_weighted_relative_eigenvalues_wsse_fit_pm20.csv
 M fermentation_model/results/fim_weighted_relative_nominal_pm20.csv
 M fermentation_model/results/fim_weighted_relative_wsse_fit_pm20.csv
 M fermentation_model/results/q_perturbation_summary_nominal_pm20.csv
 M fermentation_model/results/q_perturbation_summary_wsse_fit_pm20.csv
 M fermentation_model/results/q_sensitivity_raw_nominal_pm20.csv
 M fermentation_model/results/q_sensitivity_raw_wsse_fit_pm20.csv
 M fermentation_model/results/q_sensitivity_solve_summary_nominal_pm20.csv
 M fermentation_model/results/q_sensitivity_solve_summary_wsse_fit_pm20.csv
 M fermentation_model/results/q_sensitivity_weighted_nominal_pm20.csv
 M fermentation_model/results/q_sensitivity_weighted_relative_nominal_pm20.csv
 M fermentation_model/results/q_sensitivity_weighted_relative_wsse_fit_pm20.csv
 M fermentation_model/results/q_sensitivity_weighted_wsse_fit_pm20.csv
?? db_20260609.zip
?? db_20260609/
?? fermentation_model/aroma_partition_unifac.py
?? fermentation_model/create_volume_constrained_doe_notebook.py
?? fermentation_model/data/mosto_natural_xthiol.xlsx
?? fermentation_model/data/mosto_sintetico_vl3.xlsx
?? fermentation_model/fermentation_aroma_campaign_doe.ipynb
?? fermentation_model/fermentation_extended_campaign_doe.ipynb
?? fermentation_model/fermentation_final_operational_doe_v2.executed.ipynb
?? fermentation_model/fermentation_final_operational_doe_v2.ipynb
?? fermentation_model/fermentation_final_operational_doe_volume_constrained.executed.ipynb
?? fermentation_model/fermentation_final_operational_doe_volume_constrained.ipynb
?? fermentation_model/fermentation_model_calibration_3_effective_reformulation.ipynb
?? fermentation_model/fermentation_new_must_data_loading.ipynb
?? fermentation_model/fermentation_new_must_glycerol_estimability_doe.executed.ipynb
?? fermentation_model/fermentation_new_must_glycerol_estimability_doe.ipynb
?? fermentation_model/fermentation_secondary_fit_capacity.executed.ipynb
?? fermentation_model/fermentation_secondary_fit_capacity.ipynb
?? fermentation_model/fermentation_secondary_joint_campaign_doe.executed.ipynb
?? fermentation_model/fermentation_secondary_joint_campaign_doe.ipynb
?? fermentation_model/fermentation_secondary_metabolite_data_review.executed.ipynb
?? fermentation_model/fermentation_secondary_metabolite_data_review.ipynb
?? fermentation_model/fermentation_secondary_v2_model_evaluation.executed.ipynb
?? fermentation_model/fermentation_secondary_v2_model_evaluation.ipynb
?? fermentation_model/new_must_data_loader.py
?? fermentation_model/results/aroma_campaign_doe/
?? fermentation_model/results/aroma_joint_campaign_doe/
?? fermentation_model/results/aroma_surrogate_regression_smoke/
?? fermentation_model/results/aroma_symbolic_unifac_trial_report.md
?? fermentation_model/results/curve_validation/
?? fermentation_model/results/deep_model_selection/
?? fermentation_model/results/design_execution_bundle_2026-06-09/
?? fermentation_model/results/doe_experiment_design/
?? fermentation_model/results/extended_campaign_doe/
?? fermentation_model/results/extended_campaign_oed_refinement/
?? fermentation_model/results/fim_weighted_relative_eigendirections_nominal_pm10.csv
?? fermentation_model/results/fim_weighted_relative_eigendirections_wsse_fit_pm10.csv
?? fermentation_model/results/fim_weighted_relative_eigenvalues_nominal_pm10.csv
?? fermentation_model/results/fim_weighted_relative_eigenvalues_wsse_fit_pm10.csv
?? fermentation_model/results/fim_weighted_relative_nominal_pm10.csv
?? fermentation_model/results/fim_weighted_relative_wsse_fit_pm10.csv
?? fermentation_model/results/final_operational_doe_v2/
?? fermentation_model/results/final_operational_doe_volume_constrained/
?? fermentation_model/results/fit_strategy_analysis/
?? fermentation_model/results/identifiability_reduction/
?? fermentation_model/results/literature_text/
?? fermentation_model/results/lot1_express_optimal_sampling/
?? fermentation_model/results/lot1_pulse_timing_mbdoe/
?? fermentation_model/results/medium_transfer_diagnostics/
?? fermentation_model/results/new_must_data_loading/
?? fermentation_model/results/new_must_glycerol_estimability_doe/
?? fermentation_model/results/new_must_glycerol_overnight_validation/
?? fermentation_model/results/operational_campaign_schedule/
?? fermentation_model/results/q_perturbation_summary_nominal_pm10.csv
?? fermentation_model/results/q_perturbation_summary_wsse_fit_pm10.csv
?? fermentation_model/results/q_sensitivity_raw_nominal_pm10.csv
?? fermentation_model/results/q_sensitivity_raw_wsse_fit_pm10.csv
?? fermentation_model/results/q_sensitivity_solve_summary_nominal_pm10.csv
?? fermentation_model/results/q_sensitivity_solve_summary_wsse_fit_pm10.csv
?? fermentation_model/results/q_sensitivity_weighted_nominal_pm10.csv
?? fermentation_model/results/q_sensitivity_weighted_relative_nominal_pm10.csv
?? fermentation_model/results/q_sensitivity_weighted_relative_wsse_fit_pm10.csv
?? fermentation_model/results/q_sensitivity_weighted_wsse_fit_pm10.csv
?? fermentation_model/results/secondary_fit_capacity/
?? fermentation_model/results/secondary_joint_campaign_doe/
?? fermentation_model/results/secondary_metabolite_data_review/
?? fermentation_model/results/secondary_v2_model_evaluation/
?? fermentation_model/results/vc_doe_nb_plots/
?? fermentation_model/run_aroma_campaign_doe.py
?? fermentation_model/run_aroma_dopt_benchmark.py
?? fermentation_model/run_aroma_eigen_analysis.py
?? fermentation_model/run_batch_prioritization.py
?? fermentation_model/run_deep_model_selection.py
?? fermentation_model/run_doe_experiment_design.py
?? fermentation_model/run_extended_campaign_doe.py
?? fermentation_model/run_extended_campaign_oed_refinement.py
?? fermentation_model/run_final_operational_doe_v2.py
?? fermentation_model/run_final_operational_doe_volume_constrained.py
?? fermentation_model/run_fit_strategy_analysis.py
?? fermentation_model/run_lot1_express_optimal_sampling.py
?? fermentation_model/run_lot1_pulse_timing_mbdoe.py
?? fermentation_model/run_medium_transfer_diagnostics.py
?? fermentation_model/run_new_must_curve_validation.py
?? fermentation_model/run_new_must_glycerol_estimability_doe.py
?? fermentation_model/run_new_must_overnight_validation.py
?? fermentation_model/run_operational_campaign_schedule.py
?? fermentation_model/run_secondary_fit_capacity.py
?? fermentation_model/run_secondary_joint_campaign_doe.py
?? fermentation_model/run_secondary_metabolite_data_review.py
?? fermentation_model/run_secondary_v2_model_evaluation.py
?? rendicion_tecnica/
?? scripts/build_A02_evidence_bundle.py
```

## Diff stat

```text
...ghted_relative_eigendirections_nominal_pm20.csv |  30 +-
 ...hted_relative_eigendirections_wsse_fit_pm20.csv |  15 +-
 ..._weighted_relative_eigenvalues_nominal_pm20.csv |  30 +-
 ...weighted_relative_eigenvalues_wsse_fit_pm20.csv |  15 +-
 .../results/fim_weighted_relative_nominal_pm20.csv |  32 +-
 .../fim_weighted_relative_wsse_fit_pm20.csv        |  17 +-
 .../q_perturbation_summary_nominal_pm20.csv        |  20 +-
 .../q_perturbation_summary_wsse_fit_pm20.csv       |  16 +-
 .../results/q_sensitivity_raw_nominal_pm20.csv     | 534 +++++++++----------
 .../results/q_sensitivity_raw_wsse_fit_pm20.csv    | 574 ++++++++++-----------
 .../q_sensitivity_solve_summary_nominal_pm20.csv   | 160 +++---
 .../q_sensitivity_solve_summary_wsse_fit_pm20.csv  |  80 +--
 .../q_sensitivity_weighted_nominal_pm20.csv        | 534 +++++++++----------
 ..._sensitivity_weighted_relative_nominal_pm20.csv | 534 +++++++++----------
 ...sensitivity_weighted_relative_wsse_fit_pm20.csv | 574 ++++++++++-----------
 .../q_sensitivity_weighted_wsse_fit_pm20.csv       | 574 ++++++++++-----------
 16 files changed, 1880 insertions(+), 1859 deletions(-)
```

## Hashes SHA256 de fuentes clave

- `fermentation_model/run_new_must_glycerol_estimability_doe.py`: `7e6cb8e44e9ea5055b06bc09be8d2fcf3d951da779d049a2ccddd97c3328ae81`
- `fermentation_model/run_aroma_campaign_doe.py`: `6af582733361433392d3b0114c1e60ba9503c430839d582f0e09c7d8e4594538`
- `fermentation_model/aroma_partition_unifac.py`: `25bc43b0bc8386fae89240b58db90af389bb39ec330318609b27460d6de18fca`
- `fermentation_model/new_must_data_loader.py`: `f13c31e95313f767efde1d2c67079abf18c694b17d51880501c9e240551dd0bf`
- `fermentation_model/data/Calibration_data_vl3.xlsx`: `34167748df07b37474546305520f3f3ad9812ffc53835a519d60dda3f6dcf84b`
- `fermentation_model/data/mosto_sintetico_vl3.xlsx`: `1bbf657740503af13facefea7a1e1cae52719e6ff4c570771b678562f170b4a2`

## Nota de interpretacion

El arbol de trabajo contiene archivos sin seguimiento y archivos modificados. Para rendicion, usar este documento junto con `00_manifest.csv`, que registra hashes de los artefactos copiados al bundle. Si se requiere congelar formalmente la evidencia, crear commit/tag o almacenar el ZIP en el gestor documental del proyecto.
