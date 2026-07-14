# Manifest A02

La version tabular con hashes esta en `00_manifest.csv`.

## Conteo por categoria

| category            | status   |   n |
|:--------------------|:---------|----:|
| anexo               | ok       |   1 |
| brechas             | ok       |   2 |
| codigo              | ok       |   6 |
| criterios           | ok       |   1 |
| datos               | ok       |   7 |
| decision            | ok       |   2 |
| figuras             | ok       |  34 |
| guia                | ok       |   1 |
| matriz_cumplimiento | ok       |   1 |
| metricas            | ok       |   2 |
| notebooks           | ok       |   3 |
| predicciones        | ok       |   1 |
| predicho_observado  | ok       |   1 |
| protocolo           | ok       |   1 |
| reportes            | ok       |   6 |
| sensibilidad        | ok       |   1 |
| tabla_anexo         | ok       |   1 |
| tablas_fuente       | ok       |  17 |
| version_codigo      | ok       |   1 |

## Primeras rutas

| category      | evidence_requirement                        | status   | bundle_path                                                                                                                    | source_path                                                                                                     |
|:--------------|:--------------------------------------------|:---------|:-------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------|
| datos         | Base independiente vendimia 2025            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/data_sources/mosto_natural_xthiol.xlsx                                   | fermentation_model/data/mosto_natural_xthiol.xlsx                                                               |
| datos         | Base independiente vendimia 2025            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/data_sources/mosto_sintetico_vl3.xlsx                                    | fermentation_model/data/mosto_sintetico_vl3.xlsx                                                                |
| datos         | Base independiente vendimia 2025            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/data_sources/Calibration_data_vl3.xlsx                                   | fermentation_model/data/Calibration_data_vl3.xlsx                                                               |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/new_must_data_loader.py                                     | fermentation_model/new_must_data_loader.py                                                                      |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/run_new_must_curve_validation.py                            | fermentation_model/run_new_must_curve_validation.py                                                             |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/run_new_must_overnight_validation.py                        | fermentation_model/run_new_must_overnight_validation.py                                                         |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/run_new_must_glycerol_estimability_doe.py                   | fermentation_model/run_new_must_glycerol_estimability_doe.py                                                    |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/run_medium_transfer_diagnostics.py                          | fermentation_model/run_medium_transfer_diagnostics.py                                                           |
| codigo        | Protocolo de validacion y version de codigo | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/code_sources/build_A02_validation_bundle.py                              | scripts/build_A02_validation_bundle.py                                                                          |
| notebooks     | Predicho-observado y calibracion            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/notebooks/fermentation_new_must_data_loading.ipynb                       | fermentation_model/fermentation_new_must_data_loading.ipynb                                                     |
| notebooks     | Predicho-observado y calibracion            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/notebooks/fermentation_new_must_glycerol_estimability_doe.executed.ipynb | fermentation_model/fermentation_new_must_glycerol_estimability_doe.executed.ipynb                               |
| notebooks     | Predicho-observado y calibracion            | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/notebooks/fermentation_new_must_glycerol_estimability_doe.ipynb          | fermentation_model/fermentation_new_must_glycerol_estimability_doe.ipynb                                        |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/curve__curve_validation_report.md                         | fermentation_model/results/curve_validation/curve_validation_report.md                                          |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/curve__curve_validation_decision_summary.md               | fermentation_model/results/curve_validation/curve_validation_decision_summary.md                                |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/medium_transfer__medium_transfer_diagnostics_report.md    | fermentation_model/results/medium_transfer_diagnostics/medium_transfer_diagnostics_report.md                    |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/overnight__overnight_validation_report.md                 | fermentation_model/results/new_must_glycerol_overnight_validation/overnight_validation_report.md                |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/best_post__post_summary.md                                | fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/post_summary.md |
| reportes      | Metricas, brechas y decision tecnica        | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_reports/nmg_doe__new_must_glycerol_estimability_doe_report.md     | fermentation_model/results/new_must_glycerol_estimability_doe/new_must_glycerol_estimability_doe_report.md      |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/data_loading__new_must_batch_summary.csv                   | fermentation_model/results/new_must_data_loading/new_must_batch_summary.csv                                     |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/data_loading__new_must_normalized_long.csv                 | fermentation_model/results/new_must_data_loading/new_must_normalized_long.csv                                   |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/data_loading__density_sugar_fit_by_medium.csv              | fermentation_model/results/new_must_data_loading/density_sugar_fit_by_medium.csv                                |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/curve__current_fit_curve_metrics.csv                       | fermentation_model/results/curve_validation/current_fit_curve_metrics.csv                                       |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/curve__current_fit_curve_metrics_by_medium_state.csv       | fermentation_model/results/curve_validation/current_fit_curve_metrics_by_medium_state.csv                       |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/curve__design_prediction_feasibility.csv                   | fermentation_model/results/curve_validation/design_prediction_feasibility.csv                                   |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/curve__design_theta_disagreement.csv                       | fermentation_model/results/curve_validation/design_theta_disagreement.csv                                       |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/curve__design_pairwise_similarity.csv                      | fermentation_model/results/curve_validation/design_pairwise_similarity.csv                                      |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/medium_transfer__historical_batch_audit.csv                | fermentation_model/results/medium_transfer_diagnostics/historical_batch_audit.csv                               |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/medium_transfer__medium_map_template.csv                   | fermentation_model/results/medium_transfer_diagnostics/medium_map_template.csv                                  |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/overnight__multistart_summary.csv                          | fermentation_model/results/new_must_glycerol_overnight_validation/multistart_summary.csv                        |
| tablas_fuente | Base, predicho-observado, metricas, brechas | ok       | rendicion_tecnica/A02_validacion_modelo_vendimia_2025/source_tables/overnight__l2_scan_summary.csv                             | fermentation_model/results/new_must_glycerol_overnight_validation/l2_scan_summary.csv                           |
