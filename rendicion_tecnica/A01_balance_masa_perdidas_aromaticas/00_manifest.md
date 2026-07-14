# Manifest A02

Este archivo resume el inventario completo del bundle. La version tabular con hashes esta en `00_manifest.csv`.

## Conteo por categoria

| category            | status   |   n |
|:--------------------|:---------|----:|
| anexo               | ok       |   1 |
| cierre_masa         | ok       |   1 |
| codigo              | ok       |   7 |
| datos               | ok       |   5 |
| ecuaciones          | ok       |   1 |
| figuras             | ok       |  36 |
| guia                | ok       |   1 |
| matriz_cumplimiento | ok       |   1 |
| metricas            | ok       |   1 |
| notebooks           | ok       |   6 |
| parametros          | ok       |   2 |
| planilla            | ok       |   1 |
| reportes            | ok       |   8 |
| sensibilidad        | ok       |   1 |
| tablas_fuente       | ok       |  19 |
| version_codigo      | ok       |   1 |

## Primeras rutas

| category      | evidence_requirement                        | status   | bundle_path                                                                                                                     | source_path                                                                                                |
|:--------------|:--------------------------------------------|:---------|:--------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------|
| datos         | Datos de calibracion                        | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/data_sources/Calibration_data_vl3.xlsx                                   | fermentation_model/data/Calibration_data_vl3.xlsx                                                          |
| datos         | Datos de calibracion                        | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/data_sources/mosto_sintetico_vl3.xlsx                                    | fermentation_model/data/mosto_sintetico_vl3.xlsx                                                           |
| datos         | Datos de calibracion                        | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/data_sources/mosto_natural_xthiol.xlsx                                   | fermentation_model/data/mosto_natural_xthiol.xlsx                                                          |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/run_new_must_glycerol_estimability_doe.py                   | fermentation_model/run_new_must_glycerol_estimability_doe.py                                               |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/run_aroma_campaign_doe.py                                   | fermentation_model/run_aroma_campaign_doe.py                                                               |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/aroma_partition_unifac.py                                   | fermentation_model/aroma_partition_unifac.py                                                               |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/new_must_data_loader.py                                     | fermentation_model/new_must_data_loader.py                                                                 |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/run_secondary_fit_capacity.py                               | fermentation_model/run_secondary_fit_capacity.py                                                           |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/run_secondary_joint_campaign_doe.py                         | fermentation_model/run_secondary_joint_campaign_doe.py                                                     |
| codigo        | Version de codigo                           | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/code_sources/build_A02_evidence_bundle.py                                | scripts/build_A02_evidence_bundle.py                                                                       |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_model_calibration.ipynb                           | fermentation_model/fermentation_model_calibration.ipynb                                                    |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_model_calibration_3_effective_reformulation.ipynb | fermentation_model/fermentation_model_calibration_3_effective_reformulation.ipynb                          |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_new_must_glycerol_estimability_doe.executed.ipynb | fermentation_model/fermentation_new_must_glycerol_estimability_doe.executed.ipynb                          |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_secondary_fit_capacity.executed.ipynb             | fermentation_model/fermentation_secondary_fit_capacity.executed.ipynb                                      |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_secondary_joint_campaign_doe.executed.ipynb       | fermentation_model/fermentation_secondary_joint_campaign_doe.executed.ipynb                                |
| notebooks     | Metricas de ajuste y analisis               | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/notebooks/fermentation_aroma_campaign_doe.ipynb                          | fermentation_model/fermentation_aroma_campaign_doe.ipynb                                                   |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/nmg_doe__new_must_glycerol_estimability_doe_report.md     | fermentation_model/results/new_must_glycerol_estimability_doe/new_must_glycerol_estimability_doe_report.md |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/sec_fit__secondary_fit_capacity_report.md                 | fermentation_model/results/secondary_fit_capacity/secondary_fit_capacity_report.md                         |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/sec_joint__secondary_joint_campaign_report.md             | fermentation_model/results/secondary_joint_campaign_doe/secondary_joint_campaign_report.md                 |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/aroma__aroma_campaign_report.md                           | fermentation_model/results/aroma_campaign_doe/aroma_campaign_report.md                                     |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/aroma_joint__aroma_campaign_report.md                     | fermentation_model/results/aroma_joint_campaign_doe/aroma_campaign_report.md                               |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/aroma_joint__aroma_eigen_analysis_report.md               | fermentation_model/results/aroma_joint_campaign_doe/aroma_eigen_analysis_report.md                         |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/aroma_joint__aroma_dopt_benchmark_report.md               | fermentation_model/results/aroma_joint_campaign_doe/aroma_dopt_benchmark_report.md                         |
| reportes      | Metricas, sensibilidad y supuestos          | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_reports/results__aroma_symbolic_unifac_trial_report.md            | fermentation_model/results/aroma_symbolic_unifac_trial_report.md                                           |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/nmg_doe__fit_summary.csv                                   | fermentation_model/results/new_must_glycerol_estimability_doe/fit_summary.csv                              |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/nmg_doe__theta_fit_table.csv                               | fermentation_model/results/new_must_glycerol_estimability_doe/theta_fit_table.csv                          |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/nmg_doe__parameter_estimability_summary.csv                | fermentation_model/results/new_must_glycerol_estimability_doe/parameter_estimability_summary.csv           |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/nmg_doe__profile_summary_combined.csv                      | fermentation_model/results/new_must_glycerol_estimability_doe/profile_summary_combined.csv                 |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/sec_fit__fit_capacity_summary.csv                          | fermentation_model/results/secondary_fit_capacity/fit_capacity_summary.csv                                 |
| tablas_fuente | Parametros, cierre, metricas y sensibilidad | ok       | rendicion_tecnica/A02_balance_masa_perdidas_aromaticas/source_tables/sec_fit__fit_capacity_state_metrics.csv                    | fermentation_model/results/secondary_fit_capacity/fit_capacity_state_metrics.csv                           |
