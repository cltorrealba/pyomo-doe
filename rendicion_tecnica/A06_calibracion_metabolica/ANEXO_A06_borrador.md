# Calibracion del submodelo metabolico del Gemelo Digital (Actividad 3.13)

## Identificacion

- ID anexo: A06
- Objetivo especifico: OE3
- Actividad unica: 3.13
- Nombre de archivo sugerido: `PI-4497_IA5_A06_Act-3.13_Calibracion-submodelo-metabolico.pdf`
- Fecha de cierre documental del bundle: 2026-06-26

## Tabla fuente del anexo

| ID anexo   | OE   | Actividad unica   | Titulo propuesto                                                         | Nombre de archivo sugerido                                    | Contenido minimo                                                                                                 | Evidencia fuente                                                                                                          |
|:-----------|:-----|:------------------|:-------------------------------------------------------------------------|:--------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------|
| A06        | OE3  | 3.13              | Calibracion del submodelo metabolico del Gemelo Digital (Actividad 3.13) | PI-4497_IA5_A06_Act-3.13_Calibracion-submodelo-metabolico.pdf | Documentar dataset, estrategia de calibracion, parametros metabolicos, metricas, incertidumbre e interpretacion. | Dataset de calibracion, parametros metabolicos, funcion objetivo, metricas, incertidumbre/identificabilidad y resultados. |

## Resumen ejecutivo

Se documenta la calibracion del submodelo metabolico del Gemelo Digital a partir de bases historicas VL3, ensayos de mosto natural y ensayos de mosto sintetico. La evidencia incluye dataset integrado, estrategia de calibracion regularizada, parametros metabolicos estimados, funcion objetivo, metricas de ajuste, diagnosticos de incertidumbre/identificabilidad y una interpretacion tecnica de uso condicionado.

## Matriz de cumplimiento

| contenido_minimo                | estado_bundle   | archivos_bundle                                                                                                                | detalle                                                                                                                 |
|:--------------------------------|:----------------|:-------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------|
| Dataset de calibracion          | cubierto        | data_sources/*.xlsx; source_tables/secondary_joint_campaign_doe__secondary_joint_input_data.csv; tables/03_dataset_resumen.csv | Incluye base historica VL3, mosto natural, mosto sintetico y dataset integrado de 27 lotes para calibracion secundaria. |
| Estrategia de calibracion       | cubierto        | tables/04_estrategia_calibracion.csv; 03_metodologia_calibracion_A06.md; source_reports/*.md                                   | Documenta ajuste regularizado, evaluacion v2, diagnostico de capacidad y validacion por multistart/perfiles.            |
| Parametros metabolicos          | cubierto        | tables/06_parametros_metabolicos.csv; source_tables/*theta*.csv; source_tables/*estimability*.csv                              | Incluye parametros del ajuste regularizado y parametros v2 de piruvato, acetaldehido, acetato y oxigeno.                |
| Funcion objetivo                | cubierto        | tables/05_funcion_objetivo.csv; code_sources/run_fit_strategy_analysis.py; code_sources/run_secondary_v2_model_evaluation.py   | Describe WSSE, penalizacion log-L2, residuos escalados por sigma y diagnosticos de robustez.                            |
| Metricas                        | cubierto        | tables/07_metricas_calibracion.csv; source_tables/*fit*.csv; source_tables/*state*.csv                                         | Incluye WSSE, SSE, residuos por estado, n observaciones, convergencia y comparacion de modelos.                         |
| Incertidumbre/identificabilidad | cubierto        | tables/08_incertidumbre_identificabilidad.csv; source_tables/*profile*.csv; source_tables/*fim*.csv                            | Incluye FIM, espectro propio, perfiles de verosimilitud y clasificacion de parametros debiles/confundidos.              |
| Resultados e interpretacion     | cubierto        | tables/09_resultados_interpretacion.csv; 04_resultados_interpretacion_A06.md; ANEXO_A06_borrador.md                            | Resume decision tecnica, uso recomendado, brechas y limites de cierre parametrico.                                      |

## Dataset de calibracion

| archivo                        | estado     | n_hojas   | filas_totales   | n_lotes   | medios                           | hojas_resumen                                                                                                                                                                                                                                                               |
|:-------------------------------|:-----------|:----------|:----------------|:----------|:---------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Calibration_data_vl3.xlsx      | disponible | 8         | 312             |           |                                  | 25026:24; 25027:32; 25085:35; 25086:35; 25150:30; 25151:30; 25170:63; 25171:63                                                                                                                                                                                              |
| mosto_natural_xthiol.xlsx      | disponible | 20        | 610             |           |                                  | 00_Resumen:20; Datos_homologados:104; LAB004:14; LAB005:14; LAB006:16; LAB007:14; LAB008:14; LAB009:13; LAB010:10; LAB011:9; LAB012:9; Perfiles_temperatura:25                                                                                                              |
| mosto_sintetico_vl3.xlsx       | disponible | 16        | 310             |           |                                  | 00_Resumen:10; Datos_homologados:86; MS007:10; MS008:10; MS009:10; MS010:10; MS011:10; MS012:10; MS013:9; MS014:9; MS015:9; MS016:9                                                                                                                                         |
| secondary_joint_input_data.csv | disponible |           | 520             | 27        | historical_vl3,natural,synthetic | time_h:520; temperature_c:506; density:496; X_viable_kg_m3:307; G_g_l:315; F_g_l:315; YAN_mg_l:308; pyruvic_acid:313; acetaldehyde:193; acetic_acid:50; DO_mg_l:32; glycerol_g_l:314; E_g_l:222; ethyl_acetate_total:50; isoamyl_acetate_total:50; ethyl_octanoate_total:50 |

## Estrategia de calibracion

| etapa                      | descripcion                                                                                | parametros_o_estados                                                      | n_lotes_o_obs                              | criterio_decision                                                                                  | archivo_fuente                                                  |
|:---------------------------|:-------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------|:-------------------------------------------|:---------------------------------------------------------------------------------------------------|:----------------------------------------------------------------|
| 1_dataset_integrado        | Integracion de datos historicos VL3, mosto natural y mosto sintetico en un dataset comun.  | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE...                                     | 27 lotes en secondary_joint_input_data.csv | Dataset integrado disponible y trazable a libros fuente.                                           | secondary_joint_campaign_doe/secondary_joint_input_data.csv     |
| 2_ajuste_regularizado_core | Comparacion de estrategias ParmEst con WSSE y penalizacion log-L2.                         | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG                           | 286 observaciones                          | Seleccion de plus_qx_iG_logL2_10 por no activar bounds y mantener objetivo penalizado competitivo. | fit_strategy_analysis/fit_strategy_summary.csv                  |
| 3_submodelo_secundario_v2  | Calibracion v2 para piruvato, acetaldehido, acetato y O2 con simulacion core fija.         | kPyrS_N,kPyrO2,kPyrDrain,kAldS_N,kAldRed,kAcAld,kAcStress                 | 27 lotes; 588 residuos                     | Comparar WSSE frente a v1 y clasificar parametros por FIM.                                         | secondary_v2_model_evaluation/fit_comparison.csv                |
| 4_capacidad_de_ajuste      | Ajustes por medio y por lote como cota superior de capacidad, no como modelo transferible. | kPyrS_N,kPyrO2,kPyrDrain,kAldS_N,kAldRed,kAcAld,kAcStress,qO2,kLaO2,O2sat | 27 lotes                                   | per-batch all-free fit is an upper bound on fit capacity, not an identifiable model                | secondary_fit_capacity/fit_capacity_summary.csv                 |
| 5_incertidumbre            | Multistart, escaneo L2, perfiles de verosimilitud y FIM para robustez/identificabilidad.   | full17; v2_reduced_o2fixed                                                | 14 multistarts; 17 perfiles full17; FIM v2 | Distinguir parametros identificables, regulares y debiles/confundidos.                             | new_must_glycerol_overnight_validation/profile_full_summary.csv |

## Funcion objetivo

| workflow                      | funcion_objetivo                               | formula_operativa                                                                             | residuos_o_estados   | ponderacion                                                              | archivo_codigo                                    |
|:------------------------------|:-----------------------------------------------|:----------------------------------------------------------------------------------------------|:---------------------|:-------------------------------------------------------------------------|:--------------------------------------------------|
| fit_strategy_analysis         | WSSE + penalizacion log-L2                     | ParmEst SSE_weighted + 0.5 * lambda * sum(log(theta/theta_ref)^2) * n_batches                 | E,F,G,N,X            | residuos ponderados por funcion SSE_weighted del notebook de calibracion | code_sources/run_fit_strategy_analysis.py         |
| secondary_v2_model_evaluation | least_squares robusto sobre residuos escalados | sum(((pred-obs)/sigma_estado)^2) con loss soft_l1 y priors log para parametros seleccionados  | Pyr,AcAld,Acetate,O2 | SIGMA: Pyr=6.0, AcAld=8.0, Acetate=0.035, O2=0.20                        | code_sources/run_secondary_v2_model_evaluation.py |
| secondary_fit_capacity        | least_squares por medio/lote                   | misma base de residuos escalados; parametros liberados por lote como diagnostico de capacidad | Pyr,AcAld,Acetate,O2 | no se usa como modelo final; solo cota superior de ajuste                | code_sources/run_secondary_fit_capacity.py        |

## Parametros metabolicos

| familia                    | parameter   | value              | lower_bound   | upper_bound   | classification                   | source                                                 | std_log_approx   | approx_95_multiplier   | active_bound   |
|:---------------------------|:------------|:-------------------|:--------------|:--------------|:---------------------------------|:-------------------------------------------------------|:-----------------|:-----------------------|:---------------|
| core_regularizado_theta    | mu0         | 0.155377808467777  |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | sN          | 18.0               |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | qN          | 0.0133101185299815 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | qXG         | 0.0755014543706534 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | qXF         | 0.0916032749352908 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | betaG0      | 0.3995798062665501 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | sG          | 0.03               |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | betaF0      | 0.2889247532908035 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | sF          | 0.03               |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | qEG         | 0.5230426651395569 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | qEF         | 0.5200605577299623 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | iG          | 0.0201896772309379 |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | iE          | 0.025              |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | Kd0         | 0.00044            |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_theta    | m0          | 0.01               |               |               | seleccionado_plus_qx_iG_logL2_10 | fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv    |                  |                        |                |
| core_regularizado_physical | mu0         | 0.155377808467777  | 0.05          | 1.0           | valor_fisico_retrotransformado   | fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv |                  |                        |                |
| core_regularizado_physical | betaG0      | 0.3995798062665501 | 0.1           | 10.0          | valor_fisico_retrotransformado   | fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv |                  |                        |                |
| core_regularizado_physical | betaF0      | 0.2889247532908035 | 0.1           | 10.0          | valor_fisico_retrotransformado   | fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv |                  |                        |                |
| core_regularizado_physical | Kn0         | 0.008632100470432  | 0.01          | 0.1           | valor_fisico_retrotransformado   | fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv |                  |                        |                |
| core_regularizado_physical | Kg0         | 13.31932687555167  | 0.01          | 10.0          | valor_fisico_retrotransformado   | fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv |                  |                        |                |

## Metricas

| workflow                      | modelo_o_estrategia           | metrica_principal    | valor              | wsse               | n_obs   | estado      | comentario                                                              |
|:------------------------------|:------------------------------|:---------------------|:-------------------|:-------------------|:--------|:------------|:------------------------------------------------------------------------|
| fit_strategy_analysis         | plus_qx_iG_unregularized      | direct_plus_l2_total | 5972.547411621428  | 5972.547411621428  | 286     | ok          | active_bounds=2; lambda=0.0                                             |
| fit_strategy_analysis         | plus_qx_unregularized         | direct_plus_l2_total | 5983.45484334653   | 5983.45484334653   | 286     | ok          | active_bounds=2; lambda=0.0                                             |
| fit_strategy_analysis         | plus_qx_iG_logL2_10           | direct_plus_l2_total | 5998.054500708604  | 5992.323835014267  | 286     | ok          | active_bounds=0; lambda=10.0                                            |
| fit_strategy_analysis         | plus_iG_unregularized         | direct_plus_l2_total | 6001.98278070728   | 6001.98278070728   | 286     | ok          | active_bounds=0; lambda=0.0                                             |
| fit_strategy_analysis         | plus_qx_iG_pso_pilot_logL2_50 | direct_plus_l2_total | 6005.649337491055  | 6002.10466303318   | 286     | ok          | active_bounds=0; lambda=50.0                                            |
| fit_strategy_analysis         | plus_qx_iG_logL2_50           | direct_plus_l2_total | 6005.649340729168  | 6002.101747040763  | 286     | ok          | active_bounds=0; lambda=50.0                                            |
| fit_strategy_analysis         | plus_iG_logL2_50              | direct_plus_l2_total | 6007.175048302548  | 6005.058199919435  | 286     | ok          | active_bounds=0; lambda=50.0                                            |
| fit_strategy_analysis         | plus_qx_iG_logL2_100          | direct_plus_l2_total | 6007.718573729761  | 6005.283058791395  | 286     | ok          | active_bounds=0; lambda=100.0                                           |
| secondary_v2_model_evaluation | secondary_v1_current          | final_wsse           | 8706.95625298179   | 8706.95625298179   | 588     | ok          | wsse_per_residual=14.807748729560869; mejora_vs_v1_pct=                 |
| secondary_v2_model_evaluation | secondary_v2_phase_o2free     | final_wsse           | 5400.341991941201  | 5400.341991941201  | 592     | ok          | wsse_per_residual=9.122199310711489; mejora_vs_v1_pct=37.976695471603   |
| secondary_v2_model_evaluation | secondary_v2_phase_o2fixed    | final_wsse           | 5756.227087268427  | 5756.227087268427  | 590     | ok          | wsse_per_residual=9.75631709706513; mejora_vs_v1_pct=33.88933032370359  |
| secondary_v2_model_evaluation | secondary_v2_reduced_o2fixed  | final_wsse           | 5756.427666344939  | 5756.427666344939  | 589     | ok          | wsse_per_residual=9.773221844388692; mejora_vs_v1_pct=33.88702665901659 |
| secondary_fit_capacity        | batch_v2_all_o2free           | wsse_per_obs         | 1.629070050640295  | 957.8931897764936  | 588     | diagnostico | Cota de capacidad de ajuste; no equivale a modelo transferible.         |
| secondary_fit_capacity        | batch_v2_reduced_o2fixed      | wsse_per_obs         | 2.2937772172759057 | 1348.7410037582324 | 588     | diagnostico | Cota de capacidad de ajuste; no equivale a modelo transferible.         |
| secondary_fit_capacity        | global_v2_reduced_o2fixed     | wsse_per_obs         | 9.785874700216558  | 5754.094323727336  | 588     | diagnostico | Cota de capacidad de ajuste; no equivale a modelo transferible.         |
| secondary_fit_capacity        | medium_v2_reduced_o2fixed     | wsse_per_obs         | 8.038957228960598  | 4726.906850628831  | 588     | diagnostico | Cota de capacidad de ajuste; no equivale a modelo transferible.         |
| overnight_multistart          | multistart_07                 | final_wsse           | 23914.053825661373 | 23914.053825661373 | 1160    | ok          | wsse_per_residual=20.615563642811527; lambda=1.0                        |
| overnight_multistart          | multistart_06                 | final_wsse           | 24775.172321819    | 24775.172321819    | 1160    | ok          | wsse_per_residual=21.357907173981896; lambda=1.0                        |
| overnight_multistart          | multistart_05                 | final_wsse           | 24864.51536680058  | 24864.51536680058  | 1160    | ok          | wsse_per_residual=21.43492704034533; lambda=1.0                         |

## Incertidumbre e identificabilidad

| workflow                  | parameter   | diagnostico             | metrica              | valor              | detalle                                               |
|:--------------------------|:------------|:------------------------|:---------------------|:-------------------|:------------------------------------------------------|
| secondary_v2_FIM          | kPyrS_N     | well_estimated          | approx_95_multiplier | 1.0772804868923502 | std_log_approx=0.0379794886628891; active_bound=False |
| secondary_v2_FIM          | kPyrO2      | well_estimated          | approx_95_multiplier | 1.2180613925056278 | std_log_approx=0.100643149170135; active_bound=False  |
| secondary_v2_FIM          | kPyrDrain   | well_estimated          | approx_95_multiplier | 1.1256804543597017 | std_log_approx=0.0604018883311731; active_bound=False |
| secondary_v2_FIM          | kAldS_N     | well_estimated          | approx_95_multiplier | 1.0260821591294889 | std_log_approx=0.0131366431973644; active_bound=False |
| secondary_v2_FIM          | kAldRed     | weak_or_confounded      | approx_95_multiplier | 20752.958470953366 | std_log_approx=5.071655149188335; active_bound=False  |
| secondary_v2_FIM          | kAcAld      | well_estimated          | approx_95_multiplier | 1.1016531228114743 | std_log_approx=0.0493938217303898; active_bound=False |
| secondary_v2_FIM          | kAcStress   | well_estimated          | approx_95_multiplier | 1.06400411514056   | std_log_approx=0.0316526829210621; active_bound=False |
| profile_likelihood_full17 | mu0         | no_identifiable_95      | max_lr_stat          | 1878.6079058815667 | left=False; right=True; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | sN          | no_identifiable_95      | max_lr_stat          | 5.950374106912932  | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | qN          | profile_identifiable_95 | max_lr_stat          | 5840.911199476752  | left=True; right=True; chi2_95=3.841458820694124      |
| profile_likelihood_full17 | qXG         | no_identifiable_95      | max_lr_stat          | 0.0                | left=False; right=False; chi2_95=3.841458820694124    |
| profile_likelihood_full17 | qXF         | no_identifiable_95      | max_lr_stat          | 17.857437932496396 | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | betaG0      | no_identifiable_95      | max_lr_stat          | 0.0                | left=False; right=False; chi2_95=3.841458820694124    |
| profile_likelihood_full17 | sG          | no_identifiable_95      | max_lr_stat          | 105.51938168677589 | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | betaF0      | no_identifiable_95      | max_lr_stat          | 424.0279806899234  | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | sF          | no_identifiable_95      | max_lr_stat          | 531.0346761036053  | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | qEG         | no_identifiable_95      | max_lr_stat          | 3934.3556397044704 | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | qEF         | no_identifiable_95      | max_lr_stat          | 4017.488352924018  | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | iG          | no_identifiable_95      | max_lr_stat          | 129.04461635522603 | left=False; right=True; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | iE          | no_identifiable_95      | max_lr_stat          | 228.24435074699431 | left=True; right=False; chi2_95=3.841458820694124     |
| profile_likelihood_full17 | Kd0         | profile_identifiable_95 | max_lr_stat          | 374.761028763769   | left=True; right=True; chi2_95=3.841458820694124      |
| profile_likelihood_full17 | m0          | no_identifiable_95      | max_lr_stat          | 0.0                | left=False; right=False; chi2_95=3.841458820694124    |
| profile_likelihood_full17 | gammaG0     | profile_identifiable_95 | max_lr_stat          | 246.70044514316032 | left=True; right=True; chi2_95=3.841458820694124      |
| profile_likelihood_full17 | gammaF0     | profile_identifiable_95 | max_lr_stat          | 14.154663436827832 | left=True; right=True; chi2_95=3.841458820694124      |

## Interpretacion

| tema               | resultado                                                                                                                                    | soporte                                                     | interpretacion                                                                                                 |
|:-------------------|:---------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|
| cierre_documental  | Cobertura 100% del contenido minimo solicitado para A06.                                                                                     | 02_matriz_cumplimiento_A06.csv; 00_manifest.csv             | El bundle documenta dataset, estrategia, parametros, funcion objetivo, metricas, incertidumbre y resultados.   |
| modelo_recomendado | plus_qx_iG_logL2_10 como candidato regularizado para ajuste core; secondary_v2_reduced_o2fixed para capa secundaria.                         | fit_strategy_summary.csv; fit_comparison.csv; metadata.json | Se privilegia estabilidad sin bounds activos y reduccion de error frente al modelo secundario v1.              |
| mejora_v2          | secondary_v1_current=8.71e+03; secondary_v2_phase_o2free=5.4e+03; secondary_v2_phase_o2fixed=5.76e+03; secondary_v2_reduced_o2fixed=5.76e+03 | tables/07_metricas_calibracion.csv                          | La estructura v2 reduce el WSSE global respecto del ajuste secundario v1, pero conserva parametros debiles.    |
| identificabilidad  | kAldRed, mu0, sN, qXG, qXF, betaG0, sG, betaF0, sF, qEG, qEF, iG...                                                                          | tables/08_incertidumbre_identificabilidad.csv               | Los parametros debiles o no identificables deben mantenerse fijos, regularizados o tratarse como sensibilidad. |
| decision_tecnica   | Calibracion tecnica documentada con uso condicionado para simulacion, DOE y decision secuencial.                                             | 04_resultados_interpretacion_A06.md                         | No conviene declarar cierre parametrico definitivo para todos los parametros metabolicos.                      |

## Cumplimiento

| criterio                       | estado   | valor               | comentario                                                                                        |
|:-------------------------------|:---------|:--------------------|:--------------------------------------------------------------------------------------------------|
| Cobertura contenido minimo A06 | cumple   | 7/7 items cubiertos | Cobertura documental completa; las limitaciones tecnicas quedan explicitadas como interpretacion. |

## Texto sugerido para informe

Se consolido evidencia de calibracion del submodelo metabolico del Gemelo Digital para la actividad OE3-3.13. El dataset integrado combina informacion historica VL3, mosto natural y mosto sintetico, y se usa para evaluar parametros asociados a estados centrales y secundarios de fermentacion. La estrategia documentada incluye ajuste regularizado mediante WSSE con penalizacion log-L2, evaluacion del submodelo secundario v2 para piruvato, acetaldehido, acetato y oxigeno, diagnostico de capacidad por lote/medio y analisis de incertidumbre por FIM, multistart y perfiles de verosimilitud. Los resultados permiten cerrar documentalmente la calibracion tecnica, dejando explicitado que ciertos parametros deben mantenerse fijos, regularizados o tratados como sensibilidad antes de declarar identificabilidad parametrica completa.
