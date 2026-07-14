# Resultados e interpretacion A06

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

## Metricas de calibracion

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

## Interpretacion tecnica

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
