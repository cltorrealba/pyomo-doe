# Validacion del modelo con datos de vendimia 2025 (Actividad 2.14)

## Identificacion

- ID anexo: A02
- Objetivo especifico: OE2
- Actividad unica: 2.14
- Nombre de archivo sugerido: `PI-4497_IA5_A02_Act-2.14_Validacion-modelo-vendimia-2025.pdf`
- Fecha de cierre documental del bundle: 2026-06-23

## Tabla fuente del anexo

| ID anexo   | OE   |   Actividad unica | Titulo propuesto                                                  | Nombre de archivo sugerido                                   | Contenido minimo                                                                                                                      | Evidencia fuente                                                                                                                 |
|:-----------|:-----|------------------:|:------------------------------------------------------------------|:-------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------|
| A02        | OE2  |              2.14 | Validacion del modelo con datos de vendimia 2025 (Actividad 2.14) | PI-4497_IA5_A02_Act-2.14_Validacion-modelo-vendimia-2025.pdf | Documentar la base independiente 2025, protocolo de validacion, comparacion predicho-observado, metricas, brechas y decision tecnica. | Base independiente vendimia 2025, predicho vs. observado, metricas de error, criterio de validacion, brechas y decision tecnica. |

## Resumen ejecutivo

Se documento la validacion del modelo con una base independiente de vendimia 2025 compuesta por fermentaciones en mosto natural y mosto sintetico. La validacion integra carga y homologacion de datos, comparacion predicho-observado por lote y estado, metricas de error, revision de factibilidad de predicciones, analisis de sensibilidad a theta, brechas tecnicas y decision de uso.

La decision tecnica propuesta es `validacion_tecnica_condicionada`: el modelo queda respaldado para planificacion y diseno secuencial, condicionado a ejecutar un primer bloque mixto, incluir una prueba natural y recalibrar antes de cerrar la campana completa.

## Matriz de cumplimiento

| contenido_minimo                 | estado_bundle   | archivos_bundle                                                                                                             | detalle                                                                                                           |
|:---------------------------------|:----------------|:----------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------|
| Base independiente vendimia 2025 | cubierto        | data_sources/*.xlsx; tables/03_base_independiente_2025_resumen_lotes.csv; tables/03_base_independiente_2025_diccionario.csv | Incluye mosto natural y sintetico 2025, normalizacion long y resumen por lote.                                    |
| Protocolo de validacion          | cubierto        | 03_protocolo_validacion_A02.md; code_sources/run_new_must_curve_validation.py                                               | Define fuentes, estados, comparacion base vs multistart, predicho-observado, factibilidad y decision.             |
| Comparacion predicho-observado   | cubierto        | tables/04_predicho_observado_metricas.csv; figures/curve_fit__*.png                                                         | Incluye metricas por lote/estado y graficos fitcmp para lotes naturales y sinteticos.                             |
| Metricas de error                | cubierto        | tables/04_metricas_error_por_medio_estado.csv; tables/04_residuales_agregados_base_vs_best.csv                              | Incluye weighted RMSE, RMSE, sesgo medio y comparacion theta base vs mejor multistart.                            |
| Criterio de validacion           | cubierto        | tables/05_criterios_validacion.csv; 04_decision_tecnica_A02.md                                                              | Criterios de base cargada, cobertura de metricas, robustez, factibilidad y transferencia natural.                 |
| Brechas                          | cubierto        | tables/06_brechas_y_acciones.csv; source_reports/medium_transfer__medium_transfer_diagnostics_report.md                     | Incluye identificabilidad practica, transferencia de medio, redundancias y necesidad de recalibracion secuencial. |
| Decision tecnica                 | cubierto        | 04_decision_tecnica_A02.md; ANEXO_A02_borrador.md                                                                           | Decision de validacion tecnica condicionada para uso secuencial y primer bloque mixto.                            |

## Base independiente vendimia 2025

Resumen por medio:

| medium    |   n_batches |   n_rows |   t_min_h |   t_max_h |   G_initial_mean |   F_initial_mean |   YAN_initial_mean |   E_final_mean |   total_N_pulse_mean |
|:----------|------------:|---------:|----------:|----------:|-----------------:|-----------------:|-------------------:|---------------:|---------------------:|
| natural   |           9 |      112 |         0 |     211   |          76.9956 |          77.1433 |            230.093 |        81.3443 |              80      |
| synthetic |          10 |       96 |         0 |     213.5 |         124.064  |         117.962  |            221.8   |       102.12   |             100      |
| all       |          19 |      208 |         0 |     213.5 |         101.768  |          98.6268 |            225.728 |        92.2788 |              90.5263 |

Cobertura observacional por medio y estado:

| medium    | state          |   n_observations |   n_batches_with_observation |
|:----------|:---------------|-----------------:|-----------------------------:|
| natural   | X_viable_kg_m3 |              100 |                            9 |
| natural   | X_dead_kg_m3   |              100 |                            9 |
| natural   | N_kg_m3        |               91 |                            9 |
| natural   | G_g_l          |               91 |                            9 |
| natural   | F_g_l          |               91 |                            9 |
| natural   | E_g_l          |               44 |                            9 |
| natural   | glycerol_g_l   |               91 |                            9 |
| synthetic | X_viable_kg_m3 |               76 |                           10 |
| synthetic | X_dead_kg_m3   |               86 |                           10 |
| synthetic | N_kg_m3        |               83 |                           10 |
| synthetic | G_g_l          |               86 |                           10 |
| synthetic | F_g_l          |               86 |                           10 |
| synthetic | E_g_l          |               40 |                           10 |
| synthetic | glycerol_g_l   |               85 |                           10 |

## Protocolo de validacion

El protocolo completo esta en `03_protocolo_validacion_A02.md`. En terminos operativos, el flujo fue:

1. Cargar y homologar bases natural/sintetica 2025.
2. Simular curvas con theta base y theta best-multistart.
3. Comparar predicho-observado por lote/estado.
4. Consolidar RMSE/weighted RMSE por medio y estado.
5. Revisar factibilidad de predicciones DOE.
6. Documentar sensibilidad, redundancia, brechas y decision tecnica.

## Comparacion predicho-observado y metricas

Metricas agregadas por medio/estado:

| medium    | state   |   n |   base_weighted_rmse |   best_weighted_rmse |   best_minus_base_weighted_rmse |
|:----------|:--------|----:|---------------------:|---------------------:|--------------------------------:|
| natural   | E       |  44 |              6.30419 |              6.08947 |                      -0.214718  |
| natural   | F       |  91 |              5.82539 |              5.64507 |                      -0.180328  |
| natural   | G       |  91 |              5.18396 |              5.12519 |                      -0.058765  |
| natural   | Gly     |  91 |              3.52569 |              3.39701 |                      -0.128683  |
| natural   | N       |  91 |              4.28451 |              4.24807 |                      -0.0364374 |
| natural   | X       | 100 |              6.29427 |              6.01379 |                      -0.280483  |
| natural   | Xd      | 100 |              1.83397 |              1.79259 |                      -0.0413817 |
| synthetic | E       |  40 |              4.01818 |              3.26776 |                      -0.750423  |
| synthetic | F       |  86 |              3.28353 |              2.35161 |                      -0.931918  |
| synthetic | G       |  86 |              4.05645 |              3.63275 |                      -0.423695  |
| synthetic | Gly     |  85 |              2.31366 |              2.02307 |                      -0.290597  |
| synthetic | N       |  83 |              7.02814 |              7.01159 |                      -0.0165447 |
| synthetic | X       |  76 |              3.87264 |              3.68324 |                      -0.189399  |
| synthetic | Xd      |  86 |              5.62692 |              5.715   |                       0.088076  |

Los graficos `figures/curve_fit__*.png` contienen la comparacion visual predicho-observado por lote.

## Factibilidad de predicciones

Vista preliminar de predicciones DOE:

| candidate                                | theta         | medium    | family             |   final_sugar_GF |   final_E |   max_X |   max_Xd |   max_Gly |   issue_count | status   |
|:-----------------------------------------|:--------------|:----------|:-------------------|-----------------:|----------:|--------:|---------:|----------:|--------------:|:---------|
| synthetic_glucose_rich_fructose_pulse    | base          | synthetic | sugar_separation   |      0.585826    |   97.0614 | 1.68626 | 0.648745 |   9.59513 |             0 | ok       |
| synthetic_glucose_rich_fructose_pulse    | multistart_07 | synthetic | sugar_separation   |      3.22351     |  101.438  | 1.61493 | 0.616821 |   9.69677 |             0 | ok       |
| synthetic_fructose_rich_glucose_pulse    | base          | synthetic | sugar_separation   |     23.0806      |   90.1601 | 1.68626 | 0.370584 |   5.01638 |             0 | ok       |
| synthetic_fructose_rich_glucose_pulse    | multistart_07 | synthetic | sugar_separation   |      8.95528     |   90.9951 | 1.61493 | 0.401277 |   4.82837 |             0 | ok       |
| synthetic_high_biomass_low_N_maintenance | base          | synthetic | maintenance        |      0.00680225  |   99.7364 | 2.19641 | 0.893161 |   6.69284 |             0 | ok       |
| synthetic_high_biomass_low_N_maintenance | multistart_07 | synthetic | maintenance        |      3.41512e-05 |   98.7877 | 2.18534 | 0.798478 |   6.69315 |             0 | ok       |
| synthetic_late_ethanol_death_probe       | base          | synthetic | death              |      3.90744     |   98.701  | 1.29344 | 0.813566 |   5.33065 |             0 | ok       |
| synthetic_late_ethanol_death_probe       | multistart_07 | synthetic | death              |      0.156213    |  100.276  | 1.24425 | 0.805949 |   5.36216 |             0 | ok       |
| synthetic_viable_biomass_step            | base          | synthetic | biomass_input      |      7.35464e-06 |   79.0338 | 2.47844 | 0.572062 |   5.96864 |             0 | ok       |
| synthetic_viable_biomass_step            | multistart_07 | synthetic | biomass_input      |      2.76158e-08 |   79.1243 | 2.41204 | 0.573686 |   5.9997  |             0 | ok       |
| synthetic_ethanol_inhibition_challenge   | base          | synthetic | ethanol_inhibition |      5.82736     |  101.926  | 1.51168 | 1.01134  |   5.61871 |             0 | ok       |
| synthetic_ethanol_inhibition_challenge   | multistart_07 | synthetic | ethanol_inhibition |      2.38101     |  103.358  | 1.45018 | 0.902971 |   5.63524 |             0 | ok       |

## Criterios de validacion

| criterio                                | umbral_o_regla                                                                            | valor_observado                            | estado      | interpretacion                                                                                                    |
|:----------------------------------------|:------------------------------------------------------------------------------------------|:-------------------------------------------|:------------|:------------------------------------------------------------------------------------------------------------------|
| Base independiente 2025 cargada         | Debe contener lotes naturales y sinteticos 2025 homologados.                              | natural=9; synthetic=10                    | cumple      | La base independiente contiene ambos dominios de validacion.                                                      |
| Robustez theta base vs mejor multistart | Mejor multistart no debe deteriorar la mayoria de combinaciones medio/estado.             | 13/14 combinaciones mejoran o empatan      | cumple      | La alternativa multistart mejora casi todos los estados/medios; el deterioro principal reportado es Xd sintetico. |
| Factibilidad de predicciones DOE        | No activar banderas duras en etanol, biomasa, biomasa muerta, glicerol o azucar residual. | issue_count_total=0                        | cumple      | Las simulaciones seleccionadas quedan en rangos plausibles.                                                       |
| Prueba de transferencia a mosto natural | El set de validacion/diseno debe incluir al menos un candidato natural.                   | n_candidatos_naturales=2                   | cumple      | Existe prueba natural para transferencia de matriz, pero debe mantenerse en el primer bloque.                     |
| Sensibilidad a theta                    | Registrar desacuerdo base-vs-best para priorizar experimentos informativos.               | max_theta_disagreement_weighted_rmse=2.168 | advertencia | Alto desacuerdo implica informacion util, pero no cierre automatico de parametros debiles.                        |
| Redundancia de trayectorias             | Pares con distancia <=0.25 y correlacion >=0.90 requieren seleccion manual.               | pares_redundantes=2                        | advertencia | No ejecutar pares redundantes en el primer bloque; elegir el que ataque mejor la brecha.                          |
| Identificabilidad practica              | Registrar parametros no identificables y tratarlos como brecha tecnica.                   | 4/17 parametros identificables al 95%      | brecha      | La validacion permite uso secuencial, no cierre definitivo de todos los parametros.                               |

## Brechas

| brecha                                                | evidencia                                                                                                             | impacto                                             | accion_recomendada                                                                             | estado                                  |
|:------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------|:-----------------------------------------------------------------------------------------------|:----------------------------------------|
| Identificabilidad practica incompleta                 | Perfiles de verosimilitud con direcciones one-sided o no identificables.                                              | No declarar cierre parametrico definitivo.          | Usar campana secuencial: ejecutar bloque 1, recalibrar, recomputar FIM/perfiles y rerankear.   | controlada_para_validacion_condicionada |
| Transferencia de medio requiere prueba natural        | El reporte de transferencia indica que la matriz natural debe mantenerse para distinguir transferencia vs estructura. | No usar solo mosto sintetico como validacion final. | Mantener `natural_cold_hot_switch` o `natural_glucose_pulse_after_growth` en el primer bloque. | controlada_con_bloque_mixto             |
| Redundancia entre algunos candidatos                  | Pares con baja distancia RMS y alta correlacion de trayectoria.                                                       | Baja eficiencia si ambos se ejecutan temprano.      | No poner ambos en el primer bloque; escoger segun parametro debil objetivo.                    | controlada_por_protocolo                |
| Deterioro menor en Xd sintetico bajo mejor multistart | Curve validation report: synthetic Xd es la unica combinacion con deterioro leve.                                     | Xd permanece observable sensible/ruidoso.           | Revisar calidad de conteos de biomasa muerta y priorizar medidas consistentes en bloque 1.     | advertencia                             |
| Par candidato redundante detectado                    | synthetic_late_ethanol_death_probe vs synthetic_ethanol_inhibition_challenge; distancia=0.203; correlacion=0.966      | Puede aportar informacion duplicada.                | Ejecutar solo uno en bloque inicial.                                                           | advertencia                             |
| Par candidato redundante detectado                    | natural_glucose_pulse_after_growth vs natural_cold_hot_switch; distancia=0.243; correlacion=0.943                     | Puede aportar informacion duplicada.                | Ejecutar solo uno en bloque inicial.                                                           | advertencia                             |

## Decision tecnica

| decision_tecnica                | alcance                                                                                       | condiciones                                                                                                                        | no_declarar                                                                                    | uso_recomendado                                                                    |
|:--------------------------------|:----------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------|
| validacion_tecnica_condicionada | Modelo aceptable para planificacion y diseno secuencial con base independiente vendimia 2025. | Ejecutar primer bloque mixto, incluir una prueba natural, recalibrar y recomputar perfiles/FIM antes de cerrar bloques siguientes. | No declarar cierre parametrico definitivo ni validacion final de todos los parametros debiles. | Anexo de rendicion tecnica A02 y soporte para decision de validacion condicionada. |

## Texto sugerido para informe

Se valido el modelo con una base independiente de vendimia 2025 compuesta por fermentaciones naturales y sinteticas. La comparacion predicho-observado fue evaluada mediante metricas de error por lote, estado y medio, complementadas con revision de factibilidad de predicciones y robustez frente a una solucion best-multistart. Las predicciones seleccionadas no activaron banderas duras de factibilidad y el mejor multistart mejora la mayoria de combinaciones medio/estado. Se identifican brechas de identificabilidad practica, transferencia de medio y redundancia de algunos candidatos; por ello, la decision tecnica es validacion condicionada para uso secuencial, con primer bloque mixto, inclusion de prueba natural, recalibracion y reranking posterior.
