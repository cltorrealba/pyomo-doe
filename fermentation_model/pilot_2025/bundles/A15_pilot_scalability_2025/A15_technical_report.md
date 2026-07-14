# A15 - Validación de escalabilidad en bodega piloto vendimia 2025

Actividad 7.18. Bundle generado el 2026-06-29 22:18.

## 1. Resumen ejecutivo

Este anexo documenta la validación de escalabilidad desde ensayos de laboratorio con mosto sintético VL3 hacia vinificaciones piloto de vendimia 2025. El foco no es solo mostrar curvas de ajuste, sino evaluar si la estructura cinética/metabólica desarrollada en laboratorio puede transferirse a mosto natural y escala piloto, y qué extensiones fueron necesarias.

Resultado principal:

- La estructura primaria de fermentación es transferible de forma parcial a buena: biomasa, azúcar, nitrógeno, etanol y glicerol comparten una arquitectura ODE común.
- La transferencia directa desde mosto sintético no es suficiente para describir todo el piloto: se requirieron extensiones para CO2, transición de oxígeno, secundarios y aromas.
- Los datos piloto 2025 permitieron validar escalabilidad operacional y revelar brechas reales: CO2 online, condensado aromático, reinóculo 25171, y exclusión de sensores CO2 25150/25151.
- El enfoque recomendado es mixto: laboratorio sintético para separar mecanismos e identificabilidad; piloto natural para validar matriz, escala y outputs aromáticos.

Modelos piloto seleccionados:

- CO2: `solubility_o2_slow_transition`.
- Secundarios: `secondary_full_chem_o2fixed`.
- Aromas: `ea_ethanol_nlimited`.

## 2. Evidencia fuente incluida

El bundle contiene:

- `A15_technical_report.md`: informe autoexplicativo.
- `README_bundle.md`: guía del contenido del paquete.
- `tables/`: tablas curadas de batches, KPIs, selección de modelos, FIM, estimabilidad, incidencias y campañas.
- `figures/`: figuras seleccionadas para anexar.
- `source_data/`: workbooks fuente compactos usados como evidencia (`mosto_sintetico_vl3.xlsx` y `Calibration_data_vl3.xlsx`).
- `source_snapshot/`: scripts fuente mínimos para trazabilidad.

## 3. Configuración experimental y datasets

### 3.1 Laboratorio sintético VL3

Fuente: `fermentation_model/data/Laboratorio 2025/mosto_sintetico_vl3.xlsx`.

La planilla incluye 10 cubadas sintéticas MS007-MS016, hojas de resumen, datos homologados, diseño CCD, flags de calidad y diccionario de mapeo. El diseño VL3 cubre principalmente variación de temperatura y composición inicial de azúcares/nutrientes, con mosto sintético como matriz controlada.

Resumen de batches sintéticos:

| batch | n_rows | t_max_h | temperature_initial_c | temperature_min_c | temperature_max_c | G_initial_g_l | F_initial_g_l | S_initial_g_l |
| ----- | ------ | ------- | --------------------- | ----------------- | ----------------- | ------------- | ------------- | ------------- |
| MS007 | 10     | 213     | 15                    | 15                | 15                | 125.3         | 120.5         | 245.9         |
| MS008 | 10     | 213     | 15                    | 15                | 15                | 124.3         | 119.6         | 243.9         |
| MS009 | 10     | 213     | 23                    | 23                | 23                | 120.9         | 114.2         | 235.1         |
| MS010 | 10     | 213.5   | 19                    | 19                | 19                | 124.7         | 113.2         | 237.9         |
| MS011 | 10     | 213.5   | 19                    | 19                | 19                | 123.6         | 119.8         | 243.3         |
| MS012 | 10     | 213.5   | 23                    | 23                | 23                | 122.6         | 114.2         | 236.8         |
| MS013 | 9      | 213     | 19                    | 19                | 19                | 125.1         | 116.2         | 241.3         |
| MS014 | 9      | 213     | 19                    | 19                | 19                | 127.2         | 126.5         | 253.7         |
| MS015 | 9      | 213     | 15                    | 15                | 15                | 121.9         | 118.5         | 240.4         |
| MS016 | 9      | 208     | 23                    | 23                | 23                | 125.2         | 116.9         | 242.1         |

### 3.2 Piloto vendimia 2025

Fuente: `fermentation_model/data/Piloto 2025/Calibration_data_vl3.xlsx`.

La planilla piloto contiene 8 vinificaciones naturales: 25026, 25027, 25085, 25086, 25150, 25151, 25170 y 25171. Incluye temperatura, densidad, biomasa viable, peso seco, glucosa, fructosa, PAN, amonio, YAN, glicerol, piruvato, acetaldehído, etanol, pulsos de nutriente, aromas totales y aromas en condensado.

Resumen de batches piloto usados por el pipeline:

| batch | n_rows | t_max_h | temperature_min_c | temperature_max_c | S_initial_g_l | S_final_g_l | YAN_initial_mg_l | E_final_g_l | n_aroma_total_obs | n_co2_sensor_points |
| ----- | ------ | ------- | ----------------- | ----------------- | ------------- | ----------- | ---------------- | ----------- | ----------------- | ------------------- |
| 25026 | 24     | 138     | 18.5              | 22                | 147.9         | 0.52        | 154.5            | 96.49       | 40                | 0                   |
| 25027 | 32     | 186     | 17                | 20.5              | 155.4         | 1.49        | 162.5            | 96.57       | 40                | 0                   |
| 25085 | 35     | 216     | 17.1              | 19.1              | 100.4         | 3.46        | 93.34            | 94.84       | 40                | 0                   |
| 25086 | 35     | 216     | 15.3              | 28.9              | 136.3         | 5.28        | 92.7             | 95.31       | 40                | 0                   |
| 25150 | 30     | 174     | 19.5              | 21.6              | 132.6         | 1.65        | 94.24            | 94.52       | 48                | 0                   |
| 25151 | 30     | 174     | 19.7              | 22                | 136.8         | 1.07        | 90.78            | 94.13       | 48                | 0                   |
| 25170 | 63     | 372     | 14                | 20                | 140.2         | 8.97        | 98.42            | 95.39       | 88                | 191                 |
| 25171 | 44     | 258     | 13.2              | 20.8              | 120.4         | 4.33        | 116.6            | 93.89       | 64                | 110                 |

### 3.3 Brechas de configuración física

El workbook piloto no codifica completamente la configuración física de bodega piloto: volumen útil, geometría, agitación, transferencia térmica ni hardware de sensor. Por tanto, el análisis de escalabilidad se basa en variables operacionales medidas y en evidencia de proceso, no en un balance físico completo de estanque.

Esto debe completarse manualmente en la bitácora oficial del proyecto si se requiere trazabilidad de hardware.

## 4. Validación de datos y homologación

Se homologaron variables entre laboratorio y piloto:

- Tiempo de fermentación.
- Temperatura operacional.
- Glucosa y fructosa.
- Azúcar total a partir de densidad.
- Biomasa viable.
- Nitrógeno asimilable.
- Etanol.
- Glicerol.
- Pulsos de nutriente.
- Estados secundarios y aromas en piloto.

La regresión densidad-azúcar mostró alta consistencia dentro de cada medio:

| group     | source   | intercept | slope | r2     | rmse_g_l | n_points |
| --------- | -------- | --------- | ----- | ------ | -------- | -------- |
| natural   | G_plus_F | -2119     | 2.128 | 0.9963 | 3.539    | 88       |
| natural   | Y15      | -2119     | 2.128 | 0.9963 | 3.539    | 88       |
| synthetic | G_plus_F | -2169     | 2.19  | 0.9937 | 7.161    | 86       |
| synthetic | Y15      | -2169     | 2.19  | 0.9937 | 7.161    | 86       |

Figura clave: `figures/fig_03_density_sugar_regression_lab.png`.

Flags/incidencias de laboratorio:

| Codigo  | Cubada | Horas | Columna | Valor | Criterio  | Accion                                                                        |
| ------- | ------ | ----- | ------- | ----- | --------- | ----------------------------------------------------------------------------- |
| MS009-2 | MS009  | 22.5  | Brix    | 45797 | Brix > 40 | Valor crudo preservado; revisar si corresponde a error de carga/fecha serial. |

Incidencias piloto:

| item                        | type                  | description                                                                      | model_action                                                                     |
| --------------------------- | --------------------- | -------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| CO2 25150/25151             | sensor_data_exclusion | Los archivos CO2 de 25150 y 25151 fueron marcados como no usables para calibr... | Excluir del ajuste CO2; mantener datos offline de composición.                   |
| 25171                       | process_incident      | Inóculo inicial con baja/no viabilidad; el proceso útil se reancla al punto p... | Usar el punto ING25-SB012-Pre reinóculo como t=0 efectivo.                       |
| Configuración física piloto | metadata_gap          | El workbook no codifica volumen de estanque, geometría, agitación ni detalle ... | Reportar como brecha documental; el análisis usa variables operacionales medi... |

## 5. Modelo base transferido desde laboratorio

El modelo primario transferido usa:

\[
x=[X,X_d,N,G,F,E,Gly]^T
\]

\[
\frac{dX}{dt}=(\mu-k_d)X+u_X
\]

\[
\frac{dX_d}{dt}=k_dX
\]

\[
\frac{dN}{dt}=-q_N f_N(T,N)X+u_N
\]

\[
\frac{dG}{dt}=-\left(q_{XG}f_N+q_{EG}f_G+m\frac{G}{G+F}\right)X+u_G
\]

\[
\frac{dF}{dt}=-\left(q_{XF}f_N+q_{EF}f_F+m\frac{F}{G+F}\right)X+u_F
\]

\[
\frac{dE}{dt}=(\beta_G f_G+\beta_F f_F)X+u_E
\]

\[
\frac{dGly}{dt}=(\gamma_G f_G+\gamma_F f_F)X
\]

Los factores cinéticos incluyen dependencia de temperatura, saturación por sustrato, inhibición por etanol e interacción glucosa/fructosa. Esta estructura fue primero evaluada en laboratorio y después reutilizada como base para piloto.

## 6. Calibración laboratorio y evaluación medio sintético/natural

La calibración de laboratorio comparó fits por medio y un fit mixto:

| fit                 | n_batches | mediums           | n_parameters | success | nfev | final_wsse | wsse_per_residual | l2_lambda |
| ------------------- | --------- | ----------------- | ------------ | ------- | ---- | ---------- | ----------------- | --------- |
| natural_reduced11   | 9         | natural           | 11           | 1       | 12   | 1.341e+04  | 22.06             | 0         |
| synthetic_reduced11 | 10        | synthetic         | 11           | 1       | 16   | 1.072e+04  | 19.79             | 0         |
| mixed_reduced11     | 19        | natural,synthetic | 11           | 1       | 16   | 2.714e+04  | 23.6              | 0         |
| mixed_full17_l2     | 19        | natural,synthetic | 17           | 1       | 9    | 2.613e+04  | 22.53             | 1         |

Interpretación:

- El ajuste sintético reducido tuvo menor WSSE/residual que el natural reducido.
- El ajuste mixto completo con L2 permitió una estructura común con más parámetros sin perder estabilidad.
- La matriz natural introduce desviaciones que no aparecen en sintético puro.
- Por tanto, el mosto sintético es útil para identificación mecanística, pero no debe ser el único prior para piloto.

Figuras:

- `figures/fig_05_lab_transfer_fit_metrics.png`
- `figures/fig_13_lab_synthetic_fit_MS007.png`
- `figures/fig_14_lab_synthetic_fit_MS009.png`

## 7. FIM e identificabilidad en laboratorio

La FIM se calculó en espacio log-paramétrico:

\[
J_{ij}\approx\frac{r_i(\theta_j e^h)-r_i(\theta_j e^{-h})}{2h}
\]

\[
F=J^TJ
\]

Resumen FIM:

| analysis                   | n_residuals | logdet | min_eigenvalue | condition_number | trace_inv | rank_1e-8 |
| -------------------------- | ----------- | ------ | -------------- | ---------------- | --------- | --------- |
| natural_current            | 608         | 73.82  | 0.03686        | 2.224e+06        | 35.08     | 17        |
| synthetic_current          | 542         | 75.26  | 0.002178       | 2.893e+07        | 474.2     | 17        |
| mixed_current_reducedprior | 1150        | 85.96  | 0.004577       | 2.868e+07        | 225.9     | 17        |
| mixed_current_full_l2      | 1150        | 91.23  | 0.02552        | 5.195e+06        | 40.69     | 17        |

Estimabilidad del modelo mixto laboratorio:

| classification     | count |
| ------------------ | ----- |
| well_estimated     | 13    |
| weak_or_confounded | 2     |
| moderate           | 2     |

Conclusión de laboratorio:

- La mezcla sintético/natural mejora el contenido de información frente a cada medio por separado.
- Persisten direcciones débiles: principalmente saturación de nitrógeno, mantenimiento y algunos términos acoplados a fructosa.
- La campaña MBDoE de laboratorio prioriza diseños sintéticos porque permiten perturbar composición inicial y separar mecanismos que en mosto natural están correlacionados.

## 8. Validación piloto y extensiones requeridas

Al transferir el modelo a piloto, se mantuvo el núcleo primario pero se extendió el simulador:

1. Capa secundaria: piruvato, acetaldehído, acetato y O2 latente.
2. Aromas: ethyl acetate, isoamyl acetate y ethyl octanoate.
3. Volatilización: partición gas-líquido Antoine+UNIFAC y stripping por flujo CO2.
4. CO2: estado de CO2 disuelto y transición macroscópica de oxígeno.

Benchmark de modelos secundarios piloto:

| model                          | data_wsse | n_parameters | bic  | active_bound_count | state_penalty | selection_score |
| ------------------------------ | --------- | ------------ | ---- | ------------------ | ------------- | --------------- |
| secondary_full_chem_o2fixed    | 1817      | 12           | 1882 | 1                  | 25            | 1907            |
| secondary_redox_o2             | 1882      | 9            | 1931 | 2                  | 50            | 1981            |
| secondary_phase_split          | 2588      | 9            | 2637 | 1                  | 25            | 2662            |
| secondary_reduced_o2fixed      | 2799      | 7            | 2837 | 1                  | 25            | 2862            |
| secondary_acetate_assimilation | 2799      | 8            | 2842 | 1                  | 25            | 2867            |

Benchmark aromático piloto:

| model                    | n_parameters | data_wsse | bic  | active_bound_count | ea_retained_relative_rmse |
| ------------------------ | ------------ | --------- | ---- | ------------------ | ------------------------- |
| ea_ethanol_nlimited      | 11           | 7370      | 7431 | 2                  | 0.576                     |
| ea_redox_acetaldehyde    | 11           | 7493      | 7555 | 2                  | 0.6099                    |
| ea_combined_parsimonious | 12           | 7490      | 7557 | 3                  | 0.6102                    |
| ea_ethanol_biomass       | 10           | 7697      | 7753 | 2                  | 0.5569                    |
| ea_ethanol_temperature   | 11           | 7697      | 7758 | 3                  | 0.5569                    |
| ea_biomass_background    | 10           | 8250      | 8305 | 2                  | 0.5366                    |
| baseline_phase           | 9            | 9449      | 9499 | 2                  | 0.5136                    |

Benchmark CO2 piloto:

| mode                          | mechanistic_class      | data_wsse | bic   | active_bound_count | selection_penalty | selection_score | selected |
| ----------------------------- | ---------------------- | --------- | ----- | ------------------ | ----------------- | --------------- | -------- |
| solubility_o2_slow_transition | o2_gated_dissolved_co2 | 124.7     | 136.1 | 0                  | 0                 | 136.1           | 1        |
| solubility_o2_qfit            | o2_gated_dissolved_co2 | 124.7     | 141.8 | 1                  | 35                | 176.8           | 0        |
| old_lag_threshold             | empirical_effective    | 163.4     | 174.9 | 1                  | 43                | 217.9           | 0        |
| solubility_o2_literature      | o2_gated_dissolved_co2 | 232.3     | 243.7 | 0                  | 0                 | 243.7           | 0        |
| solubility_scaled             | dissolved_co2          | 241.5     | 252.9 | 0                  | 0                 | 252.9           | 0        |
| solubility_fixed              | dissolved_co2          | 254.6     | 260.3 | 0                  | 0                 | 260.3           | 0        |
| instant                       | rate_proxy             | 240.8     | 246.5 | 0                  | 25                | 271.5           | 0        |

## 9. KPIs y resultados piloto

Métricas finales del modelo calibrado piloto:

| group     | state          | species         | n   | relative_rmse | relative_bias | corr    |
| --------- | -------------- | --------------- | --- | ------------- | ------------- | ------- |
| core      | E              |                 | 129 | 0.1781        | 0.01296       | 0.9092  |
| core      | F              |                 | 129 | 0.5689        | -0.2664       | 0.87    |
| core      | G              |                 | 129 | 0.5388        | 0.03342       | 0.9049  |
| core      | Gly            |                 | 129 | 0.1084        | -0.000422     | 0.9103  |
| core      | N              |                 | 125 | 0.9442        | -0.2906       | 0.6739  |
| core      | X              |                 | 122 | 0.338         | -0.09685      | 0.6451  |
| secondary | AcAld          |                 | 105 | 0.1128        | -0.02193      | 0.6353  |
| secondary | Pyr            |                 | 127 | 0.3309        | -0.01179      | 0.4308  |
| aroma     |                | ethyl_acetate   | 17  | 2.483         | 1.714         | 0.5031  |
| aroma     |                | ethyl_acetate   | 17  | 0.5692        | -0.4287       | 0.02023 |
| aroma     |                | ethyl_acetate   | 51  | 0.4552        | -0.2648       | 0.627   |
| aroma     |                | ethyl_octanoate | 44  | 1.047         | -0.7051       | 0.4994  |
| aroma     |                | ethyl_octanoate | 44  | 0.6279        | -0.01445      | 0.3338  |
| aroma     |                | ethyl_octanoate | 51  | 0.5605        | -0.1897       | 0.627   |
| aroma     |                | isoamyl_acetate | 44  | 0.9022        | -0.5369       | 0.4696  |
| aroma     |                | isoamyl_acetate | 44  | 0.3523        | -0.1105       | 0.6426  |
| aroma     |                | isoamyl_acetate | 51  | 0.39          | -0.1537       | 0.7746  |
| co2       | CO2_flow_L_min |                 | 191 | 0.779         | -0.4388       | 0.8429  |
| co2       | CO2_flow_L_min |                 | 110 | 0.435         | -0.02541      | 0.9312  |

Figuras:

- `figures/fig_06_pilot_fit_metrics.png`
- `figures/fig_15_pilot_final_fit_25170.png`
- `figures/fig_16_pilot_final_fit_25171.png`
- `figures/fig_17_pilot_co2_benchmark_25170.png`
- `figures/fig_18_pilot_co2_benchmark_25171.png`

Interpretación de escalabilidad:

- Glicerol y etanol muestran transferencia razonable de estructura primaria.
- Azúcares y nitrógeno son más sensibles a matriz natural, medición y eventos operacionales.
- CO2 online no es explicable por la cinética primaria sola: requiere solubilidad, transición de O2 y liberación gas-líquido.
- Aromas requieren estructura propia y datos piloto; no son extrapolables desde mosto sintético sin medición aromática.

## 10. Diseño experimental y escalabilidad

Campaña MBDoE laboratorio:

| campaign_order | candidate                                | family              | medium    | temperature_segments | score | rationale                                                                        |
| -------------- | ---------------------------------------- | ------------------- | --------- | -------------------- | ----- | -------------------------------------------------------------------------------- |
| 1              | synthetic_glucose_rich_fructose_pulse    | sugar_separation    | synthetic | 17, 20, 23, 20       | 66.69 | High glucose plus fructose pulse separates glucose and fructose uptake/yield ... |
| 2              | synthetic_fructose_rich_glucose_pulse    | sugar_separation    | synthetic | 17, 20, 23, 20       | 70.34 | Fructose-rich run probes iG and fructose kinetic directions.                     |
| 3              | synthetic_high_biomass_low_N_maintenance | maintenance         | synthetic | 16, 18, 20, 18       | 73    | High biomass with low N creates low-growth sugar consumption windows for m0.     |
| 4              | synthetic_late_ethanol_death_probe       | death               | synthetic | 20, 25, 25, 22       | 74.54 | Late ethanol stress plus Xd observation targets Kd0/iE.                          |
| 5              | synthetic_viable_biomass_step            | biomass_input       | synthetic | 18, 22, 22, 18       | 75.59 | Known viable biomass addition tests rate proportionality to X and supports yi... |
| 6              | synthetic_ethanol_inhibition_challenge   | ethanol_inhibition  | synthetic | 18, 22, 25, 22       | 76.52 | Initial ethanol decouples ethanol inhibition from ethanol generated by fermen... |
| 7              | synthetic_high_sugar_reference           | synthetic_baseline  | synthetic | 18, 20, 20, 18       | 77.39 | Synthetic high-sugar operating point close to the current synthetic dataset.     |
| 8              | synthetic_low_yan_ladder                 | nitrogen_saturation | synthetic | 15, 18, 22, 22       | 78.18 | Low-YAN ladder targets sN and qN separation.                                     |
| 9              | natural_glucose_pulse_after_growth       | natural_sugar_pulse | natural   | 18, 20, 22, 19       | 78.8  | Natural must with glucose perturbation to test sugar-transfer validity.          |

Campaña MBDoE piloto natural:

| campaign_order | candidate                              | family                 | medium  | temperature_segments | N_pulses_kg_m3       | score |
| -------------- | -------------------------------------- | ---------------------- | ------- | -------------------- | -------------------- | ----- |
| 1              | natural_pilot_cold_to_warm_earlyN      | temperature_N          | natural | 14, 16, 22, 20       | 30h:0.045            | 72.97 |
| 2              | natural_pilot_EA_warm_early_noN        | EA_temperature_Nstress | natural | 24, 23, 19, 17       |                      | 82.94 |
| 3              | natural_pilot_EA_cold_retention_noN    | EA_retention_reference | natural | 13, 13, 16, 18       |                      | 89.43 |
| 4              | natural_pilot_EA_cold_warm_Nsplit      | EA_temperature_Nsplit  | natural | 13, 18, 23, 19       | 32h:0.025; 80h:0.035 | 94.5  |
| 5              | natural_pilot_warm_to_cool_noN         | temperature            | natural | 22, 22, 17, 16       |                      | 98.38 |
| 6              | natural_pilot_low_temp_aroma_retention | aroma_retention        | natural | 13, 15, 16, 16       | 54h:0.035            | 101.6 |

Lectura técnica:

- Laboratorio sintético maximiza separación de mecanismos por manipulación de composición inicial, azúcar, biomasa y etanol.
- Piloto natural prioriza temperatura, N y condiciones compatibles con mosto real.
- La escalabilidad se valida combinando ambas escalas: sintético para identificabilidad y piloto para validez de matriz/proceso.

## 11. Matriz de decisión de transferibilidad

| layer            | lab_synthetic_evidence                                                          | pilot_evidence                                                                   | transferability          | action                                                                     |
| ---------------- | ------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- | ------------------------ | -------------------------------------------------------------------------- |
| primary_kinetics | CCD synthetic VL3 excites sugar, temperature and nitrogen directions.           | Pilot natural must reproduces primary trends but with matrix and scale deviat... | partial_to_good          | Use shared structure; keep recalibrated pilot priors for MPCC.             |
| glycerol         | Glycerol state added and fitted in mixed lab dataset.                           | Pilot final relative RMSE for Gly is low compared with other primary states.     | good                     | Retain glycerol state in scale-up model.                                   |
| nitrogen         | Synthetic low/high YAN informs qN and sN but sN remains weak.                   | YAN residuals remain high and natural must has remanent YAN behavior.            | limited                  | Use stronger priors; consider ammonium/PAN split in future.                |
| CO2              | Not present in synthetic VL3 model stage.                                       | Online pilot CO2 required dissolved CO2 and O2-transition structure.             | requires_pilot_extension | Use pilot-specific CO2 model for MPCC constraints.                         |
| aromas           | Not directly calibrated in synthetic VL3 stage.                                 | Pilot aroma/condensate data enabled EA/IAA/EO synthesis-loss layer.              | requires_pilot_extension | Use pilot aroma structure; keep ethyl acetate uncertainty explicit.        |
| DOE              | Synthetic designs ranked high for separating sugar and weak kinetic directions. | Pilot designs favor natural-must temperature and N perturbations.                | complementary            | Use synthetic for mechanism separation; pilot for scale/matrix validation. |

## 12. Conclusiones

1. La estructura primaria desarrollada con mosto sintético VL3 es útil y parcialmente transferible a escala piloto.
2. La transferencia directa no es suficiente para los outputs de bodega piloto: CO2, O2, aromas y condensado requieren extensiones fenomenológicas.
3. La validación piloto 2025 entregó evidencia real de escalabilidad operacional: múltiples vinificaciones, tratamientos térmicos/nutricionales, datos aromáticos, CO2 online y eventos operativos.
4. El modelo escalado debe usarse con parámetros piloto recalibrados y no con parámetros sintéticos puros.
5. La campaña MBDoE sugiere continuar con un approach mixto: sintético para excitar direcciones poco identificables y piloto natural para validar transferencia.

## 13. Pasos a seguir

- Completar bitácora física de bodega piloto: volumen, geometría, sensores, configuración térmica y protocolo de muestreo.
- Ejecutar nuevos ensayos piloto con CO2 online confiable y condensado final.
- Recalibrar el modelo integrado tras los ensayos MBDoE.
- Separar YAN en amonio/PAN si se busca explicar remanentes de nitrógeno.
- Definir qué parámetros quedan libres en MPCC y cuáles se fijan por robustez.
