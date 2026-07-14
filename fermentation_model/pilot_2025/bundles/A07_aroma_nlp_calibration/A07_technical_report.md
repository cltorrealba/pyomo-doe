# A07 - Optimización NLP para calibración de síntesis aromática

Actividad 3.14. Bundle generado el 2026-06-29 18:46.

## 1. Resumen ejecutivo

Se construyó y documentó un pipeline determinista para calibrar y validar un modelo dinámico de fermentación alcohólica con síntesis y pérdida de aromas. El objetivo técnico fue transformar datos piloto de mosto natural, sensores de CO2 y mediciones aromáticas en un modelo ODE calibrable que pueda transferirse posteriormente como restricción dinámica dentro del MPCC.

La estructura vigente seleccionada es:

- Fermentación primaria con biomasa viable/muerta, nitrógeno asimilable, glucosa, fructosa, etanol y glicerol.
- Capa secundaria con piruvato, acetaldehído, acetato y oxígeno latente.
- Aromas: ethyl acetate, isoamyl acetate y ethyl octanoate, separados en retenido líquido y fracción volatilizada/condensada.
- CO2: modelo seleccionado `solubility_o2_slow_transition`, con transición macroscópica de O2 y acumulación/liberación de CO2 disuelto.

Modelos seleccionados:

- Aromas: `ea_ethanol_nlimited`.
- Secundarios: `secondary_full_chem_o2fixed`.
- CO2: `solubility_o2_slow_transition`.

El benchmark de CO2 seleccionó `solubility_o2_slow_transition` con score `136.0676470728151` y WSSE `124.65342654331737`. La selección aromática favoreció una estructura donde ethyl acetate depende de etanol, biomasa y limitación de nitrógeno, no solo de una tasa de consumo de azúcar.

## 2. Evidencia fuente incluida

El bundle contiene:

- `A07_technical_report.md`: este informe técnico.
- `tables/`: tablas curadas con selección de modelos, parámetros, métricas de ajuste, estimabilidad y campaña MBDoE.
- `figures/`: figuras seleccionadas para anexar directamente.
- `source_snapshot/`: copia compacta de scripts/notebook fuente que implementan el simulador, calibración, benchmarking, FIM y MBDoE.

No se incluyen todos los CSV o figuras intermedias para evitar un paquete excesivo. Los resultados completos siguen en `fermentation_model/pilot_2025/results/`.

## 3. Datos y validación de lectura

Se procesaron planillas piloto de mosto natural y sensores CO2 online. La curación relevante fue:

- CO2 de 25150 y 25151 se excluye.
- 25170 se usa como set de CO2 más confiable.
- 25171 se reancla a un nuevo tiempo cero en el punto post-reinóculo.
- Aromas `xxx_total` se interpretan como retenido líquido más condensado equivalente.
- Aromas `xxx_condensado` se interpretan como fracción volatilizada acumulada, expresada en concentración equivalente de mosto.

La lectura se verificó gráficamente en `figures/fig_11_*` y `figures/fig_12_*`.

## 4. Modelo dinámico de simulación

### 4.1 Estados primarios

El núcleo de fermentación usa:

\[
x = [X, X_d, N, G, F, E, Gly]^T
\]

donde `X` es biomasa viable, `X_d` biomasa muerta, `N` nitrógeno asimilable, `G` glucosa, `F` fructosa, `E` etanol y `Gly` glicerol.

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

Los factores `f_N`, `f_G` y `f_F` combinan limitación tipo Monod, corrección térmica Arrhenius, inhibición por etanol e inhibición de fructosa por glucosa. Los pulsos operacionales entran como funciones gaussianas suavizadas:

\[
u_j(t)=\sum_k \frac{a_{j,k}}{\sqrt{\pi}w}\exp\left[-\left(\frac{t-t_k}{w}\right)^2\right]
\]

### 4.2 Estados secundarios

La capa secundaria seleccionada (`secondary_full_chem_o2fixed`) representa especies que explican desviaciones de aromas y metabolismo lateral:

\[
z=[Pyr, AcAld, Acetate, O_2]^T
\]

La estructura usa términos de producción por consumo de azúcar, conversión piruvato-acetaldehído, reducción/oxidación de acetaldehído, generación/asimilación de acetato y efecto latente de O2. El benchmark de estructuras secundarias se resume en `figures/fig_02_secondary_model_selection.png` y `tables/02_secondary_model_selection.csv`.

### 4.3 CO2 disuelto y flujo gaseoso

El CO2 producido se vincula estequiométricamente a producción de etanol:

\[
r_{CO2,base}=\frac{44.01}{2\cdot 46.07}r_E
\]

En el modelo seleccionado, la producción efectiva se modula por una transición macroscópica de O2:

\[
r_{CO2,eff}=
r_{CO2,base}\left[f_C+(1-f_C)\phi_{ana}(O_2)\right]+r_{CO2,resp}
\]

\[
\phi_{ana}(O_2)=\frac{K_{ana}^n}{K_{ana}^n+O_2^n}
\]

El CO2 disuelto se representa con:

\[
\frac{dC_{CO2,L}}{dt}=r_{CO2,eff}-r_{CO2,gas}
\]

\[
r_{CO2,gas}=k_{rel}\max(C_{CO2,L}-C^*_{CO2},0)
\]

La saturación efectiva:

\[
C^*_{CO2}=
s_{CO2}1.69\exp[-0.032(T-20)]\exp(0.0016E)\exp[-0.0012(G+F)]
\]

Parámetros CO2 seleccionados:

| parameter         | value  |
| ----------------- | ------ |
| nan               | 0      |
| kCO2_release_h    | 0.2669 |
| CO2sat_scale      | 0.4066 |
| O2_qmax_mg_gdw_h  | 0.15   |
| O2_K_mg_l         | 0.25   |
| O2_ana_K_mg_l     | 0.75   |
| O2_ana_hill       | 2      |
| O2_crabtree_floor | 0.08   |
| O2sat_scale       | 1      |
| O2_kLa_h          | 0      |
| O2_25171_fraction | 0.05   |

### 4.4 Síntesis y volatilización aromática

Para cada aroma `i`:

\[
\frac{dA_i^{liq}}{dt}=r_i^{prod}-r_i^{loss}
\]

\[
\frac{dA_i^{cond}}{dt}=r_i^{loss}
\]

\[
A_i^{total}=A_i^{liq}+A_i^{cond}
\]

La producción base separa fase de crecimiento y fase estacionaria:

\[
r_i^{phase}=
\left(k_{i,g}\phi_N+k_{i,s}(1-\phi_N)\right)q_S
\]

\[
\phi_N=\frac{N}{N+K_N}
\]

Para ethyl acetate, el modelo seleccionado agrega una fuente dependiente de etanol/biomasa y limitación de nitrógeno:

\[
r_{EA}=
r_{EA}^{phase}
+k_{EA,XE}X\frac{E}{K_E+E}
+k_{EA,XE,Nlim}X\frac{E}{K_E+E}\frac{K_N}{K_N+N}
\]

La pérdida aromática no se trata como un `k_vol` libre genérico. Se usa:

\[
r_i^{loss}=
\alpha_{i,loss}K_i^{LG}(T,E,G,F)q_{gas}A_i^{liq}
\]

donde `K_i^{LG}` proviene de presión de vapor Antoine y coeficientes de actividad UNIFAC con mezcla agua-etanol-azúcar equivalente.

## 5. Formulación NLP de calibración

La calibración se formula como un problema de mínimos cuadrados no lineal con restricciones de caja:

\[
\min_{\theta} \sum_{j=1}^M \rho\left(
\frac{\hat y_j(\theta)-y_j}{\sigma_j}
\right)^2
\]

sujeto a:

\[
\dot x=f(x,z,A,C_{CO2},u,T,\theta)
\]

\[
\theta_L \leq \theta \leq \theta_U
\]

\[
x(t)\geq 0,\quad z(t)\geq 0,\quad A(t)\geq 0
\]

La optimización se ejecutó en espacio logarítmico para parámetros positivos:

\[
\eta=\log(\theta)
\]

El solver numérico fue `scipy.optimize.least_squares` con:

- método `trf` (trust-region reflective),
- restricciones de caja,
- pérdida robusta `soft_l1`,
- `f_scale=2.0`,
- multistart para modelos estructurales y CO2,
- penalización adicional por parámetros en borde durante selección de modelo.

La selección estructural usó:

\[
AICc = n\log(WSSE/n)+2p+\frac{2p(p+1)}{n-p-1}
\]

\[
BIC = n\log(WSSE/n)+p\log(n)
\]

\[
Score = BIC + P_{bounds} + P_{adequacy}
\]

## 6. Benchmarking de modelos

### 6.1 Aroma

La comparación de modelos de ethyl acetate mostró que la estructura `ea_ethanol_nlimited` entrega el mejor compromiso entre ajuste y complejidad. La evidencia principal está en `figures/fig_01_aroma_model_selection.png`.

| model                    | n_parameters | data_wsse | bic  | active_bound_count | ea_retained_relative_rmse |
| ------------------------ | ------------ | --------- | ---- | ------------------ | ------------------------- |
| ea_ethanol_nlimited      | 11           | 7370      | 7431 | 2                  | 0.576                     |
| ea_redox_acetaldehyde    | 11           | 7493      | 7555 | 2                  | 0.6099                    |
| ea_combined_parsimonious | 12           | 7490      | 7557 | 3                  | 0.6102                    |
| ea_ethanol_biomass       | 10           | 7697      | 7753 | 2                  | 0.5569                    |
| ea_ethanol_temperature   | 11           | 7697      | 7758 | 3                  | 0.5569                    |
| ea_biomass_background    | 10           | 8250      | 8305 | 2                  | 0.5366                    |
| baseline_phase           | 9            | 9449      | 9499 | 2                  | 0.5136                    |

### 6.2 Secundarios

| model                          | data_wsse | n_parameters | bic  | active_bound_count | state_penalty | selection_score |
| ------------------------------ | --------- | ------------ | ---- | ------------------ | ------------- | --------------- |
| secondary_full_chem_o2fixed    | 1817      | 12           | 1882 | 1                  | 25            | 1907            |
| secondary_redox_o2             | 1882      | 9            | 1931 | 2                  | 50            | 1981            |
| secondary_phase_split          | 2588      | 9            | 2637 | 1                  | 25            | 2662            |
| secondary_reduced_o2fixed      | 2799      | 7            | 2837 | 1                  | 25            | 2862            |
| secondary_acetate_assimilation | 2799      | 8            | 2842 | 1                  | 25            | 2867            |

### 6.3 CO2

| mode                          | mechanistic_class      | data_wsse | bic   | active_bound_count | selection_penalty | selection_score | selected |
| ----------------------------- | ---------------------- | --------- | ----- | ------------------ | ----------------- | --------------- | -------- |
| solubility_o2_slow_transition | o2_gated_dissolved_co2 | 124.7     | 136.1 | 0                  | 0                 | 136.1           | 1        |
| solubility_o2_qfit            | o2_gated_dissolved_co2 | 124.7     | 141.8 | 1                  | 35                | 176.8           | 0        |
| old_lag_threshold             | empirical_effective    | 163.4     | 174.9 | 1                  | 43                | 217.9           | 0        |
| solubility_o2_literature      | o2_gated_dissolved_co2 | 232.3     | 243.7 | 0                  | 0                 | 243.7           | 0        |
| solubility_scaled             | dissolved_co2          | 241.5     | 252.9 | 0                  | 0                 | 252.9           | 0        |
| solubility_fixed              | dissolved_co2          | 254.6     | 260.3 | 0                  | 0                 | 260.3           | 0        |
| instant                       | rate_proxy             | 240.8     | 246.5 | 0                  | 25                | 271.5           | 0        |

La estructura `old_lag_threshold` ajusta razonablemente, pero queda penalizada por ser empírica y por tocar borde de parámetro. El modelo `solubility_o2_slow_transition` fue seleccionado porque mejora la interpretación fenomenológica del retraso: transición O2 + acumulación de CO2 disuelto + liberación gaseosa.

## 7. Convergencia y robustez numérica

### CO2

| mode                          | start | success | status | nfev | data_wsse | active_bound_count |
| ----------------------------- | ----- | ------- | ------ | ---- | --------- | ------------------ |
| solubility_o2_slow_transition | 0     | 1       | 2      | 16   | 124.7     | 0                  |
| solubility_o2_slow_transition | 1     | 1       | 3      | 2    | 511.5     | 0                  |

### Ajuste integrado de aromas

| mode                          | start | success | status | nfev | initial_objective_wsse | final_objective_wsse | active_bound_count |
| ----------------------------- | ----- | ------- | ------ | ---- | ---------------------- | -------------------- | ------------------ |
| solubility_o2_slow_transition | 0     | 0       | 0      | 1    | 6572                   | 6572                 | 2                  |

Interpretación: el bloque CO2 seleccionado converge satisfactoriamente. El ajuste integrado aromático actual quedó como reconciliación/validación con parámetros heredados de la selección aromática y CO2 fijo; el campo `success=False` debe leerse como evidencia de que esta pasada no fue una reoptimización global larga. Para una versión final cerrada del anexo, conviene ejecutar un multistart integrado más largo si se quiere declarar convergencia completa de todos los parámetros aromáticos bajo el CO2 seleccionado.

## 8. Resultados de ajuste

Las métricas finales se reportan como RMSE relativo, sesgo relativo y correlación por estado o pool observado.

| group     | state          | species         | n   | relative_rmse | relative_bias | corr    |
| --------- | -------------- | --------------- | --- | ------------- | ------------- | ------- |
| core      | E              | nan             | 129 | 0.1781        | 0.01296       | 0.9092  |
| core      | F              | nan             | 129 | 0.5689        | -0.2664       | 0.87    |
| core      | G              | nan             | 129 | 0.5388        | 0.03342       | 0.9049  |
| core      | Gly            | nan             | 129 | 0.1084        | -0.000422     | 0.9103  |
| core      | N              | nan             | 125 | 0.9442        | -0.2906       | 0.6739  |
| core      | X              | nan             | 122 | 0.338         | -0.09685      | 0.6451  |
| secondary | AcAld          | nan             | 105 | 0.1128        | -0.02193      | 0.6353  |
| secondary | Pyr            | nan             | 127 | 0.3309        | -0.01179      | 0.4308  |
| aroma     | nan            | ethyl_acetate   | 17  | 2.483         | 1.714         | 0.5031  |
| aroma     | nan            | ethyl_acetate   | 17  | 0.5692        | -0.4287       | 0.02023 |
| aroma     | nan            | ethyl_acetate   | 51  | 0.4552        | -0.2648       | 0.627   |
| aroma     | nan            | ethyl_octanoate | 44  | 1.047         | -0.7051       | 0.4994  |
| aroma     | nan            | ethyl_octanoate | 44  | 0.6279        | -0.01445      | 0.3338  |
| aroma     | nan            | ethyl_octanoate | 51  | 0.5605        | -0.1897       | 0.627   |
| aroma     | nan            | isoamyl_acetate | 44  | 0.9022        | -0.5369       | 0.4696  |
| aroma     | nan            | isoamyl_acetate | 44  | 0.3523        | -0.1105       | 0.6426  |
| aroma     | nan            | isoamyl_acetate | 51  | 0.39          | -0.1537       | 0.7746  |
| co2       | CO2_flow_L_min | nan             | 191 | 0.779         | -0.4388       | 0.8429  |
| co2       | CO2_flow_L_min | nan             | 110 | 0.435         | -0.02541      | 0.9312  |

Figuras principales:

- `figures/fig_04_final_fit_relative_rmse.png`
- `figures/fig_15_final_fit_representative_25170.png`
- `figures/fig_16_final_fit_representative_25171.png`

Resultados generales:

- Estados primarios como etanol y glicerol quedan bien representados.
- Glucosa/fructosa y nitrógeno presentan errores relativos mayores, esperables por heterogeneidad de mosto natural y medición puntual.
- Isoamyl acetate total y retained quedan en rango útil para modelamiento MPCC; condensado mantiene incertidumbre.
- Ethyl acetate sigue siendo el aroma más desafiante, por eso se seleccionó una estructura enriquecida con etanol, biomasa y limitación de N.
- CO2 mejora al incorporar transición O2/solubilidad, especialmente para el inicio efectivo de fermentación.

## 9. Estimabilidad y FIM

La FIM se calculó por diferencias finitas en espacio log-paramétrico:

\[
J_{ij}\approx
\frac{r_i(\theta_j e^h)-r_i(\theta_j e^{-h})}{2h}
\]

\[
F=J^TJ
\]

La covarianza aproximada:

\[
\Sigma_{\log\theta}\approx F^{-1}
\]

Clasificación con datos actuales:

| classification     | count |
| ------------------ | ----- |
| well_estimated     | 22    |
| weak_or_confounded | 9     |
| moderate           | 3     |

Después de la campaña MBDoE global:

| classification      | count |
| ------------------- | ----- |
| well_estimated      | 26    |
| weak_or_confounded  | 7     |
| weak_but_actionable | 1     |

Para el bloque CO2/O2, la campaña propuesta transforma ambos parámetros CO2 principales en bien estimados:

| parameter      | theta  | std_log_approx | approx_95_multiplier | classification |
| -------------- | ------ | -------------- | -------------------- | -------------- |
| kCO2_release_h | 0.2669 | 0.126          | 1.28                 | well_estimated |
| CO2sat_scale   | 0.4066 | 0.08031        | 1.17                 | well_estimated |

Figuras:

- `figures/fig_05_global_estimability_current.png`
- `figures/fig_06_global_estimability_after_campaign.png`
- `figures/fig_07_co2_estimability_current.png`
- `figures/fig_08_co2_estimability_after_campaign.png`
- `figures/fig_17_*` a `fig_20_*` para espectros de eigenvalores.

## 10. MBDoE

El diseño experimental usa FIM aditiva:

\[
F_{total}=F_{data}+\sum_k F_{candidate,k}
\]

El criterio usado es híbrido:

\[
\Phi_{hybrid}=
\log\det(F)
-2\left|\log\left(\frac{\lambda_{min}}{\lambda_{max}}\right)\right|
-0.05\log(\operatorname{tr}(F^{-1}))
\]

Este criterio conserva D-optimalidad como base, pero penaliza direcciones débiles.

Campaña global seleccionada:

| campaign_order | candidate                              | family                 | temperature_segments | N_pulses_kg_m3       | score |
| -------------- | -------------------------------------- | ---------------------- | -------------------- | -------------------- | ----- |
| 1              | natural_pilot_cold_to_warm_earlyN      | temperature_N          | 14, 16, 22, 20       | 30h:0.045            | 72.97 |
| 2              | natural_pilot_EA_warm_early_noN        | EA_temperature_Nstress | 24, 23, 19, 17       | nan                  | 82.94 |
| 3              | natural_pilot_EA_cold_retention_noN    | EA_retention_reference | 13, 13, 16, 18       | nan                  | 89.43 |
| 4              | natural_pilot_EA_cold_warm_Nsplit      | EA_temperature_Nsplit  | 13, 18, 23, 19       | 32h:0.025; 80h:0.035 | 94.5  |
| 5              | natural_pilot_warm_to_cool_noN         | temperature            | 22, 22, 17, 16       | nan                  | 98.38 |
| 6              | natural_pilot_low_temp_aroma_retention | aroma_retention        | 13, 15, 16, 16       | 54h:0.035            | 101.6 |

Campaña enfocada CO2/O2 seleccionada:

| campaign_order | candidate                                   | family                 | temperature_segments | N_pulses_kg_m3 | score |
| -------------- | ------------------------------------------- | ---------------------- | -------------------- | -------------- | ----- |
| 1              | natural_pilot_noN_dynamic_temperature       | temperature            | 15, 22, 18, 22       | nan            | 6.05  |
| 2              | natural_pilot_CO2_warm_start_cool_retention | CO2_strip_retention    | 23, 22, 16, 15       | nan            | 6.667 |
| 3              | natural_pilot_warm_to_cool_noN              | temperature            | 22, 22, 17, 16       | nan            | 6.961 |
| 4              | natural_pilot_EA_warm_early_noN             | EA_temperature_Nstress | 24, 23, 19, 17       | nan            | 7.168 |
| 5              | natural_pilot_CO2_Npulse_release_probe      | CO2_N_temperature      | 18, 18, 23, 20       | 50h:0.045      | 7.359 |
| 6              | natural_pilot_cold_to_warm_earlyN           | temperature_N          | 14, 16, 22, 20       | 30h:0.045      | 7.541 |

Las figuras `figures/fig_09_*` y `figures/fig_10_*` muestran los perfiles de temperatura y pulsos de N.

## 11. Conclusiones

1. Se dispone de un simulador ODE documentado y calibrado parcialmente que conecta fermentación primaria, secundarios, CO2 y aromas.
2. La síntesis aromática no queda adecuadamente representada por una cinética de fase simple para ethyl acetate; la estructura seleccionada incorpora dependencia de etanol/biomasa y limitación de N.
3. La volatilización se modela por partición gas-líquido y stripping impulsado por CO2, no como una pérdida empírica aislada.
4. El retraso de CO2 se explica mejor con transición O2 + solubilidad de CO2 que con un lag arbitrario.
5. La FIM muestra parámetros bien estimables y parámetros aún confundidos; el MBDoE propuesto apunta a mejorar direcciones débiles.
6. Para uso MPCC inmediato, conviene transferir como robustos los parámetros bien estimados y usar bounds/priors fuertes en parámetros confounded.

## 12. Pasos siguientes

- Ejecutar un multistart integrado largo para cerrar formalmente convergencia NLP de aromas bajo el CO2 seleccionado.
- Recolectar CO2 online y condensado final en nuevos ensayos diseñados.
- Mantener ethyl acetate, isoamyl acetate y ethyl octanoate como salidas MPCC, pero tratar ethyl acetate con mayor incertidumbre.
- Revisar si el MPCC debe usar todos los parámetros libres o una versión reducida con parámetros poco estimables fijados.
- Usar la campaña MBDoE como base para priorizar nuevos experimentos con mosto natural.
