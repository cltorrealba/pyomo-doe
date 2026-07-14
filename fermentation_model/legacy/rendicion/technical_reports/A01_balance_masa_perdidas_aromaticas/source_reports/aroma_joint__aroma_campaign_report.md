# Aroma campaign DOE

This report evaluates model-based candidate experiments for aroma and fermentation parameters.

The gas-liquid loss model is fixed by UNIFAC-derived equilibrium coefficients. The target vector is:

`mu0`, `sN`, `qN`, `qXG`, `qXF`, `betaG0`, `sG`, `betaF0`, `sF`, `qEG`, `qEF`, `iG`, `iE`, `Kd0`, `m0`, `k_EA_growth`, `k_EA_stationary`, `k_IAA_growth`, `k_IAA_stationary`, `k_EO_growth`, `k_EO_stationary`

Liquid aroma observations are included at the liquid sampling interval. A single terminal condensate observation per aroma is included as a final accumulated gas-loss constraint.

## Selected Campaign

|   campaign_order | candidate                           | family           |   horizon_h | temperature_c   |   campaign_logdet |   campaign_min_relative_eigenvalue |   campaign_condition_number |   mean_var_reduction |   worst_var_reduction | rationale                                                                                                       |
|-----------------:|:------------------------------------|:-----------------|------------:|:----------------|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:----------------------------------------------------------------------------------------------------------------|
|                1 | fructose_rich_iG_probe              | sugar_initial    |         192 | 18, 20, 24, 22  |           115.984 |                        3.97464e-06 |                    251595   |             0.841626 |              0.30235  | Fructose-rich must plus glucose pulse separates fructose capacity from glucose inhibition iG.                   |
|                2 | cold_synthesis_hot_stripping        | aroma_partition  |         192 | 15, 16, 24, 25  |           129.769 |                        3.7578e-06  |                    266113   |             0.905601 |              0.479869 | Low-temperature synthesis window followed by high-temperature stripping challenge.                              |
|                3 | low_temp_growth_separation_plus_co2 | legacy_plus      |         168 | 15, 17, 20, 22  |           135.931 |                        8.83091e-06 |                    113239   |             0.931126 |              0.541207 | Previous robust N/T design, now measured with Xd and online CO2.                                                |
|                4 | combined_stress_long_horizon        | integrated       |         240 | 15, 22, 25, 18  |           143.208 |                        7.53078e-06 |                    132788   |             0.942065 |              0.566126 | Broad stress/input design for integrated all-parameter information.                                             |
|                5 | yan_saturation_scan                 | saturation       |         144 | 15, 18, 22, 22  |           145.211 |                        8.59293e-06 |                    116375   |             0.95764  |              0.712291 | Designed N ladder through low/intermediate YAN targets sN/qN separation.                                        |
|                6 | glucose_rich_growth_yield           | sugar_initial    |         168 | 18, 22, 24, 20  |           149.394 |                        9.92551e-06 |                    100750   |             0.965031 |              0.730844 | High glucose must isolates glucose growth/fermentation terms and tests fructose response after glucose history. |
|                7 | high_biomass_low_N_maintenance      | maintenance      |         168 | 16, 18, 20, 18  |           152.425 |                        1.17842e-05 |                     84859.7 |             0.968769 |              0.734643 | High biomass and low N create low-growth sugar consumption windows for m0.                                      |
|                8 | ethanol_initial_challenge           | ethanol          |         192 | 18, 22, 25, 22  |           154.452 |                        1.19934e-05 |                     83379.2 |             0.971282 |              0.754084 | Initial ethanol decouples ethanol inhibition iE from ethanol produced by fermentation.                          |
|                9 | glucose_pulse_after_N_depletion     | pulse_separation |         192 | 17, 20, 24, 22  |           155.98  |                        1.29597e-05 |                     77162   |             0.973679 |              0.768944 | Glucose pulse after likely N limitation separates fermentation/maintenance from growth uptake.                  |

## Candidate Ranking

| candidate                           | family           |   horizon_h |   combined_logdet |   combined_min_relative_eigenvalue |   combined_condition_number |   mean_var_reduction |   worst_var_reduction | status   | error                                                                                               |
|:------------------------------------|:-----------------|------------:|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:---------|:----------------------------------------------------------------------------------------------------|
| fructose_rich_iG_probe              | sugar_initial    |         192 |          115.984  |                        3.97464e-06 |                      251595 |             0.841626 |            0.30235    | ok       |                                                                                                     |
| cold_synthesis_hot_stripping        | aroma_partition  |         192 |          113.268  |                        2.65968e-06 |                      375985 |             0.819833 |            0.397601   | ok       |                                                                                                     |
| glucose_rich_growth_yield           | sugar_initial    |         168 |          113.709  |                        3.00202e-06 |                      333109 |             0.837542 |            0.339418   | ok       |                                                                                                     |
| combined_stress_long_horizon        | integrated       |         240 |          114.294  |                        3.61613e-06 |                      276538 |             0.789556 |            0.162312   | ok       |                                                                                                     |
| ethanol_partition_late_pulse        | aroma_partition  |         216 |          109.435  |                        2.92504e-06 |                      341875 |             0.756664 |            0.15218    | ok       |                                                                                                     |
| low_temp_growth_separation_plus_co2 | legacy_plus      |         168 |          103.081  |                        3.15671e-06 |                      316786 |             0.690489 |            0.406587   | ok       |                                                                                                     |
| ethanol_initial_challenge           | ethanol          |         192 |          107.064  |                        3.32112e-06 |                      301103 |             0.73816  |            0.126686   | ok       |                                                                                                     |
| viable_biomass_step                 | biomass_input    |         168 |          101.872  |                        1.6429e-06  |                      608680 |             0.647366 |            0.295954   | ok       |                                                                                                     |
| glucose_pulse_after_N_depletion     | pulse_separation |         192 |           92.2315 |                        1.63573e-06 |                      611349 |             0.549856 |            0.0129036  | ok       |                                                                                                     |
| high_biomass_low_N_maintenance      | maintenance      |         168 |           92.9609 |                        1.31022e-06 |                      763233 |             0.483228 |            0.00401706 | ok       |                                                                                                     |
| yan_saturation_scan                 | saturation       |         144 |           86.2307 |                        1.43768e-06 |                      695566 |             0.458202 |            0.0323474  | ok       |                                                                                                     |
| hot_synthesis_cold_retention        | aroma_partition  |         192 |          nan      |                      nan           |                         nan |           nan        |          nan          | failed   | RuntimeError: Model from experiment did not solve appropriately. Make sure the model is well-posed. |
| co2_stripping_without_growth        | aroma_partition  |         192 |          nan      |                      nan           |                         nan |           nan        |          nan          | failed   | RuntimeError: Model from experiment did not solve appropriately. Make sure the model is well-posed. |

## Parameter Variance Reductions

| parameter        | group        |   campaign_var_ratio |   campaign_var_reduction |
|:-----------------|:-------------|---------------------:|-------------------------:|
| mu0              | fermentation |          0.0467238   |                 0.953276 |
| sN               | fermentation |          0.231056    |                 0.768944 |
| qN               | fermentation |          0.0685908   |                 0.931409 |
| qXG              | fermentation |          0.0188627   |                 0.981137 |
| qXF              | fermentation |          0.0177801   |                 0.98222  |
| betaG0           | fermentation |          0.0117314   |                 0.988269 |
| sG               | fermentation |          0.0227014   |                 0.977299 |
| betaF0           | fermentation |          0.00316916  |                 0.996831 |
| sF               | fermentation |          0.0442392   |                 0.955761 |
| qEG              | fermentation |          0.0173885   |                 0.982611 |
| qEF              | fermentation |          0.00293252  |                 0.997067 |
| iG               | fermentation |          0.00222419  |                 0.997776 |
| iE               | fermentation |          0.0160583   |                 0.983942 |
| Kd0              | fermentation |          0.000959127 |                 0.999041 |
| m0               | fermentation |          0.0158841   |                 0.984116 |
| k_EA_growth      | synthesis    |          0.00963876  |                 0.990361 |
| k_EA_stationary  | synthesis    |          4.32661e-05 |                 0.999957 |
| k_IAA_growth     | synthesis    |          0.00725713  |                 0.992743 |
| k_IAA_stationary | synthesis    |          1.7575e-05  |                 0.999982 |
| k_EO_growth      | synthesis    |          0.0154458   |                 0.984554 |
| k_EO_stationary  | synthesis    |          3.8195e-05  |                 0.999962 |

## Interpretation

- The terminal condensate improves the mass closure but does not provide time-resolved gas-loss dynamics.
- Partition parameters are fixed by default to avoid confounding gas-liquid equilibrium with trap efficiency.
- Temperature switches are intentionally included to separate synthesis windows from stripping windows.
- Nitrogen ladders and sugar pulses are retained because the synthesis model is phase-dependent through nitrogen limitation and sugar uptake.