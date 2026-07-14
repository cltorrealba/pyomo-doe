# Aroma campaign DOE

This report evaluates model-based candidate experiments for aroma synthesis parameters.

The gas-liquid loss model is fixed by literature-scale equilibrium coefficients. The default target vector estimates synthesis only:

`k_EA_growth`, `k_EA_stationary`, `k_IAA_growth`, `k_IAA_stationary`, `k_EO_growth`, `k_EO_stationary`

Liquid aroma observations are included at the liquid sampling interval. A single terminal condensate observation per aroma is included as a final accumulated gas-loss constraint.

## Selected Campaign

|   campaign_order | candidate                           | family          |   horizon_h | temperature_c   |   campaign_logdet |   campaign_min_relative_eigenvalue |   campaign_condition_number |   mean_var_reduction |   worst_var_reduction | rationale                                                                                                       |
|-----------------:|:------------------------------------|:----------------|------------:|:----------------|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:----------------------------------------------------------------------------------------------------------------|
|                1 | cold_synthesis_hot_stripping        | aroma_partition |         192 | 15, 16, 24, 25  |           29.3691 |                        0.000372638 |                     2683.57 |             0.970138 |              0.915764 | Low-temperature synthesis window followed by high-temperature stripping challenge.                              |
|                2 | low_temp_growth_separation_plus_co2 | legacy_plus     |         168 | 15, 17, 20, 22  |           32.8685 |                        0.000461995 |                     2164.53 |             0.98511  |              0.957398 | Previous robust N/T design, now measured with Xd and online CO2.                                                |
|                3 | hot_synthesis_cold_retention        | aroma_partition |         192 | 24, 25, 17, 15  |           35.9532 |                        0.000357839 |                     2794.55 |             0.9899   |              0.970962 | High-temperature production window followed by cold retention phase.                                            |
|                4 | fructose_rich_iG_probe              | sugar_initial   |         192 | 18, 20, 24, 22  |           37.7139 |                        0.00033004  |                     3029.94 |             0.992162 |              0.977432 | Fructose-rich must plus glucose pulse separates fructose capacity from glucose inhibition iG.                   |
|                5 | high_biomass_low_N_maintenance      | maintenance     |         168 | 16, 18, 20, 18  |           38.5558 |                        0.000357961 |                     2793.6  |             0.993465 |              0.981128 | High biomass and low N create low-growth sugar consumption windows for m0.                                      |
|                6 | glucose_rich_growth_yield           | sugar_initial   |         168 | 18, 22, 24, 20  |           39.459  |                        0.000343582 |                     2910.52 |             0.994263 |              0.983426 | High glucose must isolates glucose growth/fermentation terms and tests fructose response after glucose history. |
|                7 | combined_stress_long_horizon        | integrated      |         240 | 15, 22, 25, 18  |           40.3503 |                        0.00030375  |                     3292.19 |             0.994743 |              0.984796 | Broad stress/input design for integrated all-parameter information.                                             |
|                8 | viable_biomass_step                 | biomass_input   |         168 | 18, 22, 22, 18  |           40.9335 |                        0.000290554 |                     3441.7  |             0.995123 |              0.985895 | Known viable biomass addition tests whether rates scale with X and improves q/yield separation.                 |
|                9 | ethanol_partition_late_pulse        | aroma_partition |         216 | 18, 20, 24, 22  |           41.464  |                        0.000271545 |                     3682.62 |             0.995382 |              0.986643 | Late ethanol pulse perturbs partition without changing aroma synthesis parameters directly.                     |

## Candidate Ranking

| candidate                           | family           |   horizon_h |   combined_logdet |   combined_min_relative_eigenvalue |   combined_condition_number |   mean_var_reduction |   worst_var_reduction | status   | error   |
|:------------------------------------|:-----------------|------------:|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:---------|:--------|
| cold_synthesis_hot_stripping        | aroma_partition  |         192 |           29.3691 |                        0.000372638 |                     2683.57 |             0.970138 |             0.915764  | ok       |         |
| fructose_rich_iG_probe              | sugar_initial    |         192 |           29.6405 |                        0.000285462 |                     3503.1  |             0.967447 |             0.907862  | ok       |         |
| low_temp_growth_separation_plus_co2 | legacy_plus      |         168 |           27.4695 |                        0.000560555 |                     1783.95 |             0.966585 |             0.905922  | ok       |         |
| glucose_rich_growth_yield           | sugar_initial    |         168 |           27.8158 |                        0.000298088 |                     3354.72 |             0.956709 |             0.878368  | ok       |         |
| viable_biomass_step                 | biomass_input    |         168 |           26.8728 |                        0.000222815 |                     4488.04 |             0.941186 |             0.83658   | ok       |         |
| combined_stress_long_horizon        | integrated       |         240 |           27.6137 |                        0.000130417 |                     7667.68 |             0.931703 |             0.812708  | ok       |         |
| hot_synthesis_cold_retention        | aroma_partition  |         192 |           27.6189 |                        0.000105138 |                     9511.31 |             0.924268 |             0.790796  | ok       |         |
| ethanol_partition_late_pulse        | aroma_partition  |         216 |           26.439  |                        0.000147378 |                     6785.27 |             0.92182  |             0.785978  | ok       |         |
| ethanol_initial_challenge           | ethanol          |         192 |           25.2322 |                        0.000155824 |                     6417.5  |             0.906447 |             0.747864  | ok       |         |
| glucose_pulse_after_N_depletion     | pulse_separation |         192 |           19.9142 |                        0.000118423 |                     8444.33 |             0.726499 |             0.34918   | ok       |         |
| yan_saturation_scan                 | saturation       |         144 |           16.445  |                        0.000186083 |                     5373.95 |             0.591875 |             0.125414  | ok       |         |
| high_biomass_low_N_maintenance      | maintenance      |         168 |           18.281  |                        7.26822e-05 |                    13758.5  |             0.50316  |             0.0040329 | ok       |         |
| co2_stripping_without_growth        | aroma_partition  |         192 |           17.3932 |                        9.94353e-05 |                    10056.8  |             0.508431 |             0.0106752 | ok       |         |

## Parameter Variance Reductions

| parameter        | group     |   campaign_var_ratio |   campaign_var_reduction |
|:-----------------|:----------|---------------------:|-------------------------:|
| k_EA_growth      | synthesis |          0.00838407  |                 0.991616 |
| k_EA_stationary  | synthesis |          3.38322e-05 |                 0.999966 |
| k_IAA_growth     | synthesis |          0.00589103  |                 0.994109 |
| k_IAA_stationary | synthesis |          1.14908e-05 |                 0.999989 |
| k_EO_growth      | synthesis |          0.0133573   |                 0.986643 |
| k_EO_stationary  | synthesis |          2.91807e-05 |                 0.999971 |

## Interpretation

- The terminal condensate improves the mass closure but does not provide time-resolved gas-loss dynamics.
- Partition parameters are fixed by default to avoid confounding gas-liquid equilibrium with trap efficiency.
- Temperature switches are intentionally included to separate synthesis windows from stripping windows.
- Nitrogen ladders and sugar pulses are retained because the synthesis model is phase-dependent through nitrogen limitation and sugar uptake.