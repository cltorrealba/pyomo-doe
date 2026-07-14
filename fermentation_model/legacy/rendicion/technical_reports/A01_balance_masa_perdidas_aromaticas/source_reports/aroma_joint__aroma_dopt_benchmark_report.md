# D-opt benchmark for joint fermentation-aroma DOE

This benchmark reuses the successful candidate FIMs from the joint fermentation-aroma run.
No new Ipopt/Pyomo DoE sensitivity solves are performed.

Compared designs:

- `hybrid_current`: selected with the hybrid score used in the current campaign.
- `dopt_greedy`: selected from scratch with pure D-optimality, `logdet(F_post)`.
- `dopt_exchange_from_hybrid`: one-swap local search initialized at `hybrid_current`, accepting only pure D-opt logdet gains.

## Summary

| design                    |   n_experiments |   logdet |   min_relative_eigenvalue |   condition_number |   trace_inv |   mean_var_reduction |   worst_var_reduction | worst_parameter   |   sN_var_reduction |   qN_var_reduction |   mu0_var_reduction |   k_EA_growth_var_reduction |   k_IAA_growth_var_reduction |   k_EO_growth_var_reduction |
|:--------------------------|----------------:|---------:|--------------------------:|-------------------:|------------:|---------------------:|----------------------:|:------------------|-------------------:|-------------------:|--------------------:|----------------------------:|-----------------------------:|----------------------------:|
| hybrid_current            |               9 |   155.98 |               1.29597e-05 |            77162   |    0.281606 |             0.973679 |              0.768944 | sN                |           0.768944 |           0.931409 |            0.953276 |                    0.990361 |                     0.992743 |                    0.984554 |
| dopt_greedy               |               9 |   156.48 |               1.17413e-05 |            85169.5 |    0.292098 |             0.966586 |              0.655153 | sN                |           0.655153 |           0.907198 |            0.935281 |                    0.990861 |                     0.993047 |                    0.985345 |
| dopt_exchange_from_hybrid |               9 |   156.48 |               1.17413e-05 |            85169.5 |    0.292098 |             0.966586 |              0.655153 | sN                |           0.655153 |           0.907198 |            0.935281 |                    0.990861 |                     0.993047 |                    0.985345 |

## Selected campaigns

| design                    |   campaign_order | candidate                           | family           |   horizon_h | temperature_c   |
|:--------------------------|-----------------:|:------------------------------------|:-----------------|------------:|:----------------|
| hybrid_current            |                1 | fructose_rich_iG_probe              | sugar_initial    |         192 | 18, 20, 24, 22  |
| hybrid_current            |                2 | cold_synthesis_hot_stripping        | aroma_partition  |         192 | 15, 16, 24, 25  |
| hybrid_current            |                3 | low_temp_growth_separation_plus_co2 | legacy_plus      |         168 | 15, 17, 20, 22  |
| hybrid_current            |                4 | combined_stress_long_horizon        | integrated       |         240 | 15, 22, 25, 18  |
| hybrid_current            |                5 | yan_saturation_scan                 | saturation       |         144 | 15, 18, 22, 22  |
| hybrid_current            |                6 | glucose_rich_growth_yield           | sugar_initial    |         168 | 18, 22, 24, 20  |
| hybrid_current            |                7 | high_biomass_low_N_maintenance      | maintenance      |         168 | 16, 18, 20, 18  |
| hybrid_current            |                8 | ethanol_initial_challenge           | ethanol          |         192 | 18, 22, 25, 22  |
| hybrid_current            |                9 | glucose_pulse_after_N_depletion     | pulse_separation |         192 | 17, 20, 24, 22  |
| dopt_greedy               |                1 | fructose_rich_iG_probe              | sugar_initial    |         192 | 18, 20, 24, 22  |
| dopt_greedy               |                2 | combined_stress_long_horizon        | integrated       |         240 | 15, 22, 25, 18  |
| dopt_greedy               |                3 | glucose_rich_growth_yield           | sugar_initial    |         168 | 18, 22, 24, 20  |
| dopt_greedy               |                4 | low_temp_growth_separation_plus_co2 | legacy_plus      |         168 | 15, 17, 20, 22  |
| dopt_greedy               |                5 | high_biomass_low_N_maintenance      | maintenance      |         168 | 16, 18, 20, 18  |
| dopt_greedy               |                6 | cold_synthesis_hot_stripping        | aroma_partition  |         192 | 15, 16, 24, 25  |
| dopt_greedy               |                7 | ethanol_initial_challenge           | ethanol          |         192 | 18, 22, 25, 22  |
| dopt_greedy               |                8 | glucose_pulse_after_N_depletion     | pulse_separation |         192 | 17, 20, 24, 22  |
| dopt_greedy               |                9 | viable_biomass_step                 | biomass_input    |         168 | 18, 22, 22, 18  |
| dopt_exchange_from_hybrid |                1 | fructose_rich_iG_probe              | sugar_initial    |         192 | 18, 20, 24, 22  |
| dopt_exchange_from_hybrid |                2 | cold_synthesis_hot_stripping        | aroma_partition  |         192 | 15, 16, 24, 25  |
| dopt_exchange_from_hybrid |                3 | low_temp_growth_separation_plus_co2 | legacy_plus      |         168 | 15, 17, 20, 22  |
| dopt_exchange_from_hybrid |                4 | combined_stress_long_horizon        | integrated       |         240 | 15, 22, 25, 18  |
| dopt_exchange_from_hybrid |                5 | viable_biomass_step                 | biomass_input    |         168 | 18, 22, 22, 18  |
| dopt_exchange_from_hybrid |                6 | glucose_rich_growth_yield           | sugar_initial    |         168 | 18, 22, 24, 20  |
| dopt_exchange_from_hybrid |                7 | high_biomass_low_N_maintenance      | maintenance      |         168 | 16, 18, 20, 18  |
| dopt_exchange_from_hybrid |                8 | ethanol_initial_challenge           | ethanol          |         192 | 18, 22, 25, 22  |
| dopt_exchange_from_hybrid |                9 | glucose_pulse_after_N_depletion     | pulse_separation |         192 | 17, 20, 24, 22  |

## Exchange log

|   iteration | swap_out            | swap_in             |   previous_logdet |   new_logdet |   gain_logdet |   out_index |
|------------:|:--------------------|:--------------------|------------------:|-------------:|--------------:|------------:|
|           1 | yan_saturation_scan | viable_biomass_step |            155.98 |       156.48 |      0.499085 |           4 |

## Weakest reductions by design

| design                    | parameter   | group        |   var_ratio |   var_reduction |
|:--------------------------|:------------|:-------------|------------:|----------------:|
| dopt_exchange_from_hybrid | sN          | fermentation |   0.344847  |        0.655153 |
| dopt_exchange_from_hybrid | qN          | fermentation |   0.0928016 |        0.907198 |
| dopt_exchange_from_hybrid | mu0         | fermentation |   0.0647193 |        0.935281 |
| dopt_exchange_from_hybrid | sF          | fermentation |   0.0408536 |        0.959146 |
| dopt_exchange_from_hybrid | sG          | fermentation |   0.020391  |        0.979609 |
| dopt_exchange_from_hybrid | qEG         | fermentation |   0.0178894 |        0.982111 |
| dopt_greedy               | sN          | fermentation |   0.344847  |        0.655153 |
| dopt_greedy               | qN          | fermentation |   0.0928016 |        0.907198 |
| dopt_greedy               | mu0         | fermentation |   0.0647193 |        0.935281 |
| dopt_greedy               | sF          | fermentation |   0.0408536 |        0.959146 |
| dopt_greedy               | sG          | fermentation |   0.020391  |        0.979609 |
| dopt_greedy               | qEG         | fermentation |   0.0178894 |        0.982111 |
| hybrid_current            | sN          | fermentation |   0.231056  |        0.768944 |
| hybrid_current            | qN          | fermentation |   0.0685908 |        0.931409 |
| hybrid_current            | mu0         | fermentation |   0.0467238 |        0.953276 |
| hybrid_current            | sF          | fermentation |   0.0442392 |        0.955761 |
| hybrid_current            | sG          | fermentation |   0.0227014 |        0.977299 |
| hybrid_current            | qXG         | fermentation |   0.0188627 |        0.981137 |
