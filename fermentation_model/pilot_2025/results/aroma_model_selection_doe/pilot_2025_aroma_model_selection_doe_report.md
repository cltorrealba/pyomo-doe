# Pilot 2025 Aroma Model Selection and MBDoE

## Selected structure

Selected model: `ea_ethanol_nlimited`.

Tests whether ethyl acetate increases when ethanol is available and nitrogen is depleted.

$$r_{EA}=r_{EA,phase}+k_{EA,XE}X\frac{E}{K_E+E}+k_{EA,XE,Nlim}X\frac{E}{K_E+E}\frac{K_N}{K_N+N}$$

## Gas-liquid equilibrium and CO2 stripping

Aroma loss is computed as `r_loss,i = alpha_i K_i(T,E) q_CO2 C_L,i`. `K_i` is evaluated directly from Antoine vapor pressure and original UNIFAC activity coefficients for a dilute aroma in a water-ethanol liquid mixture. The CO2 stripping proportionality follows the Mouret-style mass-balance assumption that volatile loss scales with fermentation CO2 flow.

Condensate entries equal to exactly zero are treated as missing/below-reporting-limit values, not as exact dynamic measurements of zero accumulated condensate. Positive condensate values are used to reconstruct the retained liquid pool and to constrain the accumulated loss state.

## Model selection table

| model                    |   n_parameters |   data_wsse |   wsse_per_data_residual |    aicc |     bic |   selection_score |   active_bound_count |   ea_retained_rmse |   ea_retained_relative_rmse |   ea_retained_relative_bias |   ea_total_rmse |   ea_total_relative_rmse |   ea_total_relative_bias |
|:-------------------------|---------------:|------------:|-------------------------:|--------:|--------:|------------------:|---------------------:|-------------------:|----------------------------:|----------------------------:|----------------:|-------------------------:|-------------------------:|
| ea_ethanol_nlimited      |             11 |     7369.84 |                  28.2369 | 7392.9  | 7431.05 |           7478.04 |                    2 |            4.77706 |                    0.575972 |                   -0.437476 |         3.61551 |                 0.476356 |                -0.285364 |
| ea_redox_acetaldehyde    |             11 |     7493.34 |                  28.7101 | 7516.4  | 7554.55 |           7603.94 |                    2 |            5.05853 |                    0.609909 |                   -0.46728  |         3.69386 |                 0.486679 |                -0.292988 |
| ea_combined_parsimonious |             12 |     7490.27 |                  28.6984 | 7515.53 | 7557.05 |           7626.45 |                    3 |            5.06076 |                    0.610178 |                   -0.467516 |         3.69565 |                 0.486915 |                -0.293037 |
| ea_ethanol_biomass       |             10 |     7696.92 |                  29.4901 | 7717.8  | 7752.56 |           7798.51 |                    2 |            4.61921 |                    0.55694  |                   -0.424286 |         3.49929 |                 0.461044 |                -0.279789 |
| ea_ethanol_temperature   |             11 |     7696.93 |                  29.4902 | 7719.99 | 7758.14 |           7824.09 |                    3 |            4.6192  |                    0.556939 |                   -0.424285 |         3.49928 |                 0.461043 |                -0.279788 |
| ea_biomass_background    |             10 |     8249.51 |                  31.6073 | 8270.39 | 8305.15 |           8349.98 |                    2 |            4.45073 |                    0.536626 |                   -0.410357 |         3.35281 |                 0.441744 |                -0.273907 |
| baseline_phase           |              9 |     9448.98 |                  36.203  | 9467.7  | 9499.06 |           9541.15 |                    2 |            4.25988 |                    0.513615 |                   -0.376019 |         3.24906 |                 0.428075 |                -0.28147  |

## Selected model estimability

| parameter        |       theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:-----------------|------------:|-----------------:|-----------------------:|:---------------|:-------------------|
| mu0              | 0.0689711   |        0.0116908 |            1.02318     | False          | well_estimated     |
| qN               | 0.0157715   |        0.0150312 |            1.0299      | False          | well_estimated     |
| betaG0           | 1.18634     |        0.042631  |            1.08715     | False          | well_estimated     |
| betaF0           | 0.31451     |        0.211288  |            1.51305     | False          | well_estimated     |
| qEG              | 1.11251     |        0.0399744 |            1.0815      | False          | well_estimated     |
| qEF              | 1.27801     |        0.0994617 |            1.21524     | False          | well_estimated     |
| iG               | 0.0139708   |        0.204639  |            1.49345     | False          | well_estimated     |
| iE               | 0.0156601   |        0.0866282 |            1.18506     | False          | well_estimated     |
| Kd0              | 0.000656095 |        0.0829246 |            1.17649     | False          | well_estimated     |
| gammaG0          | 0.0703099   |        0.0937404 |            1.20169     | False          | well_estimated     |
| gammaF0          | 0.00339444  |        2.00669   |           51.0659      | False          | weak_or_confounded |
| kPyrS_N          | 1.18435     |        0.0709523 |            1.1492      | False          | well_estimated     |
| kPyrO2           | 0.152368    |        0.292818  |            1.77522     | False          | well_estimated     |
| kPyrDrain        | 0.00222636  |        0.151242  |            1.34505     | False          | well_estimated     |
| kAldS_N          | 3.92542     |        0.0353663 |            1.07178     | False          | well_estimated     |
| kAldRed          | 1.00049e-05 |       15.3857    |            1.24895e+13 | True           | weak_or_confounded |
| kAcAld           | 0.00267744  |        0.0425138 |            1.0869      | False          | well_estimated     |
| kAcStress        | 0.00261587  |      363.611     |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_growth      | 0.0001      |       90.2487    |            1.05765e+17 | True           | weak_or_confounded |
| k_EA_stationary  | 0.00594688  |        0.455476  |            2.44179     | False          | moderate           |
| k_IAA_growth     | 0.053515    |        0.138806  |            1.31267     | False          | well_estimated     |
| k_IAA_stationary | 0.0400843   |        0.0592003 |            1.12303     | False          | well_estimated     |
| k_EO_growth      | 0.000898158 |        0.11657   |            1.25669     | False          | well_estimated     |
| k_EO_stationary  | 0.000149817 |        0.180195  |            1.42359     | False          | well_estimated     |
| alpha_EA_loss    | 0.2         |        0.0483356 |            1.09937     | True           | weak_or_confounded |
| alpha_IAA_loss   | 0.660648    |        0.0450883 |            1.0924      | False          | well_estimated     |
| alpha_EO_loss    | 1.47819     |        0.0506068 |            1.10428     | False          | well_estimated     |
| k_EA_XE          | 1.05635e-05 |      362.609     |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_XE_Nlim     | 0.0417864   |        0.106084  |            1.23112     | False          | well_estimated     |

## Weak FIM directions

| analysis         |   weak_direction |   eigenvalue | dominant_parameters                                                                             | dominant_abs_loadings                                  |
|:-----------------|-----------------:|-------------:|:------------------------------------------------------------------------------------------------|:-------------------------------------------------------|
| selected_current |                1 | -6.37405e-14 | kAcStress, k_EA_XE, k_EA_growth, k_EA_XE_Nlim, k_EA_stationary, kAldRed, gammaF0, alpha_EA_loss | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| selected_current |                2 |  1.88364e-08 | k_EA_XE, k_EA_growth, k_EA_XE_Nlim, kAldRed, k_EA_stationary, alpha_EA_loss, gammaF0, betaG0    | 0.998, 0.057, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| selected_current |                3 |  0.000121548 | k_EA_growth, k_EA_XE, k_EA_stationary, kAldRed, gammaF0, alpha_EA_loss, k_EA_XE_Nlim, betaF0    | 0.998, 0.057, 0.003, 0.001, 0.000, 0.000, 0.000, 0.000 |
| selected_current |                4 |  0.0042169   | kAldRed, kAcAld, gammaF0, k_EA_growth, iG, betaF0, kAldS_N, qEF                                 | 1.000, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001, 0.000 |
| selected_current |                5 |  0.247756    | gammaF0, gammaG0, betaF0, iG, qEF, iE, qEG, k_IAA_growth                                        | 0.999, 0.041, 0.021, 0.015, 0.009, 0.008, 0.004, 0.002 |
| selected_current |                6 |  6.8383      | k_EA_stationary, k_EA_XE_Nlim, alpha_EA_loss, Kd0, kPyrO2, kPyrDrain, iG, betaF0                | 0.989, 0.122, 0.076, 0.020, 0.019, 0.010, 0.005, 0.004 |
| selected_current |                7 |  9.1124      | kPyrO2, kPyrDrain, kPyrS_N, iG, betaF0, k_EA_stationary, qEF, iE                                | 0.882, 0.443, 0.153, 0.033, 0.027, 0.022, 0.016, 0.012 |
| selected_current |                8 | 12.0174      | betaF0, iG, qEF, iE, gammaG0, qEG, kPyrO2, gammaF0                                              | 0.650, 0.627, 0.327, 0.230, 0.108, 0.092, 0.042, 0.024 |

## Candidate ranking

| candidate                             | family                 |   combined_logdet |   combined_min_relative_eigenvalue |   aroma_mean_var_reduction |   aroma_worst_var_reduction | N_pulses_kg_m3       | rationale                                                                                                            |
|:--------------------------------------|:-----------------------|------------------:|-----------------------------------:|---------------------------:|----------------------------:|:---------------------|:---------------------------------------------------------------------------------------------------------------------|
| natural_pilot_cold_to_warm_earlyN     | temperature_N          |           139.125 |                        1.59216e-12 |                   0.240712 |                   0.0985989 | 30h:0.045            | Cold start followed by warm transition and early nitrogen pulse.                                                     |
| natural_pilot_midN_temperature_step   | temperature_N          |           138.996 |                        1.56743e-12 |                   0.252513 |                   0.0890028 | 54h:0.04             | Temperature step with mid-growth nitrogen perturbation.                                                              |
| natural_pilot_EA_cold_retention_noN   | EA_retention_reference |           138.743 |                        1.65226e-12 |                   0.258715 |                   0.111257  |                      | Low-temperature no-pulse comparator to decouple synthesis from CO2 stripping.                                        |
| natural_pilot_noN_dynamic_temperature | temperature            |           138.587 |                        1.70244e-12 |                   0.264407 |                   0.107629  |                      | Temperature-only perturbation for settings where nutrient action is constrained.                                     |
| natural_pilot_two_step_N_ladder       | N_timing               |           138.441 |                        1.50894e-12 |                   0.247371 |                   0.0882186 | 30h:0.03; 78h:0.035  | Two smaller nitrogen pulses to separate early growth and later metabolic response.                                   |
| natural_pilot_EA_warm_early_noN       | EA_temperature_Nstress |           138.683 |                        3.25551e-13 |                   0.269924 |                   0.0697817 |                      | Warm early natural fermentation without nutrient pulse to excite ethanol-biomass and nitrogen-limited EA directions. |
| natural_pilot_high_rate_strip         | co2_aroma              |           138.855 |                        2.73916e-13 |                   0.259047 |                   0.0953    | 30h:0.03             | High-rate natural fermentation to excite CO2 stripping and aroma loss directions.                                    |
| natural_pilot_EA_warm_early_lateN     | EA_lateN               |           138.532 |                        3.13525e-13 |                   0.257254 |                   0.0732566 | 98h:0.045            | Warm early fermentation followed by late nitrogen to separate stationary and N-limited aroma terms.                  |
| natural_pilot_warm_to_cool_noN        | temperature            |           137.995 |                        3.28179e-13 |                   0.267494 |                   0.0787818 |                      | Warm early phase to excite growth and CO2, then cool aroma-retention phase.                                          |
| natural_pilot_reference_20C           | reference              |           138.059 |                        3.04159e-13 |                   0.268911 |                   0.112586  |                      | Natural must warmer reference to increase rate and CO2 information.                                                  |
| natural_pilot_EA_cold_warm_Nsplit     | EA_temperature_Nsplit  |           138.945 |                        1.87771e-13 |                   0.244655 |                   0.108549  | 32h:0.025; 80h:0.035 | Cold start then warm acceleration with split nitrogen to test temperature and N response.                            |
| natural_pilot_reference_18C           | reference              |           137.637 |                        2.8285e-13  |                   0.261412 |                   0.112162  |                      | Natural must reference at moderate temperature.                                                                      |

## Selected campaign

|   campaign_order | candidate                             | family                 |   campaign_logdet |   campaign_min_relative_eigenvalue |   aroma_mean_var_reduction |   aroma_worst_var_reduction | N_pulses_kg_m3      | rationale                                                                                                            |
|-----------------:|:--------------------------------------|:-----------------------|------------------:|-----------------------------------:|---------------------------:|----------------------------:|:--------------------|:---------------------------------------------------------------------------------------------------------------------|
|                1 | natural_pilot_cold_to_warm_earlyN     | temperature_N          |           139.125 |                        1.59216e-12 |                   0.240712 |                   0.0985989 | 30h:0.045           | Cold start followed by warm transition and early nitrogen pulse.                                                     |
|                2 | natural_pilot_EA_cold_retention_noN   | EA_retention_reference |           145.786 |                        2.73032e-12 |                   0.367345 |                   0.216952  |                     | Low-temperature no-pulse comparator to decouple synthesis from CO2 stripping.                                        |
|                3 | natural_pilot_midN_temperature_step   | temperature_N          |           150.579 |                        3.58948e-12 |                   0.445493 |                   0.273756  | 54h:0.04            | Temperature step with mid-growth nitrogen perturbation.                                                              |
|                4 | natural_pilot_EA_warm_early_noN       | EA_temperature_Nstress |           154.73  |                        3.66231e-12 |                   0.517194 |                   0.325571  |                     | Warm early natural fermentation without nutrient pulse to excite ethanol-biomass and nitrogen-limited EA directions. |
|                5 | natural_pilot_noN_dynamic_temperature | temperature            |           157.762 |                        4.4772e-12  |                   0.556375 |                   0.35611   |                     | Temperature-only perturbation for settings where nutrient action is constrained.                                     |
|                6 | natural_pilot_two_step_N_ladder       | N_timing               |           160.368 |                        5.1105e-12  |                   0.586926 |                   0.388646  | 30h:0.03; 78h:0.035 | Two smaller nitrogen pulses to separate early growth and later metabolic response.                                   |