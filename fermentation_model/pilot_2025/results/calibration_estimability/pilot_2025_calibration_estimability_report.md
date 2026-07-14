# Pilot 2025 calibration and estimability

## Scope

This run reuses the reduced extended fermentation model for pilot-scale natural-must data: primary fermentation, glycerol, secondary v2 reduced chemistry, liquid aroma synthesis, and aroma volatilization to condenser. The future design space is restricted to natural must, temperature setpoints, and nutrient additions.

## Aroma mass-balance interpretation

`*_total` is interpreted as retained-in-wine plus accumulated condenser equivalent concentration. `*_condensate` is interpreted as the accumulated condenser equivalent concentration. When both are observed, the retained wine observation is reconstructed as `total - condensate` and fitted against the liquid aroma state. The condenser observation is fitted against the cumulative volatilized/captured state, making `alpha_*_loss` estimable as effective volatilization/capture coefficients.

## CO2 curation

|   batch |   raw_rows | used_for_calibration   | reason                                                 |   activation_time_h |   curated_rows |
|--------:|-----------:|:-----------------------|:-------------------------------------------------------|--------------------:|---------------:|
|   25150 |         56 | False                  | excluded: CO2 file not usable for this batch           |                 nan |              0 |
|   25151 |         56 | False                  | excluded: CO2 file not usable for this batch           |                 nan |              0 |
|   25170 |      22877 | True                   | kept                                                   |                   0 |          22837 |
|   25171 |      20000 | True                   | trimmed to operational CO2/process activation override |                 100 |          13963 |

## Calibration stages

| fit                                |   n_batches | mediums   | parameters                                             |   n_parameters | success   |   status | message                                    |   nfev |     initial_wsse |   final_wsse |   n_residuals |   dof |   wsse_per_residual |   wsse_per_dof |   l2_lambda | l2_parameters                                          |
|:-----------------------------------|------------:|:----------|:-------------------------------------------------------|---------------:|:----------|---------:|:-------------------------------------------|-------:|-----------------:|-------------:|--------------:|------:|--------------------:|---------------:|------------:|:-------------------------------------------------------|
| pilot_core_l2_multistart_00        |           8 | natural   | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |             11 | True      |        3 | `xtol` termination condition is satisfied. |      7 |  22183.6         |     22046.3  |           786 |   775 |             28.0487 |        28.4468 |        0.35 | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |
| pilot_core_l2_multistart_01        |           8 | natural   | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |             11 | True      |        3 | `xtol` termination condition is satisfied. |     11 |  78166.1         |     15675.6  |           786 |   775 |             19.9435 |        20.2265 |        0.35 | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |
| pilot_secondary_v2_reduced_o2fixed |         nan | nan       | nan                                                    |            nan | True      |        2 | `ftol` termination condition is satisfied. |     12 |   2963.01        |      2799.05 |           237 |   nan |             11.8103 |       nan      |      nan    | nan                                                    |
| pilot_aroma_multistart_00          |         nan | nan       | nan                                                    |            nan | True      |        2 | `ftol` termination condition is satisfied. |     10 |      1.37719e+06 |      4076.15 |           297 |   nan |             13.7244 |       nan      |      nan    | nan                                                    |
| pilot_aroma_multistart_01          |         nan | nan       | nan                                                    |            nan | True      |        2 | `ftol` termination condition is satisfied. |     15 |      2.44064e+06 |      4075.73 |           297 |   nan |             13.723  |       nan      |      nan    | nan                                                    |
| pilot_aroma_multistart_02          |         nan | nan       | nan                                                    |            nan | True      |        2 | `ftol` termination condition is satisfied. |     19 | 850682           |      4074.88 |           297 |   nan |             13.7201 |       nan      |      nan    | nan                                                    |

## Practical estimability from current pilot data

| parameter        |       theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:-----------------|------------:|-----------------:|-----------------------:|:---------------|:-------------------|
| mu0              | 0.0689711   |        0.0117835 |            1.02336     | False          | well_estimated     |
| qN               | 0.0157715   |        0.0150878 |            1.03001     | False          | well_estimated     |
| betaG0           | 1.18634     |        0.0431613 |            1.08828     | False          | well_estimated     |
| betaF0           | 0.31451     |        0.211931  |            1.51495     | False          | well_estimated     |
| qEG              | 1.11251     |        0.0403553 |            1.08231     | False          | well_estimated     |
| qEF              | 1.27801     |        0.0995529 |            1.21546     | False          | well_estimated     |
| iG               | 0.0139708   |        0.20465   |            1.49349     | False          | well_estimated     |
| iE               | 0.0156601   |        0.086965  |            1.18584     | False          | well_estimated     |
| Kd0              | 0.000656095 |        0.0835033 |            1.17782     | False          | well_estimated     |
| gammaG0          | 0.0703099   |        0.0938226 |            1.20189     | False          | well_estimated     |
| gammaF0          | 0.00339444  |        2.00828   |           51.2249      | False          | weak_or_confounded |
| kPyrS_N          | 1.18435     |        0.070993  |            1.14929     | False          | well_estimated     |
| kPyrO2           | 0.152368    |        0.292822  |            1.77523     | False          | well_estimated     |
| kPyrDrain        | 0.00222636  |        0.151243  |            1.34506     | False          | well_estimated     |
| kAldS_N          | 3.92542     |        0.0354686 |            1.07199     | False          | well_estimated     |
| kAldRed          | 1.00049e-05 |       15.3865    |            1.25092e+13 | True           | weak_or_confounded |
| kAcAld           | 0.00267744  |        0.0425154 |            1.0869      | False          | well_estimated     |
| kAcStress        | 0.00261587  |      360.858     |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_growth      | 0.000776841 |       11.3458    |            4.54709e+09 | False          | weak_or_confounded |
| k_EA_stationary  | 0.00461963  |        0.60064   |            3.24545     | False          | moderate           |
| k_IAA_growth     | 0.0493014   |        0.150535  |            1.34319     | False          | well_estimated     |
| k_IAA_stationary | 0.0410836   |        0.0577506 |            1.11985     | False          | well_estimated     |
| k_EO_growth      | 0.000839698 |        0.124645  |            1.27673     | False          | well_estimated     |
| k_EO_stationary  | 0.000166548 |        0.16272   |            1.37566     | False          | well_estimated     |
| alpha_EA_loss    | 0.2         |        0.0960622 |            1.20717     | True           | weak_or_confounded |
| alpha_IAA_loss   | 0.765044    |        0.0447403 |            1.09165     | False          | well_estimated     |
| alpha_EO_loss    | 2.20532     |        0.0496831 |            1.10228     | False          | well_estimated     |

## Weak FIM directions

| analysis      |   weak_direction |   eigenvalue | dominant_parameters                                                                 | dominant_abs_loadings                                  |
|:--------------|-----------------:|-------------:|:------------------------------------------------------------------------------------|:-------------------------------------------------------|
| current_pilot |                1 |  3.46649e-14 | kAcStress, kAldRed, k_EA_growth, gammaF0, k_EA_stationary, iG, betaF0, k_IAA_growth | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| current_pilot |                2 |  0.00421626  | kAldRed, k_EA_growth, kAcAld, gammaF0, iG, betaF0, kAldS_N, qEF                     | 1.000, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.000 |
| current_pilot |                3 |  0.00774741  | k_EA_growth, k_EA_stationary, alpha_EA_loss, kAldRed, betaF0, gammaF0, iG, qEF      | 0.999, 0.041, 0.006, 0.002, 0.000, 0.000, 0.000, 0.000 |
| current_pilot |                4 |  0.247321    | gammaF0, gammaG0, betaF0, iG, qEF, iE, qEG, k_IAA_growth                            | 0.999, 0.041, 0.022, 0.015, 0.009, 0.008, 0.004, 0.002 |
| current_pilot |                5 |  6.82966     | k_EA_stationary, alpha_EA_loss, k_EA_growth, betaF0, iG, betaG0, qEF, qEG           | 0.990, 0.136, 0.040, 0.007, 0.004, 0.003, 0.001, 0.001 |
| current_pilot |                6 |  9.11086     | kPyrO2, kPyrDrain, kPyrS_N, iG, betaF0, qEF, iE, Kd0                                | 0.882, 0.443, 0.153, 0.033, 0.027, 0.016, 0.012, 0.008 |
| current_pilot |                7 | 12.0031      | betaF0, iG, qEF, iE, gammaG0, qEG, kPyrO2, k_IAA_growth                             | 0.651, 0.626, 0.327, 0.230, 0.108, 0.092, 0.041, 0.026 |
| current_pilot |                8 | 24.8073      | k_EO_stationary, k_EO_growth, alpha_EO_loss, iG, k_IAA_growth, betaF0, iE, qEF      | 0.794, 0.588, 0.139, 0.051, 0.028, 0.018, 0.014, 0.010 |

## Natural-must candidate ranking

| candidate                              | family          |   combined_logdet |   combined_min_relative_eigenvalue |   secondary_aroma_mean_var_reduction |   secondary_aroma_worst_var_reduction | N_pulses_kg_m3      | rationale                                                                                |
|:---------------------------------------|:----------------|------------------:|-----------------------------------:|-------------------------------------:|--------------------------------------:|:--------------------|:-----------------------------------------------------------------------------------------|
| natural_pilot_cold_to_warm_earlyN      | temperature_N   |           152.372 |                        3.98542e-08 |                             0.37509  |                             0.0731657 | 30h:0.045           | Cold start followed by warm transition and early nitrogen pulse.                         |
| natural_pilot_high_rate_strip          | co2_aroma       |           152.235 |                        4.0607e-08  |                             0.398745 |                             0.0948003 | 30h:0.03            | High-rate natural fermentation to excite CO2 stripping and aroma loss directions.        |
| natural_pilot_midN_temperature_step    | temperature_N   |           152.011 |                        3.99436e-08 |                             0.381691 |                             0.0776579 | 54h:0.04            | Temperature step with mid-growth nitrogen perturbation.                                  |
| natural_pilot_noN_dynamic_temperature  | temperature     |           151.5   |                        4.17699e-08 |                             0.417104 |                             0.103834  |                     | Temperature-only perturbation for settings where nutrient action is constrained.         |
| natural_pilot_two_step_N_ladder        | N_timing        |           151.423 |                        4.08467e-08 |                             0.36018  |                             0.0726209 | 30h:0.03; 78h:0.035 | Two smaller nitrogen pulses to separate early growth and later metabolic response.       |
| natural_pilot_reference_20C            | reference       |           151.34  |                        4.19355e-08 |                             0.415013 |                             0.10167   |                     | Natural must warmer reference to increase rate and CO2 information.                      |
| natural_pilot_low_temp_aroma_retention | aroma_retention |           151.173 |                        4.01933e-08 |                             0.374357 |                             0.0791033 | 54h:0.035           | Cold profile to contrast aroma retention against high-rate stripping.                    |
| natural_pilot_warm_to_cool_noN         | temperature     |           151.171 |                        4.25753e-08 |                             0.409833 |                             0.0774378 |                     | Warm early phase to excite growth and CO2, then cool aroma-retention phase.              |
| natural_pilot_lateN_stationary_probe   | N_timing        |           150.939 |                        4.11886e-08 |                             0.364919 |                             0.0711015 | 80h:0.05            | Late nitrogen addition to test stationary/growth split in secondary and aroma formation. |
| natural_pilot_reference_18C            | reference       |           150.872 |                        4.15166e-08 |                             0.408245 |                             0.0962585 |                     | Natural must reference at moderate temperature.                                          |

## Selected campaign, hybrid criterion

| candidate                              | family          |   campaign_order |   campaign_logdet |   campaign_min_relative_eigenvalue | N_pulses_kg_m3   | rationale                                                                         |
|:---------------------------------------|:----------------|-----------------:|------------------:|-----------------------------------:|:-----------------|:----------------------------------------------------------------------------------|
| natural_pilot_cold_to_warm_earlyN      | temperature_N   |                1 |           152.372 |                        3.98542e-08 | 30h:0.045        | Cold start followed by warm transition and early nitrogen pulse.                  |
| natural_pilot_high_rate_strip          | co2_aroma       |                2 |           159.157 |                        3.86574e-08 | 30h:0.03         | High-rate natural fermentation to excite CO2 stripping and aroma loss directions. |
| natural_pilot_low_temp_aroma_retention | aroma_retention |                3 |           163.76  |                        3.70571e-08 | 54h:0.035        | Cold profile to contrast aroma retention against high-rate stripping.             |
| natural_pilot_noN_dynamic_temperature  | temperature     |                4 |           167.326 |                        3.60169e-08 |                  | Temperature-only perturbation for settings where nutrient action is constrained.  |
| natural_pilot_midN_temperature_step    | temperature_N   |                5 |           170.368 |                        3.42174e-08 | 54h:0.04         | Temperature step with mid-growth nitrogen perturbation.                           |
| natural_pilot_warm_to_cool_noN         | temperature     |                6 |           172.907 |                        3.36376e-08 |                  | Warm early phase to excite growth and CO2, then cool aroma-retention phase.       |

## Interpretation notes

- `alpha_*_loss` parameters are now estimated because condenser-equivalent aroma observations are available.
- CO2 is used as a rate-shape signal with a per-batch scale factor, not as an absolute concentration measurement.
- Candidate designs do not include glucose/fructose/ethanol/biomass injections because this pilot setting is treated as natural-must constrained.
