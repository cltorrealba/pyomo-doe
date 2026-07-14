# Estimability: historical data vs historical + Lot 1

## What was compared

Two information sets were evaluated with the same model and parameter set:

1. `historical_only`: the normalized historical natural/synthetic must database.
2. `historical_plus_lot1`: historical data plus the executed Lot 1 batches reconstructed from `fermentation_lot1_data_preview.executed.ipynb`.

The extended FIM uses observed core states and secondary states where available. Lot 1 online CO2 is included as an information-time block because the current model state is cumulative CO2, whereas the sensor reports gas flow; this should be interpreted as information availability, not as a final calibrated CO2 likelihood.

## FIM Metrics

|    logdet |   min_eigenvalue |   max_eigenvalue |   min_relative_eigenvalue |   condition_number |        trace_inv |   rank_1e-8 | case                 | block                                  |   n_parameters |
|----------:|-----------------:|-----------------:|--------------------------:|-------------------:|-----------------:|------------:|:---------------------|:---------------------------------------|---------------:|
|   76.6992 |       0.00119233 |           153664 |               7.75931e-09 |        1.28877e+08 |    820.826       |          10 | historical_only      | core_profile_parameters                |             11 |
| -158.628  |       0          |           153664 |               0           |        1e+12       |      1.36501e+06 |          10 | historical_only      | extended_observed_plus_co2_information |             26 |
|   82.6987 |       0.0158271  |           209168 |               7.5667e-08  |        1.32158e+07 |     63.0557      |          11 | historical_plus_lot1 | core_profile_parameters                |             11 |
|  -21.869  |       0          |           213252 |               0           |        1e+12       | 574054           |          17 | historical_plus_lot1 | extended_observed_plus_co2_information |             26 |

## Core Fit Candidate Audit

| case                 | seed                  |   final_wsse |   nfev | success   |
|:---------------------|:----------------------|-------------:|-------:|:----------|
| historical_only      | seed_default          |      23882.8 |      7 | True      |
| historical_only      | qN_up_25pct           |      23859.6 |      8 | True      |
| historical_only      | Kd0_up_2x             |      23895.4 |      8 | True      |
| historical_only      | qN_up_25pct_Kd0_up_2x |      23901.7 |     10 | True      |
| historical_plus_lot1 | seed_default          |      39076.6 |      7 | True      |
| historical_plus_lot1 | qN_up_25pct           |      36256.5 |      9 | True      |
| historical_plus_lot1 | Kd0_up_2x             |      38977.7 |     10 | True      |
| historical_plus_lot1 | qN_up_25pct_Kd0_up_2x |      36447.5 |      9 | True      |

## Largest Approximate Estimability Improvements

Negative `std_log_delta_plus_minus_old` means the parameter became more estimable after adding Lot 1.

| parameter        |   historical_only |   historical_plus_lot1 |   std_log_delta_plus_minus_old |   std_log_ratio_plus_over_old |
|:-----------------|------------------:|-----------------------:|-------------------------------:|------------------------------:|
| kAcStress        |           301.571 |              0.0296597 |                      -301.541  |                   9.83508e-05 |
| kAldS_N          |           301.571 |              0.0308341 |                      -301.54   |                   0.000102245 |
| kPyrS_N          |           301.571 |              0.0694011 |                      -301.501  |                   0.000230132 |
| kPyrDrain        |           301.571 |              0.162318  |                      -301.409  |                   0.000538242 |
| kPyrO2           |           301.571 |              0.166151  |                      -301.405  |                   0.000550953 |
| kAcAld           |           301.571 |              0.182029  |                      -301.389  |                   0.000603604 |
| k_EA_stationary  |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |
| k_EA_growth      |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |
| k_IAA_growth     |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |
| k_EO_stationary  |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |
| k_IAA_stationary |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |
| alpha_IAA_loss   |           301.571 |            252.541     |                       -49.0304 |                   0.837417    |

## Profile Likelihood Summary

| case                 | parameter   |   n_success |   max_lr_stat | crosses_left_95   | crosses_right_95   | profile_identifiable_95   |   chi2_95_threshold |
|:---------------------|:------------|------------:|--------------:|:------------------|:-------------------|:--------------------------|--------------------:|
| historical_only      | Kd0         |           4 |       462.483 | True              | True               | True                      |             3.84146 |
| historical_only      | qN          |           4 |      5131.63  | True              | True               | True                      |             3.84146 |
| historical_plus_lot1 | Kd0         |           4 |       418.092 | True              | True               | True                      |             3.84146 |
| historical_plus_lot1 | qN          |           4 |     13367.3   | True              | True               | True                      |             3.84146 |

## Weak Directions After Adding Lot 1

| case                 |   weak_direction |   eigenvalue | dominant_parameters                                               | dominant_abs_loadings                           |
|:---------------------|-----------------:|-------------:|:------------------------------------------------------------------|:------------------------------------------------|
| historical_plus_lot1 |                1 |            0 | k_EA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_plus_lot1 |                2 |            0 | k_EA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG  | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_plus_lot1 |                3 |            0 | k_IAA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG     | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_plus_lot1 |                4 |            0 | k_IAA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_plus_lot1 |                5 |            0 | k_EO_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_plus_lot1 |                6 |            0 | k_EO_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG  | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |

## Generated figures

- `figures/estimability_std_log_comparison.png`
- `figures/profile_likelihood_historical_only.png`
- `figures/profile_likelihood_historical_plus_lot1.png`
