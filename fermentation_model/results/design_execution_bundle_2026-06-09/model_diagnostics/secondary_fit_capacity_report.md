# Secondary fit-capacity diagnostic

## Purpose

This diagnostic asks whether the secondary-state data can be fitted at all when parameter estimability and parameter transferability are deliberately relaxed.

The comparison separates three explanations for poor global curves:

- If per-batch fits are good but global fits are poor, the equation form has enough local flexibility but the parameters are not transferable across batches/media.
- If per-medium fits improve but per-batch fits improve much more, medium-specific effects are present but do not explain all variability.
- If even per-batch fits are poor, the secondary-state structure, core forcing, or observations are inconsistent with the ODE form.

## Fit summary

| fit_label                                     | success   |   status | message                                                 |   nfev |   n_parameters_fit |   n_batches |   n_residuals |   initial_wsse |   final_wsse |   wsse_per_residual |
|:----------------------------------------------|:----------|---------:|:--------------------------------------------------------|-------:|-------------------:|------------:|--------------:|---------------:|-------------:|--------------------:|
| global_v2_reduced_o2fixed                     | True      |      nan | loaded from v2 reduced evaluation                       |    nan |                  7 |          27 |           588 |       nan      |  5754.09     |           9.78587   |
| medium_v2_reduced_o2fixed_historical_vl3      | True      |        2 | `ftol` termination condition is satisfied.              |     25 |                  7 |           8 |           251 |      2503.15   |  2418.63     |           9.63597   |
| medium_v2_reduced_o2fixed_natural             | True      |        2 | `ftol` termination condition is satisfied.              |      9 |                  7 |           9 |           253 |      2647.61   |  1949.82     |           7.70681   |
| medium_v2_reduced_o2fixed_synthetic           | True      |        2 | `ftol` termination condition is satisfied.              |     16 |                  7 |          10 |            87 |       603.339  |   365.77     |           4.20426   |
| batch_v2_reduced_o2fixed_historical_vl3_25026 | True      |        2 | `ftol` termination condition is satisfied.              |     44 |                  7 |           1 |            11 |       123.343  |     6.33316  |           0.575742  |
| batch_v2_all_o2free_historical_vl3_25026      | True      |        2 | `ftol` termination condition is satisfied.              |     46 |                 10 |           1 |            11 |       123.343  |     4.75219  |           0.432017  |
| batch_v2_reduced_o2fixed_historical_vl3_25027 | True      |        2 | `ftol` termination condition is satisfied.              |     22 |                  7 |           1 |            13 |       260.455  |    61.4148   |           4.72422   |
| batch_v2_all_o2free_historical_vl3_25027      | False     |        0 | The maximum number of function evaluations is exceeded. |     80 |                 10 |           1 |            13 |       260.455  |    32.8603   |           2.52771   |
| batch_v2_reduced_o2fixed_historical_vl3_25085 | True      |        2 | `ftol` termination condition is satisfied.              |     19 |                  7 |           1 |            34 |       153.521  |    26.8887   |           0.790846  |
| batch_v2_all_o2free_historical_vl3_25085      | True      |        2 | `ftol` termination condition is satisfied.              |     42 |                 10 |           1 |            34 |       153.521  |    25.7716   |           0.757988  |
| batch_v2_reduced_o2fixed_historical_vl3_25086 | True      |        2 | `ftol` termination condition is satisfied.              |     18 |                  7 |           1 |            34 |       242.33   |    33.4781   |           0.98465   |
| batch_v2_all_o2free_historical_vl3_25086      | True      |        2 | `ftol` termination condition is satisfied.              |     60 |                 10 |           1 |            34 |       242.33   |    21.0809   |           0.620025  |
| batch_v2_reduced_o2fixed_historical_vl3_25150 | True      |        2 | `ftol` termination condition is satisfied.              |     20 |                  7 |           1 |            27 |       458.673  |    37.2281   |           1.37882   |
| batch_v2_all_o2free_historical_vl3_25150      | True      |        2 | `ftol` termination condition is satisfied.              |     30 |                 10 |           1 |            27 |       458.673  |    26.702    |           0.988961  |
| batch_v2_reduced_o2fixed_historical_vl3_25151 | True      |        2 | `ftol` termination condition is satisfied.              |     15 |                  7 |           1 |            27 |       420.844  |    36.6519   |           1.35748   |
| batch_v2_all_o2free_historical_vl3_25151      | True      |        2 | `ftol` termination condition is satisfied.              |     18 |                 10 |           1 |            27 |       420.844  |    29.6363   |           1.09764   |
| batch_v2_reduced_o2fixed_historical_vl3_25170 | True      |        2 | `ftol` termination condition is satisfied.              |     47 |                  7 |           1 |            52 |       481.606  |   113.309    |           2.17901   |
| batch_v2_all_o2free_historical_vl3_25170      | True      |        2 | `ftol` termination condition is satisfied.              |     36 |                 10 |           1 |            52 |       481.606  |   108.754    |           2.09143   |
| batch_v2_reduced_o2fixed_historical_vl3_25171 | True      |        2 | `ftol` termination condition is satisfied.              |     17 |                  7 |           1 |            52 |       362.375  |   196.648    |           3.7817    |
| batch_v2_all_o2free_historical_vl3_25171      | True      |        2 | `ftol` termination condition is satisfied.              |     58 |                 10 |           1 |            52 |       362.375  |   183.695    |           3.5326    |
| batch_v2_reduced_o2fixed_natural_LAB004       | True      |        2 | `ftol` termination condition is satisfied.              |     22 |                  7 |           1 |            26 |       177.626  |    25.3529   |           0.975113  |
| batch_v2_all_o2free_natural_LAB004            | True      |        2 | `ftol` termination condition is satisfied.              |     49 |                 10 |           1 |            26 |       177.626  |    22.9799   |           0.883843  |
| batch_v2_reduced_o2fixed_natural_LAB005       | True      |        2 | `ftol` termination condition is satisfied.              |     19 |                  7 |           1 |            26 |       139.714  |    27.8493   |           1.07113   |
| batch_v2_all_o2free_natural_LAB005            | True      |        2 | `ftol` termination condition is satisfied.              |     63 |                 10 |           1 |            26 |       139.714  |    19.4424   |           0.747783  |
| batch_v2_reduced_o2fixed_natural_LAB006       | True      |        2 | `ftol` termination condition is satisfied.              |     20 |                  7 |           1 |            30 |       532.484  |   207.793    |           6.92643   |
| batch_v2_all_o2free_natural_LAB006            | True      |        2 | `ftol` termination condition is satisfied.              |     23 |                 10 |           1 |            30 |       532.484  |    91.7701   |           3.059     |
| batch_v2_reduced_o2fixed_natural_LAB007       | True      |        2 | `ftol` termination condition is satisfied.              |     19 |                  7 |           1 |            50 |       577.107  |   141.888    |           2.83776   |
| batch_v2_all_o2free_natural_LAB007            | True      |        2 | `ftol` termination condition is satisfied.              |     26 |                 10 |           1 |            50 |       577.107  |   124.248    |           2.48496   |
| batch_v2_reduced_o2fixed_natural_LAB008       | True      |        2 | `ftol` termination condition is satisfied.              |     19 |                  7 |           1 |            50 |       550.577  |   119.842    |           2.39684   |
| batch_v2_all_o2free_natural_LAB008            | True      |        2 | `ftol` termination condition is satisfied.              |     27 |                 10 |           1 |            50 |       550.577  |    99.6285   |           1.99257   |
| batch_v2_reduced_o2fixed_natural_LAB009       | True      |        2 | `ftol` termination condition is satisfied.              |     13 |                  7 |           1 |            46 |       507.559  |   176.788    |           3.84322   |
| batch_v2_all_o2free_natural_LAB009            | True      |        2 | `ftol` termination condition is satisfied.              |     20 |                 10 |           1 |            46 |       507.559  |   107.081    |           2.32784   |
| batch_v2_reduced_o2fixed_natural_LAB010       | True      |        2 | `ftol` termination condition is satisfied.              |     25 |                  7 |           1 |             8 |        84.1109 |     1.0148   |           0.12685   |
| batch_v2_all_o2free_natural_LAB010            | True      |        2 | `ftol` termination condition is satisfied.              |     40 |                 10 |           1 |             8 |        84.1109 |     0.962384 |           0.120298  |
| batch_v2_reduced_o2fixed_natural_LAB011       | True      |        2 | `ftol` termination condition is satisfied.              |     13 |                  7 |           1 |             8 |        57.6728 |     1.67007  |           0.208759  |
| batch_v2_all_o2free_natural_LAB011            | False     |        0 | The maximum number of function evaluations is exceeded. |     80 |                 10 |           1 |             8 |        57.6728 |     1.36024  |           0.170031  |
| batch_v2_reduced_o2fixed_natural_LAB012       | True      |        2 | `ftol` termination condition is satisfied.              |     56 |                  7 |           1 |             8 |        20.7566 |     1.04974  |           0.131218  |
| batch_v2_all_o2free_natural_LAB012            | True      |        3 | `xtol` termination condition is satisfied.              |     35 |                 10 |           1 |             8 |        20.7566 |     1.05275  |           0.131593  |
| batch_v2_reduced_o2fixed_synthetic_MS007      | False     |        0 | The maximum number of function evaluations is exceeded. |     80 |                  7 |           1 |             9 |        22.0375 |     1.9897   |           0.221078  |
| batch_v2_all_o2free_synthetic_MS007           | True      |        2 | `ftol` termination condition is satisfied.              |     43 |                 10 |           1 |             9 |        22.0375 |     0.875958 |           0.0973287 |
| batch_v2_reduced_o2fixed_synthetic_MS008      | True      |        2 | `ftol` termination condition is satisfied.              |     76 |                  7 |           1 |             9 |        18.6227 |     8.72045  |           0.968939  |
| batch_v2_all_o2free_synthetic_MS008           | True      |        2 | `ftol` termination condition is satisfied.              |     32 |                 10 |           1 |             9 |        18.6227 |     4.21334  |           0.468149  |
| batch_v2_reduced_o2fixed_synthetic_MS009      | True      |        2 | `ftol` termination condition is satisfied.              |     22 |                  7 |           1 |             9 |        81.1085 |    21.6624   |           2.40694   |
| batch_v2_all_o2free_synthetic_MS009           | True      |        2 | `ftol` termination condition is satisfied.              |     62 |                 10 |           1 |             9 |        81.1085 |     5.41141  |           0.601267  |
| batch_v2_reduced_o2fixed_synthetic_MS010      | True      |        2 | `ftol` termination condition is satisfied.              |     79 |                  7 |           1 |             9 |        41.0092 |    15.4001   |           1.71112   |
| batch_v2_all_o2free_synthetic_MS010           | True      |        2 | `ftol` termination condition is satisfied.              |     28 |                 10 |           1 |             9 |        41.0092 |    10.8823   |           1.20914   |
| batch_v2_reduced_o2fixed_synthetic_MS011      | True      |        2 | `ftol` termination condition is satisfied.              |     18 |                  7 |           1 |             9 |       114.141  |    16.6873   |           1.85414   |
| batch_v2_all_o2free_synthetic_MS011           | True      |        2 | `ftol` termination condition is satisfied.              |     29 |                 10 |           1 |             9 |       114.141  |    10.405    |           1.15611   |
| batch_v2_reduced_o2fixed_synthetic_MS012      | True      |        2 | `ftol` termination condition is satisfied.              |     17 |                  7 |           1 |             9 |        96.3182 |    11.0952   |           1.2328    |
| batch_v2_all_o2free_synthetic_MS012           | True      |        2 | `ftol` termination condition is satisfied.              |     15 |                 10 |           1 |             9 |        96.3182 |    10.4527   |           1.16141   |
| batch_v2_reduced_o2fixed_synthetic_MS013      | False     |        0 | The maximum number of function evaluations is exceeded. |     80 |                  7 |           1 |             8 |        39.1279 |     5.94224  |           0.742779  |
| batch_v2_all_o2free_synthetic_MS013           | True      |        2 | `ftol` termination condition is satisfied.              |     13 |                 10 |           1 |             8 |        39.1279 |     2.8904   |           0.3613    |
| batch_v2_reduced_o2fixed_synthetic_MS014      | True      |        2 | `ftol` termination condition is satisfied.              |     72 |                  7 |           1 |             8 |        46.9589 |    13.8903   |           1.73629   |
| batch_v2_all_o2free_synthetic_MS014           | True      |        2 | `ftol` termination condition is satisfied.              |     30 |                 10 |           1 |             8 |        46.9589 |     2.7081   |           0.338513  |
| batch_v2_reduced_o2fixed_synthetic_MS015      | True      |        2 | `ftol` termination condition is satisfied.              |     13 |                  7 |           1 |             8 |        30.5565 |    12.1959   |           1.52449   |
| batch_v2_all_o2free_synthetic_MS015           | True      |        2 | `ftol` termination condition is satisfied.              |     22 |                 10 |           1 |             8 |        30.5565 |     6.99283  |           0.874104  |
| batch_v2_reduced_o2fixed_synthetic_MS016      | True      |        2 | `ftol` termination condition is satisfied.              |     67 |                  7 |           1 |             8 |       113.459  |    25.958    |           3.24475   |
| batch_v2_all_o2free_synthetic_MS016           | False     |        0 | The maximum number of function evaluations is exceeded. |     80 |                 10 |           1 |             8 |       113.459  |     1.28422  |           0.160528  |

## Aggregated state metrics

| fit_label                 | state   |   n_obs |      wsse |   wsse_per_obs |   rmse_weighted_mean |   mae_weighted_mean |
|:--------------------------|:--------|--------:|----------:|---------------:|---------------------:|--------------------:|
| batch_v2_all_o2free       | AcAld   |     193 |  425.565  |        2.205   |           11.2486    |           7.897     |
| batch_v2_all_o2free       | Acetate |      50 |   59.1174 |        1.18235 |            0.0360319 |           0.0277258 |
| batch_v2_all_o2free       | O2      |      32 |  110.458  |        3.45182 |            0.370944  |           0.320975  |
| batch_v2_all_o2free       | Pyr     |     313 |  362.752  |        1.15895 |            5.89605   |           4.5402    |
| batch_v2_reduced_o2fixed  | AcAld   |     193 |  685.19   |        3.55021 |           13.6912    |           9.94555   |
| batch_v2_reduced_o2fixed  | Acetate |      50 |   99.2743 |        1.98549 |            0.0460072 |           0.035426  |
| batch_v2_reduced_o2fixed  | O2      |      32 |   98.1294 |        3.06654 |            0.345615  |           0.297618  |
| batch_v2_reduced_o2fixed  | Pyr     |     313 |  466.147  |        1.48929 |            6.84414   |           5.25661   |
| global_v2_reduced_o2fixed | AcAld   |     193 | 2686.89   |       13.9217  |           28.4624    |          23.6152    |
| global_v2_reduced_o2fixed | Acetate |      50 |  677.116  |       13.5423  |            0.127543  |           0.10496   |
| global_v2_reduced_o2fixed | O2      |      32 |   98.1294 |        3.06654 |            0.345615  |           0.297618  |
| global_v2_reduced_o2fixed | Pyr     |     313 | 2291.95   |        7.32254 |           15.4318    |          12.1285    |
| medium_v2_reduced_o2fixed | AcAld   |     193 | 2308.03   |       11.9587  |           25.9157    |          21.2225    |
| medium_v2_reduced_o2fixed | Acetate |      50 |  675.396  |       13.5079  |            0.127667  |           0.105373  |
| medium_v2_reduced_o2fixed | O2      |      32 |   98.1294 |        3.06654 |            0.345615  |           0.297618  |
| medium_v2_reduced_o2fixed | Pyr     |     313 | 1645.35   |        5.25672 |           12.7009    |           9.76718   |
| batch_v2_all_o2free       | ALL     |     588 |  957.893  |        1.62907 |          nan         |         nan         |
| batch_v2_reduced_o2fixed  | ALL     |     588 | 1348.74   |        2.29378 |          nan         |         nan         |
| global_v2_reduced_o2fixed | ALL     |     588 | 5754.09   |        9.78587 |          nan         |         nan         |
| medium_v2_reduced_o2fixed | ALL     |     588 | 4726.91   |        8.03896 |          nan         |         nan         |

## Parameter spread in per-batch all-free fits

| fit_label           | parameter   |   n |         min |      median |       max |   fold_range |   robust_cv |
|:--------------------|:------------|----:|------------:|------------:|----------:|-------------:|------------:|
| batch_v2_all_o2free | kPyrS_N     |  27 | 0.0001      | 0.00529661  |  2.43854  |    24385.4   |  195.061    |
| batch_v2_all_o2free | kPyrO2      |  27 | 0.0001      | 1.18057     |  9.98283  |    99828.3   |    1.52214  |
| batch_v2_all_o2free | kPyrDrain   |  27 | 0.00101903  | 0.00735251  |  0.148551 |      145.777 |    1.17107  |
| batch_v2_all_o2free | kAldS_N     |  27 | 0.000238777 | 2.16841     | 10        |    41880.1   |    1.47747  |
| batch_v2_all_o2free | kAldRed     |  27 | 1e-05       | 2.78464e-05 |  1.97695  |   197695     |   23.3802   |
| batch_v2_all_o2free | kAcAld      |  27 | 4.26744e-05 | 0.00411502  |  1.69359  |    39686.3   |   13.4771   |
| batch_v2_all_o2free | kAcStress   |  27 | 1e-06       | 0.00261587  |  0.468247 |   468247     |    0.565412 |
| batch_v2_all_o2free | qO2         |  27 | 0.000103869 | 0.0989653   |  0.359725 |     3463.25  |    0.938347 |
| batch_v2_all_o2free | kLaO2       |  27 | 1e-05       | 0.000590423 |  0.10834  |    10834     |   14.1681   |
| batch_v2_all_o2free | O2sat       |  27 | 0.5         | 2.28436     |  8        |       16     |    3.25545  |

## Interpretation guide

Treat per-batch all-free fits as an upper bound on fit capacity, not as a calibratable model. Large parameter fold-ranges mean the data can be matched only by sacrificing parameter transferability and identifiability.