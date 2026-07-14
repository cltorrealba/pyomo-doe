# Deep model selection report

Decision rule:

1. Reject candidates with failed solves or estimated parameters at active bounds as structural models.
2. Among stable candidates, require full-rank weighted-relative FIM/Q diagnostics.
3. Require two-sided 95% profile-likelihood crossings for practical identifiability.
4. Use log-L2/Bayesian candidates only as regularized curve-fit scenarios unless their added parameters pass unpenalized profiles.

## Fit and information criteria

| candidate                     |   n_estimated |   direct_WSSE_total |   AIC_like |   BIC_like |   n_active_bounds | solve_ok   |
|:------------------------------|--------------:|--------------------:|-----------:|-----------:|------------------:|:-----------|
| plus_qx_iG_unregularized      |             9 |             5972.55 |    11963.1 |    11996   |                 2 | True       |
| plus_qx_unregularized         |             8 |             5983.45 |    11982.9 |    12012.2 |                 2 | True       |
| plus_qx_iG_logL2_10           |             9 |             5992.32 |    12002.6 |    12035.6 |                 0 | True       |
| plus_iG_unregularized         |             7 |             6001.98 |    12018   |    12043.6 |                 0 | True       |
| plus_iG_logL2_50              |             7 |             6005.06 |    12024.1 |    12049.7 |                 0 | True       |
| plus_qx_iG_logL2_50           |             9 |             6002.1  |    12022.2 |    12055.1 |                 0 | True       |
| plus_qx_iG_pso_pilot_logL2_50 |             9 |             6002.1  |    12022.2 |    12055.1 |                 0 | True       |
| reduced6_local                |             6 |             6010.75 |    12033.5 |    12055.4 |                 0 | True       |
| plus_qx_logL2_50              |             8 |             6008.09 |    12032.2 |    12061.4 |                 0 | True       |
| plus_qx_iG_logL2_100          |             9 |             6005.28 |    12028.6 |    12061.5 |                 0 | True       |

## FIM/Q sensitivity summary

| candidate                |   n_requested_parameters |   n_used_parameters | used_parameters                                 |   dropped_parameters |   fim_condition_number |   fim_min_relative_eigenvalue |   fim_near_null_directions | fim_full_rank_by_threshold   | weakest_direction_parameters   | weakest_direction_abs_loadings    |   q_failed_solves |
|:-------------------------|-------------------------:|--------------------:|:------------------------------------------------|---------------------:|-----------------------:|------------------------------:|---------------------------:|:-----------------------------|:-------------------------------|:----------------------------------|------------------:|
| reduced6_local           |                        6 |                   6 | mu0, qN, betaG0, betaF0, qEG, qEF               |                  nan |          1118.31       |                   0.000894208 |                          0 | True                         | betaF0, qEF, qEG, qN, mu0      | 0.839, 0.413, 0.305, 0.138, 0.113 |                 0 |
| plus_iG_unregularized    |                        7 |                   7 | mu0, qN, betaG0, betaF0, qEG, qEF, iG           |                  nan |          8381.22       |                   0.000119314 |                          0 | True                         | iG, qEF, betaF0, betaG0, qEG   | 0.902, 0.333, 0.201, 0.166, 0.082 |                 0 |
| plus_qx_iG_unregularized |                        9 |                   9 | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |                  nan |             4.8454e+07 |                   9.06169e-09 |                          1 | False                        | qXF, qXG, iG, betaG0, qEF      | 0.929, 0.370, 0.013, 0.003, 0.002 |                 0 |
| plus_qx_iG_logL2_10      |                        9 |                   9 | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |                  nan |        370243          |                   2.70093e-06 |                          0 | True                         | qXF, qXG, iG, qEF, betaG0      | 0.806, 0.554, 0.197, 0.039, 0.037 |                 0 |

## Profile-likelihood summary

|    | profiled_theta   |   n_success |   theta_hat |   min_profile_objective |   max_lr_stat | crosses_left   | crosses_right   |   chi2_threshold | candidate                | objective_kind        | profile_identifiable_2sided   |
|---:|:-----------------|------------:|------------:|------------------------:|--------------:|:---------------|:----------------|-----------------:|:-------------------------|:----------------------|:------------------------------|
|  0 | betaF0           |           9 |   0.283496  |                 6001.98 |       512.425 | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  1 | betaG0           |           9 |   0.39504   |                 6001.98 |       846.987 | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  2 | iG               |           8 |   0.0188579 |                 6001.98 |       239.484 | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  3 | mu0              |           9 |   0.15683   |                 6001.98 |      1012.74  | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  4 | qEF              |          10 |   0.502258  |                 6001.98 |    231940     | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  5 | qEG              |           9 |   0.513865  |                 6001.98 |    213054     | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  6 | qN               |          10 |   0.0134154 |                 6001.98 |    169809     | True           | True            |          3.84146 | plus_iG_unregularized    | SSE_weighted          | True                          |
|  7 | iG               |           8 |   0.0174101 |                 5972.55 |       351.714 | True           | True            |          3.84146 | plus_qx_iG_unregularized | SSE_weighted          | True                          |
|  8 | qXF              |           8 |   0.005     |                 5972.55 |       115.756 | False          | True            |          3.84146 | plus_qx_iG_unregularized | SSE_weighted          | False                         |
|  9 | qXG              |           8 |   0.005     |                 5972.55 |       224.618 | False          | True            |          3.84146 | plus_qx_iG_unregularized | SSE_weighted          | False                         |
| 10 | iG               |           8 |   0.0201897 |                 5998.05 |       309.999 | True           | True            |          3.84146 | plus_qx_iG_logL2_10      | SSE_weighted_logL2_10 | True                          |
| 11 | qXF              |          10 |   0.0916033 |                 5998.05 |      1350.03  | True           | True            |          3.84146 | plus_qx_iG_logL2_10      | SSE_weighted_logL2_10 | True                          |
| 12 | qXG              |          10 |   0.0755015 |                 5998.05 |       352.373 | True           | True            |          3.84146 | plus_qx_iG_logL2_10      | SSE_weighted_logL2_10 | True                          |

## Bayesian/Laplace summary

|        | candidate                |   theta_hat |   reference_prior_median |   posterior_median_laplace |   posterior_log_sd |   posterior_95_lower |   posterior_95_upper |   posterior_var_over_weak_prior_var |   l2_prior_precision | data_dominated_vs_weak_prior   |
|:-------|:-------------------------|------------:|-------------------------:|---------------------------:|-------------------:|---------------------:|---------------------:|------------------------------------:|---------------------:|:-------------------------------|
| mu0    | reduced6_local           |   0.15886   |                0.158859  |                  0.15886   |          0.0278687 |            0.150415  |            0.167778  |                         0.000776666 |                    0 | True                           |
| qN     | reduced6_local           |   0.0135462 |                0.0135461 |                  0.0135462 |          0.0262804 |            0.0128661 |            0.0142622 |                         0.000690662 |                    0 | True                           |
| betaG0 | reduced6_local           |   0.409581  |                0.409583  |                  0.409581  |          0.0425131 |            0.376835  |            0.445171  |                         0.00180736  |                    0 | True                           |
| betaF0 | reduced6_local           |   0.309421  |                0.309414  |                  0.309421  |          0.0737489 |            0.267777  |            0.357541  |                         0.0054389   |                    0 | True                           |
| qEG    | reduced6_local           |   0.51583   |                0.515827  |                  0.51583   |          0.0369666 |            0.479778  |            0.554592  |                         0.00136653  |                    0 | True                           |
| qEF    | reduced6_local           |   0.579228  |                0.579225  |                  0.579228  |          0.042184  |            0.533263  |            0.629154  |                         0.00177949  |                    0 | True                           |
| mu0    | plus_iG_unregularized    |   0.15683   |                0.158859  |                  0.156792  |          0.0283548 |            0.148316  |            0.165752  |                         0.000803994 |                    0 | True                           |
| qN     | plus_iG_unregularized    |   0.0134154 |                0.0135461 |                  0.0134203 |          0.0266056 |            0.0127384 |            0.0141387 |                         0.000707856 |                    0 | True                           |
| betaG0 | plus_iG_unregularized    |   0.39504   |                0.409583  |                  0.396357  |          0.0586738 |            0.353299  |            0.444663  |                         0.00344261  |                    0 | True                           |
| betaF0 | plus_iG_unregularized    |   0.283496  |                0.309414  |                  0.28469   |          0.0836652 |            0.241632  |            0.33542   |                         0.00699986  |                    0 | True                           |
| qEG    | plus_iG_unregularized    |   0.513865  |                0.515827  |                  0.514708  |          0.0393912 |            0.476464  |            0.556021  |                         0.00155167  |                    0 | True                           |
| qEF    | plus_iG_unregularized    |   0.502258  |                0.579225  |                  0.505719  |          0.0823714 |            0.430322  |            0.594327  |                         0.00678505  |                    0 | True                           |
| iG     | plus_iG_unregularized    |   0.0188579 |                0.0270362 |                  0.0192124 |          0.207557  |            0.012791  |            0.0288574 |                         0.04308     |                    0 | True                           |
| mu0    | plus_qx_iG_unregularized |   0.150517  |                0.158859  |                  0.150464  |          0.0294807 |            0.142016  |            0.159414  |                         0.000869113 |                    0 | True                           |
| qN     | plus_qx_iG_unregularized |   0.0130038 |                0.0135461 |                  0.0129902 |          0.0272843 |            0.0123138 |            0.0137038 |                         0.000744433 |                    0 | True                           |
| betaG0 | plus_qx_iG_unregularized |   0.39507   |                0.409583  |                  0.401544  |          0.059307  |            0.357479  |            0.451042  |                         0.00351733  |                    0 | True                           |
| betaF0 | plus_qx_iG_unregularized |   0.287272  |                0.309414  |                  0.287357  |          0.0819506 |            0.244717  |            0.337427  |                         0.0067159   |                    0 | True                           |
| qEG    | plus_qx_iG_unregularized |   0.539325  |                0.515827  |                  0.53844   |          0.0395427 |            0.498285  |            0.581831  |                         0.00156363  |                    0 | True                           |
| qEF    | plus_qx_iG_unregularized |   0.517297  |                0.579225  |                  0.524168  |          0.0796054 |            0.448445  |            0.612677  |                         0.00633701  |                    0 | True                           |
| qXG    | plus_qx_iG_unregularized |   0.005     |                0.1125    |                  0.111312  |          0.99852   |            0.0157248 |            0.787951  |                         0.997042    |                    0 | False                          |
| qXF    | plus_qx_iG_unregularized |   0.005     |                0.1125    |                  0.112334  |          0.999072  |            0.015852  |            0.796047  |                         0.998144    |                    0 | False                          |
| iG     | plus_qx_iG_unregularized |   0.0174101 |                0.0270362 |                  0.0186191 |          0.20459   |            0.0124683 |            0.0278041 |                         0.0418569   |                    0 | True                           |
| mu0    | plus_qx_iG_logL2_10      |   0.155378  |                0.158859  |                  0.15524   |          0.0286291 |            0.146769  |            0.1642    |                         0.000819623 |                    0 | True                           |
| qN     | plus_qx_iG_logL2_10      |   0.0133101 |                0.0135461 |                  0.0133046 |          0.0267056 |            0.0126261 |            0.0140195 |                         0.000713191 |                    0 | True                           |
| betaG0 | plus_qx_iG_logL2_10      |   0.39958   |                0.409583  |                  0.417116  |          0.0500849 |            0.378115  |            0.46014   |                         0.00250849  |                    0 | True                           |
| betaF0 | plus_qx_iG_logL2_10      |   0.288925  |                0.309414  |                  0.295898  |          0.0780823 |            0.253908  |            0.344831  |                         0.00609685  |                    0 | True                           |
| qEG    | plus_qx_iG_logL2_10      |   0.523043  |                0.515827  |                  0.524245  |          0.0378714 |            0.48674   |            0.564639  |                         0.00143424  |                    0 | True                           |
| qEF    | plus_qx_iG_logL2_10      |   0.520061  |                0.579225  |                  0.554079  |          0.0609246 |            0.491713  |            0.624355  |                         0.0037118   |                    0 | True                           |
| qXG    | plus_qx_iG_logL2_10      |   0.0755015 |                0.1125    |                  0.111911  |          0.154878  |            0.0826108 |            0.151603  |                         0.0239873   |                   40 | True                           |
| qXF    | plus_qx_iG_logL2_10      |   0.0916033 |                0.1125    |                  0.1143    |          0.153304  |            0.084635  |            0.154362  |                         0.023502    |                   40 | True                           |
| iG     | plus_qx_iG_logL2_10      |   0.0201897 |                0.0270362 |                  0.0249596 |          0.125676  |            0.0195101 |            0.0319312 |                         0.0157943   |                   40 | True                           |

## Current recommendation

Use `plus_iG_unregularized` as the recommended structurally/practically identifiable model.
`plus_iG_unregularized` improves fit relative to `reduced6_local`, remains full-rank in the weighted-relative FIM/Q diagnostic, and passes two-sided profile-likelihood checks for all seven estimated parameters.
Use `plus_qx_iG_logL2_10` only as a Bayesian regularized curve-fit sensitivity candidate; `qXG` and `qXF` do not pass the unregularized structural/practical identifiability filters.