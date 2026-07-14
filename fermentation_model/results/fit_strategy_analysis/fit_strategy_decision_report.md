# Fit strategy analysis

Recommended regularized curve-fit candidate: `plus_qx_iG_logL2_10`.

Best log-L2 candidate with no active estimated-parameter bounds; its penalized objective is still below the reduced-model WSSE.

- Best free fit: `plus_qx_iG_unregularized`.
- Best stable direct-WSSE fit: `plus_qx_iG_logL2_10`.
- Best stable log-L2 fit: `plus_qx_iG_logL2_10`.

Selection metric for the regularized candidate: direct weighted SSE on the data plus the explicit log-L2 penalty.
Fits with active estimated-parameter bounds are reported as sensitivity fits, not as recommended estimable models.

## Recommended follow-up profile set

Run profile-likelihood diagnostics for `plus_qx_iG_logL2_10` and compare it against `reduced6_local`.
If the released parameters do not cross the 95% profile threshold without relying on the prior penalty, keep them fixed for DOE and use this candidate only as a Bayesian regularized curve-fit scenario.

Recommended candidate metrics:

|                            | value                                           |
|:---------------------------|:------------------------------------------------|
| estimated_parameters       | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |
| n_estimated                | 9                                               |
| l2_lambda                  | 10.0                                            |
| l2_parameters              | qXG, qXF, iG                                    |
| uses_pso                   | False                                           |
| pso_seed                   | nan                                             |
| ParmEst_objective_total    | 5998.0545004005935                              |
| direct_WSSE_total          | 5992.323835014267                               |
| direct_SSE_total           | 18518.288322120992                              |
| l2_penalty_total           | 5.730665694336276                               |
| direct_plus_l2_total       | 5998.054500708604                               |
| n_observations             | 286                                             |
| n_active_bounds            | 0                                               |
| active_bounds              | nan                                             |
| runtime_s                  | 5.758035399951041                               |
| solve_ok                   | True                                            |
| pso_best_direct_WSSE       | nan                                             |
| pso_evaluations            | nan                                             |
| pso_successful_evaluations | nan                                             |
| pso_failed_evaluations     | nan                                             |

## Top strategies

| strategy                      | estimated_parameters                            |   n_estimated |   l2_lambda | l2_parameters   | uses_pso   |   pso_seed |   ParmEst_objective_total |   direct_WSSE_total |   direct_SSE_total |   l2_penalty_total |   direct_plus_l2_total |   n_observations |   n_active_bounds | active_bounds        |   runtime_s | solve_ok   |   pso_best_direct_WSSE |   pso_evaluations |   pso_successful_evaluations |   pso_failed_evaluations |
|:------------------------------|:------------------------------------------------|--------------:|------------:|:----------------|:-----------|-----------:|--------------------------:|--------------------:|-------------------:|-------------------:|-----------------------:|-----------------:|------------------:|:---------------------|------------:|:-----------|-----------------------:|------------------:|-----------------------------:|-------------------------:|
| plus_qx_iG_unregularized      | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |             9 |           0 | nan             | False      |        nan |                   5972.55 |             5972.55 |            18673   |            0       |                5972.55 |              286 |                 2 | qXG:lower; qXF:lower |     4.95739 | True       |                 nan    |               nan |                          nan |                      nan |
| plus_qx_unregularized         | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF     |             8 |           0 | nan             | False      |        nan |                   5983.45 |             5983.45 |            18898.3 |            0       |                5983.45 |              286 |                 2 | qXG:lower; qXF:lower |     4.3847  | True       |                 nan    |               nan |                          nan |                      nan |
| plus_qx_iG_logL2_10           | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |             9 |          10 | qXG, qXF, iG    | False      |        nan |                   5998.05 |             5992.32 |            18518.3 |            5.73067 |                5998.05 |              286 |                 0 | nan                  |     5.75804 | True       |                 nan    |               nan |                          nan |                      nan |
| plus_iG_unregularized         | mu0, qN, betaG0, betaF0, qEG, qEF, iG           |             7 |           0 | nan             | False      |        nan |                   6001.98 |             6001.98 |            18459.3 |            0       |                6001.98 |              286 |                 0 | nan                  |     4.67737 | True       |                 nan    |               nan |                          nan |                      nan |
| plus_qx_iG_pso_pilot_logL2_50 | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |             9 |          50 | qXG, qXF, iG    | True       |        321 |                   6005.65 |             6002.1  |            18544.8 |            3.54467 |                6005.65 |              286 |                 0 | nan                  |   340.745   | True       |                6010.75 |                71 |                           61 |                       10 |
| plus_qx_iG_logL2_50           | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |             9 |          50 | qXG, qXF, iG    | False      |        nan |                   6005.65 |             6002.1  |            18544.9 |            3.54759 |                6005.65 |              286 |                 0 | nan                  |     5.4753  | True       |                 nan    |               nan |                          nan |                      nan |
| plus_iG_logL2_50              | mu0, qN, betaG0, betaF0, qEG, qEF, iG           |             7 |          50 | iG              | False      |        nan |                   6007.18 |             6005.06 |            18538.1 |            2.11685 |                6007.18 |              286 |                 0 | nan                  |     4.73414 | True       |                 nan    |               nan |                          nan |                      nan |
| plus_qx_iG_logL2_100          | mu0, qN, betaG0, betaF0, qEG, qEF, qXG, qXF, iG |             9 |         100 | qXG, qXF, iG    | False      |        nan |                   6007.72 |             6005.28 |            18565   |            2.43551 |                6007.72 |              286 |                 0 | nan                  |     5.46099 | True       |                 nan    |               nan |                          nan |                      nan |

## State residuals for selected strategies

|                                   |   SSE_raw_sum |   WSSE_raw_sum |
|:----------------------------------|--------------:|---------------:|
| ('reduced6_local', 'E')           |  7625.57      |       4287.64  |
| ('reduced6_local', 'F')           |  5305.82      |        174.53  |
| ('reduced6_local', 'G')           |  5536.92      |        211.369 |
| ('reduced6_local', 'N')           |     0.0929371 |       1063.11  |
| ('reduced6_local', 'X')           |   137.05      |        274.099 |
| ('plus_qx_iG_unregularized', 'E') |  7567.15      |       4254.79  |
| ('plus_qx_iG_unregularized', 'F') |  5464.05      |        179.735 |
| ('plus_qx_iG_unregularized', 'G') |  5509.98      |        210.34  |
| ('plus_qx_iG_unregularized', 'N') |     0.0930418 |       1064.31  |
| ('plus_qx_iG_unregularized', 'X') |   131.686     |        263.371 |
| ('plus_qx_iG_logL2_10', 'E')      |  7602.7       |       4274.78  |
| ('plus_qx_iG_logL2_10', 'F')      |  5263.4       |        173.135 |
| ('plus_qx_iG_logL2_10', 'G')      |  5516.99      |        210.608 |
| ('plus_qx_iG_logL2_10', 'N')      |     0.0929785 |       1063.58  |
| ('plus_qx_iG_logL2_10', 'X')      |   135.106     |        270.213 |
