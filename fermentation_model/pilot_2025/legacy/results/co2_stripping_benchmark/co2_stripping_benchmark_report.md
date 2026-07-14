# CO2 stripping and sugar-aware partition benchmark

## CO2 model ranking

| mode          |   data_wsse |   n_data_residuals |   n_parameters |   active_bound_count |     bic | params                                                                 | description                                                                                                 | selected_by_bic   |
|:--------------|------------:|-------------------:|---------------:|---------------------:|--------:|:-----------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------|:------------------|
| lag_threshold |     222.432 |                308 |              2 |                    0 | 233.892 | {'tau_h': 36.38236217778723, 'threshold_fraction': 0.8419628866863644} | Thresholded production-rate transform followed by first-order gas-release lag.                              | True              |
| lag_power     |     226.758 |                308 |              2 |                    1 | 238.219 | {'tau_h': 61.46629257205183, 'power_p': 2.499999999998967}             | Nonlinear production-rate transform followed by first-order gas-release lag.                                | False             |
| threshold     |     244.736 |                308 |              1 |                    0 | 250.466 | {'threshold_fraction': 0.49809756401738725}                            | Gas flow activates only above a fraction of the median production rate.                                     | False             |
| lag           |     246.625 |                308 |              1 |                    0 | 252.355 | {'tau_h': 25.853734739141235}                                          | First-order gas-release/headspace lag: tau dq/dt = q_prod - q_out.                                          | False             |
| instant       |     258.528 |                308 |              0 |                    0 | 264.258 | {}                                                                     | Current model: stripping gas flow is directly proportional to instantaneous ethanol-derived CO2 production. | False             |
| power         |     259.298 |                308 |              1 |                    0 | 265.028 | {'power_p': 1.0945677001893517}                                        | Nonlinear effective gas flow: q_out = q_ref (q_prod/q_ref)^p.                                               | False             |

## CO2 metrics by batch

| mode          |   batch |   n |   scale_factor |     rmse |       mae |        bias |   relative_rmse |   relative_bias |     corr |
|:--------------|--------:|----:|---------------:|---------:|----------:|------------:|----------------:|----------------:|---------:|
| instant       |   25170 | 191 |      1.93952   | 0.383478 | 0.333546  | -0.173784   |        1.24288  |     -0.563244   | 0.328081 |
| instant       |   25171 | 117 |      4.02669   | 0.174973 | 0.140931  |  0.0543823  |        0.391057 |      0.121542   | 0.959637 |
| lag           |   25170 | 191 |      2.35844   | 0.301231 | 0.26328   | -0.0968018  |        0.976307 |     -0.313741   | 0.576159 |
| lag           |   25171 | 117 |      3.67373   | 0.331794 | 0.281106  |  0.063058   |        0.741541 |      0.140931   | 0.820426 |
| power         |   25170 | 191 |      1.38954   | 0.389279 | 0.337835  | -0.180823   |        1.26168  |     -0.586057   | 0.310256 |
| power         |   25171 | 117 |      3.81998   | 0.155256 | 0.121598  |  0.0471689  |        0.34699  |      0.10542    | 0.966873 |
| threshold     |   25170 | 191 |      1.95871   | 0.38618  | 0.335792  | -0.177589   |        1.25163  |     -0.575578   | 0.321775 |
| threshold     |   25171 | 117 |      5.67301   | 0.113053 | 0.0745183 | -0.00199873 |        0.252668 |     -0.00446706 | 0.976235 |
| lag_power     |   25170 | 191 |      0.0201956 | 0.225279 | 0.193064  | -0.0429125  |        0.730144 |     -0.139082   | 0.754277 |
| lag_power     |   25171 | 117 |      1.48455   | 0.393152 | 0.313951  |  0.028425   |        0.878674 |      0.0635283  | 0.668171 |
| lag_threshold |   25170 | 191 |      2.53236   | 0.28108  | 0.24165   | -0.0765512  |        0.910997 |     -0.248107   | 0.619702 |
| lag_threshold |   25171 | 117 |      7.57214   | 0.322354 | 0.236703  |  0.0124141  |        0.720446 |      0.0277448  | 0.787176 |

## Aroma fixed-parameter benchmark

| co2_mode      | partition_mode              |   aroma_data_wsse |   aroma_n_residuals |   co2_data_wsse |   combined_wsse | co2_params                                                             |
|:--------------|:----------------------------|------------------:|--------------------:|----------------:|----------------:|:-----------------------------------------------------------------------|
| threshold     | water_ethanol_GF_as_glucose |           6573.28 |                 366 |         244.736 |         6818.02 | {'threshold_fraction': 0.49809756401738725}                            |
| lag_threshold | water_ethanol_GF_as_glucose |           6661.31 |                 366 |         222.432 |         6883.74 | {'tau_h': 36.38236217778723, 'threshold_fraction': 0.8419628866863644} |
| lag_threshold | water_ethanol               |           6757.4  |                 366 |         222.432 |         6979.83 | {'tau_h': 36.38236217778723, 'threshold_fraction': 0.8419628866863644} |
| threshold     | water_ethanol               |           6814.8  |                 366 |         244.736 |         7059.54 | {'threshold_fraction': 0.49809756401738725}                            |
| instant       | water_ethanol_GF_as_glucose |           8431.48 |                 366 |         258.528 |         8690.01 | {}                                                                     |
| instant       | water_ethanol               |           9272.37 |                 366 |         258.528 |         9530.9  | {}                                                                     |
| lag           | water_ethanol_GF_as_glucose |           9666.29 |                 366 |         246.625 |         9912.92 | {'tau_h': 25.853734739141235}                                          |
| lag           | water_ethanol               |          10477    |                 366 |         246.625 |        10723.6  | {'tau_h': 25.853734739141235}                                          |
| lag_power     | water_ethanol_GF_as_glucose |         304907    |                 366 |         226.758 |       305134    | {'tau_h': 61.46629257205183, 'power_p': 2.499999999998967}             |
| lag_power     | water_ethanol               |         320711    |                 366 |         226.758 |       320938    | {'tau_h': 61.46629257205183, 'power_p': 2.499999999998967}             |

Selected CO2 mode by CO2 BIC: `lag_threshold`.

Best fixed aroma/partition combination in this run: `lag_threshold` with `water_ethanol_GF_as_glucose` if it also ranks best in the aroma table; otherwise inspect the table before adopting it.
