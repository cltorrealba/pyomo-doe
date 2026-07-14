# Overnight validation technical summary

## Status

- `overnight_stderr.log` is empty; the long run completed without reported errors.
- The Pyomo.DoE checks succeeded for all 9 selected experiments using the reduced 11-parameter labeled model.

## L2 scan

| fit                  |   final_wsse |   wsse_per_residual |   l2_lambda |   nfev |
|:---------------------|-------------:|--------------------:|------------:|-------:|
| full17_l2_lambda_0.1 |      26108.5 |             22.5073 |         0.1 |      8 |
| full17_l2_lambda_1   |      26094.4 |             22.4951 |         1   |      7 |
| full17_l2_lambda_10  |      26092.8 |             22.4938 |        10   |      7 |
| full17_l2_lambda_50  |      26101.8 |             22.5016 |        50   |      8 |

## Multistart ranking

| fit           |   final_wsse |   wsse_per_residual |   nfev |
|:--------------|-------------:|--------------------:|-------:|
| multistart_07 |      23914.1 |             20.6156 |     16 |
| multistart_06 |      24775.2 |             21.3579 |     10 |
| multistart_05 |      24864.5 |             21.4349 |     14 |
| multistart_04 |      25046.9 |             21.5922 |     18 |
| multistart_10 |      25317.3 |             21.8253 |     16 |
| multistart_13 |      25511.1 |             21.9923 |     11 |
| multistart_11 |      25710.1 |             22.1639 |     19 |
| multistart_14 |      25734.4 |             22.1848 |     14 |
| multistart_08 |      26077.2 |             22.4804 |     10 |
| multistart_01 |      26094.4 |             22.4951 |      7 |

## Base vs best multistart theta

| parameter   |   base_theta |   best_theta |   ratio_best_base |
|:------------|-------------:|-------------:|------------------:|
| mu0         |  0.061443    |  0.0601697   |          0.979277 |
| sN          | 22.8223      |  8.75145     |          0.38346  |
| qN          |  0.0140776   |  0.0146092   |          1.03776  |
| qXG         |  0.0830335   |  0.081491    |          0.981423 |
| qXF         |  0.147088    |  0.0709327   |          0.482245 |
| betaG0      |  0.429577    |  0.779854    |          1.8154   |
| sG          |  0.0382373   |  0.0974551   |          2.54869  |
| betaF0      |  0.286346    |  0.510897    |          1.78419  |
| sF          |  0.0509012   |  0.178295    |          3.50277  |
| qEG         |  0.998681    |  1.53794     |          1.53997  |
| qEF         |  0.600556    |  1.26271     |          2.10257  |
| iG          |  0.00330299  |  0.00736668  |          2.23031  |
| iE          |  0.00536305  |  0.0258618   |          4.82223  |
| Kd0         |  0.000772526 |  0.000669074 |          0.866087 |
| m0          |  0.00451924  |  0.0183386   |          4.05789  |
| gammaG0     |  0.0584596   |  0.0975071   |          1.66794  |
| gammaF0     |  0.00312368  |  0.00138924  |          0.444746 |

## Full-nuisance profile likelihood

| parameter   | profile_identifiable_95   | crosses_left_95   | crosses_right_95   |   max_lr_stat |
|:------------|:--------------------------|:------------------|:-------------------|--------------:|
| mu0         | False                     | False             | True               |    1878.61    |
| sN          | False                     | True              | False              |       5.95037 |
| qN          | True                      | True              | True               |    5840.91    |
| qXG         | False                     | False             | False              |       0       |
| qXF         | False                     | True              | False              |      17.8574  |
| betaG0      | False                     | False             | False              |       0       |
| sG          | False                     | True              | False              |     105.519   |
| betaF0      | False                     | True              | False              |     424.028   |
| sF          | False                     | True              | False              |     531.035   |
| qEG         | False                     | True              | False              |    3934.36    |
| qEF         | False                     | True              | False              |    4017.49    |
| iG          | False                     | False             | True               |     129.045   |
| iE          | False                     | True              | False              |     228.244   |
| Kd0         | True                      | True              | True               |     374.761   |
| m0          | False                     | False             | False              |       0       |
| gammaG0     | True                      | True              | True               |     246.7     |
| gammaF0     | True                      | True              | True               |      14.1547  |

## Sampling policy benchmark

| policy       |   final_logdet |   final_min_relative_eigenvalue |   weak_mean_var_reduction |   weak_worst_var_reduction |
|:-------------|---------------:|--------------------------------:|--------------------------:|---------------------------:|
| front_loaded |        107.997 |                     5.10734e-07 |                  0.702303 |                   0.573201 |
| balanced     |        113.422 |                     5.81858e-07 |                  0.783129 |                   0.625768 |
| two_per_day  |        106.458 |                     4.9779e-07  |                  0.665751 |                   0.455179 |

## Base-theta selected campaign

|   campaign_order | candidate                                | medium    |   weak_mean_var_reduction |   weak_worst_var_reduction |
|-----------------:|:-----------------------------------------|:----------|--------------------------:|---------------------------:|
|                1 | synthetic_glucose_rich_fructose_pulse    | synthetic |                  0.312653 |                  0.0662948 |
|                2 | synthetic_fructose_rich_glucose_pulse    | synthetic |                  0.461598 |                  0.143794  |
|                3 | synthetic_high_biomass_low_N_maintenance | synthetic |                  0.549665 |                  0.160757  |
|                4 | synthetic_late_ethanol_death_probe       | synthetic |                  0.612382 |                  0.400003  |
|                5 | synthetic_viable_biomass_step            | synthetic |                  0.633988 |                  0.416001  |
|                6 | synthetic_ethanol_inhibition_challenge   | synthetic |                  0.655239 |                  0.418289  |
|                7 | synthetic_high_sugar_reference           | synthetic |                  0.669314 |                  0.424793  |
|                8 | synthetic_low_yan_ladder                 | synthetic |                  0.693156 |                  0.568027  |
|                9 | natural_glucose_pulse_after_growth       | natural   |                  0.702338 |                  0.573238  |

## Best-multistart selected campaign

|   campaign_order | candidate                                | medium    |   weak_mean_var_reduction |   weak_worst_var_reduction |
|-----------------:|:-----------------------------------------|:----------|--------------------------:|---------------------------:|
|                1 | synthetic_fructose_rich_glucose_pulse    | synthetic |                  0.318565 |                  0.0989876 |
|                2 | synthetic_glucose_rich_fructose_pulse    | synthetic |                  0.471565 |                  0.172418  |
|                3 | synthetic_high_biomass_low_N_maintenance | synthetic |                  0.557785 |                  0.192613  |
|                4 | synthetic_late_ethanol_death_probe       | synthetic |                  0.623506 |                  0.448382  |
|                5 | synthetic_ethanol_inhibition_challenge   | synthetic |                  0.651033 |                  0.450952  |
|                6 | synthetic_high_sugar_reference           | synthetic |                  0.6677   |                  0.457584  |
|                7 | synthetic_viable_biomass_step            | synthetic |                  0.684868 |                  0.478745  |
|                8 | synthetic_low_yan_ladder                 | synthetic |                  0.70512  |                  0.568474  |
|                9 | natural_cold_hot_switch                  | natural   |                  0.714017 |                  0.595593  |

## Pyomo.DoE selected-candidate checks

| candidate                                | status   |   pyomo_logdet |   pyomo_min_relative_eigenvalue |   pyomo_condition_number |
|:-----------------------------------------|:---------|---------------:|--------------------------------:|-------------------------:|
| synthetic_glucose_rich_fructose_pulse    | ok       |        52.3023 |                     1.31147e-06 |                 762501   |
| synthetic_fructose_rich_glucose_pulse    | ok       |        53.6033 |                     2.23324e-05 |                  44778   |
| synthetic_high_biomass_low_N_maintenance | ok       |        48.8194 |                     4.6289e-05  |                  21603.4 |
| synthetic_late_ethanol_death_probe       | ok       |        53.3124 |                     2.19598e-05 |                  45537.7 |
| synthetic_viable_biomass_step            | ok       |        54.7325 |                     3.11616e-05 |                  32090.8 |
| synthetic_ethanol_inhibition_challenge   | ok       |        54.1014 |                     1.99907e-05 |                  50023.3 |
| synthetic_high_sugar_reference           | ok       |        58.571  |                     1.31738e-05 |                  75908.3 |
| synthetic_low_yan_ladder                 | ok       |        45.1942 |                     1.65391e-05 |                  60462.7 |
| natural_glucose_pulse_after_growth       | ok       |        50.0832 |                     4.29991e-06 |                 232563   |