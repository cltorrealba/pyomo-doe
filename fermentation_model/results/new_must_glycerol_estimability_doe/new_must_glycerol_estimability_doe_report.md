# New must glycerol estimability and DOE report

## Scope

This run uses the new natural and synthetic must databases, excludes CO2 and aroma states for now, and extends the kinetic state vector with dead biomass and glycerol.

State vector: `X, Xd, N, G, F, E, Gly`.

New glycerol equations:

$$\frac{dGly}{dt} = \left(\gamma_{G,0}\,\phi_G + \gamma_{F,0}\,\phi_F\right)X$$

where `phi_G` and `phi_F` are the same temperature, substrate saturation, glucose/fructose interaction, and ethanol inhibition factors used by the ethanol-production terms.

## Calibration summary

| fit                 |   n_batches | mediums           | parameters                                                                 |   n_parameters | success   |   status | message                                    |   nfev |   initial_wsse |   final_wsse |   n_residuals |   dof |   wsse_per_residual |   wsse_per_dof |   l2_lambda | l2_parameters                              |
|:--------------------|------------:|:------------------|:---------------------------------------------------------------------------|---------------:|:----------|---------:|:-------------------------------------------|-------:|---------------:|-------------:|--------------:|------:|--------------------:|---------------:|------------:|:-------------------------------------------|
| natural_reduced11   |           9 | natural           | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0                     |             11 | True      |        3 | `xtol` termination condition is satisfied. |     12 |        54080.6 |      13409.5 |           608 |   597 |             22.0551 |        22.4615 |           0 | nan                                        |
| synthetic_reduced11 |          10 | synthetic         | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0                     |             11 | True      |        3 | `xtol` termination condition is satisfied. |     16 |        54936.9 |      10723.8 |           542 |   531 |             19.7856 |        20.1955 |           0 | nan                                        |
| mixed_reduced11     |          19 | natural,synthetic | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0                     |             11 | True      |        3 | `xtol` termination condition is satisfied. |     16 |       109017   |      27138.9 |          1150 |  1139 |             23.5991 |        23.827  |           0 | nan                                        |
| mixed_full17_l2     |          19 | natural,synthetic | mu0,sN,qN,qXG,qXF,betaG0,sG,betaF0,sF,qEG,qEF,iG,iE,Kd0,m0,gammaG0,gammaF0 |             17 | True      |        3 | `xtol` termination condition is satisfied. |      9 |        27167.6 |      26132.3 |          1160 |  1143 |             22.5278 |        22.8629 |           1 | sN,qXG,qXF,sG,sF,iE,Kd0,m0,gammaG0,gammaF0 |

Lowest normalized weighted objective in this run: `synthetic_reduced11` with WSSE/residual `19.8`.

## Parameter estimability

| analysis                   | parameter   |        theta |   std_log_approx |   approx_95_multiplier |     fim_diag | active_bound   | classification     |
|:---------------------------|:------------|-------------:|-----------------:|-----------------------:|-------------:|:---------------|:-------------------|
| mixed_current_full_l2      | gammaF0     |  0.00312368  |        0.571337  |            3.0643      |    89.0197   | False          | moderate           |
| mixed_current_full_l2      | qXG         |  0.0830335   |        0.578937  |            3.11029     |    32.5086   | False          | moderate           |
| mixed_current_full_l2      | sN          | 22.8223      |        0.822136  |            5.00975     |     9.30466  | False          | weak_or_confounded |
| mixed_current_full_l2      | m0          |  0.00451924  |        6.25346   |       210402           |     4.54168  | False          | weak_or_confounded |
| mixed_current_full_l2      | mu0         |  0.061443    |        0.0170082 |            1.0339      | 68030.9      | False          | well_estimated     |
| mixed_current_full_l2      | qN          |  0.0140776   |        0.0172734 |            1.03444     | 39646.4      | False          | well_estimated     |
| mixed_current_full_l2      | qEG         |  0.998681    |        0.0305334 |            1.06167     | 54755.4      | False          | well_estimated     |
| mixed_current_full_l2      | gammaG0     |  0.0584596   |        0.0390349 |            1.07951     | 21052.8      | False          | well_estimated     |
| mixed_current_full_l2      | Kd0         |  0.000772526 |        0.041284  |            1.08428     |   785.489    | False          | well_estimated     |
| mixed_current_full_l2      | betaG0      |  0.429577    |        0.0431798 |            1.08832     | 22136.8      | False          | well_estimated     |
| mixed_current_full_l2      | qEF         |  0.600556    |        0.0812645 |            1.17266     | 46782.1      | False          | well_estimated     |
| mixed_current_full_l2      | betaF0      |  0.286346    |        0.0876544 |            1.18744     | 16120        | False          | well_estimated     |
| mixed_current_full_l2      | sG          |  0.0382373   |        0.113229  |            1.24848     |   847.263    | False          | well_estimated     |
| mixed_current_full_l2      | iE          |  0.00536305  |        0.120901  |            1.2674      |  2722.6      | False          | well_estimated     |
| mixed_current_full_l2      | sF          |  0.0509012   |        0.192923  |            1.45955     |   420.486    | False          | well_estimated     |
| mixed_current_full_l2      | iG          |  0.00330299  |        0.225759  |            1.55658     |   697.046    | False          | well_estimated     |
| mixed_current_full_l2      | qXF         |  0.147088    |        0.338948  |            1.9432      |   146.066    | False          | well_estimated     |
| mixed_current_reducedprior | qXG         |  0.1125      |        0.467242  |            2.49876     |    41.6236   | False          | moderate           |
| mixed_current_reducedprior | qXF         |  0.1125      |        0.487296  |            2.59893     |    54.0783   | False          | moderate           |
| mixed_current_reducedprior | sN          | 18           |        0.778109  |            4.59557     |    13.3653   | False          | weak_or_confounded |
| mixed_current_reducedprior | m0          |  0.01        |        2.66614   |          185.98        |    20.4241   | False          | weak_or_confounded |
| mixed_current_reducedprior | gammaF0     |  0.000158468 |       14.7495    |            3.58914e+12 |     0.203272 | False          | weak_or_confounded |
| mixed_current_reducedprior | mu0         |  0.073746    |        0.0223142 |            1.04471     | 64605.1      | False          | well_estimated     |
| mixed_current_reducedprior | qN          |  0.0173737   |        0.022483  |            1.04505     | 41859.9      | False          | well_estimated     |
| mixed_current_reducedprior | qEG         |  0.941025    |        0.0302216 |            1.06102     | 53027.5      | False          | well_estimated     |
| mixed_current_reducedprior | gammaG0     |  0.0607084   |        0.0429509 |            1.08783     | 24258.7      | False          | well_estimated     |
| mixed_current_reducedprior | Kd0         |  0.000888893 |        0.0442301 |            1.09056     |   681.671    | False          | well_estimated     |
| mixed_current_reducedprior | betaG0      |  0.40601     |        0.0501722 |            1.10334     | 20332.3      | False          | well_estimated     |
| mixed_current_reducedprior | qEF         |  0.636993    |        0.0846682 |            1.18051     | 43477.7      | False          | well_estimated     |
| mixed_current_reducedprior | sG          |  0.03        |        0.10516   |            1.22889     |  1242.47     | False          | well_estimated     |
| mixed_current_reducedprior | iE          |  0.00589815  |        0.109221  |            1.23871     |  2926.39     | False          | well_estimated     |
| mixed_current_reducedprior | betaF0      |  0.288099    |        0.110199  |            1.24109     | 13965.9      | False          | well_estimated     |
| mixed_current_reducedprior | sF          |  0.03        |        0.157677  |            1.36213     |   875.789    | False          | well_estimated     |
| mixed_current_reducedprior | iG          |  0.00305606  |        0.267754  |            1.69012     |   697.437    | False          | well_estimated     |
| natural_current            | sN          | 18           |        1.29125   |           12.5642      |     7.81023  | False          | weak_or_confounded |
| natural_current            | qXF         |  0.1125      |        1.37711   |           14.8668      |    12.2974   | False          | weak_or_confounded |
| natural_current            | gammaF0     |  0.00505987  |        1.37974   |           14.9439      |    31.7759   | False          | weak_or_confounded |
| natural_current            | qXG         |  0.1125      |        1.51765   |           19.5815      |    11.6115   | False          | weak_or_confounded |
| natural_current            | m0          |  0.01        |        5.1911    |        26227.3         |     9.33992  | False          | weak_or_confounded |
| natural_current            | qN          |  0.0289682   |        0.0457167 |            1.09374     | 18127.4      | False          | well_estimated     |
| natural_current            | mu0         |  0.10091     |        0.0458004 |            1.09392     | 25536        | False          | well_estimated     |
| natural_current            | gammaG0     |  0.0751174   |        0.0845291 |            1.18019     |  8089.05     | False          | well_estimated     |
| natural_current            | Kd0         |  0.0037451   |        0.0845458 |            1.18023     |   449.474    | False          | well_estimated     |
| natural_current            | qEG         |  1.17707     |        0.0888974 |            1.19034     | 32123.3      | False          | well_estimated     |
| natural_current            | betaG0      |  0.502511    |        0.128106  |            1.28542     | 22812.6      | False          | well_estimated     |
| natural_current            | qEF         |  1.04431     |        0.135985  |            1.30543     | 22249.3      | False          | well_estimated     |
| natural_current            | betaF0      |  0.487841    |        0.180151  |            1.42347     | 20813.3      | False          | well_estimated     |
| natural_current            | sF          |  0.03        |        0.195403  |            1.46666     |  1060.99     | False          | well_estimated     |
| natural_current            | sG          |  0.03        |        0.229029  |            1.56658     |   997.28     | False          | well_estimated     |
| natural_current            | iE          |  0.00789219  |        0.250945  |            1.63534     |  1876.34     | False          | well_estimated     |
| natural_current            | iG          |  0.0068611   |        0.336524  |            1.93399     |   594.362    | False          | well_estimated     |
| synthetic_current          | iG          |  0.00309291  |        0.504967  |            2.69052     |   474.12     | False          | moderate           |
| synthetic_current          | qXF         |  0.1125      |        0.729108  |            4.17474     |    91.0581   | False          | moderate           |
| synthetic_current          | mu0         |  0.050014    |        0.0610516 |            1.12711     | 11212.9      | True           | weak_or_confounded |
| synthetic_current          | qXG         |  0.1125      |        0.7802    |            4.61445     |    65.4131   | False          | weak_or_confounded |
| synthetic_current          | sN          | 18           |        1.89962   |           41.3992      |     5.62999  | False          | weak_or_confounded |
| synthetic_current          | m0          |  0.01        |        3.9228    |         2183.5         |    17.904    | False          | weak_or_confounded |
| synthetic_current          | gammaF0     |  0.000222185 |       21.2972    |            1.05765e+17 |     0.181926 | False          | weak_or_confounded |
| synthetic_current          | qN          |  0.00863748  |        0.032014  |            1.06476     | 21341.3      | False          | well_estimated     |
| synthetic_current          | Kd0         |  0.000908616 |        0.0386537 |            1.07871     |  1184.18     | False          | well_estimated     |
| synthetic_current          | gammaG0     |  0.0711298   |        0.0595076 |            1.12371     | 14799.4      | False          | well_estimated     |
| synthetic_current          | qEG         |  1.12466     |        0.0684748 |            1.14363     | 32721.9      | False          | well_estimated     |
| synthetic_current          | betaG0      |  0.554456    |        0.08321   |            1.17714     | 16735.8      | False          | well_estimated     |
| synthetic_current          | iE          |  0.0141443   |        0.106007  |            1.23093     |  5471.07     | False          | well_estimated     |
| synthetic_current          | qEF         |  0.735091    |        0.135259  |            1.30357     | 26904.2      | False          | well_estimated     |
| synthetic_current          | sG          |  0.03        |        0.195976  |            1.46831     |   923.08     | False          | well_estimated     |
| synthetic_current          | betaF0      |  0.304715    |        0.215029  |            1.52418     |  7165.68     | False          | well_estimated     |
| synthetic_current          | sF          |  0.03        |        0.238079  |            1.59462     |   499.387    | False          | well_estimated     |

## Profile likelihood check

When `profile_set` is present, `core` profiles refit the reduced kinetic nuisance set, while `weak` profiles use a more flexible nuisance set around the weak directions. The weak set is the conservative diagnostic for practical identifiability.

| parameter   |   n_success |   max_lr_stat | crosses_left_95   | crosses_right_95   | profile_identifiable_95   |   chi2_95_threshold | profile_set   |
|:------------|------------:|--------------:|:------------------|:-------------------|:--------------------------|--------------------:|:--------------|
| betaF0      |           6 |   1062.63     | True              | True               | True                      |             3.84146 | core          |
| betaG0      |           6 |    407.687    | True              | False              | False                     |             3.84146 | core          |
| iE          |           6 |    423.221    | True              | True               | True                      |             3.84146 | core          |
| iG          |           6 |    193.24     | True              | True               | True                      |             3.84146 | core          |
| mu0         |           6 |   3051.08     | True              | True               | True                      |             3.84146 | core          |
| qEF         |           6 |   5272.56     | True              | True               | True                      |             3.84146 | core          |
| qEG         |           6 |   6923.71     | True              | True               | True                      |             3.84146 | core          |
| qN          |           6 |   5911.77     | True              | True               | True                      |             3.84146 | core          |
| Kd0         |           6 |    374.761    | True              | True               | True                      |             3.84146 | weak          |
| gammaF0     |           6 |     18.0824   | True              | True               | True                      |             3.84146 | weak          |
| iG          |           6 |    129.045    | False             | True               | False                     |             3.84146 | weak          |
| m0          |           6 |      0.363627 | False             | False              | False                     |             3.84146 | weak          |
| qXF         |           6 |     17.8574   | True              | False              | False                     |             3.84146 | weak          |
| qXG         |           6 |      0        | False             | False              | False                     |             3.84146 | weak          |
| sF          |           6 |    531.035    | True              | False              | False                     |             3.84146 | weak          |
| sN          |           6 |      5.95037  | True              | False              | False                     |             3.84146 | weak          |

## DOE candidate ranking

| candidate                                | family              | medium    |   combined_logdet |   combined_min_relative_eigenvalue |   combined_condition_number |   weak_mean_var_reduction |   weak_worst_var_reduction |
|:-----------------------------------------|:--------------------|:----------|------------------:|-----------------------------------:|----------------------------:|--------------------------:|---------------------------:|
| synthetic_fructose_rich_glucose_pulse    | sugar_separation    | synthetic |           96.5933 |                        2.42064e-07 |                 4.13115e+06 |                 0.320321  |                 0.0808521  |
| synthetic_glucose_rich_fructose_pulse    | sugar_separation    | synthetic |           96.2053 |                        4.47907e-07 |                 2.23261e+06 |                 0.312653  |                 0.0662948  |
| synthetic_high_biomass_low_N_maintenance | maintenance         | synthetic |           94.5007 |                        3.15029e-07 |                 3.17431e+06 |                 0.245885  |                 0.0241819  |
| synthetic_late_ethanol_death_probe       | death               | synthetic |           93.5786 |                        2.02631e-07 |                 4.93508e+06 |                 0.173266  |                 0.0561066  |
| synthetic_viable_biomass_step            | biomass_input       | synthetic |           93.1374 |                        2.07452e-07 |                 4.82039e+06 |                 0.124881  |                 0.0267492  |
| synthetic_ethanol_inhibition_challenge   | ethanol_inhibition  | synthetic |           93.0878 |                        2.02375e-07 |                 4.94132e+06 |                 0.125932  |                 0.00635496 |
| synthetic_high_sugar_reference           | synthetic_baseline  | synthetic |           93.0832 |                        1.80787e-07 |                 5.53137e+06 |                 0.0877831 |                 0.0144999  |
| natural_glucose_pulse_after_growth       | natural_sugar_pulse | natural   |           92.9861 |                        2.42419e-07 |                 4.12508e+06 |                 0.128471  |                 0.0227744  |
| synthetic_low_yan_ladder                 | nitrogen_saturation | synthetic |           92.9012 |                        2.01292e-07 |                 4.9679e+06  |                 0.126468  |                 0.0260417  |
| natural_cold_hot_switch                  | natural_temperature | natural   |           92.5941 |                        1.95199e-07 |                 5.12297e+06 |                 0.0787423 |                 0.0258684  |
| synthetic_natural_like_control           | synthetic_baseline  | synthetic |           92.5355 |                        2.03071e-07 |                 4.92439e+06 |                 0.0944213 |                 0.0186531  |
| natural_control_18C                      | natural_baseline    | natural   |           92.519  |                        2.0288e-07  |                 4.92902e+06 |                 0.0938984 |                 0.0204371  |

## Selected campaign: hybrid D/E/A score

|   campaign_order | candidate                                | family              | medium    |   campaign_logdet |   campaign_min_relative_eigenvalue |   weak_mean_var_reduction |   weak_worst_var_reduction | rationale                                                                                    |
|-----------------:|:-----------------------------------------|:--------------------|:----------|------------------:|-----------------------------------:|--------------------------:|---------------------------:|:---------------------------------------------------------------------------------------------|
|                1 | synthetic_glucose_rich_fructose_pulse    | sugar_separation    | synthetic |           96.2053 |                        4.47907e-07 |                  0.312653 |                  0.0662948 | High glucose plus fructose pulse separates glucose and fructose uptake/yield directions.     |
|                2 | synthetic_fructose_rich_glucose_pulse    | sugar_separation    | synthetic |           99.9721 |                        4.22318e-07 |                  0.461598 |                  0.143794  | Fructose-rich run probes iG and fructose kinetic directions.                                 |
|                3 | synthetic_high_biomass_low_N_maintenance | maintenance         | synthetic |          102.173  |                        5.25994e-07 |                  0.549665 |                  0.160757  | High biomass with low N creates low-growth sugar consumption windows for m0.                 |
|                4 | synthetic_late_ethanol_death_probe       | death               | synthetic |          103.701  |                        5.25428e-07 |                  0.612382 |                  0.400003  | Late ethanol stress plus Xd observation targets Kd0/iE.                                      |
|                5 | synthetic_viable_biomass_step            | biomass_input       | synthetic |          104.688  |                        5.42482e-07 |                  0.633988 |                  0.416001  | Known viable biomass addition tests rate proportionality to X and supports yield separation. |
|                6 | synthetic_ethanol_inhibition_challenge   | ethanol_inhibition  | synthetic |          105.626  |                        5.39145e-07 |                  0.655239 |                  0.418289  | Initial ethanol decouples ethanol inhibition from ethanol generated by fermentation.         |
|                7 | synthetic_high_sugar_reference           | synthetic_baseline  | synthetic |          106.579  |                        5.16005e-07 |                  0.669314 |                  0.424793  | Synthetic high-sugar operating point close to the current synthetic dataset.                 |
|                8 | synthetic_low_yan_ladder                 | nitrogen_saturation | synthetic |          107.363  |                        5.15104e-07 |                  0.693156 |                  0.568027  | Low-YAN ladder targets sN and qN separation.                                                 |
|                9 | natural_glucose_pulse_after_growth       | natural_sugar_pulse | natural   |          107.999  |                        5.12034e-07 |                  0.702338 |                  0.573238  | Natural must with glucose perturbation to test sugar-transfer validity.                      |

## Pure D-optimal benchmark

|   campaign_order | candidate                                | family              | medium    |   campaign_logdet |   campaign_min_relative_eigenvalue |   weak_mean_var_reduction |   weak_worst_var_reduction |
|-----------------:|:-----------------------------------------|:--------------------|:----------|------------------:|-----------------------------------:|--------------------------:|---------------------------:|
|                1 | synthetic_fructose_rich_glucose_pulse    | sugar_separation    | synthetic |           96.5933 |                        2.42064e-07 |                  0.320321 |                  0.0808521 |
|                2 | synthetic_glucose_rich_fructose_pulse    | sugar_separation    | synthetic |           99.9721 |                        4.22318e-07 |                  0.461598 |                  0.143794  |
|                3 | synthetic_high_biomass_low_N_maintenance | maintenance         | synthetic |          102.173  |                        5.25994e-07 |                  0.549665 |                  0.160757  |
|                4 | synthetic_late_ethanol_death_probe       | death               | synthetic |          103.701  |                        5.25428e-07 |                  0.612382 |                  0.400003  |
|                5 | synthetic_high_sugar_reference           | synthetic_baseline  | synthetic |          104.747  |                        4.99391e-07 |                  0.62986  |                  0.407828  |
|                6 | synthetic_ethanol_inhibition_challenge   | ethanol_inhibition  | synthetic |          105.699  |                        5.00291e-07 |                  0.652097 |                  0.410096  |
|                7 | synthetic_viable_biomass_step            | biomass_input       | synthetic |          106.579  |                        5.16005e-07 |                  0.669314 |                  0.424793  |
|                8 | synthetic_low_yan_ladder                 | nitrogen_saturation | synthetic |          107.363  |                        5.15104e-07 |                  0.693156 |                  0.568027  |
|                9 | natural_glucose_pulse_after_growth       | natural_sugar_pulse | natural   |          107.999  |                        5.12034e-07 |                  0.702338 |                  0.573238  |

## Pyomo.DoE check

| candidate                             | status   | parameters                                             |   pyomo_logdet |   pyomo_min_eigenvalue |   pyomo_max_eigenvalue |   pyomo_min_relative_eigenvalue |   pyomo_condition_number |   pyomo_trace |   pyomo_trace_inv |   pyomo_rank_1e-8 |
|:--------------------------------------|:---------|:-------------------------------------------------------|---------------:|-----------------------:|-----------------------:|--------------------------------:|-------------------------:|--------------:|------------------:|------------------:|
| synthetic_glucose_rich_fructose_pulse | ok       | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |        52.3023 |              0.0189604 |                14457.3 |                     1.31147e-06 |                   762501 |       30132.7 |           54.4592 |                11 |

## Operational interpretation

- Sampling policy used for future designs: Monday-Friday only, 10:00-16:00, four samples per day during the first three process weekdays and two samples per process weekday afterwards.
- Pulse/adition times were snapped to feasible working-window times; no night or weekend manual action is required in the proposed candidate set.
- Temperature setpoints are assumed automatically programmable.
- The selected mixed campaign should be interpreted as model-based screening, not as final protocol approval; the first block of three fermentations should be used to update the prior before committing to the remaining blocks.

## Main conclusions to verify experimentally

- Natural and synthetic musts excite different regions: natural contributes matrix-transfer relevance and lower-sugar dynamics; synthetic contributes stronger sugar and pulse perturbations.
- A mixed approach is preferred unless the medium-transfer residuals show that one medium has a systematic model mismatch that the current structure cannot absorb.
- Glycerol observations directly support `gammaG0` and `gammaF0`, but these parameters remain coupled to the glucose/fructose fermentation factors; sugar-composition perturbations are therefore required.
- `Kd0`, `iE`, and `m0` need targeted stress or low-growth windows; ordinary control fermentations are insufficient.