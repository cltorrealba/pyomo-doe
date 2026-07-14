# Batch prioritization

The selected nine-fermentation campaign is executed as three batches of three fermentations.
The first batch is chosen for adaptive learning, not only for maximum global logdet.

## Recommended batches

|   batch | batch_experiments                                                                            |   n_experiments |   logdet |   mean_var_reduction |   worst_var_reduction |   sN_var_reduction |   condition_number |
|--------:|:---------------------------------------------------------------------------------------------|----------------:|---------:|---------------------:|----------------------:|-------------------:|-------------------:|
|       1 | fructose_rich_iG_probe | low_temp_growth_separation_plus_co2 | yan_saturation_scan           |               3 |  129.709 |             0.927167 |              0.65019  |           0.65019  |            86443.8 |
|       2 | cold_synthesis_hot_stripping | combined_stress_long_horizon | glucose_rich_growth_yield      |               6 |  149.394 |             0.965031 |              0.730844 |           0.730844 |           100750   |
|       3 | high_biomass_low_N_maintenance | ethanol_initial_challenge | glucose_pulse_after_N_depletion |               9 |  155.98  |             0.973679 |              0.768944 |           0.768944 |            77162   |

## Top first-batch combinations by hybrid score

| experiments                                                                                    |   n_experiments |   hybrid_score |   logdet |   min_relative_eigenvalue |   condition_number |   mean_var_reduction |   worst_var_reduction |   sN_var_reduction |   qN_var_reduction |   mu0_var_reduction |   k_EO_growth_var_reduction |
|:-----------------------------------------------------------------------------------------------|----------------:|---------------:|---------:|--------------------------:|-------------------:|---------------------:|----------------------:|-------------------:|-------------------:|--------------------:|----------------------------:|
| fructose_rich_iG_probe | low_temp_growth_separation_plus_co2 | combined_stress_long_horizon    |               3 |        134.399 |  138.55  |               7.77149e-06 |             128675 |             0.930722 |              0.503613 |           0.503613 |           0.811942 |            0.854709 |                    0.967694 |
| fructose_rich_iG_probe | combined_stress_long_horizon | glucose_rich_growth_yield              |               3 |        134.268 |  139.453 |               6.91329e-06 |             144649 |             0.931759 |              0.463103 |           0.463103 |           0.801308 |            0.853634 |                    0.956982 |
| fructose_rich_iG_probe | cold_synthesis_hot_stripping | low_temp_growth_separation_plus_co2    |               3 |        132.792 |  135.931 |               8.83091e-06 |             113239 |             0.931126 |              0.541207 |           0.541207 |           0.827575 |            0.869741 |                    0.969741 |
| fructose_rich_iG_probe | cold_synthesis_hot_stripping | glucose_rich_growth_yield              |               3 |        132.778 |  137.438 |               5.85667e-06 |             170745 |             0.933166 |              0.505224 |           0.505224 |           0.827459 |            0.869222 |                    0.964229 |
| fructose_rich_iG_probe | cold_synthesis_hot_stripping | combined_stress_long_horizon           |               3 |        132.754 |  138.174 |               3.82728e-06 |             261282 |             0.92499  |              0.513852 |           0.513852 |           0.802226 |            0.856157 |                    0.963489 |
| fructose_rich_iG_probe | glucose_rich_growth_yield | ethanol_initial_challenge                 |               3 |        131.659 |  136.386 |               6.79646e-06 |             147135 |             0.929839 |              0.48866  |           0.48866  |           0.802923 |            0.851037 |                    0.953449 |
| fructose_rich_iG_probe | cold_synthesis_hot_stripping | ethanol_initial_challenge              |               3 |        131.637 |  135.657 |               6.1992e-06  |             161311 |             0.924574 |              0.535817 |           0.535817 |           0.799769 |            0.842329 |                    0.960713 |
| fructose_rich_iG_probe | low_temp_growth_separation_plus_co2 | glucose_rich_growth_yield       |               3 |        131.564 |  135.449 |               7.90834e-06 |             126449 |             0.929423 |              0.515767 |           0.515767 |           0.829471 |            0.874755 |                    0.966834 |
| low_temp_growth_separation_plus_co2 | combined_stress_long_horizon | glucose_rich_growth_yield |               3 |        131.408 |  135.283 |               9.50867e-06 |             105167 |             0.929529 |              0.497819 |           0.497819 |           0.810904 |            0.85278  |                    0.964709 |
| fructose_rich_iG_probe | combined_stress_long_horizon | yan_saturation_scan                    |               3 |        130.962 |  134.414 |               5.15437e-06 |             194010 |             0.925281 |              0.582312 |           0.582312 |           0.797478 |            0.8462   |                    0.941174 |

## Rationale

The maximum-hybrid-score three-experiment batch is globally informative, but it leaves the weakest nitrogen saturation direction less protected.
For an adaptive first batch, the recommended set includes `yan_saturation_scan` because `sN/qN` is the most important practical bottleneck before committing to the remaining six fermentations.
