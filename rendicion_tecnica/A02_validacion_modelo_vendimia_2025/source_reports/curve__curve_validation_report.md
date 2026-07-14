# Curve-level validation for the new-must glycerol DOE

## Purpose

This validation checks whether the calibration and DOE conclusions survive curve-level inspection. It compares the base theta against the best overnight multistart theta, simulates the selected DOE candidates, checks feasibility flags, and quantifies trajectory redundancy.

## Current-data fit comparison

| medium    | state   |   base_weighted_rmse |   best_weighted_rmse |   best_minus_base_weighted_rmse |
|:----------|:--------|---------------------:|---------------------:|--------------------------------:|
| natural   | E       |              6.30419 |              6.08947 |                      -0.214718  |
| natural   | F       |              5.82539 |              5.64507 |                      -0.180328  |
| natural   | G       |              5.18396 |              5.12519 |                      -0.058765  |
| natural   | Gly     |              3.52569 |              3.39701 |                      -0.128683  |
| natural   | N       |              4.28451 |              4.24807 |                      -0.0364374 |
| natural   | X       |              6.29427 |              6.01379 |                      -0.280483  |
| natural   | Xd      |              1.83397 |              1.79259 |                      -0.0413817 |
| synthetic | E       |              4.01818 |              3.26776 |                      -0.750423  |
| synthetic | F       |              3.28353 |              2.35161 |                      -0.931918  |
| synthetic | G       |              4.05645 |              3.63275 |                      -0.423695  |
| synthetic | Gly     |              2.31366 |              2.02307 |                      -0.290597  |
| synthetic | N       |              7.02814 |              7.01159 |                      -0.0165447 |
| synthetic | X       |              3.87264 |              3.68324 |                      -0.189399  |
| synthetic | Xd      |              5.62692 |              5.715   |                       0.088076  |

Negative `best_minus_base_weighted_rmse` means the best multistart improves that state/medium relative to the base theta.

## Predictive feasibility flags

No hard feasibility flags were triggered by the selected candidate simulations.

## Base-vs-best theta prediction disagreement

| candidate                                | medium    | family              |   theta_disagreement_weighted_rmse |
|:-----------------------------------------|:----------|:--------------------|-----------------------------------:|
| synthetic_fructose_rich_glucose_pulse    | synthetic | sugar_separation    |                           2.16834  |
| synthetic_late_ethanol_death_probe       | synthetic | death               |                           1.54572  |
| synthetic_low_yan_ladder                 | synthetic | nitrogen_saturation |                           1.21178  |
| synthetic_viable_biomass_step            | synthetic | biomass_input       |                           1.20851  |
| natural_glucose_pulse_after_growth       | natural   | natural_sugar_pulse |                           1.13708  |
| synthetic_glucose_rich_fructose_pulse    | synthetic | sugar_separation    |                           1.0792   |
| natural_cold_hot_switch                  | natural   | natural_temperature |                           1.0249   |
| synthetic_high_sugar_reference           | synthetic | synthetic_baseline  |                           0.872993 |
| synthetic_ethanol_inhibition_challenge   | synthetic | ethanol_inhibition  |                           0.792025 |
| synthetic_high_biomass_low_N_maintenance | synthetic | maintenance         |                           0.609699 |

High disagreement means the experiment is informative but also sensitive to the current local optimum.

## Most similar selected designs

| candidate_a                              | candidate_b                            |   scaled_rms_distance |   trajectory_correlation |
|:-----------------------------------------|:---------------------------------------|----------------------:|-------------------------:|
| synthetic_late_ethanol_death_probe       | synthetic_ethanol_inhibition_challenge |              0.203275 |                 0.965515 |
| natural_glucose_pulse_after_growth       | natural_cold_hot_switch                |              0.242676 |                 0.942738 |
| synthetic_high_sugar_reference           | natural_glucose_pulse_after_growth     |              0.348428 |                 0.91444  |
| synthetic_viable_biomass_step            | natural_cold_hot_switch                |              0.351647 |                 0.924946 |
| synthetic_viable_biomass_step            | natural_glucose_pulse_after_growth     |              0.364766 |                 0.893313 |
| synthetic_ethanol_inhibition_challenge   | natural_glucose_pulse_after_growth     |              0.372514 |                 0.84698  |
| synthetic_late_ethanol_death_probe       | natural_glucose_pulse_after_growth     |              0.381594 |                 0.839247 |
| synthetic_high_sugar_reference           | synthetic_low_yan_ladder               |              0.392754 |                 0.92874  |
| synthetic_glucose_rich_fructose_pulse    | natural_glucose_pulse_after_growth     |              0.401434 |                 0.887127 |
| synthetic_ethanol_inhibition_challenge   | synthetic_high_sugar_reference         |              0.409665 |                 0.854541 |
| synthetic_low_yan_ladder                 | natural_glucose_pulse_after_growth     |              0.43218  |                 0.787245 |
| synthetic_high_sugar_reference           | natural_cold_hot_switch                |              0.433121 |                 0.884875 |
| synthetic_low_yan_ladder                 | natural_cold_hot_switch                |              0.445248 |                 0.745038 |
| synthetic_viable_biomass_step            | synthetic_ethanol_inhibition_challenge |              0.449918 |                 0.82349  |
| synthetic_fructose_rich_glucose_pulse    | synthetic_high_sugar_reference         |              0.456606 |                 0.833328 |
| synthetic_glucose_rich_fructose_pulse    | synthetic_high_sugar_reference         |              0.45726  |                 0.833709 |
| synthetic_viable_biomass_step            | synthetic_high_sugar_reference         |              0.458637 |                 0.832108 |
| synthetic_high_biomass_low_N_maintenance | synthetic_late_ethanol_death_probe     |              0.467315 |                 0.829362 |
| synthetic_ethanol_inhibition_challenge   | natural_cold_hot_switch                |              0.472377 |                 0.755707 |
| synthetic_late_ethanol_death_probe       | natural_cold_hot_switch                |              0.476762 |                 0.725388 |

Low distance and high correlation indicate potential redundancy.

## Interpretation

- The design set should be kept only if the predictive curves are plausible and the first block contains at least one natural-must transfer test.
- If an experiment triggers feasibility flags under either theta, treat it as a candidate for manual review rather than automatic execution.
- If two candidates have very low pairwise distance, choose the one that better targets the weak profile-likelihood direction or is easier operationally.