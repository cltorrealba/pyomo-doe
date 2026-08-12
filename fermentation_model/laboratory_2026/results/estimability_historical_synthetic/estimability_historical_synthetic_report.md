# Historical synthetic estimability analysis

## Scope

This is an independent calibration and estimability analysis using only the historical `synthetic` medium subset. It is not a row-filtered view of the joint historical fit.

No matching continuous process-sensor archive was supplied for these synthetic historical batches.

The core kinetic block is re-estimated. Secondary and aroma parameters are evaluated conditionally at the shared prior values, matching the logic of `fermentation_estimability_old_vs_lot1`.

## Core fit

| fit                       |   n_batches | mediums   | parameters                                             |   n_parameters | success   |   status | message                                    |   nfev |   initial_wsse |   final_wsse |   n_residuals |   dof |   wsse_per_residual |   wsse_per_dof |   l2_lambda |   l2_parameters | best_seed                  |
|:--------------------------|------------:|:----------|:-------------------------------------------------------|---------------:|:----------|---------:|:-------------------------------------------|-------:|---------------:|-------------:|--------------:|------:|--------------------:|---------------:|------------:|----------------:|:---------------------------|
| historical_synthetic_only |          10 | synthetic | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |             11 | True      |        3 | `xtol` termination condition is satisfied. |      7 |        9304.15 |      9302.47 |           542 |   531 |             17.1632 |        17.5188 |           0 |             nan | profile_center_Kd0_release |

## Measurement support

| state                 |   n_observations |   n_batches_with_observations | used_in_core_fit   | used_in_extended_fim   |
|:----------------------|-----------------:|------------------------------:|:-------------------|:-----------------------|
| X                     |               76 |                            10 | True               | True                   |
| Xd                    |               86 |                            10 | True               | True                   |
| N                     |               83 |                            10 | True               | True                   |
| G                     |               86 |                            10 | True               | True                   |
| F                     |               86 |                            10 | True               | True                   |
| E                     |               40 |                            10 | True               | True                   |
| Gly                   |               85 |                            10 | True               | True                   |
| Pyr                   |               86 |                            10 | False              | True                   |
| AcAld                 |                0 |                             0 | False              | True                   |
| Acetate               |                0 |                             0 | False              | True                   |
| O2                    |                0 |                             0 | False              | True                   |
| ethyl_acetate_total   |                0 |                             0 | False              | False                  |
| isoamyl_acetate_total |                0 |                             0 | False              | False                  |
| ethyl_octanoate_total |                0 |                             0 | False              | False                  |

## FIM metrics

|    logdet |   min_eigenvalue |   max_eigenvalue |   min_relative_eigenvalue |   condition_number |        trace_inv |   rank_1e-8 | case                      | block                   |   n_parameters |
|----------:|-----------------:|-----------------:|--------------------------:|-------------------:|-----------------:|------------:|:--------------------------|:------------------------|---------------:|
|   68.0108 |      4.27165e-05 |          93217.6 |               4.58245e-10 |        2.18224e+09 |  17340.9         |          10 | historical_synthetic_only | core_observed           |             11 |
|   39.9815 |      0           |          95667   |               0           |        1e+12       | 302668           |          13 | historical_synthetic_only | core_secondary_observed |             17 |
| -105.48   |      0           |          95667   |               0           |        1e+12       |      1.75732e+06 |          13 | historical_synthetic_only | full_26_candidate       |             26 |

## Estimability classes

| classification     |   n_parameters |
|:-------------------|---------------:|
| weak_or_confounded |             14 |
| well_estimated     |             12 |

### Strongest local directions by marginal parameter uncertainty

| parameter   |      theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:------------|-----------:|-----------------:|-----------------------:|:---------------|:-------------------|
| qN          | 0.00878089 |       0.00886928 |                1.01754 | False          | well_estimated     |
| mu0         | 0.05       |       0.014305   |                1.02843 | True           | weak_or_confounded |
| Kd0         | 0.0007204  |       0.0381908  |                1.07773 | False          | well_estimated     |
| kPyrS_N     | 0.850862   |       0.0405073  |                1.08263 | False          | well_estimated     |
| qEG         | 2.3569     |       0.0439204  |                1.0899  | False          | well_estimated     |
| gammaG0     | 0.150049   |       0.051346   |                1.10588 | False          | well_estimated     |
| iE          | 0.0679131  |       0.058059   |                1.12052 | False          | well_estimated     |
| betaG0      | 1.14899    |       0.0581092  |                1.12063 | False          | well_estimated     |
| kPyrDrain   | 0.00359103 |       0.0708443  |                1.14896 | False          | well_estimated     |
| qEF         | 2.86071    |       0.0930672  |                1.20011 | False          | well_estimated     |

### Weakest or confounded marginal directions

| parameter        |      theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:-----------------|-----------:|-----------------:|-----------------------:|:---------------|:-------------------|
| alpha_EO_loss    | 1          |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_growth      | 0.1        |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| kAcAld           | 0.00139231 |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| kAldS_N          | 2.35989    |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_stationary  | 0.28       |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_IAA_growth     | 0.003      |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_IAA_stationary | 0.012      |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_EO_growth      | 0.0008     |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| k_EO_stationary  | 0.003      |          380.488 |            1.05765e+17 | False          | weak_or_confounded |
| alpha_EA_loss    | 1          |          380.488 |            1.05765e+17 | False          | weak_or_confounded |

## Profile likelihood screening

| case                      | parameter   |   n_success |   max_lr_stat | crosses_left_95   | crosses_right_95   | profile_identifiable_95   |   chi2_95_threshold |
|:--------------------------|:------------|------------:|--------------:|:------------------|:-------------------|:--------------------------|--------------------:|
| historical_synthetic_only | Kd0         |           4 |       399.652 | True              | True               | True                      |             3.84146 |
| historical_synthetic_only | qN          |           4 |      7670.13  | True              | True               | True                      |             3.84146 |

## Weak eigen-directions

| case                      |   weak_direction |   eigenvalue | dominant_parameters                                               | dominant_abs_loadings                           |
|:--------------------------|-----------------:|-------------:|:------------------------------------------------------------------|:------------------------------------------------|
| historical_synthetic_only |                1 |            0 | kAldS_N, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG          | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                2 |            0 | kAcAld, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG           | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                3 |            0 | kAcStress, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG        | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                4 |            0 | k_EA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                5 |            0 | k_EA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG  | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                6 |            0 | k_IAA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG     | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                7 |            0 | k_IAA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_synthetic_only |                8 |            0 | k_EO_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |

## Interpretation limits

- FIM classifications are local to the fitted/prior parameter point and use the same log-parameter finite-difference convention as the reference notebook.
- Aroma parameters without an aroma residual block remain unsupported by these notebooks even if aroma columns exist in the source workbook.
- Online CO2 is a gas-flow measurement. Its values are not treated as direct observations of the current cumulative production state until the gas-liquid observation model is calibrated.
- The secondary-state curves are conditional checks, not a new secondary-parameter calibration.
