# Historical natural estimability analysis

## Scope

This is an independent calibration and estimability analysis using only the historical `natural` medium subset. It is not a row-filtered view of the joint historical fit.

Actual temperature records replace the nominal temperature input. Valid CO2 sensor times enter the extended FIM as an online measurement schedule, while CO2 values remain diagnostic.

The core kinetic block is re-estimated. Secondary and aroma parameters are evaluated conditionally at the shared prior values, matching the logic of `fermentation_estimability_old_vs_lot1`.

## Core fit

| fit                     |   n_batches | mediums   | parameters                                             |   n_parameters | success   |   status | message                                    |   nfev |   initial_wsse |   final_wsse |   n_residuals |   dof |   wsse_per_residual |   wsse_per_dof |   l2_lambda |   l2_parameters | best_seed                  |
|:------------------------|------------:|:----------|:-------------------------------------------------------|---------------:|:----------|---------:|:-------------------------------------------|-------:|---------------:|-------------:|--------------:|------:|--------------------:|---------------:|------------:|----------------:|:---------------------------|
| historical_natural_only |           9 | natural   | mu0,qN,betaG0,betaF0,qEG,qEF,iG,iE,Kd0,gammaG0,gammaF0 |             11 | True      |        3 | `xtol` termination condition is satisfied. |      8 |        11672.1 |      11672.1 |           608 |   597 |             19.1975 |        19.5512 |           0 |             nan | profile_center_Kd0_release |

## Measurement support

| state                 |   n_observations |   n_batches_with_observations | used_in_core_fit   | used_in_extended_fim   |
|:----------------------|-----------------:|------------------------------:|:-------------------|:-----------------------|
| X                     |              100 |                             9 | True               | True                   |
| Xd                    |              100 |                             9 | True               | True                   |
| N                     |               91 |                             9 | True               | True                   |
| G                     |               91 |                             9 | True               | True                   |
| F                     |               91 |                             9 | True               | True                   |
| E                     |               44 |                             9 | True               | True                   |
| Gly                   |               91 |                             9 | True               | True                   |
| Pyr                   |               91 |                             9 | False              | True                   |
| AcAld                 |               79 |                             6 | False              | True                   |
| Acetate               |               50 |                             6 | False              | True                   |
| O2                    |               32 |                             3 | False              | True                   |
| ethyl_acetate_total   |                0 |                             0 | False              | False                  |
| isoamyl_acetate_total |                0 |                             0 | False              | False                  |
| ethyl_octanoate_total |                0 |                             0 | False              | False                  |

## FIM metrics

|   logdet |   min_eigenvalue |   max_eigenvalue |   min_relative_eigenvalue |   condition_number |       trace_inv |   rank_1e-8 | case                    | block                                     |   n_parameters |
|---------:|-----------------:|-----------------:|--------------------------:|-------------------:|----------------:|------------:|:------------------------|:------------------------------------------|---------------:|
|  54.7875 |      5.0272e-06  |          64485.6 |               7.79584e-11 |        1.28274e+10 | 64931.4         |          10 | historical_natural_only | core_observed                             |             11 |
|  90.2885 |      5.04339e-06 |          75247.7 |               6.70239e-11 |        1.492e+10   | 72582.2         |          16 | historical_natural_only | core_secondary_observed_plus_co2_schedule |             17 |
| -57.3339 |      0           |          75247.7 |               0           |        1e+12       |     1.66894e+06 |          16 | historical_natural_only | full_26_candidate_plus_co2_schedule       |             26 |

## Estimability classes

| classification     |   n_parameters |
|:-------------------|---------------:|
| moderate           |              2 |
| weak_or_confounded |             10 |
| well_estimated     |             14 |

### Strongest local directions by marginal parameter uncertainty

| parameter   |      theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification   |
|:------------|-----------:|-----------------:|-----------------------:|:---------------|:-----------------|
| mu0         | 0.0787279  |        0.0127642 |                1.02533 | False          | well_estimated   |
| qN          | 0.0210804  |        0.0146247 |                1.02908 | False          | well_estimated   |
| kAldS_N     | 2.35989    |        0.0299964 |                1.06056 | False          | well_estimated   |
| kAcStress   | 0.00261587 |        0.0343451 |                1.06963 | False          | well_estimated   |
| qEG         | 1.71281    |        0.065427  |                1.13682 | False          | well_estimated   |
| gammaG0     | 0.115559   |        0.0725893 |                1.15289 | False          | well_estimated   |
| kPyrS_N     | 0.850862   |        0.0763775 |                1.16149 | False          | well_estimated   |
| betaG0      | 0.428878   |        0.110838  |                1.24265 | False          | well_estimated   |
| iE          | 0.0401679  |        0.112305  |                1.24622 | False          | well_estimated   |
| kPyrDrain   | 0.00359103 |        0.145282  |                1.32943 | False          | well_estimated   |

### Weakest or confounded marginal directions

| parameter        |       theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:-----------------|------------:|-----------------:|-----------------------:|:---------------|:-------------------|
| alpha_EO_loss    | 1           |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_growth      | 0.1         |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| alpha_IAA_loss   | 1           |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| alpha_EA_loss    | 1           |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_EO_stationary  | 0.003       |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_EO_growth      | 0.0008      |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_IAA_stationary | 0.012       |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_IAA_growth     | 0.003       |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| k_EA_stationary  | 0.28        |          418.457 |            1.05765e+17 | False          | weak_or_confounded |
| gammaF0          | 0.000100002 |          304.937 |            1.05765e+17 | True           | weak_or_confounded |

## Profile likelihood screening

| case                    | parameter   |   n_success |   max_lr_stat | crosses_left_95   | crosses_right_95   | profile_identifiable_95   |   chi2_95_threshold |
|:------------------------|:------------|------------:|--------------:|:------------------|:-------------------|:--------------------------|--------------------:|
| historical_natural_only | Kd0         |           4 |       6.69301 | False             | True               | False                     |             3.84146 |
| historical_natural_only | qN          |           4 |    2195.54    | True              | True               | True                      |             3.84146 |

## Weak eigen-directions

| case                    |   weak_direction |   eigenvalue | dominant_parameters                                               | dominant_abs_loadings                           |
|:------------------------|-----------------:|-------------:|:------------------------------------------------------------------|:------------------------------------------------|
| historical_natural_only |                1 |            0 | k_EA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                2 |            0 | k_EA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG  | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                3 |            0 | k_IAA_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG     | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                4 |            0 | k_IAA_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                5 |            0 | k_EO_growth, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG      | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                6 |            0 | k_EO_stationary, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG  | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                7 |            0 | alpha_EA_loss, alpha_EO_loss, kPyrS_N, qN, betaG0, betaF0, qEG    | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |
| historical_natural_only |                8 |            0 | alpha_IAA_loss, alpha_EO_loss, qN, betaG0, betaF0, qEG, qEF       | 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 |

## Interpretation limits

- FIM classifications are local to the fitted/prior parameter point and use the same log-parameter finite-difference convention as the reference notebook.
- Aroma parameters without an aroma residual block remain unsupported by these notebooks even if aroma columns exist in the source workbook.
- Online CO2 is a gas-flow measurement. Its values are not treated as direct observations of the current cumulative production state until the gas-liquid observation model is calibrated.
- The secondary-state curves are conditional checks, not a new secondary-parameter calibration.
