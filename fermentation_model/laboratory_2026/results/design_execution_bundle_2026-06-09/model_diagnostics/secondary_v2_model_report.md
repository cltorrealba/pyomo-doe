# Secondary v2 model evaluation

## Purpose

This run tests whether the poor current fits for pyruvate, acetaldehyde, acetate, and DO are mainly due to the secondary-state structure. The core fermentation simulation is kept fixed and only the secondary layer is changed.

## v2 structure

The v2 layer adds nitrogen-phase and oxygen gates to the apparent production and drainage terms:

$$h_N=\frac{N}{N+K_N},\qquad h_{stat}=1-h_N,$$

$$g_{O_2}=\frac{O_2}{O_2+K_{O_2}},\qquad h_{ana}=\frac{K_{O_2}}{O_2+K_{O_2}}.$$

Pyruvate production is split between N-rich and stationary phases plus an early oxygen-linked source. Acetaldehyde production is split similarly and includes pyruvate and oxygen-linked contributions. Acetate keeps the acetaldehyde oxidation, stress, and assimilation terms. DO keeps passive transfer plus cellular consumption.

## Fit summary

| model                        | success   |   status | message                                    |   nfev |   initial_wsse |   final_wsse |   n_residuals |   wsse_per_residual |
|:-----------------------------|:----------|---------:|:-------------------------------------------|-------:|---------------:|-------------:|--------------:|--------------------:|
| secondary_v1_current         | True      |      nan | loaded from current secondary_joint fit    |    nan |          nan   |      8706.96 |           588 |            14.8077  |
| secondary_v2_phase_o2free    | True      |        2 | `ftol` termination condition is satisfied. |     16 |        10286.2 |      5400.34 |           592 |             9.1222  |
| secondary_v2_phase_o2fixed   | True      |        2 | `ftol` termination condition is satisfied. |     12 |        10286.2 |      5756.23 |           590 |             9.75632 |
| secondary_v2_reduced_o2fixed | True      |        2 | `ftol` termination condition is satisfied. |      8 |        16950.4 |      5756.43 |           589 |             9.77322 |

## Per-state comparison

| model                        | state   |   n_obs |      wsse |   wsse_per_obs |      rmse |       mae |   obs_median |   obs_range |
|:-----------------------------|:--------|--------:|----------:|---------------:|----------:|----------:|-------------:|------------:|
| secondary_v1_current         | Pyr     |     313 | 3065.92   |        9.79528 | 18.7785   | 14.7476   |       35     |       98    |
| secondary_v1_current         | AcAld   |     193 | 4855.45   |       25.1578  | 40.126    | 32.099    |      177     |      288    |
| secondary_v1_current         | Acetate |      50 |  687.448  |       13.749   |  0.129779 |  0.105739 |        0.3   |        0.39 |
| secondary_v1_current         | O2      |      32 |   98.1294 |        3.06654 |  0.350231 |  0.297618 |        0.615 |        1.35 |
| secondary_v2_phase_o2free    | Pyr     |     313 | 1861.1    |        5.94601 | 14.6307   | 11.1364   |       35     |       98    |
| secondary_v2_phase_o2free    | AcAld   |     193 | 2644.43   |       13.7017  | 29.6127   | 23.5435   |      177     |      288    |
| secondary_v2_phase_o2free    | Acetate |      50 |  673.245  |       13.4649  |  0.128431 |  0.104563 |        0.3   |        0.39 |
| secondary_v2_phase_o2free    | O2      |      32 |  188.851  |        5.90159 |  0.485864 |  0.42388  |        0.615 |        1.35 |
| secondary_v2_phase_o2fixed   | Pyr     |     313 | 2292.76   |        7.3251  | 16.239    | 12.1255   |       35     |       98    |
| secondary_v2_phase_o2fixed   | AcAld   |     193 | 2685.57   |       13.9148  | 29.8421   | 23.5987   |      177     |      288    |
| secondary_v2_phase_o2fixed   | Acetate |      50 |  677.581  |       13.5516  |  0.128844 |  0.105026 |        0.3   |        0.39 |
| secondary_v2_phase_o2fixed   | O2      |      32 |   98.1294 |        3.06654 |  0.350231 |  0.297618 |        0.615 |        1.35 |
| secondary_v2_reduced_o2fixed | Pyr     |     313 | 2291.95   |        7.32254 | 16.2361   | 12.1285   |       35     |       98    |
| secondary_v2_reduced_o2fixed | AcAld   |     193 | 2686.89   |       13.9217  | 29.8495   | 23.6152   |      177     |      288    |
| secondary_v2_reduced_o2fixed | Acetate |      50 |  677.116  |       13.5423  |  0.1288   |  0.10496  |        0.3   |        0.39 |
| secondary_v2_reduced_o2fixed | O2      |      32 |   98.1294 |        3.06654 |  0.350231 |  0.297618 |        0.615 |        1.35 |

## v2 weak/confounded parameters

| parameter   |      theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:------------|-----------:|-----------------:|-----------------------:|:---------------|:-------------------|
| kAldRed     | 1.0183e-05 |          5.07166 |                  20753 | False          | weak_or_confounded |

## Full v2 estimability

| parameter   |      theta |   std_log_approx |   approx_95_multiplier | active_bound   | classification     |
|:------------|-----------:|-----------------:|-----------------------:|:---------------|:-------------------|
| kPyrS_N     | 0.850862   |        0.0379795 |                1.07728 | False          | well_estimated     |
| kPyrO2      | 0.273189   |        0.100643  |                1.21806 | False          | well_estimated     |
| kPyrDrain   | 0.00359103 |        0.0604019 |                1.12568 | False          | well_estimated     |
| kAldS_N     | 2.35989    |        0.0131366 |                1.02608 | False          | well_estimated     |
| kAldRed     | 1.0183e-05 |        5.07166   |            20753       | False          | weak_or_confounded |
| kAcAld      | 0.00139231 |        0.0493938 |                1.10165 | False          | well_estimated     |
| kAcStress   | 0.00261587 |        0.0316527 |                1.064   | False          | well_estimated     |

## Interpretation

- v2 should be kept only if it improves state-wise errors enough to justify the extra parameters.
- Parameters that remain weak after v2 should be fixed or strongly regularized before using the model in DOE/MPCC.
- If v2 improves pyruvate/acetaldehyde but leaves acetate assimilation weak, this supports fixing acetate assimilation rather than forcing DOE to identify it.