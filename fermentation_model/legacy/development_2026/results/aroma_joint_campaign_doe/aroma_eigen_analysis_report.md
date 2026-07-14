# Eigenvalue analysis for joint aroma campaign

This analysis decomposes the posterior FIM in scaled parameter coordinates.
Small eigenvalues correspond to weakly informed linear combinations of parameters.

## Design summary

| design                    |   n_directions |   min_relative_eigenvalue |   condition_number |   numeric_rank |   effective_rank |
|:--------------------------|---------------:|--------------------------:|-------------------:|---------------:|-----------------:|
| hybrid_current            |             21 |               1.29597e-05 |            77162   |             21 |          2.88591 |
| dopt_exchange_from_hybrid |             21 |               1.17413e-05 |            85169.5 |             21 |          2.88871 |

## Dominant weak-direction loadings

| design                    |   direction_rank_weak_to_strong |   eigenvalue |   relative_eigenvalue | parameter       | group        |    loading |   abs_loading |   loading_sq |
|:--------------------------|--------------------------------:|-------------:|----------------------:|:----------------|:-------------|-----------:|--------------:|-------------:|
| dopt_exchange_from_hybrid |                               1 |      15.2102 |           1.17413e-05 | sN              | fermentation |  0.689278  |     0.689278  |   0.475104   |
| dopt_exchange_from_hybrid |                               1 |      15.2102 |           1.17413e-05 | k_EO_growth     | synthesis    | -0.68165   |     0.68165   |   0.464646   |
| dopt_exchange_from_hybrid |                               1 |      15.2102 |           1.17413e-05 | k_IAA_growth    | synthesis    | -0.177919  |     0.177919  |   0.0316553  |
| dopt_exchange_from_hybrid |                               1 |      15.2102 |           1.17413e-05 | k_EA_growth     | synthesis    | -0.14773   |     0.14773   |   0.0218242  |
| dopt_exchange_from_hybrid |                               1 |      15.2102 |           1.17413e-05 | qXF             | fermentation |  0.0573703 |     0.0573703 |   0.00329135 |
| dopt_exchange_from_hybrid |                               2 |      17.4293 |           1.34543e-05 | m0              | fermentation |  0.937301  |     0.937301  |   0.878533   |
| dopt_exchange_from_hybrid |                               2 |      17.4293 |           1.34543e-05 | k_EO_growth     | synthesis    |  0.235342  |     0.235342  |   0.055386   |
| dopt_exchange_from_hybrid |                               2 |      17.4293 |           1.34543e-05 | sN              | fermentation |  0.209776  |     0.209776  |   0.0440061  |
| dopt_exchange_from_hybrid |                               2 |      17.4293 |           1.34543e-05 | qXG             | fermentation |  0.093118  |     0.093118  |   0.00867097 |
| dopt_exchange_from_hybrid |                               2 |      17.4293 |           1.34543e-05 | qXF             | fermentation |  0.0930974 |     0.0930974 |   0.00866712 |
| dopt_exchange_from_hybrid |                               3 |      18.6837 |           1.44226e-05 | k_EO_growth     | synthesis    |  0.669524  |     0.669524  |   0.448263   |
| dopt_exchange_from_hybrid |                               3 |      18.6837 |           1.44226e-05 | sN              | fermentation |  0.64595   |     0.64595   |   0.417251   |
| dopt_exchange_from_hybrid |                               3 |      18.6837 |           1.44226e-05 | m0              | fermentation | -0.281601  |     0.281601  |   0.0792989  |
| dopt_exchange_from_hybrid |                               3 |      18.6837 |           1.44226e-05 | qXG             | fermentation | -0.186171  |     0.186171  |   0.0346598  |
| dopt_exchange_from_hybrid |                               3 |      18.6837 |           1.44226e-05 | qXF             | fermentation | -0.106772  |     0.106772  |   0.0114002  |
| dopt_exchange_from_hybrid |                               4 |      26.8751 |           2.07459e-05 | k_EA_growth     | synthesis    |  0.93433   |     0.93433   |   0.872972   |
| dopt_exchange_from_hybrid |                               4 |      26.8751 |           2.07459e-05 | k_IAA_growth    | synthesis    |  0.197588  |     0.197588  |   0.0390409  |
| dopt_exchange_from_hybrid |                               4 |      26.8751 |           2.07459e-05 | qXF             | fermentation | -0.174937  |     0.174937  |   0.0306029  |
| dopt_exchange_from_hybrid |                               4 |      26.8751 |           2.07459e-05 | qXG             | fermentation | -0.139009  |     0.139009  |   0.0193236  |
| dopt_exchange_from_hybrid |                               4 |      26.8751 |           2.07459e-05 | sN              | fermentation |  0.136411  |     0.136411  |   0.0186081  |
| dopt_exchange_from_hybrid |                               5 |      32.7968 |           2.53171e-05 | qXF             | fermentation |  0.752851  |     0.752851  |   0.566785   |
| dopt_exchange_from_hybrid |                               5 |      32.7968 |           2.53171e-05 | qXG             | fermentation | -0.640783  |     0.640783  |   0.410602   |
| dopt_exchange_from_hybrid |                               5 |      32.7968 |           2.53171e-05 | iG              | fermentation |  0.0797726 |     0.0797726 |   0.00636367 |
| dopt_exchange_from_hybrid |                               5 |      32.7968 |           2.53171e-05 | sN              | fermentation | -0.072427  |     0.072427  |   0.00524568 |
| dopt_exchange_from_hybrid |                               5 |      32.7968 |           2.53171e-05 | k_EA_growth     | synthesis    |  0.0694311 |     0.0694311 |   0.00482068 |
| hybrid_current            |                               1 |      15.7683 |           1.29597e-05 | k_EO_growth     | synthesis    |  0.964714  |     0.964714  |   0.930673   |
| hybrid_current            |                               1 |      15.7683 |           1.29597e-05 | sN              | fermentation | -0.145134  |     0.145134  |   0.0210638  |
| hybrid_current            |                               1 |      15.7683 |           1.29597e-05 | k_IAA_growth    | synthesis    |  0.11402   |     0.11402   |   0.0130006  |
| hybrid_current            |                               1 |      15.7683 |           1.29597e-05 | k_EA_growth     | synthesis    |  0.113331  |     0.113331  |   0.0128438  |
| hybrid_current            |                               1 |      15.7683 |           1.29597e-05 | qXF             | fermentation | -0.0958235 |     0.0958235 |   0.00918215 |
| hybrid_current            |                               2 |      17.7417 |           1.45817e-05 | m0              | fermentation |  0.968738  |     0.968738  |   0.938454   |
| hybrid_current            |                               2 |      17.7417 |           1.45817e-05 | qXG             | fermentation |  0.129624  |     0.129624  |   0.0168024  |
| hybrid_current            |                               2 |      17.7417 |           1.45817e-05 | k_EO_growth     | synthesis    |  0.120403  |     0.120403  |   0.014497   |
| hybrid_current            |                               2 |      17.7417 |           1.45817e-05 | qXF             | fermentation |  0.112499  |     0.112499  |   0.0126559  |
| hybrid_current            |                               2 |      17.7417 |           1.45817e-05 | sN              | fermentation |  0.111666  |     0.111666  |   0.0124692  |
| hybrid_current            |                               3 |      24.3747 |           2.00333e-05 | sN              | fermentation |  0.684501  |     0.684501  |   0.468542   |
| hybrid_current            |                               3 |      24.3747 |           2.00333e-05 | k_EA_growth     | synthesis    | -0.616138  |     0.616138  |   0.379626   |
| hybrid_current            |                               3 |      24.3747 |           2.00333e-05 | k_IAA_growth    | synthesis    | -0.267615  |     0.267615  |   0.0716176  |
| hybrid_current            |                               3 |      24.3747 |           2.00333e-05 | k_EO_growth     | synthesis    |  0.19828   |     0.19828   |   0.0393152  |
| hybrid_current            |                               3 |      24.3747 |           2.00333e-05 | qXG             | fermentation | -0.125848  |     0.125848  |   0.0158377  |
| hybrid_current            |                               4 |      26.0191 |           2.13848e-05 | k_EA_growth     | synthesis    |  0.716475  |     0.716475  |   0.513337   |
| hybrid_current            |                               4 |      26.0191 |           2.13848e-05 | sN              | fermentation |  0.625662  |     0.625662  |   0.391453   |
| hybrid_current            |                               4 |      26.0191 |           2.13848e-05 | qXG             | fermentation | -0.278166  |     0.278166  |   0.0773762  |
| hybrid_current            |                               4 |      26.0191 |           2.13848e-05 | qXF             | fermentation | -0.10101   |     0.10101   |   0.010203   |
| hybrid_current            |                               4 |      26.0191 |           2.13848e-05 | k_EA_stationary | synthesis    | -0.0457961 |     0.0457961 |   0.00209729 |
| hybrid_current            |                               5 |      31.1383 |           2.55922e-05 | qXF             | fermentation |  0.755369  |     0.755369  |   0.570583   |
| hybrid_current            |                               5 |      31.1383 |           2.55922e-05 | qXG             | fermentation | -0.621084  |     0.621084  |   0.385745   |
| hybrid_current            |                               5 |      31.1383 |           2.55922e-05 | sN              | fermentation | -0.182861  |     0.182861  |   0.033438   |
| hybrid_current            |                               5 |      31.1383 |           2.55922e-05 | iG              | fermentation |  0.0752627 |     0.0752627 |   0.00566448 |
| hybrid_current            |                               5 |      31.1383 |           2.55922e-05 | sF              | fermentation |  0.0402718 |     0.0402718 |   0.00162182 |

## Interpretation

- All tested posterior FIMs are numerically full rank under the `lambda_max * 1e-12` threshold.
- The hybrid design has the better weakest relative eigenvalue and lower condition number.
- The weakest directions are mixed fermentation-aroma combinations, not isolated single parameters.
- The recurring weak components are `sN`, `m0`, `qXG/qXF`, and growth-phase aroma yields, especially `k_EO_growth`.
- Stationary aroma yields are not among the weakest eigendirection loadings, so the remaining aroma weakness is mostly tied to growth/uptake coupling.
- Because partition parameters were not included, this analysis should not be read as evidence that partition corrections are identifiable.
