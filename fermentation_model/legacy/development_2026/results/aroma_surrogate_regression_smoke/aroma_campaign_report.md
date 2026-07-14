# Aroma campaign DOE

This report evaluates model-based candidate experiments for aroma and fermentation parameters.

The gas-liquid loss model is fixed by UNIFAC-derived equilibrium coefficients. The target vector is:

`k_EA_growth`, `k_EA_stationary`, `k_IAA_growth`, `k_IAA_stationary`, `k_EO_growth`, `k_EO_stationary`

Liquid aroma observations are included at the liquid sampling interval. A single terminal condensate observation per aroma is included as a final accumulated gas-loss constraint.

## Selected Campaign

|   campaign_order | candidate                 | family        |   horizon_h | temperature_c   |   campaign_logdet |   campaign_min_relative_eigenvalue |   campaign_condition_number |   mean_var_reduction |   worst_var_reduction | rationale                                                                                                       |
|-----------------:|:--------------------------|:--------------|------------:|:----------------|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:----------------------------------------------------------------------------------------------------------------|
|                1 | glucose_rich_growth_yield | sugar_initial |         168 | 18, 22, 24, 20  |           27.8158 |                        0.000298088 |                     3354.72 |             0.956709 |              0.878368 | High glucose must isolates glucose growth/fermentation terms and tests fructose response after glucose history. |

## Candidate Ranking

| candidate                 | family        |   horizon_h |   combined_logdet |   combined_min_relative_eigenvalue |   combined_condition_number |   mean_var_reduction |   worst_var_reduction | status   | error   |
|:--------------------------|:--------------|------------:|------------------:|-----------------------------------:|----------------------------:|---------------------:|----------------------:|:---------|:--------|
| glucose_rich_growth_yield | sugar_initial |         168 |           27.8158 |                        0.000298088 |                     3354.72 |             0.956709 |              0.878368 | ok       |         |

## Parameter Variance Reductions

| parameter        | group     |   campaign_var_ratio |   campaign_var_reduction |
|:-----------------|:----------|---------------------:|-------------------------:|
| k_EA_growth      | synthesis |          0.0799245   |                 0.920076 |
| k_EA_stationary  | synthesis |          0.000398536 |                 0.999601 |
| k_IAA_growth     | synthesis |          0.057326    |                 0.942674 |
| k_IAA_stationary | synthesis |          0.000137197 |                 0.999863 |
| k_EO_growth      | synthesis |          0.121632    |                 0.878368 |
| k_EO_stationary  | synthesis |          0.000330876 |                 0.999669 |

## Interpretation

- The terminal condensate improves the mass closure but does not provide time-resolved gas-loss dynamics.
- Partition parameters are fixed by default to avoid confounding gas-liquid equilibrium with trap efficiency.
- Temperature switches are intentionally included to separate synthesis windows from stripping windows.
- Nitrogen ladders and sugar pulses are retained because the synthesis model is phase-dependent through nitrogen limitation and sugar uptake.