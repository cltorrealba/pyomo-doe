# Water-ethanol vs water-ethanol-sugar partition benchmark

Inherited aroma model: `ea_ethanol_nlimited`.

The sugar-aware mode uses original UNIFAC with water, ethanol, and total `G+F` treated as glucose-equivalent. Antoine vapor pressure is unchanged. Only the gas-liquid partition coefficient changes:

$$K_i=\frac{\gamma_i^{UNIFAC}(T,x_{water},x_{ethanol},x_{sugar})P_i^{sat,Antoine}(T)}{RTC_{tot,L}}.$$

## Partition coefficient effect over pilot trajectories

| species         |   n_process_points |      min |   median_K_ratio |   max |
|:----------------|-------------------:|---------:|-----------------:|------:|
| ethyl_acetate   |                128 | 0.786116 |         0.958402 |     1 |
| ethyl_octanoate |                128 | 0.457411 |         0.869233 |     1 |
| isoamyl_acetate |                128 | 0.607109 |         0.9148   |     1 |

## Fixed-parameter prediction comparison

| partition_mode              | start       | success   |   data_wsse |   n_data_residuals |   active_bound_count |
|:----------------------------|:------------|:----------|------------:|-------------------:|---------------------:|
| water_ethanol               | fixed_theta | True      |     9272.37 |                366 |                    2 |
| water_ethanol_GF_as_glucose | fixed_theta | True      |     8431.48 |                366 |                    2 |

## Refit aroma-parameter comparison

| partition_mode              | start           | success   |   status | message                                       |   nfev |   initial_objective_wsse |   final_objective_wsse |   data_wsse |   n_data_residuals |   active_bound_count |   fit_selection_score |
|:----------------------------|:----------------|:----------|---------:|:----------------------------------------------|-------:|-------------------------:|-----------------------:|------------:|-------------------:|---------------------:|----------------------:|
| water_ethanol               | fixed_reference | True      |        0 | selected global theta before aroma-only refit |      0 |                  9272.37 |                9272.37 |     9272.37 |                366 |                    2 |               9312.37 |
| water_ethanol_GF_as_glucose | fixed_reference | True      |        0 | selected global theta before aroma-only refit |      0 |                  8431.48 |                8431.48 |     8431.48 |                366 |                    2 |               8471.48 |
