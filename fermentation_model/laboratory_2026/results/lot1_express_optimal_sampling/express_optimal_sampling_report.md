# Lot 1 express optimal sampling

Greedy local FIM sampling refinement for Lot 1 after real inoculation at 2026-06-15 15:30.

Budgets per 216 h fermentation: 14 `full` samples and 6 `small` samples. The initial extended `full` sample at t=0, online CO2, and terminal condensate are treated as fixed. Greedy selection fills the remaining liquid sample budget on the 10:00/12:00/14:00/16:00 operational grid with at most 4 liquid samples per fermentation per day.

## Summary

| scenario             |   n_selected |   combined_logdet |   combined_trace_inv |   new_param_mean_var_reduction |   new_param_worst_var_reduction |   focus_core_mean_var_reduction |   focus_secondary_mean_var_reduction |   focus_aroma_mean_var_reduction |   mean_post_pulse_liquid_0_6h |   min_post_pulse_liquid_0_6h |   mean_post_pulse_full_0_6h |   min_post_pulse_full_0_6h |
|:---------------------|-------------:|------------------:|---------------------:|-------------------------------:|--------------------------------:|--------------------------------:|-------------------------------------:|---------------------------------:|------------------------------:|-----------------------------:|----------------------------:|---------------------------:|
| late_16_near_nominal |           60 |           176.593 |              1.54481 |                       0.958437 |                        0.829029 |                         0.58478 |                                    1 |                         0.9492   |                           0   |                            0 |                         0   |                          0 |
| midday_12_response   |           60 |           176.514 |              1.56341 |                       0.957806 |                        0.826474 |                         0.5813  |                                    1 |                         0.948429 |                           1.2 |                            0 |                         0.6 |                          0 |