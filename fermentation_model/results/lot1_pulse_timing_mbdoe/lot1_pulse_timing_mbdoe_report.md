# Lot 1 pulse-timing MBDoE benchmark

Real inoculation start assumed: `2026-06-15 15:30:00`.

Known actual-initial override used in this quick benchmark: `synthetic_lit_SM410_18C_highN_ester` N = 0.320 kg/m3 (320 mg/L). Other initials remain nominal until measured values are loaded.

Scenarios:

- `early_10_old_clock`: pulses at 10:00 on the target calendar day, equivalent to keeping the old 08:00-clock schedule after a 15:30 inoculation.
- `midday_12_response`: pulses at 12:00 on the target calendar day to allow +2 h and +4 h response samples.
- `late_16_near_nominal`: pulses at 16:00 on the target calendar day, closest to the original process-time target without leaving the work window.

## Summary

| scenario             |   lot1_new_logdet |   combined_logdet |   combined_trace_inv |   new_param_mean_var_reduction |   new_param_worst_var_reduction |   focus_core_mean_var_reduction |   focus_secondary_mean_var_reduction |   focus_aroma_mean_var_reduction |   mean_post_pulse_liquid_0_6h |   min_post_pulse_liquid_0_6h |   mean_post_pulse_full_0_6h |   min_post_pulse_full_0_6h |
|:---------------------|------------------:|------------------:|---------------------:|-------------------------------:|--------------------------------:|--------------------------------:|-------------------------------------:|---------------------------------:|------------------------------:|-----------------------------:|----------------------------:|---------------------------:|
| late_16_near_nominal |           155.032 |           169.293 |              2.99673 |                       0.956793 |                        0.821486 |                        0.523753 |                                    1 |                         0.947192 |                             0 |                            0 |                         0   |                          0 |
| midday_12_response   |           154.731 |           169.029 |              3.10147 |                       0.95637  |                        0.820117 |                        0.513363 |                                    1 |                         0.946675 |                             2 |                            2 |                         1.2 |                          1 |
| early_10_old_clock   |           154.693 |           168.991 |              3.14482 |                       0.956151 |                        0.81937  |                        0.509293 |                                    1 |                         0.946407 |                             3 |                            3 |                         1.8 |                          1 |

## Interpretation

Higher `combined_logdet` and lower `combined_trace_inv` are better. Post-pulse counts help detect designs that cannot observe the immediate response to a pulse.