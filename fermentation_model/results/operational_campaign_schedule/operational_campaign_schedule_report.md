# Operational campaign schedule assessment

Batch mode: `recommended`.
Reference start: `2026-06-15 08:00:00`.
Temperature mode: `automatic`.
Automated dosing channel: `E`.

Operational assumptions:

- Lab days: Monday-Friday.
- Manual sampling slots: 09:00, 11:00, 14:00, and 16:00.
- Manual action slots for non-automated pulses: 09:00, 11:00, 14:00, and 16:00.
- At most one dosing channel is automated. The recommended/default automated channel is `E`.
- Temperature setpoint changes are assumed automatically programmed unless `temperature_mode = manual`.
- Initial formulation/inoculation actions at `t = 0 h` are allowed as setup actions.
- Sampling is front-loaded with up to four samples per fermentation per weekday through approximately 96 h, then one morning sample per weekday plus a final feasible terminal sample.
- Online CO2 is not constrained by staffing; this schedule covers manual sampling and manual actions.

## Daily manual load

| date       |   manual_samples |   manual_pulses |   manual_temperature_changes |   total_manual_events |
|:-----------|-----------------:|----------------:|-----------------------------:|----------------------:|
| 2026-06-15 |               12 |               0 |                            0 |                    12 |
| 2026-06-16 |               12 |               2 |                            0 |                    14 |
| 2026-06-17 |               12 |               3 |                            0 |                    15 |
| 2026-06-18 |               12 |               2 |                            0 |                    14 |
| 2026-06-19 |                6 |               0 |                            0 |                     6 |
| 2026-06-22 |                3 |               0 |                            0 |                     3 |
| 2026-06-23 |                1 |               0 |                            0 |                     1 |
| 2026-06-29 |               12 |               0 |                            0 |                    12 |
| 2026-06-30 |               12 |               1 |                            0 |                    13 |
| 2026-07-01 |               12 |               2 |                            0 |                    14 |
| 2026-07-02 |               12 |               2 |                            0 |                    14 |
| 2026-07-03 |                6 |               2 |                            0 |                     8 |
| 2026-07-06 |                3 |               0 |                            0 |                     3 |
| 2026-07-07 |                2 |               0 |                            0 |                     2 |
| 2026-07-08 |                1 |               0 |                            0 |                     1 |
| 2026-07-09 |                1 |               0 |                            0 |                     1 |
| 2026-07-13 |               12 |               0 |                            0 |                    12 |
| 2026-07-14 |               12 |               0 |                            0 |                    12 |
| 2026-07-15 |               12 |               2 |                            0 |                    14 |
| 2026-07-16 |               12 |               2 |                            0 |                    14 |
| 2026-07-17 |                6 |               1 |                            0 |                     7 |
| 2026-07-20 |                3 |               0 |                            0 |                     3 |
| 2026-07-21 |                2 |               0 |                            0 |                     2 |

## Automated dosing channel comparison

| automated_channel   |   large_manual_pulse_shifts |   max_abs_manual_pulse_shift_h |   max_total_manual_events_per_day |   max_manual_samples_per_day |   max_manual_pulses_per_day |
|:--------------------|----------------------------:|-------------------------------:|----------------------------------:|-----------------------------:|----------------------------:|
| none                |                           1 |                             16 |                                15 |                           12 |                           3 |
| N                   |                           1 |                             16 |                                14 |                           12 |                           3 |
| G                   |                           1 |                             16 |                                15 |                           12 |                           3 |
| F                   |                           1 |                             16 |                                15 |                           12 |                           3 |
| E                   |                           0 |                              8 |                                15 |                           12 |                           3 |

## Batch assignment

|   batch | candidate                           | start_datetime      |   horizon_h | temperature_c   |
|--------:|:------------------------------------|:--------------------|------------:|:----------------|
|       1 | fructose_rich_iG_probe              | 2026-06-15 08:00:00 |         192 | 18, 20, 24, 22  |
|       1 | low_temp_growth_separation_plus_co2 | 2026-06-15 08:00:00 |         168 | 15, 17, 20, 22  |
|       1 | yan_saturation_scan                 | 2026-06-15 08:00:00 |         144 | 15, 18, 22, 22  |
|       2 | cold_synthesis_hot_stripping        | 2026-06-29 08:00:00 |         192 | 15, 16, 24, 25  |
|       2 | combined_stress_long_horizon        | 2026-06-29 08:00:00 |         240 | 15, 22, 25, 18  |
|       2 | glucose_rich_growth_yield           | 2026-06-29 08:00:00 |         168 | 18, 22, 24, 20  |
|       3 | high_biomass_low_N_maintenance      | 2026-07-13 08:00:00 |         168 | 16, 18, 20, 18  |
|       3 | ethanol_initial_challenge           | 2026-07-13 08:00:00 |         192 | 18, 22, 25, 22  |
|       3 | glucose_pulse_after_N_depletion     | 2026-07-13 08:00:00 |         192 | 17, 20, 24, 22  |

## Pulse adjustments requiring reconsideration

No pulse shift exceeded the large-shift threshold.

## Manual temperature changes requiring reconsideration

No manual temperature-change shift exceeded the large-shift threshold.

## Design-level sample counts

| candidate                           |   nominal_horizon_h |   operational_final_h |   n_manual_samples |   n_final_samples |
|:------------------------------------|--------------------:|----------------------:|-------------------:|------------------:|
| cold_synthesis_hot_stripping        |                 192 |                   193 |                 20 |                 1 |
| combined_stress_long_horizon        |                 240 |                   241 |                 22 |                 1 |
| ethanol_initial_challenge           |                 192 |                   193 |                 20 |                 1 |
| fructose_rich_iG_probe              |                 192 |                   193 |                 20 |                 1 |
| glucose_pulse_after_N_depletion     |                 192 |                   193 |                 20 |                 1 |
| glucose_rich_growth_yield           |                 168 |                   169 |                 19 |                 1 |
| high_biomass_low_N_maintenance      |                 168 |                   169 |                 19 |                 1 |
| low_temp_growth_separation_plus_co2 |                 168 |                   169 |                 19 |                 1 |
| yan_saturation_scan                 |                 144 |                   169 |                 19 |                 1 |

## Recommended DOE handling

The current Pyomo DoE workflow uses uniform relative sampling intervals. Before freezing the campaign, the next technical step is to pass these explicit operational sample/action times into the FIM model and recompute the posterior FIM.

The acceptance criterion should be: preserve the hybrid campaign unless the operational calendar reduces the weakest variance reduction or minimum relative eigenvalue materially. If the degradation is large, redesign start dates, pulse times, or candidate selection under the calendar constraints.
