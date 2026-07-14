# Pure D-opt benchmark after executed Lot 1

Base information matrix:

`FIM_base = FIM_prior + FIM(executed Lot 1)`

The executed Lot 1 FIM uses actual initial states, actual pulse doses/times, actual liquid sampling times and online CO2 from `fermentation_lot1_data_preview.executed.ipynb`. Lot 1 aromas are not assumed in this conservative run.

Objective used here:

`maximize logdet(FIM_base + sum(FIM_candidate))`

No hybrid penalty, no trace inverse penalty, and no explicit eigenvalue robustness term were used for selecting the D-opt set.

## D-opt result

Pure D-opt selected:

1. `synthetic_glucose_rich_fructose_pulse`
2. `synthetic_lit_SM410_24C_highN_strip`
3. `synthetic_high_sugar_reference`
4. `synthetic_viable_biomass_step`
5. `natural_glucose_pulse_after_growth`
6. `synthetic_low_yan_ladder`

The D-opt exchange run starting from the current campaign converged to the same set, just with `synthetic_high_sugar_reference` and `synthetic_viable_biomass_step` swapped in order. Since FIM addition is commutative, this is the same design set.

## Comparison against current remaining campaign

| Scenario | logdet | hybrid score | min eigenvalue | trace_inv | worst new-param var-red | Kd0 var-red | qN var-red | qEG var-red |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Pure D-opt / exchange | 178.987 | 151.815 | 0.394 | 2.771 | 0.447 | 0.240 | 0.253 | 0.461 |
| Current F4-F9 | 178.949 | 151.751 | 0.405 | 2.711 | 0.450 | 0.245 | 0.292 | 0.427 |

Numerical reading:

- D-opt improves logdet by only `0.038` over the current campaign.
- The current campaign is better in min eigenvalue, trace inverse, worst new-parameter variance reduction, Kd0 variance reduction and qN variance reduction.
- D-opt mainly improves qEG variance reduction, from `0.427` to `0.461`.
- The swap proposed by pure D-opt is: replace `synthetic_fast_CO2_aroma_strip_highN` with `synthetic_low_yan_ladder`.

## Lot 2 fixed, D-opt for Lot 3 only

If F4-F6 are fixed as the next lot, pure D-opt selects the following Lot 3:

1. `synthetic_high_sugar_reference`
2. `natural_glucose_pulse_after_growth`
3. `synthetic_low_yan_ladder`

This is again the same replacement: `synthetic_low_yan_ladder` instead of `synthetic_fast_CO2_aroma_strip_highN`.

## Interpretation

The pure D-opt behavior is stable but not compelling enough to automatically change the campaign. It chooses the low-YAN ladder because it increases total determinant slightly, mainly through core/nitrogen kinetic directions. However, the current design has better robustness metrics and keeps the dedicated CO2/aroma stripping experiment.

My technical reading is:

- If the objective is strictly D-opt on the current model/FIM, use `synthetic_low_yan_ladder` as F9.
- If aroma/CO2 stripping and condensate balance remain important scientific endpoints, keep `synthetic_fast_CO2_aroma_strip_highN`; the D-opt sacrifice is numerically tiny.
- Lot 2 should not change under either criterion.

## Files

- `selected_after_actual_lot1_dopt_free6.csv`
- `selected_after_actual_lot1_dopt_exchange_from_current.csv`
- `selected_after_actual_lot1_dopt_exchange_from_greedy.csv`
- `selected_lot3_dopt_after_actual_lot1_and_fixed_lot2.csv`
- `scenario_metrics_after_actual_lot1.csv`
