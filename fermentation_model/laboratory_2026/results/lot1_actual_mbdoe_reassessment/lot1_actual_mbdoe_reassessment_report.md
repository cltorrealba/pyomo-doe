# Lot 1 actual MBDoE reassessment

This report re-evaluates the remaining campaign after Lot 1 using the same FIM workflow used for the operational DOE:

- residual vector from simulated outputs,
- log-parameter finite differences,
- FIM accumulation,
- D-opt criterion `logdet(FIM)`,
- hybrid score `logdet(FIM) + 2 log(min_relative_eigenvalue) - 0.05 log(trace_inv)`,
- variance reduction from `cov_prior` to `cov_post`.

The conservative case here uses the channels actually available in `fermentation_lot1_data_preview.executed.ipynb`: core states, secondary metabolites and CO2. Aroma information was not assumed for Lot 1.

## Lot 1: planned vs executed

| Basis | Hybrid score | logdet | min eig | trace inv | worst new-param var reduction | Kd0 var reduction | qN var reduction | qEG var reduction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Planned Lot 1 | 143.943 | 170.612 | 0.355 | 3.187 | 0.819 | 0.126 | 0.431 | 0.426 |
| Actual Lot 1 observed channels | 90.419 | 117.485 | 0.271 | 14.869 | ~0.000 | 0.157 | 0.300 | 0.304 |

Interpretation: the executed Lot 1 is informative for core/secondary kinetics, but it does not carry the full information assumed by the original campaign because the current processed dataset does not include aroma observations. This makes the remaining campaign still important.

## Candidate ranking after actual Lot 1

Top candidates by hybrid score, conditioned on `prior + FIM(actual Lot 1)`:

| Rank | Candidate | Hybrid score | logdet | min eig | trace inv | worst new-param var reduction |
|---:|---|---:|---:|---:|---:|---:|
| 1 | synthetic_glucose_rich_fructose_pulse | 130.335 | 157.133 | 0.335 | 3.856 | 0.117 |
| 2 | synthetic_lit_SM410_24C_highN_strip | 129.101 | 156.286 | 0.281 | 4.408 | 0.079 |
| 3 | synthetic_high_sugar_reference | 128.577 | 155.646 | 0.296 | 4.406 | 0.119 |
| 4 | synthetic_viable_biomass_step | 127.321 | 154.382 | 0.283 | 4.822 | 0.216 |
| 5 | natural_glucose_pulse_after_growth | 127.152 | 154.267 | 0.280 | 4.733 | 0.126 |
| 6 | synthetic_fast_CO2_aroma_strip_highN | 126.769 | 153.820 | 0.297 | 4.439 | 0.054 |

This supports keeping Lot 2: F4 and F5 are rank 1-2, and F6 is rank 4. The only candidate that outranks F6 is F7, but F7 is a Lot 3 validation/baseline run and does not replace the biomass perturbation role.

## Greedy selections after actual Lot 1

Hybrid and pure D-opt selected the same free six-experiment set:

1. synthetic_glucose_rich_fructose_pulse
2. synthetic_lit_SM410_24C_highN_strip
3. synthetic_high_sugar_reference
4. synthetic_viable_biomass_step
5. natural_glucose_pulse_after_growth
6. synthetic_low_yan_ladder

This differs from the original remaining design only in the last slot: `synthetic_low_yan_ladder` replaces `synthetic_fast_CO2_aroma_strip_highN`. However, the numerical gap is tiny and the original design remains slightly better in trace inverse, minimum eigenvalue, mean/worst variance reduction and Kd0 reduction.

## Scenario benchmark

| Scenario | Hybrid | logdet | min eig | trace inv | mean new-param var reduction | worst new-param var reduction | Kd0 var reduction |
|---|---:|---:|---:|---:|---:|---:|---:|
| Free greedy 6 | 151.815 | 178.987 | 0.394 | 2.771 | 0.887 | 0.447 | 0.240 |
| Lot2 fixed + greedy Lot3 | 151.815 | 178.987 | 0.394 | 2.771 | 0.887 | 0.447 | 0.240 |
| Original F4-F9 | 151.751 | 178.949 | 0.405 | 2.711 | 0.889 | 0.450 | 0.245 |
| Replace F7 with natural aroma-matrix | 151.171 | 178.330 | 0.396 | 2.778 | 0.886 | 0.444 | 0.235 |
| Replace F7 with death probe | 151.083 | 178.255 | 0.395 | 2.783 | 0.884 | 0.440 | 0.299 |

Key numerical reading:

- The original F4-F9 design is within 0.064 hybrid-score units of the free greedy optimum, a negligible difference at this scale.
- Original F4-F9 has better minimum eigenvalue (0.405 vs 0.394), better trace inverse (2.711 vs 2.771), better worst new-parameter variance reduction (0.450 vs 0.447), and better Kd0 variance reduction (0.245 vs 0.240) than the free greedy set.
- The death-probe replacement improves Kd0 variance reduction from 0.245 to 0.299, but worsens global design metrics. This is a targeted trade-off, not a globally better DOE.

## Weak directions

After actual Lot 1, the weakest directions are aroma parameters first, because aroma data were not counted in Lot 1. Among core kinetic parameters, the weakest are:

- Kd0: 0.157 variance reduction after Lot 1.
- betaG0: 0.251.
- qN: 0.300.
- qEG: 0.304.

After executing Lot 2 as planned (F4-F6), these improve to:

- Kd0: 0.321.
- qN: 0.413.
- mu0: 0.484.
- qEG: 0.493.

After the original full remaining campaign F4-F9:

- Kd0: 0.364.
- qN: 0.505.
- mu0: 0.530.
- qEG: 0.601.
- alpha_EA_loss: 0.906.

## Decision

The numerical recommendation is:

1. Keep Lot 2 unchanged: F4, F5 and F6 are strongly supported numerically.
2. Do not replace F7 yet. Contrary to my first qualitative concern, F7 remains selected by the conditioned greedy run.
3. Keep F8 because natural-must transfer remains necessary and it ranks high.
4. Keep F9 if aroma/CO2 stripping remains a campaign objective. A pure local D-opt/hybrid run can replace it with `synthetic_low_yan_ladder`, but the improvement is negligible and the scientific loss for aroma/condensate balance is not justified.
5. Use the death probe only if estimating Kd0 becomes the explicit priority; it improves Kd0 but worsens the campaign globally.
