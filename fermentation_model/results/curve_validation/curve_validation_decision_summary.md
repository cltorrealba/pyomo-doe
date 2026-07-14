# Curve-validation decision summary

## Main result

The selected DOE candidates are dynamically plausible under both the base theta and the best overnight multistart theta. No candidate triggered hard feasibility flags for ethanol, biomass, dead biomass, glycerol, or final residual sugar.

The best multistart improves the current-data curves in almost every medium/state combination. The only small deterioration is synthetic `Xd`, which remains one of the noisiest/most model-sensitive observables.

## Current-data curve fit

The largest curve-level improvements from the best multistart occur in synthetic must:

- `F`: weighted RMSE improves by about `0.93`.
- `E`: weighted RMSE improves by about `0.75`.
- `G`: weighted RMSE improves by about `0.42`.
- `Gly`: weighted RMSE improves by about `0.29`.

Natural must also improves, but more modestly. This supports using the best multistart as a serious alternative prior, but not as an automatic final parameter set because it moved weak/confounded parameters substantially.

## Predictive plausibility

All selected design predictions remain in plausible ranges:

- Ethanol stays around `67-103 g/L` across the selected simulations.
- Viable biomass stays below `2.5 kg/m3`.
- Dead biomass stays below about `1.0 kg/m3`.
- Glycerol stays below about `10 g/L`.
- Final residual sugar is acceptable for the designed perturbation experiments.

The fructose-rich/glucose-pulse and glucose-rich/fructose-pulse designs show clear separation and are worth keeping.

## Redundancy

Two pairs are notably similar:

- `synthetic_late_ethanol_death_probe` vs `synthetic_ethanol_inhibition_challenge`
  - scaled RMS distance: about `0.20`
  - trajectory correlation: about `0.97`
  - Interpretation: do not put both in the first block. Keep both only if running the full 9-fermentation campaign.

- `natural_glucose_pulse_after_growth` vs `natural_cold_hot_switch`
  - scaled RMS distance: about `0.24`
  - trajectory correlation: about `0.94`
  - Interpretation: choose one natural transfer test in the first block. The balanced-sampling and best-multistart analysis favor `natural_cold_hot_switch`.

## Recommended first block of 3

Use a mixed first block to avoid committing all early information to synthetic must:

1. `synthetic_glucose_rich_fructose_pulse`
2. `synthetic_fructose_rich_glucose_pulse`
3. `natural_cold_hot_switch`

Rationale:

- The first two experiments are robustly selected across theta assumptions and directly target G/F separation plus glycerol formation.
- The natural experiment tests transfer to real must under temperature excitation.
- This block avoids redundant ethanol-stress experiments until the prior is updated.

## Full campaign if 9 fermentations are feasible

Keep this as the current full-campaign candidate set:

1. `synthetic_glucose_rich_fructose_pulse`
2. `synthetic_fructose_rich_glucose_pulse`
3. `synthetic_high_biomass_low_N_maintenance`
4. `synthetic_late_ethanol_death_probe`
5. `synthetic_viable_biomass_step`
6. `synthetic_ethanol_inhibition_challenge`
7. `synthetic_high_sugar_reference`
8. `synthetic_low_yan_ladder`
9. `natural_cold_hot_switch`

If the campaign must be reduced below 9, remove one of the two ethanol-stress designs first. Prefer retaining `synthetic_late_ethanol_death_probe` if the priority is `Kd0`/dead-biomass behavior, and retaining `synthetic_ethanol_inhibition_challenge` if the priority is ethanol inhibition `iE`.

## Sampling recommendation

Use the balanced sampling policy when possible: 10:00, 12:00, 14:00, and 16:00 on working days. The overnight benchmark showed balanced sampling outperformed front-loaded and two-per-day sampling in FIM quality and weak-direction variance reduction.

## Remaining caution

The model still has practical-identifiability weakness in `m0`, `qXG`, and several one-sided profile directions. The campaign should be sequential:

1. Run the first block of 3.
2. Refit and recompute profile likelihood/FIM.
3. Re-rank the remaining candidates before running blocks 2 and 3.
