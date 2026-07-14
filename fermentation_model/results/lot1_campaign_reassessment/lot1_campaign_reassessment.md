# Lot 1 campaign reassessment

Generated from `fermentation_lot1_data_preview.executed.ipynb` processed outputs and the selected DOE bundle in `db_20260609/rec`.

## Scope and limitation

This is an adaptive reassessment using the real Lot 1 data, actual input log, and the previous Pyomo-DoE/FIM campaign metrics. The current runtime could not access the Python/Jupyter environment, so I did not rerun a full Pyomo sensitivity solve here. The recommendations below combine: actual Lot 1 composition/deviation, observed signal quality, and the existing D-opt/hybrid FIM diagnostics.

## Lot 1 execution signals

| Process | Candidate | Key observation | Consequence |
|---|---|---|---|
| F1 | synthetic_lit_SM410_18C_highN_ester | Actual YAN 312 mg/L vs 410 target; initial viable X 0.105 g/L vs 0.421 nominal; CO2 peak 9.23 sccm at 47.8 h. | Useful high-N/isothermal signal, but it does not fully replace a high-N target at warm strip conditions. |
| F2 | synthetic_high_biomass_low_N_maintenance | Low N was achieved (37 mg/L); initial X was 1.095 g/L vs 2.000 nominal; last sugar 26.8 g/L at 383.8 h. | Very valuable slow/maintenance trajectory, but actual biomass dose was lower than design. |
| F3 | synthetic_fructose_rich_glucose_pulse | Fructose-rich condition achieved, but N pulse was 0.015 vs planned 0.05; CO2 peak 13.68 sccm. | Useful sugar-separation/high-temperature signal, but the N-rescue part was under-excited. |

## Nominal DOE contribution still relevant

After the first three nominal experiments, weakest directions were: Kd0 (0.13 var-red), qEG (0.43 var-red), qN (0.43 var-red), mu0 (0.44 var-red), gammaG0 (0.48 var-red), betaG0 (0.49 var-red), iE (0.51 var-red), gammaF0 (0.56 var-red).

The remaining nominal design continues improving these directions. The largest remaining marginal FIM gain is still F4 (logdet +4.62, trace-inv -0.49), followed by F5 and F6. Gains after F6 are smaller but still useful for validation, natural transfer, and aroma/CO2 loss.

## Recommendation for immediate preparation

1. Keep Lot 2 unchanged: F4, F5, F6 should be prepared and executed next.
2. Add stricter execution QC before Lot 2: initial YAN check/correction, inoculum/biomass standardization, and exact pulse dose logging.
3. Reconsider only Lot 3 after Lot 2 data. The main adaptive change I would consider is replacing F7 (synthetic_high_sugar_reference) with a death/Kd0 probe if Kd0 is still intended to be estimated. If Kd0 will be fixed or strongly regularized, keep F7 or replace it with a second natural-matrix experiment depending on wine-transfer priority.

## Experiment-by-experiment decision

| Order | Candidate | Decision | Priority | Main reason |
|---:|---|---|---|---|
| 4 | synthetic_glucose_rich_fructose_pulse | KEEP | very_high | F3 gave fructose-rich information, but the planned N pulse was underdosed to 30% of target. The glucose-rich mirror is still needed to separate G/F uptake, yields and inhibition directions. |
| 5 | synthetic_lit_SM410_24C_highN_strip | KEEP_WITH_QC | very_high | F1 was high-N but actual YAN was 312 mg/L vs 410 mg/L target and was isothermal 18 C. This run is the main high-N warm/high-CO2 strip condition for aroma and CO2 loss. |
| 6 | synthetic_viable_biomass_step | KEEP_WITH_QC | high | Lot 1 measured initial viable X was far below nominal in F1/F3 and about 55% of nominal in F2. A deliberate biomass excitation remains valuable for X proportionality and yield separation. |
| 7 | synthetic_high_sugar_reference | CONDITIONAL | medium_low | After F1/F3 and planned F4/F5, another synthetic high-sugar reference is less critical than death/Kd0 or natural matrix transfer. |
| 8 | natural_glucose_pulse_after_growth | KEEP | very_high | Lot 1 is fully synthetic. At least one natural-must perturbation is required to test transferability to wine-relevant matrix behavior. |
| 9 | synthetic_fast_CO2_aroma_strip_highN | KEEP | high | F5 and F9 are related but not identical; this run targets fast CO2/aroma stripping and terminal condensate information. |

## Operational notes before Lot 2

- Use actual event times and actual measured initial states in calibration; do not force nominal design states.
- For high-N recipes, do not assume recipe equals measured YAN. F1 suggests a systematic shortfall risk.
- For biomass, Oculyze-derived viable X at t0 is the state that should enter the model. The nominal initial X was not achieved in Lot 1.
- F2 should be retained as an informative slow/low-N trajectory, not treated as failed data.
- If a fermentation is called dry by density but Y15 still reports residual fructose, keep the chemical value in the model and record the operational dry criterion separately.

## Files

- `lot1_actual_vs_design.csv`
- `lot1_signal_features.csv`
- `doe_incremental_nominal.csv`
- `remaining_experiment_reassessment.csv`
