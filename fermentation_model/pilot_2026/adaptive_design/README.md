# Pilot 2026 adaptive MBDoE workspace

This directory owns the campaign-specific model adapter, calibration, gates and
3-2-1 campaign state. It does not issue physical temperature or nutrition
profiles.

## Reproducible entry points

Run from the repository root in `environment-maintainer.yml`:

```bash
python fermentation_model/pilot_2026/adaptive_design/resolve_historical_audit.py
python fermentation_model/pilot_2026/adaptive_design/build_model_dataset.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_model_ready_qc_notebook.py
python fermentation_model/pilot_2026/adaptive_design/run_hierarchical_calibration.py
python fermentation_model/pilot_2026/adaptive_design/run_aroma_calibration.py
python fermentation_model/pilot_2026/adaptive_design/build_joint_ensemble.py
python fermentation_model/pilot_2026/adaptive_design/qualify_hybrid_engine.py
python fermentation_model/pilot_2026/adaptive_design/fit_temperature_actuator.py
```

The corrected Wave-1 steps require explicit immutable sources; they never
discover a source by directory ordering:

```bash
python fermentation_model/pilot_2026/adaptive_design/validate_wave1_mbdoe_adapter.py --source-ensemble-run fermentation_model/pilot_2026/results/adaptive_design_2026/joint_ensemble/20260717T162323Z_6b9889 --source-engine-run fermentation_model/pilot_2026/results/adaptive_design_2026/hybrid_engine_qualification/20260717T164729Z_c7a6d5 --source-aroma-run fermentation_model/pilot_2026/results/adaptive_design_2026/aroma_calibration/20260717T161527Z_e927a5 --source-actuator-run fermentation_model/pilot_2026/results/adaptive_design_2026/temperature_actuator/20260717T170749Z_95efe0
python fermentation_model/pilot_2026/adaptive_design/run_wave1_hybrid_search.py --source-adapter-run fermentation_model/pilot_2026/results/adaptive_design_2026/wave1_mbdoe_adapter_v2/20260718T012212Z_7afa2a --source-ensemble-run fermentation_model/pilot_2026/results/adaptive_design_2026/joint_ensemble/20260717T162323Z_6b9889 --source-aroma-run fermentation_model/pilot_2026/results/adaptive_design_2026/aroma_calibration/20260717T161527Z_e927a5 --source-actuator-run fermentation_model/pilot_2026/results/adaptive_design_2026/temperature_actuator/20260717T170749Z_95efe0
python fermentation_model/pilot_2026/adaptive_design/optimize_wave1_sampling_and_plots.py --source-search-run fermentation_model/pilot_2026/results/adaptive_design_2026/wave1_hybrid_search_v2/20260718T020952Z_70e343 --source-ensemble-run fermentation_model/pilot_2026/results/adaptive_design_2026/joint_ensemble/20260717T162323Z_6b9889 --source-aroma-run fermentation_model/pilot_2026/results/adaptive_design_2026/aroma_calibration/20260717T161527Z_e927a5
```

The historical 53-error audit is not rebaselined: it remains 52 raw worktree
size divergences plus one ambiguous generated-results owner. The raw-data
manifest and protected source files are unchanged. Immutable text artifacts
also record their Windows checkout bytes and accept canonical LF verification
only when that canonical digest exactly matches the declared source digest.

## Scientific contract

- Model likelihood tables contain active fermentation only. Cooling remains in
  integration QC and is excluded from kinetic likelihoods.
- The mass-derived 90.435 mg/L pulse remains auditable in integration outputs.
  The owner-approved kinetic basis is fixed at 80 mg/L for each later
  historical pulse; no YAN delivery fraction is estimated. The time-zero pulse
  is already represented by measured initial YAN.
- Oculyze cells are mapped to viable/dead biomass through an estimated nuisance
  conversion with a bounded prior, not a confirmed physical constant.
- Observation errors have explicit absolute and relative terms. Nominal scales
  are frozen while sensitivities are calculated; state-specific model
  discrepancy propagates into posterior/design uncertainty.
- Wine and condensate censoring operators and scales are explicit. Aroma
  parameters outside the reduced Stage-B calibration remain broad in the joint
  ensemble.
- CO2 is fit only for Lots 2-3. Effective sample size is recomputed per run from
  final standardized residuals before residual-weighted refitting.
- Pulses are causal state jumps. At coincident events, the declared
  sample-before-action or action-before-sample convention is applied exactly.

## Calibration and ensemble status

The fixed-YAN reduced calibration in
`baseline_calibration/20260717T153612Z_b5b67b` is computational `PASS` and
release `PASS_CONDITIONAL`. Ten primary fits converged, the selected basin was
replicated, broad starts retained alternative minima, and the residual-weighted
CO2 refit improved its objective. It is a reduced prior, not final model
validation.

The aroma calibration in
`aroma_calibration/20260717T161527Z_e927a5` is `PASS_CONDITIONAL`. Stationary
formation is identified for the three priority aromas; ethyl-acetate loss and
growth-associated formation remain weak and therefore broad.

The deterministic 64-member joint ensemble in
`joint_ensemble/20260717T162323Z_6b9889` passed its gate and represents all ten
successful primary multistart centres. Point-estimate-only design is prohibited.
The Windows environment passed the real Pyomo/IPOPT benchmark in
`hybrid_engine_qualification/20260717T164729Z_c7a6d5`; SciPy is not a final
local solver.

## Corrected Wave-1 outcome

The adapter qualification in
`wave1_mbdoe_adapter_v2/20260718T012212Z_7afa2a` is `PASS`. It freezes nominal
observation scales during finite differencing, passes three sensitivity-step
and three integration-grid checks, and agrees with an independent Pyomo-DOE
linear FIM to relative error 1.71e-15.

The multiseed Sobol-PSO plus trust-region IPOPT search in
`wave1_hybrid_search_v2/20260718T020952Z_70e343` is `PASS_CONDITIONAL`. All
three seeds produced finite candidates, the accepted IPOPT proposal improved
the actual 64-member objective, and all nominal actions precede biological
drying. Ranking stability and all encoded actuator scenarios did not pass.

The full-ensemble sampling run in
`wave1_sampling_v2/20260718T021159Z_70e343` improved the robust information
score from 20.520 to 22.172, with 63 paired wins and one loss relative to the
preliminary schedule. Its final gate is nevertheless `FAIL`: capture
interval/loading/change rules and manual sampling conflict rules are not
approved. Nutrition product combinations are reported only as exact endpoint
translations; no product mix is selected.

The historical regression references
`wave1_hybrid_search/20260717T225010Z_9d52be` and
`wave1_sampling/20260717T231040Z_9d52be` remain immutable. All current results
are computational candidates only. There is no physical temperature/nutrition
profile, tank assignment, or execution authorization.
