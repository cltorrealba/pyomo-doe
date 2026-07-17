# Pilot 2026 adaptive MBDoE workspace

This directory owns the campaign-specific model adapter, calibration, gates and
3-2-1 campaign state.  It does not yet issue temperature or nutrition profiles.

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
python fermentation_model/pilot_2026/adaptive_design/validate_wave1_mbdoe_adapter.py
python fermentation_model/pilot_2026/adaptive_design/run_wave1_hybrid_search.py
python fermentation_model/pilot_2026/adaptive_design/optimize_wave1_sampling_and_plots.py
```

The historical 53-error artifact is not rebaselined.  It is classified as 52
raw worktree size divergences plus one generated aggregate result directory.
A clean committed tree passes the current hash audit with the original raw
manifest byte-for-byte unchanged.

## Scientific contract

- Model tables contain active fermentation only. Cooling remains in integration
  QC and is formally absent from every kinetic likelihood.
- The mass-derived pulse of 90.434783 mg/L remains in integration outputs for
  auditability.  The owner-approved kinetic input is fixed at 80 mg/L for each
  later historical pulse; no YAN delivery fraction is estimated.  The time-zero
  pulse is already represented by measured initial YAN.
- Oculyze cells are mapped to viable/dead biomass through an estimated nuisance
  conversion with bounded prior. It is not presented as a confirmed physical
  constant.
- Primary observation errors have explicit absolute and relative terms.  The
  calibration also records state-specific empirical model-discrepancy
  multipliers; these must propagate into posterior/design uncertainty.
- Wine and condensate left-censoring operators and scales are explicit. Aroma
  kinetic parameters are outside the reduced Stage-B calibration.
- CO2 is fit only for Lots 2-3. Effective sample size is recomputed per run from
  final standardized residuals, followed by a residual-weighted refit.

## Current outcome

The owner-approved fixed-YAN run under
`results/adaptive_design_2026/baseline_calibration/20260717T153612Z_b5b67b`
has computational verdict `PASS` and release verdict `PASS_CONDITIONAL`:

- 10 primary fits converged: five broad starts and five selected-basin replicas;
- three selected-basin replicas reproduce the optimum within 5%;
- broad starts reveal alternative minima, which remain in the ensemble;
- primary objective improvement is 60.38%;
- all five residual-weighted CO2 starts reproduce one basin;
- CO2 objective improvement is 83.70%;
- residual CO2 ESS totals 469.53 over six runs;
- no estimated parameter is bound-active.

The result is a reduced baseline prior, not final model validation. Owner review
is complete.

The aroma calibration under
`results/adaptive_design_2026/aroma_calibration/20260717T161527Z_e927a5` has
verdict `PASS_CONDITIONAL`.  Formation during the stationary-associated regime is
identified for all three priority aromas.  Effective loss is separately identified
for ethyl octanoate and isoamyl acetate, but not for ethyl acetate because all 33
condensate observations are left-censored.  Growth-associated formation is weak
for all three analytes.  These directions must remain broad in the robust MBDoE
ensemble and are explicit Wave-1 information targets.

Independent Ultra audit and PSO-to-IPOPT engine qualification remain release
conditions. No existing Pilot 2025 result is an official prior, no profile is
approved, and physical execution is not authorized.

## Phase C status

The deterministic 64-member joint ensemble under
`results/adaptive_design_2026/joint_ensemble/20260717T162323Z_6b9889` passed its
gate. It represents all ten successful primary multistart centres and inflates
only the aroma directions declared weak by the bilateral profile-likelihood
gate. Point-estimate-only design is prohibited.

The hybrid-engine qualification under
`results/adaptive_design_2026/hybrid_engine_qualification/20260717T162820Z_c7a6d5`
is deliberately `FAIL` in the current container: PSO is reproducible and its
SciPy diagnostic refinement solves the bounded Rosenbrock benchmark, but no
Pyomo/IPOPT or cyipopt interface is installed. SciPy is not accepted as a final
substitute. Run the qualification command above in the approved environment
where IPOPT is available; candidate generation stays locked until that gate and
the real MBDoE objective adapter both pass.

The approved Windows environment subsequently passed the real IPOPT benchmark
in `hybrid_engine_qualification/20260717T164729Z_c7a6d5`.  The Wave-1 adapter
then passed after adding a first-order setpoint-to-Sonda1 actuator model fitted
from all nine Pilot 2026 runs (global step-response tau = 0.318 h).

The current preliminary pair under
`wave1_hybrid_search/20260717T171902Z_9d52be` was evaluated against all 64 joint
ensemble members. It reaches median information gain 7.89, lower-decile gain
5.22, 100% completion by the three-week horizon and maximum terminal residual
sugar 3.46 g/L. Its ten-sample-per-process schedules and review figures are in
`wave1_sampling/20260717T172453Z_9d52be`. The executable sampling gate is
`PASS_CONDITIONAL`: samples are restricted to weekdays 09:00–17:00 and each has
at least 95% probability of preceding biological drying.

These are deliberately high-excitation computational candidates, not physical
instructions. The same Wave-1 search must still be rerun in the approved IPOPT
environment so the local trust-region surrogate is solved and revalidated
against the actual MBDoE objective. Owner plot review and independent Ultra
audit remain mandatory before tank assignment or physical release.
