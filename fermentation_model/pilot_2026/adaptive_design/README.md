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

The result is a reduced baseline prior, not final model validation. Owner review,
is complete. Independent Ultra audit and aroma calibration remain release conditions.
No existing Pilot 2025 result is an official prior, no profile is approved, and
physical execution is not authorized.
