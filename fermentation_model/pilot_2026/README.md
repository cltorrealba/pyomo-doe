# Pilot 2026 campaign

## Scope

Nine pilot fermentations were executed as three lots of three reactors using
natural Sauvignon blanc must, constant nutrition and dynamic temperature
profiles. The reactor/run mapping is confirmed in
`../campaigns/experiments.csv` and `data_integration_config.json`.

## Authoritative integration

Run from the repository root:

```bash
cd fermentation_model/pilot_2026
npm install
cd ../..
python fermentation_model/pilot_2026/run_data_integration.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_data_loading_qc_notebook.py
```

The runner reads immutable files under `../data/Piloto 2026/Lotes1a3/raw/`
and regenerates `results/data_integration_2026/`. It never modifies source
workbooks, ZIP archives, MDB files or the original `Resultados_GC.csv` export.
The executed notebook provides the owner-facing visual review. The next-stage
implementation brief is `ADAPTIVE_HYBRID_DESIGN_CODEX_PROMPT.md`.

## Confirmed scientific rules

- The sampling window is the inclusive interval between the first and last
  primary-result sample in the master workbook.
- Lot 2 dates displayed as June through December are corrected to 6–12 April.
- MassView `fmeasure` is already normalized and is retained in Ln/min without
  an additional temperature/pressure conversion.
- Raw values exactly equal to 2.621 are invalid and are removed before
  same-second averaging.
- Leading/trailing zeros are valid; zeros between the first and last positive
  active-process observation are invalid.
- Missing acquisition seconds are filled with zero only at the two edges and
  are explicitly marked artificial. Internal gaps are not imputed.
- A controller `IMPOSTA_CTA` event with mode 2 and a 9–11 °C band marks the
  start of cooling. Cooling/postprocess observations are excluded from dynamic
  calibration.
- Protocol C is 16→18→21 °C. The 21 °C stage was not observed inside the
  calibration window for 26159 and is not represented as an executed input.
- GC numeric results describe the 1:1000 diluted vial. The MIX concentration is
  the reported numeric value multiplied by 1000; wine MEF results remain on
  their undiluted basis. NQ and below-LOQ values are retained as left-censored
  observations.
- Condensate A/B volumes come from `09_Aromas_Cond`. Lot 1 condensates were not
  prepared or analyzed and remain outside the GC model scope.
- The 39 pilot wine GC samples comprise six initial baselines without
  condensate and 33 samples paired one-to-one with the 33 MIX samples.
- MIX observations are accumulated masses over collection intervals, computed
  from corrected MIX concentration and recovered A+B volume.
- Lot 3 operational events are reconstructed by shifting the owner-confirmed
  matching Lot 1 schedule to each Lot 3 process start.
- Eight `-0.3 % v/v` Alcolyzer formula artifacts with blank raw readings are
  set to zero in processed data; the original value and correction flag remain.
- Lot 1 CO2 is retained for QC plots but excluded from model calibration because
  the owner confirmed that those three signals are unreliable.
- `Sonda1` is the authoritative measured temperature. Controller setpoint is
  retained as the commanded input.
- The mass-derived audit calculation for a listed `116 + 46 g` pulse is
  90.435 mg/L YAN at 230 L. It remains in the integration evidence, but it does
  not replace the owner-approved 80 mg/L-per-pulse historical calibration basis.

## Calibration gate

The integration verdict is `processed_qc_ready_conditional_for_calibration`.
Before parameter estimation or optimal design, review the executed loading-QC
notebook, complete the independent verification using
`results/data_integration_2026/ULTRA_VERIFICATION_PROMPT.md`, and document the
factory normal-reference temperature and pressure configured in MassView
(20 °C is currently recorded only as the owner's estimate).

## Model-ready adapter and adaptive campaign gate

Build one immutable model-dataset run and execute its reader-facing QC notebook:

```bash
python fermentation_model/pilot_2026/adaptive_design/build_model_dataset.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_model_ready_qc_notebook.py
```

The adapter writes only under
`results/adaptive_design_2026/model_dataset/<run-id>/`. It keeps observation
operators and units explicit, excludes cooling from the kinetic tables, prevents
total/individual sugar double counting, retains GC censoring and capture
intervals, excludes Lot 1 CO2 by assertion, and applies deterministic
10-minute CO2 bins. Signal-based weights are provisional; calibration replaces
them with residual-based effective-sample weights. Density-triggered
second pulses carry data-derived timing intervals. The carbon balance is
explicitly diagnostic-only, not a closed elemental balance.

Evaluate the calibration gate with:

```bash
python fermentation_model/pilot_2026/adaptive_design/run_hierarchical_calibration.py
```

The reconstructed 2026-07-17 reduced calibration executes five broad primary
multistarts, five selected-basin replicas, and two-pass CO2 fitting with ESS
computed from final standardized residuals. Its computational verdict is
`PASS` and its release verdict is `PASS_CONDITIONAL`. The result uses the fixed
owner-approved historical YAN basis, parameterizes the Oculyze cell-to-biomass
conversion, estimates model-discrepancy scales, and preserves alternative broad
minima. It is suitable as a conditional reduced prior for the computational
MBDoE stage, not as a final physical schedule.

The corrected Wave-1 adapter and multiseed PSO-to-IPOPT search have been
qualified computationally. The subsequent 64-member sampling qualification is
`FAIL`: capture limits and operational conflict rules remain unapproved.
`adaptive_design/design_constraints.json` records those values as `null` and
fail-closed. No executable profile or tank assignment has been issued.
