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
```

The runner reads immutable files under `../data/Piloto 2026/Lotes1a3/raw/`
and regenerates `results/data_integration_2026/`. It never modifies source
workbooks, ZIP archives, MDB files or the original `Resultados_GC.csv` export.

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
  the reported numeric value multiplied by 1000. NQ remains categorical.
- Condensate A/B volumes come from `09_Aromas_Cond`. Lot 1 condensates were not
  prepared or analyzed and remain outside the GC model scope.

## Calibration gate

The integration verdict is `processed_qc_ready_conditional_for_calibration`.
Before parameter estimation or optimal design, complete the independent Ultra
verification using `results/data_integration_2026/ULTRA_VERIFICATION_PROMPT.md`
and confirm the factory normal-reference temperature and pressure configured in
MassView (20 °C is currently recorded only as the owner's estimate).
