# Independent Ultra audit — fixed-YAN Pilot 2026 calibration

Act as an independent scientific-data and scientific-code auditor. Work in
read-only mode on branch `ctorrealba_fermentation` at or after the commit that
contains calibration run `20260717T153612Z_b5b67b`. Do not modify files,
rebaseline manifests, or approve physical profiles.

Re-run:

```bash
python fermentation_model/pilot_2026/run_data_integration.py
python fermentation_model/pilot_2026/adaptive_design/build_model_dataset.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_model_ready_qc_notebook.py
python fermentation_model/pilot_2026/adaptive_design/run_hierarchical_calibration.py
python -m unittest discover -s fermentation_model/tests -v
python fermentation_model/tools/campaign_audit.py --check-hashes
```

Confirm independently:

1. The owner approved the loading/QC notebook and the active-fermentation versus
   cooling separation on 2026-07-17.
2. The integration tables retain the mass-derived 90.434783 mg/L YAN value for
   auditability, but calibration fixes every later historical pulse at exactly
   80 mg/L. `yan_delivery_fraction` must not be estimated or appear in the
   covariance/FIM.
3. Oculyze cell-to-biomass scale remains an estimated nuisance parameter.
4. Cooling and postprocess rows are absent from every kinetic likelihood.
5. Pilot Lot 1 CO2 is absent from calibration; Lots 2-3 use residual-derived
   AR(1) effective sample sizes.
6. The fixed-YAN calibration estimates 14 parameters, has no active parameter
   bounds, and records alternative primary minima rather than hiding them.
7. Recomputed results are consistent with computational `PASS`, approximately
   60.38% primary objective improvement, 83.70% CO2 improvement, and total CO2
   ESS near 469.53. Explain any numerical differences.
8. The MBDoE prior must include the multistart ensemble, covariance, empirical
   discrepancy scales, and nuisance marginalization; it must not use only the
   best local vector.
9. Current Pilot 2026 aroma evidence is reported accurately: wine has 39 values
   for each priority analyte; condensate ethyl acetate is entirely left-censored,
   while isoamyl acetate and ethyl octanoate contain observed and censored data.
10. No aroma design or physical schedule is released before the aroma-parameter
    gate and hybrid-engine qualification pass.

Return `PASS`, `PASS CONDITIONAL`, or `FAIL`, followed by a reproducible finding
matrix with severity, file/row evidence, scientific impact, and required action.
