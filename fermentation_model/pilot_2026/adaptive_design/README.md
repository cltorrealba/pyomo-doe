# Pilot 2026 adaptive MBDoE workspace

This directory owns campaign-specific adapters, gate configuration and campaign
state for the 3-2-1 workflow. It does not currently contain an executable design
runner because the Phase B calibration gate failed and the reviewed design
constraints intentionally contain unresolved `null` safety bounds.

## Validated entry points

```bash
python fermentation_model/pilot_2026/adaptive_design/build_model_dataset.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_model_ready_qc_notebook.py
python fermentation_model/pilot_2026/adaptive_design/run_hierarchical_calibration.py
```

The first two commands pass against the current integrated data. The calibration
command exits with status 2 after writing an immutable `FAIL` gate run. It does
not fit parameters. This is deliberate: `calibration_config.json` lists the
scientific and owner prerequisites, and `design_constraints.json` lists every
missing operational bound that must fail closed.

No existing Pilot 2025 integrated result is accepted as an official prior. Any
future Pilot 2025 contribution must be regenerated coherently and restricted to
weak physical-structure support.

## Current boundary

- Phase A model-ready adapter: `PASS`.
- Phase B hierarchical calibration: `FAIL`; no fit executed.
- Phases C-E: not implemented or executed after the mandatory stop gate.
- Physical execution: not authorized; no executable schedule issued.
