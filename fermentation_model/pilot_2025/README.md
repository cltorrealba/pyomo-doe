# Pilot 2025 workspace

## Current integrated workflow

The authoritative pilot-scale workflow is:

- Runner: `run_pilot_2025_co2_solubility_integrated_doe.py`
- Source notebook: `pilot_2025_co2_solubility_integrated_doe.ipynb`
- Executed notebook: `pilot_2025_co2_solubility_integrated_doe.executed.ipynb`
- Final results: `results/co2_solubility_integrated_doe/`

The executed notebook is the current scientific narrative. It includes data
curation, the complete model equations, CO2-solubility and macro-O2 benchmark,
integrated fit, eigenvalue/FIM diagnostics and the selected natural-must MBDoE.

Run from the repository root:

```powershell
python fermentation_model\pilot_2025\run_pilot_2025_co2_solubility_integrated_doe.py
```

## Active dependencies

`support/` contains helper modules and prior-stage runners imported by the
current workflow. They are active code dependencies, even though their own
notebooks are no longer recommended endpoints.

The current runner also imports these shared modules from the parent directory:

- `run_new_must_glycerol_estimability_doe.py`
- `run_secondary_metabolite_data_review.py`
- `run_secondary_joint_campaign_doe.py`
- `run_secondary_v2_model_evaluation.py`
- `aroma_partition_unifac.py`

The following result families provide the calibrated starting point:

- `results/calibration_estimability/`
- `results/aroma_model_selection_doe/`
- `results/global_state_model_selection_doe/`
- `results/global_sugar/`

`global_sugar/` is the sugar-partition rerun retained under its historical
name. `global_state_model_selection_doe/` supplies the selected secondary-model
metadata. Both are therefore still active inputs.

## Bundles and legacy

- `bundles/`: frozen evidence packages with source snapshots and selected
  result tables.
- `legacy/notebooks/`: superseded exploratory and executed notebooks.
- `legacy/scripts/`: runners replaced by the integrated workflow.
- `legacy/results/`: outputs from superseded branches.

Do not move `support/` or the four active prior-result directories without
updating and retesting the current integrated runner.
