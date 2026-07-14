# Fermentation modeling workspace

This directory contains the fermentation-model calibration, validation and
model-based design of experiments (MBDoE) work. The organization below marks
which artifacts are current, which are active dependencies, and which are kept
only for traceability.

## Authoritative workflows

### Pilot 2025: current integrated model

The current pilot-scale endpoint is:

- Runner: `pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`
- Source notebook: `pilot_2025/pilot_2025_co2_solubility_integrated_doe.ipynb`
- Executed notebook: `pilot_2025/pilot_2025_co2_solubility_integrated_doe.executed.ipynb`
- Results: `pilot_2025/results/co2_solubility_integrated_doe/`

This workflow combines the reduced primary fermentation model, glycerol,
secondary metabolites, aromas, liquid-gas partition, dissolved CO2, macro-O2,
estimability diagnostics and a natural-must MBDoE campaign.

Run from the repository root:

```powershell
python fermentation_model\pilot_2025\run_pilot_2025_co2_solubility_integrated_doe.py
```

For a structural/calibration smoke run without recomputing FIM and DOE:

```powershell
python fermentation_model\pilot_2025\run_pilot_2025_co2_solubility_integrated_doe.py --skip-fim --skip-doe
```

### Laboratory campaign: operational design

The latest volume-constrained laboratory campaign is represented by:

- `fermentation_final_operational_doe_volume_constrained.ipynb`
- `fermentation_final_operational_doe_volume_constrained.executed.ipynb`
- `run_final_operational_doe_volume_constrained.py`
- `results/final_operational_doe_volume_constrained/`

`fermentation_final_operational_doe_v2.*` and
`run_final_operational_doe_v2.py` remain in the root because the
volume-constrained and Lot 1 workflows use them as an active baseline.

### Executed Lot 1 follow-up

The current Lot 1 data and estimability follow-up is represented by:

- `fermentation_lot1_data_preview.ipynb`
- `fermentation_estimability_old_vs_lot1.ipynb`
- `run_lot1_actual_mbdoe_reassessment.py`
- `run_lot1_pulse_timing_mbdoe.py`
- `run_lot1_express_optimal_sampling.py`
- `run_estimability_old_vs_lot1.py`

## Active shared model modules

These root-level files are dependencies of the current pilot or laboratory
workflows and must not be treated as obsolete:

- `new_must_data_loader.py`: normalized synthetic/natural must ingestion.
- `run_new_must_glycerol_estimability_doe.py`: primary model and glycerol base.
- `run_secondary_metabolite_data_review.py`: secondary-state data preparation.
- `run_secondary_joint_campaign_doe.py`: joint secondary/aroma formulation.
- `run_secondary_v2_model_evaluation.py`: selected secondary-model utilities.
- `aroma_partition_unifac.py`: liquid-gas partition calculations.

The pilot-specific dependency layer is in `pilot_2025/support/`.

## Data

Raw workbooks and sensor exports remain under `data/`, separated by source and
year:

- `data/Laboratorio 2025/`
- `data/Laboratorio 2026/`
- `data/Piloto 2025/`
- `data/Piloto 2026/`

Raw data should not be moved into `results/` or modified by analysis scripts.
Normalized tables and fitted parameters belong in the corresponding results
directory.

## Results and legacy material

See `results/README.md` for the status of every result family. Some historical
results intentionally remain in place because current evidence-bundle scripts
read them by path.

Superseded notebooks and runners were moved to
`legacy/development_2026/`. They are retained as a frozen development record,
not as recommended entry points. See `legacy/README.md` and
`legacy/development_2026/ARCHIVE_MANIFEST.md`.

For the complete dependency and ownership map, see `REPOSITORY_MAP.md`.
