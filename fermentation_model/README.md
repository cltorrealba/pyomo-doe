# Fermentation modeling workspace

This directory contains the four fermentation campaigns, their shared dynamic
model and the model-based design of experiments (MBDoE) workflows.

## Start here

1. Read `REPOSITORY_MAP.md` to locate a campaign, experiment, runner or result.
2. Use `campaigns/experiments.csv` as the master fermentation registry.
3. Use `campaigns/workflows.csv` to find the authoritative Python entry point
   and its output directory.
4. Read `REPRODUCIBILITY.md` before reporting new parameter estimates, FIM
   comparisons or profile-likelihood results.

## Four campaigns

| Campaign | Workspace | Raw data | State |
| --- | --- | --- | --- |
| Laboratory 2025 | `laboratory_2025/` | `data/Laboratorio 2025/` | complete CCD |
| Pilot 2025 | `pilot_2025/` | `data/Piloto 2025/` | complete paired vintage |
| Laboratory 2026 | `laboratory_2026/` | `data/Laboratorio 2026/` | active sequential MBDoE |
| Pilot 2026 | `pilot_2026/` | `data/Piloto 2026/` | data integration |

Raw files are immutable. Processed tables, parameters, figures and FIM outputs
belong to a workflow-owned results directory.

## Code and result ownership

- `shared/`: reusable data loading, fermentation/glycerol, secondary-metabolite,
  aroma and partition code. Shared priors and outputs are in `shared/results/`.
- `laboratory_2026/`: operational design, Lot 1 processing, estimability and
  sequential design. Outputs are in `laboratory_2026/results/`.
- `pilot_2025/`: pilot-specific integrated model and support layer. Outputs are
  in `pilot_2025/results/`.
- `pilot_2026/`: future processed pilot-2026 datasets and analyses.
- `legacy/`: frozen superseded science and administrative rendition artifacts.

The `fermentation_model/` root intentionally contains no active `.py` files and
no generic `results/` directory. This makes ownership visible from the path.

## Authoritative current endpoints

### Laboratory 2026 selected design

- Runner:
  `laboratory_2026/run_final_operational_doe_volume_constrained.py`
- Notebook:
  `laboratory_2026/notebooks/fermentation_final_operational_doe_volume_constrained.ipynb`
- Results:
  `laboratory_2026/results/final_operational_doe_volume_constrained/`
- Frozen execution handoff:
  `laboratory_2026/results/design_execution_bundle_2026-06-09/`

Run from the repository root:

```powershell
python fermentation_model\laboratory_2026\run_final_operational_doe_volume_constrained.py
```

### Pilot 2025 integrated model

- Runner: `pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py`
- Notebook: `pilot_2025/pilot_2025_co2_solubility_integrated_doe.ipynb`
- Results: `pilot_2025/results/co2_solubility_integrated_doe/`

```powershell
python fermentation_model\pilot_2025\run_pilot_2025_co2_solubility_integrated_doe.py
```

## Structural audit

From the repository root:

```powershell
python fermentation_model\tools\campaign_audit.py --check-hashes
python -m unittest discover -s fermentation_model\tests -v
```

The audit rejects active Python files or ambiguous results placed back in the
`fermentation_model/` root.
