# Repository map

Status date: 2026-07-16.

This is the canonical navigation map for the fermentation work. Every active
artifact now has one explicit owner: the shared model layer or one experimental
campaign. There is no generic `fermentation_model/results/` directory and no
active Python module in the `fermentation_model/` root.

## Directory tree

```text
fermentation_model/
├── README.md                  orientation and entry points
├── REPOSITORY_MAP.md          this ownership/dependency map
├── REPRODUCIBILITY.md         environment and run-manifest rules
├── campaigns/                 machine-readable campaign/experiment registry
├── config/                    solver and run-manifest configuration
├── data/                      immutable experimental sources
├── shared/                    reusable model package and shared results
├── laboratory_2025/           Laboratory 2025 campaign metadata
├── pilot_2025/                Pilot 2025 code, notebooks and results
├── laboratory_2026/           Laboratory 2026 code, notebooks and results
├── pilot_2026/                Pilot 2026 data-integration workspace
├── legacy/                    frozen superseded work and rendition bundles
├── tools/                     repository audit and run-context capture
└── tests/                     structural and provenance checks
```

## Where to start

| Question | Start at | Then follow |
| --- | --- | --- |
| What campaigns exist? | `campaigns/campaigns.csv` | campaign `README.md` |
| Where is one fermentation? | `campaigns/experiments.csv` | its `raw_location` |
| Which runner is authoritative? | `campaigns/workflows.csv` | `path` and `authoritative_output` |
| Are raw files unchanged? | `campaigns/raw_data_manifest.csv` | `tools/campaign_audit.py --check-hashes` |
| How do I reproduce a result? | `REPRODUCIBILITY.md` | run config and `run_manifest.json` |
| Is an old artifact still current? | nearest `README.md` | otherwise treat `legacy/` as frozen |

## Ownership by layer

### Immutable data

```text
data/Laboratorio 2025/
data/Piloto 2025/
data/Laboratorio 2026/
data/Piloto 2026/
```

Analysis code never writes into `data/`. Derived normalized tables belong to
the results directory of the workflow that created them.

### Shared model

```text
shared/
├── paths.py
├── new_must_data_loader.py
├── run_new_must_glycerol_estimability_doe.py
├── run_new_must_overnight_validation.py
├── run_secondary_metabolite_data_review.py
├── run_secondary_joint_campaign_doe.py
├── run_secondary_v2_model_evaluation.py
├── aroma_partition_unifac.py
├── notebooks/
└── results/
```

`shared/paths.py` is the single path contract. Shared modules must import
canonical directories from it instead of deriving data or result ownership from
their own file location.

### Laboratory 2026

```text
laboratory_2026/
├── run_final_operational_doe_v2.py
├── run_final_operational_doe_volume_constrained.py
├── run_estimability_old_vs_lot1.py
├── run_lot1_actual_mbdoe_reassessment.py
├── run_lot1_pulse_timing_mbdoe.py
├── run_lot1_express_optimal_sampling.py
├── create_*.py
├── notebooks/
└── results/
```

The authoritative design handoff is
`laboratory_2026/results/design_execution_bundle_2026-06-09/`. Sequential Lot 1
processing and design updates live beside it under `laboratory_2026/results/`.

### Pilot campaigns

Pilot 2025 owns its integrated model, support code, notebooks and results under
`pilot_2025/`. Pilot 2026 has nine confirmed reactor/run mappings and owns the
lossless integration runner, QC tables, figures and executed loading-QC notebook
under `pilot_2026/`. Its corrected data contract includes paired wine/MIX aroma
observations, reconstructed Lot 3 operations and explicit CO2/model masks. The
processed dataset is conditional for calibration pending owner notebook review,
independent verification and confirmation of the MassView factory normal
reference.

### Frozen history

`legacy/development_2026/` contains superseded runners, notebooks and their
matching result families. `legacy/rendicion/` and `pilot_2025/bundles/` contain
administrative AXX deliverables. Frozen source snapshots retain their original
paths intentionally and are not imported by active workflows.

## Active dependency chain

```text
data/Laboratorio 2025 + data/Laboratorio 2026
                         │
                         ▼
             shared/new_must_data_loader.py
                         │
                         ▼
 shared/run_new_must_glycerol_estimability_doe.py
                         │
             ┌───────────┴───────────┐
             ▼                       ▼
 secondary metabolite review   shared validation
             │
             ▼
 secondary joint/aroma model
             │
             ▼
 secondary v2 model evaluation
             │
       ┌─────┴────────────────┐
       ▼                      ▼
 laboratory_2026 runners   pilot_2025 runners
       │                      │
       ▼                      ▼
 campaign-owned results   campaign-owned results
```

## Rules for future additions

1. Put reusable model code in `shared/`; put campaign-specific code in its
   campaign workspace.
2. Write outputs only to `shared/results/` or `<campaign>/results/`.
3. Add every experiment to `campaigns/experiments.csv`.
4. Add every active runner to `campaigns/workflows.csv`.
5. Regenerate the raw-data manifest only after intentionally adding raw files.
6. Move superseded code and its outputs together into `legacy/`.
7. Do not import code from notebooks, executed notebooks, bundles or `legacy/`.
8. Run `python fermentation_model/tools/campaign_audit.py --check-hashes`
   before committing a scientific result.
