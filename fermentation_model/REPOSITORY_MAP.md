# Repository map

Status date: 2026-07-14.

## Current dependency chain

```text
data/Laboratorio 2025-2026
        |
        v
new_must_data_loader.py
        |
        v
run_new_must_glycerol_estimability_doe.py
        |
        +--> run_secondary_metabolite_data_review.py
        |             |
        |             v
        +--> run_secondary_joint_campaign_doe.py
                      |
                      v
            run_secondary_v2_model_evaluation.py

data/Piloto 2025 + CO2 sensors
        |
        v
pilot_2025/support/pilot_2025_data_loader.py
        |
        +--> calibration_estimability prior
        +--> aroma_model_selection_doe prior
        +--> global model-selection prior
        |
        v
pilot_2025/run_pilot_2025_co2_solubility_integrated_doe.py
        |
        +--> dissolved CO2 and macro-O2 benchmark
        +--> integrated calibration
        +--> FIM/eigenvalue diagnostics
        +--> natural-must MBDoE
        |
        v
pilot_2025/results/co2_solubility_integrated_doe/
```

## Artifact status

| Area | Status | Authoritative location |
| --- | --- | --- |
| Pilot integrated model | Current | `pilot_2025/` |
| Pilot final results | Current | `pilot_2025/results/co2_solubility_integrated_doe/` |
| Pilot calibration/aroma/global priors | Active dependency | `pilot_2025/results/` |
| Laboratory volume-constrained DOE | Current | root notebooks and `results/final_operational_doe_volume_constrained/` |
| Lot 1 follow-up | Current | root Lot 1 notebooks/runners and `results/lot1_*` |
| Primary/glycerol model | Active shared dependency | root Python modules and `results/new_must_*` |
| Secondary/aroma model | Active shared dependency | root Python modules and `results/secondary_*` |
| Initial calibration and early DOE iterations | Archived | `legacy/development_2026/` |
| Pilot superseded branches | Archived | `pilot_2025/legacy/` |
| Evidence bundles | Frozen deliverables | `pilot_2025/bundles/` and selected root result bundles |

## Why some old-looking results remain in `results/`

The report and evidence-bundle generators still consume selected earlier
outputs, including curve validation, aroma campaign, secondary-model and
operational-design tables. Moving those directories would silently break the
evidence chain. They therefore remain in place and are labelled in
`results/README.md` instead of being relocated.

## Rules for future work

1. Add new pilot iterations under `pilot_2025/`; do not place pilot notebooks in
   the fermentation-model root.
2. Use one result directory per runner and keep its name equal to the workflow
   name.
3. Keep one source notebook and, when needed, one `.executed.ipynb` beside it.
4. When a workflow is superseded, move its notebooks and runner to `legacy/`
   only after checking imports with `rg`.
5. Do not delete prior parameter tables or FIM outputs used by a current runner
   or evidence bundle.
6. Record the authoritative runner, notebook and result directory in the
   nearest `README.md`.
