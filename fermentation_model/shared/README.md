# Shared model layer

This is the reusable Python package imported by the laboratory and pilot
campaigns. Campaign folders must not contain copied working versions of these
modules.

## Modules

- `paths.py`: canonical data, workspace and result locations.
- `new_must_data_loader.py`: normalized synthetic/natural-must ingestion.
- `run_new_must_glycerol_estimability_doe.py`: primary fermentation and
  glycerol calibration/estimability layer.
- `run_new_must_overnight_validation.py`: retained validation support.
- `run_secondary_metabolite_data_review.py`: secondary-state preparation.
- `run_secondary_joint_campaign_doe.py`: joint secondary/aroma formulation.
- `run_secondary_v2_model_evaluation.py`: selected secondary-model utilities.
- `aroma_partition_unifac.py`: liquid-gas partition calculations.

Import shared modules explicitly:

```python
from shared import run_new_must_glycerol_estimability_doe as base
from shared import run_secondary_joint_campaign_doe as joint
```

## Results

`results/` contains only outputs owned by the shared model layer:

- `new_must_data_loading/`
- `identifiability_reduction/`
- `new_must_glycerol_estimability_doe/`
- `new_must_glycerol_overnight_validation/`
- `secondary_metabolite_data_review/`
- `secondary_joint_campaign_doe/`
- `secondary_v2_model_evaluation/`

Campaign outputs must not be written here. They belong under the corresponding
campaign workspace.

The preferred estimation family is local IPOPT calibration followed by weighted
FIM/eigen-analysis and profile likelihood. PSO and broad multistart artifacts
are frozen under `../legacy/development_2026/results/`.
