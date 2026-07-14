# Campaign registry

This directory is the navigation and provenance layer for the four experimental
campaigns. It does not duplicate raw data or model code.

- `campaigns.csv` defines the campaign-level scope and authoritative paths.
- `experiments.csv` is the master experiment registry. Unknown or incomplete
  metadata is recorded explicitly instead of being inferred silently.
- `workflows.csv` classifies the current Python entry points and separates
  shared model code from campaign-specific runners.
- `raw_data_manifest.csv` records file size and SHA-256 for the current raw-data
  snapshot. Regenerate it only when intentional raw inputs are added.

The campaign workspaces are:

- `../laboratory_2025/`
- `../pilot_2025/`
- `../laboratory_2026/`
- `../pilot_2026/`

Shared kinetic, secondary-metabolite and aroma code remains in the
`fermentation_model/` root during the transition. See `../shared/README.md`.

## Status vocabulary

- `complete`: campaign execution and core data ingestion are complete.
- `active`: experiments are being executed or incorporated.
- `data_integration`: experimental work exists, but the repository dataset is
  still being assembled or homologated.
- `planned_or_in_progress`: protocol exists but the raw experiment folder is not
  yet present in Git.
- `data_present_scope_unconfirmed`: source data exists, but its exact relationship
  to the campaign still needs confirmation.

Run the structural audit from the repository root:

```powershell
python fermentation_model\tools\campaign_audit.py
```

Regenerate the raw-data checksum manifest after deliberately adding source
files:

```powershell
python fermentation_model\tools\campaign_audit.py --write-raw-manifest
```
