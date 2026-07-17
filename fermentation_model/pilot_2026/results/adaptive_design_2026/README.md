# Pilot 2026 adaptive-design results

Every child run directory is immutable. Earlier attempts are retained rather
than overwritten so changes in the adapter contract remain auditable.

## Model-dataset run ledger

| Run | Status | Interpretation |
| --- | --- | --- |
| `20260717T001101.548263Z_ced2c885b6` | FAIL | Windows path-length failure; failure manifest added, no deletion. |
| `20260717T001228Z_ced2c8` | PASS, superseded | Initial model-ready pass before the stricter censored-value projection. |
| `20260717T001546Z_ced2c8` | PASS, superseded | Censoring fixed; CO2 time origin still referenced the first valid sensor minute. |
| `20260717T002514Z_ced2c8` | PASS, current | Censoring fixed and CO2 time referenced to the authoritative sampling-window start. |

## Calibration-gate run ledger

| Run | Status | Interpretation |
| --- | --- | --- |
| `20260717T002128Z_62fd66` | FAIL, superseded | Gate evaluated against the superseded adapter run. No fit executed. |
| `20260717T002546Z_62fd66` | FAIL, current | Gate evaluated against current adapter run. No fit executed. |

The current calibration gate reports ten blockers, 53 raw-audit errors and one
raw-audit warning. Pyomo-DOE and IPOPT are available. No directory here contains
an executable fermentation schedule.
