# Pilot 2026 adaptive-design results

Every child run directory is immutable. Earlier attempts are retained rather
than overwritten so changes in the adapter contract remain auditable.

## Model-dataset run ledger

| Run | Status | Interpretation |
| --- | --- | --- |
| `20260717T001101.548263Z_ced2c885b6` | FAIL | Windows path-length failure; failure manifest added, no deletion. |
| `20260717T001228Z_ced2c8` | PASS, superseded | Initial model-ready pass before the stricter censored-value projection. |
| `20260717T001546Z_ced2c8` | PASS, superseded | Censoring fixed; CO2 time origin still referenced the first valid sensor minute. |
| `20260717T002514Z_ced2c8` | PASS, superseded | Censoring fixed and CO2 time referenced to the authoritative sampling-window start. |
| `20260717T125004Z_b03413` | PASS, current | Active-only kinetic tables; cooling and the one out-of-process event are formally excluded. |

## Calibration-gate run ledger

| Run | Status | Interpretation |
| --- | --- | --- |
| `20260717T002128Z_62fd66` | FAIL, superseded | Gate evaluated against the superseded adapter run. No fit executed. |
| `20260717T002546Z_62fd66` | FAIL, superseded | Historical gate with 53 worktree errors. No fit executed. |
| `20260717T125012Z_c17f74` | FAIL, superseded | Dependency preflight failure before optimization; no scientific result. |
| `20260717T125123Z_c17f74` | FAIL, superseded | First executed fit; broad primary basin was not replicated. |
| `20260717T125412Z_0743f8` | PASS CONDITIONAL, superseded | Selected basin and residual CO2 fit replicated; error-scale contract still incomplete. |
| `20260717T125732Z_0743f8` | PASS CONDITIONAL, current | Validated reduced fit with empirical model-discrepancy scales and residual-based CO2 ESS. |

The current audit has zero errors and retains one unrelated warning for the
missing Laboratory 2026 Lot 3 raw folder. The historical 53 errors are
classified without changing raw data or the manifest. The current calibration
estimates 15 primary, YAN/Oculyze and CO2 parameters; five selected-basin
replicas and five residual-weighted CO2 starts converge. Broad starts retain
alternative primary minima, so the result is a conditional reduced prior. No
directory here contains an executable fermentation schedule.
