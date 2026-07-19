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

The repository audit still reports the historical 53 errors: 52 raw size
divergences and one ambiguous generated-results owner. They are classified
without changing raw data or the manifest. The current calibration
estimates 15 primary, YAN/Oculyze and CO2 parameters; five selected-basin
replicas and five residual-weighted CO2 starts converge. Broad starts retain
alternative primary minima, so the result is a conditional reduced prior. No
directory here contains an executable fermentation schedule.

## Corrected Wave-1 qualification ledger

| Run | Status | Interpretation |
| --- | --- | --- |
| `wave1_mbdoe_adapter_requalification/20260718T011943Z_7afa2a` | FAIL, retained | Windows path-length failure recorded before qualification. |
| `wave1_mbdoe_adapter_v2/20260718T012044Z_7afa2a` | FAIL, retained | Second path-length failure recorded before the extended-path fix. |
| `wave1_mbdoe_adapter_v2/20260718T012212Z_7afa2a` | PASS, superseded | Frozen-scale FIM, causal pulse semantics, finite-difference/grid checks and independent Pyomo-DOE check pass. |
| `wave1_hybrid_search_v2/20260718T014738Z_70e343` | REJECTED, retained | Search was superseded after full sampling exposed actions after biological drying. |
| `wave1_sampling_v2/20260718T015547Z_70e343` | FAIL, retained | Correctly rejected the superseded search and unresolved capture/conflict rules. |
| `wave1_hybrid_search_v2/20260718T020952Z_70e343` | PASS CONDITIONAL, superseded | Three Sobol-PSO seeds and accepted trust-region IPOPT refinement; ranking stability and one tracking-error actuator scenario fail. |
| `wave1_sampling_v2/20260718T021159Z_70e343` | FAIL, superseded | Full 64-member optimization improves information but capture and conflict approvals remained unresolved. |
| `wave1_mbdoe_adapter_v3/20260718T061321Z_b41bbd` | PASS, current | Eight-case policy/member/scenario scope, corrected FIM, t=0 baseline, nine capture intervals and independent direct-versus-cached 9x9 FIM all pass. |
| `wave1_hybrid_search_v2/20260718T074703Z_b41bbd` | FAIL, current | Five-seed minimum search and complete approved actuator envelope run; convergence, cross-seed finalist stability, plateau and local-refinement qualification fail closed. |
| `wave1_sampling_v2/20260718T075319Z_b41bbd` | FAIL, retained | Complete candidate; visual QA found title/subtitle overlap in three multipanel figures. |
| `wave1_sampling_v2/20260718T075725Z_b41bbd` | FAIL, retained | Title overlap fixed; visual QA found top sample labels entering the agenda title area. |
| `wave1_sampling_v2/20260718T080054Z_b41bbd` | FAIL, current | Operational sampling/capture/randomization checks pass and all ten figures pass visual QA; source search and missing nutrition product compositions remain fail-closed. |
| `wave1_final_search/20260718T211306Z_cf5179` | PASS, current | Resumable five-seed multifidelity search; 155 checkpoints, real 64-member local acceptance, corrected eligible-set selection, practical convergence, full actuator envelope and exact post-search FIM pass. |
| `wave1_final_sampling/20260719T042832Z_cf5179` | PASS, current | Optimized ten-sample schedules and nine capture intervals pass the full 64-member tail guardrails; all ten watermarked figures pass automated and manual visual QA. |
| `wave1_final_execution_package/20260719T043206Z_cf5179` | FAIL, retained | Package files were created, but strict JSON serialization rejected NumPy boolean audit values before a gate or manifest was emitted. |
| `wave1_final_execution_package/20260719T043339Z_cf5179` | PASS, superseded | Operational checks passed, but independent review found that Windows `MAX_PATH` enumeration omitted package files from the run manifest. |
| `wave1_final_execution_package/20260719T043613Z_cf5179` | PASS, current | Operational translation and automated audit pass; all 23 package files and 26 declared run outputs are hashed. Physical authorization remains pending owner approval. |

The immutable historical regression references are
`wave1_hybrid_search/20260717T225010Z_9d52be` and
`wave1_sampling/20260717T231040Z_9d52be`. The corrected runs do not overwrite
them. The campaign state has no tank assignments, no executable profile and no
physical release. The tank CSV freezes the owner-approved logical mapping for
review only and is never copied into `tank_assignments`.
