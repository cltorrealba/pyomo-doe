# Prompt for Codex — adaptive hybrid MBDoE for Pilot 2026

Work in the `ctorrealba_fermentation` branch of `cltorrealba/pyomo-doe`, starting
from the current HEAD. Act as a scientific-computing engineer with expertise in
fermentation kinetics, parameter estimation, Pyomo, Pyomo-DOE, IPOPT and
model-based design of experiments.

The objective is to implement and execute the computational plan for a 3–2–1
adaptive pilot campaign. Use a hybrid optimizer: a global Particle Swarm
Optimization (PSO) layer followed by local refinement and verification with
IPOPT/Pyomo-DOE. Do not predefine the final thermal profiles: profile topology,
setpoint values and legal changepoints must be design variables. Existing A/B/C
profiles may be seeds and replay cases only.

## Working rules

1. Inspect `git status` first. Preserve all unrelated/user-owned changes,
   especially any modified Laboratory 2026 CSV. Never reset or overwrite them.
2. Treat everything under `fermentation_model/data/` as immutable.
3. Put campaign-specific adapters and workflows under `pilot_2026/`; put only
   genuinely reusable model/optimization components under `shared/`.
4. Write results only to campaign-owned result directories. Each run must have
   an immutable directory containing configuration, source hashes, code hash,
   Git commit, environment, random seeds, solver status and convergence data.
5. Use deterministic seeds and add automated tests. Do not claim runtime
   validation from syntax-only checks.
6. Do not use the existing integrated Pilot 2025 result as an official prior.
   It did not converge and its result directory contains artifacts from
   different executions. Regenerate any 2025 contribution coherently.
7. Do not push to GitHub until explicitly authorized. Make small, intentional
   local commits and report their hashes.

## Authoritative Pilot 2026 data contract

Re-run and verify:

```bash
python fermentation_model/pilot_2026/run_data_integration.py
python fermentation_model/pilot_2026/notebooks/build_and_execute_data_loading_qc_notebook.py
python -m unittest discover -s fermentation_model/tests -v
python fermentation_model/tools/campaign_audit.py --check-hashes
```

Read:

- `pilot_2026/data_integration_config.json`
- `pilot_2026/results/data_integration_2026/qc_summary.json`
- `pilot_2026/notebooks/pilot_2026_data_loading_qc.ipynb`
- `pilot_2026/results/data_integration_2026/primary_results_qc.csv`
- `temperature_controller_qc.csv`
- `co2_minute_qc.csv`
- `operational_events_qc.csv`
- `gc_results_long_qc.csv`
- `gc_mix_wine_pairs_long_qc.csv`

Enforce these confirmed rules:

- Nine fermentation windows run from the first to last primary sample.
- `Sonda1` is the executed temperature; setpoint is the commanded input.
- Protocol C is 16→18→21 °C, but an unexecuted stage after biological drying
  contributes no information.
- MassView is in Ln/min and must not be normalized a second time.
- Lot 1 CO2 (`26134–26136`) is QC-only and excluded from calibration.
- Lot 2–3 CO2 may be used with invalid/artificial-zero masks and a correlated
  error or effective-sample-size model; millions of seconds must not dominate
  the likelihood or FIM.
- There are 39 wine GC samples: six initial baselines and 33 paired 1:1 with
  33 MIX samples.
- MIX concentration is diluted-vial concentration ×1,000. Wine uses factor 1.
- Each MIX is an interval accumulation. Use captured mass
  `C_MIX × (V_A + V_B)` and its collection interval.
- NQ and below-LOQ are left-censored; never replace them with zero.
- Eight flagged `−0.3 % v/v` ethanol formula artifacts are processed as zero,
  while the original values remain auditable.
- Lot 3 operational events are owner-confirmed relative reconstructions from
  the matching Lot 1 runs.
- Second nutrient-pulse timestamps are density-triggered proxies unless a more
  exact record is found; represent this timing uncertainty explicitly.

## Phase A — model-ready adapter

Implement a tested Pilot 2026 adapter to the shared model contract. It must:

1. Load primary states, measured temperature, setpoint, usable CO2, operational
   events, wine aromas and interval condensate masses.
2. Keep units and observation operators explicit.
3. Separate active fermentation, cooling and postprocess.
4. Prevent double counting of total glucose+fructose and individual sugars.
5. Support censored GC likelihoods and interval-integrated aroma predictions.
6. Downsample CO2 reproducibly or model autocorrelation/effective sample size.
7. Carry lot, tank, campaign, yeast, storage, capture-system and sensor metadata.
8. Add a carbon-balance diagnostic and fail if an excluded run re-enters the
   calibration dataset.

Deliver a model-ready QC notebook and unit tests before fitting parameters.

## Phase B — regenerate the calibration prior

Build a hierarchical/multicampaign calibration:

- Pilot 2026 with X5 is the target dataset and dominates target kinetics.
- Laboratory 2026 may support mechanistic/kinetic identification after its own
  QC, with matrix and scale effects.
- Pilot 2025 may inform physical structure and weak priors only. Keep
  yeast-dependent kinetics campaign/yeast-specific unless pooling is justified
  quantitatively.
- Share physically justified thermodynamic or partition parameters.
- Treat tank, lot, scale, capture efficiency and sensor gain as nuisance/random
  effects where supported.

Fit in stages:

1. Primary fermentation states.
2. CO2 release/transfer using Pilot 2026 Lots 2–3.
3. Wine and interval-condensate aromas with censored likelihood.
4. Limited joint polish.

Use log-parameters, multiple starts/global basin exploration and IPOPT for local
fitting. Do not estimate all available parameters simultaneously. Use
estimability, profile likelihood and correlation diagnostics to choose a
reduced target set. Preserve the full posterior/parameter ensemble and
covariance, not only one local optimum. Use Schur complements or equivalent
marginalization for nuisance parameters in the design FIM.

Stop and report a failed calibration gate if multistarts do not converge, if
target profiles remain bound-dominated, or if structural residuals are large.

## Phase C — hybrid PSO–IPOPT design engine

Implement a finite but non-predefined design representation:

- Piecewise-constant temperature setpoint with an optimized number of active
  plateaus, optimized setpoint levels and optimized changepoint slots.
- Up to the operationally approved number of nutrient pulses, with optimized
  legal time and masses of organic product and DAP/FDA.
- Optional sampling/capture schedule only if the current analytical budget is
  explicitly configured.

Hard constraints:

- Temperature 15–25 °C.
- Minimum thermal-segment duration 12 h.
- Manual actions Monday–Friday, 09:00–17:00, America/Santiago. Until the owner
  approves otherwise, apply the same legal-window rule to setpoint changes.
- YAN conversion:
  `YAN mg/L = (100*m_org_g + 200*m_DAP_g)/V_L`.
- Put maximum dose per event, total YAN, number of pulses and latest allowable
  pulse in a reviewed configuration file. Do not invent missing safety bounds;
  fail closed before producing executable profiles if they are absent.
- Penalize or prohibit setpoint changes/pulses predicted to occur after drying.
- Use an identified setpoint→Sonda1 actuator model when predicting executed
  temperature; do not assume instantaneous perfect tracking.
- Add feasibility and minimum probability-of-completion constraints.

Global layer:

1. Implement or use a declared, tested PSO with Sobol-initialized particles,
   deterministic seeds, bounded velocities and explicit repair/projection for
   calendar and discrete topology constraints.
2. Seed part of the swarm with executed A/B/C profiles and the agreed anchor,
   but allow the swarm to leave those profiles.
3. Cache simulations by design hash and parallelize safely.
4. Run several independent seeds and report convergence/diversity, not only the
   best particle.

Local layer:

1. For each leading PSO topology/calendar solution, fix its discrete choices
   and refine continuous setpoints, doses and locally admissible times with
   IPOPT.
2. Validate the resulting design independently with Pyomo-DOE.
3. Compare finite-difference and Pyomo-DOE FIMs, scaling and positive
   semidefiniteness.
4. Reject candidates whose results depend strongly on finite-difference step,
   mesh or a single posterior sample.

Use a robust/pseudo-Bayesian objective over the posterior ensemble. Combine
expected D-optimality with an E-optimal or worst-quantile guardrail, then
subtract penalties for infeasibility, redundancy, sensor/capture risk and
failure to dry. Freeze the observation-error scale at the nominal prediction
for a standard FIM unless parameter-dependent variance is modeled formally.

## Phase D — qualify the engine

Before recommending experiments, pass:

- Replay of executed 2025/2026 profiles.
- Synthetic parameter-recovery tests.
- Leave-one-run/lot/tank-out predictive validation.
- PSO repeatability across seeds and population sizes.
- IPOPT multistart convergence.
- FIM agreement and finite-difference stability.
- Time-grid convergence.
- Calendar, temperature, minimum-duration and nutrition-unit tests.
- Predicted versus realized information for existing experiments.
- No silent fallback in the aroma thermodynamics model. Resolve the current
  conflicting trap-efficiency constants and make any fallback explicit.

## Phase E — adaptive campaign 3–2–1

### Wave 1: three fermentations

1. Anchor: operational wine-like reference, initially 18 °C isothermal with the
   approved standard nutrition.
2. Adaptive design 1: information-rich excitation selected by the hybrid
   optimizer.
3. Adaptive design 2: complementary excitation selected jointly with design 1,
   not simply the second-best individual candidate.

Optimize designs 2 and 3 as a pair. Rotate/randomize profiles among TK31, TK32
and TK33 to reduce tank/profile confounding. Record initial composition and
storage history for every reserved aliquot.

After Wave 1, refit with executed temperature/events and compare predicted with
realized information, posterior contraction, weakest eigenvector and residual
adequacy. If realized information is less than 50% of predicted or systematic
residuals exceed the preregistered threshold, change Wave 2 from parameter
estimation to model discrimination.

### Wave 2: two fermentations

Select the pair jointly using the updated posterior. One design should target
the weakest remaining parameter direction; the other should provide thermal,
nutritional, aroma or model-discrimination complementarity. Include a bridge
segment with the anchor to help distinguish wave/storage effects.

### Wave 3: one fermentation

Freeze model, posterior and policy before observing the run. Score prediction
interval coverage for primary states, CO2, drying time, wine aroma and
interval/cumulative captured aroma before refitting. This single run is an
out-of-sample trajectory confirmation, not replicated efficacy evidence. If
the Wave 2 gate fails, convert it explicitly into another adaptive
model-discrimination run instead of calling it confirmation.

## Required repository outputs

Create a clear campaign-owned structure, for example:

```text
pilot_2026/adaptive_design/
  build_model_dataset.py
  run_hierarchical_calibration.py
  design_space.py
  pso_optimizer.py
  run_adaptive_mbdoe.py
  run_realized_information.py
  design_constraints.json
  campaign_state.json
pilot_2026/results/adaptive_design_2026/
  baseline_calibration/
  engine_qualification/
  wave1_design/
  wave1_update/
  wave2_design/
  wave2_update/
  wave3_confirmation/
```

Update campaign/workflow registries and the repository map. Add notebooks or
reports that show parameter identifiability, posterior uncertainty, PSO
convergence, candidate diversity, feasibility, optimized profiles, expected
information and tank randomization.

## Final handoff

At each gate, report:

1. PASS, PASS CONDITIONAL or FAIL.
2. Exact commands executed and runtime status.
3. Data/code hashes and local commit.
4. Scientific findings and limitations.
5. Owner decisions still needed.
6. Whether profiles are exploratory or approved for physical execution.

Do not issue an executable fermentation schedule until all safety/operational
bounds are explicit, the calibration gate passes and the owner approves the
final candidate plots.
