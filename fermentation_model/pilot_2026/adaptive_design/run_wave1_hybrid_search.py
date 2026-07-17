from __future__ import annotations

import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.hybrid_optimizer import particle_swarm  # noqa: E402
from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    allowed_nutrition_times,
    anchor_policy,
    evaluate_campaign,
    load_json,
    policy_from_vector,
    prior_precision,
    representative_members,
    vector_bounds,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def _policy_vector(policy: DesignPolicy, config: dict) -> np.ndarray:
    slots = int(config["future_process"]["optimized_temperature_slots"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    allowed = allowed_nutrition_times(config)
    times = np.zeros(pulses, dtype=float)
    amounts = np.zeros(pulses, dtype=float)
    for index, (time, amount) in enumerate(policy.nutrition_mg_yan_l[:pulses]):
        times[index] = int(np.argmin(np.abs(allowed - float(time))))
        amounts[index] = float(amount)
    return np.concatenate([np.asarray(policy.temperature_c[:slots]), times, amounts])


def _seed_pairs(config: dict) -> np.ndarray:
    anchor = anchor_policy(config)
    slots = int(config["future_process"]["optimized_temperature_slots"])
    half = slots // 2
    policies = [
        anchor,
        DesignPolicy("warm_cool", tuple([23.0] * half + [16.0] * (slots - half)), ((24.0, 70.0), (72.0, 70.0))),
        DesignPolicy("cool_warm", tuple([16.0] * half + [23.0] * (slots - half)), ((48.0, 80.0),)),
        DesignPolicy("moderate_ramp", tuple(np.linspace(17.0, 22.0, slots)), ((48.0, 80.0),)),
    ]
    vectors = [_policy_vector(policy, config) for policy in policies]
    return np.vstack(
        [
            np.concatenate([vectors[0], vectors[0]]),
            np.concatenate([vectors[1], vectors[2]]),
            np.concatenate([vectors[2], vectors[1]]),
            np.concatenate([vectors[0], vectors[3]]),
        ]
    )


def _split_pair(values: np.ndarray, config: dict, prefix: str) -> tuple[DesignPolicy, DesignPolicy]:
    single = len(vector_bounds(config))
    return (
        policy_from_vector(f"{prefix}_A", values[:single], config),
        policy_from_vector(f"{prefix}_B", values[single:], config),
    )


def _rows_for_policy(policy: DesignPolicy) -> list[dict]:
    rows = [
        {
            "policy": policy.name,
            "action": "temperature_setpoint",
            "time_h": 12.0 * index,
            "value": float(value),
            "unit": "degC",
        }
        for index, value in enumerate(policy.temperature_c)
    ]
    rows.extend(
        {
            "policy": policy.name,
            "action": "nutrition",
            "time_h": float(time),
            "value": float(amount),
            "unit": "mgYAN/L",
        }
        for time, amount in policy.nutrition_mg_yan_l
    )
    return rows


def _ipopt_surrogate_refine(
    x: np.ndarray,
    bounds: np.ndarray,
    objective,
    config: dict,
) -> tuple[np.ndarray, dict]:
    executable = shutil.which("ipopt")
    try:
        import pyomo.environ as pyo
    except ImportError:
        return x, {"executed": False, "reason": "pyomo_not_available"}
    if not executable:
        return x, {"executed": False, "reason": "ipopt_not_on_path"}
    single = len(vector_bounds(config))
    slots = int(config["future_process"]["optimized_temperature_slots"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    local_indices = []
    steps = []
    for offset in (0, single):
        local_indices.extend(range(offset, offset + slots))
        steps.extend([float(config["search"]["local_temperature_step_c"])] * slots)
        amount_start = offset + slots + pulses
        local_indices.extend(range(amount_start, amount_start + pulses))
        steps.extend([float(config["search"]["local_yan_step_mg_l"])] * pulses)
    f0 = float(objective(x))
    gradient, curvature = [], []
    for index, step in zip(local_indices, steps):
        plus, minus = x.copy(), x.copy()
        plus[index] = min(plus[index] + step, bounds[index, 1])
        minus[index] = max(minus[index] - step, bounds[index, 0])
        hp = plus[index] - x[index]
        hm = x[index] - minus[index]
        fp, fm = float(objective(plus)), float(objective(minus))
        if hp <= 0.0 or hm <= 0.0:
            gradient.append(0.0)
            curvature.append(1e-3)
            continue
        gradient.append((fp - fm) / (hp + hm))
        symmetric_h = 0.5 * (hp + hm)
        curvature.append(max((fp - 2.0 * f0 + fm) / (symmetric_h**2), 1e-4))
    model = pyo.ConcreteModel()
    model.local_index = pyo.RangeSet(0, len(local_indices) - 1)
    trust = float(config["search"]["local_trust_fraction"])
    local_bounds = []
    for index in local_indices:
        radius = trust * (bounds[index, 1] - bounds[index, 0])
        local_bounds.append(
            (max(bounds[index, 0], x[index] - radius), min(bounds[index, 1], x[index] + radius))
        )
    model.z = pyo.Var(
        model.local_index,
        bounds=lambda _m, i: local_bounds[int(i)],
        initialize=lambda _m, i: float(x[local_indices[int(i)]]),
    )
    model.objective = pyo.Objective(
        expr=f0
        + sum(
            gradient[i] * (model.z[i] - float(x[local_indices[i]]))
            + 0.5 * curvature[i] * (model.z[i] - float(x[local_indices[i]])) ** 2
            for i in range(len(local_indices))
        )
    )
    result = pyo.SolverFactory("ipopt", executable=executable).solve(
        model, tee=False, options={"tol": 1e-8, "max_iter": 500}
    )
    candidate = x.copy()
    for local, index in enumerate(local_indices):
        candidate[index] = float(pyo.value(model.z[local]))
    actual = float(objective(candidate))
    accepted = actual <= f0
    return (
        candidate if accepted else x,
        {
            "executed": True,
            "interface": "pyomo.SolverFactory(ipopt)",
            "termination": str(result.solver.termination_condition),
            "surrogate_objective": float(pyo.value(model.objective)),
            "actual_objective_before": f0,
            "actual_objective_after": actual,
            "accepted_after_actual_revalidation": bool(accepted),
            "dimensions_refined": len(local_indices),
            "pulse_slot_indices_fixed": True,
        },
    )


def main() -> None:
    config = load_json(CONFIG_PATH)
    aroma_config = load_json(AROMA_CONFIG_PATH)
    ensemble_run = REPOSITORY_DIR / config["source_contract"]["joint_ensemble_run"]
    adapter_root = RESULT_ROOT / "wave1_mbdoe_adapter"
    adapter_runs = sorted(adapter_root.glob("*/adapter_gate.json"))
    if not adapter_runs or load_json(adapter_runs[-1])["verdict"] != "PASS":
        raise RuntimeError("Wave-1 MBDoE adapter gate is not PASS")
    ensemble = pd.read_csv(ensemble_run / "joint_parameter_ensemble.csv")
    partitions, partition_provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    all_representative = representative_members(ensemble, config)
    search_representative = all_representative[: int(config["search"]["representative_members"])]
    anchor = anchor_policy(config)
    single_bounds = vector_bounds(config)
    bounds = np.vstack([single_bounds, single_bounds])
    cache: dict[tuple, float] = {}

    def objective(values: np.ndarray) -> float:
        values = np.asarray(values, dtype=float)
        slots = int(config["future_process"]["optimized_temperature_slots"])
        pulses = int(config["nutrition"]["maximum_pulses"])
        rounded = values.copy()
        for offset in (0, len(single_bounds)):
            rounded[offset : offset + slots] = np.round(
                rounded[offset : offset + slots], int(config["search"]["cache_temperature_decimals"])
            )
            rounded[offset + slots : offset + slots + pulses] = np.round(
                rounded[offset + slots : offset + slots + pulses]
            )
            rounded[offset + slots + pulses : offset + len(single_bounds)] = np.round(
                rounded[offset + slots + pulses : offset + len(single_bounds)],
                int(config["search"]["cache_yan_decimals"]),
            )
        key = tuple(rounded.tolist())
        if key not in cache:
            pair = _split_pair(rounded, config, "candidate")
            score, _ = evaluate_campaign(
                (anchor, *pair), ensemble, search_representative, prior, config, partitions
            )
            cache[key] = -float(score)
        return cache[key]

    search = config["search"]
    pso = particle_swarm(
        objective,
        bounds,
        particles=int(search["particles"]),
        iterations=int(search["iterations"]),
        seed=int(config["seed"]),
        inertia=float(search["inertia"]),
        cognitive=float(search["cognitive"]),
        social=float(search["social"]),
        initial_positions=_seed_pairs(config),
    )
    refined_x, local = _ipopt_surrogate_refine(pso.x, bounds, objective, config)
    pair = _split_pair(refined_x, config, "wave1")
    robust_score, robust_evaluations = evaluate_campaign(
        (anchor, *pair), ensemble, list(range(len(ensemble))), prior, config, partitions
    )
    completion_probability = float(np.mean([row.completion for row in robust_evaluations]))
    robust_gains = np.asarray([row.information_gain for row in robust_evaluations])
    checks = {
        "adapter_gate_pass": True,
        "pso_objective_finite": math.isfinite(float(pso.fun)),
        "pso_improved_seeded_population": bool(pso.history[-1] < pso.history[0] - 1e-9),
        "full_ensemble_evaluated": len(robust_evaluations) == 64,
        "full_ensemble_completion_probability_at_least_95pct": completion_probability
        >= float(config["completion"]["minimum_probability"]),
        "all_information_gains_finite": bool(np.isfinite(robust_gains).all()),
        "ipopt_local_refinement_executed": bool(local["executed"]),
        "actual_objective_revalidated_after_local_step": bool(
            local.get("accepted_after_actual_revalidation", False)
        ),
    }
    computational_pass = all(
        value for name, value in checks.items() if name != "actual_objective_revalidated_after_local_step"
    )
    gate = {
        "gate": "phase_D_wave1_hybrid_search",
        "verdict": "PASS" if computational_pass else "PASS_CONDITIONAL" if checks["pso_objective_finite"] else "FAIL",
        "checks": checks,
        "pso": {
            "objective": float(pso.fun),
            "evaluations": int(pso.evaluations),
            "cache_entries": len(cache),
        },
        "ipopt_local": local,
        "full_ensemble": {
            "robust_information_score": robust_score,
            "median_information_gain": float(np.median(robust_gains)),
            "lower_decile_information_gain": float(np.quantile(robust_gains, 0.1)),
            "completion_probability": completion_probability,
            "maximum_residual_sugar_g_l": float(max(row.residual_sugar_g_l for row in robust_evaluations)),
        },
        "profiles_for_physical_execution": False,
        "next_action": "optimize ten-sample schedules and generate owner-review plots",
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_hybrid_search", config)
    history_path = run_dir / "pso_history.csv"
    policy_path = run_dir / "candidate_policy_actions.csv"
    scenarios_path = run_dir / "full_ensemble_validation.csv"
    gate_path = run_dir / "hybrid_search_gate.json"
    vector_path = run_dir / "candidate_vector.csv"
    config_path = run_dir / "wave1_mbdoe_config.json"
    pd.DataFrame({"iteration": np.arange(len(pso.history)), "best_objective": pso.history}).to_csv(history_path, index=False)
    pd.DataFrame(sum((_rows_for_policy(policy) for policy in (anchor, *pair)), [])).to_csv(policy_path, index=False)
    pd.DataFrame(
        {
            "ensemble_member": [row.member for row in robust_evaluations],
            "information_gain": [row.information_gain for row in robust_evaluations],
            "completion": [row.completion for row in robust_evaluations],
            "residual_sugar_g_l": [row.residual_sugar_g_l for row in robust_evaluations],
        }
    ).to_csv(scenarios_path, index=False)
    pd.DataFrame({"variable_index": np.arange(len(refined_x)), "value": refined_x}).to_csv(vector_path, index=False)
    write_json(gate_path, gate)
    write_json(config_path, config)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_hybrid_search",
        config=config,
        sources={
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "adapter_gate": adapter_runs[-1],
            "wave1_config": CONFIG_PATH,
        },
        code_paths=[Path(__file__), ADAPTIVE_DIR / "pilot_mbdoe_adapter.py", ADAPTIVE_DIR / "hybrid_optimizer.py"],
        random_seeds=[int(config["seed"])],
        status="completed" if gate["verdict"] == "PASS" else "conditional_or_failed",
        convergence={"pso": gate["pso"], "ipopt_local": local},
        gate=gate,
        outputs=[history_path, policy_path, scenarios_path, gate_path, vector_path, config_path],
    )
    write_json(run_dir / "partition_surrogate_provenance.json", partition_provenance)
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if gate["verdict"] == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
