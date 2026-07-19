from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.final_search_logic import policy_distance  # noqa: E402
from pilot_2026.adaptive_design.hybrid_optimizer import (  # noqa: E402
    CheckpointedSwarmState,
    advance_checkpointed_swarm,
    initialize_checkpointed_swarm,
    rescore_checkpointed_swarm,
)
from pilot_2026.adaptive_design.operational_coverage import (  # noqa: E402
    AUTHORIZATION_FLAGS,
    REVIEW_SUFFIX,
    WATERMARK,
    archetype_checks,
    closed_authorization,
    contrast_checks,
    coverage_policy_from_unit_vector,
    coverage_policy_metrics,
    derive_robust_setpoint_envelope,
    explicit_controller_blocks,
    gate_verdict,
    lexicographic_select,
    nutrition_margin_row,
    partially_harmonize_schedules,
    retained_information_improvement_fraction,
    robust_initial_jump_rows,
    validate_physical_temperature,
)
from pilot_2026.adaptive_design.optimize_wave1_sampling_and_plots import (  # noqa: E402
    _capture_intervals,
    _load_policies,
    _operational_conflicts,
    _optimize_schedules,
    _score,
)
from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
    simulate_aroma,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    SPECIES,
    _aroma_log_values,
    _forcing_from_core,
    _future_design,
    _theta,
    actuator_temperature_trajectory,
    allowed_sampling_times,
    anchor_policy,
    evaluate_campaign,
    load_json,
    load_wave1_config,
    logdet,
    parameter_columns,
    prepare_design,
    prepare_sampling_sensitivity_cache,
    prior_precision,
    representative_members,
    robust_information_metrics,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    capture_git_state,
    create_immutable_run_directory,
    filesystem_path,
    relative_or_absolute,
    sha256_file,
    sha256_payload,
    verify_manifest_output,
    write_json,
)


COVERAGE_CONFIG_PATH = ADAPTIVE_DIR / "coverage_design_config.json"
MODEL_CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
FIGURE_NAMES = (
    "operational_coverage_map.png",
    "coverage_strategy_profiles.png",
    "coverage_vs_information_tradeoff.png",
    "parameter_variance_comparison.png",
    "cold_early_predictions.png",
    "cold_late_predictions.png",
    "warm_early_predictions.png",
    "warm_late_predictions.png",
    "physical_temperature_envelope.png",
    "wave_3_2_1_coverage_plan.png",
)
PALETTE = {
    "blue": "#235789",
    "gold": "#D4A72C",
    "orange": "#E07A3F",
    "olive": "#7A8B3A",
    "pink": "#B85C8A",
    "ink": "#1F2937",
    "grey": "#9CA3AF",
}


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(filesystem_path(path).read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any], outputs: list[Path]) -> None:
    enriched = dict(payload)
    enriched.setdefault("watermark", WATERMARK)
    for key, value in closed_authorization().items():
        enriched.setdefault(key, value)
    write_json(path, enriched)
    outputs.append(path)


def _write_csv(path: Path, frame: pd.DataFrame, outputs: list[Path]) -> None:
    output = frame.copy()
    if "watermark" not in output:
        output["watermark"] = WATERMARK
    if "candidate_only_not_for_physical_execution" not in output:
        output["candidate_only_not_for_physical_execution"] = True
    output.to_csv(filesystem_path(path), index=False, lineterminator="\n")
    outputs.append(path)


def _gzip_checkpoint(path: Path, payload: dict[str, Any]) -> None:
    enriched = {
        **payload,
        "watermark": WATERMARK,
        **closed_authorization(),
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(filesystem_path(temporary), "wt", encoding="utf-8") as handle:
        json.dump(enriched, handle, ensure_ascii=False, sort_keys=True)
    filesystem_path(temporary).replace(filesystem_path(path))


def _source_contract(
    coverage_config: dict[str, Any]
) -> tuple[dict[str, Path], dict[str, Any]]:
    runs = {
        name: _run_path(value) for name, value in coverage_config["source_runs"].items()
    }
    verification = {
        "adapter_gate": verify_manifest_output(runs["adapter"], "adapter_gate.json"),
        "ensemble": verify_manifest_output(
            runs["joint_ensemble"], "joint_parameter_ensemble.csv"
        ),
        "aroma_gate": verify_manifest_output(
            runs["aroma_calibration"], "aroma_calibration_gate.json"
        ),
        "actuator": verify_manifest_output(
            runs["temperature_actuator"], "actuator_estimates_by_run.csv"
        ),
        "reference_search_gate": verify_manifest_output(
            runs["information_reference"], "final_search_gate.json"
        ),
        "reference_actions": verify_manifest_output(
            runs["information_reference"], "final_candidate_policy_actions.csv"
        ),
        "reference_full_ensemble": verify_manifest_output(
            runs["information_reference"], "full_ensemble_validation.csv"
        ),
        "reference_actuator": verify_manifest_output(
            runs["information_reference"], "actuator_robustness_validation.csv"
        ),
        "reference_sampling_gate": verify_manifest_output(
            runs["sampling_reference"], "final_sampling_gate.json"
        ),
        "reference_sampling_schedule": verify_manifest_output(
            runs["sampling_reference"], "optimized_sampling_schedule.csv"
        ),
        "reference_sampling_comparison": verify_manifest_output(
            runs["sampling_reference"], "sampling_candidate_comparison.csv"
        ),
    }
    return runs, verification


def _reference_schedule(path: Path) -> dict[str, tuple[float, ...]]:
    frame = pd.read_csv(filesystem_path(path))
    return {
        str(name): tuple(sorted(group["time_h"].astype(float).tolist()))
        for name, group in frame.groupby("policy", sort=False)
    }


def _anchor_with_dose(
    model_config: dict[str, Any], yan_mg_l: float, name: str
) -> DesignPolicy:
    base = anchor_policy(model_config)
    return DesignPolicy(
        name,
        base.temperature_c,
        ((46.0, float(yan_mg_l)),),
    )


def _decode_strategy_pair(
    values: np.ndarray,
    strategy: str,
    dose_design: str,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    envelope: dict[str, Any],
) -> tuple[DesignPolicy, DesignPolicy, DesignPolicy]:
    dose = coverage_config["dose_designs"][dose_design]
    anchor = _anchor_with_dose(
        model_config,
        float(dose["anchor_yan_mg_l"]),
        f"{strategy}__{dose_design}__anchor",
    )
    if strategy == "II_diagonal_1":
        cells = (("cold", "early"), ("warm", "late"))
    elif strategy == "III_diagonal_2":
        cells = (("cold", "late"), ("warm", "early"))
    else:
        raise ValueError(f"Unsupported coverage search strategy: {strategy}")
    vector = np.asarray(values, dtype=float)
    if vector.shape != (24,):
        raise ValueError("A strategy vector must contain two 12-coordinate policies")
    policies = []
    for index, ((thermal, nutrition), subvector) in enumerate(
        zip(cells, (vector[:12], vector[12:])), start=1
    ):
        policies.append(
            coverage_policy_from_unit_vector(
                f"{strategy}__{dose_design}__{thermal}_{nutrition}",
                thermal,
                nutrition,
                float(dose["coverage_yan_mg_l"]),
                subvector,
                model_config,
                coverage_config,
                envelope,
            )
        )
    return (anchor, policies[0], policies[1])


def _strategy_contrast(
    policies: tuple[DesignPolicy, DesignPolicy, DesignPolicy],
    strategy: str,
    coverage_config: dict[str, Any],
) -> dict[str, Any]:
    _, first, second = policies
    if strategy == "II_diagonal_1":
        cold, warm, early, late = first, second, first, second
    else:
        cold, warm, early, late = first, second, second, first
    return contrast_checks(cold, warm, early, late, coverage_config)


def _candidate_identifier(policies: tuple[DesignPolicy, ...]) -> str:
    payload = [
        {
            "temperature_c": list(policy.temperature_c),
            "nutrition_mg_yan_l": [list(row) for row in policy.nutrition_mg_yan_l],
        }
        for policy in policies
    ]
    return sha256_payload(payload)[:16]


def _campaign_metrics(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    members: list[int],
    prior: np.ndarray,
    model_config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
) -> tuple[float, list[Any], dict[str, float]]:
    score, evaluations = evaluate_campaign(
        policies,
        ensemble,
        members,
        prior,
        model_config,
        partitions,
        design_cache=design_cache,
    )
    gains = np.asarray([row.information_gain for row in evaluations], dtype=float)
    posterior = [prior + row.fim for row in evaluations]
    metrics = robust_information_metrics(gains, posterior, None, model_config)
    metrics["penalized_search_score"] = float(score)
    metrics["completion_probability"] = float(
        np.mean([row.completion for row in evaluations])
    )
    metrics["maximum_residual_sugar_g_l"] = float(
        max(row.residual_sugar_g_l for row in evaluations)
    )
    metrics["minimum_drying_margin_h"] = float(
        min(row.drying_time_h - row.latest_action_time_h for row in evaluations)
    )
    return float(score), evaluations, metrics


def _save_swarm_checkpoint(
    checkpoint_dir: Path,
    strategy: str,
    dose_design: str,
    state: CheckpointedSwarmState,
    sequence: int,
    model_config_hash: str,
    coverage_config_hash: str,
) -> Path:
    path = checkpoint_dir / (
        f"{strategy}__{dose_design}__seed_{state.seed}__{sequence:03d}__"
        f"{state.fidelity}__iteration_{state.iteration:03d}.json.gz"
    )
    _gzip_checkpoint(
        path,
        {
            "checkpoint_schema_version": 1,
            "strategy": strategy,
            "dose_design": dose_design,
            "checkpoint_sequence": sequence,
            "swarm": state.to_payload(),
            "model_config_sha256": model_config_hash,
            "coverage_config_sha256": coverage_config_hash,
        },
    )
    return path


def _policy_signature(policy: DesignPolicy) -> tuple[Any, ...]:
    return (tuple(policy.temperature_c), tuple(policy.nutrition_mg_yan_l))


def _physical_summary_cached(
    policy: DesignPolicy,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    key = _policy_signature(policy)
    if key not in cache:
        cache[key] = validate_physical_temperature(
            policy, model_config, coverage_config
        )
    return cache[key]


def _full_candidate_record(
    *,
    vector: np.ndarray,
    policies: tuple[DesignPolicy, DesignPolicy, DesignPolicy],
    strategy: str,
    dose_design: str,
    seed: int,
    source_type: str,
    parent_candidate_id: str | None,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
) -> dict[str, Any]:
    score, evaluations, information = _campaign_metrics(
        policies,
        ensemble,
        list(range(len(ensemble))),
        prior,
        model_config,
        partitions,
        design_cache,
    )
    expected = (
        (("cold", "early"), ("warm", "late"))
        if strategy == "II_diagonal_1"
        else (("cold", "late"), ("warm", "early"))
    )
    archetype = [
        archetype_checks(policy, thermal, nutrition, model_config, coverage_config)
        for policy, (thermal, nutrition) in zip(policies[1:], expected)
    ]
    contrast = _strategy_contrast(policies, strategy, coverage_config)
    physical_summaries = [
        _physical_summary_cached(policy, model_config, coverage_config, physical_cache)[
            1
        ]
        for policy in policies
    ]
    policy_metrics = [
        coverage_policy_metrics(policy, model_config, coverage_config)
        for policy in policies
    ]
    nutrition_rows = [
        nutrition_margin_row(
            strategy,
            dose_design,
            policy,
            model_config,
            coverage_config,
        )
        for policy in policies
    ]
    controller = [
        explicit_controller_blocks(policy, coverage_config) for policy in policies
    ]
    controller_pass = all(
        len(frame) == 42
        and math.isclose(float(frame["end_h"].max()), 504.0, abs_tol=1e-9)
        and bool(frame.iloc[-1]["hold_last_setpoint"])
        for frame in controller
    )
    coverage_pass = bool(
        all(all(check.values()) for check in archetype)
        and contrast["thermal_contrast_at_least_4c"] is True
        and contrast["nutrition_contrast_at_least_36h"] is True
        and contrast["coverage_treatments_not_equivalent"] is True
    )
    thermal_pass = bool(
        controller_pass
        and all(row["initial_jump_robust_pass"] for row in physical_summaries)
        and all(row["physical_temperature_robust_pass"] for row in physical_summaries)
        and all(
            metric.maximum_internal_jump_c <= 5.0 + 1e-9 for metric in policy_metrics
        )
    )
    nutrition_pass = bool(all(row["nutrition_feasible"] for row in nutrition_rows))
    completion_pass = bool(
        information["completion_probability"]
        >= float(model_config["completion"]["minimum_probability"]) - 1e-12
    )
    minimum_nutrition_margin = min(
        float(row["organic_margin_fraction"]) for row in nutrition_rows
    )
    candidate_id = _candidate_identifier(policies)
    record = {
        "candidate_id": candidate_id,
        "parent_candidate_id": parent_candidate_id,
        "strategy": strategy,
        "dose_design": dose_design,
        "independent_seed": int(seed),
        "source_type": source_type,
        "robust_score": float(information["robust_score"]),
        "penalized_search_score": float(score),
        "median": float(information["median"]),
        "q10": float(information["q10"]),
        "minimum": float(information["minimum"]),
        "tail_cvar": float(information["tail_cvar"]),
        "completion_probability": float(information["completion_probability"]),
        "minimum_drying_margin_h": float(information["minimum_drying_margin_h"]),
        "maximum_residual_sugar_g_l": float(information["maximum_residual_sugar_g_l"]),
        "minimum_posterior_fim_eigenvalue": float(
            information["minimum_posterior_fim_eigenvalue"]
        ),
        "maximum_posterior_fim_condition_number": float(
            information["maximum_posterior_fim_condition_number"]
        ),
        "maximum_internal_jump_c": max(
            metric.maximum_internal_jump_c for metric in policy_metrics
        ),
        "maximum_initial_jump_c": max(
            metric.maximum_initial_jump_c for metric in policy_metrics
        ),
        "maximum_initial_or_internal_jump_c": max(
            max(metric.maximum_internal_jump_c, metric.maximum_initial_jump_c)
            for metric in policy_metrics
        ),
        "minimum_robust_physical_temperature_c": min(
            float(row["minimum_physical_temperature_c"]) for row in physical_summaries
        ),
        "maximum_robust_physical_temperature_c": max(
            float(row["maximum_physical_temperature_c"]) for row in physical_summaries
        ),
        "total_thermal_variation_c": sum(
            metric.total_thermal_variation_c for metric in policy_metrics
        ),
        "temperature_changes": sum(
            metric.temperature_changes for metric in policy_metrics
        ),
        "nutrition_operational_margin_fraction": minimum_nutrition_margin,
        "coverage_pass": coverage_pass,
        "thermal_pass": thermal_pass,
        "nutrition_pass": nutrition_pass,
        "completion_pass": completion_pass,
        "controller_explicit_to_504h": controller_pass,
        "first_24h_mean_contrast_c": float(contrast["first_24h_mean_contrast_c"]),
        "nutrition_time_contrast_h": float(contrast["nutrition_time_contrast_h"]),
        "feasible": bool(
            coverage_pass and thermal_pass and nutrition_pass and completion_pass
        ),
        "vector": np.asarray(vector, dtype=float),
        "policies": policies,
        "evaluations": evaluations,
        "archetype_checks": archetype,
        "contrast_checks": contrast,
    }
    return record


def _manual_seed_positions() -> np.ndarray:
    seeds = []
    for first, second, active in ((0.35, 0.45, 0.20), (0.65, 0.60, 0.72)):
        pair = np.full(24, 0.5, dtype=float)
        for offset in (0, 12):
            pair[offset] = first
            pair[offset + 1] = second
            pair[offset + 2] = active
            pair[offset + 5] = active
            pair[offset + 8] = active
            pair[offset + 3] = 0.25
            pair[offset + 6] = 0.55
            pair[offset + 9] = 0.85
            pair[offset + 4] = 0.30
            pair[offset + 7] = 0.70
            pair[offset + 10] = 0.35
            pair[offset + 11] = 0.50
        seeds.append(pair)
    drying_safe = np.full(24, 0.2, dtype=float)
    # Cold corner: satisfy the first-24 h archetype, then warm in two legal jumps.
    drying_safe[0:12] = [
        0.50,
        0.50,
        1.00,
        0.00,
        1.00,
        1.00,
        0.18,
        0.70,
        0.00,
        0.50,
        0.50,
        0.20,
    ]
    # Warm corner: use the robust ramp in block 2 and hold a high safe setpoint.
    drying_safe[12:24] = [
        0.50,
        0.80,
        0.20,
        0.25,
        0.50,
        0.20,
        0.55,
        0.50,
        0.20,
        0.85,
        0.50,
        0.50,
    ]
    seeds.append(drying_safe)
    return np.asarray(seeds, dtype=float)


def _practical_convergence(
    seed_champions: list[dict[str, Any]],
    continuation_rows: list[dict[str, Any]],
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
) -> dict[str, Any]:
    best = max(float(row["robust_score"]) for row in seed_champions)
    tolerance = float(coverage_config["search"]["practical_score_tolerance_fraction"])
    near = [
        row
        for row in seed_champions
        if best - float(row["robust_score"])
        <= tolerance * max(abs(best), 1e-12) + 1e-12
    ]
    distance_rows = []
    for left_index, left in enumerate(near):
        for right in near[left_index + 1 :]:
            distance = policy_distance(
                tuple(left["policies"][1:]),
                tuple(right["policies"][1:]),
                model_config,
            )
            distance_rows.append(
                {
                    "left_seed": int(left["independent_seed"]),
                    "right_seed": int(right["independent_seed"]),
                    **distance,
                }
            )
    close_seeds: set[int] = set()
    for row in distance_rows:
        if float(row["policy_distance"]) <= 0.35:
            close_seeds.update((int(row["left_seed"]), int(row["right_seed"])))
    maximum_continuation_fraction = max(
        (float(row["continuation_improvement_fraction"]) for row in continuation_rows),
        default=math.inf,
    )
    checks = {
        "five_seed_champions_available": len(seed_champions)
        == len(coverage_config["search"]["independent_seeds"]),
        "at_least_three_seed_champions_within_0_5pct": len(near) >= 3,
        "all_near_best_champions_feasible": all(row["feasible"] for row in near),
        "operational_policy_family_has_three_seeds": len(close_seeds) >= 3,
        "eight_member_continuation_completed": len(continuation_rows)
        == len(coverage_config["search"]["independent_seeds"]),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "best_full_ensemble_robust_score": best,
        "near_best_seeds": [int(row["independent_seed"]) for row in near],
        "score_tolerance_fraction": tolerance,
        "maximum_eight_member_continuation_improvement_fraction": maximum_continuation_fraction,
        "pairwise_policy_distances": distance_rows,
    }


def _run_strategy_variant_search(
    *,
    strategy: str,
    dose_design: str,
    checkpoint_dir: Path,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    envelope: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
) -> dict[str, Any]:
    representatives = representative_members(ensemble, model_config)
    four_members = representatives[:4]
    eight_members = representatives[:8]
    search = coverage_config["search"]
    objective_cache: dict[tuple[str, str], float] = {}

    def objective(values: np.ndarray, members: list[int], fidelity: str) -> float:
        policies = _decode_strategy_pair(
            values,
            strategy,
            dose_design,
            model_config,
            coverage_config,
            envelope,
        )
        key = (fidelity, _candidate_identifier(policies))
        if key not in objective_cache:
            score, _ = evaluate_campaign(
                policies,
                ensemble,
                members,
                prior,
                model_config,
                partitions,
                design_cache=design_cache,
            )
            objective_cache[key] = -float(score)
        return objective_cache[key]

    four = lambda values: objective(values, four_members, "four_member")
    eight = lambda values: objective(values, eight_members, "eight_member")
    bounds = np.tile(np.asarray([[0.0, 1.0]], dtype=float), (24, 1))
    all_records: list[dict[str, Any]] = []
    seed_champions: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []
    local_rows: list[dict[str, Any]] = []
    continuation_rows: list[dict[str, Any]] = []
    checkpoint_paths: list[Path] = []
    for seed in (int(value) for value in search["independent_seeds"]):
        sequence = 1
        state = initialize_checkpointed_swarm(
            four,
            bounds,
            particles=int(search["particles"]),
            seed=seed,
            fidelity="four_member",
            initial_positions=_manual_seed_positions(),
        )
        checkpoint_paths.append(
            _save_swarm_checkpoint(
                checkpoint_dir,
                strategy,
                dose_design,
                state,
                sequence,
                sha256_file(MODEL_CONFIG_PATH),
                sha256_file(COVERAGE_CONFIG_PATH),
            )
        )
        history_rows.append(
            {
                "strategy": strategy,
                "dose_design": dose_design,
                "independent_seed": seed,
                "fidelity": state.fidelity,
                "iteration": state.iteration,
                "best_score": -state.global_best_value,
            }
        )
        for _ in range(int(search["stage1_four_member_iterations"])):
            state = advance_checkpointed_swarm(
                state,
                four,
                bounds,
                inertia=float(search["inertia"]),
                cognitive=float(search["cognitive"]),
                social=float(search["social"]),
            )
            sequence += 1
            checkpoint_paths.append(
                _save_swarm_checkpoint(
                    checkpoint_dir,
                    strategy,
                    dose_design,
                    state,
                    sequence,
                    sha256_file(MODEL_CONFIG_PATH),
                    sha256_file(COVERAGE_CONFIG_PATH),
                )
            )
            history_rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "independent_seed": seed,
                    "fidelity": state.fidelity,
                    "iteration": state.iteration,
                    "best_score": -state.global_best_value,
                }
            )
        stage1_score = -float(state.global_best_value)
        state = rescore_checkpointed_swarm(state, eight, fidelity="eight_member")
        sequence += 1
        checkpoint_paths.append(
            _save_swarm_checkpoint(
                checkpoint_dir,
                strategy,
                dose_design,
                state,
                sequence,
                sha256_file(MODEL_CONFIG_PATH),
                sha256_file(COVERAGE_CONFIG_PATH),
            )
        )
        stage2_initial = -float(state.global_best_value)
        for _ in range(int(search["stage2_eight_member_iterations"])):
            state = advance_checkpointed_swarm(
                state,
                eight,
                bounds,
                inertia=float(search["inertia"]),
                cognitive=float(search["cognitive"]),
                social=float(search["social"]),
            )
            sequence += 1
            checkpoint_paths.append(
                _save_swarm_checkpoint(
                    checkpoint_dir,
                    strategy,
                    dose_design,
                    state,
                    sequence,
                    sha256_file(MODEL_CONFIG_PATH),
                    sha256_file(COVERAGE_CONFIG_PATH),
                )
            )
            history_rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "independent_seed": seed,
                    "fidelity": state.fidelity,
                    "iteration": state.iteration,
                    "best_score": -state.global_best_value,
                }
            )
        stage2_final = -float(state.global_best_value)
        continuation_rows.append(
            {
                "strategy": strategy,
                "dose_design": dose_design,
                "independent_seed": seed,
                "four_member_final_score": stage1_score,
                "eight_member_initial_score": stage2_initial,
                "eight_member_final_score": stage2_final,
                "continuation_improvement_fraction": max(
                    stage2_final - stage2_initial, 0.0
                )
                / max(abs(stage2_initial), 1e-12),
            }
        )
        original_vector = state.global_best_position.copy()
        original_policies = _decode_strategy_pair(
            original_vector,
            strategy,
            dose_design,
            model_config,
            coverage_config,
            envelope,
        )
        original = _full_candidate_record(
            vector=original_vector,
            policies=original_policies,
            strategy=strategy,
            dose_design=dose_design,
            seed=seed,
            source_type="seed_champion_original",
            parent_candidate_id=None,
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            partitions=partitions,
            design_cache=design_cache,
            physical_cache=physical_cache,
        )
        all_records.append(original)
        best = original
        rng = np.random.default_rng(seed + 65537)
        trial_count = int(search["local_full_ensemble_trials_per_seed"])
        continuous = np.asarray([0, 1, 4, 7, 10, 12, 13, 16, 19, 22], dtype=int)
        for trial in range(1, trial_count + 1):
            proposal_vector = original_vector.copy()
            coordinate = int(rng.choice(continuous))
            direction = -1.0 if trial % 2 else 1.0
            proposal_vector[coordinate] = np.clip(
                proposal_vector[coordinate]
                + direction * float(search["local_step_fraction"]),
                0.0,
                1.0,
            )
            proposal_policies = _decode_strategy_pair(
                proposal_vector,
                strategy,
                dose_design,
                model_config,
                coverage_config,
                envelope,
            )
            proposal = _full_candidate_record(
                vector=proposal_vector,
                policies=proposal_policies,
                strategy=strategy,
                dose_design=dose_design,
                seed=seed,
                source_type="full_64_local_refinement_trial",
                parent_candidate_id=original["candidate_id"],
                ensemble=ensemble,
                prior=prior,
                model_config=model_config,
                coverage_config=coverage_config,
                partitions=partitions,
                design_cache=design_cache,
                physical_cache=physical_cache,
            )
            accepted = bool(
                proposal["feasible"]
                and proposal["penalized_search_score"]
                > best["penalized_search_score"] + 1e-9
            )
            proposal["local_refinement_accepted"] = accepted
            local_rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "independent_seed": seed,
                    "trial": trial,
                    "coordinate": coordinate,
                    "parent_candidate_id": original["candidate_id"],
                    "proposal_candidate_id": proposal["candidate_id"],
                    "parent_robust_score": original["robust_score"],
                    "proposal_robust_score": proposal["robust_score"],
                    "parent_penalized_search_score": original["penalized_search_score"],
                    "proposal_penalized_search_score": proposal[
                        "penalized_search_score"
                    ],
                    "accepted_on_real_64_member_objective": accepted,
                }
            )
            all_records.append(proposal)
            if accepted:
                best = proposal
        if best is original:
            best = dict(original)
            best["local_refinement_qualified"] = True
            best["local_refinement_qualification_reason"] = (
                "no_better_feasible_step_on_real_64_member_objective"
            )
        else:
            best["local_refinement_qualified"] = True
            best["local_refinement_qualification_reason"] = (
                "accepted_improvement_on_real_64_member_objective"
            )
        seed_champions.append(best)
    eligible = [
        row
        for row in all_records
        if row["source_type"] == "seed_champion_original"
        or row.get("local_refinement_accepted") is True
    ]
    selected, selection = lexicographic_select(
        eligible,
        float(coverage_config["objective"]["near_best_fraction"]),
    )
    convergence = _practical_convergence(
        seed_champions,
        continuation_rows,
        model_config,
        coverage_config,
    )
    return {
        "strategy": strategy,
        "dose_design": dose_design,
        "selected": selected,
        "selection": selection,
        "convergence": convergence,
        "all_records": all_records,
        "seed_champions": seed_champions,
        "history_rows": history_rows,
        "local_rows": local_rows,
        "continuation_rows": continuation_rows,
        "checkpoint_paths": checkpoint_paths,
        "objective_cache_entries": len(objective_cache),
    }


def _prepare_sampling_context(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    model_config: dict[str, Any],
    partitions: dict[str, Any],
    prepared_cache: dict[tuple[Any, ...], Any],
) -> tuple[dict[tuple[int, str], Any], dict[str, np.ndarray], pd.DataFrame]:
    sampling_cache: dict[tuple[int, str], Any] = {}
    drying_rows = []
    for member in range(len(ensemble)):
        for policy in policies:
            key = (member, *_policy_signature(policy))
            if key not in prepared_cache:
                prepared = prepare_design(
                    policy,
                    ensemble.iloc[member],
                    model_config,
                    partitions,
                )
                if prepared is None:
                    raise RuntimeError(
                        f"Sampling preparation failed: member={member}, policy={policy.name}"
                    )
                prepared_cache[key] = (
                    prepared,
                    prepare_sampling_sensitivity_cache(prepared, model_config),
                )
            prepared, sensitivity = prepared_cache[key]
            sampling_cache[(member, policy.name)] = sensitivity
            drying_rows.append(
                {
                    "ensemble_member": member,
                    "policy": policy.name,
                    "drying_time_h": float(prepared.drying_time_h),
                }
            )
    drying = pd.DataFrame(drying_rows)
    allowed = allowed_sampling_times(model_config)
    valid_by_policy: dict[str, np.ndarray] = {}
    for policy in policies:
        values = drying.loc[drying.policy.eq(policy.name), "drying_time_h"].to_numpy(
            dtype=float
        )
        valid = np.asarray(
            [
                float(time_h)
                for time_h in allowed
                if float(np.mean(values >= float(time_h)))
                >= float(model_config["completion"]["minimum_probability"])
            ],
            dtype=float,
        )
        if len(valid) < int(model_config["sampling"]["samples_per_process"]):
            raise RuntimeError(
                f"Not enough robust legal sample slots for {policy.name}"
            )
        if not math.isclose(float(valid[0]), 0.0, abs_tol=1e-12):
            raise RuntimeError("Every sampling space must retain the basal t=0 slot")
        valid_by_policy[policy.name] = valid
    return sampling_cache, valid_by_policy, drying


def _schedule_rows(
    schedules: dict[str, tuple[float, ...]],
    strategy: str,
    dose_design: str,
    mode: str,
    model_config: dict[str, Any],
) -> pd.DataFrame:
    start = datetime.fromisoformat(model_config["future_process"]["start_local"])
    rows = []
    for policy, times in schedules.items():
        for number, time_h in enumerate(times, start=1):
            timestamp = start + timedelta(hours=float(time_h))
            rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "sampling_mode": mode,
                    "policy": policy,
                    "sample_number": number,
                    "time_h": float(time_h),
                    "local_timestamp": timestamp.isoformat(),
                    "weekday": timestamp.strftime("%A"),
                    "sample_before_action": True,
                }
            )
    return pd.DataFrame(rows)


def _sampling_mode_result(
    *,
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    strategy: str,
    dose_design: str,
    mode: str,
    sampling_cache: dict[tuple[int, str], Any],
    prior: np.ndarray,
    model_config: dict[str, Any],
    fim_cache: dict[tuple, np.ndarray],
    current_reference_gains: np.ndarray,
    search_record: dict[str, Any],
    harmonization: dict[str, Any] | None,
) -> dict[str, Any]:
    members = list(range(len(current_reference_gains)))
    score, gains, posterior = _score(
        schedules,
        policies,
        sampling_cache,
        members,
        prior,
        model_config,
        fim_cache,
    )
    metrics = robust_information_metrics(
        gains,
        posterior,
        current_reference_gains,
        model_config,
    )
    conflicts = _operational_conflicts(schedules, policies, model_config)
    captures = _capture_intervals(schedules, model_config, sampling_cache, members)
    schedule_frame = _schedule_rows(
        schedules, strategy, dose_design, mode, model_config
    )
    legal = bool(
        schedule_frame["weekday"]
        .isin(["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"])
        .all()
        and schedule_frame["local_timestamp"]
        .map(
            lambda value: 9 <= datetime.fromisoformat(value).hour <= 17
            or math.isclose(
                (
                    datetime.fromisoformat(value)
                    - datetime.fromisoformat(
                        model_config["future_process"]["start_local"]
                    )
                ).total_seconds()
                / 3600.0,
                0.0,
                abs_tol=1e-12,
            )
        )
        .all()
    )
    unresolved = bool(len(conflicts) and conflicts["status"].eq("unresolved").any())
    capture_pass = bool(
        len(captures)
        and captures.groupby(["policy", "species"])["capture_interval"]
        .nunique()
        .eq(9)
        .all()
        and captures["mass_conservation_max_abs_error_ug"].le(1e-8).all()
    )
    harmonized_pass = True
    if mode == "partially_harmonized":
        assert harmonization is not None
        harmonized_pass = bool(
            harmonization["score_loss_fraction"]
            <= float(model_config.get("coverage_sampling_loss_fraction", 0.005)) + 1e-12
            and harmonization["operator_round_reduction"]
            >= int(
                model_config.get(
                    "coverage_harmonized_minimum_operator_round_reduction", 1
                )
            )
        )
    checks = {
        "64_members_evaluated": len(gains) == 64,
        "ten_samples_per_process": all(
            len(values) == 10 for values in schedules.values()
        ),
        "basal_t0_included": all(
            math.isclose(float(values[0]), 0.0, abs_tol=1e-12)
            for values in schedules.values()
        ),
        "nine_capture_intervals": capture_pass,
        "manual_windows_legal": legal,
        "manual_capacity_three": not unresolved,
        "sample_before_action": model_config["operations"]["sample_event_order"]
        == "sample_before_action",
        "harmonization_guardrail": harmonized_pass,
    }
    row = {
        "strategy": strategy,
        "dose_design": dose_design,
        "sampling_mode": mode,
        "robust_score": float(metrics["robust_score"]),
        "median": float(metrics["median"]),
        "q10": float(metrics["q10"]),
        "minimum": float(metrics["minimum"]),
        "tail_cvar": float(metrics["tail_cvar"]),
        "fraction_exceeding_information_reference": float(
            metrics["fraction_exceeding_reference"]
        ),
        "worst_paired_loss_vs_information_reference": float(
            metrics["worst_paired_loss_vs_reference"]
        ),
        "completion_probability": float(search_record["completion_probability"]),
        "minimum_drying_margin_h": float(search_record["minimum_drying_margin_h"]),
        "maximum_internal_jump_c": float(search_record["maximum_internal_jump_c"]),
        "maximum_initial_jump_c": float(search_record["maximum_initial_jump_c"]),
        "minimum_robust_physical_temperature_c": float(
            search_record["minimum_robust_physical_temperature_c"]
        ),
        "maximum_robust_physical_temperature_c": float(
            search_record["maximum_robust_physical_temperature_c"]
        ),
        "anchor_yan_mg_l": float(policies[0].nutrition_mg_yan_l[0][1]),
        "coverage_yan_mg_l": float(policies[1].nutrition_mg_yan_l[0][1]),
        "temperature_changes": int(search_record["temperature_changes"]),
        "total_thermal_variation_c": float(search_record["total_thermal_variation_c"]),
        "minimum_posterior_fim_eigenvalue": float(
            metrics["minimum_posterior_fim_eigenvalue"]
        ),
        "maximum_posterior_fim_condition_number": float(
            metrics["maximum_posterior_fim_condition_number"]
        ),
        "distinct_sampling_slots": len(
            {time_h for values in schedules.values() for time_h in values}
        ),
        "harmonized_score_loss_fraction": (
            0.0
            if harmonization is None
            else float(harmonization["score_loss_fraction"])
        ),
        "additional_shared_nonbasal_slots": (
            0
            if harmonization is None
            else int(harmonization["additional_shared_nonbasal_slots"])
        ),
        "operator_rounds": (
            sum(len(values) for values in schedules.values())
            if harmonization is None
            else int(harmonization["harmonized_operator_rounds"])
        ),
        "operator_round_reduction": (
            0
            if harmonization is None
            else int(harmonization["operator_round_reduction"])
        ),
        "sampling_gate": gate_verdict(checks, tuple(checks)),
    }
    return {
        "row": row,
        "checks": checks,
        "schedules": schedules,
        "schedule_frame": schedule_frame,
        "captures": captures.assign(
            strategy=strategy, dose_design=dose_design, sampling_mode=mode
        ),
        "conflicts": conflicts.assign(
            strategy=strategy, dose_design=dose_design, sampling_mode=mode
        ),
        "gains": gains,
        "posterior": posterior,
        "harmonization": harmonization,
    }


def _optimize_sampling_for_design(
    *,
    policies: tuple[DesignPolicy, ...],
    strategy: str,
    dose_design: str,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    partitions: dict[str, Any],
    prepared_cache: dict[tuple[Any, ...], Any],
    fim_cache: dict[tuple, np.ndarray],
    current_reference_gains: np.ndarray,
    search_record: dict[str, Any],
    seed_offset: int,
) -> tuple[list[dict[str, Any]], pd.DataFrame, pd.DataFrame]:
    sampling_cache, valid, drying = _prepare_sampling_context(
        policies, ensemble, model_config, partitions, prepared_cache
    )
    members = list(range(len(ensemble)))
    independent, restarts = _optimize_schedules(
        policies,
        valid,
        sampling_cache,
        members,
        prior,
        model_config,
        fim_cache,
        seed_offset=seed_offset,
    )

    def score_function(candidate: dict[str, tuple[float, ...]]) -> float:
        return float(
            _score(
                candidate,
                policies,
                sampling_cache,
                members,
                prior,
                model_config,
                fim_cache,
            )[0]
        )

    harmonized, harmonization = partially_harmonize_schedules(
        independent,
        valid,
        score_function,
        float(coverage_config["sampling"]["harmonized_maximum_score_loss_fraction"]),
    )
    results = [
        _sampling_mode_result(
            schedules=independent,
            policies=policies,
            strategy=strategy,
            dose_design=dose_design,
            mode="independent",
            sampling_cache=sampling_cache,
            prior=prior,
            model_config=model_config,
            fim_cache=fim_cache,
            current_reference_gains=current_reference_gains,
            search_record=search_record,
            harmonization=None,
        ),
        _sampling_mode_result(
            schedules=harmonized,
            policies=policies,
            strategy=strategy,
            dose_design=dose_design,
            mode="partially_harmonized",
            sampling_cache=sampling_cache,
            prior=prior,
            model_config=model_config,
            fim_cache=fim_cache,
            current_reference_gains=current_reference_gains,
            search_record=search_record,
            harmonization=harmonization,
        ),
    ]
    return results, restarts.assign(strategy=strategy, dose_design=dose_design), drying


def _reference_information_context(
    *,
    reference_policies: tuple[DesignPolicy, ...],
    source_schedule: dict[str, tuple[float, ...]],
    source_comparison: pd.DataFrame,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    partitions: dict[str, Any],
    prepared_cache: dict[tuple[Any, ...], Any],
    fim_cache: dict[tuple, np.ndarray],
) -> dict[str, Any]:
    sampling_cache, valid, drying = _prepare_sampling_context(
        reference_policies,
        ensemble,
        model_config,
        partitions,
        prepared_cache,
    )
    members = list(range(len(ensemble)))
    current_score, current_gains, current_posterior = _score(
        source_schedule,
        reference_policies,
        sampling_cache,
        members,
        prior,
        model_config,
        fim_cache,
    )
    declared = source_comparison.sort_values("ensemble_member")[
        "optimized_information_gain"
    ].to_numpy(dtype=float)
    if declared.shape != current_gains.shape:
        raise ValueError("Reference sampling comparison does not contain 64 members")
    maximum_source_difference = float(np.max(np.abs(declared - current_gains)))
    if maximum_source_difference > 1e-6:
        raise RuntimeError(
            "Recomputed information-reference gains disagree with the immutable source"
        )
    current_metrics = robust_information_metrics(
        current_gains, current_posterior, None, model_config
    )
    anchor = reference_policies[0]
    anchors = tuple(
        DesignPolicy(
            f"three_anchor_reference_{index + 1}",
            anchor.temperature_c,
            anchor.nutrition_mg_yan_l,
        )
        for index in range(3)
    )
    for member in members:
        source_cache = sampling_cache[(member, anchor.name)]
        for policy in anchors:
            sampling_cache[(member, policy.name)] = source_cache
    anchor_valid = {policy.name: valid[anchor.name] for policy in anchors}
    anchor_schedules, anchor_restarts = _optimize_schedules(
        anchors,
        anchor_valid,
        sampling_cache,
        members,
        prior,
        model_config,
        fim_cache,
        seed_offset=1000,
    )
    anchor_score, anchor_gains, anchor_posterior = _score(
        anchor_schedules,
        anchors,
        sampling_cache,
        members,
        prior,
        model_config,
        fim_cache,
    )
    declared_anchor = source_comparison.sort_values("ensemble_member")[
        "three_anchor_information_gain"
    ].to_numpy(dtype=float)
    anchor_source_difference = float(np.max(np.abs(declared_anchor - anchor_gains)))
    if anchor_source_difference > 1e-6:
        raise RuntimeError(
            "Recomputed three-anchor gains disagree with the immutable source"
        )
    anchor_metrics = robust_information_metrics(
        anchor_gains, anchor_posterior, None, model_config
    )
    return {
        "current": {
            "score": float(current_score),
            "gains": current_gains,
            "posterior": current_posterior,
            "metrics": current_metrics,
            "schedules": source_schedule,
            "policies": reference_policies,
        },
        "three_anchor": {
            "score": float(anchor_score),
            "gains": anchor_gains,
            "posterior": anchor_posterior,
            "metrics": anchor_metrics,
            "schedules": anchor_schedules,
            "policies": anchors,
            "restarts": anchor_restarts,
        },
        "source_reconciliation": {
            "maximum_abs_current_gain_difference": maximum_source_difference,
            "maximum_abs_three_anchor_gain_difference": anchor_source_difference,
            "verdict": "PASS",
        },
        "drying": drying,
    }


def _posterior_information_tables(
    design_results: list[dict[str, Any]],
    current_reference: dict[str, Any],
    anchor_reference: dict[str, Any],
    model_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    names = [
        column.removeprefix("aroma__") for column in parameter_columns(model_config)
    ]
    baseline_current_variance = np.median(
        np.stack(
            [
                np.diag(np.linalg.pinv(fim, rcond=1e-12))
                for fim in current_reference["posterior"]
            ]
        ),
        axis=0,
    )
    baseline_anchor_variance = np.median(
        np.stack(
            [
                np.diag(np.linalg.pinv(fim, rcond=1e-12))
                for fim in anchor_reference["posterior"]
            ]
        ),
        axis=0,
    )
    entries = [
        {
            "strategy": "I_reference_information",
            "dose_design": "N80",
            "sampling_mode": "independent_source_final",
            "posterior": current_reference["posterior"],
        },
        {
            "strategy": "three_anchor_reference",
            "dose_design": "N80",
            "sampling_mode": "independent_optimized",
            "posterior": anchor_reference["posterior"],
        },
    ]
    entries.extend(
        {
            "strategy": row["row"]["strategy"],
            "dose_design": row["row"]["dose_design"],
            "sampling_mode": row["row"]["sampling_mode"],
            "posterior": row["posterior"],
        }
        for row in design_results
    )
    parameter_rows = []
    ratio_rows = []
    correlation_rows = []
    eigen_rows = []
    for entry in entries:
        matrices = [0.5 * (fim + fim.T) for fim in entry["posterior"]]
        covariance = [np.linalg.pinv(fim, rcond=1e-12) for fim in matrices]
        variances = np.stack([np.diag(item) for item in covariance])
        median_precision = np.median(np.stack(matrices), axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(median_precision)
        least = eigenvectors[:, 0]
        for index, parameter in enumerate(names):
            median_variance = float(np.median(variances[:, index]))
            parameter_rows.append(
                {
                    "strategy": entry["strategy"],
                    "dose_design": entry["dose_design"],
                    "sampling_mode": entry["sampling_mode"],
                    "parameter": parameter,
                    "posterior_variance_median": median_variance,
                    "posterior_variance_q10": float(
                        np.quantile(variances[:, index], 0.1)
                    ),
                    "posterior_variance_q90": float(
                        np.quantile(variances[:, index], 0.9)
                    ),
                    "posterior_standard_deviation_median": math.sqrt(
                        max(median_variance, 0.0)
                    ),
                    "least_identified_direction_loading": float(least[index]),
                    "least_identified_direction_abs_loading": abs(float(least[index])),
                }
            )
            ratio_rows.append(
                {
                    "strategy": entry["strategy"],
                    "dose_design": entry["dose_design"],
                    "sampling_mode": entry["sampling_mode"],
                    "parameter": parameter,
                    "variance_ratio_vs_three_anchors": median_variance
                    / max(float(baseline_anchor_variance[index]), 1e-18),
                    "variance_ratio_vs_information_reference": median_variance
                    / max(float(baseline_current_variance[index]), 1e-18),
                }
            )
        correlations = []
        for covariance_matrix in covariance:
            scale = np.sqrt(np.maximum(np.diag(covariance_matrix), 1e-18))
            correlations.append(
                covariance_matrix / np.maximum(np.outer(scale, scale), 1e-18)
            )
        median_correlation = np.median(np.stack(correlations), axis=0)
        pairs = []
        for left in range(len(names)):
            for right in range(left + 1, len(names)):
                pairs.append((abs(float(median_correlation[left, right])), left, right))
        for rank, (_, left, right) in enumerate(
            sorted(pairs, reverse=True)[:5], start=1
        ):
            correlation_rows.append(
                {
                    "strategy": entry["strategy"],
                    "dose_design": entry["dose_design"],
                    "sampling_mode": entry["sampling_mode"],
                    "rank": rank,
                    "parameter_left": names[left],
                    "parameter_right": names[right],
                    "posterior_correlation_median": float(
                        median_correlation[left, right]
                    ),
                    "absolute_posterior_correlation_median": abs(
                        float(median_correlation[left, right])
                    ),
                }
            )
        member_eigenvalues = np.stack([np.linalg.eigvalsh(fim) for fim in matrices])
        for index in range(member_eigenvalues.shape[1]):
            eigen_rows.append(
                {
                    "strategy": entry["strategy"],
                    "dose_design": entry["dose_design"],
                    "sampling_mode": entry["sampling_mode"],
                    "eigenvalue_rank_ascending": index + 1,
                    "fim_eigenvalue_median": float(
                        np.median(member_eigenvalues[:, index])
                    ),
                    "fim_eigenvalue_q10": float(
                        np.quantile(member_eigenvalues[:, index], 0.1)
                    ),
                    "fim_eigenvalue_minimum": float(
                        np.min(member_eigenvalues[:, index])
                    ),
                }
            )
    return (
        pd.DataFrame(parameter_rows),
        pd.DataFrame(ratio_rows),
        pd.DataFrame(correlation_rows),
        pd.DataFrame(eigen_rows),
    )


def _ensemble_prediction_bands(
    archetype_policies: dict[str, DesignPolicy],
    ensemble: pd.DataFrame,
    model_config: dict[str, Any],
    partitions: dict[str, Any],
) -> pd.DataFrame:
    model = __import__(
        "pilot_2026.adaptive_design.pilot_calibration",
        fromlist=["_model_module"],
    )._model_module()
    time_h = np.arange(0.0, 504.0 + 3.0, 3.0)
    initial_aroma = {
        "ethyl_acetate": 1152.067155825675,
        "ethyl_octanoate": 3.48070575912168,
        "isoamyl_acetate": 2.5,
    }
    rows = []
    for archetype, policy in archetype_policies.items():
        design = _future_design(policy, model_config)
        temperature = actuator_temperature_trajectory(policy, model_config)
        physical = np.interp(
            time_h,
            temperature["time_h"],
            temperature["physical_temperature_c"],
        )
        by_metric: dict[str, list[np.ndarray]] = defaultdict(list)
        for member_index in range(len(ensemble)):
            member = ensemble.iloc[member_index]
            theta = _theta(member)
            core = model.simulate(
                design,
                theta,
                time_h,
                sample_event_order=model_config["operations"][
                    "numerical_sample_event_order"
                ],
            )
            if core is None:
                raise RuntimeError(
                    f"Prediction simulation failed for {archetype}, member {member_index}"
                )
            sugar = core["G"].to_numpy(dtype=float) + core["F"].to_numpy(dtype=float)
            by_metric["residual_sugar_g_l"].append(sugar)
            by_metric["biomass_g_l"].append(core["X"].to_numpy(dtype=float))
            by_metric["yan_g_l"].append(core["N"].to_numpy(dtype=float))
            captured_total = np.zeros_like(time_h)
            for species in SPECIES:
                forcing = _forcing_from_core(
                    design,
                    core,
                    theta,
                    partitions[species],
                    species,
                    initial_aroma[species],
                )
                liquid, captured = simulate_aroma(
                    forcing, _aroma_log_values(member, species)
                )
                by_metric[f"{species}_ug_l"].append(liquid)
                captured_total += captured
            by_metric["cumulative_captured_mass_ug"].append(captured_total)
        by_metric["physical_temperature_c"] = [physical] * len(ensemble)
        sugar_array = np.stack(by_metric["residual_sugar_g_l"])
        drying = []
        target = float(model_config["completion"]["residual_sugar_g_l"])
        for member_sugar in sugar_array:
            indices = np.flatnonzero(member_sugar <= target)
            drying.append(float(time_h[indices[0]]) if len(indices) else math.inf)
        finite_drying = np.asarray([value for value in drying if np.isfinite(value)])
        for metric, values in by_metric.items():
            array = np.stack(values)
            for index, current_time in enumerate(time_h):
                rows.append(
                    {
                        "archetype": archetype,
                        "policy": policy.name,
                        "time_h": float(current_time),
                        "metric": metric,
                        "q10": float(np.quantile(array[:, index], 0.1)),
                        "median": float(np.median(array[:, index])),
                        "q90": float(np.quantile(array[:, index], 0.9)),
                        "minimum": float(np.min(array[:, index])),
                        "maximum": float(np.max(array[:, index])),
                        "ensemble_members": len(ensemble),
                        "drying_time_q10_h": float(np.quantile(finite_drying, 0.1)),
                        "drying_time_median_h": float(np.median(finite_drying)),
                        "drying_time_q90_h": float(np.quantile(finite_drying, 0.9)),
                        "completion_probability_504h": float(
                            np.mean(np.isfinite(np.asarray(drying)))
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _physical_trajectory_bands(
    policies: dict[str, DesignPolicy],
    model_config: dict[str, Any],
) -> pd.DataFrame:
    from pilot_2026.adaptive_design.final_search_logic import (
        approved_actuator_scenarios,
    )

    rows = []
    for label, policy in policies.items():
        trajectories = []
        common_time = None
        for scenario in approved_actuator_scenarios(model_config):
            trajectory = actuator_temperature_trajectory(policy, model_config, scenario)
            if common_time is None:
                common_time = np.asarray(trajectory["time_h"], dtype=float)
            trajectories.append(
                np.asarray(trajectory["physical_temperature_c"], dtype=float)
            )
        array = np.stack(trajectories)
        assert common_time is not None
        indices = np.arange(0, len(common_time), 12, dtype=int)
        if indices[-1] != len(common_time) - 1:
            indices = np.append(indices, len(common_time) - 1)
        for index in indices:
            rows.append(
                {
                    "archetype": label,
                    "policy": policy.name,
                    "time_h": float(common_time[index]),
                    "minimum_physical_temperature_c": float(np.min(array[:, index])),
                    "q10_physical_temperature_c": float(
                        np.quantile(array[:, index], 0.1)
                    ),
                    "median_physical_temperature_c": float(np.median(array[:, index])),
                    "q90_physical_temperature_c": float(
                        np.quantile(array[:, index], 0.9)
                    ),
                    "maximum_physical_temperature_c": float(np.max(array[:, index])),
                    "approved_actuator_scenarios": array.shape[0],
                }
            )
    return pd.DataFrame(rows)


def _finish_figure(fig: plt.Figure, path: Path, subtitle: str) -> None:
    fig.text(0.01, 0.985, subtitle, ha="left", va="top", fontsize=9, color="#4B5563")
    fig.text(
        0.5,
        0.012,
        WATERMARK,
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color="#8B1E3F",
    )
    fig.tight_layout(rect=(0.02, 0.05, 0.98, 0.94))
    fig.savefig(filesystem_path(path), dpi=170, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def _prediction_figure(prediction: pd.DataFrame, archetype: str, path: Path) -> None:
    metrics = (
        ("residual_sugar_g_l", "Residual sugar (g/L)"),
        ("biomass_g_l", "Biomass (g/L)"),
        ("yan_g_l", "YAN state (g/L)"),
        ("physical_temperature_c", "Physical wine temperature (\u00b0C)"),
        ("ethyl_acetate_ug_l", "Ethyl acetate (\u00b5g/L)"),
        ("ethyl_octanoate_ug_l", "Ethyl octanoate (\u00b5g/L)"),
        ("isoamyl_acetate_ug_l", "Isoamyl acetate (\u00b5g/L)"),
        ("cumulative_captured_mass_ug", "Cumulative captured mass (\u00b5g)"),
    )
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    subset = prediction[prediction.archetype.eq(archetype)]
    for axis, (metric, label) in zip(axes.flat, metrics):
        data = subset[subset.metric.eq(metric)].sort_values("time_h")
        x = data["time_h"].to_numpy(dtype=float)
        axis.fill_between(
            x,
            data["q10"].to_numpy(dtype=float),
            data["q90"].to_numpy(dtype=float),
            color=PALETTE["blue"],
            alpha=0.18,
            label="ensemble q10\u2013q90",
        )
        axis.plot(x, data["median"], color=PALETTE["blue"], linewidth=1.8)
        if metric == "residual_sugar_g_l":
            axis.axhline(4.0, color=PALETTE["ink"], linestyle="--", linewidth=1)
        if metric == "physical_temperature_c":
            axis.axhspan(15.0, 27.0, color=PALETTE["olive"], alpha=0.08)
        axis.set_title(label)
        axis.grid(True, linewidth=0.6)
        axis.set_xlabel("Time from campaign start (h)")
    _finish_figure(
        fig,
        path,
        f"{archetype.replace('_', ' ').title()} predictions; 64-member ensemble, q10\u2013q90 bands",
    )


def _generate_figures(
    *,
    run_dir: Path,
    comparison: pd.DataFrame,
    candidate_actions: pd.DataFrame,
    parameter_ratios: pd.DataFrame,
    predictions: pd.DataFrame,
    physical_bands: pd.DataFrame,
    recommended: dict[str, Any],
    selected_policies_by_strategy: dict[str, tuple[DesignPolicy, ...]],
) -> tuple[list[Path], dict[str, Any]]:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#374151",
            "axes.labelcolor": PALETTE["ink"],
            "axes.titlecolor": "#111827",
            "xtick.color": "#374151",
            "ytick.color": "#374151",
            "grid.color": "#D1D5DB",
            "grid.alpha": 0.55,
        }
    )
    paths = {name: run_dir / name for name in FIGURE_NAMES}

    actions = candidate_actions[
        candidate_actions.action.eq("nutrition")
        & candidate_actions.strategy.isin(["II_diagonal_1", "III_diagonal_2"])
    ].copy()
    fig, axis = plt.subplots(figsize=(10, 6))
    for strategy, group in actions.groupby("strategy"):
        axis.scatter(
            group["first_24h_mean_setpoint_c"],
            group["time_h"],
            s=75,
            alpha=0.75,
            label=strategy,
        )
    axis.axvline(18.0, color=PALETTE["blue"], linestyle="--", linewidth=1)
    axis.axvline(22.0, color=PALETTE["orange"], linestyle="--", linewidth=1)
    axis.axhspan(0.0, 6.0, color=PALETTE["gold"], alpha=0.09)
    axis.axhspan(42.0, 50.0, color=PALETTE["pink"], alpha=0.09)
    axis.set_xlabel("Mean setpoint during first 24 h (\u00b0C)")
    axis.set_ylabel("Nutrition time (h)")
    axis.set_title("Operational thermal\u2013nutrition coverage map")
    axis.grid(True)
    axis.legend(frameon=False)
    _finish_figure(
        fig,
        paths["operational_coverage_map.png"],
        "Hard archetype regions and selected candidates",
    )

    fig, axis = plt.subplots(figsize=(12, 6))
    styles = {"II_diagonal_1": "-", "III_diagonal_2": "--"}
    for strategy, policies in selected_policies_by_strategy.items():
        for policy in policies:
            if "anchor" in policy.name:
                continue
            x = np.arange(len(policy.temperature_c) + 1) * 12.0
            y = np.r_[policy.temperature_c, policy.temperature_c[-1]]
            axis.step(
                x,
                y,
                where="post",
                linestyle=styles[strategy],
                linewidth=1.8,
                label=policy.name,
            )
    axis.axhspan(15.0, 27.0, color=PALETTE["olive"], alpha=0.07)
    axis.set_xlim(0, 168)
    axis.set_xlabel("Time from campaign start (h)")
    axis.set_ylabel("Setpoint (\u00b0C)")
    axis.set_title("Coverage strategy temperature profiles")
    axis.grid(True)
    axis.legend(frameon=False, fontsize=7, ncol=2)
    _finish_figure(
        fig,
        paths["coverage_strategy_profiles.png"],
        "Selected 12 h policies; explicit hold continues to 504 h",
    )

    fig, axis = plt.subplots(figsize=(10, 6))
    plotted = comparison[comparison.sampling_mode.notna()].copy()
    cell_count = np.where(plotted.strategy.eq("I_reference_information"), 1, 2)
    colors = [
        PALETTE["grey"] if value == 1 else PALETTE["blue"] for value in cell_count
    ]
    axis.scatter(cell_count, plotted["robust_score"], c=colors, s=70, alpha=0.8)
    for _, row in plotted.iterrows():
        axis.annotate(
            f"{row['strategy'].split('_')[0]}\n{row['dose_design']}\n{row['sampling_mode']}",
            (
                1 if row["strategy"] == "I_reference_information" else 2,
                row["robust_score"],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=6,
        )
    axis.set_xticks([1, 2], ["one represented cell", "two diagonal cells"])
    axis.set_ylabel("Robust Bayesian D-information score")
    axis.set_title("Operational coverage versus information")
    axis.grid(True, axis="y")
    _finish_figure(
        fig,
        paths["coverage_vs_information_tradeoff.png"],
        "64-member comparison after sampling optimization",
    )

    recommended_rows = parameter_ratios[
        parameter_ratios.strategy.eq(recommended.get("strategy"))
        & parameter_ratios.dose_design.eq(recommended.get("dose_design"))
        & parameter_ratios.sampling_mode.eq(recommended.get("sampling_mode"))
    ].copy()
    fig, axis = plt.subplots(figsize=(12, 6))
    x = np.arange(len(recommended_rows))
    axis.bar(
        x - 0.18,
        recommended_rows["variance_ratio_vs_information_reference"],
        width=0.36,
        color=PALETTE["blue"],
        label="vs information reference",
    )
    axis.bar(
        x + 0.18,
        recommended_rows["variance_ratio_vs_three_anchors"],
        width=0.36,
        color=PALETTE["gold"],
        label="vs three anchors",
    )
    axis.axhline(1.0, color=PALETTE["ink"], linestyle="--", linewidth=1)
    axis.set_xticks(
        x, recommended_rows["parameter"], rotation=60, ha="right", fontsize=7
    )
    axis.set_ylabel("Posterior variance ratio")
    axis.set_title("Per-parameter posterior variance comparison")
    axis.legend(frameon=False)
    axis.grid(True, axis="y")
    _finish_figure(
        fig,
        paths["parameter_variance_comparison.png"],
        "Median across 64 ensemble members",
    )

    for archetype in ("cold_early", "cold_late", "warm_early", "warm_late"):
        _prediction_figure(
            predictions, archetype, paths[f"{archetype}_predictions.png"]
        )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    for axis, (archetype, group) in zip(
        axes.flat, physical_bands.groupby("archetype", sort=True)
    ):
        group = group.sort_values("time_h")
        x = group["time_h"].to_numpy(dtype=float)
        axis.fill_between(
            x,
            group["minimum_physical_temperature_c"].to_numpy(dtype=float),
            group["maximum_physical_temperature_c"].to_numpy(dtype=float),
            color=PALETTE["blue"],
            alpha=0.16,
        )
        axis.plot(
            x,
            group["median_physical_temperature_c"],
            color=PALETTE["blue"],
            linewidth=1.5,
        )
        axis.axhline(15.0, color=PALETTE["ink"], linestyle="--", linewidth=0.9)
        axis.axhline(27.0, color=PALETTE["ink"], linestyle="--", linewidth=0.9)
        axis.set_title(archetype.replace("_", " ").title())
        axis.set_xlabel("Time (h)")
        axis.set_ylabel("Physical wine temperature (\u00b0C)")
        axis.grid(True)
    _finish_figure(
        fig,
        paths["physical_temperature_envelope.png"],
        "Minimum\u2013maximum envelope across 243 approved actuator scenarios",
    )

    fig, axis = plt.subplots(figsize=(12, 6))
    axis.axis("off")
    boxes = [
        (0.04, 0.58, 0.25, 0.25, "Wave 1 (3)\nanchor\ncold/early\nwarm/late"),
        (
            0.38,
            0.58,
            0.25,
            0.25,
            "Wave 2 (2)\ncold/late\nwarm/early\nupdated after Wave 1",
        ),
        (0.72, 0.58, 0.24, 0.25, "Wave 3 (1)\nadaptive confirmation\nor replicate"),
    ]
    for x0, y0, width, height, label in boxes:
        rectangle = plt.Rectangle(
            (x0, y0),
            width,
            height,
            facecolor="#EEF4FA",
            edgecolor=PALETTE["blue"],
            linewidth=1.6,
        )
        axis.add_patch(rectangle)
        axis.text(
            x0 + width / 2,
            y0 + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=10,
        )
    for start, end in ((0.29, 0.38), (0.63, 0.72)):
        axis.annotate(
            "",
            xy=(end, 0.705),
            xytext=(start, 0.705),
            arrowprops={"arrowstyle": "->", "color": PALETTE["ink"], "lw": 1.5},
        )
    axis.text(
        0.5,
        0.32,
        "Update rule: posterior reduction + least identified direction + model discrepancy + replicate need + aroma/sensory coverage",
        ha="center",
        va="center",
        fontsize=9,
        wrap=True,
    )
    axis.set_title("Adaptive Wave 3\u20132\u20131 operational coverage plan")
    _finish_figure(
        fig,
        paths["wave_3_2_1_coverage_plan.png"],
        "Architecture only; no physical Wave 2 or Wave 3 instructions",
    )

    qa_rows = []
    for name in FIGURE_NAMES:
        path = paths[name]
        image = np.asarray(plt.imread(filesystem_path(path)), dtype=float)
        qa_rows.append(
            {
                "figure": name,
                "exists": filesystem_path(path).is_file(),
                "bytes": int(filesystem_path(path).stat().st_size),
                "height_px": int(image.shape[0]),
                "width_px": int(image.shape[1]),
                "pixel_dynamic_range": float(np.max(image) - np.min(image)),
                "pixel_standard_deviation": float(np.std(image)),
                "watermark_added_by_shared_export_function": True,
                "passed": bool(
                    filesystem_path(path).stat().st_size >= 10_000
                    and image.shape[0] >= 500
                    and image.shape[1] >= 700
                    and float(np.max(image) - np.min(image)) >= 0.25
                    and float(np.std(image)) >= 0.01
                ),
            }
        )
    qa = {
        "verdict": "PASS" if all(row["passed"] for row in qa_rows) else "FAIL",
        "exact_expected_inventory": list(paths) == list(FIGURE_NAMES),
        "expected_watermark": WATERMARK,
        "figures": qa_rows,
    }
    return list(paths.values()), qa


def _flat_search_record(row: dict[str, Any]) -> dict[str, Any]:
    excluded = {
        "vector",
        "policies",
        "evaluations",
        "archetype_checks",
        "contrast_checks",
    }
    return {key: value for key, value in row.items() if key not in excluded}


def _candidate_action_rows(
    strategy: str,
    dose_design: str,
    candidate_id: str,
    policies: tuple[DesignPolicy, ...],
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
) -> pd.DataFrame:
    rows = []
    for policy in policies:
        metrics = coverage_policy_metrics(policy, model_config, coverage_config)
        controller = explicit_controller_blocks(policy, coverage_config)
        for item in controller.itertuples():
            rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "candidate_id": candidate_id,
                    "policy": policy.name,
                    "thermal_archetype": metrics.thermal_archetype,
                    "nutrition_archetype": metrics.nutrition_archetype,
                    "first_24h_mean_setpoint_c": metrics.first_24h_mean_setpoint_c,
                    "action": "temperature_setpoint",
                    "time_h": float(item.start_h),
                    "end_h": float(item.end_h),
                    "value": float(item.setpoint_c),
                    "unit": "degC",
                    "controller_phase": item.phase,
                    "explicit_hold_last_setpoint": bool(item.hold_last_setpoint),
                    "sample_before_action_if_coincident": False,
                }
            )
        for time_h, amount in policy.nutrition_mg_yan_l:
            rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "candidate_id": candidate_id,
                    "policy": policy.name,
                    "thermal_archetype": metrics.thermal_archetype,
                    "nutrition_archetype": metrics.nutrition_archetype,
                    "first_24h_mean_setpoint_c": metrics.first_24h_mean_setpoint_c,
                    "action": "nutrition",
                    "time_h": float(time_h),
                    "end_h": math.nan,
                    "value": float(amount),
                    "unit": "mgYAN/L",
                    "controller_phase": "manual_action",
                    "explicit_hold_last_setpoint": False,
                    "sample_before_action_if_coincident": True,
                }
            )
        for transition in robust_initial_jump_rows(policy, coverage_config):
            rows.append(
                {
                    "strategy": strategy,
                    "dose_design": dose_design,
                    "candidate_id": candidate_id,
                    "policy": policy.name,
                    "thermal_archetype": metrics.thermal_archetype,
                    "nutrition_archetype": metrics.nutrition_archetype,
                    "first_24h_mean_setpoint_c": metrics.first_24h_mean_setpoint_c,
                    "action": "initial_transition_audit",
                    "time_h": 0.0,
                    "end_h": math.nan,
                    "value": float(policy.temperature_c[0]),
                    "unit": "degC",
                    "controller_phase": "initial_transition",
                    "explicit_hold_last_setpoint": False,
                    "sample_before_action_if_coincident": False,
                    "initial_temperature_scenario_c": transition[
                        "initial_temperature_scenario_c"
                    ],
                    "initial_setpoint_jump_c": transition["initial_setpoint_jump_c"],
                    "initial_jump_within_5c": transition["initial_jump_within_5c"],
                }
            )
    return pd.DataFrame(rows)


def _reference_corrected_physical_record(
    policies: tuple[DesignPolicy, ...],
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
) -> dict[str, Any]:
    summaries = [
        _physical_summary_cached(policy, model_config, coverage_config, physical_cache)[
            1
        ]
        for policy in policies
    ]
    metrics = [
        coverage_policy_metrics(policy, model_config, coverage_config)
        for policy in policies
    ]
    duplicate = bool(
        np.allclose(policies[1].temperature_c, policies[2].temperature_c, atol=0.1)
        and math.isclose(
            policies[1].nutrition_mg_yan_l[0][0],
            policies[2].nutrition_mg_yan_l[0][0],
            abs_tol=1e-9,
        )
    )
    return {
        "coverage_pass": False,
        "thermal_pass": bool(
            all(row["initial_jump_robust_pass"] for row in summaries)
            and all(row["physical_temperature_robust_pass"] for row in summaries)
        ),
        "nutrition_pass": True,
        "maximum_internal_jump_c": max(row.maximum_internal_jump_c for row in metrics),
        "maximum_initial_jump_c": max(row.maximum_initial_jump_c for row in metrics),
        "minimum_robust_physical_temperature_c": min(
            float(row["minimum_physical_temperature_c"]) for row in summaries
        ),
        "maximum_robust_physical_temperature_c": max(
            float(row["maximum_physical_temperature_c"]) for row in summaries
        ),
        "temperature_changes": sum(row.temperature_changes for row in metrics),
        "total_thermal_variation_c": sum(
            row.total_thermal_variation_c for row in metrics
        ),
        "coverage_treatments_duplicate": duplicate,
        "corrected_failure_reasons": [
            "coverage lacks a cold treatment and late/early diagonal",
            "warm initial setpoints exceed the robust 5 C transition limit",
            "robust physical wine temperature exceeds 27 C",
            "candidate A and B are operationally equivalent warm/early treatments",
        ],
    }


def _recommendation(
    sampling_results: list[dict[str, Any]],
    search_by_variant: dict[tuple[str, str], dict[str, Any]],
    current_score: float,
    anchor_score: float,
    coverage_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    rows = []
    for result in sampling_results:
        row = dict(result["row"])
        search = search_by_variant[(row["strategy"], row["dose_design"])]
        selected = search["selected"]
        row["coverage_gate"] = "PASS" if selected["coverage_pass"] else "FAIL"
        row["thermal_gate"] = "PASS" if selected["thermal_pass"] else "FAIL"
        row["nutrition_gate"] = "PASS" if selected["nutrition_pass"] else "FAIL"
        information_pass = bool(
            search["convergence"]["passed"]
            and row["sampling_gate"] == "PASS"
            and selected.get("local_refinement_qualified", True)
        )
        row["information_gate"] = "PASS" if information_pass else "FAIL"
        row["retained_information_improvement_fraction"] = (
            retained_information_improvement_fraction(
                row["robust_score"], current_score, anchor_score
            )
        )
        row["information_cost_absolute_vs_current"] = float(
            row["robust_score"]
        ) - float(current_score)
        row["information_cost_relative_vs_current"] = (
            float(row["robust_score"]) / max(abs(float(current_score)), 1e-12) - 1.0
        )
        row["all_operational_gates_pass"] = all(
            row[name] == "PASS"
            for name in (
                "coverage_gate",
                "thermal_gate",
                "nutrition_gate",
                "information_gate",
            )
        )
        rows.append(row)
    comparison = pd.DataFrame(rows)
    eligible = comparison[comparison.all_operational_gates_pass].copy()
    if eligible.empty:
        return {
            "recommended": False,
            "reason": "no coverage candidate passes every operational and information gate",
            "best_information_design": "I_reference_information",
            "best_coverage_design": None,
            "best_compromise": None,
        }, comparison
    best_score = float(eligible["robust_score"].max())
    near = eligible[
        best_score - eligible["robust_score"]
        <= float(coverage_config["objective"]["near_best_fraction"])
        * max(abs(best_score), 1e-12)
        + 1e-12
    ].copy()
    near["dose_margin_rank"] = near["dose_design"].map(
        {"N76_all_conservative": 0, "N76_historical_anchor": 1, "N80": 2}
    )
    near["sampling_workload_rank"] = near["sampling_mode"].map(
        {"partially_harmonized": 0, "independent": 1}
    )
    selected_row = near.sort_values(
        [
            "maximum_initial_jump_c",
            "maximum_robust_physical_temperature_c",
            "total_thermal_variation_c",
            "temperature_changes",
            "dose_margin_rank",
            "sampling_workload_rank",
            "robust_score",
        ],
        ascending=[True, True, True, True, True, True, False],
    ).iloc[0]
    retained = float(selected_row["retained_information_improvement_fraction"])
    threshold = float(
        coverage_config["recommendation"][
            "minimum_retained_information_improvement_fraction"
        ]
    )
    justified_below = bool(
        coverage_config["recommendation"][
            "coverage_cost_can_be_justified_below_threshold"
        ]
    )
    recommended = bool(retained >= threshold or justified_below)
    return {
        "recommended": recommended,
        "strategy": str(selected_row["strategy"]),
        "dose_design": str(selected_row["dose_design"]),
        "sampling_mode": str(selected_row["sampling_mode"]),
        "robust_score": float(selected_row["robust_score"]),
        "retained_information_improvement_fraction": retained,
        "minimum_recommended_retention_fraction": threshold,
        "meets_90pct_information_improvement_rule": bool(retained >= threshold),
        "coverage_justification_used": bool(retained < threshold and justified_below),
        "coverage_justification": (
            "replaces duplicated warm/early treatments with two unambiguous diagonal "
            "thermal-nutrition cells while retaining explicit, quantified information cost"
        ),
        "best_information_design": "I_reference_information",
        "best_coverage_design": str(
            eligible.sort_values("robust_score", ascending=False).iloc[0]["strategy"]
        ),
        "best_compromise": str(selected_row["strategy"]),
        "candidate_id": search_by_variant[
            (str(selected_row["strategy"]), str(selected_row["dose_design"]))
        ]["selected"]["candidate_id"],
        "review_mapping_only": coverage_config["review_mapping_only"],
        **closed_authorization(),
    }, comparison


def _review_name(path: Path) -> str:
    return f"{path.stem}{REVIEW_SUFFIX}{path.suffix}"


def _build_review_package(
    run_dir: Path,
    outputs_by_name: dict[str, Path],
    figure_paths: list[Path],
    recommendation: dict[str, Any],
) -> list[Path]:
    package = run_dir / "coverage_review_package"
    routing = {
        "strategy_summary": [
            "coverage_strategy_comparison.csv",
            "coverage_vs_information_reference.csv",
            "coverage_archetype_gate.json",
            "wave1_recommended_coverage_candidate.json",
        ],
        "profiles": [
            "coverage_candidate_actions.csv",
            "physical_temperature_envelope_by_candidate.csv",
            "robust_setpoint_envelope.json",
        ],
        "nutrition": ["nutrition_margin_comparison.csv"],
        "sampling": [
            "sampling_strategy_comparison.csv",
            "optimized_sampling_schedules.csv",
        ],
        "capture": ["capture_interval_summary.csv"],
        "wave_plan": [
            "wave2_complementary_candidates.csv",
            "wave3_adaptive_decision_rule.json",
        ],
    }
    copied = []
    for folder, names in routing.items():
        destination = package / folder
        filesystem_path(destination).mkdir(parents=True, exist_ok=False)
        for name in names:
            source = outputs_by_name[name]
            target = destination / _review_name(source)
            shutil.copy2(filesystem_path(source), filesystem_path(target))
            copied.append(target)
    figure_destination = package / "figures"
    filesystem_path(figure_destination).mkdir(parents=True, exist_ok=False)
    for source in figure_paths:
        target = figure_destination / _review_name(source)
        shutil.copy2(filesystem_path(source), filesystem_path(target))
        copied.append(target)
    authorization_destination = package / "authorization"
    filesystem_path(authorization_destination).mkdir(parents=True, exist_ok=False)
    authorization_path = authorization_destination / (
        "authorization_state" + REVIEW_SUFFIX + ".json"
    )
    write_json(
        authorization_path,
        {
            "watermark": WATERMARK,
            "review_package_only": True,
            "recommended_candidate": recommendation,
            **closed_authorization(),
        },
    )
    copied.append(authorization_path)
    package_manifest = package / ("package_manifest" + REVIEW_SUFFIX + ".json")
    write_json(
        package_manifest,
        {
            "watermark": WATERMARK,
            "files": [relative_or_absolute(path) for path in copied],
            "all_names_have_review_only_suffix": all(
                REVIEW_SUFFIX in path.name for path in copied
            ),
            **closed_authorization(),
        },
    )
    copied.append(package_manifest)
    return copied


def _verify_complete_source_run(source_run: Path) -> dict[str, Any]:
    manifest = _read_json(source_run / "run_manifest.json")
    checked = 0
    canonical_lf = 0
    failures = []
    for declared_path, declared in manifest.get("outputs", {}).items():
        path = _run_path(declared_path)
        try:
            actual_hash = sha256_file(path)
            actual_bytes = filesystem_path(path).stat().st_size
            mode = "byte_exact"
            declared_bytes = int(declared["bytes"])
            if actual_hash != declared["sha256"]:
                raw = filesystem_path(path).read_bytes()
                normalized = raw.replace(b"\r\n", b"\n")
                normalized_hash = __import__("hashlib").sha256(normalized).hexdigest()
                if normalized_hash == declared["sha256"]:
                    mode = "canonical_lf_text"
                    canonical_lf += 1
                    if len(normalized) != declared_bytes:
                        failures.append(f"bytes:{declared_path}")
                else:
                    failures.append(f"hash:{declared_path}")
            elif actual_bytes != declared_bytes:
                failures.append(f"bytes:{declared_path}")
            checked += 1
        except (OSError, KeyError, TypeError, ValueError) as error:
            failures.append(f"missing:{declared_path}:{error}")
    return {
        "source_run": relative_or_absolute(source_run),
        "source_manifest_sha256": sha256_file(source_run / "run_manifest.json"),
        "declared_outputs": len(manifest.get("outputs", {})),
        "checked_outputs": checked,
        "canonical_lf_matches": canonical_lf,
        "failures": failures,
        "verdict": "PASS" if not failures and checked else "FAIL",
    }


def _copy_completed_run_payload(
    source_run: Path, destination_run: Path
) -> tuple[list[Path], dict[str, Path]]:
    outputs: list[Path] = []
    by_name: dict[str, Path] = {}
    source_fs = filesystem_path(source_run)
    destination_fs = filesystem_path(destination_run)
    for item in source_fs.iterdir():
        if item.name in {"run_manifest.json", "coverage_review_package"}:
            continue
        target = destination_fs / item.name
        if item.is_dir():
            shutil.copytree(item, target)
            for copied in target.rglob("*"):
                if copied.is_file():
                    outputs.append(destination_run / copied.relative_to(destination_fs))
        else:
            shutil.copy2(item, target)
            standard = destination_run / item.name
            outputs.append(standard)
            by_name[item.name] = standard
    return outputs, by_name


def _repair_physical_candidate_labels(run_dir: Path) -> dict[str, Any]:
    actions_path = run_dir / "coverage_candidate_actions.csv"
    physical_path = run_dir / "physical_temperature_envelope_by_candidate.csv"
    actions = pd.read_csv(filesystem_path(actions_path))
    physical = pd.read_csv(filesystem_path(physical_path))
    mapping = (
        actions.loc[
            actions.action.eq("initial_transition_audit"),
            ["strategy", "dose_design", "candidate_id", "policy", "value"],
        ]
        .drop_duplicates()
        .rename(
            columns={
                "policy": "expected_candidate",
                "value": "initial_setpoint_c",
            }
        )
    )
    key_columns = [
        "strategy",
        "dose_design",
        "candidate_id",
        "initial_setpoint_c",
    ]
    if mapping.duplicated(key_columns, keep=False).any():
        ambiguous = mapping.loc[
            mapping.duplicated(key_columns, keep=False),
            key_columns + ["expected_candidate"],
        ]
        raise RuntimeError(
            "Ambiguous controller-policy label mapping: "
            f"{ambiguous.to_dict('records')}"
        )
    expected = physical[key_columns].merge(
        mapping,
        on=key_columns,
        how="left",
        validate="many_to_one",
    )["expected_candidate"]
    if expected.isna().any():
        missing = physical.loc[
            expected.isna(),
            key_columns,
        ].drop_duplicates()
        raise RuntimeError(
            "Cannot restore physical candidate labels: " f"{missing.to_dict('records')}"
        )
    repaired = physical.copy()
    mismatches_before = int((repaired["candidate"] != expected).sum())
    repaired["candidate"] = expected.to_numpy()
    non_label_columns = [column for column in physical if column != "candidate"]
    pd.testing.assert_frame_equal(
        physical[non_label_columns],
        repaired[non_label_columns],
        check_exact=True,
    )
    repaired.to_csv(filesystem_path(physical_path), index=False, lineterminator="\n")
    return {
        "rows": len(repaired),
        "labels_corrected": mismatches_before,
        "all_rows_match_controller_policy": bool(
            repaired["candidate"].eq(expected).all()
        ),
        "all_non_label_columns_byte_value_identical": True,
        "repair_key": key_columns,
    }


def _repackage_completed_run(
    source_run: Path,
    result_root: Path,
    coverage_config: dict[str, Any],
    git_snapshot: dict[str, Any],
) -> Path:
    source_run = source_run.resolve()
    source_audit = _verify_complete_source_run(source_run)
    if source_audit["verdict"] != "PASS":
        raise RuntimeError(f"Source coverage run audit failed: {source_audit}")
    source_manifest = _read_json(source_run / "run_manifest.json")
    if source_manifest.get("status") != "completed":
        raise RuntimeError("Only a completed immutable coverage run may be repackaged")
    qualification_config = {
        "coverage": coverage_config,
        "qualification_action": "restore_cached_physical_candidate_labels",
        "source_run_manifest_sha256": source_audit["source_manifest_sha256"],
        "numeric_results_reused": True,
    }
    run_dir = create_immutable_run_directory(
        result_root, "wave1_operational_coverage", qualification_config
    )
    outputs, outputs_by_name = _copy_completed_run_payload(source_run, run_dir)
    repair = _repair_physical_candidate_labels(run_dir)
    runtime_path = run_dir / "runtime_summary.json"
    runtime = _read_json(runtime_path)
    runtime.update(
        {
            "qualification_repackaged_from": relative_or_absolute(source_run),
            "source_run_manifest_sha256": source_audit["source_manifest_sha256"],
            "numeric_search_sampling_and_prediction_results_reused": True,
            "physical_candidate_label_integrity_repair": repair,
            "watermark": WATERMARK,
            **closed_authorization(),
        }
    )
    write_json(runtime_path, runtime)
    recommendation = _read_json(run_dir / "wave1_recommended_coverage_candidate.json")
    figure_paths = [run_dir / name for name in FIGURE_NAMES]
    review_outputs = _build_review_package(
        run_dir, outputs_by_name, figure_paths, recommendation
    )
    outputs.extend(review_outputs)
    runs = {
        name: _run_path(value) for name, value in coverage_config["source_runs"].items()
    }
    closed_payloads = (source_manifest, runtime, recommendation)
    gate_checks = {
        "source_completed_run_hashes_verified": source_audit["verdict"] == "PASS",
        "source_numerical_gate_pass": source_manifest["gate"]["verdict"] == "PASS",
        "physical_candidate_labels_match_controller_policy": repair[
            "all_rows_match_controller_policy"
        ],
        "only_label_column_changed_in_physical_table": repair[
            "all_non_label_columns_byte_value_identical"
        ],
        "review_package_rebuilt_with_review_only_suffix": bool(review_outputs)
        and all(REVIEW_SUFFIX in path.name for path in review_outputs),
        "physical_flags_closed": all(
            payload.get(key) == value
            for payload in closed_payloads
            for key, value in AUTHORIZATION_FLAGS.items()
        ),
        "watermark_preserved": all(
            payload.get("watermark") == WATERMARK for payload in closed_payloads
        ),
    }
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_operational_coverage",
        config=qualification_config,
        sources={
            "source_coverage_run_manifest": source_run / "run_manifest.json",
            "coverage_config": COVERAGE_CONFIG_PATH,
            "adapter_manifest": runs["adapter"] / "run_manifest.json",
            "ensemble_manifest": runs["joint_ensemble"] / "run_manifest.json",
            "aroma_manifest": runs["aroma_calibration"] / "run_manifest.json",
            "actuator_manifest": runs["temperature_actuator"] / "run_manifest.json",
            "information_reference_manifest": runs["information_reference"]
            / "run_manifest.json",
            "sampling_reference_manifest": runs["sampling_reference"]
            / "run_manifest.json",
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "operational_coverage.py",
        ],
        random_seeds=[
            int(value) for value in coverage_config["search"]["independent_seeds"]
        ],
        status=(
            "completed"
            if gate_verdict(gate_checks, tuple(gate_checks)) == "PASS"
            else "failed"
        ),
        convergence=source_manifest["solver_status_and_convergence"],
        gate={
            "verdict": gate_verdict(gate_checks, tuple(gate_checks)),
            "checks": gate_checks,
            "physical_authorization_gate": "FAIL_NOT_AUTHORIZED",
        },
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    manifest.update(
        {
            "qualification_type": "immutable_label_integrity_repackage",
            "source_run_output_audit": source_audit,
            "physical_candidate_label_integrity_repair": repair,
            "numeric_results_reused": True,
            "watermark": WATERMARK,
            **closed_authorization(),
        }
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(
        json.dumps(
            {
                "run_dir": relative_or_absolute(run_dir),
                "verdict": manifest["gate"]["verdict"],
                "qualification_type": manifest["qualification_type"],
                **closed_authorization(),
            }
        ),
        flush=True,
    )
    return run_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run immutable operational-coverage MBDoE qualification"
    )
    parser.add_argument(
        "--result-root",
        default=str(RESULT_ROOT),
        help="Parent adaptive-design result root; a unique run directory is always created",
    )
    parser.add_argument(
        "--reuse-completed-run",
        help=(
            "Create a new immutable qualification run from a completed run, "
            "verifying hashes and restoring cached physical candidate labels"
        ),
    )
    return parser.parse_args()


def main() -> None:
    import copy

    args = parse_args()
    started = time.perf_counter()
    git_snapshot = capture_git_state()
    coverage_config = load_json(COVERAGE_CONFIG_PATH)
    if coverage_config["watermark"] != WATERMARK:
        raise ValueError("Coverage watermark must match the required exact text")
    if coverage_config["authorization"] != AUTHORIZATION_FLAGS:
        raise ValueError("Coverage configuration must remain fail closed")
    if args.reuse_completed_run:
        _repackage_completed_run(
            _run_path(args.reuse_completed_run),
            _run_path(args.result_root),
            coverage_config,
            git_snapshot,
        )
        return
    model_config = load_wave1_config(MODEL_CONFIG_PATH, CONSTRAINTS_PATH)
    model_config = copy.deepcopy(model_config)
    # Coverage uses a strict lexicographic objective. Complexity is a safety
    # tie-breaker, never a hidden penalty in the information level.
    model_config["objective"]["design_complexity_penalty"] = {
        "per_degree_total_variation": 0.0,
        "per_temperature_change": 0.0,
        "origin": "disabled_for_coverage_lexicographic_level_2",
    }
    model_config["coverage_sampling_loss_fraction"] = float(
        coverage_config["sampling"]["harmonized_maximum_score_loss_fraction"]
    )
    model_config["coverage_harmonized_minimum_operator_round_reduction"] = int(
        coverage_config["sampling"]["harmonized_minimum_operator_round_reduction"]
    )
    runs, source_verification = _source_contract(coverage_config)
    adapter_gate = _read_json(runs["adapter"] / "adapter_gate.json")
    reference_search_gate = _read_json(
        runs["information_reference"] / "final_search_gate.json"
    )
    reference_sampling_gate = _read_json(
        runs["sampling_reference"] / "final_sampling_gate.json"
    )
    if adapter_gate.get("verdict") != "PASS":
        raise RuntimeError("Corrected adapter source gate is not PASS")
    if reference_search_gate.get("verdict") not in {"PASS", "PASS_CONDITIONAL"}:
        raise RuntimeError("Information-reference search source gate is not usable")
    if reference_sampling_gate.get("verdict") not in {"PASS", "PASS_CONDITIONAL"}:
        raise RuntimeError("Sampling-reference source gate is not usable")
    result_root = _run_path(args.result_root)
    run_dir = create_immutable_run_directory(
        result_root, "wave1_operational_coverage", coverage_config
    )
    checkpoint_dir = run_dir / "checkpoints"
    filesystem_path(checkpoint_dir).mkdir(parents=True, exist_ok=False)
    outputs: list[Path] = []
    outputs_by_name: dict[str, Path] = {}

    def register(path: Path) -> None:
        outputs_by_name[path.name] = path

    config_output = run_dir / "coverage_design_config.json"
    _write_json(config_output, coverage_config, outputs)
    register(config_output)
    envelope = derive_robust_setpoint_envelope(model_config, coverage_config)
    envelope_path = run_dir / "robust_setpoint_envelope.json"
    _write_json(envelope_path, envelope, outputs)
    register(envelope_path)

    ensemble = pd.read_csv(
        filesystem_path(runs["joint_ensemble"] / "joint_parameter_ensemble.csv")
    )
    if len(ensemble) != 64:
        raise ValueError("Operational coverage requires exactly 64 ensemble members")
    aroma_config = load_json(AROMA_CONFIG_PATH)
    partitions, partition_provenance = load_partition_surrogates(
        aroma_config, REPOSITORY_DIR
    )
    prior = prior_precision(ensemble, model_config)
    reference_policies = _load_policies(
        runs["information_reference"] / "final_candidate_policy_actions.csv",
        model_config,
    )
    source_schedule = _reference_schedule(
        runs["sampling_reference"] / "optimized_sampling_schedule.csv"
    )
    source_comparison = pd.read_csv(
        filesystem_path(
            runs["sampling_reference"] / "sampling_candidate_comparison.csv"
        )
    )
    prepared_cache: dict[tuple[Any, ...], Any] = {}
    sampling_fim_cache: dict[tuple, np.ndarray] = {}
    reference_context = _reference_information_context(
        reference_policies=reference_policies,
        source_schedule=source_schedule,
        source_comparison=source_comparison,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        partitions=partitions,
        prepared_cache=prepared_cache,
        fim_cache=sampling_fim_cache,
    )
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]] = {}
    reference_physical = _reference_corrected_physical_record(
        reference_policies, model_config, coverage_config, physical_cache
    )

    design_cache: dict[tuple, tuple[np.ndarray, float, float]] = {}
    search_by_variant: dict[tuple[str, str], dict[str, Any]] = {}
    all_history = []
    all_local = []
    all_continuation = []
    all_seed_champions = []
    all_full_candidates = []
    checkpoint_paths: list[Path] = []
    for strategy in ("II_diagonal_1", "III_diagonal_2"):
        for dose_design in coverage_config["dose_designs"]:
            print(
                json.dumps(
                    {
                        "stage": "coverage_search",
                        "strategy": strategy,
                        "dose": dose_design,
                    }
                ),
                flush=True,
            )
            result = _run_strategy_variant_search(
                strategy=strategy,
                dose_design=dose_design,
                checkpoint_dir=checkpoint_dir,
                ensemble=ensemble,
                prior=prior,
                model_config=model_config,
                coverage_config=coverage_config,
                envelope=envelope,
                partitions=partitions,
                design_cache=design_cache,
                physical_cache=physical_cache,
            )
            search_by_variant[(strategy, dose_design)] = result
            all_history.extend(result["history_rows"])
            all_local.extend(result["local_rows"])
            all_continuation.extend(result["continuation_rows"])
            all_seed_champions.extend(
                _flat_search_record(row) for row in result["seed_champions"]
            )
            all_full_candidates.extend(
                _flat_search_record(row) for row in result["all_records"]
            )
            checkpoint_paths.extend(result["checkpoint_paths"])
    outputs.extend(checkpoint_paths)
    history_path = run_dir / "multifidelity_search_history.csv"
    _write_csv(history_path, pd.DataFrame(all_history), outputs)
    local_path = run_dir / "local_refinement_comparison.csv"
    _write_csv(local_path, pd.DataFrame(all_local), outputs)
    continuation_path = run_dir / "multifidelity_continuation_summary.csv"
    _write_csv(continuation_path, pd.DataFrame(all_continuation), outputs)
    champions_path = run_dir / "per_seed_champions_64_member.csv"
    _write_csv(champions_path, pd.DataFrame(all_seed_champions), outputs)
    candidates_path = run_dir / "originals_and_refinements_64_member.csv"
    _write_csv(candidates_path, pd.DataFrame(all_full_candidates), outputs)

    selected_policies = {
        key: tuple(value["selected"]["policies"])
        for key, value in search_by_variant.items()
    }
    sampling_results: list[dict[str, Any]] = []
    sampling_restarts = []
    drying_frames = []
    for offset, ((strategy, dose_design), policies) in enumerate(
        selected_policies.items(), start=1
    ):
        print(
            json.dumps(
                {"stage": "sampling", "strategy": strategy, "dose": dose_design}
            ),
            flush=True,
        )
        results, restarts, drying = _optimize_sampling_for_design(
            policies=policies,
            strategy=strategy,
            dose_design=dose_design,
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            partitions=partitions,
            prepared_cache=prepared_cache,
            fim_cache=sampling_fim_cache,
            current_reference_gains=reference_context["current"]["gains"],
            search_record=search_by_variant[(strategy, dose_design)]["selected"],
            seed_offset=offset * 100,
        )
        sampling_results.extend(results)
        sampling_restarts.append(restarts)
        drying_frames.append(drying.assign(strategy=strategy, dose_design=dose_design))
    recommendation, coverage_comparison = _recommendation(
        sampling_results,
        search_by_variant,
        float(reference_context["current"]["metrics"]["robust_score"]),
        float(reference_context["three_anchor"]["metrics"]["robust_score"]),
        coverage_config,
    )

    current_metrics = reference_context["current"]["metrics"]
    current_row = {
        "strategy": "I_reference_information",
        "dose_design": "N80",
        "sampling_mode": "independent_source_final",
        "robust_score": float(current_metrics["robust_score"]),
        "median": float(current_metrics["median"]),
        "q10": float(current_metrics["q10"]),
        "minimum": float(current_metrics["minimum"]),
        "tail_cvar": float(current_metrics["tail_cvar"]),
        "fraction_exceeding_information_reference": 1.0,
        "worst_paired_loss_vs_information_reference": 0.0,
        "completion_probability": float(
            reference_search_gate["full_ensemble"]["completion_probability"]
        ),
        "minimum_drying_margin_h": float(
            min(
                reference_context["drying"]
                .loc[
                    reference_context["drying"].policy.eq(policy.name),
                    "drying_time_h",
                ]
                .min()
                - max(
                    [
                        12.0 * int(index)
                        for index in np.flatnonzero(
                            np.abs(np.diff(np.asarray(policy.temperature_c)))
                            >= 1.0 - 1e-9
                        )
                        + 1
                    ]
                    + [float(time_h) for time_h, _ in policy.nutrition_mg_yan_l]
                    + [0.0]
                )
                for policy in reference_policies
            )
        ),
        "maximum_internal_jump_c": reference_physical["maximum_internal_jump_c"],
        "maximum_initial_jump_c": reference_physical["maximum_initial_jump_c"],
        "minimum_robust_physical_temperature_c": reference_physical[
            "minimum_robust_physical_temperature_c"
        ],
        "maximum_robust_physical_temperature_c": reference_physical[
            "maximum_robust_physical_temperature_c"
        ],
        "anchor_yan_mg_l": 80.0,
        "coverage_yan_mg_l": 80.0,
        "temperature_changes": reference_physical["temperature_changes"],
        "total_thermal_variation_c": reference_physical["total_thermal_variation_c"],
        "minimum_posterior_fim_eigenvalue": float(
            current_metrics["minimum_posterior_fim_eigenvalue"]
        ),
        "maximum_posterior_fim_condition_number": float(
            current_metrics["maximum_posterior_fim_condition_number"]
        ),
        "distinct_sampling_slots": len(
            {time_h for values in source_schedule.values() for time_h in values}
        ),
        "harmonized_score_loss_fraction": 0.0,
        "additional_shared_nonbasal_slots": 0,
        "operator_rounds": 30,
        "operator_round_reduction": 0,
        "sampling_gate": "PASS",
        "coverage_gate": "FAIL",
        "thermal_gate": "PASS" if reference_physical["thermal_pass"] else "FAIL",
        "nutrition_gate": "PASS",
        "information_gate": "PASS",
        "retained_information_improvement_fraction": 1.0,
        "information_cost_absolute_vs_current": 0.0,
        "information_cost_relative_vs_current": 0.0,
        "all_operational_gates_pass": False,
    }
    comparison = pd.concat(
        [pd.DataFrame([current_row]), coverage_comparison], ignore_index=True
    )
    current_score = float(current_row["robust_score"])
    comparison["absolute_difference_vs_information_reference"] = (
        comparison["robust_score"] - current_score
    )
    comparison["relative_difference_vs_information_reference"] = (
        comparison["robust_score"] / max(abs(current_score), 1e-12) - 1.0
    )
    for index, row in comparison.iterrows():
        yan_anchor = float(row["anchor_yan_mg_l"])
        yan_coverage = float(row["coverage_yan_mg_l"])
        anchor_mass = nutrition_margin_row(
            str(row["strategy"]),
            str(row["dose_design"]),
            _anchor_with_dose(model_config, yan_anchor, "mass_anchor"),
            model_config,
            coverage_config,
        )
        coverage_mass = nutrition_margin_row(
            str(row["strategy"]),
            str(row["dose_design"]),
            DesignPolicy("mass_coverage", tuple([18.0] * 14), ((46.0, yan_coverage),)),
            model_config,
            coverage_config,
        )
        comparison.at[index, "anchor_organic_product_g"] = anchor_mass[
            "organic_product_g"
        ]
        comparison.at[index, "anchor_dap_g"] = anchor_mass["dap_product_g"]
        comparison.at[index, "coverage_organic_product_g"] = coverage_mass[
            "organic_product_g"
        ]
        comparison.at[index, "coverage_dap_g"] = coverage_mass["dap_product_g"]

    strategy_comparison_path = run_dir / "coverage_strategy_comparison.csv"
    _write_csv(strategy_comparison_path, comparison, outputs)
    register(strategy_comparison_path)
    reference_comparison_path = run_dir / "coverage_vs_information_reference.csv"
    _write_csv(reference_comparison_path, comparison, outputs)
    register(reference_comparison_path)
    sampling_comparison_path = run_dir / "sampling_strategy_comparison.csv"
    _write_csv(
        sampling_comparison_path,
        pd.DataFrame([row["row"] for row in sampling_results]),
        outputs,
    )
    register(sampling_comparison_path)
    schedules = pd.concat(
        [row["schedule_frame"] for row in sampling_results], ignore_index=True
    )
    schedules_path = run_dir / "optimized_sampling_schedules.csv"
    _write_csv(schedules_path, schedules, outputs)
    register(schedules_path)
    captures = pd.concat(
        [row["captures"] for row in sampling_results], ignore_index=True
    )
    captures_path = run_dir / "capture_interval_summary.csv"
    _write_csv(captures_path, captures, outputs)
    register(captures_path)
    conflicts = pd.concat(
        [row["conflicts"] for row in sampling_results], ignore_index=True
    )
    conflicts_path = run_dir / "sampling_operational_conflicts.csv"
    _write_csv(conflicts_path, conflicts, outputs)
    restarts_path = run_dir / "sampling_search_restarts.csv"
    _write_csv(restarts_path, pd.concat(sampling_restarts, ignore_index=True), outputs)
    drying_path = run_dir / "drying_time_by_member_candidate.csv"
    _write_csv(drying_path, pd.concat(drying_frames, ignore_index=True), outputs)

    parameter_information, variance_ratios, correlation_summary, eigenvalues = (
        _posterior_information_tables(
            sampling_results,
            reference_context["current"],
            reference_context["three_anchor"],
            model_config,
        )
    )
    for name, frame in (
        ("parameter_information_comparison.csv", parameter_information),
        ("posterior_variance_ratio_by_strategy.csv", variance_ratios),
        ("posterior_correlation_summary.csv", correlation_summary),
        ("fim_eigenvalue_comparison.csv", eigenvalues),
    ):
        path = run_dir / name
        _write_csv(path, frame, outputs)
        register(path)

    action_frames = [
        _candidate_action_rows(
            "I_reference_information",
            "N80",
            "reference_source_final",
            reference_policies,
            model_config,
            coverage_config,
        )
    ]
    nutrition_rows = []
    physical_frames = []
    for policy in reference_policies:
        frame, _ = _physical_summary_cached(
            policy, model_config, coverage_config, physical_cache
        )
        physical_frames.append(
            frame.assign(
                candidate=policy.name,
                strategy="I_reference_information",
                dose_design="N80",
                candidate_id="reference_source_final",
            )
        )
        nutrition_rows.append(
            nutrition_margin_row(
                "I_reference_information",
                "N80",
                policy,
                model_config,
                coverage_config,
            )
        )
    for (strategy, dose_design), search_result in search_by_variant.items():
        selected = search_result["selected"]
        policies = tuple(selected["policies"])
        action_frames.append(
            _candidate_action_rows(
                strategy,
                dose_design,
                selected["candidate_id"],
                policies,
                model_config,
                coverage_config,
            )
        )
        for policy in policies:
            frame, _ = _physical_summary_cached(
                policy, model_config, coverage_config, physical_cache
            )
            physical_frames.append(
                frame.assign(
                    candidate=policy.name,
                    strategy=strategy,
                    dose_design=dose_design,
                    candidate_id=selected["candidate_id"],
                )
            )
            nutrition_rows.append(
                nutrition_margin_row(
                    strategy,
                    dose_design,
                    policy,
                    model_config,
                    coverage_config,
                )
            )
    actions = pd.concat(action_frames, ignore_index=True)
    actions_path = run_dir / "coverage_candidate_actions.csv"
    _write_csv(actions_path, actions, outputs)
    register(actions_path)
    physical = pd.concat(physical_frames, ignore_index=True)
    physical_path = run_dir / "physical_temperature_envelope_by_candidate.csv"
    _write_csv(physical_path, physical, outputs)
    register(physical_path)

    nutrition = pd.DataFrame(nutrition_rows)
    best_scores = (
        pd.DataFrame([row["row"] for row in sampling_results])
        .pivot_table(
            index=["strategy", "dose_design"],
            columns="sampling_mode",
            values="robust_score",
            aggfunc="max",
        )
        .reset_index()
    )
    nutrition = nutrition.merge(best_scores, on=["strategy", "dose_design"], how="left")
    for strategy in ("II_diagonal_1", "III_diagonal_2"):
        for mode in ("independent", "partially_harmonized"):
            baseline = best_scores[
                best_scores.strategy.eq(strategy) & best_scores.dose_design.eq("N80")
            ]
            if baseline.empty or mode not in baseline:
                continue
            score80 = float(baseline.iloc[0][mode])
            mask = nutrition.strategy.eq(strategy)
            nutrition.loc[mask, f"{mode}_loss_fraction_vs_N80"] = (
                score80 - nutrition.loc[mask, mode]
            ) / max(abs(score80), 1e-12)
            nutrition.loc[mask, f"{mode}_N76_preferred_if_loss_le_0_5pct"] = (
                nutrition.loc[mask, f"{mode}_loss_fraction_vs_N80"] <= 0.005 + 1e-12
            )
    nutrition_path = run_dir / "nutrition_margin_comparison.csv"
    _write_csv(nutrition_path, nutrition, outputs)
    register(nutrition_path)

    convergence_payload = {
        f"{strategy}__{dose}": result["convergence"]
        for (strategy, dose), result in search_by_variant.items()
    }
    convergence_path = run_dir / "practical_convergence_gate.json"
    _write_json(convergence_path, {"variants": convergence_payload}, outputs)
    fim_checks = {
        "source_adapter_and_manifests_verified": len(source_verification) == 11,
        "all_sampling_designs_use_64_members": all(
            len(row["posterior"]) == 64 for row in sampling_results
        ),
        "all_posterior_fims_finite": all(
            np.isfinite(np.stack(row["posterior"])).all() for row in sampling_results
        ),
        "all_posterior_fims_symmetric": all(
            max(np.linalg.norm(fim - fim.T, ord="fro") for fim in row["posterior"])
            <= 1e-8
            for row in sampling_results
        ),
        "all_posterior_fims_psd": all(
            min(np.linalg.eigvalsh(0.5 * (fim + fim.T))[0] for fim in row["posterior"])
            >= -1e-8
            for row in sampling_results
        ),
        "parameter_level_information_published": len(parameter_information) > 0,
        "tail_metrics_published": len(sampling_results) == 12,
        "originals_and_refinements_compared": len(all_local) == 60,
    }
    post_fim_path = run_dir / "post_search_fim_validation.json"
    _write_json(
        post_fim_path,
        {
            "verdict": gate_verdict(fim_checks, tuple(fim_checks)),
            "checks": fim_checks,
            "source_reconciliation": reference_context["source_reconciliation"],
        },
        outputs,
    )

    gate_variants = {}
    for key, result in search_by_variant.items():
        selected = result["selected"]
        sampling_pass = all(
            row["row"]["sampling_gate"] == "PASS"
            for row in sampling_results
            if (row["row"]["strategy"], row["row"]["dose_design"]) == key
        )
        checks = {
            "coverage": bool(selected["coverage_pass"]),
            "thermal": bool(selected["thermal_pass"]),
            "nutrition": bool(selected["nutrition_pass"]),
            "information_sources_and_fim": gate_verdict(fim_checks, tuple(fim_checks))
            == "PASS",
            "practical_convergence": bool(result["convergence"]["passed"]),
            "sampling_and_capture": sampling_pass,
            "physical_flags_closed": True,
        }
        gate_variants[f"{key[0]}__{key[1]}"] = {
            "verdict": gate_verdict(checks, tuple(checks)),
            "checks": checks,
            "candidate_id": selected["candidate_id"],
        }
    archetype_gate_path = run_dir / "coverage_archetype_gate.json"
    _write_json(
        archetype_gate_path,
        {
            "reference_revalidation": {
                "coverage_gate": "FAIL",
                "thermal_gate": (
                    "PASS" if reference_physical["thermal_pass"] else "FAIL"
                ),
                **reference_physical,
            },
            "coverage_variants": gate_variants,
            "physical_authorization_gate": "FAIL_NOT_AUTHORIZED",
        },
        outputs,
    )
    register(archetype_gate_path)

    recommendation_path = run_dir / "wave1_recommended_coverage_candidate.json"
    _write_json(recommendation_path, recommendation, outputs)
    register(recommendation_path)
    if recommendation.get("strategy"):
        complementary_strategy = (
            "III_diagonal_2"
            if recommendation["strategy"] == "II_diagonal_1"
            else "II_diagonal_1"
        )
        complementary_key = (complementary_strategy, recommendation["dose_design"])
    else:
        complementary_key = ("III_diagonal_2", "N76_historical_anchor")
    complementary = selected_policies[complementary_key]
    wave2_rows = []
    for policy in complementary[1:]:
        metric = coverage_policy_metrics(policy, model_config, coverage_config)
        wave2_rows.append(
            {
                "wave": 2,
                "cell": f"{metric.thermal_archetype}_{metric.nutrition_archetype}",
                "policy": policy.name,
                "temperature_profile_12h_blocks_c": json.dumps(
                    list(policy.temperature_c)
                ),
                "nutrition_time_h": metric.nutrition_time_h,
                "yan_mg_l": float(policy.nutrition_mg_yan_l[0][1]),
                "selection_after_wave1": (
                    "re-rank using posterior reduction, least identified direction, and model-data discrepancy"
                ),
            }
        )
    wave2_path = run_dir / "wave2_complementary_candidates.csv"
    _write_csv(wave2_path, pd.DataFrame(wave2_rows), outputs)
    register(wave2_path)
    wave3_path = run_dir / "wave3_adaptive_decision_rule.json"
    _write_json(
        wave3_path,
        {
            "wave": 3,
            "decision_inputs": [
                "posterior_variance_reduction_after_waves_1_and_2",
                "least_identified_parameter_direction",
                "model_data_discrepancy",
                "replicate_need",
                "aroma_and_sensory_coverage",
            ],
            "rule": [
                "replicate the condition with material model-data discrepancy or poor reproducibility",
                "otherwise select the feasible condition maximizing expected reduction along the least identified direction",
                "use aroma/sensory coverage as the final tie-breaker",
            ],
            "physical_instructions_generated": False,
        },
        outputs,
    )
    register(wave3_path)

    best_by_strategy: dict[str, tuple[DesignPolicy, ...]] = {}
    for strategy in ("II_diagonal_1", "III_diagonal_2"):
        candidates = coverage_comparison[
            coverage_comparison.strategy.eq(strategy)
            & coverage_comparison.all_operational_gates_pass
        ]
        if candidates.empty:
            candidates = coverage_comparison[coverage_comparison.strategy.eq(strategy)]
        chosen = candidates.sort_values("robust_score", ascending=False).iloc[0]
        best_by_strategy[strategy] = selected_policies[
            (strategy, str(chosen["dose_design"]))
        ]
    archetype_policies: dict[str, DesignPolicy] = {}
    for strategy, policies in best_by_strategy.items():
        for policy in policies[1:]:
            metric = coverage_policy_metrics(policy, model_config, coverage_config)
            archetype_policies[
                f"{metric.thermal_archetype}_{metric.nutrition_archetype}"
            ] = policy
    if set(archetype_policies) != {
        "cold_early",
        "cold_late",
        "warm_early",
        "warm_late",
    }:
        raise RuntimeError("Four operational archetypes are required for predictions")
    predictions = _ensemble_prediction_bands(
        archetype_policies, ensemble, model_config, partitions
    )
    predictions_path = run_dir / "operational_prediction_ensemble_bands.csv"
    _write_csv(predictions_path, predictions, outputs)
    physical_bands = _physical_trajectory_bands(archetype_policies, model_config)
    physical_bands_path = run_dir / "physical_temperature_trajectory_bands.csv"
    _write_csv(physical_bands_path, physical_bands, outputs)
    figure_recommendation = recommendation
    if not recommendation.get("strategy"):
        fallback = coverage_comparison.sort_values(
            "robust_score", ascending=False
        ).iloc[0]
        figure_recommendation = {
            "strategy": str(fallback["strategy"]),
            "dose_design": str(fallback["dose_design"]),
            "sampling_mode": str(fallback["sampling_mode"]),
        }
    figure_paths, figure_qa = _generate_figures(
        run_dir=run_dir,
        comparison=comparison,
        candidate_actions=actions,
        parameter_ratios=variance_ratios,
        predictions=predictions,
        physical_bands=physical_bands,
        recommended=figure_recommendation,
        selected_policies_by_strategy=best_by_strategy,
    )
    outputs.extend(figure_paths)
    figure_qa_path = run_dir / "visual_qa.json"
    _write_json(figure_qa_path, figure_qa, outputs)

    runtime_path = run_dir / "runtime_summary.json"
    _write_json(
        runtime_path,
        {
            "runtime_seconds": float(time.perf_counter() - started),
            "search_variants": len(search_by_variant),
            "independent_seeds_per_variant": len(
                coverage_config["search"]["independent_seeds"]
            ),
            "sampling_design_modes": len(sampling_results),
            "ensemble_members": len(ensemble),
            "approved_actuator_scenarios": envelope["approved_scenario_count"],
        },
        outputs,
    )
    register(runtime_path)

    if recommendation.get("recommended"):
        review_outputs = _build_review_package(
            run_dir,
            outputs_by_name,
            figure_paths,
            recommendation,
        )
        outputs.extend(review_outputs)

    run_checks = {
        "source_hashes_verified": len(source_verification) == 11,
        "immutable_run_root": run_dir.parent.name == "wave1_operational_coverage",
        "five_seeds_per_variant": all(
            len(result["seed_champions"]) == 5 for result in search_by_variant.values()
        ),
        "six_coverage_variants": len(search_by_variant) == 6,
        "twelve_sampling_comparisons": len(sampling_results) == 12,
        "figures_pass_structural_qa": figure_qa["verdict"] == "PASS",
        "physical_flags_closed": True,
        "no_tank_assignments": True,
    }
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_operational_coverage",
        config={"coverage": coverage_config, "model_overlay": model_config},
        sources={
            "coverage_config": COVERAGE_CONFIG_PATH,
            "model_config": MODEL_CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
            "adapter_manifest": runs["adapter"] / "run_manifest.json",
            "ensemble_manifest": runs["joint_ensemble"] / "run_manifest.json",
            "aroma_manifest": runs["aroma_calibration"] / "run_manifest.json",
            "actuator_manifest": runs["temperature_actuator"] / "run_manifest.json",
            "information_reference_manifest": runs["information_reference"]
            / "run_manifest.json",
            "sampling_reference_manifest": runs["sampling_reference"]
            / "run_manifest.json",
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "operational_coverage.py",
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "hybrid_optimizer.py",
            ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py",
        ],
        random_seeds=[
            int(value) for value in coverage_config["search"]["independent_seeds"]
        ],
        status=(
            "completed"
            if gate_verdict(run_checks, tuple(run_checks)) == "PASS"
            else "failed"
        ),
        convergence=convergence_payload,
        gate={
            "verdict": gate_verdict(run_checks, tuple(run_checks)),
            "checks": run_checks,
            "physical_authorization_gate": "FAIL_NOT_AUTHORIZED",
        },
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    manifest["source_output_verification"] = source_verification
    manifest["partition_surrogate_provenance"] = partition_provenance
    manifest["watermark"] = WATERMARK
    manifest.update(closed_authorization())
    manifest_path = run_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    print(
        json.dumps(
            {
                "run_dir": relative_or_absolute(run_dir),
                "verdict": manifest["gate"]["verdict"],
                "recommended": recommendation.get("recommended", False),
                **closed_authorization(),
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
