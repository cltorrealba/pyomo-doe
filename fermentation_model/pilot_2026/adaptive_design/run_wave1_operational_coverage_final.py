from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

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

from pilot_2026.adaptive_design.final_search_logic import (  # noqa: E402
    approved_actuator_scenarios,
    policy_distance,
)
from pilot_2026.adaptive_design.hybrid_optimizer import (  # noqa: E402
    CheckpointedSwarmState,
    advance_checkpointed_swarm,
    initialize_checkpointed_swarm,
    rescore_checkpointed_swarm,
)
from pilot_2026.adaptive_design.local_refinement import (  # noqa: E402
    sequential_ipopt_refine,
)
from pilot_2026.adaptive_design.operational_coverage import (  # noqa: E402
    AUTHORIZATION_FLAGS,
    REVIEW_SUFFIX,
    WATERMARK,
    coverage_policy_metrics,
    derive_robust_setpoint_envelope,
    drying_margin_by_policy_member,
    explicit_controller_blocks,
    feasible_information_policy_from_unit_vector,
    gate_verdict,
    lexicographic_select,
    nutrition_margin_row,
    validate_physical_temperature,
)
from pilot_2026.adaptive_design.optimize_wave1_sampling_and_plots import (  # noqa: E402
    _load_policies,
    _optimize_schedules,
    _score,
)
from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    design_fim_with_drying,
    evaluate_campaign,
    fim_from_prepared,
    load_json,
    load_wave1_config,
    parameter_columns,
    prepare_design,
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
from pilot_2026.adaptive_design.run_wave1_operational_coverage import (  # noqa: E402
    _anchor_with_dose,
    _candidate_action_rows,
    _decode_strategy_pair,
    _flat_search_record,
    _full_candidate_record,
    _manual_seed_positions,
    _optimize_sampling_for_design,
    _physical_summary_cached,
    _posterior_information_tables,
    _prepare_sampling_context,
    _reference_information_context,
    _reference_schedule,
    _source_contract,
)


FINAL_CONFIG_PATH = ADAPTIVE_DIR / "coverage_final_design_config.json"
COVERAGE_CONFIG_PATH = ADAPTIVE_DIR / "coverage_design_config.json"
MODEL_CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
BALANCED_STRATEGY = "II_diagonal_1"
BALANCED_DOSE = "N76_all_conservative"
BENCHMARK_STRATEGY = "I_feasible_information_only_N76"
DEPENDENCY_PATHS = (
    ADAPTIVE_DIR / "operational_coverage.py",
    ADAPTIVE_DIR / "run_wave1_operational_coverage.py",
    Path(__file__),
    ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
    ADAPTIVE_DIR / "hybrid_optimizer.py",
    ADAPTIVE_DIR / "local_refinement.py",
    ADAPTIVE_DIR / "final_search_logic.py",
    ADAPTIVE_DIR / "pilot_aroma_calibration.py",
    ADAPTIVE_DIR / "pilot_calibration.py",
    ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py",
    ADAPTIVE_DIR / "run_artifacts.py",
    FINAL_CONFIG_PATH,
    COVERAGE_CONFIG_PATH,
    MODEL_CONFIG_PATH,
    CONSTRAINTS_PATH,
    AROMA_CONFIG_PATH,
)
PALETTE = {
    "blue": "#235789",
    "gold": "#D4A72C",
    "orange": "#E07A3F",
    "olive": "#7A8B3A",
    "ink": "#1F2937",
    "grey": "#9CA3AF",
}


def _run_path(value: str | Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _json(path: Path) -> dict[str, Any]:
    return json.loads(filesystem_path(path).read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any], outputs: list[Path]) -> None:
    write_json(path, payload)
    outputs.append(path)


def _write_csv(path: Path, frame: pd.DataFrame, outputs: list[Path]) -> None:
    output = frame.copy()
    if "watermark" not in output:
        output["watermark"] = WATERMARK
    if "candidate_only_not_for_physical_execution" not in output:
        output["candidate_only_not_for_physical_execution"] = True
    output.to_csv(filesystem_path(path), index=False, lineterminator="\n")
    outputs.append(path)


def _dependency_manifest() -> dict[str, Any]:
    files = []
    for path in DEPENDENCY_PATHS:
        if not filesystem_path(path).is_file():
            raise FileNotFoundError(f"Scientific dependency is missing: {path}")
        files.append(
            {
                "path": relative_or_absolute(path),
                "sha256": sha256_file(path),
                "bytes": filesystem_path(path).stat().st_size,
            }
        )
    digest = sha256_payload(
        [{"path": row["path"], "sha256": row["sha256"]} for row in files]
    )
    return {
        "schema_version": 1,
        "scientific_dependency_sha256": digest,
        "dependencies": files,
        "transitive_dependency_count": len(files),
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }


def _model_config(final_config: dict[str, Any]) -> dict[str, Any]:
    model = copy.deepcopy(load_wave1_config(MODEL_CONFIG_PATH, CONSTRAINTS_PATH))
    model["objective"]["design_complexity_penalty"] = {
        "per_degree_total_variation": 0.0,
        "per_temperature_change": 0.0,
        "origin": "disabled_for_final_lexicographic_information_level",
    }
    owner = final_config["owner_selection"]
    model["anchor"] = {
        "temperature_c": 18.0,
        "nutrition_pulses_mg_yan_l": [
            [float(owner["anchor_nutrition_h"]), float(owner["yan_mg_l"])]
        ],
        "status": "operational_anchor_18C_N76_review_candidate",
    }
    required_margin = float(final_config["required_minimum_action_to_drying_margin_h"])
    model["operations"]["minimum_action_to_drying_margin_h"] = required_margin
    model["coverage_sampling_loss_fraction"] = float(
        final_config["sampling"]["harmonized_maximum_score_loss_fraction"]
    )
    model["coverage_harmonized_minimum_operator_round_reduction"] = 0
    return model


def _rename_final_policies(policies: tuple[DesignPolicy, ...]) -> tuple[DesignPolicy, ...]:
    names = (
        "operational_anchor_18C_N76",
        "cold_early_N76",
        "warm_late_N76",
    )
    return tuple(
        DesignPolicy(name, policy.temperature_c, policy.nutrition_mg_yan_l)
        for name, policy in zip(names, policies)
    )


def _policy_hash(policies: tuple[DesignPolicy, ...]) -> str:
    return sha256_payload(
        [
            {
                "name": policy.name,
                "temperature_c": list(policy.temperature_c),
                "nutrition_mg_yan_l": [list(row) for row in policy.nutrition_mg_yan_l],
            }
            for policy in policies
        ]
    )


def _candidate_id(policies: tuple[DesignPolicy, ...]) -> str:
    return _policy_hash(policies)[:16]


def _read_gzip(path: Path) -> dict[str, Any]:
    with gzip.open(filesystem_path(path), "rt", encoding="utf-8") as stream:
        return json.load(stream)


def _write_gzip(path: Path, payload: dict[str, Any]) -> None:
    filesystem_path(path.parent).mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(filesystem_path(temporary), "wt", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, sort_keys=True)
    filesystem_path(temporary).replace(filesystem_path(path))


def _warm_start_positions(final_config: dict[str, Any]) -> np.ndarray:
    positions = [row for row in _manual_seed_positions()]
    source = _run_path(final_config["balanced_search"]["warm_start_source_run"])
    checkpoint_dir = source / "checkpoints"
    latest: dict[int, tuple[int, np.ndarray]] = {}
    if filesystem_path(checkpoint_dir).is_dir():
        for path in filesystem_path(checkpoint_dir).glob(
            "II_diagonal_1__N76_all_conservative__seed_*.json.gz"
        ):
            payload = _read_gzip(Path(path))
            swarm = payload["swarm"]
            seed = int(swarm["seed"])
            sequence = int(payload["checkpoint_sequence"])
            vector = np.asarray(swarm["global_best_position"], dtype=float)
            if vector.shape == (24,) and (
                seed not in latest or sequence > latest[seed][0]
            ):
                latest[seed] = (sequence, vector)
    positions.extend(value[1] for _, value in sorted(latest.items()))
    unique: dict[str, np.ndarray] = {}
    for position in positions:
        key = sha256_payload(np.round(np.asarray(position), 12).tolist())
        unique.setdefault(key, np.asarray(position, dtype=float))
    return np.stack(list(unique.values()))


def _objective_factory(
    *,
    kind: str,
    members: list[int],
    fidelity: str,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    envelope: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    objective_cache: dict[str, float],
) -> Callable[[np.ndarray], float]:
    anchor = _anchor_with_dose(
        model_config, 76.0, "operational_anchor_18C_N76"
    )

    def decode(vector: np.ndarray) -> tuple[DesignPolicy, ...]:
        if kind == "balanced":
            return _decode_strategy_pair(
                vector,
                BALANCED_STRATEGY,
                BALANCED_DOSE,
                model_config,
                coverage_config,
                envelope,
            )
        return (
            anchor,
            feasible_information_policy_from_unit_vector(
                f"{BENCHMARK_STRATEGY}__A",
                76.0,
                vector[:12],
                model_config,
                coverage_config,
                envelope,
            ),
            feasible_information_policy_from_unit_vector(
                f"{BENCHMARK_STRATEGY}__B",
                76.0,
                vector[12:],
                model_config,
                coverage_config,
                envelope,
            ),
        )

    def objective(vector: np.ndarray) -> float:
        policies = decode(vector)
        key = f"{fidelity}|{_candidate_id(tuple(policies))}"
        if key not in objective_cache:
            score, _ = evaluate_campaign(
                tuple(policies),
                ensemble,
                members,
                prior,
                model_config,
                partitions,
                design_cache=design_cache,
            )
            objective_cache[key] = -float(score)
        return float(objective_cache[key])

    return objective


def _checkpoint_path(
    checkpoint_dir: Path, kind: str, seed: int, sequence: int, phase: str, iteration: int
) -> Path:
    prefix = "b" if kind == "balanced" else "i"
    safe_phase = phase.replace("_member", "m").replace("continuation", "cont")
    return checkpoint_dir / (
        f"{prefix}_s{seed}_{sequence:03d}_{safe_phase}_i{iteration:03d}.json.gz"
    )


def _save_checkpoint(
    *,
    checkpoint_dir: Path,
    kind: str,
    state: CheckpointedSwarmState,
    sequence: int,
    phase: str,
    dependency_hash: str,
    config_hash: str,
    objective_cache: dict[str, float],
    base_vector: np.ndarray | None,
    continuation_windows: list[dict[str, Any]],
    complete: bool,
) -> Path:
    path = _checkpoint_path(
        checkpoint_dir, kind, state.seed, sequence, phase, state.iteration
    )
    candidate_id = sha256_payload(np.asarray(state.global_best_position).tolist())[:16]
    _write_gzip(
        path,
        {
            "checkpoint_schema_version": 2,
            "search_kind": kind,
            "checkpoint_sequence": int(sequence),
            "phase": phase,
            "complete": bool(complete),
            "saved_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "scientific_dependency_sha256": dependency_hash,
            "final_config_sha256": config_hash,
            "swarm": state.to_payload(),
            "objective_cache": objective_cache,
            "base_eight_member_vector": (
                None if base_vector is None else np.asarray(base_vector).tolist()
            ),
            "candidate_by_seed": {
                str(state.seed): {
                    "candidate_id": candidate_id,
                    "vector": np.asarray(state.global_best_position).tolist(),
                }
            },
            "continuation_windows": continuation_windows,
            "watermark": WATERMARK,
            **AUTHORIZATION_FLAGS,
        },
    )
    return path


def _latest_valid_checkpoint(
    checkpoint_dir: Path,
    kind: str,
    seed: int,
    dependency_hash: str,
    config_hash: str,
) -> tuple[Path, dict[str, Any]] | None:
    prefix = "b" if kind == "balanced" else "i"
    candidates = []
    for path in filesystem_path(checkpoint_dir).glob(f"{prefix}_s{seed}_*.json.gz"):
        try:
            payload = _read_gzip(Path(path))
        except (OSError, json.JSONDecodeError, KeyError, ValueError):
            continue
        if (
            payload.get("scientific_dependency_sha256") == dependency_hash
            and payload.get("final_config_sha256") == config_hash
            and payload.get("search_kind") == kind
            and int(payload["swarm"]["seed"]) == int(seed)
        ):
            candidates.append((int(payload["checkpoint_sequence"]), Path(path), payload))
    if not candidates:
        return None
    _, path, payload = max(candidates, key=lambda row: row[0])
    return path, payload


def _run_swarm_seed(
    *,
    kind: str,
    seed: int,
    search_config: dict[str, Any],
    checkpoint_dir: Path,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    envelope: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    dependency_hash: str,
    config_hash: str,
    warm_positions: np.ndarray,
    stage_started: float,
) -> dict[str, Any]:
    representatives = representative_members(ensemble, model_config)
    four_members = representatives[:4]
    eight_members = representatives[:8]
    objective_cache: dict[str, float] = {}
    bounds = np.tile(np.asarray([[0.0, 1.0]], dtype=float), (24, 1))
    latest = _latest_valid_checkpoint(
        checkpoint_dir, kind, seed, dependency_hash, config_hash
    )
    base_vector: np.ndarray | None = None
    continuation_windows: list[dict[str, Any]] = []
    checkpoint_paths: list[Path] = []
    if latest is None:
        four = _objective_factory(
            kind=kind,
            members=four_members,
            fidelity="four_member",
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            envelope=envelope,
            partitions=partitions,
            design_cache=design_cache,
            objective_cache=objective_cache,
        )
        seeds = warm_positions[: int(search_config["particles"])]
        state = initialize_checkpointed_swarm(
            four,
            bounds,
            particles=int(search_config["particles"]),
            seed=seed,
            fidelity="four_member",
            initial_positions=seeds,
        )
        sequence = 1
        phase = "four_member"
        checkpoint_paths.append(
            _save_checkpoint(
                checkpoint_dir=checkpoint_dir,
                kind=kind,
                state=state,
                sequence=sequence,
                phase=phase,
                dependency_hash=dependency_hash,
                config_hash=config_hash,
                objective_cache=objective_cache,
                base_vector=None,
                continuation_windows=[],
                complete=False,
            )
        )
    else:
        _path, payload = latest
        state = CheckpointedSwarmState.from_payload(payload["swarm"])
        sequence = int(payload["checkpoint_sequence"])
        phase = str(payload["phase"])
        objective_cache.update(
            {str(key): float(value) for key, value in payload["objective_cache"].items()}
        )
        if payload.get("base_eight_member_vector") is not None:
            base_vector = np.asarray(payload["base_eight_member_vector"], dtype=float)
        continuation_windows = list(payload.get("continuation_windows", []))
        if payload.get("complete") is True:
            return {
                "state": state,
                "base_vector": base_vector,
                "continuation_windows": continuation_windows,
                "checkpoint_paths": [],
                "resumed_from": relative_or_absolute(_path),
                "objective_cache_entries": len(objective_cache),
            }

    four = _objective_factory(
        kind=kind,
        members=four_members,
        fidelity="four_member",
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        objective_cache=objective_cache,
    )
    eight = _objective_factory(
        kind=kind,
        members=eight_members,
        fidelity="eight_member",
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        objective_cache=objective_cache,
    )
    four_target = int(search_config["four_member_iterations"])
    eight_target = four_target + int(search_config["eight_member_iterations"])
    if state.fidelity == "four_member":
        while state.iteration < four_target:
            state = advance_checkpointed_swarm(
                state,
                four,
                bounds,
                inertia=float(model_config["search"]["inertia"]),
                cognitive=float(model_config["search"]["cognitive"]),
                social=float(model_config["search"]["social"]),
            )
            sequence += 1
            checkpoint_paths.append(
                _save_checkpoint(
                    checkpoint_dir=checkpoint_dir,
                    kind=kind,
                    state=state,
                    sequence=sequence,
                    phase="four_member",
                    dependency_hash=dependency_hash,
                    config_hash=config_hash,
                    objective_cache=objective_cache,
                    base_vector=None,
                    continuation_windows=continuation_windows,
                    complete=False,
                )
            )
        sequence += 1
        checkpoint_paths.append(
            _save_checkpoint(
                checkpoint_dir=checkpoint_dir,
                kind=kind,
                state=state,
                sequence=sequence,
                phase="before_fidelity_change",
                dependency_hash=dependency_hash,
                config_hash=config_hash,
                objective_cache=objective_cache,
                base_vector=None,
                continuation_windows=continuation_windows,
                complete=False,
            )
        )
        state = rescore_checkpointed_swarm(state, eight, fidelity="eight_member")
        sequence += 1
        checkpoint_paths.append(
            _save_checkpoint(
                checkpoint_dir=checkpoint_dir,
                kind=kind,
                state=state,
                sequence=sequence,
                phase="after_fidelity_change",
                dependency_hash=dependency_hash,
                config_hash=config_hash,
                objective_cache=objective_cache,
                base_vector=None,
                continuation_windows=continuation_windows,
                complete=False,
            )
        )
    while state.iteration < eight_target:
        state = advance_checkpointed_swarm(
            state,
            eight,
            bounds,
            inertia=float(model_config["search"]["inertia"]),
            cognitive=float(model_config["search"]["cognitive"]),
            social=float(model_config["search"]["social"]),
        )
        sequence += 1
        checkpoint_paths.append(
            _save_checkpoint(
                checkpoint_dir=checkpoint_dir,
                kind=kind,
                state=state,
                sequence=sequence,
                phase="eight_member",
                dependency_hash=dependency_hash,
                config_hash=config_hash,
                objective_cache=objective_cache,
                base_vector=None,
                continuation_windows=continuation_windows,
                complete=False,
            )
        )
    if base_vector is None:
        base_vector = state.global_best_position.copy()
    if kind == "balanced":
        window = int(search_config["continuation_window_iterations"])
        tolerance = float(
            search_config["continuation_improvement_tolerance_fraction"]
        )
        maximum_windows = int(search_config["maximum_continuation_windows"])
        while len(continuation_windows) < maximum_windows:
            before = -float(state.global_best_value)
            start_iteration = state.iteration
            for _ in range(window):
                state = advance_checkpointed_swarm(
                    state,
                    eight,
                    bounds,
                    inertia=float(model_config["search"]["inertia"]),
                    cognitive=float(model_config["search"]["cognitive"]),
                    social=float(model_config["search"]["social"]),
                )
                sequence += 1
                checkpoint_paths.append(
                    _save_checkpoint(
                        checkpoint_dir=checkpoint_dir,
                        kind=kind,
                        state=state,
                        sequence=sequence,
                        phase="continuation",
                        dependency_hash=dependency_hash,
                        config_hash=config_hash,
                        objective_cache=objective_cache,
                        base_vector=base_vector,
                        continuation_windows=continuation_windows,
                        complete=False,
                    )
                )
            after = -float(state.global_best_value)
            improvement = max(after - before, 0.0) / max(abs(before), 1e-12)
            continuation_windows.append(
                {
                    "seed": seed,
                    "window": len(continuation_windows) + 1,
                    "start_iteration": start_iteration,
                    "end_iteration": state.iteration,
                    "score_before": before,
                    "score_after": after,
                    "continuation_improvement_fraction": improvement,
                }
            )
            if improvement <= tolerance + 1e-12:
                break
            if time.perf_counter() - stage_started >= float(
                search_config["budget_hours"]
            ) * 3600.0:
                break
    sequence += 1
    checkpoint_paths.append(
        _save_checkpoint(
            checkpoint_dir=checkpoint_dir,
            kind=kind,
            state=state,
            sequence=sequence,
            phase="complete",
            dependency_hash=dependency_hash,
            config_hash=config_hash,
            objective_cache=objective_cache,
            base_vector=base_vector,
            continuation_windows=continuation_windows,
            complete=True,
        )
    )
    return {
        "state": state,
        "base_vector": base_vector,
        "continuation_windows": continuation_windows,
        "checkpoint_paths": checkpoint_paths,
        "resumed_from": None if latest is None else relative_or_absolute(latest[0]),
        "objective_cache_entries": len(objective_cache),
    }


def _benchmark_record(
    *,
    vector: np.ndarray,
    policies: tuple[DesignPolicy, ...],
    seed: int,
    source_type: str,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
) -> dict[str, Any]:
    score, evaluations = evaluate_campaign(
        policies,
        ensemble,
        list(range(len(ensemble))),
        prior,
        model_config,
        partitions,
        design_cache=design_cache,
    )
    gains = np.asarray([row.information_gain for row in evaluations], dtype=float)
    posterior = [prior + row.fim for row in evaluations]
    metrics = robust_information_metrics(gains, posterior, None, model_config)
    drying = drying_margin_by_policy_member(
        policies,
        ensemble,
        model_config,
        partitions,
        design_cache=design_cache,
    )
    physical = [
        _physical_summary_cached(policy, model_config, coverage_config, physical_cache)[1]
        for policy in policies
    ]
    policy_metrics = [
        coverage_policy_metrics(policy, model_config, coverage_config)
        for policy in policies
    ]
    nutrition = [
        nutrition_margin_row(
            BENCHMARK_STRATEGY,
            BALANCED_DOSE,
            policy,
            model_config,
            coverage_config,
        )
        for policy in policies
    ]
    controller_pass = all(
        len(explicit_controller_blocks(policy, coverage_config)) == 42
        for policy in policies
    )
    thermal_pass = bool(
        controller_pass
        and all(row["initial_jump_robust_pass"] for row in physical)
        and all(row["physical_temperature_robust_pass"] for row in physical)
        and all(metric.maximum_internal_jump_c <= 5.0 + 1e-9 for metric in policy_metrics)
    )
    nutrition_pass = bool(all(row["nutrition_feasible"] for row in nutrition))
    completion_probability = float(np.mean([row.completion for row in evaluations]))
    completion_pass = completion_probability >= float(
        model_config["completion"]["minimum_probability"]
    ) - 1e-12
    minimum_margin = float(drying["action_to_drying_margin_h"].min())
    required_margin = float(
        model_config["operations"]["minimum_action_to_drying_margin_h"]
    )
    drying_pass = minimum_margin >= required_margin - 1e-12
    candidate_id = _candidate_id(tuple(policies))
    return {
        "candidate_id": candidate_id,
        "parent_candidate_id": None,
        "strategy": BENCHMARK_STRATEGY,
        "dose_design": BALANCED_DOSE,
        "independent_seed": int(seed),
        "source_type": source_type,
        "robust_score": float(metrics["robust_score"]),
        "penalized_search_score": float(score),
        "median": float(metrics["median"]),
        "q10": float(metrics["q10"]),
        "minimum": float(metrics["minimum"]),
        "tail_cvar": float(metrics["tail_cvar"]),
        "completion_probability": completion_probability,
        "minimum_drying_margin_h": minimum_margin,
        "maximum_residual_sugar_g_l": float(
            max(row.residual_sugar_g_l for row in evaluations)
        ),
        "minimum_posterior_fim_eigenvalue": float(
            metrics["minimum_posterior_fim_eigenvalue"]
        ),
        "maximum_posterior_fim_condition_number": float(
            metrics["maximum_posterior_fim_condition_number"]
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
            float(row["minimum_physical_temperature_c"]) for row in physical
        ),
        "maximum_robust_physical_temperature_c": max(
            float(row["maximum_physical_temperature_c"]) for row in physical
        ),
        "total_thermal_variation_c": sum(
            metric.total_thermal_variation_c for metric in policy_metrics
        ),
        "temperature_changes": sum(metric.temperature_changes for metric in policy_metrics),
        "nutrition_operational_margin_fraction": min(
            float(row["organic_margin_fraction"]) for row in nutrition
        ),
        "coverage_pass": True,
        "coverage_required": False,
        "thermal_pass": thermal_pass,
        "nutrition_pass": nutrition_pass,
        "completion_pass": completion_pass,
        "drying_margin_pass": drying_pass,
        "controller_explicit_to_504h": controller_pass,
        "first_24h_mean_contrast_c": math.nan,
        "nutrition_time_contrast_h": math.nan,
        "feasible": bool(
            thermal_pass
            and nutrition_pass
            and completion_pass
            and drying_pass
        ),
        "vector": np.asarray(vector, dtype=float),
        "policies": policies,
        "evaluations": evaluations,
        "drying_rows": drying,
    }


def _balanced_practical_convergence(
    champions: list[dict[str, Any]],
    seed_results: list[dict[str, Any]],
    model_config: dict[str, Any],
    final_config: dict[str, Any],
) -> dict[str, Any]:
    search = final_config["balanced_search"]
    best = max(float(row["robust_score"]) for row in champions)
    tolerance = float(search["near_best_fraction"])
    near = [
        row
        for row in champions
        if best - float(row["robust_score"])
        <= tolerance * max(abs(best), 1e-12) + 1e-12
    ]
    distances = []
    close_seeds: set[int] = set()
    for left_index, left in enumerate(near):
        for right in near[left_index + 1 :]:
            distance = policy_distance(
                tuple(left["policies"][1:]),
                tuple(right["policies"][1:]),
                model_config,
            )
            row = {
                "left_seed": int(left["independent_seed"]),
                "right_seed": int(right["independent_seed"]),
                **distance,
            }
            distances.append(row)
            if float(distance["policy_distance"]) <= float(
                search["policy_family_distance_threshold"]
            ):
                close_seeds.update((row["left_seed"], row["right_seed"]))
    last_windows = [
        result["continuation_windows"][-1]
        for result in seed_results
        if result["continuation_windows"]
    ]
    maximum_improvement = max(
        (
            float(row["continuation_improvement_fraction"])
            for row in last_windows
        ),
        default=math.inf,
    )
    checks = {
        "five_seed_champions_available": len(champions) == 5,
        "at_least_three_seed_champions_within_0_5pct": len(near)
        >= int(search["minimum_near_best_seed_champions"]),
        "all_near_best_champions_feasible": bool(near)
        and all(row["feasible"] for row in near),
        "operational_policy_family_has_three_seeds": len(close_seeds) >= 3,
        "eight_member_continuation_completed": len(last_windows) == 5,
        "continuation_improvement_at_most_0_1pct": maximum_improvement
        <= float(search["continuation_improvement_tolerance_fraction"]) + 1e-12,
    }
    return {
        "verdict": gate_verdict(checks, tuple(checks)),
        "checks": checks,
        "best_full_ensemble_robust_score": best,
        "near_best_seeds": [int(row["independent_seed"]) for row in near],
        "maximum_continuation_improvement_fraction": maximum_improvement,
        "required_maximum_continuation_improvement_fraction": float(
            search["continuation_improvement_tolerance_fraction"]
        ),
        "latest_continuation_window_by_seed": last_windows,
        "pairwise_policy_distances": distances,
    }


def _balanced_search(
    *,
    run_dir: Path,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    final_config: dict[str, Any],
    envelope: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
    dependency_hash: str,
    config_hash: str,
) -> dict[str, Any]:
    checkpoint_dir = run_dir / "balanced_search_checkpoints"
    filesystem_path(checkpoint_dir).mkdir(parents=True, exist_ok=True)
    search = final_config["balanced_search"]
    warm_positions = _warm_start_positions(final_config)
    stage_started = time.perf_counter()
    seed_results = []
    all_checkpoint_paths: list[Path] = []
    for seed in (int(value) for value in search["seeds"]):
        print(json.dumps({"stage": "balanced_search", "seed": seed}), flush=True)
        result = _run_swarm_seed(
            kind="balanced",
            seed=seed,
            search_config=search,
            checkpoint_dir=checkpoint_dir,
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            envelope=envelope,
            partitions=partitions,
            design_cache=design_cache,
            dependency_hash=dependency_hash,
            config_hash=config_hash,
            warm_positions=warm_positions,
            stage_started=stage_started,
        )
        seed_results.append(result)
        all_checkpoint_paths.extend(result["checkpoint_paths"])
    records: list[dict[str, Any]] = []
    champions: list[dict[str, Any]] = []
    for seed, result in zip(search["seeds"], seed_results):
        vectors = (
            ("original_seed_champion", result["base_vector"]),
            ("continued_seed_champion", result["state"].global_best_position),
        )
        per_seed = []
        seen: set[str] = set()
        for source_type, vector in vectors:
            policies = _decode_strategy_pair(
                np.asarray(vector),
                BALANCED_STRATEGY,
                BALANCED_DOSE,
                model_config,
                coverage_config,
                envelope,
            )
            record = _full_candidate_record(
                vector=np.asarray(vector),
                policies=policies,
                strategy=BALANCED_STRATEGY,
                dose_design=BALANCED_DOSE,
                seed=int(seed),
                source_type=source_type,
                parent_candidate_id=None,
                ensemble=ensemble,
                prior=prior,
                model_config=model_config,
                coverage_config=coverage_config,
                partitions=partitions,
                design_cache=design_cache,
                physical_cache=physical_cache,
            )
            record["eligible_for_final_selection"] = False
            key = f"{source_type}|{record['candidate_id']}"
            if key not in seen:
                records.append(record)
                per_seed.append(record)
                seen.add(key)
        champion = max(
            (row for row in per_seed if row["feasible"]),
            key=lambda row: float(row["robust_score"]),
            default=max(per_seed, key=lambda row: float(row["robust_score"])),
        )
        champions.append(champion)
    convergence = _balanced_practical_convergence(
        champions, seed_results, model_config, final_config
    )
    unique_candidates: dict[str, dict[str, Any]] = {}
    for row in records:
        current = unique_candidates.get(row["candidate_id"])
        if current is None or float(row["robust_score"]) > float(current["robust_score"]):
            unique_candidates[row["candidate_id"]] = row
    full_ranked = sorted(
        unique_candidates.values(),
        key=lambda row: (not bool(row["feasible"]), -float(row["robust_score"])),
    )
    top_for_refinement = full_ranked[: int(search["local_refinement_top_k"])]
    local_trace: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    local_records: list[dict[str, Any]] = []
    continuous = np.asarray([0, 1, 4, 7, 10, 12, 13, 16, 19, 22], dtype=int)
    bounds = np.tile(np.asarray([[0.0, 1.0]], dtype=float), (24, 1))
    reduced_cache: dict[str, float] = {}
    full_cache: dict[str, float] = {}
    reduced = _objective_factory(
        kind="balanced",
        members=representative_members(ensemble, model_config)[:8],
        fidelity="local_reduced_eight_member",
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        objective_cache=reduced_cache,
    )
    full = _objective_factory(
        kind="balanced",
        members=list(range(len(ensemble))),
        fidelity="local_full_64_member",
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        objective_cache=full_cache,
    )
    refinement_dir = run_dir / "balanced_search_checkpoints" / "local_refinement"
    filesystem_path(refinement_dir).mkdir(parents=True, exist_ok=True)
    for rank, parent in enumerate(top_for_refinement, start=1):
        before_path = refinement_dir / f"rank_{rank:02d}_before_ipopt.json.gz"
        _write_gzip(
            before_path,
            {
                "phase": "before_local_refinement",
                "parent_candidate_id": parent["candidate_id"],
                "vector": np.asarray(parent["vector"]).tolist(),
                "scientific_dependency_sha256": dependency_hash,
                "final_config_sha256": config_hash,
                "watermark": WATERMARK,
                **AUTHORIZATION_FLAGS,
            },
        )
        all_checkpoint_paths.append(before_path)
        result = sequential_ipopt_refine(
            np.asarray(parent["vector"], dtype=float),
            bounds,
            reduced,
            continuous,
            model_config,
            candidate_id=str(parent["candidate_id"]),
            validation_objective=full,
        )
        for row in result.trace:
            local_trace.append(
                {
                    "refinement_rank": rank,
                    "parent_candidate_id": parent["candidate_id"],
                    **row,
                }
            )
        rejected.extend(
            {
                "refinement_rank": rank,
                "parent_candidate_id": parent["candidate_id"],
                **row,
            }
            for row in result.rejected_candidates
        )
        policies = _decode_strategy_pair(
            result.x,
            BALANCED_STRATEGY,
            BALANCED_DOSE,
            model_config,
            coverage_config,
            envelope,
        )
        record = _full_candidate_record(
            vector=result.x,
            policies=policies,
            strategy=BALANCED_STRATEGY,
            dose_design=BALANCED_DOSE,
            seed=int(parent["independent_seed"]),
            source_type="ipopt_local_refinement",
            parent_candidate_id=parent["candidate_id"],
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            partitions=partitions,
            design_cache=design_cache,
            physical_cache=physical_cache,
        )
        record["local_refinement_qualified"] = bool(result.qualified)
        record["local_refinement_accepted"] = bool(result.accepted_improvement)
        record["local_refinement_state"] = result.state
        record["eligible_for_final_selection"] = bool(result.qualified and record["feasible"])
        local_records.append(record)
        after_path = refinement_dir / f"rank_{rank:02d}_after_ipopt.json.gz"
        _write_gzip(
            after_path,
            {
                "phase": "after_local_refinement",
                "parent_candidate_id": parent["candidate_id"],
                "result_candidate_id": record["candidate_id"],
                "vector": np.asarray(result.x).tolist(),
                "state": result.state,
                "qualified": bool(result.qualified),
                "accepted_improvement": bool(result.accepted_improvement),
                "trace": list(result.trace),
                "scientific_dependency_sha256": dependency_hash,
                "final_config_sha256": config_hash,
                "watermark": WATERMARK,
                **AUTHORIZATION_FLAGS,
            },
        )
        all_checkpoint_paths.append(after_path)
    records.extend(local_records)
    eligible = [row for row in records if row.get("eligible_for_final_selection") is True]
    if not eligible:
        raise RuntimeError("No IPOPT-qualified balanced candidate is eligible")
    selected, selection = lexicographic_select(
        eligible, float(search["near_best_fraction"])
    )
    convergence["checks"]["ipopt_executed_for_top_three"] = bool(
        len(local_records) == 3
        and all(row.get("local_refinement_state") != "not_executed" for row in local_records)
    )
    convergence["checks"]["selected_candidate_locally_qualified"] = bool(
        selected.get("local_refinement_qualified") is True
    )
    convergence["verdict"] = gate_verdict(
        convergence["checks"], tuple(convergence["checks"])
    )
    return {
        "selected": selected,
        "selection": selection,
        "records": records,
        "champions": champions,
        "seed_results": seed_results,
        "convergence": convergence,
        "local_trace": pd.DataFrame(local_trace),
        "rejected": rejected,
        "checkpoint_paths": all_checkpoint_paths,
        "runtime_seconds": time.perf_counter() - stage_started,
    }


def _benchmark_search(
    *,
    run_dir: Path,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    final_config: dict[str, Any],
    envelope: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]],
    dependency_hash: str,
    config_hash: str,
) -> dict[str, Any]:
    search = final_config["feasible_information_benchmark"]
    checkpoint_dir = run_dir / "feasible_information_search_checkpoints"
    filesystem_path(checkpoint_dir).mkdir(parents=True, exist_ok=True)
    stage_started = time.perf_counter()
    warm = _warm_start_positions(final_config)
    seed_results = []
    checkpoint_paths = []
    for seed in (int(value) for value in search["seeds"]):
        print(json.dumps({"stage": "feasible_information_benchmark", "seed": seed}), flush=True)
        result = _run_swarm_seed(
            kind="benchmark",
            seed=seed,
            search_config=search,
            checkpoint_dir=checkpoint_dir,
            ensemble=ensemble,
            prior=prior,
            model_config=model_config,
            coverage_config=coverage_config,
            envelope=envelope,
            partitions=partitions,
            design_cache=design_cache,
            dependency_hash=dependency_hash,
            config_hash=config_hash,
            warm_positions=warm,
            stage_started=stage_started,
        )
        seed_results.append(result)
        checkpoint_paths.extend(result["checkpoint_paths"])
    anchor = _anchor_with_dose(model_config, 76.0, "operational_anchor_18C_N76")
    records = []
    for seed, result in zip(search["seeds"], seed_results):
        vector = np.asarray(result["state"].global_best_position, dtype=float)
        policies = (
            anchor,
            feasible_information_policy_from_unit_vector(
                f"{BENCHMARK_STRATEGY}__A",
                76.0,
                vector[:12],
                model_config,
                coverage_config,
                envelope,
            ),
            feasible_information_policy_from_unit_vector(
                f"{BENCHMARK_STRATEGY}__B",
                76.0,
                vector[12:],
                model_config,
                coverage_config,
                envelope,
            ),
        )
        records.append(
            _benchmark_record(
                vector=vector,
                policies=policies,
                seed=int(seed),
                source_type="benchmark_seed_champion",
                ensemble=ensemble,
                prior=prior,
                model_config=model_config,
                coverage_config=coverage_config,
                partitions=partitions,
                design_cache=design_cache,
                physical_cache=physical_cache,
            )
        )
    feasible = [row for row in records if row["feasible"]]
    pool = feasible or records
    selected = max(pool, key=lambda row: float(row["robust_score"]))
    elapsed = time.perf_counter() - stage_started
    budget_exhausted = elapsed >= float(search["budget_hours"]) * 3600.0
    label = (
        str(search["timeout_label"])
        if budget_exhausted
        else "feasible_information_benchmark"
    )
    return {
        "selected": selected,
        "records": records,
        "seed_results": seed_results,
        "checkpoint_paths": checkpoint_paths,
        "benchmark_status": label,
        "budget_exhausted": budget_exhausted,
        "runtime_seconds": elapsed,
    }


def _select_sampling_result(
    results: list[dict[str, Any]], *, prefer_harmonized: bool
) -> dict[str, Any]:
    independent = next(row for row in results if row["row"]["sampling_mode"] == "independent")
    harmonized = next(
        row for row in results if row["row"]["sampling_mode"] == "partially_harmonized"
    )
    if prefer_harmonized:
        loss = float(harmonized["row"]["harmonized_score_loss_fraction"])
        if loss <= 0.005 + 1e-12:
            return harmonized
    return independent


def _sampling_and_information(
    *,
    balanced: dict[str, Any],
    benchmark: dict[str, Any],
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    model_config: dict[str, Any],
    coverage_config: dict[str, Any],
    partitions: dict[str, Any],
) -> dict[str, Any]:
    prepared_cache: dict[tuple[Any, ...], Any] = {}
    fim_cache: dict[tuple, np.ndarray] = {}
    benchmark_policies = tuple(benchmark["selected"]["policies"])
    benchmark_reference_gains = np.asarray(
        [row.information_gain for row in benchmark["selected"]["evaluations"]],
        dtype=float,
    )
    benchmark_results, benchmark_restarts, _ = _optimize_sampling_for_design(
        policies=benchmark_policies,
        strategy=BENCHMARK_STRATEGY,
        dose_design=BALANCED_DOSE,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        partitions=partitions,
        prepared_cache=prepared_cache,
        fim_cache=fim_cache,
        current_reference_gains=benchmark_reference_gains,
        search_record=benchmark["selected"],
        seed_offset=5000,
    )
    benchmark_selected = _select_sampling_result(
        benchmark_results, prefer_harmonized=False
    )
    balanced_policies = _rename_final_policies(tuple(balanced["selected"]["policies"]))
    balanced_record = dict(balanced["selected"])
    balanced_record["policies"] = balanced_policies
    balanced_results, balanced_restarts, _ = _optimize_sampling_for_design(
        policies=balanced_policies,
        strategy=BALANCED_STRATEGY,
        dose_design=BALANCED_DOSE,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        partitions=partitions,
        prepared_cache=prepared_cache,
        fim_cache=fim_cache,
        current_reference_gains=np.asarray(benchmark_selected["gains"], dtype=float),
        search_record=balanced_record,
        seed_offset=7000,
    )
    balanced_selected = _select_sampling_result(
        balanced_results, prefer_harmonized=True
    )
    anchor = balanced_policies[0]
    anchors = tuple(
        DesignPolicy(
            f"three_anchor_N76_{index + 1}",
            anchor.temperature_c,
            anchor.nutrition_mg_yan_l,
        )
        for index in range(3)
    )
    sampling_cache, valid, _drying = _prepare_sampling_context(
        anchors, ensemble, model_config, partitions, prepared_cache
    )
    for member in range(len(ensemble)):
        source = sampling_cache[(member, anchors[0].name)]
        for policy in anchors[1:]:
            sampling_cache[(member, policy.name)] = source
            valid[policy.name] = valid[anchors[0].name]
    anchor_schedules, anchor_restarts = _optimize_schedules(
        anchors,
        valid,
        sampling_cache,
        list(range(len(ensemble))),
        prior,
        model_config,
        fim_cache,
        seed_offset=9000,
    )
    anchor_score, anchor_gains, anchor_posterior = _score(
        anchor_schedules,
        anchors,
        sampling_cache,
        list(range(len(ensemble))),
        prior,
        model_config,
        fim_cache,
    )
    anchor_metrics = robust_information_metrics(
        anchor_gains, anchor_posterior, None, model_config
    )
    parameter, ratios, correlations, eigen = _posterior_information_tables(
        [balanced_selected],
        {
            "posterior": benchmark_selected["posterior"],
            "metrics": robust_information_metrics(
                benchmark_selected["gains"],
                benchmark_selected["posterior"],
                None,
                model_config,
            ),
        },
        {
            "posterior": anchor_posterior,
            "metrics": anchor_metrics,
        },
        model_config,
    )
    parameter.loc[
        parameter.strategy.eq("I_reference_information"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = [BENCHMARK_STRATEGY, BALANCED_DOSE, benchmark_selected["row"]["sampling_mode"]]
    parameter.loc[
        parameter.strategy.eq("three_anchor_reference"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = ["three_anchor_N76", BALANCED_DOSE, "independent_optimized"]
    ratios.loc[
        ratios.strategy.eq("I_reference_information"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = [BENCHMARK_STRATEGY, BALANCED_DOSE, benchmark_selected["row"]["sampling_mode"]]
    ratios.loc[
        ratios.strategy.eq("three_anchor_reference"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = ["three_anchor_N76", BALANCED_DOSE, "independent_optimized"]
    ratios["variance_ratio_vs_feasible_information_benchmark"] = ratios[
        "variance_ratio_vs_information_reference"
    ]
    correlations.loc[
        correlations.strategy.eq("I_reference_information"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = [BENCHMARK_STRATEGY, BALANCED_DOSE, benchmark_selected["row"]["sampling_mode"]]
    correlations.loc[
        correlations.strategy.eq("three_anchor_reference"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = ["three_anchor_N76", BALANCED_DOSE, "independent_optimized"]
    eigen.loc[
        eigen.strategy.eq("I_reference_information"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = [BENCHMARK_STRATEGY, BALANCED_DOSE, benchmark_selected["row"]["sampling_mode"]]
    eigen.loc[
        eigen.strategy.eq("three_anchor_reference"),
        ["strategy", "dose_design", "sampling_mode"],
    ] = ["three_anchor_N76", BALANCED_DOSE, "independent_optimized"]
    return {
        "balanced_policies": balanced_policies,
        "balanced_results": balanced_results,
        "balanced_selected": balanced_selected,
        "benchmark_results": benchmark_results,
        "benchmark_selected": benchmark_selected,
        "all_results": [*balanced_results, *benchmark_results],
        "restarts": pd.concat(
            [balanced_restarts, benchmark_restarts, anchor_restarts],
            ignore_index=True,
        ),
        "three_anchor": {
            "policies": anchors,
            "schedules": anchor_schedules,
            "score": float(anchor_score),
            "gains": anchor_gains,
            "posterior": anchor_posterior,
            "metrics": anchor_metrics,
        },
        "parameter_information": parameter,
        "variance_ratios": ratios,
        "correlations": correlations,
        "eigenvalues": eigen,
        "prepared_cache": prepared_cache,
        "fim_cache": fim_cache,
    }


def _relative_fim_difference(left: np.ndarray, right: np.ndarray) -> float:
    return float(
        np.linalg.norm(left - right, ord="fro")
        / max(np.linalg.norm(right, ord="fro"), 1e-12)
    )


def _post_search_fim_validation(
    *,
    policies: tuple[DesignPolicy, ...],
    schedules: dict[str, tuple[float, ...]],
    ensemble: pd.DataFrame,
    evaluations: list[Any],
    drying: pd.DataFrame,
    model_config: dict[str, Any],
    final_config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
) -> dict[str, Any]:
    representative = representative_members(ensemble, model_config)
    gains = np.asarray([row.information_gain for row in evaluations], dtype=float)
    member_ids = np.asarray([int(row.member) for row in evaluations], dtype=int)
    critical_drying = int(
        drying.sort_values("action_to_drying_margin_h").iloc[0]["ensemble_member"]
    )
    cases = {
        "central_member": int(representative[0]),
        "lowest_information_member": int(member_ids[np.argmin(gains)]),
        "highest_information_member": int(member_ids[np.argmax(gains)]),
        "worst_gain_member": int(member_ids[np.argmin(gains)]),
        "critical_drying_member": critical_drying,
        "extreme_actuator_member": int(representative[0]),
        "selected_candidate_exact": int(representative[0]),
    }
    scenarios = approved_actuator_scenarios(model_config)
    extreme = max(
        scenarios,
        key=lambda row: (
            abs(float(row["tracking_error_c"])),
            abs(float(row["probe_bias_c"])),
            abs(float(row["initial_temperature_offset_c"])),
            float(row["tau_h"]),
        ),
    )
    case_rows: list[dict[str, Any]] = []
    fd_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    direct_rows: list[dict[str, Any]] = []
    fd_steps = [float(value) for value in final_config["fim_validation"]["finite_difference_log_steps"]]
    grid_steps = [float(value) for value in final_config["fim_validation"]["time_grid_steps_h"]]
    fd_threshold = float(model_config["fim_validation"]["maximum_relative_step_fim_difference"])
    grid_threshold = float(model_config["fim_validation"]["maximum_relative_grid_fim_difference"])
    psd_tolerance = float(model_config["fim_validation"]["minimum_eigenvalue_tolerance"])
    for policy in policies:
        sample_times = np.asarray(schedules[policy.name], dtype=float)
        for case_label, member_index in cases.items():
            scenario = extreme if case_label == "extreme_actuator_member" else None
            member = ensemble.iloc[member_index]
            prepared_by_grid = {
                grid: prepare_design(
                    policy,
                    member,
                    model_config,
                    partitions,
                    simulation_grid_step_h=grid,
                    actuator_scenario=scenario,
                )
                for grid in grid_steps
            }
            if any(value is None for value in prepared_by_grid.values()):
                case_rows.append(
                    {
                        "policy": policy.name,
                        "member_case": case_label,
                        "ensemble_member": member_index,
                        "actuator_extreme": scenario is not None,
                        "passed": False,
                        "reason": "simulation_failed",
                    }
                )
                continue
            nominal_prepared = prepared_by_grid[2.0]
            nominal_fim = fim_from_prepared(
                nominal_prepared,
                sample_times,
                model_config,
                finite_difference_log_step=0.02,
            )
            fd_differences = []
            for step in fd_steps:
                fim = fim_from_prepared(
                    nominal_prepared,
                    sample_times,
                    model_config,
                    finite_difference_log_step=step,
                )
                difference = _relative_fim_difference(fim, nominal_fim)
                fd_differences.append(difference)
                fd_rows.append(
                    {
                        "policy": policy.name,
                        "member_case": case_label,
                        "ensemble_member": member_index,
                        "finite_difference_log_step": step,
                        "relative_fim_difference_vs_0_02": difference,
                    }
                )
            grid_differences = []
            for grid, prepared in prepared_by_grid.items():
                fim = fim_from_prepared(
                    prepared,
                    sample_times,
                    model_config,
                    finite_difference_log_step=0.02,
                )
                difference = _relative_fim_difference(fim, nominal_fim)
                grid_differences.append(difference)
                grid_rows.append(
                    {
                        "policy": policy.name,
                        "member_case": case_label,
                        "ensemble_member": member_index,
                        "time_grid_step_h": grid,
                        "relative_fim_difference_vs_2h": difference,
                    }
                )
            direct, _residual, _drying = design_fim_with_drying(
                policy,
                member,
                model_config,
                partitions,
                actuator_scenario=scenario,
            )
            cache_key = (
                int(member.get("ensemble_member", member_index)),
                tuple(policy.temperature_c),
                tuple(policy.nutrition_mg_yan_l),
                tuple(sorted((scenario or {}).items())),
            )
            cached = design_cache.get(cache_key)
            if cached is None:
                design_cache[cache_key] = (direct, float(_residual), float(_drying))
                cached_fim = direct
                cache_origin = "populated_from_direct_validation"
            else:
                cached_fim = cached[0]
                cache_origin = "search_or_sampling_cache"
            direct_difference = _relative_fim_difference(direct, cached_fim)
            direct_rows.append(
                {
                    "policy": policy.name,
                    "member_case": case_label,
                    "ensemble_member": member_index,
                    "cache_origin": cache_origin,
                    "relative_direct_vs_cached_fim_difference": direct_difference,
                }
            )
            symmetric = bool(np.allclose(nominal_fim, nominal_fim.T, atol=1e-8, rtol=0.0))
            eigenvalues = np.linalg.eigvalsh(0.5 * (nominal_fim + nominal_fim.T))
            positive = eigenvalues[eigenvalues > 1e-12]
            condition = float(eigenvalues[-1] / positive[0]) if len(positive) else math.inf
            checks = {
                "finite": bool(np.isfinite(nominal_fim).all()),
                "symmetric": symmetric,
                "psd": bool(eigenvalues[0] >= psd_tolerance),
                "fd_stable": max(fd_differences) <= fd_threshold + 1e-12,
                "grid_stable": max(grid_differences) <= grid_threshold + 1e-12,
                "direct_matches_cache": direct_difference <= 1e-10,
            }
            case_rows.append(
                {
                    "policy": policy.name,
                    "member_case": case_label,
                    "ensemble_member": member_index,
                    "actuator_extreme": scenario is not None,
                    "minimum_eigenvalue": float(eigenvalues[0]),
                    "maximum_eigenvalue": float(eigenvalues[-1]),
                    "condition_number": condition,
                    "maximum_fd_relative_difference": max(fd_differences),
                    "maximum_grid_relative_difference": max(grid_differences),
                    "direct_vs_cached_relative_difference": direct_difference,
                    **checks,
                    "passed": all(checks.values()),
                    "reason": None,
                }
            )
    case_frame = pd.DataFrame(case_rows)
    fd_frame = pd.DataFrame(fd_rows)
    grid_frame = pd.DataFrame(grid_rows)
    direct_frame = pd.DataFrame(direct_rows)
    checks = {
        "exact_selected_policies_covered": set(case_frame["policy"])
        == {policy.name for policy in policies},
        "all_required_member_cases_covered": set(case_frame["member_case"])
        == set(final_config["fim_validation"]["member_cases"]),
        "fd_steps_0_01_0_02_0_04_covered": set(fd_frame["finite_difference_log_step"])
        == {0.01, 0.02, 0.04},
        "grids_4_2_1h_covered": set(grid_frame["time_grid_step_h"])
        == {1.0, 2.0, 4.0},
        "all_fims_finite_symmetric_psd": bool(
            case_frame[["finite", "symmetric", "psd"]].all(axis=None)
        ),
        "all_fd_stability_checks_pass": bool(case_frame["fd_stable"].all()),
        "all_grid_stability_checks_pass": bool(case_frame["grid_stable"].all()),
        "all_direct_vs_cached_checks_pass": bool(
            case_frame["direct_matches_cache"].all()
        ),
        "all_cases_pass": bool(case_frame["passed"].all()),
    }
    gate = {
        "verdict": gate_verdict(checks, tuple(checks)),
        "checks": checks,
        "thresholds": {
            "maximum_relative_step_fim_difference": fd_threshold,
            "maximum_relative_grid_fim_difference": grid_threshold,
            "minimum_eigenvalue_tolerance": psd_tolerance,
            "maximum_direct_vs_cached_relative_difference": 1e-10,
        },
        "extreme_actuator_scenario": extreme,
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }
    return {
        "gate": gate,
        "cases": case_frame,
        "fd": fd_frame,
        "grid": grid_frame,
        "direct": direct_frame,
    }


def _tag_frame(frame: pd.DataFrame, state: dict[str, Any]) -> pd.DataFrame:
    tagged = frame.copy()
    tagged["qualification_candidate_id"] = state["candidate_id"]
    tagged["canonical_policy_hash"] = state["canonical_policy_hash"]
    tagged["final_verdict"] = state["final_verdict"]
    for gate in ("coverage", "thermal", "information", "nutrition", "sampling"):
        tagged[f"qualification_{gate}_gate"] = state["gates"][gate]["verdict"]
    return tagged


def _artifact_payload(
    state: dict[str, Any], payload: dict[str, Any]
) -> dict[str, Any]:
    return {
        **payload,
        "candidate_id": state["candidate_id"],
        "canonical_policy_hash": state["canonical_policy_hash"],
        "final_verdict": state["final_verdict"],
        "gate_verdicts": {
            key: value["verdict"]
            for key, value in state["gates"].items()
            if isinstance(value, dict) and "verdict" in value
        },
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }


def _final_state(
    *,
    balanced: dict[str, Any],
    benchmark: dict[str, Any],
    sampling: dict[str, Any],
    drying: pd.DataFrame,
    physical: pd.DataFrame,
    fim_validation: dict[str, Any],
    final_config: dict[str, Any],
    model_config: dict[str, Any],
    historical_metrics: dict[str, Any],
) -> dict[str, Any]:
    selected = balanced["selected"]
    policies = sampling["balanced_policies"]
    sampling_selected = sampling["balanced_selected"]
    sampling_row = sampling_selected["row"]
    benchmark_row = sampling["benchmark_selected"]["row"]
    three_anchor = sampling["three_anchor"]
    policy_metrics = [
        coverage_policy_metrics(policy, model_config, load_json(COVERAGE_CONFIG_PATH))
        for policy in policies
    ]
    cold, warm = policy_metrics[1], policy_metrics[2]
    coverage_checks = {
        "cold_early_classified_correctly": cold.thermal_archetype == "cold"
        and cold.nutrition_archetype == "early",
        "warm_late_classified_correctly": warm.thermal_archetype == "warm"
        and warm.nutrition_archetype == "late",
        "thermal_contrast_at_least_4c": warm.first_24h_mean_setpoint_c
        - cold.first_24h_mean_setpoint_c
        >= 4.0 - 1e-9,
        "nutrition_contrast_at_least_36h": warm.nutrition_time_h
        - cold.nutrition_time_h
        >= 36.0 - 1e-9,
        "treatments_not_equivalent": bool(selected["coverage_pass"]),
    }
    final_physical = physical[physical["strategy"].eq(BALANCED_STRATEGY)]
    controller = [explicit_controller_blocks(policy, load_json(COVERAGE_CONFIG_PATH)) for policy in policies]
    thermal_checks = {
        "maximum_robust_initial_jump_at_most_5c": float(
            final_physical["initial_setpoint_jump_c"].max()
        )
        <= 5.0 + 1e-9,
        "maximum_internal_jump_at_most_5c": max(
            metric.maximum_internal_jump_c for metric in policy_metrics
        )
        <= 5.0 + 1e-9,
        "all_effective_changes_at_least_1c": all(
            np.all(
                np.abs(np.diff(np.asarray(policy.temperature_c)))[
                    np.abs(np.diff(np.asarray(policy.temperature_c))) > 1e-12
                ]
                >= 1.0 - 1e-9
            )
            for policy in policies
        ),
        "physical_temperature_within_15_27c": bool(
            final_physical["physical_temperature_within_15_27c"].all()
        ),
        "exactly_243_scenarios_per_policy": bool(
            final_physical.groupby("candidate").size().eq(243).all()
        ),
        "controller_explicit_to_504h": all(
            len(frame) == 42
            and math.isclose(float(frame["end_h"].max()), 504.0, abs_tol=1e-9)
            and bool(frame.iloc[-1]["hold_last_setpoint"])
            for frame in controller
        ),
    }
    nutrition_rows = [
        nutrition_margin_row(
            BALANCED_STRATEGY,
            BALANCED_DOSE,
            policy,
            model_config,
            load_json(COVERAGE_CONFIG_PATH),
        )
        for policy in policies
    ]
    owner = final_config["owner_selection"]
    nutrition_checks = {
        "n76_in_all_three_tanks": all(
            math.isclose(
                float(policy.nutrition_mg_yan_l[0][1]), 76.0, abs_tol=1e-9
            )
            for policy in policies
        ),
        "masses_reconstruct_76_mgYAN_L": all(
            math.isclose(float(row["organic_product_g"]), 87.4, abs_tol=1e-9)
            and math.isclose(float(row["dap_product_g"]), 43.7, abs_tol=1e-9)
            and math.isclose(float(row["reconstructed_yan_mg_l"]), 76.0, abs_tol=1e-9)
            for row in nutrition_rows
        ),
        "net_yan_split_50_50": all(row["net_yan_split_50_50"] for row in nutrition_rows),
        "organic_limit_and_margin_met": all(
            row["nutrition_feasible"]
            and math.isclose(float(row["organic_margin_fraction"]), 0.05, abs_tol=1e-9)
            for row in nutrition_rows
        ),
    }
    search_checks = dict(balanced["convergence"]["checks"])
    sampling_checks = {
        **sampling_selected["checks"],
        "minimum_margin_24h_by_policy_member": bool(
            drying["drying_margin_pass"].all()
            and drying["action_to_drying_margin_h"].min() >= 24.0 - 1e-9
        ),
        "exactly_3_by_64_drying_rows": len(drying) == 3 * 64,
        "harmonized_loss_at_most_0_5pct": float(
            sampling_row["harmonized_score_loss_fraction"]
        )
        <= 0.005 + 1e-12,
        "operator_rounds_at_most_12": int(sampling_row["operator_rounds"]) <= 12,
        "selected_sampling_is_partially_harmonized": sampling_row["sampling_mode"]
        == "partially_harmonized",
    }
    gates = {
        "coverage": {
            "verdict": gate_verdict(coverage_checks, tuple(coverage_checks)),
            "checks": coverage_checks,
        },
        "thermal": {
            "verdict": gate_verdict(thermal_checks, tuple(thermal_checks)),
            "checks": thermal_checks,
        },
        "nutrition": {
            "verdict": gate_verdict(nutrition_checks, tuple(nutrition_checks)),
            "checks": nutrition_checks,
        },
        "search": {
            "verdict": gate_verdict(search_checks, tuple(search_checks)),
            "checks": search_checks,
        },
        "information": fim_validation["gate"],
        "sampling": {
            "verdict": gate_verdict(sampling_checks, tuple(sampling_checks)),
            "checks": sampling_checks,
        },
        "consistency": {
            "verdict": "PASS",
            "checks": {
                "canonical_state_defined_before_artifact_generation": True,
                "cross_artifact_audit_required": True,
            },
        },
        "physical_authorization": {
            "verdict": "FAIL / PENDING OWNER APPROVAL",
            "checks": {
                "profiles_for_physical_execution": False,
                "physical_execution_authorized": False,
                "executable_schedule_issued": False,
                "tank_assignments_empty": True,
            },
        },
    }
    mandatory = ("coverage", "thermal", "nutrition", "search", "information", "sampling", "consistency")
    scientific_pass = all(gates[name]["verdict"] == "PASS" for name in mandatory)
    verdict = (
        "PASS — READY FOR OWNER PHYSICAL RELEASE REVIEW"
        if scientific_pass
        else "FAIL — NOT READY FOR EXECUTION"
    )
    balanced_score = float(sampling_row["robust_score"])
    benchmark_score = float(benchmark_row["robust_score"])
    anchor_score = float(three_anchor["metrics"]["robust_score"])
    coverage_cost = benchmark_score - balanced_score
    coverage_cost_fraction = coverage_cost / max(abs(benchmark_score), 1e-12)
    retained = (balanced_score - anchor_score) / max(
        benchmark_score - anchor_score, 1e-12
    )
    profiles = []
    for policy, metrics in zip(policies, policy_metrics):
        profiles.append(
            {
                "policy": policy.name,
                "temperature_12h_blocks_to_168h_c": list(policy.temperature_c),
                "controller_blocks_to_504h": explicit_controller_blocks(
                    policy, load_json(COVERAGE_CONFIG_PATH)
                ).to_dict(orient="records"),
                "nutrition_time_h": metrics.nutrition_time_h,
                "nutrition_timestamp_local": (
                    datetime.fromisoformat(final_config["campaign_start_local"])
                    + timedelta(hours=metrics.nutrition_time_h)
                ).isoformat(),
                "yan_mg_l": 76.0,
                "first_24h_mean_setpoint_c": metrics.first_24h_mean_setpoint_c,
                "maximum_initial_jump_c": metrics.maximum_initial_jump_c,
                "maximum_internal_jump_c": metrics.maximum_internal_jump_c,
                "total_thermal_variation_c": metrics.total_thermal_variation_c,
                "temperature_changes": metrics.temperature_changes,
            }
        )
    return {
        "schema_version": 1,
        "candidate_id": selected["candidate_id"],
        "canonical_policy_hash": _policy_hash(policies),
        "final_verdict": verdict,
        "scientific_qualification_pass": scientific_pass,
        "strategy": BALANCED_STRATEGY,
        "dose_design": BALANCED_DOSE,
        "logical_treatments": [policy.name for policy in policies],
        "proposed_tank_mapping_for_review": owner[
            "proposed_tank_mapping_for_review"
        ],
        "mapping_authorized": False,
        "scores": {
            "balanced_search_robust_score": float(selected["robust_score"]),
            "balanced_final_sampling_robust_score": balanced_score,
            "feasible_information_benchmark_robust_score": benchmark_score,
            "three_anchor_N76_robust_score": anchor_score,
            "historical_infeasible_information_upper_bound_robust_score": float(
                historical_metrics["robust_score"]
            ),
            "coverage_cost_absolute": coverage_cost,
            "coverage_cost_fraction": coverage_cost_fraction,
            "retained_information_improvement_fraction_vs_three_anchor": retained,
        },
        "metrics": {
            "median": float(sampling_row["median"]),
            "q10": float(sampling_row["q10"]),
            "minimum": float(sampling_row["minimum"]),
            "tail_cvar": float(sampling_row["tail_cvar"]),
            "completion_probability": float(selected["completion_probability"]),
            "minimum_action_to_drying_margin_h": float(
                drying["action_to_drying_margin_h"].min()
            ),
            "maximum_initial_jump_c": float(
                final_physical["initial_setpoint_jump_c"].max()
            ),
            "maximum_internal_jump_c": max(
                metric.maximum_internal_jump_c for metric in policy_metrics
            ),
            "minimum_robust_physical_temperature_c": float(
                final_physical["minimum_physical_temperature_c"].min()
            ),
            "maximum_robust_physical_temperature_c": float(
                final_physical["maximum_physical_temperature_c"].max()
            ),
            "total_thermal_variation_c": sum(
                metric.total_thermal_variation_c for metric in policy_metrics
            ),
            "temperature_changes": sum(metric.temperature_changes for metric in policy_metrics),
        },
        "gates": gates,
        "dose": {
            "yan_mg_l": 76.0,
            "springferm_xtrem_g_per_230_l": 87.4,
            "dap_g_per_230_l": 43.7,
            "net_yan_fraction_springferm": 0.5,
            "net_yan_fraction_dap": 0.5,
            "organic_limit_fraction": 0.95,
        },
        "profiles": profiles,
        "sampling": {
            "mode": sampling_row["sampling_mode"],
            "robust_score": balanced_score,
            "harmonized_score_loss_fraction": float(
                sampling_row["harmonized_score_loss_fraction"]
            ),
            "operator_rounds": int(sampling_row["operator_rounds"]),
            "schedules_h": {
                key: list(value)
                for key, value in sampling_selected["schedules"].items()
            },
            "sample_event_order": "sample_before_action",
        },
        "drying": {
            "definition": "drying_time(policy, member) - latest_active_action(policy)",
            "row_count": len(drying),
            "minimum_margin_h": float(drying["action_to_drying_margin_h"].min()),
            "required_margin_h": 24.0,
        },
        "benchmark": {
            "strategy": BENCHMARK_STRATEGY,
            "candidate_id": benchmark["selected"]["candidate_id"],
            "status": benchmark["benchmark_status"],
            "robust_score": benchmark_score,
            "budget_exhausted": bool(benchmark["budget_exhausted"]),
        },
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }


def _figure_setup() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.edgecolor": PALETTE["ink"],
            "axes.labelcolor": PALETTE["ink"],
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def _finish_figure(fig: plt.Figure, path: Path, subtitle: str) -> None:
    fig.text(0.01, 0.012, subtitle, fontsize=8, color=PALETTE["ink"])
    fig.text(
        0.99,
        0.012,
        WATERMARK,
        fontsize=7,
        color=PALETTE["orange"],
        ha="right",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.98))
    fig.savefig(filesystem_path(path), dpi=180, bbox_inches="tight")
    plt.close(fig)


def _generate_figures(
    *,
    run_dir: Path,
    state: dict[str, Any],
    comparison: pd.DataFrame,
    drying: pd.DataFrame,
    physical: pd.DataFrame,
    variance_ratios: pd.DataFrame,
    eigenvalues: pd.DataFrame,
    sampling_schedule: pd.DataFrame,
) -> tuple[list[Path], dict[str, Any]]:
    _figure_setup()
    paths: list[Path] = []
    fig, axis = plt.subplots(figsize=(9.5, 5.5))
    plot = comparison[comparison["comparator_type"].ne("historical_upper_bound")]
    axis.barh(
        plot["strategy"],
        plot["robust_score"],
        color=[PALETTE["blue"], PALETTE["gold"]],
        edgecolor=PALETTE["ink"],
    )
    for index, value in enumerate(plot["robust_score"]):
        axis.text(float(value), index, f" {float(value):.3f}", va="center")
    axis.set_xlabel("Robust Bayesian D-information score")
    axis.set_title("Operational coverage versus feasible information benchmark")
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    path = run_dir / "coverage_vs_feasible_information.png"
    _finish_figure(fig, path, "64-member final sampling comparison; exact values shown")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(9.5, 5.5))
    for profile, color in zip(state["profiles"], (PALETTE["grey"], PALETTE["blue"], PALETTE["orange"])):
        values = profile["temperature_12h_blocks_to_168h_c"]
        axis.step(
            np.arange(len(values) + 1) * 12.0,
            [*values, values[-1]],
            where="post",
            label=profile["policy"],
            color=color,
            linewidth=2.2,
        )
    axis.set_xlabel("Hours from inoculation")
    axis.set_ylabel("Setpoint (°C)")
    axis.set_ylim(15, 27)
    axis.set_title("Final Wave 1 temperature profiles")
    axis.legend(frameon=False, loc="best")
    axis.grid(color="#E5E7EB", linewidth=0.8)
    path = run_dir / "final_temperature_profiles.png"
    _finish_figure(fig, path, "12 h blocks to 168 h; controller holds the last setpoint to 504 h")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(8.5, 6.0))
    points = [(row["nutrition_time_h"], row["first_24h_mean_setpoint_c"], row["policy"]) for row in state["profiles"]]
    for (x, y, label), color in zip(points, (PALETTE["grey"], PALETTE["blue"], PALETTE["orange"])):
        axis.scatter(x, y, s=90, color=color, edgecolor=PALETTE["ink"], zorder=3)
        axis.annotate(label, (x, y), xytext=(7, 7), textcoords="offset points", fontsize=9)
    axis.axvspan(0, 6, color=PALETTE["blue"], alpha=0.08)
    axis.axvspan(42, 50, color=PALETTE["orange"], alpha=0.08)
    axis.axhline(18, color=PALETTE["grey"], linestyle="--")
    axis.axhline(22, color=PALETTE["grey"], linestyle="--")
    axis.set_xlabel("Nutrition time (h)")
    axis.set_ylabel("First-24 h mean setpoint (°C)")
    axis.set_title("Wave 1 operational coverage map")
    axis.grid(color="#E5E7EB", linewidth=0.8)
    path = run_dir / "operational_coverage_map.png"
    _finish_figure(fig, path, "Logical treatment geometry; tank mapping remains unauthorised")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(9.5, 5.5))
    final_physical = physical[physical["strategy"].eq(BALANCED_STRATEGY)]
    summary = final_physical.groupby("candidate").agg(
        minimum=("minimum_physical_temperature_c", "min"),
        maximum=("maximum_physical_temperature_c", "max"),
    )
    y = np.arange(len(summary))
    axis.hlines(y, summary["minimum"], summary["maximum"], color=PALETTE["blue"], linewidth=5)
    axis.scatter(summary["minimum"], y, color=PALETTE["blue"], edgecolor=PALETTE["ink"], label="Minimum")
    axis.scatter(summary["maximum"], y, color=PALETTE["gold"], edgecolor=PALETTE["ink"], label="Maximum")
    axis.axvline(15, color=PALETTE["ink"], linestyle="--")
    axis.axvline(27, color=PALETTE["ink"], linestyle="--")
    axis.set_yticks(y, summary.index)
    axis.set_xlabel("Robust physical temperature (°C)")
    axis.set_title("Physical temperature envelope across actuator scenarios")
    axis.legend(frameon=False)
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    path = run_dir / "physical_temperature_envelope.png"
    _finish_figure(fig, path, "243 approved actuator scenarios per policy")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(10.5, 6.0))
    balanced_ratios = variance_ratios[
        variance_ratios["strategy"].eq(BALANCED_STRATEGY)
    ].sort_values("variance_ratio_vs_feasible_information_benchmark")
    axis.barh(
        balanced_ratios["parameter"],
        balanced_ratios["variance_ratio_vs_feasible_information_benchmark"],
        color=PALETTE["olive"],
        edgecolor=PALETTE["ink"],
    )
    axis.axvline(1.0, color=PALETTE["ink"], linestyle="--")
    axis.set_xlabel("Posterior variance ratio vs feasible information benchmark")
    axis.set_title("Parameter-level information tradeoff")
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    path = run_dir / "parameter_information_comparison.png"
    _finish_figure(fig, path, "Values above 1 indicate higher posterior variance in the balanced design")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(10.0, 5.5))
    for index, (policy, group) in enumerate(sampling_schedule.groupby("policy", sort=False)):
        axis.scatter(group["time_h"], np.full(len(group), index), s=45, label=policy)
    axis.set_yticks(range(sampling_schedule["policy"].nunique()), sampling_schedule["policy"].drop_duplicates())
    axis.set_xlabel("Hours from inoculation")
    axis.set_title("Selected partially harmonized sampling calendar")
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    path = run_dir / "wave1_sampling_calendar.png"
    _finish_figure(fig, path, "Ten samples per tank; basal at t=0; sample_before_action")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(9.5, 5.5))
    ordered = drying.sort_values("action_to_drying_margin_h")
    axis.plot(
        np.arange(len(ordered)),
        ordered["action_to_drying_margin_h"],
        color=PALETTE["blue"],
        linewidth=1.8,
    )
    axis.axhline(24.0, color=PALETTE["orange"], linestyle="--", label="Required 24 h")
    axis.set_xlabel("Policy-member record, ordered")
    axis.set_ylabel("Action-to-drying margin (h)")
    axis.set_title("Drying margin by policy and ensemble member")
    axis.legend(frameon=False)
    axis.grid(color="#E5E7EB", linewidth=0.8)
    path = run_dir / "drying_margin_by_policy_member.png"
    _finish_figure(fig, path, "Each margin pairs drying and latest active action from the same policy")
    paths.append(path)

    fig, axis = plt.subplots(figsize=(9.5, 5.5))
    first_eigen = eigenvalues[eigenvalues["eigenvalue_rank_ascending"].eq(1)]
    axis.barh(
        first_eigen["strategy"],
        first_eigen["fim_eigenvalue_median"],
        color=[PALETTE["grey"], PALETTE["gold"], PALETTE["blue"]][: len(first_eigen)],
        edgecolor=PALETTE["ink"],
    )
    axis.set_xlabel("Median smallest posterior FIM eigenvalue")
    axis.set_title("Least-identified FIM direction by strategy")
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    path = run_dir / "fim_eigenvalue_comparison.png"
    _finish_figure(fig, path, "64-member posterior FIM comparison")
    paths.append(path)

    qa_rows = []
    for path in paths:
        image = plt.imread(filesystem_path(path))
        checks = {
            "file_nonempty": filesystem_path(path).stat().st_size > 10_000,
            "width_at_least_1200px": int(image.shape[1]) >= 1200,
            "height_at_least_700px": int(image.shape[0]) >= 700,
            "finite_pixels": bool(np.isfinite(image).all()),
            "not_blank": float(np.std(image)) > 0.01,
        }
        qa_rows.append(
            {
                "figure": path.name,
                "width_px": int(image.shape[1]),
                "height_px": int(image.shape[0]),
                "checks": checks,
                "verdict": gate_verdict(checks, tuple(checks)),
            }
        )
    qa = {
        "verdict": "PASS" if all(row["verdict"] == "PASS" for row in qa_rows) else "FAIL",
        "figures": qa_rows,
        "candidate_id": state["candidate_id"],
        "canonical_policy_hash": state["canonical_policy_hash"],
        "final_verdict": state["final_verdict"],
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }
    return paths, qa


def _ics_datetime(value: datetime) -> str:
    return value.strftime("%Y%m%dT%H%M%S")


def _ics_event(
    *, uid: str, start: datetime, summary: str, description: str
) -> list[str]:
    return [
        "BEGIN:VEVENT",
        f"UID:{uid}",
        f"DTSTAMP:{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}",
        f"DTSTART;TZID=America/Santiago:{_ics_datetime(start)}",
        f"SUMMARY:{summary}",
        f"DESCRIPTION:{description}",
        "END:VEVENT",
    ]


def _build_protocol_package(
    *,
    run_dir: Path,
    state: dict[str, Any],
    sampling_selected: dict[str, Any],
    capture_selected: pd.DataFrame,
    drying: pd.DataFrame,
    figure_paths: list[Path],
) -> list[Path]:
    package = run_dir / "final_protocol_review_package"
    subdirs = (
        "controller",
        "calendar",
        "nutrition",
        "sampling",
        "capture",
        "tanks",
        "checklists",
        "figures",
        "authorization",
    )
    for name in subdirs:
        filesystem_path(package / name).mkdir(parents=True, exist_ok=True)
    files: list[Path] = []
    profile_by_name = {row["policy"]: row for row in state["profiles"]}
    mapping = state["proposed_tank_mapping_for_review"]
    reverse = {policy: tank for tank, policy in mapping.items()}
    controller_rows = []
    controller_gui = {
        "candidate_id": state["candidate_id"],
        "canonical_policy_hash": state["canonical_policy_hash"],
        "review_only": True,
        "controllers": [],
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }
    for policy, profile in profile_by_name.items():
        tank = reverse[policy]
        rows = []
        for block in profile["controller_blocks_to_504h"]:
            row = {
                "tank_review_mapping": tank,
                "policy": policy,
                "block_index": block["block_index"],
                "start_h": block["start_h"],
                "end_h": block["end_h"],
                "setpoint_c": block["setpoint_c"],
                "phase": block["phase"],
                "hold_last_setpoint": block["hold_last_setpoint"],
            }
            rows.append(row)
            controller_rows.append(row)
        controller_gui["controllers"].append(
            {
                "tank_review_mapping": tank,
                "policy": policy,
                "blocks": rows,
                "checksum_sha256": sha256_payload(rows),
            }
        )
    controller_csv = package / "controller" / f"controller_profiles{REVIEW_SUFFIX}.csv"
    pd.DataFrame(controller_rows).assign(
        watermark=WATERMARK,
        candidate_only_not_for_physical_execution=True,
    ).to_csv(filesystem_path(controller_csv), index=False, lineterminator="\n")
    files.append(controller_csv)
    controller_json = package / "controller" / f"controller_profiles_gui{REVIEW_SUFFIX}.json"
    write_json(controller_json, controller_gui)
    files.append(controller_json)
    checksum_path = package / "controller" / f"controller_checksums{REVIEW_SUFFIX}.json"
    write_json(
        checksum_path,
        {
            "candidate_id": state["candidate_id"],
            "canonical_policy_hash": state["canonical_policy_hash"],
            "files": {
                controller_csv.name: sha256_file(controller_csv),
                controller_json.name: sha256_file(controller_json),
            },
            "watermark": WATERMARK,
            **AUTHORIZATION_FLAGS,
        },
    )
    files.append(checksum_path)
    schedule = sampling_selected["schedule_frame"].copy()
    schedule["tank_review_mapping"] = schedule["policy"].map(reverse)
    schedule["candidate_id"] = state["candidate_id"]
    schedule["canonical_policy_hash"] = state["canonical_policy_hash"]
    schedule["watermark"] = WATERMARK
    schedule["candidate_only_not_for_physical_execution"] = True
    sampling_path = package / "sampling" / f"sampling_calendar{REVIEW_SUFFIX}.csv"
    schedule.to_csv(filesystem_path(sampling_path), index=False, lineterminator="\n")
    files.append(sampling_path)
    captures = capture_selected.copy()
    captures["tank_review_mapping"] = captures["policy"].map(reverse)
    captures["candidate_id"] = state["candidate_id"]
    captures["canonical_policy_hash"] = state["canonical_policy_hash"]
    captures["watermark"] = WATERMARK
    captures["candidate_only_not_for_physical_execution"] = True
    capture_path = package / "capture" / f"capture_intervals{REVIEW_SUFFIX}.csv"
    captures.to_csv(filesystem_path(capture_path), index=False, lineterminator="\n")
    files.append(capture_path)
    nutrition_rows = []
    campaign_start = datetime.fromisoformat(state["profiles"][0]["nutrition_timestamp_local"]) - timedelta(
        hours=state["profiles"][0]["nutrition_time_h"]
    )
    for profile in state["profiles"]:
        nutrition_rows.append(
            {
                "tank_review_mapping": reverse[profile["policy"]],
                "policy": profile["policy"],
                "nutrition_time_h": profile["nutrition_time_h"],
                "nutrition_timestamp_local": profile["nutrition_timestamp_local"],
                "sequence_if_sample_coincident": "sample_then_nutrition",
                "springferm_xtrem_g": 87.4,
                "dap_g": 43.7,
                "yan_mg_l": 76.0,
                "cumulative_springferm_xtrem_g": 87.4,
                "cumulative_dap_g": 43.7,
                "organic_limit_margin_fraction": 0.05,
                "watermark": WATERMARK,
                "candidate_only_not_for_physical_execution": True,
            }
        )
    nutrition_path = package / "nutrition" / f"nutrition_plan{REVIEW_SUFFIX}.csv"
    pd.DataFrame(nutrition_rows).to_csv(
        filesystem_path(nutrition_path), index=False, lineterminator="\n"
    )
    files.append(nutrition_path)
    tank_path = package / "tanks" / f"proposed_logical_tank_mapping{REVIEW_SUFFIX}.json"
    write_json(
        tank_path,
        {
            "candidate_id": state["candidate_id"],
            "canonical_policy_hash": state["canonical_policy_hash"],
            "proposed_tank_mapping_for_review": mapping,
            "mapping_authorized": False,
            "watermark": WATERMARK,
            **AUTHORIZATION_FLAGS,
        },
    )
    files.append(tank_path)
    for tank, policy in mapping.items():
        lines = [
            "BEGIN:VCALENDAR",
            "VERSION:2.0",
            "PRODID:-//pyomo-doe//Wave1 Review Only//ES",
            "CALSCALE:GREGORIAN",
            "METHOD:PUBLISH",
        ]
        lines.extend(
            _ics_event(
                uid=f"{state['candidate_id']}-{tank}-basal",
                start=campaign_start,
                summary=f"{tank} basal e inoculación (REVISIÓN)",
                description=f"{WATERMARK}; sample_before_action",
            )
        )
        tank_samples = schedule[schedule["policy"].eq(policy)]
        for _, sample in tank_samples.iterrows():
            when = datetime.fromisoformat(str(sample["local_timestamp"]))
            lines.extend(
                _ics_event(
                    uid=f"{state['candidate_id']}-{tank}-sample-{int(sample['sample_number'])}",
                    start=when,
                    summary=f"{tank} muestra {int(sample['sample_number'])} (REVISIÓN)",
                    description=f"{WATERMARK}; sample_before_action",
                )
            )
        profile = profile_by_name[policy]
        nutrition_when = datetime.fromisoformat(profile["nutrition_timestamp_local"])
        lines.extend(
            _ics_event(
                uid=f"{state['candidate_id']}-{tank}-nutrition",
                start=nutrition_when,
                summary=f"{tank} nutrición N76 (REVISIÓN)",
                description=f"{WATERMARK}; 87.4 g SpringFerm Xtrem + 43.7 g DAP; sample_then_nutrition if coincident",
            )
        )
        prior_setpoint = None
        for block in profile["controller_blocks_to_504h"]:
            setpoint = float(block["setpoint_c"])
            if prior_setpoint is None or not math.isclose(setpoint, prior_setpoint, abs_tol=1e-12):
                when = campaign_start + timedelta(hours=float(block["start_h"]))
                lines.extend(
                    _ics_event(
                        uid=f"{state['candidate_id']}-{tank}-temp-{int(block['block_index'])}",
                        start=when,
                        summary=f"{tank} cambio automático informativo a {setpoint:.2f} C",
                        description=f"{WATERMARK}; evento informativo no ejecutable",
                    )
                )
            prior_setpoint = setpoint
        lines.extend(
            _ics_event(
                uid=f"{state['candidate_id']}-{tank}-capture-start",
                start=campaign_start,
                summary=f"{tank} inicio de captura (REVISIÓN)",
                description=WATERMARK,
            )
        )
        lines.extend(
            _ics_event(
                uid=f"{state['candidate_id']}-{tank}-close",
                start=campaign_start + timedelta(hours=504.0),
                summary=f"{tank} cierre de horizonte de controlador (REVISIÓN)",
                description=WATERMARK,
            )
        )
        lines.append("END:VCALENDAR")
        ics_path = package / "calendar" / f"{tank}_calendar{REVIEW_SUFFIX}.ics"
        filesystem_path(ics_path).write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")
        files.append(ics_path)
        checklist = package / "checklists" / f"{tank}_checklist{REVIEW_SUFFIX}.md"
        checklist_text = "\n".join(
            [
                f"# {tank} — {policy} (REVIEW ONLY)",
                "",
                WATERMARK,
                "",
                "- [ ] Aprobación física explícita del owner (aún pendiente)",
                "- [ ] Carga del perfil y verificación de checksum",
                "- [ ] Inoculación y muestra basal",
                "- [ ] Inicio de captura",
                "- [ ] Nutrición; tomar muestra antes si coincide",
                "- [ ] Muestreos en ventanas legales",
                "- [ ] Vigilancia de temperatura y registro de desvíos",
                "- [ ] Cierre de captura y proceso",
                "- [ ] Firma operador: ____________________",
                "- [ ] Firma revisor: ____________________",
                "",
                "Este checklist no es vigente hasta la aprobación física del owner.",
            ]
        )
        filesystem_path(checklist).write_text(checklist_text + "\n", encoding="utf-8")
        files.append(checklist)
    drying_package = package / "calendar" / f"drying_margin_by_policy_member{REVIEW_SUFFIX}.csv"
    drying.assign(
        candidate_id=state["candidate_id"],
        canonical_policy_hash=state["canonical_policy_hash"],
        watermark=WATERMARK,
        candidate_only_not_for_physical_execution=True,
    ).to_csv(filesystem_path(drying_package), index=False, lineterminator="\n")
    files.append(drying_package)
    for figure in figure_paths:
        destination = package / "figures" / f"{figure.stem}{REVIEW_SUFFIX}{figure.suffix}"
        shutil.copy2(filesystem_path(figure), filesystem_path(destination))
        files.append(destination)
    authorization_path = package / "authorization" / f"authorization_state{REVIEW_SUFFIX}.json"
    write_json(
        authorization_path,
        {
            "candidate_id": state["candidate_id"],
            "canonical_policy_hash": state["canonical_policy_hash"],
            "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
            "mapping_authorized": False,
            "watermark": WATERMARK,
            **AUTHORIZATION_FLAGS,
        },
    )
    files.append(authorization_path)
    manifest_path = package / f"package_manifest{REVIEW_SUFFIX}.json"
    manifest = {
        "candidate_id": state["candidate_id"],
        "canonical_policy_hash": state["canonical_policy_hash"],
        "final_verdict": state["final_verdict"],
        "review_only": True,
        "files": [
            {
                "path": str(path.relative_to(package)).replace("\\", "/"),
                "sha256": sha256_file(path),
                "bytes": filesystem_path(path).stat().st_size,
            }
            for path in sorted(files, key=lambda item: str(item))
        ],
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }
    write_json(manifest_path, manifest)
    files.append(manifest_path)
    return files


def _cross_artifact_consistency_audit(
    *, run_dir: Path, state: dict[str, Any], package_manifest: Path | None
) -> dict[str, Any]:
    expected = {
        "candidate_id": state["candidate_id"],
        "canonical_policy_hash": state["canonical_policy_hash"],
        "final_verdict": state["final_verdict"],
        "drying_margin_h": state["metrics"]["minimum_action_to_drying_margin_h"],
        "sampling_gate": state["gates"]["sampling"]["verdict"],
        "coverage_gate": state["gates"]["coverage"]["verdict"],
        "thermal_gate": state["gates"]["thermal"]["verdict"],
        "information_gate": state["gates"]["information"]["verdict"],
        "nutrition_gate": state["gates"]["nutrition"]["verdict"],
    }
    checks: dict[str, bool] = {}
    state_disk = _json(run_dir / "qualification_state.json")
    checks["qualification_state_candidate_id"] = state_disk["candidate_id"] == expected["candidate_id"]
    checks["qualification_state_policy_hash"] = state_disk["canonical_policy_hash"] == expected["canonical_policy_hash"]
    checks["qualification_state_verdict"] = state_disk["final_verdict"] == expected["final_verdict"]
    gate_summary = _json(run_dir / "final_gate_summary.json")
    checks["gate_summary_candidate_id"] = gate_summary["candidate_id"] == expected["candidate_id"]
    checks["gate_summary_policy_hash"] = gate_summary["canonical_policy_hash"] == expected["canonical_policy_hash"]
    checks["gate_summary_verdict"] = gate_summary["final_verdict"] == expected["final_verdict"]
    checks["all_gate_verdicts_match"] = all(
        gate_summary["gates"][name]["verdict"] == state["gates"][name]["verdict"]
        for name in ("coverage", "thermal", "information", "nutrition", "sampling")
    )
    provenance = _json(run_dir / "selected_candidate_provenance.json")
    checks["provenance_candidate_id"] = provenance["candidate_id"] == expected["candidate_id"]
    checks["provenance_policy_hash"] = provenance["canonical_policy_hash"] == expected["canonical_policy_hash"]
    drying = pd.read_csv(filesystem_path(run_dir / "drying_margin_by_policy_member.csv"))
    checks["drying_candidate_id"] = drying["qualification_candidate_id"].eq(expected["candidate_id"]).all()
    checks["drying_policy_hash"] = drying["canonical_policy_hash"].eq(expected["canonical_policy_hash"]).all()
    checks["drying_margin_matches"] = math.isclose(
        float(drying["action_to_drying_margin_h"].min()),
        float(expected["drying_margin_h"]),
        abs_tol=1e-9,
    )
    sampling = pd.read_csv(filesystem_path(run_dir / "sampling_strategy_comparison.csv"))
    selected_sampling = sampling[
        sampling["strategy"].eq(BALANCED_STRATEGY)
        & sampling["sampling_mode"].eq(state["sampling"]["mode"])
    ]
    checks["sampling_candidate_metadata"] = bool(
        len(selected_sampling) == 1
        and selected_sampling["qualification_candidate_id"].eq(expected["candidate_id"]).all()
        and selected_sampling["canonical_policy_hash"].eq(expected["canonical_policy_hash"]).all()
    )
    visual = _json(run_dir / "visual_qa.json")
    checks["visual_qa_candidate_id"] = visual["candidate_id"] == expected["candidate_id"]
    checks["visual_qa_policy_hash"] = visual["canonical_policy_hash"] == expected["canonical_policy_hash"]
    checks["visual_qa_pass"] = visual["verdict"] == "PASS"
    if package_manifest is not None:
        package = _json(package_manifest)
        checks["package_candidate_id"] = package["candidate_id"] == expected["candidate_id"]
        checks["package_policy_hash"] = package["canonical_policy_hash"] == expected["canonical_policy_hash"]
        checks["package_verdict"] = package["final_verdict"] == expected["final_verdict"]
    else:
        checks["package_omitted_only_on_failed_scientific_gate"] = not bool(
            state["scientific_qualification_pass"]
        )
    verdict = gate_verdict(checks, tuple(checks))
    return {
        "verdict": verdict,
        "checks": checks,
        "expected": expected,
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }


def _flat_record(row: dict[str, Any]) -> dict[str, Any]:
    excluded = {
        "vector",
        "policies",
        "evaluations",
        "archetype_checks",
        "contrast_checks",
        "drying_rows",
    }
    output = {key: value for key, value in row.items() if key not in excluded}
    if "vector" in row:
        output["vector_sha256"] = sha256_payload(
            np.asarray(row["vector"], dtype=float).tolist()
        )
    return output


def _preflight() -> dict[str, Any]:
    final_config = load_json(FINAL_CONFIG_PATH)
    coverage_config = load_json(COVERAGE_CONFIG_PATH)
    dependencies = _dependency_manifest()
    git = capture_git_state()
    import pyomo.environ as pyo

    solver = pyo.SolverFactory("ipopt")
    warm = _warm_start_positions(final_config)
    expected_start = "9ecc162668deb8f21b86d0188bddd754dfc99b59"
    start_is_ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", expected_start, str(git.get("commit"))],
        cwd=REPOSITORY_DIR,
        check=False,
        capture_output=True,
    ).returncode == 0
    checks = {
        "branch_is_wave1_operational_coverage": git.get("branch")
        == "wave1_operational_coverage",
        "expected_remote_start_is_ancestor": start_is_ancestor,
        "working_tree_clean": git.get("dirty") is False,
        "watermark_exact": final_config["watermark"] == WATERMARK
        and coverage_config["watermark"] == WATERMARK,
        "authorization_fail_closed": final_config["authorization"]
        == AUTHORIZATION_FLAGS,
        "owner_strategy_fixed": final_config["owner_selection"]["strategy"]
        == BALANCED_STRATEGY,
        "n76_fixed": float(final_config["owner_selection"]["yan_mg_l"]) == 76.0,
        "balanced_minimum_budget_declared": len(
            final_config["balanced_search"]["seeds"]
        )
        == 5
        and int(final_config["balanced_search"]["particles"]) >= 24
        and int(final_config["balanced_search"]["four_member_iterations"]) >= 10
        and int(final_config["balanced_search"]["eight_member_iterations"]) >= 8,
        "benchmark_minimum_budget_declared": len(
            final_config["feasible_information_benchmark"]["seeds"]
        )
        == 3
        and int(final_config["feasible_information_benchmark"]["particles"]) >= 16
        and int(
            final_config["feasible_information_benchmark"]["four_member_iterations"]
        )
        >= 6
        and int(
            final_config["feasible_information_benchmark"]["eight_member_iterations"]
        )
        >= 4,
        "continuation_threshold_is_0_1pct": math.isclose(
            float(
                final_config["balanced_search"][
                    "continuation_improvement_tolerance_fraction"
                ]
            ),
            0.001,
            abs_tol=0.0,
        ),
        "drying_margin_is_24h": math.isclose(
            float(final_config["required_minimum_action_to_drying_margin_h"]),
            24.0,
            abs_tol=0.0,
        ),
        "ipopt_available": bool(solver.available(exception_flag=False)),
        "warm_starts_available": len(warm) >= 3 and warm.shape[1] == 24,
        "transitive_dependencies_hashed": dependencies[
            "transitive_dependency_count"
        ]
        == len(DEPENDENCY_PATHS),
    }
    return {
        "verdict": gate_verdict(checks, tuple(checks)),
        "checks": checks,
        "git": git,
        "scientific_dependency_sha256": dependencies[
            "scientific_dependency_sha256"
        ],
        "warm_start_count": len(warm),
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the single final Wave 1 operational-coverage qualification"
    )
    parser.add_argument(
        "--result-root",
        default=str(RESULT_ROOT),
        help="Adaptive-design result root; the final stage and one immutable run are created",
    )
    parser.add_argument(
        "--resume-run",
        help="Resume the exact same incomplete final run from validated checkpoints",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Validate code, configuration, sources, warm starts, and IPOPT without creating results",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.preflight_only:
        result = _preflight()
        print(json.dumps(result, ensure_ascii=False, sort_keys=True), flush=True)
        if result["verdict"] != "PASS":
            raise SystemExit(2)
        return
    started = time.perf_counter()
    stage_times: dict[str, float] = {}
    final_config = load_json(FINAL_CONFIG_PATH)
    coverage_config = load_json(COVERAGE_CONFIG_PATH)
    if final_config["watermark"] != WATERMARK:
        raise ValueError("Final watermark differs from the exact required text")
    if final_config["authorization"] != AUTHORIZATION_FLAGS:
        raise ValueError("Final configuration must remain fail closed")
    dependencies = _dependency_manifest()
    dependency_hash = str(dependencies["scientific_dependency_sha256"])
    config_hash = sha256_file(FINAL_CONFIG_PATH)
    git_snapshot = capture_git_state()
    if args.resume_run:
        run_dir = _run_path(args.resume_run)
        if run_dir.parent.name != final_config["result_stage"]:
            raise ValueError("Resume path is not a final operational-coverage run")
        if filesystem_path(run_dir / "run_manifest.json").is_file():
            manifest = _json(run_dir / "run_manifest.json")
            if manifest.get("status") == "completed":
                print(
                    json.dumps(
                        {
                            "run_dir": relative_or_absolute(run_dir),
                            "status": "already_completed_no_restart",
                            "final_verdict": manifest.get("final_verdict"),
                            **AUTHORIZATION_FLAGS,
                        }
                    ),
                    flush=True,
                )
                return
    else:
        if git_snapshot.get("dirty") is not False:
            raise RuntimeError("The official run requires a clean Git worktree")
        run_dir = create_immutable_run_directory(
            _run_path(args.result_root), final_config["result_stage"], final_config
        )
    outputs: list[Path] = []
    final_config_output = run_dir / "coverage_final_design_config.json"
    if not filesystem_path(final_config_output).exists():
        _write_json(final_config_output, final_config, outputs)
    else:
        if sha256_file(final_config_output) != sha256_file(FINAL_CONFIG_PATH):
            raise RuntimeError("Resume final configuration hash mismatch")
        outputs.append(final_config_output)
    dependency_output = run_dir / "scientific_dependency_manifest.json"
    if filesystem_path(dependency_output).exists():
        prior_dependency = _json(dependency_output)
        if prior_dependency.get("scientific_dependency_sha256") != dependency_hash:
            raise RuntimeError("Resume scientific dependency hash mismatch")
    model_config = _model_config(final_config)
    source_started = time.perf_counter()
    runs, source_verification = _source_contract(coverage_config)
    dependencies["verified_source_outputs"] = source_verification
    dependencies["all_source_output_hashes_verified"] = all(
        value.get("declared_output_sha256")
        in {
            value.get("working_tree_output_sha256"),
            value.get("canonical_lf_output_sha256"),
        }
        for value in source_verification.values()
    )
    if not dependencies["all_source_output_hashes_verified"]:
        raise RuntimeError("At least one immutable source output hash failed")
    if not filesystem_path(dependency_output).exists():
        _write_json(dependency_output, dependencies, outputs)
    else:
        outputs.append(dependency_output)
    envelope = derive_robust_setpoint_envelope(model_config, coverage_config)
    envelope_path = run_dir / "robust_setpoint_envelope.json"
    if not filesystem_path(envelope_path).exists():
        _write_json(envelope_path, envelope, outputs)
    else:
        outputs.append(envelope_path)
    ensemble = pd.read_csv(
        filesystem_path(runs["joint_ensemble"] / "joint_parameter_ensemble.csv")
    )
    if len(ensemble) != 64:
        raise ValueError("Final qualification requires exactly 64 ensemble members")
    aroma_config = load_json(AROMA_CONFIG_PATH)
    partitions, partition_provenance = load_partition_surrogates(
        aroma_config, REPOSITORY_DIR
    )
    prior = prior_precision(ensemble, model_config)
    stage_times["preflight_and_sources_seconds"] = time.perf_counter() - source_started
    design_cache: dict[tuple, tuple[np.ndarray, float, float]] = {}
    physical_cache: dict[tuple[Any, ...], tuple[pd.DataFrame, dict[str, Any]]] = {}

    balanced_started = time.perf_counter()
    balanced = _balanced_search(
        run_dir=run_dir,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        final_config=final_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        physical_cache=physical_cache,
        dependency_hash=dependency_hash,
        config_hash=config_hash,
    )
    stage_times["balanced_search_and_ipopt_seconds"] = time.perf_counter() - balanced_started
    outputs.extend(balanced["checkpoint_paths"])

    benchmark_started = time.perf_counter()
    benchmark = _benchmark_search(
        run_dir=run_dir,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        final_config=final_config,
        envelope=envelope,
        partitions=partitions,
        design_cache=design_cache,
        physical_cache=physical_cache,
        dependency_hash=dependency_hash,
        config_hash=config_hash,
    )
    stage_times["feasible_information_benchmark_seconds"] = time.perf_counter() - benchmark_started
    outputs.extend(benchmark["checkpoint_paths"])

    post_started = time.perf_counter()
    sampling = _sampling_and_information(
        balanced=balanced,
        benchmark=benchmark,
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        coverage_config=coverage_config,
        partitions=partitions,
    )
    final_policies = sampling["balanced_policies"]
    drying = drying_margin_by_policy_member(
        final_policies,
        ensemble,
        model_config,
        partitions,
        design_cache=design_cache,
    )
    drying["strategy"] = BALANCED_STRATEGY
    drying["dose_design"] = BALANCED_DOSE
    drying["candidate_id"] = balanced["selected"]["candidate_id"]
    physical_frames = []
    for strategy, candidate_id, policies in (
        (BALANCED_STRATEGY, balanced["selected"]["candidate_id"], final_policies),
        (BENCHMARK_STRATEGY, benchmark["selected"]["candidate_id"], tuple(benchmark["selected"]["policies"])),
    ):
        for policy in policies:
            frame, _summary = _physical_summary_cached(
                policy, model_config, coverage_config, physical_cache
            )
            physical_frames.append(
                frame.assign(
                    strategy=strategy,
                    dose_design=BALANCED_DOSE,
                    candidate_id=candidate_id,
                )
            )
    physical = pd.concat(physical_frames, ignore_index=True)
    fim_validation = _post_search_fim_validation(
        policies=final_policies,
        schedules=sampling["balanced_selected"]["schedules"],
        ensemble=ensemble,
        evaluations=balanced["selected"]["evaluations"],
        drying=drying,
        model_config=model_config,
        final_config=final_config,
        partitions=partitions,
        design_cache=design_cache,
    )
    reference_policies = _load_policies(
        runs["information_reference"] / "final_candidate_policy_actions.csv",
        model_config,
    )
    historical_context = _reference_information_context(
        reference_policies=reference_policies,
        source_schedule=_reference_schedule(
            runs["sampling_reference"] / "optimized_sampling_schedule.csv"
        ),
        source_comparison=pd.read_csv(
            filesystem_path(
                runs["sampling_reference"] / "sampling_candidate_comparison.csv"
            )
        ),
        ensemble=ensemble,
        prior=prior,
        model_config=model_config,
        partitions=partitions,
        prepared_cache={},
        fim_cache={},
    )
    state = _final_state(
        balanced=balanced,
        benchmark=benchmark,
        sampling=sampling,
        drying=drying,
        physical=physical,
        fim_validation=fim_validation,
        final_config=final_config,
        model_config=model_config,
        historical_metrics=historical_context["current"]["metrics"],
    )
    qualification_path = run_dir / "qualification_state.json"
    _write_json(qualification_path, state, outputs)

    champions_path = run_dir / "balanced_per_seed_champions_64_member.csv"
    _write_csv(
        champions_path,
        _tag_frame(pd.DataFrame([_flat_record(row) for row in balanced["champions"]]), state),
        outputs,
    )
    convergence_path = run_dir / "balanced_practical_convergence_gate.json"
    _write_json(
        convergence_path,
        _artifact_payload(state, balanced["convergence"]),
        outputs,
    )
    trace_path = run_dir / "balanced_local_refinement_trace.csv"
    tagged_trace = _tag_frame(balanced["local_trace"], state)
    _write_csv(trace_path, tagged_trace, outputs)
    legacy_trace_path = run_dir / "local_refinement_trace.csv"
    _write_csv(legacy_trace_path, tagged_trace, outputs)
    rejected_path = run_dir / "local_refinement_rejected_candidates.json"
    _write_json(
        rejected_path,
        _artifact_payload(state, {"rejected_candidates": balanced["rejected"]}),
        outputs,
    )
    eligible_frame = pd.DataFrame([_flat_record(row) for row in balanced["records"]])
    eligible_path = run_dir / "balanced_eligible_candidate_comparison.csv"
    _write_csv(eligible_path, _tag_frame(eligible_frame, state), outputs)
    complete_eligible_path = run_dir / "eligible_candidate_comparison.csv"
    _write_csv(complete_eligible_path, _tag_frame(eligible_frame, state), outputs)
    provenance_path = run_dir / "selected_candidate_provenance.json"
    _write_json(
        provenance_path,
        _artifact_payload(
            state,
            {
                "selected_source_type": balanced["selected"]["source_type"],
                "parent_candidate_id": balanced["selected"].get("parent_candidate_id"),
                "independent_seed": int(balanced["selected"]["independent_seed"]),
                "local_refinement_qualified": bool(
                    balanced["selected"].get("local_refinement_qualified")
                ),
                "local_refinement_state": balanced["selected"].get(
                    "local_refinement_state"
                ),
                "selection_rule": balanced["selection"],
                "compared_candidate_ids": sorted(
                    {str(row["candidate_id"]) for row in balanced["records"]}
                ),
            },
        ),
        outputs,
    )
    benchmark_path = run_dir / "feasible_information_benchmark.csv"
    _write_csv(
        benchmark_path,
        _tag_frame(
            pd.DataFrame([_flat_record(row) for row in benchmark["records"]]), state
        ),
        outputs,
    )
    balanced_row = sampling["balanced_selected"]["row"]
    benchmark_row = sampling["benchmark_selected"]["row"]
    comparison = pd.DataFrame(
        [
            {
                "strategy": BALANCED_STRATEGY,
                "comparator_type": "balanced_final",
                **{
                    key: balanced_row[key]
                    for key in ("robust_score", "median", "q10", "minimum", "tail_cvar")
                },
                "feasible": True,
                "benchmark_status": "not_applicable",
            },
            {
                "strategy": BENCHMARK_STRATEGY,
                "comparator_type": "feasible_information_benchmark",
                **{
                    key: benchmark_row[key]
                    for key in ("robust_score", "median", "q10", "minimum", "tail_cvar")
                },
                "feasible": bool(benchmark["selected"]["feasible"]),
                "benchmark_status": benchmark["benchmark_status"],
            },
            {
                "strategy": final_config["historical_comparator_label"],
                "comparator_type": "historical_upper_bound",
                **{
                    key: historical_context["current"]["metrics"][key]
                    for key in ("robust_score", "median", "q10", "minimum", "tail_cvar")
                },
                "feasible": False,
                "benchmark_status": "historical_infeasible_comparator_only",
            },
        ]
    )
    comparison["absolute_score_cost_vs_feasible_benchmark"] = (
        float(benchmark_row["robust_score"]) - comparison["robust_score"]
    )
    comparison["score_cost_fraction_vs_feasible_benchmark"] = comparison[
        "absolute_score_cost_vs_feasible_benchmark"
    ] / max(abs(float(benchmark_row["robust_score"])), 1e-12)
    comparison_path = run_dir / "coverage_vs_feasible_information.csv"
    _write_csv(comparison_path, _tag_frame(comparison, state), outputs)
    strategy_comparison_path = run_dir / "coverage_strategy_comparison.csv"
    _write_csv(strategy_comparison_path, _tag_frame(comparison, state), outputs)
    drying_path = run_dir / "drying_margin_by_policy_member.csv"
    _write_csv(drying_path, _tag_frame(drying, state), outputs)
    physical_path = run_dir / "physical_temperature_envelope_by_candidate.csv"
    _write_csv(physical_path, _tag_frame(physical, state), outputs)
    for name, frame in (
        ("post_search_fim_validation_cases.csv", fim_validation["cases"]),
        ("post_search_fd_stability.csv", fim_validation["fd"]),
        ("post_search_grid_stability.csv", fim_validation["grid"]),
        ("post_search_direct_vs_cached_fim.csv", fim_validation["direct"]),
    ):
        _write_csv(run_dir / name, _tag_frame(frame, state), outputs)
    _write_json(
        run_dir / "post_search_fim_validation.json",
        _artifact_payload(state, fim_validation["gate"]),
        outputs,
    )
    schedules = pd.concat(
        [row["schedule_frame"] for row in sampling["all_results"]],
        ignore_index=True,
    )
    _write_csv(
        run_dir / "optimized_sampling_schedules.csv",
        _tag_frame(schedules, state),
        outputs,
    )
    sampling_comparison = pd.DataFrame(
        [row["row"] for row in sampling["all_results"]]
    )
    _write_csv(
        run_dir / "sampling_strategy_comparison.csv",
        _tag_frame(sampling_comparison, state),
        outputs,
    )
    captures = pd.concat(
        [row["captures"] for row in sampling["all_results"]], ignore_index=True
    )
    _write_csv(
        run_dir / "capture_interval_summary.csv",
        _tag_frame(captures, state),
        outputs,
    )
    for name, frame in (
        ("parameter_information_comparison.csv", sampling["parameter_information"]),
        ("posterior_variance_ratio_by_strategy.csv", sampling["variance_ratios"]),
        ("posterior_correlation_summary.csv", sampling["correlations"]),
        ("fim_eigenvalue_comparison.csv", sampling["eigenvalues"]),
    ):
        _write_csv(run_dir / name, _tag_frame(frame, state), outputs)
    _write_csv(
        run_dir / "sampling_search_restarts.csv",
        _tag_frame(sampling["restarts"], state),
        outputs,
    )
    action_frames = []
    for strategy, candidate_id, policies in (
        (BALANCED_STRATEGY, state["candidate_id"], final_policies),
        (BENCHMARK_STRATEGY, benchmark["selected"]["candidate_id"], tuple(benchmark["selected"]["policies"])),
    ):
        action_frames.append(
            _candidate_action_rows(
                strategy,
                BALANCED_DOSE,
                candidate_id,
                policies,
                model_config,
                coverage_config,
            )
        )
    _write_csv(
        run_dir / "coverage_candidate_actions.csv",
        _tag_frame(pd.concat(action_frames, ignore_index=True), state),
        outputs,
    )
    gate_summary_path = run_dir / "final_gate_summary.json"
    _write_json(
        gate_summary_path,
        {
            "candidate_id": state["candidate_id"],
            "canonical_policy_hash": state["canonical_policy_hash"],
            "final_verdict": state["final_verdict"],
            "gates": state["gates"],
            "watermark": WATERMARK,
            **AUTHORIZATION_FLAGS,
        },
        outputs,
    )
    selected_schedule = sampling["balanced_selected"]["schedule_frame"]
    figure_paths, visual_qa = _generate_figures(
        run_dir=run_dir,
        state=state,
        comparison=comparison,
        drying=drying,
        physical=physical,
        variance_ratios=sampling["variance_ratios"],
        eigenvalues=sampling["eigenvalues"],
        sampling_schedule=selected_schedule,
    )
    outputs.extend(figure_paths)
    visual_path = run_dir / "visual_qa.json"
    _write_json(visual_path, visual_qa, outputs)
    package_files: list[Path] = []
    package_manifest = None
    if state["scientific_qualification_pass"] and visual_qa["verdict"] == "PASS":
        package_files = _build_protocol_package(
            run_dir=run_dir,
            state=state,
            sampling_selected=sampling["balanced_selected"],
            capture_selected=sampling["balanced_selected"]["captures"],
            drying=drying,
            figure_paths=figure_paths,
        )
        outputs.extend(package_files)
        package_manifest = next(
            path for path in package_files if path.name.startswith("package_manifest")
        )
    audit = _cross_artifact_consistency_audit(
        run_dir=run_dir, state=state, package_manifest=package_manifest
    )
    audit_path = run_dir / "cross_artifact_consistency_audit.json"
    _write_json(audit_path, audit, outputs)
    stage_times["sampling_fim_figures_package_seconds"] = time.perf_counter() - post_started
    total_runtime = time.perf_counter() - started
    runtime = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_seconds": total_runtime,
        "runtime_hours": total_runtime / 3600.0,
        "maximum_total_hours": float(final_config["runtime"]["maximum_total_hours"]),
        "within_total_budget": total_runtime
        <= float(final_config["runtime"]["maximum_total_hours"]) * 3600.0,
        "stage_runtime_seconds": stage_times,
        "balanced_checkpoint_count": len(
            list(filesystem_path(run_dir / "balanced_search_checkpoints").rglob("*.json.gz"))
        ),
        "benchmark_checkpoint_count": len(
            list(filesystem_path(run_dir / "feasible_information_search_checkpoints").rglob("*.json.gz"))
        ),
        "balanced_resumed_seed_count": sum(
            row["resumed_from"] is not None for row in balanced["seed_results"]
        ),
        "benchmark_resumed_seed_count": sum(
            row["resumed_from"] is not None for row in benchmark["seed_results"]
        ),
        "watermark": WATERMARK,
        **AUTHORIZATION_FLAGS,
    }
    runtime_path = run_dir / "runtime_summary.json"
    _write_json(runtime_path, runtime, outputs)
    run_checks = {
        "scientific_qualification_pass": bool(state["scientific_qualification_pass"]),
        "cross_artifact_consistency_pass": audit["verdict"] == "PASS",
        "visual_qa_pass": visual_qa["verdict"] == "PASS",
        "source_hashes_verified": dependencies["all_source_output_hashes_verified"] is True,
        "runtime_within_8h": runtime["within_total_budget"] is True,
        "physical_flags_closed": AUTHORIZATION_FLAGS["physical_execution_authorized"] is False,
        "tank_assignments_empty": AUTHORIZATION_FLAGS["tank_assignments"] == [],
    }
    manifest_status = "completed" if gate_verdict(run_checks, tuple(run_checks)) == "PASS" else "failed"
    unique_outputs = []
    seen_outputs: set[str] = set()
    for path in outputs:
        key = str(path.resolve())
        if key not in seen_outputs:
            unique_outputs.append(path)
            seen_outputs.add(key)
    manifest = build_manifest(
        run_dir=run_dir,
        stage=final_config["result_stage"],
        config={
            "final": final_config,
            "coverage": coverage_config,
            "model_overlay": model_config,
        },
        sources={
            "final_config": FINAL_CONFIG_PATH,
            "coverage_config": COVERAGE_CONFIG_PATH,
            "model_config": MODEL_CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
            "adapter_manifest": runs["adapter"] / "run_manifest.json",
            "ensemble_manifest": runs["joint_ensemble"] / "run_manifest.json",
            "aroma_manifest": runs["aroma_calibration"] / "run_manifest.json",
            "actuator_manifest": runs["temperature_actuator"] / "run_manifest.json",
            "information_reference_manifest": runs["information_reference"] / "run_manifest.json",
            "sampling_reference_manifest": runs["sampling_reference"] / "run_manifest.json",
        },
        code_paths=list(DEPENDENCY_PATHS),
        random_seeds=[
            *[int(value) for value in final_config["balanced_search"]["seeds"]],
            *[
                int(value)
                for value in final_config["feasible_information_benchmark"]["seeds"]
            ],
        ],
        status=manifest_status,
        convergence=balanced["convergence"],
        gate={
            "verdict": gate_verdict(run_checks, tuple(run_checks)),
            "checks": run_checks,
            "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
        },
        outputs=unique_outputs,
        git_snapshot=git_snapshot,
    )
    manifest["final_verdict"] = state["final_verdict"]
    manifest["candidate_id"] = state["candidate_id"]
    manifest["canonical_policy_hash"] = state["canonical_policy_hash"]
    manifest["source_output_verification"] = source_verification
    manifest["partition_surrogate_provenance"] = partition_provenance
    manifest["scientific_dependency_sha256"] = dependency_hash
    manifest["watermark"] = WATERMARK
    manifest.update(AUTHORIZATION_FLAGS)
    manifest_path = run_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    print(
        json.dumps(
            {
                "run_dir": relative_or_absolute(run_dir),
                "manifest_status": manifest_status,
                "manifest_gate": manifest["gate"]["verdict"],
                "final_verdict": state["final_verdict"],
                "candidate_id": state["candidate_id"],
                "canonical_policy_hash": state["canonical_policy_hash"],
                **AUTHORIZATION_FLAGS,
            },
            ensure_ascii=False,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
