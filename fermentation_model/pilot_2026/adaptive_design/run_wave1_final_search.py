from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
import time
from itertools import combinations
from pathlib import Path
from typing import Any

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
    candidate_operational_metrics,
    policy_distance,
    practical_convergence,
    select_final_candidate,
    selected_candidate_local_qualified,
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
from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    anchor_policy,
    fim_from_prepared,
    load_json,
    load_wave1_config,
    prepare_design,
    prior_precision,
    representative_members,
    robust_information_metrics,
    vector_bounds,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    capture_git_state,
    create_immutable_run_directory,
    filesystem_path,
    sha256_file,
    sha256_payload,
    verify_manifest_output,
    write_json,
    write_json_atomic,
)
from pilot_2026.adaptive_design.run_wave1_hybrid_search import (  # noqa: E402
    _continuous_pair_indices,
    _decode_pair,
    _evaluation_metrics,
    _rows_for_policy,
    _seed_pairs,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _gzip_json(path: Path, payload: dict[str, Any]) -> None:
    filesystem_path(path.parent).mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(filesystem_path(temporary), "wt", encoding="utf-8") as stream:
        json.dump(payload, stream, sort_keys=True, separators=(",", ":"))
    filesystem_path(temporary).replace(filesystem_path(path))


def _read_gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(filesystem_path(path), "rt", encoding="utf-8") as stream:
        return json.load(stream)


def _relative_fim_difference(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right, ord="fro") / max(np.linalg.norm(right, ord="fro"), 1e-12))


def _candidate_id(canonical: np.ndarray) -> str:
    return sha256_payload(np.asarray(canonical, dtype=float).tolist())[:16]


def _common_completed_continuation_windows(
    states: dict[int, CheckpointedSwarmState], base_iteration: int, window: int
) -> int:
    """Return only extension windows completed by every independent seed."""

    if window <= 0:
        raise ValueError("Continuation window must be positive")
    return min(
        max((state.iteration - base_iteration) // window, 0)
        for state in states.values()
    )


def _history_from_checkpoints(
    checkpoint_dir: Path, maximum_iterations: int, stage2_iterations: int
) -> list[dict[str, Any]]:
    """Reconstruct the complete multifidelity trajectory for resumed runs."""

    rows: list[dict[str, Any]] = []
    by_seed: dict[int, list[dict[str, Any]]] = {}
    for path in checkpoint_dir.glob("*.json.gz"):
        payload = _read_gzip_json(path)
        swarm = payload["swarm"]
        seed = int(swarm["seed"])
        by_seed.setdefault(seed, []).append(
            {
                "checkpoint_sequence": int(payload["checkpoint_sequence"]),
                "iteration": int(swarm["iteration"]),
                "fidelity": str(swarm["fidelity"]),
                "global_best_value": float(swarm["global_best_value"]),
            }
        )
    for seed, checkpoints in sorted(by_seed.items()):
        previous_by_fidelity: dict[str, float] = {}
        for item in sorted(checkpoints, key=lambda row: row["checkpoint_sequence"]):
            fidelity = str(item["fidelity"])
            objective = float(item["global_best_value"])
            previous = previous_by_fidelity.get(fidelity, objective)
            iteration = int(item["iteration"])
            if fidelity == "four_member":
                stage = "stage1_four_member_exploration"
            elif iteration <= maximum_iterations:
                stage = "stage2_eight_member_rescore"
            elif iteration <= maximum_iterations + stage2_iterations:
                stage = "stage2_eight_member_continuation"
            else:
                stage = "stage5_eight_member_extension"
            rows.append(
                {
                    "independent_seed": seed,
                    "stage": stage,
                    "iteration": iteration,
                    "robust_score": -objective,
                    "improvement": max(previous - objective, 0.0),
                }
            )
            previous_by_fidelity[fidelity] = objective
    return rows


def _full_candidate(
    candidate: dict[str, Any],
    *,
    anchor: DesignPolicy,
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
    reference_evaluations: list | None = None,
) -> dict[str, Any]:
    from pilot_2026.adaptive_design.pilot_mbdoe_adapter import evaluate_campaign

    policies = (anchor, *candidate["pair"])
    score, evaluations = evaluate_campaign(
        policies,
        ensemble,
        list(range(len(ensemble))),
        prior,
        config,
        partitions,
        design_cache=design_cache,
    )
    metrics = _evaluation_metrics(evaluations, prior, config, reference=reference_evaluations)
    margin_h = float(config["operations"]["minimum_action_to_drying_margin_h"])
    margin_probability = float(
        np.mean(
            [
                row.minimum_action_to_drying_margin_h >= margin_h - 1e-9
                for row in evaluations
            ]
        )
    )
    minimum_margin = float(
        min(row.minimum_action_to_drying_margin_h for row in evaluations)
    )
    operational = candidate_operational_metrics(policies, config, margin_h=minimum_margin)
    feasible = bool(
        metrics["completion_probability"] >= float(config["completion"]["minimum_probability"])
        and margin_probability >= float(config["completion"]["minimum_probability"])
        and operational["temperature_feasible"]
        and operational["nutrition_feasible"]
    )
    return {
        **candidate,
        "policies": policies,
        "full_score": float(score),
        "full_objective": -float(score),
        "full_evaluations": evaluations,
        "full_metrics": metrics,
        "action_margin_probability": margin_probability,
        "feasible": feasible,
        **operational,
    }


def _flat_candidate(row: dict[str, Any]) -> dict[str, Any]:
    metrics = row.get("full_metrics", {})
    return {
        "candidate_id": row["candidate_id"],
        "canonical_policy_hash": row["canonical_hash"],
        "source_type": row["source_type"],
        "source_types": ";".join(sorted(row.get("source_types", {row["source_type"]}))),
        "source_seed": row.get("source_seed"),
        "source_seeds": ";".join(str(value) for value in sorted(row.get("source_seeds", set()))),
        "parent_candidate_id": row.get("parent_candidate_id"),
        "four_member_robust_score": -float(row.get("four_objective", math.nan)),
        "eight_member_robust_score": -float(row.get("eight_objective", math.nan)),
        "full_ensemble_robust_score": float(row["full_score"]),
        "full_objective": float(row["full_objective"]),
        "full_improvement": row.get("full_improvement"),
        "feasible": bool(row["feasible"]),
        "action_margin_probability": float(row["action_margin_probability"]),
        "minimum_action_margin_h": float(row["minimum_action_margin_h"]),
        "maximum_temperature_jump_c": float(row["maximum_temperature_jump_c"]),
        "total_thermal_variation_c": float(row["total_thermal_variation_c"]),
        "temperature_changes": int(row["temperature_changes"]),
        "total_yan_mg_l": float(row["total_yan_mg_l"]),
        **{key: value for key, value in metrics.items()},
    }


def _actuator_validation(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    prior: np.ndarray,
    config: dict[str, Any],
    partitions: dict[str, Any],
    design_cache: dict[tuple, tuple[np.ndarray, float, float]],
) -> tuple[pd.DataFrame, bool]:
    from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (
        actuator_temperature_trajectory,
        evaluate_campaign,
    )

    margin_h = float(config["operations"]["minimum_action_to_drying_margin_h"])
    required = float(config["completion"]["minimum_probability"])
    rows = []
    scenarios = approved_actuator_scenarios(config)
    for index, scenario in enumerate(scenarios, start=1):
        numeric = {key: float(value) for key, value in scenario.items() if key != "scenario"}
        score, evaluations = evaluate_campaign(
            policies,
            ensemble,
            list(range(len(ensemble))),
            prior,
            config,
            partitions,
            actuator_scenario=numeric,
            design_cache=design_cache,
        )
        metrics = _evaluation_metrics(evaluations, prior, config)
        margin_probability = float(
            np.mean(
                [
                    row.minimum_action_to_drying_margin_h >= margin_h - 1e-9
                    for row in evaluations
                ]
            )
        )
        trajectories = [actuator_temperature_trajectory(policy, config, numeric) for policy in policies]
        physical_min = min(float(np.min(row["physical_temperature_c"])) for row in trajectories)
        physical_max = max(float(np.max(row["physical_temperature_c"])) for row in trajectories)
        probe_min = min(float(np.min(row["probe_reading_c"])) for row in trajectories)
        probe_max = max(float(np.max(row["probe_reading_c"])) for row in trajectories)
        feasible = bool(
            metrics["completion_probability"] >= required and margin_probability >= required
        )
        rows.append(
            {
                **scenario,
                "ensemble_members": len(ensemble),
                "robust_score": float(score),
                **metrics,
                "action_margin_probability": margin_probability,
                "minimum_physical_temperature_c": physical_min,
                "maximum_physical_temperature_c": physical_max,
                "minimum_probe_reading_c": probe_min,
                "maximum_probe_reading_c": probe_max,
                "probe_bias_model": config["future_process"]["temperature_actuator"][
                    "probe_bias_control_approximation"
                ],
                "feasible": feasible,
            }
        )
        print(
            json.dumps(
                {
                    "stage": "actuator_envelope",
                    "scenario": index,
                    "total": len(scenarios),
                    "feasible": feasible,
                }
            ),
            flush=True,
        )
    frame = pd.DataFrame(rows)
    complete = bool(
        len(frame) == 243
        and frame["ensemble_members"].eq(64).all()
        and frame["feasible"].all()
    )
    return frame, complete


def _post_search_fim(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    selected_evaluations: list,
    actuator: pd.DataFrame,
    config: dict[str, Any],
    partitions: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    representative = representative_members(ensemble, config)
    gains = np.asarray([row.information_gain for row in selected_evaluations], dtype=float)
    drying = np.asarray([row.drying_time_h for row in selected_evaluations], dtype=float)
    member_cases = {
        "central_member": int(representative[0]),
        "lowest_information_member": int(np.argmin(gains)),
        "highest_information_member": int(np.argmax(gains)),
        "critical_drying_member": int(np.argmin(drying)),
        "worst_gain_member": int(np.argmin(gains)),
    }
    worst_actuator = actuator.sort_values("robust_score").iloc[0]
    extreme_scenario = {
        key: float(worst_actuator[key])
        for key in (
            "tau_h",
            "tracking_error_c",
            "initial_temperature_offset_c",
            "probe_bias_c",
            "command_delay_h",
        )
    }
    cases = []
    fd_rows = []
    grid_rows = []
    case_index = 0
    for policy in policies:
        for member_label, member_index in member_cases.items():
            case_index += 1
            scenario = extreme_scenario if member_label == "critical_drying_member" else None
            prepared_by_grid = {
                grid: prepare_design(
                    policy,
                    ensemble.iloc[member_index],
                    config,
                    partitions,
                    simulation_grid_step_h=float(grid),
                    actuator_scenario=scenario,
                )
                for grid in config["fim_validation"]["time_grid_steps_h"]
            }
            if any(value is None for value in prepared_by_grid.values()):
                cases.append(
                    {
                        "case_id": case_index,
                        "policy": policy.name,
                        "member_case": member_label,
                        "ensemble_member": member_index,
                        "actuator_extreme": scenario is not None,
                        "passed": False,
                        "reason": "simulation_failed",
                    }
                )
                continue
            nominal_grid = float(config["fim_validation"]["nominal_time_grid_step_h"])
            nominal_prepared = prepared_by_grid[nominal_grid]
            nominal_step = float(config["objective"]["finite_difference_log_step"])
            nominal_fim = fim_from_prepared(
                nominal_prepared,
                np.asarray(config["sampling"]["preliminary_times_h"], dtype=float),
                config,
                finite_difference_log_step=nominal_step,
            )
            case_fd = []
            for step in config["fim_validation"]["log_steps"]:
                fim = fim_from_prepared(
                    nominal_prepared,
                    np.asarray(config["sampling"]["preliminary_times_h"], dtype=float),
                    config,
                    finite_difference_log_step=float(step),
                )
                difference = _relative_fim_difference(fim, nominal_fim)
                case_fd.append(difference)
                fd_rows.append(
                    {
                        "case_id": case_index,
                        "policy": policy.name,
                        "member_case": member_label,
                        "ensemble_member": member_index,
                        "finite_difference_log_step": float(step),
                        "relative_fim_difference_vs_0_02": difference,
                        "fim_logdet": float(np.linalg.slogdet(fim + np.eye(9) * 1e-12)[1]),
                    }
                )
            case_grid = []
            for grid, prepared in prepared_by_grid.items():
                fim = fim_from_prepared(
                    prepared,
                    np.asarray(config["sampling"]["preliminary_times_h"], dtype=float),
                    config,
                    finite_difference_log_step=nominal_step,
                )
                difference = _relative_fim_difference(fim, nominal_fim)
                case_grid.append(difference)
                grid_rows.append(
                    {
                        "case_id": case_index,
                        "policy": policy.name,
                        "member_case": member_label,
                        "ensemble_member": member_index,
                        "time_grid_step_h": float(grid),
                        "relative_fim_difference_vs_2h": difference,
                        "fim_logdet": float(np.linalg.slogdet(fim + np.eye(9) * 1e-12)[1]),
                    }
                )
            fd_pass = max(case_fd) <= float(
                config["fim_validation"]["maximum_relative_step_fim_difference"]
            )
            grid_pass = max(case_grid) <= float(
                config["fim_validation"]["maximum_relative_grid_fim_difference"]
            )
            cases.append(
                {
                    "case_id": case_index,
                    "policy": policy.name,
                    "member_case": member_label,
                    "ensemble_member": member_index,
                    "actuator_extreme": scenario is not None,
                    "maximum_fd_relative_difference": max(case_fd),
                    "maximum_grid_relative_difference": max(case_grid),
                    "fd_pass": fd_pass,
                    "grid_pass": grid_pass,
                    "passed": bool(fd_pass and grid_pass),
                    "reason": None,
                }
            )
    cases_frame = pd.DataFrame(cases)
    checks = {
        "exact_final_anchor_A_B_covered": set(cases_frame["policy"]) == {
            policy.name for policy in policies
        },
        "fd_steps_0_01_0_02_0_04_covered": set(pd.DataFrame(fd_rows)["finite_difference_log_step"])
        == {0.01, 0.02, 0.04},
        "grids_4_2_1h_covered": set(pd.DataFrame(grid_rows)["time_grid_step_h"])
        == {1.0, 2.0, 4.0},
        "critical_member_cases_covered": set(cases_frame["member_case"])
        == set(member_cases),
        "approved_actuator_extreme_covered": bool(cases_frame["actuator_extreme"].any()),
        "all_cases_pass": bool(cases_frame["passed"].all()),
    }
    gate = {
        "gate": "post_search_fim_validation",
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "thresholds": {
            "maximum_relative_step_fim_difference": config["fim_validation"][
                "maximum_relative_step_fim_difference"
            ],
            "maximum_relative_grid_fim_difference": config["fim_validation"][
                "maximum_relative_grid_fim_difference"
            ],
        },
        "worst_approved_actuator_scenario": extreme_scenario,
    }
    return gate, cases_frame, pd.DataFrame(fd_rows), pd.DataFrame(grid_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run final checkpointed multifidelity Wave-1 search")
    parser.add_argument("--source-adapter-run", required=True)
    parser.add_argument("--source-ensemble-run", required=True)
    parser.add_argument("--source-aroma-run", required=True)
    parser.add_argument("--source-actuator-run", required=True)
    parser.add_argument("--resume-run")
    parser.add_argument("--maximum-wall-clock-hours", type=float)
    return parser.parse_args()


def main() -> None:
    from pilot_2026.adaptive_design.pilot_mbdoe_adapter import evaluate_campaign

    args = parse_args()
    total_started = time.perf_counter()
    search_started = time.perf_counter()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    aroma_config = load_json(AROMA_CONFIG_PATH)
    adapter_run = _run_path(args.source_adapter_run)
    ensemble_run = _run_path(args.source_ensemble_run)
    aroma_run = _run_path(args.source_aroma_run)
    actuator_run = _run_path(args.source_actuator_run)
    source_verification = {
        "adapter": verify_manifest_output(adapter_run, "adapter_gate.json"),
        "ensemble": verify_manifest_output(ensemble_run, "joint_parameter_ensemble.csv"),
        "aroma": verify_manifest_output(aroma_run, "aroma_calibration_gate.json"),
        "actuator": verify_manifest_output(actuator_run, "actuator_estimates_by_run.csv"),
    }
    adapter_gate = load_json(adapter_run / "adapter_gate.json")
    if adapter_gate["verdict"] != "PASS":
        raise RuntimeError("Source corrected adapter gate is not PASS")
    ensemble = pd.read_csv(filesystem_path(ensemble_run / "joint_parameter_ensemble.csv"))
    if len(ensemble) != int(config["search"]["full_ensemble_members"]):
        raise ValueError("The final search requires exactly the approved 64-member ensemble")
    partitions, partition_provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    representatives = representative_members(ensemble, config)
    four_members = representatives[: int(config["search"]["stage1_members"])]
    eight_members = representatives[: int(config["search"]["stage2_members"])]
    full_members = list(range(len(ensemble)))
    anchor = anchor_policy(config)
    single_bounds = vector_bounds(config)
    bounds = np.vstack([single_bounds, single_bounds])
    design_cache: dict[tuple, tuple[np.ndarray, float, float]] = {}
    objective_cache: dict[str, float] = {}

    def objective(values: np.ndarray, members: list[int]) -> float:
        pair, canonical, _ = _decode_pair(values, config, "objective")
        key = f"{','.join(str(value) for value in members)}|{sha256_payload(canonical.tolist())}"
        if key not in objective_cache:
            score, _ = evaluate_campaign(
                (anchor, *pair),
                ensemble,
                members,
                prior,
                config,
                partitions,
                design_cache=design_cache,
            )
            objective_cache[key] = -float(score)
        return objective_cache[key]

    four_objective = lambda values: objective(values, four_members)
    eight_objective = lambda values: objective(values, eight_members)
    full_objective = lambda values: objective(values, full_members)
    if args.resume_run:
        run_dir = _run_path(args.resume_run)
        if not filesystem_path(run_dir).is_dir():
            raise FileNotFoundError(run_dir)
    else:
        run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_final_search", config)
    checkpoint_dir = run_dir / "checkpoints"
    filesystem_path(checkpoint_dir).mkdir(parents=True, exist_ok=True)
    config_hashes = {
        "wave1_mbdoe_config_sha256": sha256_file(CONFIG_PATH),
        "design_constraints_sha256": sha256_file(CONSTRAINTS_PATH),
        "source_adapter_manifest_sha256": sha256_file(adapter_run / "run_manifest.json"),
        "source_ensemble_manifest_sha256": sha256_file(ensemble_run / "run_manifest.json"),
    }
    max_hours = float(
        args.maximum_wall_clock_hours
        if args.maximum_wall_clock_hours is not None
        else config["search"]["maximum_wall_clock_h"]
    )
    if max_hours <= 0.0 or max_hours > 8.0 + 1e-12:
        raise ValueError("Final search wall-clock limit must be in (0, 8] hours")
    existing_checkpoints = list(checkpoint_dir.glob("*.json.gz"))
    prior_search_elapsed = (
        max(filesystem_path(path).stat().st_mtime for path in existing_checkpoints)
        - filesystem_path(run_dir).stat().st_ctime
        if args.resume_run and existing_checkpoints
        else 0.0
    )
    prior_search_elapsed = max(float(prior_search_elapsed), 0.0)
    remaining_search_seconds = max(max_hours * 3600.0 - prior_search_elapsed, 0.0)
    deadline = search_started + remaining_search_seconds
    seed_profiles = _seed_pairs(config)
    allowed_seed_count = max(
        1,
        min(
            len(seed_profiles),
            int(
                round(
                    float(config["search"]["seed_profile_fraction"])
                    * int(config["search"]["particles"])
                )
            ),
        ),
    )
    seed_profiles = seed_profiles[:allowed_seed_count]
    states: dict[int, CheckpointedSwarmState] = {}
    checkpoint_paths: list[Path] = []
    history_rows: list[dict[str, Any]] = (
        _history_from_checkpoints(
            checkpoint_dir,
            int(config["search"]["maximum_iterations"]),
            int(config["search"]["stage2_continuation_iterations"]),
        )
        if args.resume_run
        else []
    )
    sequence_by_seed: dict[int, int] = {}

    def save_checkpoint(state: CheckpointedSwarmState, candidates: list[str]) -> None:
        sequence = sequence_by_seed.get(state.seed, 0) + 1
        sequence_by_seed[state.seed] = sequence
        path = checkpoint_dir / (
            f"seed_{state.seed}_{sequence:04d}_{state.fidelity}_iteration_{state.iteration:04d}.json.gz"
        )
        payload = {
            "checkpoint_schema_version": 1,
            "checkpoint_sequence": sequence,
            "swarm": state.to_payload(),
            "evaluation_cache": objective_cache,
            "candidate_ids": candidates,
            "config_hashes": config_hashes,
            "numpy_version": np.__version__,
            "cumulative_search_elapsed_seconds": prior_search_elapsed
            + float(time.perf_counter() - search_started),
        }
        _gzip_json(path, payload)
        checkpoint_paths.append(path)
        print(
            json.dumps(
                {
                    "checkpoint": path.name,
                    "seed": state.seed,
                    "fidelity": state.fidelity,
                    "iteration": state.iteration,
                    "best_score": -state.global_best_value,
                }
            ),
            flush=True,
        )

    def latest_checkpoint(seed: int) -> Path | None:
        candidates = sorted(checkpoint_dir.glob(f"seed_{seed}_*.json.gz"))
        return candidates[-1] if candidates else None

    search_timed_out = False
    for seed in (int(value) for value in config["independent_seeds"]):
        latest = latest_checkpoint(seed) if args.resume_run else None
        if latest is not None:
            payload = _read_gzip_json(latest)
            if payload["config_hashes"] != config_hashes:
                raise RuntimeError("Checkpoint configuration/source hashes do not match")
            state = CheckpointedSwarmState.from_payload(payload["swarm"])
            sequence_by_seed[seed] = int(payload["checkpoint_sequence"])
            objective_cache.update(
                {str(key): float(value) for key, value in payload["evaluation_cache"].items()}
            )
            sequence_by_seed[seed] = int(payload["checkpoint_sequence"])
            stage1_iterations = int(config["search"]["maximum_iterations"])
            for history_index, score in enumerate(state.history[1:], start=1):
                if history_index <= stage1_iterations:
                    stage = "stage1_four_member"
                    iteration = history_index
                elif history_index == stage1_iterations + 1:
                    continue  # fidelity rescore at the same swarm iteration
                else:
                    stage = "stage2_eight_member_continuation"
                    iteration = stage1_iterations + (
                        history_index - stage1_iterations - 1
                    )
                history_rows.append(
                    {
                        "independent_seed": seed,
                        "stage": stage,
                        "iteration": iteration,
                        "robust_score": -float(score),
                        "improvement": float(state.improvement_history[history_index]),
                    }
                )
        else:
            state = initialize_checkpointed_swarm(
                four_objective,
                bounds,
                particles=int(config["search"]["particles"]),
                seed=seed,
                fidelity="four_member",
                initial_positions=seed_profiles,
            )
            save_checkpoint(state, [_candidate_id(_decode_pair(state.global_best_position, config, "checkpoint")[1])])
        while (
            state.fidelity == "four_member"
            and state.iteration < int(config["search"]["maximum_iterations"])
        ):
            if time.perf_counter() >= deadline:
                search_timed_out = True
                break
            previous = state.global_best_value
            state = advance_checkpointed_swarm(
                state,
                four_objective,
                bounds,
                inertia=float(config["search"]["inertia"]),
                cognitive=float(config["search"]["cognitive"]),
                social=float(config["search"]["social"]),
            )
            history_rows.append(
                {
                    "independent_seed": seed,
                    "stage": "stage1_four_member",
                    "iteration": state.iteration,
                    "robust_score": -state.global_best_value,
                    "improvement": max(previous - state.global_best_value, 0.0),
                }
            )
            save_checkpoint(state, [_candidate_id(_decode_pair(state.global_best_position, config, "checkpoint")[1])])
        if search_timed_out:
            states[seed] = state
            break
        if state.fidelity == "four_member":
            state = rescore_checkpointed_swarm(
                state, eight_objective, fidelity="eight_member"
            )
            save_checkpoint(state, [_candidate_id(_decode_pair(state.global_best_position, config, "checkpoint")[1])])
        target = int(config["search"]["maximum_iterations"]) + int(
            config["search"]["stage2_continuation_iterations"]
        )
        while state.iteration < target:
            if time.perf_counter() >= deadline:
                search_timed_out = True
                break
            previous = state.global_best_value
            state = advance_checkpointed_swarm(
                state,
                eight_objective,
                bounds,
                inertia=float(config["search"]["inertia"]),
                cognitive=float(config["search"]["cognitive"]),
                social=float(config["search"]["social"]),
            )
            history_rows.append(
                {
                    "independent_seed": seed,
                    "stage": "stage2_eight_member_continuation",
                    "iteration": state.iteration,
                    "robust_score": -state.global_best_value,
                    "improvement": max(previous - state.global_best_value, 0.0),
                }
            )
            save_checkpoint(state, [_candidate_id(_decode_pair(state.global_best_position, config, "checkpoint")[1])])
        states[seed] = state
        if search_timed_out:
            break

    completed_seed_count = len(states)
    continuation_window = int(config["search"]["stage5_continuation_window"])
    continuation_improvement_fraction = math.inf
    if completed_seed_count == len(config["independent_seeds"]) and not search_timed_out:
        base_iteration = int(config["search"]["maximum_iterations"]) + int(
            config["search"]["stage2_continuation_iterations"]
        )

        def checkpoint_global_best(seed: int, iteration: int) -> float:
            state = states[seed]
            if state.iteration == iteration:
                return float(state.global_best_value)
            candidates = sorted(
                checkpoint_dir.glob(
                    f"seed_{seed}_*_{state.fidelity}_iteration_{iteration:04d}.json.gz"
                )
            )
            if not candidates:
                raise RuntimeError(
                    f"Missing checkpoint for seed {seed}, iteration {iteration}"
                )
            payload = _read_gzip_json(candidates[-1])
            return float(payload["swarm"]["global_best_value"])

        completed_windows = _common_completed_continuation_windows(
            states, base_iteration, continuation_window
        )
        if completed_windows:
            window_end = base_iteration + completed_windows * continuation_window
            before = min(
                checkpoint_global_best(seed, window_end - continuation_window)
                for seed in states
            )
            after = min(
                checkpoint_global_best(seed, window_end) for seed in states
            )
            continuation_improvement_fraction = max(before - after, 0.0) / max(
                abs(after), 1e-12
            )
        tolerance = float(
            config["search"]["continuation_improvement_tolerance_fraction"]
        )
        while (
            time.perf_counter() < deadline
            and continuation_improvement_fraction > tolerance
        ):
            target_iteration = base_iteration + (completed_windows + 1) * continuation_window
            for seed, state in list(states.items()):
                while state.iteration < target_iteration:
                    if time.perf_counter() >= deadline:
                        search_timed_out = True
                        break
                    previous = state.global_best_value
                    state = advance_checkpointed_swarm(
                        state,
                        eight_objective,
                        bounds,
                        inertia=float(config["search"]["inertia"]),
                        cognitive=float(config["search"]["cognitive"]),
                        social=float(config["search"]["social"]),
                    )
                    history_rows.append(
                        {
                            "independent_seed": seed,
                            "stage": "stage5_eight_member_extension",
                            "iteration": state.iteration,
                            "robust_score": -state.global_best_value,
                            "improvement": max(previous - state.global_best_value, 0.0),
                        }
                    )
                    save_checkpoint(
                        state,
                        [
                            _candidate_id(
                                _decode_pair(
                                    state.global_best_position, config, "checkpoint"
                                )[1]
                            )
                        ],
                    )
                states[seed] = state
                if search_timed_out:
                    break
            if search_timed_out:
                break
            before = min(
                checkpoint_global_best(seed, target_iteration - continuation_window)
                for seed in states
            )
            after = min(
                checkpoint_global_best(seed, target_iteration) for seed in states
            )
            continuation_improvement_fraction = max(before - after, 0.0) / max(
                abs(after), 1e-12
            )
            completed_windows += 1

    # Preserve at least one canonical champion from every seed before global ranking.
    candidate_pool: dict[str, dict[str, Any]] = {}

    def register(values: np.ndarray, source_type: str, source_seed: int) -> dict[str, Any]:
        pair, canonical, repairs = _decode_pair(values, config, f"seed_{source_seed}")
        identifier = _candidate_id(canonical)
        candidate = candidate_pool.setdefault(
            identifier,
            {
                "candidate_id": identifier,
                "canonical_hash": sha256_payload(canonical.tolist()),
                "canonical": canonical,
                "pair": pair,
                "repairs": repairs,
                "source_type": source_type,
                "source_types": set(),
                "source_seed": source_seed,
                "source_seeds": set(),
            },
        )
        candidate["source_types"].add(source_type)
        candidate["source_seeds"].add(source_seed)
        return candidate

    champions = []
    for seed, state in states.items():
        champion = register(state.global_best_position, "per_seed_champion", seed)
        champion["four_objective"] = four_objective(champion["canonical"])
        champion["eight_objective"] = eight_objective(champion["canonical"])
        champions.append(champion)
        for rank_index in np.argsort(state.personal_best_values)[: max(2, int(config["search"]["top_k_candidates"]))]:
            candidate = register(
                state.personal_best_positions[int(rank_index)],
                "multifidelity_continuation",
                seed,
            )
            candidate["four_objective"] = four_objective(candidate["canonical"])
            candidate["eight_objective"] = eight_objective(candidate["canonical"])
    candidates = list(candidate_pool.values())
    candidates.sort(key=lambda row: float(row["eight_objective"]))
    for candidate in candidates[: int(config["search"]["top_k_candidates"])]:
        candidate["source_types"].add("pso_original_top_k")
        candidate["source_type"] = "pso_original_top_k"
    # Every seed champion plus every global top-K original is full-ensemble eligible.
    eligible_ids = {
        row["candidate_id"] for row in champions
    } | {
        row["candidate_id"]
        for row in candidates[: int(config["search"]["top_k_candidates"])]
    }
    anchor_score, anchor_evaluations = evaluate_campaign(
        (anchor, anchor, anchor),
        ensemble,
        full_members,
        prior,
        config,
        partitions,
        design_cache=design_cache,
    )
    full_by_id: dict[str, dict[str, Any]] = {}
    for identifier in eligible_ids:
        full_by_id[identifier] = _full_candidate(
            candidate_pool[identifier],
            anchor=anchor,
            ensemble=ensemble,
            prior=prior,
            config=config,
            partitions=partitions,
            design_cache=design_cache,
            reference_evaluations=anchor_evaluations,
        )
    per_seed_eight_rows = []
    per_seed_full_rows = []
    for champion in champions:
        full = full_by_id[champion["candidate_id"]]
        per_seed_eight_rows.append(
            {
                "independent_seed": champion["source_seed"],
                "candidate_id": champion["candidate_id"],
                "canonical_policy_hash": champion["canonical_hash"],
                "eight_member_robust_score": -float(champion["eight_objective"]),
            }
        )
        per_seed_full_rows.append(
            {
                "independent_seed": champion["source_seed"],
                "candidate_id": champion["candidate_id"],
                "canonical_policy_hash": champion["canonical_hash"],
                "full_ensemble_robust_score": float(full["full_score"]),
                "completion_probability": float(full["full_metrics"]["completion_probability"]),
                "action_margin_probability": float(full["action_margin_probability"]),
                "feasible": bool(full["feasible"]),
            }
        )
    distance_rows = []
    champion_by_seed = {int(row["source_seed"]): row for row in champions}
    for left_seed, right_seed in combinations(sorted(champion_by_seed), 2):
        components = policy_distance(
            champion_by_seed[left_seed]["pair"],
            champion_by_seed[right_seed]["pair"],
            config,
        )
        distance_rows.append(
            {"left_seed": left_seed, "right_seed": right_seed, **components}
        )
    convergence = practical_convergence(
        per_seed_full_rows,
        distance_rows,
        continuation_improvement_fraction,
        config,
    )

    # Full-ensemble acceptance is performed at every trust-region proposal.
    local_rows = []
    rejected_rows = []
    eligible = list(full_by_id.values())
    local_results = []
    for local_rank, parent in enumerate(
        sorted(eligible, key=lambda row: row["full_objective"])[
            : int(config["search"]["top_k_local_refinement"])
        ],
        start=1,
    ):
        result = sequential_ipopt_refine(
            parent["canonical"],
            bounds,
            eight_objective,
            _continuous_pair_indices(parent["canonical"], config),
            config,
            candidate_id=f"full_top_{local_rank}",
            validation_objective=full_objective,
        )
        refined_pair, refined_canonical, refined_repairs = _decode_pair(
            result.x, config, f"refined_{local_rank}"
        )
        refined_objective = full_objective(refined_canonical)
        improvement = float(parent["full_objective"] - refined_objective)
        accepted = bool(
            result.accepted_improvement
            and improvement > float(config["search"]["improvement_tolerance"])
        )
        local_results.append(
            {
                "parent_candidate_id": parent["candidate_id"],
                "state": result.state,
                "accepted": accepted,
            }
        )
        local_rows.extend(result.trace)
        rejected_rows.extend(result.rejected_candidates)
        if accepted:
            refined = {
                "candidate_id": _candidate_id(refined_canonical),
                "canonical_hash": sha256_payload(refined_canonical.tolist()),
                "canonical": refined_canonical,
                "pair": refined_pair,
                "repairs": refined_repairs,
                "source_type": "accepted_refinement",
                "source_types": {"accepted_refinement"},
                "source_seed": parent.get("source_seed"),
                "source_seeds": set(parent.get("source_seeds", set())),
                "parent_candidate_id": parent["candidate_id"],
                "four_objective": four_objective(refined_canonical),
                "eight_objective": eight_objective(refined_canonical),
                "full_improvement": improvement,
            }
            refined_full = _full_candidate(
                refined,
                anchor=anchor,
                ensemble=ensemble,
                prior=prior,
                config=config,
                partitions=partitions,
                design_cache=design_cache,
                reference_evaluations=anchor_evaluations,
            )
            eligible.append(refined_full)
    selected, selection_reason = select_final_candidate(eligible, config)
    local_qualified, local_reason = selected_candidate_local_qualified(
        selected, eligible, float(config["search"]["improvement_tolerance"])
    )
    local_multistart_consistency = len({row["state"] for row in local_results}) <= 1
    selected_pair = (
        DesignPolicy("candidate_A", selected["pair"][0].temperature_c, selected["pair"][0].nutrition_mg_yan_l),
        DesignPolicy("candidate_B", selected["pair"][1].temperature_c, selected["pair"][1].nutrition_mg_yan_l),
    )
    selected_policies = (anchor, *selected_pair)
    selected_score, selected_evaluations = evaluate_campaign(
        selected_policies,
        ensemble,
        full_members,
        prior,
        config,
        partitions,
        design_cache=design_cache,
    )
    selected_metrics = _evaluation_metrics(
        selected_evaluations, prior, config, reference=anchor_evaluations
    )
    search_runtime = prior_search_elapsed + float(time.perf_counter() - search_started)
    actuator_frame, actuator_pass = _actuator_validation(
        selected_policies, ensemble, prior, config, partitions, design_cache
    )
    post_fim_gate, post_fim_cases, post_fd, post_grid = _post_search_fim(
        selected_policies,
        ensemble,
        selected_evaluations,
        actuator_frame,
        config,
        partitions,
    )
    checks = {
        "adapter_fim_pass": adapter_gate["verdict"] == "PASS",
        "sources_verified": len(source_verification) == 4,
        "five_seeds_completed": completed_seed_count == 5,
        "three_seed_champions_within_0_5pct": bool(
            convergence.get("checks", {}).get(
                "at_least_three_seed_champions_within_0_5pct", False
            )
        ),
        "continuation_improvement_at_most_0_1pct": bool(
            convergence.get("checks", {}).get(
                "continuation_improvement_at_most_0_1pct", False
            )
        ),
        "operational_policy_family_converged": bool(
            convergence.get("checks", {}).get(
                "operationally_equivalent_policy_family", False
            )
        ),
        "top_k_revalidated_on_64_members": sum(
            "pso_original_top_k" in row["source_types"] for row in eligible
        )
        >= int(config["search"]["top_k_candidates"]),
        "originals_and_refinements_compared": bool(
            any(row["source_type"] == "pso_original_top_k" for row in eligible)
        ),
        "selected_candidate_local_qualified": local_qualified,
        "post_search_fim_pass": post_fim_gate["verdict"] == "PASS",
        "completion_probability_at_least_95pct": selected_metrics["completion_probability"]
        >= float(config["completion"]["minimum_probability"]),
        "actions_respect_24h_drying_margin": selected["action_margin_probability"]
        >= float(config["completion"]["minimum_probability"]),
        "temperature_range_and_jumps_feasible": bool(selected["temperature_feasible"]),
        "nutrition_physically_feasible": bool(selected["nutrition_feasible"]),
        "approved_actuator_envelope_64_members_pass": actuator_pass,
        "checkpoint_every_iteration_generated": bool(
            list(checkpoint_dir.glob("*.json.gz"))
        ),
        "search_wall_clock_within_8h": search_runtime <= 8.0 * 3600.0 + 1.0,
        "profiles_for_physical_execution_false": True,
        "physical_execution_authorized_false": True,
        "executable_schedule_issued_false": True,
        "tank_assignments_empty": True,
    }
    critical = tuple(checks)
    verdict = "PASS" if all(checks[name] for name in critical) else "FAIL"
    if verdict == "PASS" and not local_multistart_consistency:
        verdict = "PASS_CONDITIONAL"
    failed = [name for name, value in checks.items() if not value]
    gate = {
        "gate": "final_wave1_multifidelity_search",
        "verdict": verdict,
        "checks": checks,
        "conditions": [f"failed_check:{name}" for name in failed]
        + ([] if local_multistart_consistency else ["local_multistart_consistency_mixed"]),
        "selected_candidate_id": selected["candidate_id"],
        "selected_candidate_local_qualification_reason": local_reason,
        "local_multistart_consistency": local_multistart_consistency,
        "selection_reason": selection_reason,
        "search_timed_out": search_timed_out,
        "maximum_wall_clock_hours": max_hours,
        "search_runtime_seconds": search_runtime,
        "full_ensemble": {
            "selected_robust_score": float(selected_score),
            "three_anchor_robust_score": float(anchor_score),
            **selected_metrics,
        },
        "post_search_fim_verdict": post_fim_gate["verdict"],
        "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
        "profiles_for_physical_execution": False,
        "physical_execution_authorized": False,
        "executable_schedule_issued": False,
        "tank_assignments": [],
    }

    paths = {
        "seed8": run_dir / "per_seed_champions_8_member.csv",
        "seed64": run_dir / "per_seed_champions_64_member.csv",
        "distance": run_dir / "cross_seed_policy_distance.csv",
        "convergence": run_dir / "practical_convergence_gate.json",
        "checkpoint": run_dir / "checkpoint_manifest.json",
        "eligible": run_dir / "eligible_candidate_comparison.csv",
        "provenance": run_dir / "selected_candidate_provenance.json",
        "actions": run_dir / "final_candidate_policy_actions.csv",
        "full": run_dir / "full_ensemble_validation.csv",
        "actuator": run_dir / "actuator_robustness_validation.csv",
        "post_gate": run_dir / "post_search_fim_validation.json",
        "post_cases": run_dir / "post_search_fim_validation_cases.csv",
        "post_fd": run_dir / "post_search_fd_stability.csv",
        "post_grid": run_dir / "post_search_grid_stability.csv",
        "local": run_dir / "local_refinement_trace.csv",
        "rejected": run_dir / "local_refinement_rejected_candidates.json",
        "history": run_dir / "multifidelity_pso_history.csv",
        "gate": run_dir / "final_search_gate.json",
        "config": run_dir / "wave1_mbdoe_config.json",
        "partition": run_dir / "partition_surrogate_provenance.json",
        "runtime": run_dir / "runtime_summary.json",
    }
    pd.DataFrame(per_seed_eight_rows).to_csv(filesystem_path(paths["seed8"]), index=False)
    pd.DataFrame(per_seed_full_rows).to_csv(filesystem_path(paths["seed64"]), index=False)
    pd.DataFrame(distance_rows).to_csv(filesystem_path(paths["distance"]), index=False)
    write_json(paths["convergence"], convergence)
    checkpoint_entries = [
        {
            "path": path.relative_to(run_dir).as_posix(),
            "sha256": sha256_file(path),
            "bytes": filesystem_path(path).stat().st_size,
        }
        for path in sorted(checkpoint_dir.glob("*.json.gz"))
    ]
    write_json(
        paths["checkpoint"],
        {
            "checkpoint_schema_version": 1,
            "checkpoint_every_iteration": True,
            "resumable": True,
            "config_hashes": config_hashes,
            "checkpoints": checkpoint_entries,
        },
    )
    pd.DataFrame([_flat_candidate(row) for row in eligible]).sort_values(
        "full_objective"
    ).to_csv(filesystem_path(paths["eligible"]), index=False)
    write_json(
        paths["provenance"],
        {
            "selected_candidate_id": selected["candidate_id"],
            "canonical_policy_hash": selected["canonical_hash"],
            "source_type": selected["source_type"],
            "source_types": sorted(selected["source_types"]),
            "source_seeds": sorted(selected["source_seeds"]),
            "parent_candidate_id": selected.get("parent_candidate_id"),
            "selection_reason": selection_reason,
            "selected_candidate_local_qualified": local_qualified,
            "selected_candidate_local_qualification_reason": local_reason,
            "local_multistart_consistency": local_multistart_consistency,
        },
    )
    pd.DataFrame(
        sum((_rows_for_policy(policy, config) for policy in selected_policies), [])
    ).to_csv(filesystem_path(paths["actions"]), index=False)
    pd.DataFrame(
        [
            {
                "ensemble_member": row.member,
                "information_gain": row.information_gain,
                "completion": row.completion,
                "residual_sugar_g_l": row.residual_sugar_g_l,
                "drying_time_h": row.drying_time_h,
                "latest_action_time_h": row.latest_action_time_h,
                "action_margin_to_drying_h": row.minimum_action_to_drying_margin_h,
                "three_anchor_information_gain": reference.information_gain,
                "paired_delta_vs_three_anchor": row.information_gain - reference.information_gain,
            }
            for row, reference in zip(selected_evaluations, anchor_evaluations)
        ]
    ).to_csv(filesystem_path(paths["full"]), index=False)
    actuator_frame.to_csv(filesystem_path(paths["actuator"]), index=False)
    write_json(paths["post_gate"], post_fim_gate)
    post_fim_cases.to_csv(filesystem_path(paths["post_cases"]), index=False)
    post_fd.to_csv(filesystem_path(paths["post_fd"]), index=False)
    post_grid.to_csv(filesystem_path(paths["post_grid"]), index=False)
    pd.DataFrame(local_rows).to_csv(filesystem_path(paths["local"]), index=False)
    write_json(paths["rejected"], rejected_rows)
    pd.DataFrame(history_rows).to_csv(filesystem_path(paths["history"]), index=False)
    write_json(paths["gate"], gate)
    write_json(paths["config"], config)
    write_json(paths["partition"], partition_provenance)
    write_json(
        paths["runtime"],
        {
            "search_runtime_seconds": search_runtime,
            "total_runtime_seconds": float(time.perf_counter() - total_started),
            "objective_cache_entries": len(objective_cache),
            "design_cache_entries": len(design_cache),
            "completed_seeds": completed_seed_count,
        },
    )
    outputs = list(paths.values()) + sorted(checkpoint_dir.glob("*.json.gz"))
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_final_multifidelity_search",
        config=config,
        sources={
            "adapter_manifest": adapter_run / "run_manifest.json",
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "aroma_calibration_manifest": aroma_run / "run_manifest.json",
            "temperature_actuator_manifest": actuator_run / "run_manifest.json",
            "wave1_config": CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "final_search_logic.py",
            ADAPTIVE_DIR / "hybrid_optimizer.py",
            ADAPTIVE_DIR / "local_refinement.py",
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[int(seed) for seed in config["independent_seeds"]],
        status="completed" if verdict != "FAIL" else "validation_failed",
        convergence=convergence,
        gate=gate,
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    write_json(run_dir / "run_manifest.json", manifest)
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "final_wave1_computational_qualification",
            "current_gate": gate["gate"],
            "gate_verdict": verdict,
            "latest_wave1_search_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "gate_blockers": failed + ["explicit_owner_physical_release_approval"],
            "physical_execution_status": "pending_owner_approval",
            "profiles_for_physical_execution": False,
            "physical_execution_authorized": False,
            "executable_schedule_issued": False,
            "tank_assignments": [],
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2), flush=True)
    if verdict == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
