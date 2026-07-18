from __future__ import annotations

import argparse
import json
import math
import sys
import time
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
from pilot_2026.adaptive_design.local_refinement import (  # noqa: E402
    sequential_ipopt_refine,
)
from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    anchor_policy,
    continuous_design_indices,
    decode_policy_vector,
    evaluate_campaign,
    effective_temperature_change_indices,
    load_json,
    load_wave1_config,
    policy_to_canonical_vector,
    prior_precision,
    representative_members,
    robust_information_metrics,
    temperature_profile_metrics,
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


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def hybrid_gate_verdict(checks: dict[str, bool], critical_checks: tuple[str, ...]) -> str:
    if not all(bool(checks.get(name, False)) for name in critical_checks):
        return "FAIL"
    return "PASS" if all(bool(value) for value in checks.values()) else "PASS_CONDITIONAL"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _policy_vector(policy: DesignPolicy, config: dict) -> np.ndarray:
    return policy_to_canonical_vector(policy, config)


def _seed_pairs(config: dict) -> np.ndarray:
    anchor = anchor_policy(config)
    slots = int(config["future_process"]["optimized_temperature_slots"])
    half = slots // 2
    policies = [
        anchor,
        DesignPolicy(
            "warm_cool",
            tuple([23.0] * half + [16.0] * (slots - half)),
            ((24.0, 70.0), (72.0, 70.0)),
        ),
        DesignPolicy(
            "cool_warm",
            tuple([16.0] * half + [23.0] * (slots - half)),
            ((48.0, 80.0),),
        ),
        DesignPolicy(
            "moderate_step",
            tuple([17.0] * half + [22.0] * (slots - half)),
            ((48.0, 80.0),),
        ),
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


def _decode_pair(
    values: np.ndarray,
    config: dict,
    prefix: str,
) -> tuple[tuple[DesignPolicy, DesignPolicy], np.ndarray, list[dict]]:
    single = len(vector_bounds(config))
    left = decode_policy_vector(f"{prefix}_A", values[:single], config)
    right = decode_policy_vector(f"{prefix}_B", values[single:], config)
    canonical = np.concatenate([left.canonical_vector, right.canonical_vector])
    repairs = [
        {"policy": left.policy.name, **row} for row in left.repair_log
    ] + [{"policy": right.policy.name, **row} for row in right.repair_log]
    return (left.policy, right.policy), canonical, repairs


def _continuous_pair_indices(canonical: np.ndarray, config: dict) -> np.ndarray:
    single = len(vector_bounds(config))
    left = continuous_design_indices(canonical[:single], config)
    right = continuous_design_indices(canonical[single:], config) + single
    return np.concatenate([left, right])


def _rows_for_policy(policy: DesignPolicy, config: dict) -> list[dict]:
    rows = [
        {
            "policy": policy.name,
            "action": "temperature_setpoint",
            "time_h": 0.0,
            "value": float(policy.temperature_c[0]),
            "unit": "degC",
            "canonical_effective_action": True,
        }
    ]
    slot_h = float(config["future_process"]["temperature_slot_h"])
    for slot in effective_temperature_change_indices(policy.temperature_c, config):
        rows.append(
            {
                "policy": policy.name,
                "action": "temperature_setpoint",
                "time_h": float(slot * slot_h),
                "value": float(policy.temperature_c[slot]),
                "unit": "degC",
                "canonical_effective_action": True,
            }
        )
    rows.extend(
        {
            "policy": policy.name,
            "action": "nutrition",
            "time_h": float(time_h),
            "value": float(amount),
            "unit": "mgYAN/L",
            "canonical_effective_action": True,
        }
        for time_h, amount in policy.nutrition_mg_yan_l
    )
    return rows


def _evaluation_metrics(evaluations, prior: np.ndarray, config: dict, reference=None) -> dict:
    gains = np.asarray([row.information_gain for row in evaluations], dtype=float)
    posterior = [prior + row.fim for row in evaluations]
    reference_gains = (
        None
        if reference is None
        else np.asarray([row.information_gain for row in reference], dtype=float)
    )
    metrics = robust_information_metrics(gains, posterior, reference_gains, config)
    metrics.update(
        {
            "completion_probability": float(np.mean([row.completion for row in evaluations])),
            "maximum_residual_sugar_g_l": float(
                max(row.residual_sugar_g_l for row in evaluations)
            ),
        }
    )
    return metrics


def _actuator_validation(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    representative: list[int],
    prior: np.ndarray,
    config: dict,
    partitions: dict,
    design_cache: dict,
) -> pd.DataFrame:
    actuator = config["future_process"]["temperature_actuator"]
    tau_values = sorted(float(value) for value in actuator["empirical_tau_h"])
    rows = []
    for index, tau_h in enumerate(tau_values):
        members = list(range(len(ensemble)))
        score, evaluations = evaluate_campaign(
            policies,
            ensemble,
            members,
            prior,
            config,
            partitions,
            actuator_scenario={"tau_h": tau_h},
            design_cache=design_cache,
        )
        metrics = _evaluation_metrics(evaluations, prior, config)
        rows.append(
            {
                "scenario": f"empirical_tau_{index + 1}",
                "tau_h": tau_h,
                "tracking_error_c": 0.0,
                "probe_bias_c": 0.0,
                "command_delay_h": 0.0,
                "ensemble_members": len(members),
                "score": score,
                **metrics,
                "feasible": metrics["completion_probability"]
                >= float(config["completion"]["minimum_probability"]),
            }
        )
    tracking = float(np.median(actuator["observed_tracking_rmse_c"]))
    median_tau = float(np.median(tau_values))
    for error in (-tracking, tracking):
        score, evaluations = evaluate_campaign(
            policies,
            ensemble,
            list(range(len(ensemble))),
            prior,
            config,
            partitions,
            actuator_scenario={"tau_h": median_tau, "tracking_error_c": error},
            design_cache=design_cache,
        )
        metrics = _evaluation_metrics(evaluations, prior, config)
        rows.append(
            {
                "scenario": f"tracking_error_{error:+.6f}C",
                "tau_h": median_tau,
                "tracking_error_c": error,
                "probe_bias_c": 0.0,
                "command_delay_h": 0.0,
                "ensemble_members": len(ensemble),
                "score": score,
                **metrics,
                "feasible": metrics["completion_probability"]
                >= float(config["completion"]["minimum_probability"]),
            }
        )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run corrected multiseed Wave-1 hybrid MBDoE search")
    parser.add_argument("--source-adapter-run", required=True)
    parser.add_argument("--source-ensemble-run", required=True)
    parser.add_argument("--source-aroma-run", required=True)
    parser.add_argument("--source-actuator-run", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_started = time.perf_counter()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    constraints = config["design_constraints"]
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
        raise RuntimeError("Corrected Wave-1 MBDoE adapter gate is not PASS")
    ensemble = pd.read_csv(ensemble_run / "joint_parameter_ensemble.csv")
    partitions, partition_provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    all_representative = representative_members(ensemble, config)
    search_representative = all_representative[: int(config["search"]["representative_members"])]
    rerank_representative = all_representative[: int(config["search"]["frequent_rerank_members"])]
    anchor = anchor_policy(config)
    single_bounds = vector_bounds(config)
    bounds = np.vstack([single_bounds, single_bounds])
    objective_cache: dict[tuple[tuple[int, ...], tuple[float, ...]], float] = {}
    design_cache: dict[tuple, tuple[np.ndarray, float, float]] = {}

    def objective_for_members(values: np.ndarray, members: list[int]) -> float:
        pair, canonical, _ = _decode_pair(values, config, "candidate")
        key = (tuple(int(value) for value in members), tuple(np.round(canonical, 12)))
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

    search_objective = lambda values: objective_for_members(values, search_representative)
    rerank_objective = lambda values: objective_for_members(values, rerank_representative)
    full_members = list(range(len(ensemble)))
    full_objective = lambda values: objective_for_members(values, full_members)
    search = config["search"]
    seed_profiles = _seed_pairs(config)
    allowed_seed_count = max(
        1,
        min(
            len(seed_profiles),
            int(round(float(search["seed_profile_fraction"]) * int(search["particles"]))),
        ),
    )
    seed_profiles = seed_profiles[:allowed_seed_count]
    pso_results = []
    summary_rows = []
    history_rows = []
    diversity_rows = []
    baseline_rows = []
    candidate_pool: dict[str, dict] = {}
    for independent_seed in config["independent_seeds"]:
        seed_started = time.perf_counter()
        result = particle_swarm(
            search_objective,
            bounds,
            particles=int(search["particles"]),
            iterations=int(search["maximum_iterations"]),
            seed=int(independent_seed),
            inertia=float(search["inertia"]),
            cognitive=float(search["cognitive"]),
            social=float(search["social"]),
            initial_positions=seed_profiles,
            stagnation_iterations=int(search["stagnation_iterations"]),
            improvement_tolerance=float(search["improvement_tolerance"]),
            restarts=int(search["restarts"]),
        )
        pso_results.append((int(independent_seed), result))
        for particle_index, (source, value, position) in enumerate(
            zip(result.initial_sources, result.initial_values, result.initial_positions)
        ):
            _, canonical, _ = _decode_pair(position, config, "initial")
            baseline_rows.append(
                {
                    "independent_seed": int(independent_seed),
                    "particle": particle_index,
                    "source": source,
                    "search_objective": float(value),
                    "search_robust_score": -float(value),
                    "canonical_policy_hash": sha256_payload(canonical.tolist()),
                }
            )
        for iteration, (best_value, improvement, positions, values) in enumerate(
            zip(
                result.history,
                result.improvement_history,
                result.position_history,
                result.value_history,
            )
        ):
            best_position = positions[int(np.argmin(values))]
            reranked = (
                rerank_objective(best_position)
                if iteration % int(search["rerank_interval"]) == 0
                or iteration == result.iterations_completed
                else math.nan
            )
            history_rows.append(
                {
                    "independent_seed": int(independent_seed),
                    "iteration": iteration,
                    "best_search_objective": float(best_value),
                    "best_search_robust_score": -float(best_value),
                    "leader_eight_member_objective": float(reranked),
                    "leader_eight_member_robust_score": -float(reranked),
                    "improvement": float(improvement),
                    "still_improving": bool(
                        improvement > float(search["improvement_tolerance"])
                    ),
                }
            )
            canonical_hashes = {
                sha256_payload(_decode_pair(position, config, "diversity")[1].tolist())
                for position in positions
            }
            diversity_rows.append(
                {
                    "independent_seed": int(independent_seed),
                    "iteration": iteration,
                    "normalized_position_diversity": float(result.diversity_history[iteration]),
                    "unique_canonical_policy_pairs": len(canonical_hashes),
                    "canonical_policy_pair_fraction": len(canonical_hashes) / len(positions),
                }
            )
        summary_rows.append(
            {
                "independent_seed": int(independent_seed),
                "best_search_objective": float(result.fun),
                "best_search_robust_score": -float(result.fun),
                "evaluations": int(result.evaluations),
                "iterations_completed": int(result.iterations_completed),
                "restarts_completed": int(result.restarts_completed),
                "converged": bool(result.converged),
                "stop_reason": result.stop_reason,
                "last_iteration_improvement": float(result.improvement_history[-1]),
                "final_position_diversity": float(result.diversity_history[-1]),
                "runtime_seconds": float(time.perf_counter() - seed_started),
            }
        )
        candidate_count = max(int(search["top_k_candidates"]) * 2, 4)
        for rank_index in np.argsort(result.personal_best_values)[:candidate_count]:
            raw = result.personal_best_positions[int(rank_index)]
            pair, canonical, repairs = _decode_pair(raw, config, f"seed{independent_seed}")
            key = sha256_payload(canonical.tolist())
            candidate = candidate_pool.setdefault(
                key,
                {
                    "canonical_hash": key,
                    "raw": raw.copy(),
                    "canonical": canonical,
                    "pair": pair,
                    "repairs": repairs,
                    "source_seed": int(independent_seed),
                    "source_seeds": set(),
                    "search_objective": search_objective(canonical),
                },
            )
            candidate["source_seeds"].add(int(independent_seed))
    candidates = list(candidate_pool.values())
    for candidate in candidates:
        candidate["eight_objective"] = rerank_objective(candidate["canonical"])
    candidates.sort(key=lambda row: row["eight_objective"])
    top_candidates = candidates[: int(search["top_k_candidates"])]
    full_candidate_rows = []
    for rank, candidate in enumerate(top_candidates, start=1):
        policies = (anchor, *candidate["pair"])
        score, evaluations = evaluate_campaign(
            policies,
            ensemble,
            full_members,
            prior,
            config,
            partitions,
            design_cache=design_cache,
        )
        candidate["full_score"] = score
        candidate["full_objective"] = -score
        candidate["full_evaluations"] = evaluations
        metrics = _evaluation_metrics(evaluations, prior, config)
        candidate["full_metrics"] = metrics
        full_candidate_rows.append(
            {
                "rank_after_eight_member_rerank": rank,
                "canonical_policy_hash": candidate["canonical_hash"],
                "source_seed": candidate["source_seed"],
                "source_seed_contributors": ";".join(
                    str(value) for value in sorted(candidate["source_seeds"])
                ),
                "search_robust_score": -candidate["search_objective"],
                "eight_member_robust_score": -candidate["eight_objective"],
                "full_ensemble_robust_score": score,
                **metrics,
            }
        )
    top_candidates.sort(key=lambda row: row["full_objective"])
    search_rank = pd.Series([row["search_objective"] for row in top_candidates]).rank()
    full_rank = pd.Series([row["full_objective"] for row in top_candidates]).rank()
    ranking_stability = float(search_rank.corr(full_rank)) if len(top_candidates) > 1 else 1.0

    local_trace_rows: list[dict] = []
    rejected_local_rows: list[dict] = []
    local_results = []
    for candidate_index, candidate in enumerate(
        top_candidates[: int(search["top_k_local_refinement"])], start=1
    ):
        continuous = _continuous_pair_indices(candidate["canonical"], config)
        result = sequential_ipopt_refine(
            candidate["canonical"],
            bounds,
            rerank_objective,
            continuous,
            config,
            candidate_id=f"top_{candidate_index}",
        )
        refined_pair, refined_canonical, refined_repairs = _decode_pair(
            result.x, config, f"refined_{candidate_index}"
        )
        refined_full_objective = full_objective(refined_canonical)
        full_improvement = candidate["full_objective"] - refined_full_objective
        accepted_full = bool(
            result.accepted_improvement
            and full_improvement > float(search["improvement_tolerance"])
        )
        qualified = bool(accepted_full or result.state == "stationary_no_improving_step")
        local_results.append(
            {
                "candidate": candidate,
                "result": result,
                "canonical": refined_canonical,
                "pair": refined_pair,
                "repairs": refined_repairs,
                "full_objective": refined_full_objective,
                "full_improvement": full_improvement,
                "accepted_full": accepted_full,
                "qualified": qualified,
                "state": (
                    "accepted_improvement"
                    if accepted_full
                    else result.state
                    if result.state == "stationary_no_improving_step"
                    else "rejected_model_mismatch"
                ),
            }
        )
        for row in result.trace:
            local_trace_rows.append(
                {
                    **row,
                    "full_ensemble_objective_before": candidate["full_objective"],
                    "full_ensemble_objective_after": refined_full_objective,
                    "full_ensemble_actual_improvement": full_improvement,
                    "accepted_after_full_ensemble_revalidation": accepted_full,
                }
            )
        rejected_local_rows.extend(result.rejected_candidates)
    accepted_refinements = [row for row in local_results if row["accepted_full"]]
    if accepted_refinements:
        selected_local = min(accepted_refinements, key=lambda row: row["full_objective"])
        selected_pair = selected_local["pair"]
        selected_canonical = selected_local["canonical"]
        selected_raw = selected_local["candidate"]["raw"]
        selected_repairs = selected_local["candidate"]["repairs"] + selected_local["repairs"]
        candidate_source = "ipopt_refined"
        selected_objective = selected_local["full_objective"]
    else:
        best_pso = top_candidates[0]
        selected_pair = best_pso["pair"]
        selected_canonical = best_pso["canonical"]
        selected_raw = best_pso["raw"]
        selected_repairs = best_pso["repairs"]
        candidate_source = "pso"
        selected_objective = best_pso["full_objective"]
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
    anchor_score, anchor_evaluations = evaluate_campaign(
        (anchor, anchor, anchor),
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
    actuator = _actuator_validation(
        selected_policies,
        ensemble,
        full_members,
        prior,
        config,
        partitions,
        design_cache,
    )
    all_local_qualified = bool(local_results) and all(row["qualified"] for row in local_results)
    local_candidate_accepted = bool(accepted_refinements)
    pso_evidence = bool(
        any(result.converged for _, result in pso_results)
        or ranking_stability >= float(search.get("minimum_ranking_spearman", 0.7))
    )
    finalist_seed_contributors = sorted(
        {
            seed
            for candidate in top_candidates
            for seed in candidate["source_seeds"]
        }
    )
    multi_seed_finalists = len(finalist_seed_contributors) >= int(
        search.get("minimum_finalist_seed_contributors", 2)
    )
    plateau_evidence = all(
        result.converged
        or float(result.improvement_history[-1]) <= float(search["improvement_tolerance"])
        for _, result in pso_results
    )
    margin_h = float(config["operations"]["minimum_action_to_drying_margin_h"])
    margin_probability = float(
        np.mean(
            [
                row.drying_time_h - row.latest_action_time_h >= margin_h - 1e-9
                for row in selected_evaluations
            ]
        )
    )
    temperature_metrics = [
        temperature_profile_metrics(policy, config) for policy in selected_policies
    ]
    actuator_configuration_complete = all(
        config["future_process"]["temperature_actuator"][name] is not None
        for name in (
            "probe_bias_scenarios_c",
            "command_delay_scenarios_h",
            "initial_temperature_uncertainty_c",
        )
    )
    checks = {
        "adapter_gate_pass": adapter_gate["verdict"] == "PASS",
        "source_hashes_verified": len(source_verification) == 4,
        "corrected_fim_stability_pass": bool(
            adapter_gate["checks"]["corrected_fim_validation_pass"]
        ),
        "pso_three_independent_seeds_completed": len(pso_results) >= 3,
        "pso_five_independent_seeds_completed": len(pso_results) >= 5,
        "pso_objectives_finite": bool(
            np.isfinite([result.fun for _, result in pso_results]).all()
        ),
        "pso_convergence_or_ranking_stability": pso_evidence,
        "pso_finalists_from_multiple_seeds_or_policy_convergence": bool(
            multi_seed_finalists or any(len(row["source_seeds"]) >= 2 for row in top_candidates)
        ),
        "pso_stagnation_or_plateau": plateau_evidence,
        "top_k_revalidated_on_64_members": len(top_candidates)
        == int(search["top_k_candidates"]),
        "actual_objective_evaluated": bool(local_trace_rows),
        "local_candidate_accepted": local_candidate_accepted,
        "local_refinement_qualified": all_local_qualified,
        "full_ensemble_completion_probability_at_least_95pct": selected_metrics[
            "completion_probability"
        ]
        >= float(config["completion"]["minimum_probability"]),
        "actuator_scenarios_evaluated": len(actuator) >= 11,
        "all_encoded_actuator_scenarios_feasible": bool(actuator["feasible"].all()),
        "all_approved_actuator_scenarios_feasible": bool(actuator["feasible"].all()),
        "canonical_policies_have_no_duplicate_pulses": all(
            len({time_h for time_h, _ in policy.nutrition_mg_yan_l})
            == len(policy.nutrition_mg_yan_l)
            for policy in selected_pair
        ),
        "actions_respect_24h_drying_margin": margin_probability
        >= float(config["completion"]["minimum_probability"]),
        "temperature_changes_respect_minimum_1C": all(
            int(metrics["subthreshold_temperature_changes"]) == 0
            for metrics in temperature_metrics
        ),
        "profiles_for_physical_execution_false": True,
    }
    critical_checks = (
        "adapter_gate_pass",
        "source_hashes_verified",
        "corrected_fim_stability_pass",
        "pso_three_independent_seeds_completed",
        "pso_five_independent_seeds_completed",
        "pso_objectives_finite",
        "pso_convergence_or_ranking_stability",
        "pso_finalists_from_multiple_seeds_or_policy_convergence",
        "pso_stagnation_or_plateau",
        "top_k_revalidated_on_64_members",
        "actual_objective_evaluated",
        "local_refinement_qualified",
        "full_ensemble_completion_probability_at_least_95pct",
        "canonical_policies_have_no_duplicate_pulses",
        "actions_respect_24h_drying_margin",
        "temperature_changes_respect_minimum_1C",
        "all_approved_actuator_scenarios_feasible",
    )
    if actuator_configuration_complete:
        critical_checks += ("all_encoded_actuator_scenarios_feasible",)
    verdict = hybrid_gate_verdict(checks, critical_checks)
    physical_blockers = list(constraints["fail_closed_fields"])
    if config["future_process"]["temperature_actuator"]["probe_bias_scenarios_c"] is None:
        physical_blockers.append("temperature_actuator.probe_bias_scenarios_c")
    if config["future_process"]["temperature_actuator"]["command_delay_scenarios_h"] is None:
        physical_blockers.append("temperature_actuator.command_delay_scenarios_h")
    if config["future_process"]["temperature_actuator"]["initial_temperature_uncertainty_c"] is None:
        physical_blockers.append("temperature_actuator.initial_temperature_uncertainty_c")
    failed_checks = [name for name, passed in checks.items() if not passed]
    gate = {
        "gate": "phase_D_wave1_hybrid_search_requalification",
        "verdict": verdict,
        "checks": checks,
        "conditions": [f"failed_check:{name}" for name in failed_checks]
        + [f"physical_release_blocker:{name}" for name in physical_blockers],
        "candidate_source": candidate_source,
        "pso": {
            "independent_seeds": list(config["independent_seeds"]),
            "particles_per_seed": int(search["particles"]),
            "maximum_iterations": int(search["maximum_iterations"]),
            "ranking_stability_spearman": ranking_stability,
            "finalist_seed_contributors": finalist_seed_contributors,
            "stagnation_or_plateau": plateau_evidence,
            "cache_entries": len(objective_cache),
        },
        "local_refinement": {
            "candidates_refined": len(local_results),
            "candidate_states": [row["state"] for row in local_results],
            "accepted_after_full_ensemble_revalidation": local_candidate_accepted,
            "qualified": all_local_qualified,
        },
        "full_ensemble": {
            "candidate_robust_information_score": selected_score,
            "three_anchor_robust_information_score": anchor_score,
            **selected_metrics,
            "action_margin_probability": margin_probability,
            "minimum_action_to_drying_margin_h": margin_h,
        },
        "nutrition_translation_blocker": any(
            constraints["nutrition"][name] is None
            for name in ("organic_product_yan_mass_fraction", "dap_yan_mass_fraction")
        ),
        "physical_release_blockers": physical_blockers,
        "tank_assignments": [],
        "profiles_for_physical_execution": False,
        "next_action": (
            "qualify local refinement and unresolved numerical checks"
            if verdict == "PASS_CONDITIONAL"
            else "repair critical computational failure"
            if verdict == "FAIL"
            else "optimize sampling while keeping physical release locked"
        ),
    }

    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_hybrid_search_v2", config)
    paths = {
        "summary": run_dir / "pso_multiseed_summary.csv",
        "history": run_dir / "pso_history_by_seed.csv",
        "diversity": run_dir / "pso_population_diversity.csv",
        "baselines": run_dir / "pso_seed_baselines.csv",
        "top": run_dir / "top_candidates_full_ensemble.csv",
        "raw": run_dir / "raw_particle_vector.csv",
        "canonical": run_dir / "canonical_candidate_vector.csv",
        "candidate_alias": run_dir / "candidate_vector.csv",
        "actions": run_dir / "candidate_policy_actions.csv",
        "repair": run_dir / "repair_log.json",
        "local": run_dir / "local_refinement_trace.csv",
        "rejected": run_dir / "local_refinement_rejected_candidates.json",
        "full": run_dir / "full_ensemble_validation.csv",
        "actuator": run_dir / "actuator_robustness_validation.csv",
        "gate": run_dir / "hybrid_search_gate.json",
        "config": run_dir / "wave1_mbdoe_config.json",
        "partition": run_dir / "partition_surrogate_provenance.json",
        "runtime": run_dir / "runtime_summary.json",
    }
    pd.DataFrame(summary_rows).to_csv(filesystem_path(paths["summary"]), index=False)
    pd.DataFrame(history_rows).to_csv(filesystem_path(paths["history"]), index=False)
    pd.DataFrame(diversity_rows).to_csv(filesystem_path(paths["diversity"]), index=False)
    pd.DataFrame(baseline_rows).to_csv(filesystem_path(paths["baselines"]), index=False)
    pd.DataFrame(full_candidate_rows).to_csv(filesystem_path(paths["top"]), index=False)
    pd.DataFrame(
        {"variable_index": np.arange(len(selected_raw)), "value": selected_raw}
    ).to_csv(filesystem_path(paths["raw"]), index=False)
    canonical_frame = pd.DataFrame(
        {"variable_index": np.arange(len(selected_canonical)), "value": selected_canonical}
    )
    canonical_frame.to_csv(filesystem_path(paths["canonical"]), index=False)
    canonical_frame.to_csv(filesystem_path(paths["candidate_alias"]), index=False)
    pd.DataFrame(
        sum((_rows_for_policy(policy, config) for policy in selected_policies), [])
    ).to_csv(filesystem_path(paths["actions"]), index=False)
    write_json(
        paths["repair"],
        {
            "raw_particle_to_canonical_repairs": selected_repairs,
            "raw_particle_hash": sha256_payload(selected_raw.tolist()),
            "canonical_candidate_hash": sha256_payload(selected_canonical.tolist()),
            "canonical_round_trip_exact": bool(
                np.array_equal(
                    selected_canonical,
                    np.concatenate(
                        [policy_to_canonical_vector(policy, config) for policy in selected_pair]
                    ),
                )
            ),
        },
    )
    pd.DataFrame(local_trace_rows).to_csv(filesystem_path(paths["local"]), index=False)
    write_json(paths["rejected"], rejected_local_rows)
    pd.DataFrame(
        [
            {
                "ensemble_member": row.member,
                "information_gain": row.information_gain,
                "completion": row.completion,
                "residual_sugar_g_l": row.residual_sugar_g_l,
                "drying_time_h": row.drying_time_h,
                "latest_action_time_h": row.latest_action_time_h,
                "action_margin_to_drying_h": row.drying_time_h - row.latest_action_time_h,
                "three_anchor_information_gain": reference.information_gain,
                "paired_delta_vs_three_anchor": row.information_gain - reference.information_gain,
            }
            for row, reference in zip(selected_evaluations, anchor_evaluations)
        ]
    ).to_csv(filesystem_path(paths["full"]), index=False)
    actuator.to_csv(filesystem_path(paths["actuator"]), index=False)
    write_json(paths["gate"], gate)
    write_json(paths["config"], config)
    write_json(paths["partition"], partition_provenance)
    write_json(
        paths["runtime"],
        {
            "total_runtime_seconds": float(time.perf_counter() - run_started),
            "per_seed": summary_rows,
        },
    )
    outputs = list(paths.values())
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_hybrid_search_requalification",
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
            ADAPTIVE_DIR / "hybrid_optimizer.py",
            ADAPTIVE_DIR / "local_refinement.py",
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "pilot_aroma_calibration.py",
            ADAPTIVE_DIR / "pilot_calibration.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            FERMENTATION_DIR / "shared" / "run_new_must_glycerol_estimability_doe.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[int(seed) for seed in config["independent_seeds"]],
        status="completed" if verdict != "FAIL" else "validation_failed",
        convergence={
            "pso_multiseed": summary_rows,
            "ranking_stability_spearman": ranking_stability,
            "local_states": [row["state"] for row in local_results],
        },
        gate=gate,
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    manifest_path = run_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    if sha256_file(manifest_path) == "":  # pragma: no cover - defensive
        raise RuntimeError("Final manifest hash could not be computed")
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "computational_wave1_requalification",
            "current_gate": gate["gate"],
            "gate_verdict": verdict,
            "latest_wave1_search_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "gate_blockers": failed_checks + physical_blockers,
            "physical_execution_status": "not_authorized",
            "executable_schedule_issued": False,
            "tank_assignments": [],
            "profiles_for_physical_execution": False,
            "wave1": {
                "status": f"computational_candidate_{verdict}_not_physically_released",
                "candidate_source": candidate_source,
                "ipopt_executed": bool(local_trace_rows),
                "ipopt_proposal_accepted": local_candidate_accepted,
                "local_refinement_qualified": all_local_qualified,
                "full_ensemble_robust_score": selected_score,
                "tank_randomization_required": True,
            },
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if verdict == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
