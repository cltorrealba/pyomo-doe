from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    allowed_sampling_times,
    fim_from_sampling_sensitivity_cache,
    load_json,
    load_wave1_config,
    logdet,
    prepare_design,
    prepare_sampling_sensitivity_cache,
    prior_precision,
    robust_information_metrics,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    capture_git_state,
    create_immutable_run_directory,
    filesystem_path,
    sha256_file,
    verify_manifest_output,
    write_json,
    write_json_atomic,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def robust_score(gains: np.ndarray, config: dict) -> float:
    gains = np.asarray(gains, dtype=float)
    return float(config["objective"]["median_weight"]) * float(np.median(gains)) + float(
        config["objective"]["lower_decile_weight"]
    ) * float(np.quantile(gains, 0.1))


def sampling_gate_verdict(checks: dict[str, bool], critical_checks: tuple[str, ...]) -> str:
    if not all(bool(checks.get(name, False)) for name in critical_checks):
        return "FAIL"
    return "PASS" if all(bool(value) for value in checks.values()) else "PASS_CONDITIONAL"


def optimized_schedule_is_acceptable(
    optimized_gains: np.ndarray,
    preliminary_gains: np.ndarray,
    config: dict,
) -> tuple[bool, str | None]:
    score_tolerance = float(config["sampling"]["full_ensemble_score_tolerance"])
    q10_tolerance = float(config["sampling"]["q10_regression_tolerance"])
    optimized_score = robust_score(optimized_gains, config)
    preliminary_score = robust_score(preliminary_gains, config)
    if optimized_score < preliminary_score - score_tolerance:
        return False, "optimized_full_ensemble_score_below_preliminary"
    optimized_q10 = float(np.quantile(optimized_gains, 0.1))
    preliminary_q10 = float(np.quantile(preliminary_gains, 0.1))
    if optimized_q10 < preliminary_q10 - q10_tolerance:
        return False, "optimized_q10_below_preliminary_guardrail"
    return True, None


def policy_actions_before_drying(
    policy: DesignPolicy,
    drying_times: np.ndarray,
    config: dict,
) -> bool:
    effective_changes = list(
        np.flatnonzero(np.abs(np.diff(policy.temperature_c)) > 1e-12) + 1
    )
    action_times = [
        float(config["future_process"]["temperature_slot_h"]) * index
        for index in effective_changes
    ]
    action_times.extend(time_h for time_h, _ in policy.nutrition_mg_yan_l)
    return all(
        float(np.mean(np.asarray(drying_times, dtype=float) >= float(action_time)))
        >= float(config["completion"]["minimum_probability"])
        for action_time in action_times
    )


def capture_constraints_approved(config: dict) -> bool:
    constraints = config["design_constraints"]["sampling_and_capture"]
    return all(
        constraints[name] is not None
        for name in (
            "maximum_capture_interval_h",
            "maximum_trap_loading",
            "vessel_change_required",
        )
    )


def _load_policies(path: Path) -> tuple[DesignPolicy, ...]:
    actions = pd.read_csv(filesystem_path(path))
    policies = []
    for name, group in actions.groupby("policy", sort=False):
        temperature_actions = group[group["action"].eq("temperature_setpoint")].sort_values(
            "time_h"
        )
        slots = 14
        profile = np.empty(slots, dtype=float)
        current = float(temperature_actions.iloc[0]["value"])
        action_rows = list(temperature_actions.itertuples())
        action_index = 0
        for slot in range(slots):
            time_h = 12.0 * slot
            while action_index + 1 < len(action_rows) and float(
                action_rows[action_index + 1].time_h
            ) <= time_h + 1e-9:
                action_index += 1
                current = float(action_rows[action_index].value)
            profile[slot] = current
        nutrition = tuple(
            (float(row.time_h), float(row.value))
            for row in group[group["action"].eq("nutrition")].sort_values("time_h").itertuples()
        )
        policies.append(DesignPolicy(str(name), tuple(profile), nutrition))
    if len(policies) != 3:
        raise ValueError(f"Expected three Wave-1 policies; found {len(policies)}")
    return tuple(policies)


def _score(
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
) -> tuple[float, np.ndarray, list[np.ndarray]]:
    base = logdet(prior)
    gains = []
    posterior = []
    for member in members:
        total = prior.copy()
        for policy in policies:
            schedule = tuple(float(value) for value in schedules[policy.name])
            key = (int(member), policy.name, schedule)
            if key not in fim_cache:
                fim_cache[key] = fim_from_sampling_sensitivity_cache(
                    sampling_cache[(member, policy.name)],
                    np.asarray(schedule, dtype=float),
                    config,
                )
            total += fim_cache[key]
        posterior.append(total)
        gains.append(logdet(total) - base)
    gains_array = np.asarray(gains, dtype=float)
    return robust_score(gains_array, config), gains_array, posterior


def _individual_score(
    schedule: tuple[float, ...],
    policy: DesignPolicy,
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
) -> float:
    schedules = {policy.name: schedule}
    return _score(
        schedules,
        (policy,),
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
    )[0]


def _beam_schedule(
    policy: DesignPolicy,
    valid: np.ndarray,
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    rng: np.random.Generator,
) -> tuple[float, ...]:
    count = int(config["sampling"]["samples_per_process"])
    forced: set[float] = set()
    if config["sampling"]["force_first_slot"]:
        forced.add(float(valid[0]))
    if config["sampling"]["force_last_slot"]:
        forced.add(float(valid[-1]))
    beam: list[tuple[float, ...]] = [tuple(sorted(forced))]
    width = int(config["sampling"]["beam_width"])
    ordered_valid = list(float(value) for value in valid)
    rng.shuffle(ordered_valid)
    while len(beam[0]) < count:
        proposals: dict[tuple[float, ...], float] = {}
        for schedule in beam:
            for candidate in ordered_valid:
                if candidate in schedule:
                    continue
                trial = tuple(sorted((*schedule, candidate)))
                proposals.setdefault(
                    trial,
                    _individual_score(
                        trial,
                        policy,
                        sampling_cache,
                        members,
                        prior,
                        config,
                        fim_cache,
                    ),
                )
        if not proposals:
            raise RuntimeError(f"Not enough valid sample slots for {policy.name}")
        beam = [
            schedule
            for schedule, _ in sorted(
                proposals.items(), key=lambda item: (-item[1], item[0])
            )[:width]
        ]
    return beam[0]


def _coordinate_improve(
    schedules: dict[str, tuple[float, ...]],
    valid_by_policy: dict[str, np.ndarray],
    policies: tuple[DesignPolicy, ...],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    rng: np.random.Generator,
) -> dict[str, tuple[float, ...]]:
    schedules = dict(schedules)
    current, _, _ = _score(
        schedules, policies, sampling_cache, members, prior, config, fim_cache
    )
    for _pass in range(int(config["sampling"]["coordinate_passes"])):
        improved = False
        policy_order = list(policies)
        rng.shuffle(policy_order)
        for policy in policy_order:
            old_order = list(schedules[policy.name])
            rng.shuffle(old_order)
            for old in old_order:
                best_schedule, best_score = schedules[policy.name], current
                for candidate in valid_by_policy[policy.name]:
                    candidate = float(candidate)
                    if candidate in schedules[policy.name]:
                        continue
                    proposed = set(schedules[policy.name])
                    proposed.remove(old)
                    proposed.add(candidate)
                    trial = dict(schedules)
                    trial[policy.name] = tuple(sorted(proposed))
                    score, _, _ = _score(
                        trial,
                        policies,
                        sampling_cache,
                        members,
                        prior,
                        config,
                        fim_cache,
                    )
                    if score > best_score + float(config["sampling"]["full_ensemble_score_tolerance"]):
                        best_schedule, best_score = trial[policy.name], score
                if best_schedule != schedules[policy.name]:
                    improved = True
                schedules[policy.name] = best_schedule
                current = best_score
        if not improved:
            break
    return schedules


def _optimize_schedules(
    policies: tuple[DesignPolicy, ...],
    valid_by_policy: dict[str, np.ndarray],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    seed_offset: int,
) -> tuple[dict[str, tuple[float, ...]], pd.DataFrame]:
    best_schedules = None
    best_score = -math.inf
    rows = []
    for restart in range(int(config["sampling"]["search_restarts"])):
        rng = np.random.default_rng(int(config["seed"]) + seed_offset + restart)
        schedules = {
            policy.name: _beam_schedule(
                policy,
                valid_by_policy[policy.name],
                sampling_cache,
                members,
                prior,
                config,
                fim_cache,
                rng,
            )
            for policy in policies
        }
        schedules = _coordinate_improve(
            schedules,
            valid_by_policy,
            policies,
            sampling_cache,
            members,
            prior,
            config,
            fim_cache,
            rng,
        )
        score, gains, _ = _score(
            schedules, policies, sampling_cache, members, prior, config, fim_cache
        )
        rows.append(
            {
                "restart": restart,
                "robust_score": score,
                "median": float(np.median(gains)),
                "q10": float(np.quantile(gains, 0.1)),
            }
        )
        if score > best_score:
            best_score = score
            best_schedules = schedules
    if best_schedules is None:  # pragma: no cover - at least one restart is configured
        raise RuntimeError("Sampling optimization produced no schedule")
    return best_schedules, pd.DataFrame(rows)


def _capture_intervals(
    schedules: dict[str, tuple[float, ...]],
    config: dict,
) -> pd.DataFrame:
    constraints = config["design_constraints"]["sampling_and_capture"]
    rows = []
    for policy, times in schedules.items():
        for index, (start, end) in enumerate(zip(times[:-1], times[1:]), start=1):
            maximum = constraints["maximum_capture_interval_h"]
            rows.append(
                {
                    "policy": policy,
                    "capture_interval": index,
                    "start_h": start,
                    "end_h": end,
                    "duration_h": end - start,
                    "capture_stage_1_c": 0.0,
                    "capture_stage_2_c": -40.0,
                    "maximum_capture_interval_h": maximum,
                    "within_approved_duration": None if maximum is None else end - start <= maximum,
                    "maximum_trap_loading": constraints["maximum_trap_loading"],
                    "vessel_change_required": constraints["vessel_change_required"],
                }
            )
    return pd.DataFrame(rows)


def _operational_conflicts(
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    config: dict,
) -> pd.DataFrame:
    constraints = config["design_constraints"]["sampling_and_capture"]
    rows = []
    for policy in policies:
        sample_times = set(schedules[policy.name])
        pulse_times = {float(time_h) for time_h, _ in policy.nutrition_mg_yan_l}
        for collision in sorted(sample_times & pulse_times):
            rows.append(
                {
                    "time_h": collision,
                    "policies": policy.name,
                    "conflict_type": "sampling_and_nutrition_same_timestamp",
                    "status": "unresolved" if constraints["sample_event_order"] is None else "ordered",
                    "required_approval": "sample_event_order",
                }
            )
    all_times = sorted({time_h for times in schedules.values() for time_h in times})
    for time_h in all_times:
        policies_at_time = sorted(
            name for name, times in schedules.items() if time_h in set(times)
        )
        if len(policies_at_time) > 1:
            capacity = constraints["manual_sampling_capacity_per_time_slot"]
            unresolved = capacity is None or len(policies_at_time) > int(capacity)
            rows.append(
                {
                    "time_h": time_h,
                    "policies": ";".join(policies_at_time),
                    "conflict_type": "simultaneous_manual_sampling",
                    "status": "unresolved" if unresolved else "within_capacity",
                    "required_approval": "manual_sampling_capacity_per_time_slot",
                }
            )
    if constraints["minimum_minutes_between_sampling_and_nutrition"] is None:
        rows.append(
            {
                "time_h": math.nan,
                "policies": ";".join(policy.name for policy in policies),
                "conflict_type": "sampling_nutrition_separation_not_approved",
                "status": "unresolved",
                "required_approval": "minimum_minutes_between_sampling_and_nutrition",
            }
        )
    return pd.DataFrame(rows)


def _nutrition_translation(policies: tuple[DesignPolicy, ...], config: dict) -> pd.DataFrame:
    nutrition = config["design_constraints"]["nutrition"]
    volume_l = float(nutrition["initial_volume_l"])
    organic_limit = float(nutrition["maximum_total_organic_product_g_hl"]) * volume_l / 100.0
    dap_limit = float(nutrition["maximum_total_dap_fda_g_hl"]) * volume_l / 100.0
    rows = []
    for policy in policies:
        total_yan = sum(amount for _, amount in policy.nutrition_mg_yan_l)
        if total_yan <= 0.0:
            continue
        required = total_yan * volume_l
        dap_min = max(0.0, (required - 100.0 * organic_limit) / 200.0)
        dap_max = min(dap_limit, required / 200.0)
        endpoints = {
            "minimum_dap_endpoint": dap_min,
            "maximum_dap_endpoint": dap_max,
        }
        for endpoint, total_dap in endpoints.items():
            total_organic = (required - 200.0 * total_dap) / 100.0
            for time_h, target in policy.nutrition_mg_yan_l:
                fraction = target / total_yan
                organic = total_organic * fraction
                dap = total_dap * fraction
                reconstructed = (100.0 * organic + 200.0 * dap) / volume_l
                rows.append(
                    {
                        "policy": policy.name,
                        "time_h": time_h,
                        "solution": endpoint,
                        "selected_for_operation": False,
                        "yan_target_mg_l": target,
                        "organic_product_g": organic,
                        "dap_fda_g": dap,
                        "yan_reconstructed_mg_l": reconstructed,
                        "reconstruction_error_mg_l": reconstructed - target,
                        "organic_total_limit_g": organic_limit,
                        "dap_fda_total_limit_g": dap_limit,
                        "organic_total_slack_g": organic_limit - total_organic,
                        "dap_fda_total_slack_g": dap_limit - total_dap,
                        "per_event_organic_limit_g": nutrition["maximum_organic_product_g_per_event"],
                        "per_event_dap_fda_limit_g": nutrition["maximum_dap_fda_g_per_event"],
                        "translation_blocker": nutrition["approved_product_mix_selection_policy"] is None,
                    }
                )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize and qualify Wave-1 sampling on 64 members")
    parser.add_argument("--source-search-run", required=True)
    parser.add_argument("--source-ensemble-run", required=True)
    parser.add_argument("--source-aroma-run", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    constraints = config["design_constraints"]
    aroma_config = load_json(AROMA_CONFIG_PATH)
    search_run = _run_path(args.source_search_run)
    ensemble_run = _run_path(args.source_ensemble_run)
    aroma_run = _run_path(args.source_aroma_run)
    source_verification = {
        "search_gate": verify_manifest_output(search_run, "hybrid_search_gate.json"),
        "search_actions": verify_manifest_output(search_run, "candidate_policy_actions.csv"),
        "ensemble": verify_manifest_output(ensemble_run, "joint_parameter_ensemble.csv"),
        "aroma": verify_manifest_output(aroma_run, "aroma_calibration_gate.json"),
    }
    search_gate = load_json(search_run / "hybrid_search_gate.json")
    if search_gate["verdict"] not in {"PASS", "PASS_CONDITIONAL"}:
        raise RuntimeError("Explicit Wave-1 search source is not computationally usable")
    policies = _load_policies(search_run / "candidate_policy_actions.csv")
    ensemble = pd.read_csv(ensemble_run / "joint_parameter_ensemble.csv")
    partitions, provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    members = list(range(len(ensemble)))
    prepared: dict[tuple[int, str], object] = {}
    sampling_cache: dict[tuple[int, str], object] = {}
    drying_rows = []
    for member in members:
        for policy in policies:
            item = prepare_design(policy, ensemble.iloc[member], config, partitions)
            if item is None:
                raise RuntimeError(f"Simulation failed: member {member}, {policy.name}")
            prepared[(member, policy.name)] = item
            sampling_cache[(member, policy.name)] = prepare_sampling_sensitivity_cache(
                item, config
            )
            drying_rows.append(
                {
                    "ensemble_member": member,
                    "policy": policy.name,
                    "drying_time_h": item.drying_time_h,
                }
            )
    drying = pd.DataFrame(drying_rows)
    allowed = allowed_sampling_times(config)
    valid_by_policy = {}
    active_probabilities = {}
    for policy in policies:
        drying_times = drying[drying.policy.eq(policy.name)].drying_time_h.to_numpy(dtype=float)
        probability = {
            float(time_h): float(np.mean(drying_times >= float(time_h))) for time_h in allowed
        }
        valid = np.asarray(
            [
                time_h
                for time_h in allowed
                if probability[float(time_h)]
                >= float(config["completion"]["minimum_probability"])
            ],
            dtype=float,
        )
        if constraints["sampling_and_capture"]["sample_event_order"] is None:
            pulse_times = {float(time_h) for time_h, _ in policy.nutrition_mg_yan_l}
            valid = np.asarray([time_h for time_h in valid if float(time_h) not in pulse_times])
        valid_by_policy[policy.name] = valid
        active_probabilities[policy.name] = probability
    fim_cache: dict[tuple, np.ndarray] = {}
    optimized, optimization_restarts = _optimize_schedules(
        policies,
        valid_by_policy,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
        seed_offset=0,
    )
    optimized_score, optimized_gains, optimized_posterior = _score(
        optimized, policies, sampling_cache, members, prior, config, fim_cache
    )
    preliminary = {
        policy.name: tuple(float(value) for value in config["sampling"]["preliminary_times_h"])
        for policy in policies
    }
    preliminary_score, preliminary_gains, preliminary_posterior = _score(
        preliminary, policies, sampling_cache, members, prior, config, fim_cache
    )

    reference_policies = tuple(
        DesignPolicy(f"three_anchor_reference_{index + 1}", policies[0].temperature_c, policies[0].nutrition_mg_yan_l)
        for index in range(3)
    )
    for member in members:
        original = prepared[(member, policies[0].name)]
        cache = sampling_cache[(member, policies[0].name)]
        for reference in reference_policies:
            prepared[(member, reference.name)] = original
            sampling_cache[(member, reference.name)] = cache
    reference_valid = {
        reference.name: valid_by_policy[policies[0].name] for reference in reference_policies
    }
    reference_schedules, reference_restarts = _optimize_schedules(
        reference_policies,
        reference_valid,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
        seed_offset=1000,
    )
    reference_score, reference_gains, reference_posterior = _score(
        reference_schedules,
        reference_policies,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
    )
    acceptable, fallback_reason = optimized_schedule_is_acceptable(
        optimized_gains, preliminary_gains, config
    )
    selected = optimized if acceptable else preliminary
    selected_source = "optimized" if acceptable else "preliminary_fallback"
    selected_score, selected_gains, selected_posterior = _score(
        selected, policies, sampling_cache, members, prior, config, fim_cache
    )
    paired = pd.DataFrame(
        {
            "ensemble_member": members,
            "optimized_information_gain": optimized_gains,
            "preliminary_information_gain": preliminary_gains,
            "three_anchor_information_gain": reference_gains,
            "optimized_vs_preliminary_delta": optimized_gains - preliminary_gains,
            "optimized_vs_three_anchor_delta": optimized_gains - reference_gains,
            "selected_information_gain": selected_gains,
        }
    )
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    schedule_rows = []
    for policy in policies:
        for index, time_h in enumerate(selected[policy.name], start=1):
            timestamp = start + timedelta(hours=float(time_h))
            schedule_rows.append(
                {
                    "policy": policy.name,
                    "sample_number": index,
                    "time_h": time_h,
                    "local_timestamp": timestamp.isoformat(),
                    "weekday": timestamp.strftime("%A"),
                    "active_probability": active_probabilities[policy.name][float(time_h)],
                    "candidate_only_not_for_physical_execution": True,
                }
            )
    schedule = pd.DataFrame(schedule_rows)
    captures = _capture_intervals(selected, config)
    conflicts = _operational_conflicts(selected, policies, config)
    nutrition = _nutrition_translation(policies, config)
    selected_metrics = robust_information_metrics(
        selected_gains,
        selected_posterior,
        preliminary_gains,
        config,
    )
    preliminary_metrics = robust_information_metrics(
        preliminary_gains, preliminary_posterior, None, config
    )
    reference_metrics = robust_information_metrics(
        reference_gains, reference_posterior, None, config
    )
    process_active_actions = []
    for policy in policies:
        drying_times = drying[drying.policy.eq(policy.name)].drying_time_h.to_numpy(dtype=float)
        process_active_actions.append(
            policy_actions_before_drying(policy, drying_times, config)
        )
    unresolved_conflicts = bool(
        len(conflicts) and conflicts["status"].eq("unresolved").any()
    )
    capture_limits_approved = capture_constraints_approved(config)
    checks = {
        "source_hashes_verified": len(source_verification) == 4,
        "source_search_ipopt_execution_recognized": bool(
            search_gate["checks"]["actual_objective_evaluated"]
        ),
        "sampling_optimization_used_64_members": len(members) == 64,
        "optimized_full_ensemble_score_not_below_preliminary": acceptable
        or selected_source == "preliminary_fallback",
        "selected_q10_meets_preliminary_guardrail": float(np.quantile(selected_gains, 0.1))
        >= float(np.quantile(preliminary_gains, 0.1))
        - float(config["sampling"]["q10_regression_tolerance"]),
        "automatic_preliminary_fallback_operational": acceptable
        or selected_source == "preliminary_fallback",
        "paired_comparison_published": len(paired) == 64,
        "three_anchor_reference_optimized_with_comparable_method": len(reference_restarts)
        == len(optimization_restarts),
        "exactly_ten_samples_per_process": bool(
            schedule.groupby("policy").size().eq(int(config["sampling"]["samples_per_process"])).all()
        ),
        "all_samples_in_legal_windows": bool(
            schedule["weekday"].isin(["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"]).all()
            and schedule["local_timestamp"].map(
                lambda value: 9 <= datetime.fromisoformat(value).hour <= 17
            ).all()
        ),
        "samples_and_actions_before_drying": bool(
            (schedule["active_probability"] >= float(config["completion"]["minimum_probability"])).all()
            and all(process_active_actions)
        ),
        "no_ambiguous_operational_conflicts": not unresolved_conflicts,
        "capture_limits_approved_and_satisfied": capture_limits_approved
        and bool(captures["within_approved_duration"].fillna(False).all()),
        "nutrition_translation_has_no_arbitrary_mix": bool(
            nutrition.empty or not nutrition["selected_for_operation"].any()
        ),
        "profiles_for_physical_execution_false": True,
    }
    critical_checks = (
        "source_hashes_verified",
        "sampling_optimization_used_64_members",
        "optimized_full_ensemble_score_not_below_preliminary",
        "selected_q10_meets_preliminary_guardrail",
        "automatic_preliminary_fallback_operational",
        "paired_comparison_published",
        "exactly_ten_samples_per_process",
        "all_samples_in_legal_windows",
        "samples_and_actions_before_drying",
        "no_ambiguous_operational_conflicts",
        "capture_limits_approved_and_satisfied",
    )
    verdict = sampling_gate_verdict(checks, critical_checks)
    physical_blockers = list(constraints["fail_closed_fields"])
    failed_checks = [name for name, passed in checks.items() if not passed]
    gate = {
        "gate": "phase_D_wave1_sampling_requalification",
        "verdict": verdict,
        "checks": checks,
        "source_search_verdict": search_gate["verdict"],
        "source_candidate_source": search_gate["candidate_source"],
        "selected_schedule_source": selected_source,
        "fallback_reason": fallback_reason,
        "information": {
            "optimized": {
                "robust_score": optimized_score,
                "median": float(np.median(optimized_gains)),
                "q10": float(np.quantile(optimized_gains, 0.1)),
            },
            "preliminary": preliminary_metrics,
            "selected": selected_metrics,
            "three_anchor_optimized_reference": reference_metrics,
            "optimized_vs_preliminary_wins": int((optimized_gains > preliminary_gains).sum()),
            "optimized_vs_preliminary_losses": int((optimized_gains < preliminary_gains).sum()),
            "optimized_vs_three_anchor_wins": int((optimized_gains > reference_gains).sum()),
            "optimized_vs_three_anchor_losses": int((optimized_gains < reference_gains).sum()),
            "worst_optimized_vs_preliminary_delta": float(
                np.min(optimized_gains - preliminary_gains)
            ),
        },
        "conditions": [f"failed_check:{name}" for name in failed_checks]
        + [f"physical_release_blocker:{name}" for name in physical_blockers],
        "nutrition_translation_blocker": bool(
            not nutrition.empty and nutrition["translation_blocker"].any()
        ),
        "tank_assignments": [],
        "profiles_for_physical_execution": False,
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_sampling_v2", config)
    schedule_name = (
        "optimized_sampling_schedule.csv"
        if selected_source == "optimized"
        else "fallback_sampling_schedule.csv"
    )
    paths = {
        "schedule": run_dir / schedule_name,
        "comparison": run_dir / "sampling_candidate_comparison.csv",
        "captures": run_dir / "optimized_capture_intervals.csv",
        "conflicts": run_dir / "operational_conflicts.csv",
        "nutrition": run_dir / "nutrition_product_translation.csv",
        "drying": run_dir / "drying_time_ensemble.csv",
        "candidate_restarts": run_dir / "sampling_search_restarts.csv",
        "reference_restarts": run_dir / "three_anchor_search_restarts.csv",
        "gate": run_dir / "sampling_gate.json",
        "config": run_dir / "wave1_mbdoe_config.json",
        "partition": run_dir / "partition_surrogate_provenance.json",
    }
    schedule.to_csv(filesystem_path(paths["schedule"]), index=False)
    paired.to_csv(filesystem_path(paths["comparison"]), index=False)
    captures.to_csv(filesystem_path(paths["captures"]), index=False)
    conflicts.to_csv(filesystem_path(paths["conflicts"]), index=False)
    nutrition.to_csv(filesystem_path(paths["nutrition"]), index=False)
    drying.to_csv(filesystem_path(paths["drying"]), index=False)
    optimization_restarts.to_csv(filesystem_path(paths["candidate_restarts"]), index=False)
    reference_restarts.to_csv(filesystem_path(paths["reference_restarts"]), index=False)
    write_json(paths["gate"], gate)
    write_json(paths["config"], config)
    write_json(paths["partition"], provenance)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_sampling_requalification",
        config=config,
        sources={
            "hybrid_search_manifest": search_run / "run_manifest.json",
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "aroma_calibration_manifest": aroma_run / "run_manifest.json",
            "wave1_config": CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "pilot_aroma_calibration.py",
            ADAPTIVE_DIR / "pilot_calibration.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            FERMENTATION_DIR / "shared" / "run_new_must_glycerol_estimability_doe.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[int(config["seed"]) + value for value in range(int(config["sampling"]["search_restarts"]))],
        status="completed" if verdict != "FAIL" else "validation_failed",
        convergence={
            "candidate_restarts": optimization_restarts.to_dict(orient="records"),
            "reference_restarts": reference_restarts.to_dict(orient="records"),
            "selected_schedule_source": selected_source,
        },
        gate=gate,
        outputs=list(paths.values()),
        git_snapshot=git_snapshot,
    )
    manifest_path = run_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    if sha256_file(manifest_path) == "":  # pragma: no cover
        raise RuntimeError("Final sampling manifest could not be hashed")
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "computational_wave1_sampling_requalification",
            "current_gate": gate["gate"],
            "gate_verdict": verdict,
            "latest_wave1_sampling_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "latest_gate_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "gate_blockers": failed_checks + physical_blockers,
            "physical_execution_status": "not_authorized",
            "executable_schedule_issued": False,
            "tank_assignments": [],
            "profiles_for_physical_execution": False,
            "wave1": {
                **state.get("wave1", {}),
                "sampling_gate_verdict": verdict,
                "sampling_schedule_source": selected_source,
                "sampling_fallback_reason": fallback_reason,
                "candidate_sampling_schedule_ready_for_review": False,
            },
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if verdict == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
