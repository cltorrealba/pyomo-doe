from __future__ import annotations

import math
from itertools import product
from typing import Any

import numpy as np

from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (
    DesignPolicy,
    effective_temperature_change_indices,
    nutrition_policy_metrics,
    temperature_profile_metrics,
)


def approved_actuator_scenarios(config: dict[str, Any]) -> list[dict[str, float | str]]:
    """Enumerate the complete owner-approved actuator envelope (9 x 3 x 3 x 3 x 1)."""

    actuator = config["future_process"]["temperature_actuator"]
    tracking = float(np.median(np.asarray(actuator["observed_tracking_rmse_c"], dtype=float)))
    tracking_values = (-tracking, 0.0, tracking)
    scenarios: list[dict[str, float | str]] = []
    for index, (tau, tracking_error, initial_offset, probe_bias, delay) in enumerate(
        product(
            actuator["empirical_tau_h"],
            tracking_values,
            actuator["initial_temperature_uncertainty_c"],
            actuator["probe_bias_scenarios_c"],
            actuator["command_delay_scenarios_h"],
        ),
        start=1,
    ):
        scenarios.append(
            {
                "scenario": f"approved_actuator_{index:03d}",
                "tau_h": float(tau),
                "tracking_error_c": float(tracking_error),
                "initial_temperature_offset_c": float(initial_offset),
                "probe_bias_c": float(probe_bias),
                "command_delay_h": float(delay),
            }
        )
    return scenarios


def candidate_operational_metrics(
    policies: tuple[DesignPolicy, ...], config: dict[str, Any], *, margin_h: float
) -> dict[str, float | int | bool]:
    thermal = [temperature_profile_metrics(policy, config) for policy in policies]
    nutrition = [nutrition_policy_metrics(policy, config) for policy in policies]
    return {
        "maximum_temperature_jump_c": max(
            float(row["maximum_temperature_jump_c"]) for row in thermal
        ),
        "total_thermal_variation_c": sum(
            float(row["total_thermal_variation_c"]) for row in thermal
        ),
        "temperature_changes": sum(int(row["temperature_changes"]) for row in thermal),
        "total_yan_mg_l": sum(float(row["total_yan_mg_l"]) for row in nutrition),
        "minimum_action_margin_h": float(margin_h),
        "temperature_feasible": all(
            float(row["maximum_temperature_jump_c"])
            <= float(config["future_process"]["maximum_temperature_jump_c"]) + 1e-9
            and int(row["subthreshold_temperature_changes"]) == 0
            for row in thermal
        ),
        "nutrition_feasible": all(bool(row["feasible"]) for row in nutrition),
    }


def select_final_candidate(
    eligible: list[dict[str, Any]], config: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Select over every eligible 64-member candidate with the approved tie-breakers."""

    feasible = [row for row in eligible if bool(row.get("feasible", False))]
    if not feasible:
        raise ValueError("No feasible full-ensemble candidate is eligible")
    if any(not np.isfinite(float(row["full_objective"])) for row in feasible):
        raise ValueError("Eligible full-ensemble objectives must be finite")
    best_objective = min(float(row["full_objective"]) for row in feasible)
    scale = max(abs(best_objective), 1e-12)
    tolerance_fraction = float(config["search"]["practical_score_tolerance_fraction"])
    near_best = [
        row
        for row in feasible
        if float(row["full_objective"]) - best_objective
        <= tolerance_fraction * scale + 1e-12
    ]
    selected = min(
        near_best,
        key=lambda row: (
            float(row["maximum_temperature_jump_c"]),
            float(row["total_thermal_variation_c"]),
            int(row["temperature_changes"]),
            float(row["total_yan_mg_l"]),
            -float(row["minimum_action_margin_h"]),
            float(row["full_objective"]),
            str(row["candidate_id"]),
        ),
    )
    reason = {
        "primary_rule": "all eligible candidates compared on real 64-member objective",
        "near_best_tolerance_fraction": tolerance_fraction,
        "near_best_candidate_ids": [str(row["candidate_id"]) for row in near_best],
        "secondary_rule": (
            "minimum maximum jump, then minimum total variation, temperature changes, "
            "total YAN, and maximum drying margin"
        ),
        "selected_candidate_id": str(selected["candidate_id"]),
        "selected_is_exact_objective_best": math.isclose(
            float(selected["full_objective"]), best_objective, rel_tol=0.0, abs_tol=1e-12
        ),
    }
    return selected, reason


def selected_candidate_local_qualified(
    selected: dict[str, Any], eligible: list[dict[str, Any]], improvement_tolerance: float
) -> tuple[bool, str]:
    """Qualify the selected policy independently of multistart consistency."""

    if selected.get("source_type") == "accepted_refinement":
        improved = float(selected.get("full_improvement", 0.0)) > float(
            improvement_tolerance
        )
        return improved, "selected_refinement_improved_real_full_ensemble_objective"
    descendants = [
        row
        for row in eligible
        if row.get("parent_candidate_id") == selected.get("candidate_id")
        and row.get("source_type") == "accepted_refinement"
    ]
    if not descendants:
        return True, "no_full_ensemble_acceptable_local_step_found"
    best_descendant = min(float(row["full_objective"]) for row in descendants)
    retained = float(selected["full_objective"]) <= best_descendant + float(
        improvement_tolerance
    )
    return retained, "better_original_correctly_retained_after_full_ensemble_comparison"


def _pulse_vectors(policy: DesignPolicy, maximum_pulses: int) -> tuple[np.ndarray, np.ndarray]:
    times = np.full(maximum_pulses, 96.0, dtype=float)
    amounts = np.zeros(maximum_pulses, dtype=float)
    for index, (time_h, amount) in enumerate(policy.nutrition_mg_yan_l[:maximum_pulses]):
        times[index] = float(time_h)
        amounts[index] = float(amount)
    return times, amounts


def policy_distance(
    left: tuple[DesignPolicy, ...], right: tuple[DesignPolicy, ...], config: dict[str, Any]
) -> dict[str, float]:
    """Distance across setpoints, changepoints, dose, pulse timing, and complexity."""

    if len(left) != len(right):
        raise ValueError("Policy families must have the same logical profiles")
    maximum_pulses = int(config["nutrition"]["maximum_pulses"])
    setpoint_terms: list[float] = []
    change_terms: list[float] = []
    dose_terms: list[float] = []
    pulse_time_terms: list[float] = []
    complexity_terms: list[float] = []
    slots = int(config["future_process"]["optimized_temperature_slots"])
    for left_policy, right_policy in zip(left, right):
        left_temperature = np.asarray(left_policy.temperature_c, dtype=float)
        right_temperature = np.asarray(right_policy.temperature_c, dtype=float)
        setpoint_terms.append(float(np.sqrt(np.mean((left_temperature - right_temperature) ** 2)) / 12.0))
        left_changes = np.zeros(slots, dtype=float)
        right_changes = np.zeros(slots, dtype=float)
        left_changes[effective_temperature_change_indices(left_temperature, config)] = 1.0
        right_changes[effective_temperature_change_indices(right_temperature, config)] = 1.0
        change_terms.append(float(np.mean(np.abs(left_changes - right_changes))))
        left_times, left_amounts = _pulse_vectors(left_policy, maximum_pulses)
        right_times, right_amounts = _pulse_vectors(right_policy, maximum_pulses)
        dose_terms.append(float(abs(left_amounts.sum() - right_amounts.sum()) / 80.0))
        active = (left_amounts > 0.0) | (right_amounts > 0.0)
        pulse_time_terms.append(
            float(np.mean(np.abs(left_times[active] - right_times[active])) / 96.0)
            if np.any(active)
            else 0.0
        )
        left_metrics = temperature_profile_metrics(left_policy, config)
        right_metrics = temperature_profile_metrics(right_policy, config)
        complexity_terms.append(
            abs(
                int(left_metrics["temperature_changes"])
                - int(right_metrics["temperature_changes"])
            )
            / max(float(config["future_process"]["maximum_temperature_changes"]), 1.0)
        )
    components = {
        "setpoint_distance": float(np.mean(setpoint_terms)),
        "changepoint_distance": float(np.mean(change_terms)),
        "dose_distance": float(np.mean(dose_terms)),
        "pulse_time_distance": float(np.mean(pulse_time_terms)),
        "complexity_distance": float(np.mean(complexity_terms)),
    }
    components["policy_distance"] = float(
        0.40 * components["setpoint_distance"]
        + 0.20 * components["changepoint_distance"]
        + 0.15 * components["dose_distance"]
        + 0.15 * components["pulse_time_distance"]
        + 0.10 * components["complexity_distance"]
    )
    return components


def practical_convergence(
    champions: list[dict[str, Any]],
    pairwise_distances: list[dict[str, Any]],
    continuation_improvement_fraction: float,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate owner-approved score and operational-family convergence."""

    if not champions:
        return {"passed": False, "checks": {"champions_available": False}}
    best_score = max(float(row["full_ensemble_robust_score"]) for row in champions)
    tolerance = float(config["search"]["practical_score_tolerance_fraction"])
    near = [
        row
        for row in champions
        if best_score - float(row["full_ensemble_robust_score"])
        <= tolerance * max(abs(best_score), 1e-12) + 1e-12
    ]
    near_seeds = {int(row["independent_seed"]) for row in near}
    threshold = float(config["search"].get("policy_family_distance_threshold", 0.35))
    adjacency = {seed: {seed} for seed in near_seeds}
    for row in pairwise_distances:
        left, right = int(row["left_seed"]), int(row["right_seed"])
        if left in near_seeds and right in near_seeds and float(row["policy_distance"]) <= threshold:
            adjacency[left].add(right)
            adjacency[right].add(left)
    families: list[list[int]] = []
    unseen = set(near_seeds)
    while unseen:
        root = unseen.pop()
        family = {root}
        frontier = [root]
        while frontier:
            current = frontier.pop()
            for neighbor in adjacency[current] - family:
                family.add(neighbor)
                unseen.discard(neighbor)
                frontier.append(neighbor)
        families.append(sorted(family))
    required = int(config["search"]["minimum_practically_converged_seed_champions"])
    checks = {
        "champions_available": len(champions) >= len(config["independent_seeds"]),
        "at_least_three_seed_champions_within_0_5pct": len(near_seeds) >= required,
        "all_near_best_champions_feasible": all(bool(row["feasible"]) for row in near),
        "continuation_improvement_at_most_0_1pct": continuation_improvement_fraction
        <= float(config["search"]["continuation_improvement_tolerance_fraction"]) + 1e-12,
        "operationally_equivalent_policy_family": max(
            (len(family) for family in families), default=0
        )
        >= required,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "best_full_ensemble_robust_score": best_score,
        "near_best_seeds": sorted(near_seeds),
        "score_tolerance_fraction": tolerance,
        "policy_family_distance_threshold": threshold,
        "policy_families": families,
        "continuation_improvement_fraction": float(continuation_improvement_fraction),
    }
