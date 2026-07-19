from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from pilot_2026.adaptive_design.final_search_logic import approved_actuator_scenarios
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (
    DesignPolicy,
    actuator_temperature_trajectory,
    effective_temperature_change_indices,
    nutrition_product_masses,
    temperature_profile_metrics,
)


WATERMARK = (
    "COMPUTATIONAL COVERAGE CANDIDATE \u2014 NOT AUTHORIZED FOR PHYSICAL EXECUTION"
)
REVIEW_SUFFIX = "_REVIEW_ONLY_NOT_AUTHORIZED"
AUTHORIZATION_FLAGS: dict[str, Any] = {
    "profiles_for_physical_execution": False,
    "physical_execution_authorized": False,
    "executable_schedule_issued": False,
    "tank_assignments": [],
}


@dataclass(frozen=True)
class CoveragePolicyMetrics:
    thermal_archetype: str
    nutrition_archetype: str
    first_24h_mean_setpoint_c: float
    nutrition_time_h: float
    maximum_internal_jump_c: float
    maximum_initial_jump_c: float
    total_thermal_variation_c: float
    temperature_changes: int


def closed_authorization() -> dict[str, Any]:
    """Return a fresh fail-closed authorization payload."""

    return {
        "profiles_for_physical_execution": False,
        "physical_execution_authorized": False,
        "executable_schedule_issued": False,
        "tank_assignments": [],
    }


def first_24h_mean_setpoint(policy: DesignPolicy) -> float:
    if len(policy.temperature_c) < 2:
        raise ValueError("A coverage policy requires at least two 12 h blocks")
    return float(np.mean(np.asarray(policy.temperature_c[:2], dtype=float)))


def classify_thermal_archetype(
    policy: DesignPolicy, coverage_config: Mapping[str, Any]
) -> str:
    """Classify only unambiguous hard-constrained cold or warm policies."""

    thermal = coverage_config["thermal"]
    first = np.asarray(policy.temperature_c[:2], dtype=float)
    first_three = np.asarray(policy.temperature_c[:3], dtype=float)
    if len(first) != 2 or len(first_three) < 3:
        raise ValueError("Thermal classification requires three 12 h blocks")
    mean = float(np.mean(first))
    cold = thermal["cold"]
    warm = thermal["warm"]
    is_cold = bool(
        mean <= float(cold["first_24h_mean_maximum_c"]) + 1e-9
        and np.any(first <= float(cold["required_block_maximum_c"]) + 1e-9)
        and np.all(first <= float(cold["first_24h_block_maximum_c"]) + 1e-9)
    )
    is_warm = bool(
        mean >= float(warm["first_24h_mean_minimum_c"]) - 1e-9
        and np.any(
            first_three[1:3] >= float(warm["required_12_to_36h_block_minimum_c"]) - 1e-9
        )
        and np.all(first >= float(warm["first_24h_block_minimum_c"]) - 1e-9)
    )
    if is_cold == is_warm:
        return "unclassified"
    return "cold" if is_cold else "warm"


def classify_nutrition_archetype(
    policy: DesignPolicy, coverage_config: Mapping[str, Any]
) -> str:
    nutrition = coverage_config["nutrition"]
    if len(policy.nutrition_mg_yan_l) != 1:
        return "unclassified"
    time_h = float(policy.nutrition_mg_yan_l[0][0])
    early = nutrition["early_window_h"]
    late = nutrition["late_window_h"]
    is_early = float(early[0]) - 1e-9 <= time_h <= float(early[1]) + 1e-9
    is_late = float(late[0]) - 1e-9 <= time_h <= float(late[1]) + 1e-9
    if is_early == is_late:
        return "unclassified"
    return "early" if is_early else "late"


def robust_initial_jump_rows(
    policy: DesignPolicy, coverage_config: Mapping[str, Any]
) -> list[dict[str, Any]]:
    first = float(policy.temperature_c[0])
    maximum = float(coverage_config["thermal"]["maximum_initial_jump_c"])
    rows = []
    for initial in coverage_config["thermal"]["initial_temperature_scenarios_c"]:
        jump = abs(first - float(initial))
        rows.append(
            {
                "policy": policy.name,
                "initial_temperature_scenario_c": float(initial),
                "initial_setpoint_c": first,
                "initial_setpoint_jump_c": jump,
                "initial_jump_within_5c": bool(jump <= maximum + 1e-9),
            }
        )
    return rows


def coverage_policy_metrics(
    policy: DesignPolicy,
    model_config: Mapping[str, Any],
    coverage_config: Mapping[str, Any],
) -> CoveragePolicyMetrics:
    thermal = temperature_profile_metrics(policy, dict(model_config))
    initial_rows = robust_initial_jump_rows(policy, coverage_config)
    nutrition_time = (
        float(policy.nutrition_mg_yan_l[0][0])
        if len(policy.nutrition_mg_yan_l) == 1
        else math.nan
    )
    return CoveragePolicyMetrics(
        thermal_archetype=classify_thermal_archetype(policy, coverage_config),
        nutrition_archetype=classify_nutrition_archetype(policy, coverage_config),
        first_24h_mean_setpoint_c=first_24h_mean_setpoint(policy),
        nutrition_time_h=nutrition_time,
        maximum_internal_jump_c=float(thermal["maximum_temperature_jump_c"]),
        maximum_initial_jump_c=max(
            float(row["initial_setpoint_jump_c"]) for row in initial_rows
        ),
        total_thermal_variation_c=float(thermal["total_thermal_variation_c"]),
        temperature_changes=int(thermal["temperature_changes"]),
    )


def archetype_checks(
    policy: DesignPolicy,
    expected_thermal: str,
    expected_nutrition: str,
    model_config: Mapping[str, Any],
    coverage_config: Mapping[str, Any],
) -> dict[str, bool]:
    metrics = coverage_policy_metrics(policy, model_config, coverage_config)
    initial_rows = robust_initial_jump_rows(policy, coverage_config)
    return {
        "thermal_archetype_exact": metrics.thermal_archetype == expected_thermal,
        "nutrition_archetype_exact": metrics.nutrition_archetype == expected_nutrition,
        "initial_jump_robust": all(
            bool(row["initial_jump_within_5c"]) for row in initial_rows
        ),
        "internal_jumps_at_most_5c": metrics.maximum_internal_jump_c
        <= float(coverage_config["thermal"]["maximum_change_c"]) + 1e-9,
        "active_internal_changes_at_least_1c": _active_changes_at_least_minimum(
            policy, coverage_config
        ),
    }


def _active_changes_at_least_minimum(
    policy: DesignPolicy, coverage_config: Mapping[str, Any]
) -> bool:
    differences = np.abs(np.diff(np.asarray(policy.temperature_c, dtype=float)))
    active = differences[differences > 1e-12]
    minimum = float(coverage_config["thermal"]["minimum_effective_change_c"])
    return bool(len(active) == 0 or np.all(active >= minimum - 1e-9))


def policies_operationally_equivalent(
    left: DesignPolicy, right: DesignPolicy, *, atol: float = 1e-9
) -> bool:
    if len(left.temperature_c) != len(right.temperature_c):
        return False
    temperature_equal = np.allclose(
        np.asarray(left.temperature_c),
        np.asarray(right.temperature_c),
        atol=atol,
        rtol=0.0,
    )
    if len(left.nutrition_mg_yan_l) != len(right.nutrition_mg_yan_l):
        return False
    nutrition_equal = all(
        math.isclose(float(lt), float(rt), abs_tol=atol, rel_tol=0.0)
        and math.isclose(float(la), float(ra), abs_tol=atol, rel_tol=0.0)
        for (lt, la), (rt, ra) in zip(left.nutrition_mg_yan_l, right.nutrition_mg_yan_l)
    )
    return bool(temperature_equal and nutrition_equal)


def contrast_checks(
    cold_policy: DesignPolicy,
    warm_policy: DesignPolicy,
    early_policy: DesignPolicy,
    late_policy: DesignPolicy,
    coverage_config: Mapping[str, Any],
) -> dict[str, bool | float]:
    thermal_difference = first_24h_mean_setpoint(warm_policy) - first_24h_mean_setpoint(
        cold_policy
    )
    early_time = float(early_policy.nutrition_mg_yan_l[0][0])
    late_time = float(late_policy.nutrition_mg_yan_l[0][0])
    nutrition_difference = late_time - early_time
    return {
        "first_24h_mean_contrast_c": thermal_difference,
        "nutrition_time_contrast_h": nutrition_difference,
        "thermal_contrast_at_least_4c": thermal_difference
        >= float(
            coverage_config["thermal"]["minimum_cold_warm_first_24h_mean_contrast_c"]
        )
        - 1e-9,
        "nutrition_contrast_at_least_36h": nutrition_difference
        >= float(coverage_config["nutrition"]["minimum_early_late_contrast_h"]) - 1e-9,
        "coverage_treatments_not_equivalent": not policies_operationally_equivalent(
            cold_policy, warm_policy
        ),
    }


def derive_robust_setpoint_envelope(
    model_config: Mapping[str, Any], coverage_config: Mapping[str, Any]
) -> dict[str, Any]:
    """Derive physical and first-block setpoint bounds from every approved scenario."""

    scenarios = approved_actuator_scenarios(dict(model_config))
    physical_min = float(coverage_config["thermal"]["physical_minimum_c"])
    physical_max = float(coverage_config["thermal"]["physical_maximum_c"])
    lower_candidates = [
        physical_min + float(row["probe_bias_c"]) - float(row["tracking_error_c"])
        for row in scenarios
    ]
    upper_candidates = [
        physical_max + float(row["probe_bias_c"]) - float(row["tracking_error_c"])
        for row in scenarios
    ]
    physical_lower = max(
        float(coverage_config["thermal"]["setpoint_minimum_c"]),
        max(lower_candidates),
    )
    physical_upper = min(
        float(coverage_config["thermal"]["setpoint_maximum_c"]),
        min(upper_candidates),
    )
    initials = np.asarray(
        coverage_config["thermal"]["initial_temperature_scenarios_c"], dtype=float
    )
    maximum_jump = float(coverage_config["thermal"]["maximum_initial_jump_c"])
    initial_lower = float(np.max(initials - maximum_jump))
    initial_upper = float(np.min(initials + maximum_jump))
    first_lower = max(physical_lower, initial_lower)
    first_upper = min(physical_upper, initial_upper)
    if physical_lower > physical_upper or first_lower > first_upper:
        raise ValueError("Approved actuator scenarios produce an empty safe envelope")
    tracking = sorted({float(row["tracking_error_c"]) for row in scenarios})
    return {
        "derivation": (
            "physical_equilibrium=setpoint-probe_bias+tracking_error; bounds are "
            "the intersection over all 243 approved actuator scenarios"
        ),
        "approved_scenario_count": len(scenarios),
        "steady_state_safe_setpoint_minimum_c": physical_lower,
        "steady_state_safe_setpoint_maximum_c": physical_upper,
        "initial_jump_only_setpoint_minimum_c": initial_lower,
        "initial_jump_only_setpoint_maximum_c": initial_upper,
        "first_block_safe_setpoint_minimum_c": first_lower,
        "first_block_safe_setpoint_maximum_c": first_upper,
        "tracking_error_scenarios_c": tracking,
        "probe_bias_scenarios_c": sorted(
            {float(row["probe_bias_c"]) for row in scenarios}
        ),
        "initial_temperature_scenarios_c": [float(value) for value in initials],
        "physical_temperature_bounds_c": [physical_min, physical_max],
        "setpoint_hardware_bounds_c": [
            float(coverage_config["thermal"]["setpoint_minimum_c"]),
            float(coverage_config["thermal"]["setpoint_maximum_c"]),
        ],
        **closed_authorization(),
        "watermark": WATERMARK,
    }


def explicit_controller_blocks(
    policy: DesignPolicy, coverage_config: Mapping[str, Any]
) -> pd.DataFrame:
    thermal = coverage_config["thermal"]
    block_h = float(thermal["block_duration_h"])
    horizon = float(thermal["controller_horizon_h"])
    total = int(round(horizon / block_h))
    if not math.isclose(total * block_h, horizon, abs_tol=1e-9):
        raise ValueError("Controller horizon must be divisible by the block duration")
    optimized = len(policy.temperature_c)
    rows = []
    for block in range(total):
        value = float(policy.temperature_c[min(block, optimized - 1)])
        rows.append(
            {
                "policy": policy.name,
                "block_index": block,
                "start_h": block * block_h,
                "end_h": (block + 1) * block_h,
                "setpoint_c": value,
                "phase": (
                    "optimized_window" if block < optimized else "explicit_hold_last"
                ),
                "hold_last_setpoint": bool(block >= optimized),
                "candidate_only_not_for_physical_execution": True,
                "watermark": WATERMARK,
            }
        )
    return pd.DataFrame(rows)


def validate_physical_temperature(
    policy: DesignPolicy,
    model_config: Mapping[str, Any],
    coverage_config: Mapping[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Evaluate all 243 actuator scenarios and audit the real initial transition."""

    lower = float(coverage_config["thermal"]["physical_minimum_c"])
    upper = float(coverage_config["thermal"]["physical_maximum_c"])
    maximum_initial = float(coverage_config["thermal"]["maximum_initial_jump_c"])
    rows = []
    global_min = math.inf
    global_max = -math.inf
    for scenario in approved_actuator_scenarios(dict(model_config)):
        trajectory = actuator_temperature_trajectory(
            policy, dict(model_config), scenario
        )
        physical = np.asarray(trajectory["physical_temperature_c"], dtype=float)
        initial_temperature = float(
            model_config["future_process"]["temperature_actuator"][
                "initial_temperature_c"
            ]
        ) + float(scenario["initial_temperature_offset_c"])
        initial_jump = abs(float(policy.temperature_c[0]) - initial_temperature)
        minimum = float(np.min(physical))
        maximum = float(np.max(physical))
        global_min = min(global_min, minimum)
        global_max = max(global_max, maximum)
        rows.append(
            {
                "candidate": policy.name,
                **scenario,
                "initial_temperature_scenario_c": initial_temperature,
                "initial_setpoint_c": float(policy.temperature_c[0]),
                "initial_setpoint_jump_c": initial_jump,
                "initial_jump_within_5c": bool(initial_jump <= maximum_initial + 1e-9),
                "minimum_physical_temperature_c": minimum,
                "maximum_physical_temperature_c": maximum,
                "physical_temperature_within_15_27c": bool(
                    minimum >= lower - 1e-9 and maximum <= upper + 1e-9
                ),
                "controller_horizon_h": float(
                    coverage_config["thermal"]["controller_horizon_h"]
                ),
                "candidate_only_not_for_physical_execution": True,
                "watermark": WATERMARK,
            }
        )
    frame = pd.DataFrame(rows)
    summary = {
        "candidate": policy.name,
        "approved_scenario_count": len(frame),
        "minimum_physical_temperature_c": global_min,
        "maximum_physical_temperature_c": global_max,
        "initial_jump_maximum_c": float(frame["initial_setpoint_jump_c"].max()),
        "initial_jump_robust_pass": bool(frame["initial_jump_within_5c"].all()),
        "physical_temperature_robust_pass": bool(
            frame["physical_temperature_within_15_27c"].all()
        ),
    }
    return frame, summary


def _nearest_free_slot(desired: int, used: set[int], minimum: int, maximum: int) -> int:
    candidates = sorted(
        range(minimum, maximum + 1), key=lambda value: (abs(value - desired), value)
    )
    return next(value for value in candidates if value not in used)


def _apply_safe_change(
    current: float, encoded: float, lower: float, upper: float
) -> float | None:
    sign = -1.0 if encoded < 0.5 else 1.0
    magnitude = 1.0 + 4.0 * abs(2.0 * encoded - 1.0)
    available = (current - lower) if sign < 0.0 else (upper - current)
    if available < 1.0 - 1e-9:
        sign *= -1.0
        available = (current - lower) if sign < 0.0 else (upper - current)
    if available < 1.0 - 1e-9:
        return None
    magnitude = min(magnitude, available, 5.0)
    return float(current + sign * magnitude)


def coverage_policy_from_unit_vector(
    name: str,
    thermal_archetype: str,
    nutrition_archetype: str,
    yan_mg_l: float,
    values: Sequence[float],
    model_config: Mapping[str, Any],
    coverage_config: Mapping[str, Any],
    robust_envelope: Mapping[str, Any],
) -> DesignPolicy:
    """Decode a unit vector directly into one hard-feasible archetype region."""

    vector = np.clip(np.asarray(values, dtype=float), 0.0, 1.0)
    if vector.shape != (12,):
        raise ValueError(
            "Coverage policy vectors must have exactly 12 unit coordinates"
        )
    slots = int(coverage_config["thermal"]["optimized_blocks"])
    lower = float(robust_envelope["steady_state_safe_setpoint_minimum_c"])
    upper = float(robust_envelope["steady_state_safe_setpoint_maximum_c"])
    first_upper = float(robust_envelope["first_block_safe_setpoint_maximum_c"])
    if thermal_archetype == "cold":
        first = lower + vector[0] * (17.0 - lower)
        maximum_delta = min(5.0, 36.0 - 2.0 * first, 20.0 - first)
        second = first + 1.0 + vector[1] * max(maximum_delta - 1.0, 0.0)
    elif thermal_archetype == "warm":
        first = 20.0 + vector[0] * (min(22.0, first_upper) - 20.0)
        second_lower = max(23.0, 44.0 - first, first + 1.0)
        second_upper = min(upper, first + 5.0)
        if second_lower > second_upper + 1e-9:
            raise ValueError("Warm archetype has no feasible second block")
        second = second_lower + vector[1] * max(second_upper - second_lower, 0.0)
    else:
        raise ValueError(f"Unknown thermal archetype: {thermal_archetype}")
    profile = np.full(slots, second, dtype=float)
    profile[0] = first
    used: set[int] = set()
    changes: list[tuple[int, float]] = []
    for index in range(3):
        offset = 2 + 3 * index
        if vector[offset] < 0.35:
            continue
        desired = 2 + int(round(vector[offset + 1] * (slots - 3)))
        slot = _nearest_free_slot(desired, used, 2, slots - 1)
        used.add(slot)
        changes.append((slot, float(vector[offset + 2])))
    current = float(second)
    for slot, encoded in sorted(changes):
        changed = _apply_safe_change(current, encoded, lower, upper)
        if changed is None:
            continue
        profile[slot:] = changed
        current = changed
    nutrition = coverage_config["nutrition"]
    if nutrition_archetype == "early":
        # Deliberately bias the search toward the preferred 2 h slot while retaining t=0.
        time_h = (
            0.0 if vector[11] >= 0.90 else float(nutrition["preferred_early_slot_h"])
        )
    elif nutrition_archetype == "late":
        allowed = list(float(value) for value in nutrition["late_slots_h"])
        chosen = min(int(math.floor(vector[11] * len(allowed))), len(allowed) - 1)
        time_h = allowed[chosen]
    else:
        raise ValueError(f"Unknown nutrition archetype: {nutrition_archetype}")
    policy = DesignPolicy(
        name,
        tuple(float(value) for value in profile),
        ((float(time_h), float(yan_mg_l)),),
    )
    expected = archetype_checks(
        policy,
        thermal_archetype,
        nutrition_archetype,
        model_config,
        coverage_config,
    )
    if not all(expected.values()):
        raise ValueError(f"Coverage decoder produced an invalid policy: {expected}")
    return policy


def nutrition_margin_row(
    strategy: str,
    dose_design: str,
    policy: DesignPolicy,
    model_config: Mapping[str, Any],
    coverage_config: Mapping[str, Any],
) -> dict[str, Any]:
    if len(policy.nutrition_mg_yan_l) != 1:
        raise ValueError("Coverage nutrition comparison requires exactly one pulse")
    yan = float(policy.nutrition_mg_yan_l[0][1])
    masses = nutrition_product_masses(yan, dict(model_config))
    limit = float(
        model_config["nutrition"]["derived_product_limits"][
            "maximum_total_organic_product_g"
        ]
    )
    organic_mass = float(masses["organic_product_g"])
    return {
        "strategy": strategy,
        "dose_design": dose_design,
        "policy": policy.name,
        "yan_total_mg_l": yan,
        "organic_product": coverage_config["nutrition"]["organic_product"],
        "organic_product_g": organic_mass,
        "dap_product_g": float(masses["dap_product_g"]),
        "organic_limit_g": limit,
        "organic_margin_g": limit - organic_mass,
        "organic_margin_fraction": (limit - organic_mass) / limit,
        "net_yan_split_50_50": bool(masses["net_yan_split_50_50"]),
        "reconstructed_yan_mg_l": float(masses["reconstructed_yan_mg_l"]),
        "nutrition_feasible": bool(
            organic_mass <= limit + 1e-9
            and float(masses["dap_product_g"])
            <= float(
                model_config["nutrition"]["derived_product_limits"][
                    "maximum_total_dap_g"
                ]
            )
            + 1e-9
        ),
        "candidate_only_not_for_physical_execution": True,
        "watermark": WATERMARK,
    }


def gate_verdict(checks: Mapping[str, Any], critical_checks: Sequence[str]) -> str:
    """Fail closed when a critical check is missing, false, or non-boolean."""

    return (
        "PASS" if all(checks.get(name) is True for name in critical_checks) else "FAIL"
    )


def lexicographic_select(
    candidates: Sequence[Mapping[str, Any]], near_best_fraction: float = 0.005
) -> tuple[dict[str, Any], dict[str, Any]]:
    feasible = [dict(row) for row in candidates if row.get("feasible") is True]
    if not feasible:
        raise ValueError("No feasible coverage candidate is available")
    best_score = max(float(row["robust_score"]) for row in feasible)
    tolerance = float(near_best_fraction) * max(abs(best_score), 1e-12)
    near = [
        row
        for row in feasible
        if best_score - float(row["robust_score"]) <= tolerance + 1e-12
    ]
    selected = min(
        near,
        key=lambda row: (
            float(row["maximum_initial_or_internal_jump_c"]),
            float(row["maximum_robust_physical_temperature_c"]),
            float(row["total_thermal_variation_c"]),
            int(row["temperature_changes"]),
            -float(row["minimum_drying_margin_h"]),
            -float(row.get("nutrition_operational_margin_fraction", 0.0)),
            -float(row["robust_score"]),
            str(row.get("candidate_id", "")),
        ),
    )
    return selected, {
        "primary_rule": "maximum robust score among physically feasible coverage candidates",
        "near_best_fraction": float(near_best_fraction),
        "best_robust_score": best_score,
        "near_best_candidate_ids": [str(row.get("candidate_id", "")) for row in near],
        "secondary_rule": (
            "minimum initial-or-internal jump, robust physical maximum, total variation, "
            "changes; maximum drying and nutrition margins"
        ),
        "selected_candidate_id": str(selected.get("candidate_id", "")),
    }


def _shared_slot_count(schedules: Mapping[str, Sequence[float]]) -> int:
    counts: dict[float, int] = {}
    for times in schedules.values():
        for value in set(float(item) for item in times):
            counts[value] = counts.get(value, 0) + 1
    return sum(1 for time_h, count in counts.items() if time_h > 0.0 and count >= 2)


def partially_harmonize_schedules(
    independent: Mapping[str, Sequence[float]],
    valid_by_policy: Mapping[str, Sequence[float]],
    score_function: Callable[[dict[str, tuple[float, ...]]], float],
    maximum_loss_fraction: float,
) -> tuple[dict[str, tuple[float, ...]], dict[str, Any]]:
    """Greedily share legal slots while retaining the independent score guardrail."""

    current = {
        name: tuple(sorted(float(value) for value in times))
        for name, times in independent.items()
    }
    base_score = float(score_function(current))
    if not np.isfinite(base_score):
        raise ValueError("Independent sampling score must be finite")
    original_union = len({time for times in current.values() for time in times})
    original_shared = _shared_slot_count(current)
    accepted_moves: list[dict[str, Any]] = []
    while True:
        targets = sorted(
            {time for times in current.values() for time in times if time > 0.0}
        )
        proposals: list[
            tuple[tuple[Any, ...], dict[str, tuple[float, ...]], dict[str, Any]]
        ] = []
        for name, times in current.items():
            valid = {float(value) for value in valid_by_policy[name]}
            for old in times:
                if old == 0.0:
                    continue
                for target in targets:
                    if target == old or target in times or target not in valid:
                        continue
                    trial = dict(current)
                    replacement = set(times)
                    replacement.remove(old)
                    replacement.add(target)
                    trial[name] = tuple(sorted(replacement))
                    score = float(score_function(trial))
                    loss_fraction = max(base_score - score, 0.0) / max(
                        abs(base_score), 1e-12
                    )
                    if loss_fraction > float(maximum_loss_fraction) + 1e-12:
                        continue
                    union = len({time for values in trial.values() for time in values})
                    shared = _shared_slot_count(trial)
                    if union >= len(
                        {time for values in current.values() for time in values}
                    ) and shared <= _shared_slot_count(current):
                        continue
                    detail = {
                        "policy": name,
                        "from_h": old,
                        "to_h": target,
                        "score": score,
                        "loss_fraction_vs_independent": loss_fraction,
                        "distinct_slots": union,
                        "additional_shared_slots": shared - original_shared,
                    }
                    rank = (union, -shared, -score, name, old, target)
                    proposals.append((rank, trial, detail))
        if not proposals:
            break
        _, current, detail = min(proposals, key=lambda row: row[0])
        accepted_moves.append(detail)
    final_score = float(score_function(current))
    final_union = len({time for times in current.values() for time in times})
    final_shared = _shared_slot_count(current)
    independent_operator_rounds = sum(len(times) for times in independent.values())
    harmonized_operator_rounds = final_union
    return current, {
        "independent_robust_score": base_score,
        "harmonized_robust_score": final_score,
        "score_loss_fraction": max(base_score - final_score, 0.0)
        / max(abs(base_score), 1e-12),
        "independent_distinct_slots": original_union,
        "harmonized_distinct_slots": final_union,
        "distinct_slot_reduction": original_union - final_union,
        "independent_shared_nonbasal_slots": original_shared,
        "harmonized_shared_nonbasal_slots": final_shared,
        "additional_shared_nonbasal_slots": final_shared - original_shared,
        "independent_operator_rounds": independent_operator_rounds,
        "harmonized_operator_rounds": harmonized_operator_rounds,
        "operator_round_reduction": independent_operator_rounds
        - harmonized_operator_rounds,
        "accepted_moves": accepted_moves,
    }


def retained_information_improvement_fraction(
    coverage_score: float, information_score: float, three_anchor_score: float
) -> float:
    denominator = float(information_score) - float(three_anchor_score)
    if denominator <= 0.0:
        return math.nan
    return (float(coverage_score) - float(three_anchor_score)) / denominator
