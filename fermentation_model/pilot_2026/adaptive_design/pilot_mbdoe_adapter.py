from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pilot_2026.adaptive_design.pilot_aroma_calibration import (
    AromaForcing,
    PARAMETER_NAMES,
    _partition_basis,
    load_partition_surrogates,
    simulate_aroma,
)
from pilot_2026.adaptive_design import pilot_calibration


SPECIES = ("ethyl_acetate", "ethyl_octanoate", "isoamyl_acetate")


@dataclass(frozen=True)
class DesignPolicy:
    name: str
    temperature_c: tuple[float, ...]
    nutrition_mg_yan_l: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class CanonicalDesign:
    policy: DesignPolicy
    raw_vector: np.ndarray
    canonical_vector: np.ndarray
    repair_log: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class ScenarioEvaluation:
    member: int
    information_gain: float
    completion: bool
    residual_sugar_g_l: float
    fim: np.ndarray
    drying_time_h: float
    latest_action_time_h: float


@dataclass(frozen=True)
class PreparedDesign:
    policy: DesignPolicy
    member: int
    forcings: dict[str, AromaForcing]
    aroma_log_values: dict[str, np.ndarray]
    residual_sugar_g_l: float
    drying_time_h: float


@dataclass(frozen=True)
class FIMComponents:
    physical_prediction: np.ndarray
    nominal_sigma: np.ndarray
    scaled_observation: np.ndarray
    sensitivity: np.ndarray
    fim: np.ndarray
    difference_methods: tuple[str, ...]


@dataclass(frozen=True)
class SamplingSensitivityCache:
    prepared: PreparedDesign
    nominal_liquid: dict[str, np.ndarray]
    nominal_captured: dict[str, np.ndarray]
    liquid_derivatives: dict[str, np.ndarray]
    captured_derivatives: dict[str, np.ndarray]
    difference_methods: dict[str, tuple[str, ...]]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_wave1_config(
    config_path: Path,
    constraints_path: Path | None = None,
) -> dict[str, Any]:
    """Load Wave-1 settings and derive every physical limit from constraints."""

    config = load_json(config_path)
    if constraints_path is None:
        declared = config.get("design_constraints_path")
        constraints_path = (
            Path(config_path).parent / "design_constraints.json"
            if declared is None
            else Path(config_path).resolve().parents[3] / declared
        )
    constraints = load_json(Path(constraints_path))
    config["design_constraints"] = constraints
    temperature = constraints["temperature"]
    nutrition = constraints["nutrition"]
    completion = constraints["completion"]
    sampling = constraints["sampling_and_capture"]
    manual = constraints["manual_action_window"]
    config["future_process"]["temperature_minimum_c"] = float(temperature["minimum_c"])
    config["future_process"]["temperature_maximum_c"] = float(temperature["maximum_c"])
    config["future_process"]["temperature_slot_h"] = float(
        temperature["minimum_segment_duration_h"]
    )
    config["nutrition"].update(
        {
            "maximum_pulses": int(nutrition["maximum_pulses"]),
            "maximum_yan_per_pulse_mg_l": float(nutrition["maximum_yan_per_pulse_mg_l"]),
            "maximum_total_yan_mg_l": float(nutrition["maximum_total_yan_mg_l"]),
            "latest_h": float(nutrition["latest_allowable_pulse_h"]),
            "manual_weekdays": list(range(len(manual["weekdays"]))),
        }
    )
    config["completion"]["residual_sugar_g_l"] = float(
        completion["drying_definition_residual_sugar_g_l"]
    )
    config["completion"]["minimum_probability"] = float(completion["minimum_probability"])
    config["sampling"]["samples_per_process"] = int(sampling["maximum_samples_per_process"])
    config["sampling"]["capture_intervals_per_process"] = int(
        sampling["maximum_capture_intervals_per_process"]
    )
    config["sampling"]["capture_stages_c"] = list(sampling["capture_stages_c"])
    config["operations"].update(
        {
            "sample_event_order": sampling["sample_event_order"],
            "minimum_minutes_between_sampling_and_nutrition": sampling[
                "minimum_minutes_between_sampling_and_nutrition"
            ],
            "manual_sampling_capacity_per_time_slot": sampling[
                "manual_sampling_capacity_per_time_slot"
            ],
        }
    )
    return config


def parameter_columns(config: dict[str, Any]) -> tuple[str, ...]:
    return tuple(f"aroma__{name}" for name in config["objective"]["priority_parameters"])


def log_parameter_matrix(ensemble: pd.DataFrame, config: dict[str, Any]) -> np.ndarray:
    columns = parameter_columns(config)
    values = ensemble.loc[:, columns].to_numpy(dtype=float)
    if np.any(values <= 0.0) or not np.isfinite(values).all():
        raise ValueError("Aroma ensemble parameters must be finite and positive")
    return np.log(values)


def representative_members(ensemble: pd.DataFrame, config: dict[str, Any]) -> list[int]:
    values = log_parameter_matrix(ensemble, config)
    primary_columns = [column for column in ensemble if column.startswith("primary__")]
    primary = np.log(ensemble[primary_columns].to_numpy(dtype=float))
    values = np.column_stack([values, primary])
    median = np.median(values, axis=0)
    scale = np.maximum(np.std(values, axis=0, ddof=1), 1e-8)
    standardized = (values - median) / scale
    target = int(config["objective"]["representative_ensemble_members"])
    selected = [int(np.argmin(np.linalg.norm(standardized, axis=1)))]
    while len(selected) < min(target, len(ensemble)):
        distance = np.min(
            np.stack(
                [np.linalg.norm(standardized - standardized[index], axis=1) for index in selected]
            ),
            axis=0,
        )
        distance[selected] = -np.inf
        selected.append(int(np.argmax(distance)))
    return selected


def prior_precision(ensemble: pd.DataFrame, config: dict[str, Any]) -> np.ndarray:
    values = log_parameter_matrix(ensemble, config)
    covariance = np.cov(values, rowvar=False)
    ridge = float(config["objective"]["prior_ridge"])
    covariance = 0.5 * (covariance + covariance.T) + np.eye(covariance.shape[0]) * ridge
    return np.linalg.pinv(covariance, rcond=1e-10)


def allowed_nutrition_times(config: dict[str, Any]) -> np.ndarray:
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    latest = float(config["nutrition"]["latest_h"])
    weekdays = set(int(value) for value in config["nutrition"]["manual_weekdays"])
    hours = [int(value) for value in config["nutrition"]["manual_hours_local"]]
    times = []
    date = start.date()
    while True:
        for hour in hours:
            current = datetime.combine(date, datetime.min.time(), tzinfo=start.tzinfo) + timedelta(hours=hour)
            relative = (current - start).total_seconds() / 3600.0
            if date.weekday() in weekdays and 0.0 <= relative <= latest:
                times.append(relative)
        if (datetime.combine(date, datetime.min.time(), tzinfo=start.tzinfo) - start).total_seconds() / 3600.0 > latest + 24.0:
            break
        date += timedelta(days=1)
    return np.asarray(sorted(set(times)), dtype=float)


def allowed_sampling_times(config: dict[str, Any]) -> np.ndarray:
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    horizon = float(config["future_process"]["information_horizon_h"])
    weekdays = set(int(value) for value in config["nutrition"]["manual_weekdays"])
    hours = [int(value) for value in config["nutrition"]["manual_hours_local"]]
    times = []
    date = start.date()
    while True:
        midnight = datetime.combine(date, datetime.min.time(), tzinfo=start.tzinfo)
        if (midnight - start).total_seconds() / 3600.0 > horizon + 24.0:
            break
        if date.weekday() in weekdays:
            for hour in hours:
                current = midnight + timedelta(hours=hour)
                relative = (current - start).total_seconds() / 3600.0
                if 0.0 <= relative <= horizon:
                    times.append(relative)
        date += timedelta(days=1)
    return np.asarray(sorted(set(times)), dtype=float)


def _nearest_unused_index(desired: int, used: set[int], size: int) -> int | None:
    candidates = sorted(range(size), key=lambda value: (abs(value - desired), value))
    return next((value for value in candidates if value not in used), None)


def normalize_pulses(
    pulse_times: np.ndarray, pulse_amounts: np.ndarray, config: dict[str, Any]
) -> tuple[tuple[float, float], ...]:
    """Canonicalize pulses without silently merging or discarding a dose."""

    allowed = allowed_nutrition_times(config)
    max_per = float(config["nutrition"]["maximum_yan_per_pulse_mg_l"])
    max_total = float(config["nutrition"]["maximum_total_yan_mg_l"])
    if len(allowed) == 0:
        raise ValueError("No legal nutrition times are configured")
    used: set[int] = set()
    items: list[tuple[float, float]] = []
    for time, amount in zip(pulse_times, pulse_amounts):
        desired = int(np.clip(round(float(time)), 0, len(allowed) - 1))
        value = float(np.clip(amount, 0.0, max_per))
        if value <= 0.0:
            continue
        chosen_index = _nearest_unused_index(desired, used, len(allowed))
        if chosen_index is None:
            raise ValueError("Duplicate nutrition times cannot be repaired without losing dose")
        used.add(chosen_index)
        items.append((float(allowed[chosen_index]), value))
    items.sort()
    total = sum(amount for _, amount in items)
    if total > max_total and total > 0.0:
        scale = max_total / total
        items = [(time, amount * scale) for time, amount in items]
    if len(items) > int(config["nutrition"]["maximum_pulses"]):
        raise ValueError("Pulse count exceeds configured maximum")
    return tuple(items)


def _design_layout(config: dict[str, Any]) -> dict[str, slice | int]:
    changes = int(config["future_process"]["maximum_temperature_changes"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    cursor = 1
    layout: dict[str, slice | int] = {"temperature_initial": 0}
    for name, length in (
        ("temperature_active", changes),
        ("temperature_slots", changes),
        ("temperature_levels", changes),
        ("nutrition_active", pulses),
        ("nutrition_times", pulses),
        ("nutrition_amounts", pulses),
    ):
        layout[name] = slice(cursor, cursor + length)
        cursor += length
    layout["size"] = cursor
    return layout


def policy_to_canonical_vector(policy: DesignPolicy, config: dict[str, Any]) -> np.ndarray:
    layout = _design_layout(config)
    vector = np.zeros(int(layout["size"]), dtype=float)
    slots = int(config["future_process"]["optimized_temperature_slots"])
    changes = int(config["future_process"]["maximum_temperature_changes"])
    profile = np.asarray(policy.temperature_c, dtype=float)
    if len(profile) != slots:
        raise ValueError(f"A canonical temperature profile must contain {slots} legal slots")
    vector[int(layout["temperature_initial"])] = float(profile[0])
    change_indices = np.flatnonzero(np.abs(np.diff(profile)) > 1e-12) + 1
    if len(change_indices) > changes:
        raise ValueError("Policy contains more effective temperature changes than configured")
    for output_index, slot_index in enumerate(change_indices):
        vector[layout["temperature_active"]][output_index] = 1.0
        vector[layout["temperature_slots"]][output_index] = float(slot_index)
        vector[layout["temperature_levels"]][output_index] = float(profile[slot_index])
    allowed = allowed_nutrition_times(config)
    if len(policy.nutrition_mg_yan_l) > int(config["nutrition"]["maximum_pulses"]):
        raise ValueError("Policy contains too many nutrition pulses")
    for index, (time_h, amount) in enumerate(sorted(policy.nutrition_mg_yan_l)):
        matches = np.flatnonzero(np.isclose(allowed, float(time_h), atol=1e-9, rtol=0.0))
        if len(matches) != 1:
            raise ValueError(f"Nutrition time {time_h} is not a unique legal slot")
        vector[layout["nutrition_active"]][index] = 1.0
        vector[layout["nutrition_times"]][index] = float(matches[0])
        vector[layout["nutrition_amounts"]][index] = float(amount)
    return vector


def decode_policy_vector(name: str, values: np.ndarray, config: dict[str, Any]) -> CanonicalDesign:
    raw = np.asarray(values, dtype=float).copy()
    layout = _design_layout(config)
    if len(raw) != int(layout["size"]):
        raise ValueError(f"Expected {layout['size']} design variables; found {len(raw)}")
    if not np.isfinite(raw).all():
        raise ValueError("Design vector must be finite")
    bounds = vector_bounds(config)
    clipped = np.clip(raw, bounds[:, 0], bounds[:, 1])
    repairs: list[dict[str, Any]] = []
    if not np.array_equal(raw, clipped):
        repairs.append({"action": "clip_to_bounds", "changed_coordinates": int(np.sum(raw != clipped))})
    slots = int(config["future_process"]["optimized_temperature_slots"])
    changes = int(config["future_process"]["maximum_temperature_changes"])
    initial = float(clipped[int(layout["temperature_initial"])])
    proposed_changes: list[tuple[int, float, int]] = []
    used_slots: set[int] = set()
    for index in range(changes):
        active = float(clipped[layout["temperature_active"]][index]) >= 0.5
        if not active:
            continue
        desired = int(np.clip(round(float(clipped[layout["temperature_slots"]][index])), 1, slots - 1))
        chosen = _nearest_unused_index(desired - 1, {value - 1 for value in used_slots}, slots - 1)
        if chosen is None:
            raise ValueError("Temperature changepoints cannot be made unique")
        chosen += 1
        if chosen != desired:
            repairs.append(
                {"action": "relocate_duplicate_temperature_changepoint", "entry": index, "from_slot": desired, "to_slot": chosen}
            )
        used_slots.add(chosen)
        proposed_changes.append((chosen, float(clipped[layout["temperature_levels"]][index]), index))
    profile = np.full(slots, initial, dtype=float)
    effective_changes: list[tuple[int, float]] = []
    current = initial
    for slot_index, level, source_index in sorted(proposed_changes):
        if math.isclose(level, current, abs_tol=1e-12, rel_tol=0.0):
            repairs.append({"action": "merge_equal_temperature_segment", "entry": source_index, "slot": slot_index})
            continue
        profile[slot_index:] = level
        current = level
        effective_changes.append((slot_index, level))

    allowed = allowed_nutrition_times(config)
    max_per = float(config["nutrition"]["maximum_yan_per_pulse_mg_l"])
    max_total = float(config["nutrition"]["maximum_total_yan_mg_l"])
    pulse_count = int(config["nutrition"]["maximum_pulses"])
    used_times: set[int] = set()
    pulse_rows: list[tuple[float, float, int]] = []
    for index in range(pulse_count):
        active = float(clipped[layout["nutrition_active"]][index]) >= 0.5
        amount = float(np.clip(clipped[layout["nutrition_amounts"]][index], 0.0, max_per))
        if not active:
            continue
        if amount <= 0.0:
            repairs.append({"action": "deactivate_zero_dose_pulse", "entry": index})
            continue
        desired = int(np.clip(round(float(clipped[layout["nutrition_times"]][index])), 0, len(allowed) - 1))
        chosen = _nearest_unused_index(desired, used_times, len(allowed))
        if chosen is None:
            raise ValueError("Duplicate nutrition times cannot be repaired without losing dose")
        if chosen != desired:
            repairs.append(
                {"action": "relocate_duplicate_nutrition_time", "entry": index, "from_index": desired, "to_index": chosen, "dose_preserved_mg_yan_l": amount}
            )
        used_times.add(chosen)
        pulse_rows.append((float(allowed[chosen]), amount, index))
    total_before = float(sum(row[1] for row in pulse_rows))
    if total_before > max_total and total_before > 0.0:
        scale = max_total / total_before
        pulse_rows = [(time_h, amount * scale, index) for time_h, amount, index in pulse_rows]
        repairs.append(
            {"action": "scale_total_yan_to_limit", "total_before_mg_yan_l": total_before, "total_after_mg_yan_l": max_total, "scale": scale}
        )
    schedule = tuple((time_h, amount) for time_h, amount, _ in sorted(pulse_rows))
    policy = DesignPolicy(name, tuple(float(value) for value in profile), schedule)
    canonical = policy_to_canonical_vector(policy, config)
    return CanonicalDesign(policy, raw, canonical, tuple(repairs))


def continuous_design_indices(canonical_vector: np.ndarray, config: dict[str, Any]) -> np.ndarray:
    """Return continuous coordinates while holding all mixed decisions fixed."""

    vector = np.asarray(canonical_vector, dtype=float)
    layout = _design_layout(config)
    if len(vector) != int(layout["size"]):
        raise ValueError("Canonical vector has incorrect length")
    indices = [int(layout["temperature_initial"])]
    for active_index, level_index in zip(
        range(layout["temperature_active"].start, layout["temperature_active"].stop),
        range(layout["temperature_levels"].start, layout["temperature_levels"].stop),
    ):
        if vector[active_index] >= 0.5:
            indices.append(level_index)
    for active_index, amount_index in zip(
        range(layout["nutrition_active"].start, layout["nutrition_active"].stop),
        range(layout["nutrition_amounts"].start, layout["nutrition_amounts"].stop),
    ):
        if vector[active_index] >= 0.5:
            indices.append(amount_index)
    return np.asarray(indices, dtype=int)


def anchor_policy(config: dict[str, Any]) -> DesignPolicy:
    slots = int(config["future_process"]["optimized_temperature_slots"])
    temperature = float(config["anchor"]["temperature_c"])
    return DesignPolicy(
        "anchor_18C",
        tuple([temperature] * slots),
        tuple((float(time), float(amount)) for time, amount in config["anchor"]["nutrition_pulses_mg_yan_l"]),
    )


def policy_from_vector(name: str, values: np.ndarray, config: dict[str, Any]) -> DesignPolicy:
    return decode_policy_vector(name, values, config).policy


def vector_bounds(config: dict[str, Any]) -> np.ndarray:
    slots = int(config["future_process"]["optimized_temperature_slots"])
    changes = int(config["future_process"]["maximum_temperature_changes"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    allowed = allowed_nutrition_times(config)
    minimum = float(config["future_process"]["temperature_minimum_c"])
    maximum = float(config["future_process"]["temperature_maximum_c"])
    return np.asarray(
        [[minimum, maximum]]
        + [[0.0, 1.0]] * changes
        + [[0.0, float(slots - 1)]] * changes
        + [[minimum, maximum]] * changes
        + [[0.0, 1.0]] * pulses
        + [[0.0, float(len(allowed) - 1)]] * pulses
        + [[0.0, config["nutrition"]["maximum_yan_per_pulse_mg_l"]]] * pulses,
        dtype=float,
    )


def temperature_profile_metrics(policy: DesignPolicy, config: dict[str, Any]) -> dict[str, float | int]:
    profile = np.asarray(policy.temperature_c, dtype=float)
    jumps = np.diff(profile)
    effective = jumps[np.abs(jumps) > 1e-12]
    slot_h = float(config["future_process"]["temperature_slot_h"])
    lower = float(config["future_process"]["temperature_minimum_c"])
    upper = float(config["future_process"]["temperature_maximum_c"])
    initial = float(config["future_process"]["temperature_actuator"]["initial_temperature_c"])
    return {
        "temperature_changes": int(len(effective)),
        "total_thermal_variation_c": float(np.sum(np.abs(effective))),
        "maximum_temperature_jump_c": float(np.max(np.abs(effective))) if len(effective) else 0.0,
        "time_at_temperature_limits_h": float(
            slot_h * np.sum(np.isclose(profile, lower) | np.isclose(profile, upper))
        ),
        "thermal_effort_proxy_c2_h": float(slot_h * np.sum((profile - initial) ** 2)),
    }


def design_complexity_penalty(policies: tuple[DesignPolicy, ...], config: dict[str, Any]) -> float:
    penalty = config["objective"]["design_complexity_penalty"]
    total = 0.0
    for policy in policies:
        metrics = temperature_profile_metrics(policy, config)
        total += float(penalty["per_temperature_change"]) * float(metrics["temperature_changes"])
        total += float(penalty["per_degree_total_variation"]) * float(
            metrics["total_thermal_variation_c"]
        )
    return float(total)


def _future_design(
    policy: DesignPolicy,
    config: dict[str, Any],
    actuator_scenario: dict[str, float] | None = None,
):
    model = pilot_calibration._model_module()
    horizon = float(config["future_process"]["maximum_horizon_h"])
    total_slots = int(round(horizon / float(config["future_process"]["temperature_slot_h"])))
    setpoints = list(policy.temperature_c)
    setpoints.extend([setpoints[-1]] * max(0, total_slots - len(setpoints)))
    actuator = config["future_process"]["temperature_actuator"]
    scenario = dict(actuator_scenario or {})
    actuator_step = float(actuator["simulation_step_h"])
    actuator_time = np.arange(0.0, horizon + actuator_step * 0.5, actuator_step)
    executed = np.empty_like(actuator_time)
    executed[0] = float(scenario.get("initial_temperature_c", actuator["initial_temperature_c"]))
    tau = float(scenario.get("tau_h", actuator["global_tau_h"]))
    command_delay_h = max(float(scenario.get("command_delay_h", 0.0)), 0.0)
    tracking_error_c = float(scenario.get("tracking_error_c", 0.0))
    probe_bias_c = float(scenario.get("probe_bias_c", 0.0))
    slot_h = float(config["future_process"]["temperature_slot_h"])
    for index in range(1, len(actuator_time)):
        command_time = max(float(actuator_time[index - 1]) - command_delay_h, 0.0)
        setpoint_index = min(int(command_time // slot_h), total_slots - 1)
        setpoint = float(setpoints[setpoint_index])
        decay = math.exp(-actuator_step / tau)
        executed[index] = setpoint + (executed[index - 1] - setpoint) * decay
    executed = executed + tracking_error_c + probe_bias_c
    pulses = tuple((time, amount / 1000.0) for time, amount in policy.nutrition_mg_yan_l)
    return model.BatchData(
        medium="natural_pilot_2026",
        batch=policy.name,
        time=actuator_time,
        temperature_c=executed,
        initials={key: float(value) for key, value in config["future_process"]["initials"].items()},
        pulses={"N": pulses, "G": tuple(), "F": tuple(), "E": tuple(), "X": tuple()},
        observations={},
    )


def _theta(member: pd.Series) -> dict[str, float]:
    model = pilot_calibration._model_module()
    theta = dict(model.DEFAULT_THETA)
    for name in tuple(theta):
        column = f"primary__{name}"
        if column in member and np.isfinite(member[column]):
            theta[name] = float(member[column])
    return theta


def _aroma_log_values(member: pd.Series, species: str) -> np.ndarray:
    return np.log(
        [float(member[f"aroma__{species}__{parameter}"]) for parameter in PARAMETER_NAMES]
    )


def _forcing_from_core(
    design,
    core: pd.DataFrame,
    theta: dict[str, float],
    partition: dict[str, float],
    species: str,
    initial: float,
) -> AromaForcing:
    model = pilot_calibration._model_module()
    growth, uptake, loss = [], [], []
    for time_h, row in core.iterrows():
        terms = model.kinetic_terms(
            theta,
            model.temperature_at(design, float(time_h)),
            float(row.X),
            float(row.N),
            float(row.G),
            float(row.F),
            float(row.E),
        )
        total = max(float(terms["sugar_total"]), 1e-8)
        maintenance = float(terms["maintenance"] * terms["maintenance_availability"])
        sugar_uptake = (
            theta["qXG"] * terms["growth_factor"]
            + theta["qEG"] * terms["glucose_ferm_factor"]
            + maintenance * float(row.G) / total
            + theta["qXF"] * terms["growth_factor"]
            + theta["qEF"] * terms["fructose_ferm_factor"]
            + maintenance * float(row.F) / total
        ) * float(row.X)
        ethanol_rate = max(float(terms["beta_g"] + terms["beta_f"]) * float(row.X), 0.0)
        co2_rate = (44.01 / (2.0 * 46.07)) * ethanol_rate
        k_part = _partition_basis(
            partition,
            model.temperature_at(design, float(time_h)),
            float(row.E),
            float(row.G + row.F),
        )
        growth.append(float(row.N) / (float(row.N) + 0.035))
        uptake.append(max(float(sugar_uptake), 0.0))
        loss.append(max(k_part * co2_rate, 0.0))
    return AromaForcing(
        experiment_id=str(getattr(design, "name", design.batch)),
        time_h=core.index.to_numpy(dtype=float),
        growth_fraction=np.asarray(growth),
        sugar_uptake_g_l_h=np.asarray(uptake),
        loss_basis_h_inv=np.asarray(loss),
        initial_concentration_ug_l=float(initial),
        volume_l=230.0,
        trap_efficiency=float(partition["trap_efficiency"]),
    )


def physical_observation_vector(
    forcing: AromaForcing,
    log_values: np.ndarray,
    sample_times: np.ndarray,
) -> np.ndarray:
    """Return physical wine concentrations followed by captured interval masses."""

    liquid, captured = simulate_aroma(forcing, log_values)
    liquid_pred = np.interp(sample_times, forcing.time_h, liquid)
    interval = []
    for start, end in zip(sample_times[:-1], sample_times[1:]):
        mass = np.interp(end, forcing.time_h, captured) - np.interp(start, forcing.time_h, captured)
        interval.append(float(mass))
    return np.concatenate([liquid_pred, np.asarray(interval, dtype=float)])


def nominal_observation_scale(
    nominal_prediction: np.ndarray,
    sample_count: int,
    config: dict[str, Any],
) -> np.ndarray:
    """Compute observational scales once at the nominal parameter centre."""

    prediction = np.asarray(nominal_prediction, dtype=float)
    wine = prediction[:sample_count]
    captured = prediction[sample_count:]
    error = config["observation_error_model"]
    wine_sigma = np.maximum(
        float(error["wine_minimum_sigma_ug_l"]),
        float(error["wine_relative_sigma"]) * np.abs(wine),
    )
    capture_sigma = np.maximum(
        float(error["condensate_minimum_sigma_ug"]),
        float(error["condensate_relative_sigma"]) * np.abs(captured),
    )
    sigma = np.concatenate([wine_sigma, capture_sigma])
    if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0.0):
        raise ValueError("Nominal observation scales must be finite and positive")
    return sigma


def scaled_observation_vector(prediction: np.ndarray, nominal_sigma: np.ndarray) -> np.ndarray:
    prediction = np.asarray(prediction, dtype=float)
    sigma = np.asarray(nominal_sigma, dtype=float)
    if prediction.shape != sigma.shape:
        raise ValueError("Prediction and nominal sigma must have identical shapes")
    return prediction / sigma


def finite_difference_sensitivity(
    predictor,
    centre: np.ndarray,
    step: float,
    nominal_sigma: np.ndarray,
    lower: np.ndarray | None = None,
    upper: np.ndarray | None = None,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Differentiate physical predictions while keeping nominal sigma frozen."""

    centre = np.asarray(centre, dtype=float)
    step = float(step)
    if step <= 0.0 or not np.isfinite(step):
        raise ValueError("Finite-difference step must be finite and positive")
    lower = np.full_like(centre, -np.inf) if lower is None else np.asarray(lower, dtype=float)
    upper = np.full_like(centre, np.inf) if upper is None else np.asarray(upper, dtype=float)
    if lower.shape != centre.shape or upper.shape != centre.shape:
        raise ValueError("Finite-difference bounds must match parameter centre")
    nominal = np.asarray(predictor(centre.copy()), dtype=float)
    sigma = np.asarray(nominal_sigma, dtype=float)
    if nominal.shape != sigma.shape:
        raise ValueError("Predictor output and nominal sigma must have identical shapes")
    columns: list[np.ndarray] = []
    methods: list[str] = []
    for index in range(len(centre)):
        can_minus = centre[index] - step >= lower[index] - 1e-14
        can_plus = centre[index] + step <= upper[index] + 1e-14
        if can_minus and can_plus:
            plus, minus = centre.copy(), centre.copy()
            plus[index] += step
            minus[index] -= step
            derivative = (np.asarray(predictor(plus)) - np.asarray(predictor(minus))) / (2.0 * step)
            method = "central"
        elif can_plus:
            h = min(step, (upper[index] - centre[index]) / 2.0)
            if not np.isfinite(h):
                h = step
            if h <= 0.0:
                raise ValueError("No feasible forward finite-difference step")
            one, two = centre.copy(), centre.copy()
            one[index] += h
            two[index] += 2.0 * h
            derivative = (
                -3.0 * nominal + 4.0 * np.asarray(predictor(one)) - np.asarray(predictor(two))
            ) / (2.0 * h)
            method = "forward_second_order"
        elif can_minus:
            h = min(step, (centre[index] - lower[index]) / 2.0)
            if not np.isfinite(h):
                h = step
            if h <= 0.0:
                raise ValueError("No feasible backward finite-difference step")
            one, two = centre.copy(), centre.copy()
            one[index] -= h
            two[index] -= 2.0 * h
            derivative = (
                3.0 * nominal - 4.0 * np.asarray(predictor(one)) + np.asarray(predictor(two))
            ) / (2.0 * h)
            method = "backward_second_order"
        else:
            raise ValueError("Parameter has no feasible finite-difference direction")
        columns.append(np.asarray(derivative, dtype=float) / sigma)
        methods.append(method)
    return np.column_stack(columns), tuple(methods)


def prepare_design(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
    *,
    simulation_grid_step_h: float | None = None,
    actuator_scenario: dict[str, float] | None = None,
) -> PreparedDesign | None:
    model = pilot_calibration._model_module()
    design = _future_design(policy, config, actuator_scenario)
    info_horizon = float(config["future_process"]["information_horizon_h"])
    sample_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    grid_step = float(
        simulation_grid_step_h
        if simulation_grid_step_h is not None
        else config["fim_validation"]["nominal_time_grid_step_h"]
    )
    event_times = np.asarray([time_h for time_h, _ in policy.nutrition_mg_yan_l], dtype=float)
    grid = np.unique(
        np.concatenate(
            [np.arange(0.0, info_horizon + grid_step * 0.5, grid_step), sample_times, event_times]
        )
    )
    theta = _theta(member)
    numerical_order = config["operations"].get("numerical_sample_event_order", "sample_before_action")
    core = model.simulate(design, theta, grid, sample_event_order=numerical_order)
    if core is None:
        return None
    drying_step = float(config["completion"]["drying_grid_step_h"])
    final_grid = np.unique(
        np.concatenate(
            [
                np.arange(
                    0.0,
                    float(config["future_process"]["maximum_horizon_h"]) + drying_step * 0.5,
                    drying_step,
                ),
                event_times,
            ]
        )
    )
    completion_core = model.simulate(
        design, theta, final_grid, sample_event_order=numerical_order
    )
    residual_sugar = (
        math.inf
        if completion_core is None
        else float(completion_core.iloc[-1]["G"] + completion_core.iloc[-1]["F"])
    )
    if completion_core is None:
        drying_time = math.inf
    else:
        sugar = completion_core["G"].to_numpy(dtype=float) + completion_core["F"].to_numpy(dtype=float)
        dry = np.flatnonzero(sugar <= float(config["completion"]["residual_sugar_g_l"]))
        if len(dry):
            upper_index = int(dry[0])
            if upper_index == 0:
                drying_time = float(completion_core.index[0])
            else:
                lower_index = upper_index - 1
                t0 = float(completion_core.index[lower_index])
                t1 = float(completion_core.index[upper_index])
                s0 = float(sugar[lower_index])
                s1 = float(sugar[upper_index])
                target = float(config["completion"]["residual_sugar_g_l"])
                fraction = 0.0 if math.isclose(s0, s1) else (s0 - target) / (s0 - s1)
                drying_time = float(t0 + np.clip(fraction, 0.0, 1.0) * (t1 - t0))
        else:
            drying_time = math.inf
    initial_aroma = {
        "ethyl_acetate": 1152.067155825675,
        "ethyl_octanoate": 3.48070575912168,
        "isoamyl_acetate": 2.5,
    }
    forcings = {}
    aroma_values = {}
    for species in SPECIES:
        forcings[species] = _forcing_from_core(
            design, core, theta, partitions[species], species, initial_aroma[species]
        )
        aroma_values[species] = _aroma_log_values(member, species)
    return PreparedDesign(
        policy,
        int(member["ensemble_member"]),
        forcings,
        aroma_values,
        residual_sugar,
        drying_time,
    )


def fim_components_from_prepared(
    prepared: PreparedDesign,
    sample_times: np.ndarray,
    config: dict[str, Any],
    *,
    finite_difference_log_step: float | None = None,
) -> FIMComponents:
    step = float(
        finite_difference_log_step
        if finite_difference_log_step is not None
        else config["objective"]["finite_difference_log_step"]
    )
    fim = np.zeros((9, 9), dtype=float)
    physical_blocks: list[np.ndarray] = []
    sigma_blocks: list[np.ndarray] = []
    scaled_blocks: list[np.ndarray] = []
    sensitivity_blocks: list[np.ndarray] = []
    methods: list[str] = []
    log_bounds = np.asarray(config["objective"]["finite_difference_log_parameter_bounds"], dtype=float)
    if log_bounds.shape != (3, 2):
        raise ValueError("finite_difference_log_parameter_bounds must be a 3-by-2 array")
    for species_index, species in enumerate(SPECIES):
        forcing = prepared.forcings[species]
        centre = prepared.aroma_log_values[species]
        predictor = lambda values, forcing=forcing: physical_observation_vector(
            forcing, values, sample_times
        )
        nominal = predictor(centre)
        sigma = nominal_observation_scale(nominal, len(sample_times), config)
        block, block_methods = finite_difference_sensitivity(
            predictor,
            centre,
            step,
            sigma,
            lower=log_bounds[:, 0],
            upper=log_bounds[:, 1],
        )
        start = 3 * species_index
        fim[start : start + 3, start : start + 3] = block.T @ block
        physical_blocks.append(nominal)
        sigma_blocks.append(sigma)
        scaled_blocks.append(scaled_observation_vector(nominal, sigma))
        sensitivity_blocks.append(block)
        methods.extend(f"{species}:{method}" for method in block_methods)
    rows_per_species = sensitivity_blocks[0].shape[0]
    sensitivity = np.zeros((rows_per_species * len(SPECIES), 9), dtype=float)
    for species_index, block in enumerate(sensitivity_blocks):
        row_start = species_index * rows_per_species
        col_start = species_index * 3
        sensitivity[row_start : row_start + rows_per_species, col_start : col_start + 3] = block
    fim = 0.5 * (fim + fim.T)
    return FIMComponents(
        physical_prediction=np.concatenate(physical_blocks),
        nominal_sigma=np.concatenate(sigma_blocks),
        scaled_observation=np.concatenate(scaled_blocks),
        sensitivity=sensitivity,
        fim=fim,
        difference_methods=tuple(methods),
    )


def fim_from_prepared(
    prepared: PreparedDesign,
    sample_times: np.ndarray,
    config: dict[str, Any],
    *,
    finite_difference_log_step: float | None = None,
) -> np.ndarray:
    return fim_components_from_prepared(
        prepared,
        sample_times,
        config,
        finite_difference_log_step=finite_difference_log_step,
    ).fim


def prepare_sampling_sensitivity_cache(
    prepared: PreparedDesign,
    config: dict[str, Any],
    *,
    finite_difference_log_step: float | None = None,
) -> SamplingSensitivityCache:
    """Precompute physical trajectory derivatives for fast discrete sampling search."""

    step = float(
        finite_difference_log_step
        if finite_difference_log_step is not None
        else config["objective"]["finite_difference_log_step"]
    )
    log_bounds = np.asarray(config["objective"]["finite_difference_log_parameter_bounds"], dtype=float)
    nominal_liquid: dict[str, np.ndarray] = {}
    nominal_captured: dict[str, np.ndarray] = {}
    liquid_derivatives: dict[str, np.ndarray] = {}
    captured_derivatives: dict[str, np.ndarray] = {}
    methods: dict[str, tuple[str, ...]] = {}
    for species in SPECIES:
        forcing = prepared.forcings[species]
        centre = prepared.aroma_log_values[species]

        def trajectory(values):
            liquid, captured = simulate_aroma(forcing, values)
            return np.concatenate([liquid, captured])

        nominal = trajectory(centre)
        derivatives, method = finite_difference_sensitivity(
            trajectory,
            centre,
            step,
            np.ones_like(nominal),
            lower=log_bounds[:, 0],
            upper=log_bounds[:, 1],
        )
        size = len(forcing.time_h)
        nominal_liquid[species] = nominal[:size]
        nominal_captured[species] = nominal[size:]
        liquid_derivatives[species] = derivatives[:size, :]
        captured_derivatives[species] = derivatives[size:, :]
        methods[species] = method
    return SamplingSensitivityCache(
        prepared,
        nominal_liquid,
        nominal_captured,
        liquid_derivatives,
        captured_derivatives,
        methods,
    )


def fim_from_sampling_sensitivity_cache(
    cache: SamplingSensitivityCache,
    sample_times: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    sample_times = np.asarray(sample_times, dtype=float)
    fim = np.zeros((9, 9), dtype=float)
    for species_index, species in enumerate(SPECIES):
        forcing = cache.prepared.forcings[species]
        time = forcing.time_h
        liquid_nominal = np.interp(sample_times, time, cache.nominal_liquid[species])
        captured_nominal = np.interp(sample_times, time, cache.nominal_captured[species])
        interval_nominal = np.diff(captured_nominal)
        physical = np.concatenate([liquid_nominal, interval_nominal])
        sigma = nominal_observation_scale(physical, len(sample_times), config)
        columns = []
        for parameter_index in range(3):
            liquid = np.interp(
                sample_times,
                time,
                cache.liquid_derivatives[species][:, parameter_index],
            )
            captured = np.interp(
                sample_times,
                time,
                cache.captured_derivatives[species][:, parameter_index],
            )
            columns.append(np.concatenate([liquid, np.diff(captured)]) / sigma)
        block = np.column_stack(columns)
        start = 3 * species_index
        fim[start : start + 3, start : start + 3] = block.T @ block
    return 0.5 * (fim + fim.T)


def design_fim(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
    *,
    actuator_scenario: dict[str, float] | None = None,
) -> tuple[np.ndarray, float]:
    fim, residual, _drying = design_fim_with_drying(
        policy,
        member,
        config,
        partitions,
        actuator_scenario=actuator_scenario,
    )
    return fim, residual


def design_fim_with_drying(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
    *,
    actuator_scenario: dict[str, float] | None = None,
) -> tuple[np.ndarray, float, float]:
    prepared = prepare_design(
        policy, member, config, partitions, actuator_scenario=actuator_scenario
    )
    if prepared is None:
        return np.zeros((9, 9), dtype=float), math.inf, -math.inf
    sample_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    return (
        fim_from_prepared(prepared, sample_times, config),
        prepared.residual_sugar_g_l,
        prepared.drying_time_h,
    )


def logdet(matrix: np.ndarray) -> float:
    sign, value = np.linalg.slogdet(0.5 * (matrix + matrix.T))
    return float(value) if sign > 0 and np.isfinite(value) else -math.inf


def evaluate_campaign(
    policies: tuple[DesignPolicy, ...],
    ensemble: pd.DataFrame,
    representative: list[int],
    prior: np.ndarray,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
    *,
    actuator_scenario: dict[str, float] | None = None,
    design_cache: dict[tuple, tuple[np.ndarray, float, float]] | None = None,
) -> tuple[float, list[ScenarioEvaluation]]:
    base_logdet = logdet(prior)
    evaluations = []
    for member_index in representative:
        member = ensemble.iloc[int(member_index)]
        total = prior.copy()
        completed = True
        max_residual = 0.0
        minimum_drying = math.inf
        latest_action = 0.0
        for policy in policies:
            scenario_key = tuple(sorted((actuator_scenario or {}).items()))
            cache_key = (
                int(member.get("ensemble_member", member_index)),
                tuple(policy.temperature_c),
                tuple(policy.nutrition_mg_yan_l),
                scenario_key,
            )
            cached = None if design_cache is None else design_cache.get(cache_key)
            if cached is None:
                cached = design_fim_with_drying(
                    policy,
                    member,
                    config,
                    partitions,
                    actuator_scenario=actuator_scenario,
                )
                if design_cache is not None:
                    design_cache[cache_key] = cached
            fim, residual, drying_time = cached
            total += fim
            max_residual = max(max_residual, residual)
            effective_change_slots = np.flatnonzero(
                np.abs(np.diff(policy.temperature_c)) > 1e-12
            ) + 1
            action_times = [
                float(config["future_process"]["temperature_slot_h"]) * int(slot)
                for slot in effective_change_slots
            ]
            action_times.extend(time_h for time_h, _ in policy.nutrition_mg_yan_l)
            policy_latest_action = max(action_times, default=0.0)
            latest_action = max(latest_action, policy_latest_action)
            minimum_drying = min(minimum_drying, drying_time)
            approved_margin = config["operations"].get("minimum_action_to_drying_margin_h")
            numerical_margin = 0.0 if approved_margin is None else float(approved_margin)
            completed = (
                completed
                and residual <= float(config["completion"]["residual_sugar_g_l"])
                and policy_latest_action <= drying_time - numerical_margin + 1e-9
            )
        evaluations.append(
            ScenarioEvaluation(
                int(member["ensemble_member"]),
                logdet(total) - base_logdet,
                completed,
                max_residual,
                total - prior,
                minimum_drying,
                latest_action,
            )
        )
    gains = np.asarray([item.information_gain for item in evaluations], dtype=float)
    completion_probability = float(np.mean([item.completion for item in evaluations]))
    score = (
        float(config["objective"]["median_weight"]) * float(np.median(gains))
        + float(config["objective"]["lower_decile_weight"]) * float(np.quantile(gains, 0.1))
    )
    shortfall = max(float(config["completion"]["minimum_probability"]) - completion_probability, 0.0)
    score -= float(config["completion"]["penalty_per_failed_fraction"]) * shortfall
    score -= design_complexity_penalty(policies, config)
    return score, evaluations


def robust_information_metrics(
    gains: np.ndarray,
    posterior_fims: list[np.ndarray],
    reference_gains: np.ndarray | None,
    config: dict[str, Any],
) -> dict[str, float]:
    gains = np.asarray(gains, dtype=float)
    tail_probability = float(config["objective"]["tail_probability"])
    tail_count = max(1, int(math.ceil(tail_probability * len(gains))))
    cvar = float(np.mean(np.sort(gains)[:tail_count]))
    minimum_eigenvalues = []
    conditions = []
    for fim in posterior_fims:
        eigenvalues = np.linalg.eigvalsh(0.5 * (fim + fim.T))
        minimum_eigenvalues.append(float(eigenvalues[0]))
        positive = eigenvalues[eigenvalues > 1e-12]
        conditions.append(float(eigenvalues[-1] / positive[0]) if len(positive) else math.inf)
    metrics = {
        "robust_score": float(config["objective"]["median_weight"]) * float(np.median(gains))
        + float(config["objective"]["lower_decile_weight"]) * float(np.quantile(gains, 0.1)),
        "median": float(np.median(gains)),
        "q10": float(np.quantile(gains, 0.1)),
        "minimum": float(np.min(gains)),
        "tail_cvar": cvar,
        "minimum_posterior_fim_eigenvalue": float(np.min(minimum_eigenvalues)),
        "maximum_posterior_fim_condition_number": float(np.max(conditions)),
    }
    if reference_gains is not None:
        reference = np.asarray(reference_gains, dtype=float)
        if reference.shape != gains.shape:
            raise ValueError("Reference gains must match candidate gains")
        delta = gains - reference
        metrics.update(
            {
                "fraction_exceeding_reference": float(np.mean(delta >= 0.0)),
                "worst_paired_loss_vs_reference": float(np.min(delta)),
            }
        )
    return metrics
