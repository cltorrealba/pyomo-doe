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
class ScenarioEvaluation:
    member: int
    information_gain: float
    completion: bool
    residual_sugar_g_l: float
    fim: np.ndarray


@dataclass(frozen=True)
class PreparedDesign:
    policy: DesignPolicy
    member: int
    forcings: dict[str, AromaForcing]
    aroma_log_values: dict[str, np.ndarray]
    residual_sugar_g_l: float
    drying_time_h: float


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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


def normalize_pulses(
    pulse_times: np.ndarray, pulse_amounts: np.ndarray, config: dict[str, Any]
) -> tuple[tuple[float, float], ...]:
    allowed = allowed_nutrition_times(config)
    max_per = float(config["nutrition"]["maximum_yan_per_pulse_mg_l"])
    max_total = float(config["nutrition"]["maximum_total_yan_mg_l"])
    threshold = float(config["nutrition"]["zero_pulse_threshold_mg_l"])
    merged: dict[float, float] = {}
    for time, amount in zip(pulse_times, pulse_amounts):
        chosen = float(allowed[int(np.clip(round(float(time)), 0, len(allowed) - 1))])
        value = float(np.clip(amount, 0.0, max_per))
        merged[chosen] = min(max_per, merged.get(chosen, 0.0) + value)
    items = sorted((time, amount) for time, amount in merged.items() if amount >= threshold)
    total = sum(amount for _, amount in items)
    if total > max_total and total > 0.0:
        scale = max_total / total
        items = [(time, amount * scale) for time, amount in items]
    return tuple(items[: int(config["nutrition"]["maximum_pulses"])])


def anchor_policy(config: dict[str, Any]) -> DesignPolicy:
    slots = int(config["future_process"]["optimized_temperature_slots"])
    temperature = float(config["anchor"]["temperature_c"])
    return DesignPolicy(
        "anchor_18C",
        tuple([temperature] * slots),
        tuple((float(time), float(amount)) for time, amount in config["anchor"]["nutrition_pulses_mg_yan_l"]),
    )


def policy_from_vector(name: str, values: np.ndarray, config: dict[str, Any]) -> DesignPolicy:
    values = np.asarray(values, dtype=float)
    slots = int(config["future_process"]["optimized_temperature_slots"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    expected = slots + 2 * pulses
    if len(values) != expected:
        raise ValueError(f"Expected {expected} design variables; found {len(values)}")
    temperature = np.clip(
        values[:slots],
        float(config["future_process"]["temperature_minimum_c"]),
        float(config["future_process"]["temperature_maximum_c"]),
    )
    schedule = normalize_pulses(
        values[slots : slots + pulses], values[slots + pulses :], config
    )
    return DesignPolicy(name, tuple(float(value) for value in temperature), schedule)


def vector_bounds(config: dict[str, Any]) -> np.ndarray:
    slots = int(config["future_process"]["optimized_temperature_slots"])
    pulses = int(config["nutrition"]["maximum_pulses"])
    allowed = allowed_nutrition_times(config)
    return np.asarray(
        [[config["future_process"]["temperature_minimum_c"], config["future_process"]["temperature_maximum_c"]]] * slots
        + [[0.0, float(len(allowed) - 1)]] * pulses
        + [[0.0, config["nutrition"]["maximum_yan_per_pulse_mg_l"]]] * pulses,
        dtype=float,
    )


def _future_design(policy: DesignPolicy, config: dict[str, Any]):
    model = pilot_calibration._model_module()
    horizon = float(config["future_process"]["maximum_horizon_h"])
    total_slots = int(round(horizon / float(config["future_process"]["temperature_slot_h"])))
    setpoints = list(policy.temperature_c)
    setpoints.extend([setpoints[-1]] * max(0, total_slots - len(setpoints)))
    actuator = config["future_process"]["temperature_actuator"]
    actuator_step = float(actuator["simulation_step_h"])
    actuator_time = np.arange(0.0, horizon + actuator_step * 0.5, actuator_step)
    executed = np.empty_like(actuator_time)
    executed[0] = float(actuator["initial_temperature_c"])
    tau = float(actuator["global_tau_h"])
    slot_h = float(config["future_process"]["temperature_slot_h"])
    for index in range(1, len(actuator_time)):
        setpoint_index = min(int(actuator_time[index - 1] // slot_h), total_slots - 1)
        setpoint = float(setpoints[setpoint_index])
        decay = math.exp(-actuator_step / tau)
        executed[index] = setpoint + (executed[index - 1] - setpoint) * decay
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


def _observation_vector(
    forcing: AromaForcing,
    log_values: np.ndarray,
    sample_times: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    liquid, captured = simulate_aroma(forcing, log_values)
    liquid_pred = np.interp(sample_times, forcing.time_h, liquid)
    wine_sigma = np.maximum(2.0, 0.2 * np.maximum(liquid_pred, 2.0))
    wine_scaled = liquid_pred / wine_sigma
    interval = []
    for start, end in zip(sample_times[:-1], sample_times[1:]):
        mass = np.interp(end, forcing.time_h, captured) - np.interp(start, forcing.time_h, captured)
        sigma = max(1.0, 0.25 * max(float(mass), 1.0))
        interval.append(float(mass) / sigma)
    return np.concatenate([wine_scaled, np.asarray(interval)])


def prepare_design(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
) -> PreparedDesign | None:
    model = pilot_calibration._model_module()
    design = _future_design(policy, config)
    info_horizon = float(config["future_process"]["information_horizon_h"])
    sample_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    grid = np.unique(
        np.concatenate([np.arange(0.0, info_horizon + 0.001, 2.0), sample_times])
    )
    theta = _theta(member)
    core = model.simulate(design, theta, grid)
    if core is None:
        return None
    final_grid = np.arange(0.0, float(config["future_process"]["maximum_horizon_h"]) + 0.001, 6.0)
    completion_core = model.simulate(design, theta, final_grid)
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
        drying_time = float(completion_core.index[int(dry[0])]) if len(dry) else math.inf
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


def fim_from_prepared(
    prepared: PreparedDesign,
    sample_times: np.ndarray,
    config: dict[str, Any],
) -> np.ndarray:
    step = float(config["objective"]["finite_difference_log_step"])
    fim = np.zeros((9, 9), dtype=float)
    for species_index, species in enumerate(SPECIES):
        forcing = prepared.forcings[species]
        centre = prepared.aroma_log_values[species]
        columns = []
        for index in range(3):
            plus, minus = centre.copy(), centre.copy()
            plus[index] += step
            minus[index] -= step
            columns.append(
                (
                    _observation_vector(forcing, plus, sample_times, config)
                    - _observation_vector(forcing, minus, sample_times, config)
                )
                / (2.0 * step)
            )
        block = np.column_stack(columns)
        start = 3 * species_index
        fim[start : start + 3, start : start + 3] = block.T @ block
    return 0.5 * (fim + fim.T)


def design_fim(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
) -> tuple[np.ndarray, float]:
    prepared = prepare_design(policy, member, config, partitions)
    if prepared is None:
        return np.zeros((9, 9), dtype=float), math.inf
    sample_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    return fim_from_prepared(prepared, sample_times, config), prepared.residual_sugar_g_l


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
) -> tuple[float, list[ScenarioEvaluation]]:
    base_logdet = logdet(prior)
    evaluations = []
    for member_index in representative:
        member = ensemble.iloc[int(member_index)]
        total = prior.copy()
        completed = True
        max_residual = 0.0
        for policy in policies:
            fim, residual = design_fim(policy, member, config, partitions)
            total += fim
            max_residual = max(max_residual, residual)
            completed = completed and residual <= float(config["completion"]["residual_sugar_g_l"])
        evaluations.append(
            ScenarioEvaluation(
                int(member["ensemble_member"]),
                logdet(total) - base_logdet,
                completed,
                max_residual,
                total - prior,
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
    return score, evaluations
