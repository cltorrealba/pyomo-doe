from __future__ import annotations

"""Nested production/observation models for the pilot-2026 aroma data.

The module deliberately keeps the published Morakul/Mouret equilibrium loss
driver fixed.  Candidate models add either a smooth post-nutrient production
term, a mass-conserving gas-line reservoir, or both.  This makes the competing
explanations testable without allowing a fitted transfer coefficient to absorb
errors in biological production.

The ``ASSUMED_COMPLETE_*`` variants encode the pilot train's engineering
contract explicitly: all aroma mass drained from the gas-line state is assigned
to the combined A+B condensate.  This is an assumption, not an independently
measured condenser efficiency, so those variants contain no fitted recovery
parameter.
"""

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import log_ndtr
from scipy.stats import qmc

from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma


BASELINE = "equilibrium_baseline"
DELAYED = "delayed_biological_response"
RESERVOIR = "line_reservoir"
COMBINED = "delayed_response_plus_reservoir"
CAPTURE = "fitted_capture_efficiency"
DELAYED_CAPTURE = "delayed_response_plus_capture_efficiency"
RESERVOIR_CAPTURE = "reservoir_plus_capture_efficiency"
FULL_CAPTURE = "delayed_reservoir_plus_capture_efficiency"
ETHANOL_CAPTURE = "ethanol_dependent_capture_efficiency"
ETHANOL_CAPTURE_RESERVOIR = "ethanol_capture_plus_line_reservoir"
ASSUMED_COMPLETE = "assumed_complete_capture_equilibrium"
ASSUMED_COMPLETE_DELAYED = "assumed_complete_capture_delayed_response"
ASSUMED_COMPLETE_RESERVOIR = "assumed_complete_capture_line_reservoir"
ASSUMED_COMPLETE_COMBINED = (
    "assumed_complete_capture_delayed_response_plus_line_reservoir"
)

ASSUMED_COMPLETE_MODELS = {
    ASSUMED_COMPLETE,
    ASSUMED_COMPLETE_DELAYED,
    ASSUMED_COMPLETE_RESERVOIR,
    ASSUMED_COMPLETE_COMBINED,
}

PARAMETER_NAMES_BY_MODEL: dict[str, tuple[str, ...]] = {
    BASELINE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
    ),
    DELAYED: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
    ),
    RESERVOIR: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "line_reservoir_tau_h",
    ),
    COMBINED: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
        "line_reservoir_tau_h",
    ),
    CAPTURE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "capture_efficiency_fraction",
    ),
    DELAYED_CAPTURE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
        "capture_efficiency_fraction",
    ),
    RESERVOIR_CAPTURE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "line_reservoir_tau_h",
        "capture_efficiency_fraction",
    ),
    FULL_CAPTURE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
        "line_reservoir_tau_h",
        "capture_efficiency_fraction",
    ),
    ETHANOL_CAPTURE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "capture_efficiency_fraction",
        "ethanol_capture_multiplier_per_10_g_l",
    ),
    ETHANOL_CAPTURE_RESERVOIR: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "capture_efficiency_fraction",
        "ethanol_capture_multiplier_per_10_g_l",
        "line_reservoir_tau_h",
    ),
    ASSUMED_COMPLETE: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
    ),
    ASSUMED_COMPLETE_DELAYED: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
    ),
    ASSUMED_COMPLETE_RESERVOIR: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "line_reservoir_tau_h",
    ),
    ASSUMED_COMPLETE_COMBINED: (
        "formation_growth_ug_per_g_sugar",
        "formation_stationary_ug_per_g_sugar",
        "pulse_formation_ug_per_g_sugar",
        "activation_peak_h",
        "line_reservoir_tau_h",
    ),
}


@dataclass(frozen=True)
class NestedSimulation:
    time_h: np.ndarray
    liquid_ug_l: np.ndarray
    captured_ug: np.ndarray
    line_inventory_ug: np.ndarray
    cumulative_production_ug: np.ndarray
    cumulative_volatilized_ug: np.ndarray
    cumulative_drained_ug: np.ndarray
    production_rate_ug_l_h: np.ndarray
    pulse_activation: np.ndarray
    capture_efficiency_fraction: np.ndarray
    relative_mass_balance_error: np.ndarray


@dataclass(frozen=True)
class NestedAromaFit:
    species: str
    model_variant: str
    parameter_names: tuple[str, ...]
    log_values: np.ndarray
    parameter_table: pd.DataFrame
    multistart_summary: pd.DataFrame
    covariance: pd.DataFrame
    wine_predictions: pd.DataFrame
    condensate_predictions: pd.DataFrame
    validation: dict[str, Any]


def pulse_activation(
    time_h: np.ndarray, pulse_time_h: float, activation_peak_h: float
) -> np.ndarray:
    """Unit-height, smooth gamma-shaped activation after a nutrient pulse."""

    time = np.asarray(time_h, dtype=float)
    tau = max(float(activation_peak_h), 1e-12)
    elapsed = time - float(pulse_time_h)
    activation = np.zeros_like(elapsed)
    active = elapsed > 0.0
    scaled = elapsed[active] / tau
    activation[active] = scaled * np.exp(1.0 - scaled)
    return activation


def _parameter_values(
    log_values: np.ndarray, model_variant: str
) -> dict[str, float]:
    if model_variant not in PARAMETER_NAMES_BY_MODEL:
        raise ValueError(f"Unknown nested aroma model {model_variant!r}")
    names = PARAMETER_NAMES_BY_MODEL[model_variant]
    values = np.exp(np.asarray(log_values, dtype=float))
    if len(values) != len(names):
        raise ValueError(
            f"{model_variant} expects {len(names)} parameters, got {len(values)}"
        )
    return dict(zip(names, values))


def simulate_nested(
    forcing: aroma.AromaForcing,
    pulse_time_h: float,
    log_values: np.ndarray,
    model_variant: str,
) -> NestedSimulation:
    """Simulate liquid aroma and captured mass with explicit mass accounting."""

    if forcing.loss_model != aroma.LOSS_MODEL_EQUILIBRIUM:
        raise ValueError("Nested validation requires the fixed equilibrium loss model")
    parameters = _parameter_values(log_values, model_variant)
    time = np.asarray(forcing.time_h, dtype=float)
    if len(time) < 2 or np.any(np.diff(time) <= 0.0):
        raise ValueError("Forcing time grid must be strictly increasing")

    delayed = model_variant in {
        DELAYED,
        COMBINED,
        DELAYED_CAPTURE,
        FULL_CAPTURE,
        ASSUMED_COMPLETE_DELAYED,
        ASSUMED_COMPLETE_COMBINED,
    }
    reservoir = model_variant in {
        RESERVOIR,
        COMBINED,
        RESERVOIR_CAPTURE,
        FULL_CAPTURE,
        ETHANOL_CAPTURE_RESERVOIR,
        ASSUMED_COMPLETE_RESERVOIR,
        ASSUMED_COMPLETE_COMBINED,
    }
    activation_peak = parameters.get("activation_peak_h", 1.0)
    activation = pulse_activation(time, pulse_time_h, activation_peak)
    pulse_yield = parameters.get("pulse_formation_ug_per_g_sugar", 0.0)
    growth_yield = parameters["formation_growth_ug_per_g_sugar"]
    stationary_yield = parameters["formation_stationary_ug_per_g_sugar"]
    phi = np.clip(np.asarray(forcing.growth_fraction, dtype=float), 0.0, 1.0)
    uptake = np.maximum(np.asarray(forcing.sugar_uptake_g_l_h, dtype=float), 0.0)
    production_rate = uptake * (
        growth_yield * phi
        + stationary_yield * (1.0 - phi)
        + (pulse_yield * activation if delayed else 0.0)
    )
    # build_forcings stores K(T,E)*Qgas/VL directly in loss_basis_h_inv for
    # the equilibrium model.  Reading that array avoids constructing a pandas
    # diagnostics table at every optimizer evaluation.
    loss = np.maximum(np.asarray(forcing.loss_basis_h_inv, dtype=float), 0.0)

    n_time = len(time)
    liquid = np.empty(n_time, dtype=float)
    captured = np.zeros(n_time, dtype=float)
    line = np.zeros(n_time, dtype=float)
    cumulative_production = np.zeros(n_time, dtype=float)
    cumulative_volatilized = np.zeros(n_time, dtype=float)
    cumulative_drained = np.zeros(n_time, dtype=float)
    relative_balance = np.zeros(n_time, dtype=float)
    liquid[0] = max(float(forcing.initial_concentration_ug_l), 0.0)
    volume = float(forcing.volume_l)
    initial_mass = liquid[0] * volume
    line_tau = parameters.get("line_reservoir_tau_h", 1.0)
    capture_efficiency_reference = (
        1.0
        if model_variant in ASSUMED_COMPLETE_MODELS
        else parameters.get(
            "capture_efficiency_fraction", float(forcing.trap_efficiency)
        )
    )
    if not 0.0 < capture_efficiency_reference <= 1.0:
        raise ValueError("Capture efficiency must be in (0, 1]")
    if model_variant in {ETHANOL_CAPTURE, ETHANOL_CAPTURE_RESERVOIR}:
        if forcing.ethanol_g_l is None:
            raise ValueError("Ethanol-dependent capture requires ethanol forcing")
        ethanol = np.asarray(forcing.ethanol_g_l, dtype=float)
        multiplier = parameters["ethanol_capture_multiplier_per_10_g_l"]
        capture_efficiency = np.clip(
            capture_efficiency_reference
            * np.exp(((ethanol - 50.0) / 10.0) * np.log(multiplier)),
            1e-12,
            1.0,
        )
    else:
        capture_efficiency = np.full(n_time, capture_efficiency_reference)

    for index in range(1, n_time):
        dt = float(time[index] - time[index - 1])
        production = float(production_rate[index - 1])
        loss_coefficient = max(float(loss[index - 1]), 0.0)
        previous = max(float(liquid[index - 1]), 0.0)
        produced_per_l = production * dt
        if loss_coefficient > 1e-12:
            decay = math.exp(-loss_coefficient * dt)
            current = previous * decay + production / loss_coefficient * (1.0 - decay)
        else:
            current = previous + produced_per_l
        current = max(float(current), 0.0)
        # Roundoff in the exact first-order update can place current a few ULPs
        # above the no-loss upper bound when production and loss are both tiny.
        # Enforce that physical bound before deriving volatilized mass.
        current = min(current, previous + produced_per_l)
        volatilized_per_l = produced_per_l - (current - previous)
        volatilized_mass = volatilized_per_l * volume

        if reservoir:
            decay_line = math.exp(-dt / line_tau)
            input_rate = volatilized_mass / dt
            end_inventory = (
                line[index - 1] * decay_line
                + input_rate * line_tau * (1.0 - decay_line)
            )
            drained_mass = max(
                line[index - 1] + volatilized_mass - end_inventory, 0.0
            )
            line[index] = max(float(end_inventory), 0.0)
        else:
            drained_mass = volatilized_mass
            line[index] = 0.0

        liquid[index] = current
        captured[index] = (
            captured[index - 1]
            + float(capture_efficiency[index - 1]) * drained_mass
        )
        cumulative_production[index] = (
            cumulative_production[index - 1] + produced_per_l * volume
        )
        cumulative_volatilized[index] = (
            cumulative_volatilized[index - 1] + volatilized_mass
        )
        cumulative_drained[index] = cumulative_drained[index - 1] + drained_mass
        lhs = initial_mass + cumulative_production[index]
        rhs = liquid[index] * volume + line[index] + cumulative_drained[index]
        relative_balance[index] = abs(lhs - rhs) / max(abs(lhs), 1.0)

    return NestedSimulation(
        time_h=time,
        liquid_ug_l=liquid,
        captured_ug=captured,
        line_inventory_ug=line,
        cumulative_production_ug=cumulative_production,
        cumulative_volatilized_ug=cumulative_volatilized,
        cumulative_drained_ug=cumulative_drained,
        production_rate_ug_l_h=production_rate,
        pulse_activation=activation,
        capture_efficiency_fraction=capture_efficiency,
        relative_mass_balance_error=relative_balance,
    )


def _sigma(value: float, relative: float, floor: float) -> float:
    return max(float(floor), float(relative) * max(abs(float(value)), float(floor)))


def prediction_tables_nested(
    forcings: dict[str, aroma.AromaForcing],
    pulse_times_h: dict[str, float],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    log_values: np.ndarray,
    model_variant: str,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, NestedSimulation]]:
    simulations = {
        run: simulate_nested(
            forcing, pulse_times_h[run], log_values, model_variant
        )
        for run, forcing in forcings.items()
    }
    error = config["error_model"]
    wine_rows: list[dict[str, Any]] = []
    for row in wine.itertuples(index=False):
        run = str(row.experiment_id)
        forcing = forcings[run]
        simulation = simulations[run]
        predicted = float(
            np.interp(float(row.time_h), forcing.time_h, simulation.liquid_ug_l)
        )
        status = str(row.model_observation_type)
        reference = (
            float(row.observed_value) if status == "observed" else float(row.upper_bound)
        )
        wine_rows.append(
            {
                "experiment_id": run,
                "sample_id": str(row.sample_id),
                "time_h": float(row.time_h),
                "status": status,
                "observed_or_upper_bound": reference,
                "predicted_ug_l": predicted,
                "sigma": _sigma(
                    reference,
                    error["wine_relative_sigma"],
                    error["wine_minimum_sigma_ug_l"],
                ),
                "initialization_point": bool(
                    np.isclose(float(row.time_h), float(wine[wine["experiment_id"].astype(str).eq(run)]["time_h"].min()))
                ),
            }
        )

    condensate_rows: list[dict[str, Any]] = []
    for row in condensate.itertuples(index=False):
        run = str(row.experiment_id)
        forcing = forcings[run]
        simulation = simulations[run]
        end = float(row.time_h)
        start = end - float(row.capture_interval_h)
        predicted = float(
            np.interp(end, forcing.time_h, simulation.captured_ug)
            - np.interp(start, forcing.time_h, simulation.captured_ug)
        )
        status = str(row.mix_model_observation_type)
        reference = (
            float(row.observed_value) if status == "observed" else float(row.upper_bound)
        )
        condensate_rows.append(
            {
                "experiment_id": run,
                "mix_id": str(row.mix_id),
                "interval_start_h": start,
                "interval_end_h": end,
                "status": status,
                "observed_or_upper_bound": reference,
                "predicted_captured_ug": max(predicted, 0.0),
                "sigma": _sigma(
                    reference,
                    error["condensate_relative_sigma"],
                    error["condensate_minimum_sigma_ug"],
                ),
            }
        )
    return pd.DataFrame(wine_rows), pd.DataFrame(condensate_rows), simulations


def residual_vector_nested(
    log_values: np.ndarray,
    forcings: dict[str, aroma.AromaForcing],
    pulse_times_h: dict[str, float],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    model_variant: str,
    config: dict[str, Any],
    default_log: np.ndarray,
    *,
    include_prior: bool = True,
) -> np.ndarray:
    wine_prediction, condensate_prediction, _ = prediction_tables_nested(
        forcings,
        pulse_times_h,
        wine,
        condensate,
        log_values,
        model_variant,
        config,
    )
    residuals: list[np.ndarray] = []
    calibration_domains = set(
        config["error_model"].get(
            "calibration_domains", ["wine", "condensate"]
        )
    )
    unknown_domains = calibration_domains - {"wine", "condensate"}
    if unknown_domains:
        raise ValueError(
            f"Unknown aroma calibration domains: {sorted(unknown_domains)}"
        )
    for domain, frame, prediction_column, score_mask in (
        (
            "wine",
            wine_prediction,
            "predicted_ug_l",
            ~wine_prediction["initialization_point"].to_numpy(dtype=bool),
        ),
        (
            "condensate",
            condensate_prediction,
            "predicted_captured_ug",
            np.ones(len(condensate_prediction), dtype=bool),
        ),
    ):
        if domain not in calibration_domains:
            continue
        observed = frame["observed_or_upper_bound"].to_numpy(dtype=float)
        predicted = frame[prediction_column].to_numpy(dtype=float)
        sigma = frame["sigma"].to_numpy(dtype=float)
        direct = frame["status"].eq("observed").to_numpy() & score_mask
        if direct.any():
            residuals.append((predicted[direct] - observed[direct]) / sigma[direct])
        censored = (~frame["status"].eq("observed").to_numpy()) & score_mask
        if censored.any():
            z = (observed[censored] - predicted[censored]) / sigma[censored]
            residuals.append(np.sqrt(np.maximum(-2.0 * log_ndtr(z), 0.0)))
    if include_prior:
        prior_sigma = float(config["optimization"]["weak_log_prior_sigma"])
        residuals.append((np.asarray(log_values) - default_log) / prior_sigma)
    return np.concatenate(residuals)


def _compile_observations(
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[list[tuple[str, float, str, float, float]], list[tuple[str, float, float, str, float, float]]]:
    """Compile observation rows once for repeated optimizer evaluations."""

    error = config["error_model"]
    calibration_domains = set(
        error.get("calibration_domains", ["wine", "condensate"])
    )
    unknown_domains = calibration_domains - {"wine", "condensate"}
    if unknown_domains:
        raise ValueError(
            f"Unknown aroma calibration domains: {sorted(unknown_domains)}"
        )
    if not calibration_domains:
        raise ValueError("At least one aroma calibration domain is required")
    weighting = str(error.get("calibration_weighting", "row_relative"))
    if weighting not in {"row_relative", "domain_mean"}:
        raise ValueError(f"Unknown calibration weighting {weighting!r}")
    wine_domain_sigma = np.nan
    if "wine" in calibration_domains:
        wine_observed = wine[wine["model_observation_type"].eq("observed")]
        wine_domain_sigma = max(
            float(error["wine_minimum_sigma_ug_l"]),
            float(error["wine_relative_sigma"])
            * float(wine_observed["observed_value"].astype(float).mean()),
        )
    condensate_domain_sigma = np.nan
    if "condensate" in calibration_domains:
        condensate_observed = condensate[
            condensate["mix_model_observation_type"].eq("observed")
        ]
        condensate_domain_sigma = max(
            float(error["condensate_minimum_sigma_ug"]),
            float(error["condensate_relative_sigma"])
            * float(condensate_observed["observed_value"].astype(float).mean()),
        )
    first_time = (
        wine.assign(experiment_id=wine["experiment_id"].astype(str))
        .groupby("experiment_id")["time_h"]
        .min()
        .to_dict()
    )
    wine_rows: list[tuple[str, float, str, float, float]] = []
    for row in wine.itertuples(index=False) if "wine" in calibration_domains else ():
        run = str(row.experiment_id)
        time_h = float(row.time_h)
        if np.isclose(time_h, float(first_time[run])):
            continue
        status = str(row.model_observation_type)
        reference = (
            float(row.observed_value) if status == "observed" else float(row.upper_bound)
        )
        wine_rows.append(
            (
                run,
                time_h,
                status,
                reference,
                wine_domain_sigma
                if weighting == "domain_mean"
                else _sigma(
                    reference,
                    error["wine_relative_sigma"],
                    error["wine_minimum_sigma_ug_l"],
                ),
            )
        )
    condensate_rows: list[tuple[str, float, float, str, float, float]] = []
    for row in (
        condensate.itertuples(index=False)
        if "condensate" in calibration_domains
        else ()
    ):
        run = str(row.experiment_id)
        end = float(row.time_h)
        start = end - float(row.capture_interval_h)
        status = str(row.mix_model_observation_type)
        reference = (
            float(row.observed_value) if status == "observed" else float(row.upper_bound)
        )
        condensate_rows.append(
            (
                run,
                start,
                end,
                status,
                reference,
                condensate_domain_sigma
                if weighting == "domain_mean"
                else _sigma(
                    reference,
                    error["condensate_relative_sigma"],
                    error["condensate_minimum_sigma_ug"],
                ),
            )
        )
    return wine_rows, condensate_rows


def _compiled_residual_vector(
    log_values: np.ndarray,
    forcings: dict[str, aroma.AromaForcing],
    pulse_times_h: dict[str, float],
    wine_rows: list[tuple[str, float, str, float, float]],
    condensate_rows: list[tuple[str, float, float, str, float, float]],
    model_variant: str,
    default_log: np.ndarray,
    prior_sigma: float,
) -> np.ndarray:
    simulations = {
        run: simulate_nested(
            forcing, pulse_times_h[run], log_values, model_variant
        )
        for run, forcing in forcings.items()
    }
    direct: list[float] = []
    censored_z: list[float] = []
    for run, time_h, status, reference, sigma in wine_rows:
        simulation = simulations[run]
        predicted = float(
            np.interp(time_h, simulation.time_h, simulation.liquid_ug_l)
        )
        if status == "observed":
            direct.append((predicted - reference) / sigma)
        else:
            censored_z.append((reference - predicted) / sigma)
    for run, start, end, status, reference, sigma in condensate_rows:
        simulation = simulations[run]
        predicted = float(
            np.interp(end, simulation.time_h, simulation.captured_ug)
            - np.interp(start, simulation.time_h, simulation.captured_ug)
        )
        if status == "observed":
            direct.append((predicted - reference) / sigma)
        else:
            censored_z.append((reference - predicted) / sigma)
    result = [np.asarray(direct, dtype=float)]
    if censored_z:
        result.append(
            np.sqrt(
                np.maximum(-2.0 * log_ndtr(np.asarray(censored_z, dtype=float)), 0.0)
            )
        )
    result.append((np.asarray(log_values) - default_log) / prior_sigma)
    return np.concatenate(result)


def fit_nested_species(
    species: str,
    model_variant: str,
    forcings: dict[str, aroma.AromaForcing],
    pulse_times_h: dict[str, float],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    config: dict[str, Any],
) -> NestedAromaFit:
    names = PARAMETER_NAMES_BY_MODEL[model_variant]
    defaults = np.asarray(
        [config["parameter_defaults"][species][name] for name in names], dtype=float
    )
    default_log = np.log(defaults)
    lower = np.log(
        [float(config["parameter_bounds"][name][0]) for name in names]
    )
    upper = np.log(
        [float(config["parameter_bounds"][name][1]) for name in names]
    )
    optimization = config["optimization"]
    n_starts = int(optimization["sobol_multistarts"])
    sampler = qmc.Sobol(
        d=len(names), scramble=True, seed=int(optimization["seed"])
    )
    positions = sampler.random_base2(
        int(math.ceil(math.log2(max(n_starts, 2))))
    )[:n_starts]
    starts = lower + positions * (upper - lower)
    starts[0] = np.clip(default_log, lower + 1e-8, upper - 1e-8)
    compiled_wine, compiled_condensate = _compile_observations(
        wine, condensate, config
    )
    prior_sigma = float(optimization["weak_log_prior_sigma"])
    objective = lambda values: _compiled_residual_vector(
        values,
        forcings,
        pulse_times_h,
        compiled_wine,
        compiled_condensate,
        model_variant,
        default_log,
        prior_sigma,
    )

    results = []
    rows: list[dict[str, Any]] = []
    for start_index, start in enumerate(starts):
        initial = objective(start)
        result = least_squares(
            objective,
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            loss=str(optimization["loss"]),
            f_scale=float(optimization["f_scale"]),
            max_nfev=int(optimization["max_nfev"]),
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
        )
        final = objective(result.x)
        results.append(result)
        rows.append(
            {
                "species": species,
                "model_variant": model_variant,
                "start_index": start_index,
                "success": bool(result.success),
                "nfev": int(result.nfev),
                "initial_weighted_sse": float(np.dot(initial, initial)),
                "final_weighted_sse": float(np.dot(final, final)),
                "robust_cost": float(2.0 * result.cost),
                **{
                    f"estimate__{name}": float(math.exp(value))
                    for name, value in zip(names, result.x)
                },
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["robust_cost", "final_weighted_sse"]
    ).reset_index(drop=True)
    best_index = int(summary.iloc[0]["start_index"])
    best = results[best_index]
    residual = objective(best.x)
    degrees_freedom = max(len(residual) - len(names), 1)
    residual_variance = float(np.dot(residual, residual) / degrees_freedom)
    covariance_array = (
        np.linalg.pinv(best.jac.T @ best.jac, rcond=1e-10) * residual_variance
    )
    covariance = pd.DataFrame(covariance_array, index=names, columns=names)
    active = np.isclose(best.x, lower, atol=1e-5) | np.isclose(
        best.x, upper, atol=1e-5
    )
    parameter_table = pd.DataFrame(
        {
            "species": species,
            "model_variant": model_variant,
            "parameter": names,
            "estimate": np.exp(best.x),
            "lower_bound": np.exp(lower),
            "upper_bound": np.exp(upper),
            "active_bound": active,
            "std_log_local": np.sqrt(
                np.maximum(np.diag(covariance_array), 0.0)
            ),
        }
    )
    wine_prediction, condensate_prediction, simulations = prediction_tables_nested(
        forcings,
        pulse_times_h,
        wine,
        condensate,
        best.x,
        model_variant,
        config,
    )
    max_balance = max(
        float(np.max(simulation.relative_mass_balance_error))
        for simulation in simulations.values()
    )
    validation = {
        "success": bool(best.success),
        "converged_multistarts": int(summary["success"].sum()),
        "robust_cost": float(2.0 * best.cost),
        "weighted_sse": float(np.dot(residual, residual)),
        "n_residuals_including_prior": int(len(residual)),
        "active_bound_fraction": float(np.mean(active)),
        "covariance_finite": bool(np.isfinite(covariance_array).all()),
        "covariance_condition": float(np.linalg.cond(best.jac.T @ best.jac)),
        "maximum_relative_mass_balance_error": max_balance,
        "calibration_domains": "+".join(
            config["error_model"].get(
                "calibration_domains", ["wine", "condensate"]
            )
        ),
    }
    return NestedAromaFit(
        species=species,
        model_variant=model_variant,
        parameter_names=names,
        log_values=best.x.copy(),
        parameter_table=parameter_table,
        multistart_summary=summary,
        covariance=covariance,
        wine_predictions=wine_prediction,
        condensate_predictions=condensate_prediction,
        validation=validation,
    )
