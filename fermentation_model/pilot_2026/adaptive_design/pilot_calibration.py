from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


KINETIC_STATES = ("X", "Xd", "N", "G", "F", "E", "Gly")
NUISANCE_PARAMETERS = (
    "oculyze_biomass_scale_kg_m3_per_million_cells_ml",
)


@dataclass(frozen=True)
class CalibrationTables:
    primary: pd.DataFrame
    temperature: pd.DataFrame
    co2: pd.DataFrame
    events: pd.DataFrame
    metadata: pd.DataFrame
    wine_aroma: pd.DataFrame
    condensate: pd.DataFrame


@dataclass(frozen=True)
class FitResult:
    parameter_names: tuple[str, ...]
    values: np.ndarray
    multistart_summary: pd.DataFrame
    predictions: pd.DataFrame
    residual_summary: pd.DataFrame
    error_scales: pd.DataFrame
    covariance: pd.DataFrame
    parameter_table: pd.DataFrame
    validation: dict[str, Any]


@dataclass(frozen=True)
class CO2FitResult:
    parameter_names: tuple[str, ...]
    values: np.ndarray
    multistart_summary: pd.DataFrame
    predictions: pd.DataFrame
    ess: pd.DataFrame
    parameter_table: pd.DataFrame
    validation: dict[str, Any]


def _model_module():
    # The shared campaign model is the single source of kinetic equations.  Its
    # module also contains Pyomo-DOE utilities, so the repository conda
    # environment (environment-maintainer.yml) must be active even though this
    # calibration uses SciPy's bounded trust-region least-squares solver.
    try:
        from shared import run_new_must_glycerol_estimability_doe as model
    except ModuleNotFoundError as error:
        if error.name == "pyomo":
            raise RuntimeError(
                "Pyomo is required by the shared kinetic-model module. Activate the "
                "repository idaes-pse environment before running calibration."
            ) from error
        raise
    return model


def load_tables(model_run: Path) -> CalibrationTables:
    model_run = Path(model_run)
    primary = pd.read_csv(model_run / "primary_observations.csv", parse_dates=["timestamp"])
    temperature = pd.read_csv(
        model_run / "temperature_inputs.csv", parse_dates=["timestamp"]
    )
    co2 = pd.read_csv(model_run / "co2_observations.csv", parse_dates=["timestamp"])
    events = pd.read_csv(
        model_run / "operational_events.csv",
        parse_dates=["timestamp", "timing_interval_start", "timing_interval_end"],
    )
    metadata = pd.read_csv(
        model_run / "run_metadata.csv",
        parse_dates=["sampling_start", "sampling_end", "active_end"],
    )
    wine = pd.read_csv(
        model_run / "wine_aroma_observations.csv", parse_dates=["sample_timestamp"]
    )
    condensate = pd.read_csv(
        model_run / "condensate_observations.csv",
        parse_dates=["timestamp", "capture_interval_start", "capture_interval_end"],
    )
    for table in (primary, temperature, co2, events, metadata, wine, condensate):
        table["experiment_id"] = (
            table["experiment_id"].astype(str).str.replace(r"\.0$", "", regex=True)
        )
    if "kinetic_include" not in temperature:
        raise ValueError("Model dataset predates formal kinetic-phase filtering")
    if not temperature["process_phase"].eq("active_process").all():
        raise ValueError("Cooling/postprocess temperature entered the kinetic table")
    if not primary["calibration_include"].astype(str).str.lower().isin({"true", "1"}).all():
        raise ValueError("Postprocess primary observations entered the kinetic table")
    if not events["calibration_include"].astype(str).str.lower().isin({"true", "1"}).all():
        raise ValueError("Out-of-process operational events entered the kinetic table")
    return CalibrationTables(primary, temperature, co2, events, metadata, wine, condensate)


def _pivot_primary(primary: pd.DataFrame) -> pd.DataFrame:
    index = ["experiment_id", "timestamp", "time_h"]
    wide = primary.pivot_table(
        index=index,
        columns="state",
        values="observed_value",
        aggfunc="first",
    ).reset_index()
    wide.columns.name = None
    return wide.sort_values(index).reset_index(drop=True)


def _bounds(config: dict[str, Any]) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    model = _model_module()
    names = tuple(config["primary_fit"]["parameters"]) + NUISANCE_PARAMETERS
    lower: list[float] = []
    upper: list[float] = []
    for name in config["primary_fit"]["parameters"]:
        lb, ub = model.PARAMETER_BOUNDS[name]
        lower.append(math.log(float(lb)))
        upper.append(math.log(float(ub)))
    ocu_lb, ocu_ub = config["oculyze"]["scale_bounds"]
    lower.append(math.log(float(ocu_lb)))
    upper.append(math.log(float(ocu_ub)))
    return names, np.asarray(lower), np.asarray(upper)


def _initial_vector(config: dict[str, Any], names: tuple[str, ...]) -> np.ndarray:
    model = _model_module()
    values: list[float] = []
    for name in names:
        if name in model.DEFAULT_THETA:
            values.append(math.log(float(model.DEFAULT_THETA[name])))
        elif name == "oculyze_biomass_scale_kg_m3_per_million_cells_ml":
            values.append(math.log(float(config["oculyze"]["scale_prior"])))
        else:  # pragma: no cover - guarded by the fixed name list
            raise KeyError(name)
    return np.asarray(values, dtype=float)


def _decode(
    log_values: np.ndarray, names: tuple[str, ...], config: dict[str, Any]
) -> tuple[dict[str, float], dict[str, float]]:
    model = _model_module()
    theta = dict(model.DEFAULT_THETA)
    decoded = {name: float(math.exp(value)) for name, value in zip(names, log_values)}
    for name in config["primary_fit"]["parameters"]:
        theta[name] = decoded[name]
    nuisance = {name: decoded[name] for name in NUISANCE_PARAMETERS}
    return theta, nuisance


def _first_finite(values: np.ndarray, default: float) -> float:
    finite = values[np.isfinite(values)]
    return float(finite[0]) if len(finite) else float(default)


def build_batches(
    tables: CalibrationTables,
    nuisance: dict[str, float],
    config: dict[str, Any],
) -> list[Any]:
    model = _model_module()
    wide = _pivot_primary(tables.primary)
    scale = nuisance["oculyze_biomass_scale_kg_m3_per_million_cells_ml"]
    concentration = pd.to_numeric(wide.get("cell_concentration"), errors="coerce")
    viability = pd.to_numeric(wide.get("cell_viability"), errors="coerce") / 100.0
    viability = viability.clip(0.0, 1.0)
    wide["X"] = concentration * viability * scale
    wide["Xd"] = concentration * (1.0 - viability) * scale
    pan = pd.to_numeric(wide.get("primary_amino_nitrogen"), errors="coerce")
    ammonium = pd.to_numeric(wide.get("ammonium_nitrogen"), errors="coerce")
    wide["N"] = (pan.fillna(0.0) + ammonium.fillna(0.0)) / 1000.0
    wide["N"] = wide["N"].where(pan.notna() | ammonium.notna())
    wide["G"] = pd.to_numeric(wide.get("glucose"), errors="coerce").clip(lower=0.0)
    wide["F"] = pd.to_numeric(wide.get("fructose"), errors="coerce").clip(lower=0.0)
    wide["E"] = (
        pd.to_numeric(wide.get("ethanol"), errors="coerce").clip(lower=0.0) * 7.89
    )
    wide["Gly"] = pd.to_numeric(wide.get("glycerol"), errors="coerce").clip(
        lower=0.0
    )

    batches: list[Any] = []
    for run, group in wide.groupby("experiment_id", sort=True):
        group = group.sort_values("time_h").drop_duplicates("time_h").copy()
        time = group["time_h"].to_numpy(dtype=float)
        if len(time) < 3 or np.any(np.diff(time) <= 0.0):
            raise ValueError(f"Invalid primary time grid for {run}")
        temp = tables.temperature[tables.temperature["experiment_id"].eq(run)].copy()
        temp = temp.sort_values("time_h").dropna(
            subset=["time_h", "executed_temperature_c"]
        )
        if temp.empty:
            raise ValueError(f"No active Sonda1 temperature for {run}")
        temperature = np.interp(
            time,
            temp["time_h"].to_numpy(dtype=float),
            temp["executed_temperature_c"].to_numpy(dtype=float),
        )
        event_rows = tables.events[
            tables.events["experiment_id"].eq(run)
            & tables.events["yan_added_mg_l"].notna()
            & (pd.to_numeric(tables.events["relative_time_h"], errors="coerce") > 1e-9)
        ]
        pulses = tuple(
            (
                float(row.relative_time_h),
                float(config["yan"]["fixed_historical_pulse_mg_l"]) / 1000.0,
            )
            for row in event_rows.itertuples()
        )
        observations = {
            state: group[state].to_numpy(dtype=float) for state in KINETIC_STATES
        }
        initials = {
            "X": _first_finite(observations["X"], 0.45),
            "Xd": _first_finite(observations["Xd"], 0.0),
            "N": _first_finite(observations["N"], 0.18),
            "G": _first_finite(observations["G"], 80.0),
            "F": _first_finite(observations["F"], 80.0),
            "E": _first_finite(observations["E"], 0.0),
            "Gly": _first_finite(observations["Gly"], 0.0),
        }
        batches.append(
            model.BatchData(
                medium="natural_pilot_2026",
                batch=run,
                time=time,
                temperature_c=temperature,
                pulses={
                    "N": pulses,
                    "G": tuple(),
                    "F": tuple(),
                    "E": tuple(),
                    "X": tuple(),
                },
                observations=observations,
                initials=initials,
            )
        )
    if len(batches) != 9:
        raise ValueError(f"Expected nine calibration batches; found {len(batches)}")
    return batches


def _sigma(state: str, observation: np.ndarray, config: dict[str, Any]) -> np.ndarray:
    error = config["observation_error"]
    floor = float(error["absolute_floor"][state])
    relative = float(error["relative_fraction"][state])
    return np.maximum(floor, relative * np.maximum(np.abs(observation), floor))


def _primary_residuals(
    log_values: np.ndarray,
    names: tuple[str, ...],
    tables: CalibrationTables,
    config: dict[str, Any],
    *,
    include_priors: bool = True,
    return_predictions: bool = False,
) -> np.ndarray | tuple[np.ndarray, pd.DataFrame]:
    model = _model_module()
    theta, nuisance = _decode(log_values, names, config)
    batches = build_batches(tables, nuisance, config)
    residuals: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    for batch in batches:
        simulation = model.simulate(batch, theta, batch.time)
        if simulation is None:
            bad = np.full(1000, 1e6, dtype=float)
            return (bad, pd.DataFrame()) if return_predictions else bad
        for state in KINETIC_STATES:
            observed = np.asarray(batch.observations[state], dtype=float)
            predicted = simulation.loc[batch.time, state].to_numpy(dtype=float)
            mask = np.isfinite(observed)
            if not mask.any():
                continue
            sigma = _sigma(state, observed[mask], config)
            standardized = (predicted[mask] - observed[mask]) / sigma
            residuals.append(standardized)
            if return_predictions:
                for time_h, obs, pred, sig, resid in zip(
                    batch.time[mask], observed[mask], predicted[mask], sigma, standardized
                ):
                    records.append(
                        {
                            "experiment_id": batch.batch,
                            "time_h": float(time_h),
                            "state": state,
                            "observed": float(obs),
                            "predicted": float(pred),
                            "sigma": float(sig),
                            "residual": float(pred - obs),
                            "standardized_residual": float(resid),
                        }
                    )
    if include_priors:
        x0 = _initial_vector(config, names)
        priors: list[float] = []
        kinetic_sigma = float(config["primary_fit"]["weak_log_prior_sigma"])
        for index, name in enumerate(names):
            if name in config["primary_fit"]["parameters"]:
                sigma = kinetic_sigma
            else:
                sigma = float(config["oculyze"]["scale_log_prior_sigma"])
            priors.append((float(log_values[index]) - float(x0[index])) / sigma)
        residuals.append(np.asarray(priors, dtype=float))
    vector = np.concatenate(residuals) if residuals else np.array([], dtype=float)
    if return_predictions:
        return vector, pd.DataFrame.from_records(records)
    return vector


def _covariance(
    jacobian: np.ndarray, residuals: np.ndarray, names: tuple[str, ...]
) -> pd.DataFrame:
    jacobian = np.asarray(jacobian, dtype=float)
    residuals = np.asarray(residuals, dtype=float)
    dof = max(len(residuals) - len(names), 1)
    variance = float(np.dot(residuals, residuals) / dof)
    covariance = np.linalg.pinv(jacobian.T @ jacobian, rcond=1e-10) * variance
    return pd.DataFrame(covariance, index=names, columns=names)


def fit_primary(
    tables: CalibrationTables, config: dict[str, Any]
) -> FitResult:
    names, lower, upper = _bounds(config)
    base = np.clip(_initial_vector(config, names), lower + 1e-9, upper - 1e-9)
    rng_scale = float(config["primary_fit"]["multistart_log_sigma"])
    rows: list[dict[str, Any]] = []
    results: list[Any] = []
    for index, seed in enumerate(config["random_seeds"]):
        rng = np.random.default_rng(int(seed))
        start = base if index == 0 else np.clip(
            base + rng.normal(0.0, rng_scale, size=len(base)),
            lower + 1e-7,
            upper - 1e-7,
        )
        initial = np.asarray(
            _primary_residuals(start, names, tables, config), dtype=float
        )
        result = least_squares(
            lambda values: _primary_residuals(values, names, tables, config),
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            loss=str(config["primary_fit"]["loss"]),
            f_scale=float(config["primary_fit"]["f_scale"]),
            max_nfev=int(config["primary_fit"]["max_nfev"]),
            ftol=2e-6,
            xtol=2e-6,
            gtol=2e-6,
        )
        final = np.asarray(
            _primary_residuals(result.x, names, tables, config), dtype=float
        )
        result_index = len(results)
        results.append(result)
        decoded_result = {
            f"estimate__{name}": float(math.exp(value))
            for name, value in zip(names, result.x)
        }
        rows.append(
            {
                "stage": "broad_global_screen",
                "result_index": result_index,
                "start_index": index,
                "seed": int(seed),
                "success": bool(result.success),
                "status": int(result.status),
                "nfev": int(result.nfev),
                "robust_objective": float(2.0 * result.cost),
                "initial_wsse": float(np.dot(initial, initial)),
                "final_wsse": float(np.dot(final, final)),
                "message": str(result.message),
                **decoded_result,
            }
        )
    global_summary = pd.DataFrame(rows).sort_values(
        ["robust_objective", "final_wsse"]
    ).reset_index(drop=True)
    selected_global = results[int(global_summary.iloc[0]["result_index"])]

    # A broad screen tests locality; a second, independently seeded screen
    # tests reproducibility of the selected basin.  Both are retained so the
    # local validation cannot hide alternative broad-screen minima.
    basin_sigma = float(config["primary_fit"]["selected_basin_replication_log_sigma"])
    for index, seed in enumerate(config["random_seeds"]):
        rng = np.random.default_rng(int(seed) + 1701)
        start = np.clip(
            selected_global.x + rng.normal(0.0, basin_sigma, size=len(base)),
            lower + 1e-7,
            upper - 1e-7,
        )
        initial = np.asarray(
            _primary_residuals(start, names, tables, config), dtype=float
        )
        result = least_squares(
            lambda values: _primary_residuals(values, names, tables, config),
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            loss=str(config["primary_fit"]["loss"]),
            f_scale=float(config["primary_fit"]["f_scale"]),
            max_nfev=int(
                config["primary_fit"]["selected_basin_replication_max_nfev"]
            ),
            ftol=2e-7,
            xtol=2e-7,
            gtol=2e-7,
        )
        final = np.asarray(
            _primary_residuals(result.x, names, tables, config), dtype=float
        )
        result_index = len(results)
        results.append(result)
        rows.append(
            {
                "stage": "selected_basin_replication",
                "result_index": result_index,
                "start_index": index,
                "seed": int(seed),
                "success": bool(result.success),
                "status": int(result.status),
                "nfev": int(result.nfev),
                "robust_objective": float(2.0 * result.cost),
                "initial_wsse": float(np.dot(initial, initial)),
                "final_wsse": float(np.dot(final, final)),
                "message": str(result.message),
                **{
                    f"estimate__{name}": float(math.exp(value))
                    for name, value in zip(names, result.x)
                },
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["robust_objective", "final_wsse"]
    ).reset_index(drop=True)
    best = results[int(summary.iloc[0]["result_index"])]
    residuals, predictions = _primary_residuals(
        best.x, names, tables, config, return_predictions=True
    )
    residuals = np.asarray(residuals, dtype=float)
    covariance = _covariance(best.jac, residuals, names)
    _, nuisance = _decode(best.x, names, config)
    decoded = {name: float(math.exp(value)) for name, value in zip(names, best.x)}
    active = np.isclose(best.x, lower, atol=1e-5) | np.isclose(best.x, upper, atol=1e-5)
    parameters = pd.DataFrame(
        {
            "parameter": names,
            "estimate": [decoded[name] for name in names],
            "log_estimate": best.x,
            "lower_bound": np.exp(lower),
            "upper_bound": np.exp(upper),
            "active_bound": active,
            "role": [
                "kinetic" if name in config["primary_fit"]["parameters"] else "nuisance"
                for name in names
            ],
        }
    )
    residual_summary = (
        predictions.groupby(["experiment_id", "state"], sort=True)
        .agg(
            n=("residual", "size"),
            rmse=("residual", lambda x: float(np.sqrt(np.mean(np.square(x))))),
            mean_bias=("residual", "mean"),
            weighted_rmse=(
                "standardized_residual",
                lambda x: float(np.sqrt(np.mean(np.square(x)))),
            ),
        )
        .reset_index()
    )
    error_scales = (
        predictions.groupby("state", sort=True)
        .agg(
            n=("standardized_residual", "size"),
            nominal_standardized_rmse=(
                "standardized_residual",
                lambda x: float(np.sqrt(np.mean(np.square(x)))),
            ),
            nominal_standardized_median=("standardized_residual", "median"),
        )
        .reset_index()
    )
    error_scales["model_discrepancy_multiplier"] = error_scales[
        "nominal_standardized_rmse"
    ].clip(lower=1.0)
    floors = config["observation_error"]["absolute_floor"]
    relative = config["observation_error"]["relative_fraction"]
    error_scales["nominal_absolute_floor"] = error_scales["state"].map(floors)
    error_scales["nominal_relative_fraction"] = error_scales["state"].map(relative)
    error_scales["effective_absolute_floor"] = (
        error_scales["nominal_absolute_floor"]
        * error_scales["model_discrepancy_multiplier"]
    )
    error_scales["effective_relative_fraction"] = (
        error_scales["nominal_relative_fraction"]
        * error_scales["model_discrepancy_multiplier"]
    )
    error_scales["postscale_standardized_rmse"] = (
        error_scales["nominal_standardized_rmse"]
        / error_scales["model_discrepancy_multiplier"]
    )
    scale_map = error_scales.set_index("state")[
        "model_discrepancy_multiplier"
    ].to_dict()
    residual_summary["model_discrepancy_multiplier"] = residual_summary["state"].map(
        scale_map
    )
    residual_summary["postscale_weighted_rmse"] = (
        residual_summary["weighted_rmse"]
        / residual_summary["model_discrepancy_multiplier"]
    )
    best_obj = float(summary.iloc[0]["robust_objective"])
    tolerance = float(config["validation"]["objective_relative_basin_tolerance"])
    converged = summary[summary["success"]]
    basin_replicates = converged[
        converged["stage"].eq("selected_basin_replication")
    ]
    near = basin_replicates[
        basin_replicates["robust_objective"]
        <= best_obj * (1.0 + tolerance) + 1e-12
    ]
    global_near = converged[
        converged["stage"].eq("broad_global_screen")
        & (converged["robust_objective"] <= best_obj * (1.0 + tolerance) + 1e-12)
    ]
    base_row = summary[
        summary["stage"].eq("broad_global_screen")
        & summary["start_index"].eq(0)
    ].iloc[0]
    improvement = 1.0 - float(base_row["final_wsse"]) / max(
        float(base_row["initial_wsse"]), 1e-12
    )
    validation = {
        "broad_multistarts": int(summary["stage"].eq("broad_global_screen").sum()),
        "selected_basin_replicates": int(
            summary["stage"].eq("selected_basin_replication").sum()
        ),
        "converged_multistarts": int(summary["success"].sum()),
        "multistarts_in_best_basin": len(near),
        "broad_multistarts_in_best_basin": len(global_near),
        "alternative_broad_minima_observed": len(global_near)
        < int(summary["stage"].eq("broad_global_screen").sum()),
        "objective_relative_basin_tolerance": tolerance,
        "primary_objective_improvement_fraction": float(improvement),
        "active_bound_fraction": float(np.mean(active)),
        "covariance_finite": bool(np.isfinite(covariance.to_numpy()).all()),
        "error_scale_model": (
            "nominal_measurement_sigma_times_state_specific_empirical_model_"
            "discrepancy_multiplier"
        ),
        "maximum_postscale_state_rmse": float(
            error_scales["postscale_standardized_rmse"].max()
        ),
        "fixed_historical_pulse_mg_l": float(
            config["yan"]["fixed_historical_pulse_mg_l"]
        ),
        "oculyze_scale": nuisance[
            "oculyze_biomass_scale_kg_m3_per_million_cells_ml"
        ],
    }
    return FitResult(
        names,
        np.asarray(best.x, dtype=float),
        summary,
        predictions,
        residual_summary,
        error_scales,
        covariance,
        parameters,
        validation,
    )


def _co2_base_predictions(
    tables: CalibrationTables,
    primary_fit: FitResult,
    config: dict[str, Any],
) -> pd.DataFrame:
    model = _model_module()
    theta, nuisance = _decode(primary_fit.values, primary_fit.parameter_names, config)
    batches = {
        batch.batch: batch for batch in build_batches(tables, nuisance, config)
    }
    reactor = tables.metadata.set_index("experiment_id")["reactor"].astype(str)
    volume = tables.metadata.set_index("experiment_id")["initial_volume_l"].astype(float)
    rows: list[pd.DataFrame] = []
    for run, group in tables.co2.groupby("experiment_id", sort=True):
        batch = batches[run]
        group = group.sort_values("time_h").copy()
        group = group[
            group["time_h"].between(float(batch.time.min()), float(batch.time.max()))
        ].copy()
        times = group["time_h"].to_numpy(dtype=float)
        output_time = np.unique(np.concatenate([batch.time, times]))
        simulation = model.simulate(batch, theta, output_time)
        if simulation is None:
            raise RuntimeError(f"CO2 support simulation failed for {run}")
        state_values = {
            state: np.interp(times, simulation.index.to_numpy(dtype=float), simulation[state])
            for state in model.STATE_NAMES
        }
        ethanol_rate = []
        for idx, time_h in enumerate(times):
            state = np.asarray(
                [state_values[name][idx] for name in model.STATE_NAMES], dtype=float
            )
            ethanol_rate.append(max(float(model.rhs(time_h, state, theta, batch)[5]), 0.0))
        ethanol_rate = np.asarray(ethanol_rate)
        # g ethanol/L/h -> mol ethanol/L/h -> Ln CO2/min for the full vessel,
        # using one mol CO2 per mol ethanol as the explicit stoichiometric base.
        base_flow = (
            ethanol_rate
            / 46.068
            * float(volume.loc[run])
            * 22.414
            / 60.0
        )
        out = group[
            ["experiment_id", "timestamp", "time_h", "observed_flow_ln_min"]
        ].copy()
        out["tank"] = reactor.loc[run].replace("-", "")
        out["base_stoichiometric_flow_ln_min"] = base_flow
        rows.append(out)
    return pd.concat(rows, ignore_index=True)


def _co2_parameter_spec(config: dict[str, Any]) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    names = (
        "co2_yield",
        "sensor_gain_TK31",
        "sensor_gain_TK32",
        "sensor_gain_TK33",
        "co2_background_Ln_min",
    )
    bounds = config["co2_fit"]["bounds"]
    lower = np.log(
        [
            bounds["co2_yield"][0],
            bounds["sensor_gain"][0],
            bounds["sensor_gain"][0],
            bounds["sensor_gain"][0],
            max(bounds["co2_background_Ln_min"][0], 1e-6),
        ]
    )
    upper = np.log(
        [
            bounds["co2_yield"][1],
            bounds["sensor_gain"][1],
            bounds["sensor_gain"][1],
            bounds["sensor_gain"][1],
            bounds["co2_background_Ln_min"][1],
        ]
    )
    return names, lower, upper


def _co2_predict(data: pd.DataFrame, values: np.ndarray) -> np.ndarray:
    decoded = np.exp(values)
    co2_yield = decoded[0]
    gains = {"TK31": decoded[1], "TK32": decoded[2], "TK33": decoded[3]}
    background = decoded[4]
    gain = data["tank"].map(gains).to_numpy(dtype=float)
    return (
        co2_yield
        * gain
        * data["base_stoichiometric_flow_ln_min"].to_numpy(dtype=float)
        + background
    )


def _co2_sigma(observed: np.ndarray, config: dict[str, Any]) -> np.ndarray:
    error = config["observation_error"]
    floor = float(error["absolute_floor"]["co2_Ln_min"])
    relative = float(error["relative_fraction"]["co2_Ln_min"])
    return np.maximum(floor, relative * np.maximum(np.abs(observed), floor))


def _co2_residuals(
    values: np.ndarray,
    data: pd.DataFrame,
    config: dict[str, Any],
    weights: dict[str, float] | None,
    *,
    include_priors: bool = True,
) -> np.ndarray:
    observed = data["observed_flow_ln_min"].to_numpy(dtype=float)
    standardized = (_co2_predict(data, values) - observed) / _co2_sigma(
        observed, config
    )
    if weights:
        row_weight = data["experiment_id"].map(weights).to_numpy(dtype=float)
        standardized = standardized * np.sqrt(row_weight)
    if include_priors:
        # Global yield and tank gains are multiplicatively confounded. Weak
        # zero-centred gain priors make the global yield the reference scale.
        standardized = np.concatenate([standardized, values[1:4] / 0.5])
    return standardized


def _residual_ess(
    data: pd.DataFrame, predictions: np.ndarray, config: dict[str, Any]
) -> pd.DataFrame:
    residual = (
        predictions - data["observed_flow_ln_min"].to_numpy(dtype=float)
    ) / _co2_sigma(data["observed_flow_ln_min"].to_numpy(dtype=float), config)
    work = data[["experiment_id", "timestamp"]].copy()
    work["standardized_residual"] = residual
    low, high = config["co2_fit"]["ess_clip"]
    rows = []
    for run, group in work.groupby("experiment_id", sort=True):
        values = group.sort_values("timestamp")["standardized_residual"].to_numpy(
            dtype=float
        )
        n = len(values)
        rho = float(np.corrcoef(values[:-1], values[1:])[0, 1]) if n >= 3 else 0.0
        if not np.isfinite(rho):
            rho = 0.0
        rho = float(np.clip(rho, float(low), float(high)))
        ess = float(np.clip(n * (1.0 - rho) / (1.0 + rho), 1.0, n))
        rows.append(
            {
                "experiment_id": run,
                "n_bins": n,
                "residual_lag1_autocorrelation": rho,
                "effective_sample_size": ess,
                "likelihood_weight_per_bin": ess / max(n, 1),
                "ess_basis": "final_standardized_fit_residuals_AR1",
            }
        )
    return pd.DataFrame(rows)


def _fit_co2_multistart(
    data: pd.DataFrame,
    config: dict[str, Any],
    weights: dict[str, float] | None,
    label: str,
) -> tuple[Any, pd.DataFrame]:
    names, lower, upper = _co2_parameter_spec(config)
    base = np.log([1.0, 1.0, 1.0, 1.0, 0.05])
    rows = []
    results = []
    for index, seed in enumerate(config["random_seeds"]):
        rng = np.random.default_rng(int(seed) + 71)
        start = base if index == 0 else base + rng.normal(0.0, 0.35, len(base))
        start = np.clip(start, lower + 1e-7, upper - 1e-7)
        initial = _co2_residuals(start, data, config, weights)
        result = least_squares(
            lambda values: _co2_residuals(values, data, config, weights),
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            loss="soft_l1",
            f_scale=1.0,
            max_nfev=int(config["co2_fit"]["max_nfev"]),
            ftol=1e-7,
            xtol=1e-7,
            gtol=1e-7,
        )
        final = _co2_residuals(result.x, data, config, weights)
        results.append(result)
        rows.append(
            {
                "stage": label,
                "start_index": index,
                "seed": int(seed),
                "success": bool(result.success),
                "nfev": int(result.nfev),
                "robust_objective": float(2.0 * result.cost),
                "initial_wsse": float(np.dot(initial, initial)),
                "final_wsse": float(np.dot(final, final)),
                "message": str(result.message),
            }
        )
    summary = pd.DataFrame(rows).sort_values(
        ["robust_objective", "final_wsse"]
    ).reset_index(drop=True)
    return results[int(summary.iloc[0]["start_index"])], summary


def fit_co2(
    tables: CalibrationTables,
    primary_fit: FitResult,
    config: dict[str, Any],
) -> CO2FitResult:
    data = _co2_base_predictions(tables, primary_fit, config)
    provisional, provisional_summary = _fit_co2_multistart(
        data, config, None, "provisional_unweighted"
    )
    provisional_prediction = _co2_predict(data, provisional.x)
    ess = _residual_ess(data, provisional_prediction, config)
    weights = ess.set_index("experiment_id")["likelihood_weight_per_bin"].to_dict()
    final, final_summary = _fit_co2_multistart(
        data, config, weights, "residual_ESS_weighted"
    )
    prediction = _co2_predict(data, final.x)
    # The required ESS is based on final residuals, then recorded without a
    # third circular refit. This makes the two-pass algorithm deterministic.
    final_ess = _residual_ess(data, prediction, config)
    final_weights = final_ess.set_index("experiment_id")[
        "likelihood_weight_per_bin"
    ].to_dict()
    observed = data["observed_flow_ln_min"].to_numpy(dtype=float)
    sigma = _co2_sigma(observed, config)
    data = data.copy()
    data["predicted_flow_ln_min"] = prediction
    data["sigma_ln_min"] = sigma
    data["residual_ln_min"] = prediction - observed
    data["standardized_residual"] = (prediction - observed) / sigma
    data["residual_ess_weight"] = data["experiment_id"].map(final_weights)
    names, lower, upper = _co2_parameter_spec(config)
    estimates = np.exp(final.x)
    parameter_table = pd.DataFrame(
        {
            "parameter": names,
            "estimate": estimates,
            "lower_bound": np.exp(lower),
            "upper_bound": np.exp(upper),
            "active_bound": np.isclose(final.x, lower, atol=1e-5)
            | np.isclose(final.x, upper, atol=1e-5),
            "role": ["physical", "nuisance", "nuisance", "nuisance", "nuisance"],
        }
    )
    weighted = np.sqrt(data["residual_ess_weight"].to_numpy(dtype=float))
    fitted_residual = (prediction - observed) / sigma * weighted
    zero_residual = observed / sigma * weighted
    improvement = 1.0 - float(np.dot(fitted_residual, fitted_residual)) / max(
        float(np.dot(zero_residual, zero_residual)), 1e-12
    )
    combined_summary = pd.concat(
        [provisional_summary, final_summary], ignore_index=True
    )
    best_obj = float(final_summary.iloc[0]["robust_objective"])
    tolerance = float(config["validation"]["objective_relative_basin_tolerance"])
    near = final_summary[
        final_summary["success"]
        & (final_summary["robust_objective"] <= best_obj * (1.0 + tolerance) + 1e-12)
    ]
    validation = {
        "co2_rows_fitted": len(data),
        "co2_runs_fitted": int(data["experiment_id"].nunique()),
        "residual_ess_total": float(final_ess["effective_sample_size"].sum()),
        "ess_available_for_every_run": bool(
            len(final_ess) == data["experiment_id"].nunique()
            and final_ess["effective_sample_size"].gt(0.0).all()
        ),
        "co2_objective_improvement_fraction": float(improvement),
        "converged_multistarts": int(final_summary["success"].sum()),
        "multistarts_in_best_basin": len(near),
        "weights_recomputed_from_residuals": True,
    }
    return CO2FitResult(
        names,
        np.asarray(final.x, dtype=float),
        combined_summary,
        data,
        final_ess,
        parameter_table,
        validation,
    )


def censoring_contract(tables: CalibrationTables, config: dict[str, Any]) -> dict[str, Any]:
    wine_censored = tables.wine_aroma["model_observation_type"].eq("left_censored")
    condensate_censored = tables.condensate["mix_model_observation_type"].eq(
        "left_censored"
    )
    return {
        "wine": {
            "observed": int((~wine_censored).sum()),
            "left_censored": int(wine_censored.sum()),
            "all_censored_upper_bounds_positive": bool(
                tables.wine_aroma.loc[wine_censored, "upper_bound"].gt(0.0).all()
            ),
            "operator": config["censoring"]["wine"],
        },
        "condensate": {
            "observed": int((~condensate_censored).sum()),
            "left_censored": int(condensate_censored.sum()),
            "all_censored_upper_bounds_positive": bool(
                tables.condensate.loc[condensate_censored, "upper_bound"].gt(0.0).all()
            ),
            "operator": config["censoring"]["condensate"],
        },
        "likelihood": config["censoring"]["likelihood"],
        "scales": {
            "relative_sigma": config["censoring"]["analyte_relative_sigma"],
            "minimum_sigma_ug_l": config["censoring"]["minimum_sigma_ug_l"],
        },
        "calibration_scope_note": config["censoring"]["status"],
    }


def validate_fit(
    primary: FitResult,
    co2: CO2FitResult,
    audit_resolution: dict[str, Any],
    censoring: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    limits = config["validation"]
    checks = {
        "historical_audit_resolved_without_rebaseline": audit_resolution["verdict"]
        == "PASS",
        "primary_multistarts_converged": primary.validation["converged_multistarts"]
        >= int(limits["minimum_converged_multistarts"]),
        "primary_basin_replicated": primary.validation["multistarts_in_best_basin"]
        >= int(limits["minimum_converged_multistarts"]),
        "primary_objective_improved": primary.validation[
            "primary_objective_improvement_fraction"
        ]
        >= float(limits["minimum_primary_objective_improvement_fraction"]),
        "target_fit_not_bound_dominated": primary.validation["active_bound_fraction"]
        <= float(limits["maximum_active_bound_fraction"]),
        "primary_error_and_model_discrepancy_scales_estimated": primary.validation[
            "maximum_postscale_state_rmse"
        ]
        <= 1.05,
        "covariance_available": primary.validation["covariance_finite"],
        "co2_multistarts_converged": co2.validation["converged_multistarts"]
        >= int(limits["minimum_converged_multistarts"]),
        "co2_basin_replicated": co2.validation["multistarts_in_best_basin"]
        >= int(limits["minimum_converged_multistarts"]),
        "co2_objective_improved": co2.validation[
            "co2_objective_improvement_fraction"
        ]
        >= float(limits["minimum_co2_objective_improvement_fraction"]),
        "co2_ESS_recomputed_from_final_residuals": co2.validation[
            "weights_recomputed_from_residuals"
        ]
        and co2.validation["ess_available_for_every_run"],
        "censoring_bounds_valid": censoring["wine"][
            "all_censored_upper_bounds_positive"
        ]
        and censoring["condensate"]["all_censored_upper_bounds_positive"],
    }
    computational_pass = all(checks.values())
    return {
        "computational_verdict": "PASS" if computational_pass else "FAIL",
        "release_verdict": "PASS_CONDITIONAL" if computational_pass else "FAIL",
        "checks": checks,
        "primary": primary.validation,
        "co2": co2.validation,
        "conditions": [
            "Independent Ultra audit remains required before physical design release.",
            "The Oculyze conversion is an estimated nuisance parameter, not a confirmed physical constant.",
            "Broad primary multistarts found alternative minima and the empirical state-specific discrepancy multipliers must propagate into the MBDoE posterior ensemble.",
            "Aroma censoring is fully specified, but aroma kinetic parameters are outside this reduced Stage-B fit.",
            "No temperature/nutrition profile is approved for physical execution by this calibration run.",
        ],
        "profiles_for_physical_execution": False,
    }


def load_config(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
