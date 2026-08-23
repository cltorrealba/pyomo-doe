from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import log_ndtr
from scipy.stats import qmc

from pilot_2026.adaptive_design import pilot_calibration


FORMATION_PARAMETER_NAMES = (
    "formation_growth_ug_per_g_sugar",
    "formation_stationary_ug_per_g_sugar",
)
PARAMETER_NAMES = FORMATION_PARAMETER_NAMES + (
    "effective_loss_scale",
)
DYNAMIC_TRANSFER_PARAMETER_NAMES = FORMATION_PARAMETER_NAMES + (
    "mass_transfer_kla_ref_h_inv",
    "ethanol_kla_multiplier_per_10_g_l",
)
LOSS_MODEL_EMPIRICAL = "empirical_mass_rate_scale"
LOSS_MODEL_DYNAMIC_TRANSFER = "dynamic_gas_liquid_transfer"
LOSS_MODEL_EQUILIBRIUM = "equilibrium_gas_liquid_partition"
CO2_MOLAR_MASS_KG_MOL = 0.0440095
IDEAL_GAS_CONSTANT_PA_M3_MOL_K = 8.314462618
REFERENCE_PRESSURE_PA = 101325.0
IDEAL_GAS_CONSTANT_J_MOL_K = 8.314462618
MORAKUL_REFERENCE_TEMPERATURE_K = 293.15


@dataclass(frozen=True)
class AromaForcing:
    experiment_id: str
    time_h: np.ndarray
    growth_fraction: np.ndarray
    sugar_uptake_g_l_h: np.ndarray
    loss_basis_h_inv: np.ndarray
    initial_concentration_ug_l: float
    volume_l: float
    trap_efficiency: float
    co2_rate_g_l_h: np.ndarray | None = None
    partition_basis_l_g: np.ndarray | None = None
    gas_turnover_h_inv: np.ndarray | None = None
    temperature_c: np.ndarray | None = None
    ethanol_g_l: np.ndarray | None = None
    total_sugar_g_l: np.ndarray | None = None
    loss_model: str = LOSS_MODEL_EMPIRICAL


@dataclass(frozen=True)
class AromaFit:
    species: str
    parameter_names: tuple[str, ...]
    log_values: np.ndarray
    parameter_table: pd.DataFrame
    multistart_summary: pd.DataFrame
    wine_predictions: pd.DataFrame
    condensate_predictions: pd.DataFrame
    covariance: pd.DataFrame
    profiles: pd.DataFrame
    validation: dict[str, Any]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_partition_surrogates(
    config: dict[str, Any], repository_dir: Path
) -> tuple[dict[str, dict[str, float]], dict[str, str]]:
    thermo = config["thermodynamics"]
    source = repository_dir / thermo["surrogate_source"]
    payload = load_json(source)
    models = payload.get(thermo["surrogate_key"])
    if not isinstance(models, dict):
        raise ValueError("Configured UNIFAC surrogate table is absent")
    required = set(config["priority_analytes"])
    if set(models) != required:
        raise ValueError(
            f"UNIFAC surrogate analytes differ: expected {sorted(required)}, "
            f"found {sorted(models)}"
        )
    expected_label = thermo["required_source_label"]
    for species, model in models.items():
        if model.get("source") != expected_label:
            raise ValueError(f"Unapproved partition source for {species}")
        if not 0.0 < float(model["trap_efficiency"]) <= 1.0:
            raise ValueError(f"Invalid trap efficiency for {species}")
    provenance = {"path": str(source.relative_to(repository_dir)), "sha256": sha256(source)}
    return models, provenance


def _dense_grid(batch: Any, wine: pd.DataFrame, condensate: pd.DataFrame) -> np.ndarray:
    horizon = max(
        float(batch.time.max()),
        float(pd.to_numeric(wine["time_h"], errors="coerce").max()),
        float(pd.to_numeric(condensate["time_h"], errors="coerce").max()),
    )
    regular = np.arange(0.0, horizon + 0.500001, 0.5)
    wine_time = pd.to_numeric(wine["time_h"], errors="coerce").dropna().to_numpy()
    end = pd.to_numeric(condensate["time_h"], errors="coerce").dropna().to_numpy()
    duration = pd.to_numeric(condensate["capture_interval_h"], errors="coerce")
    start = (pd.to_numeric(condensate["time_h"], errors="coerce") - duration).dropna().to_numpy()
    return np.unique(np.clip(np.concatenate([regular, wine_time, start, end]), 0.0, horizon))


def _initial_concentration(wine: pd.DataFrame) -> float:
    first = wine.sort_values("time_h").iloc[0]
    if str(first["model_observation_type"]) == "observed":
        return max(float(first["observed_value"]), 0.0)
    return max(0.5 * float(first["upper_bound"]), 0.0)


def _partition_basis(model: dict[str, Any], temp: float, ethanol: float, sugar: float) -> float:
    if model.get("model_kind") == "morakul_ethanol_temperature":
        temp_k = float(temp) + 273.15
        enthalpy_kj_mol = float(model["F3_kj_mol"]) + float(
            model["F4_kj_l_mol_g"]
        ) * float(ethanol)
        log_k = (
            float(model["F1"])
            + float(model["F2_l_g"]) * float(ethanol)
            - enthalpy_kj_mol
            # F3/F4 are reported in kJ/mol while the published expression uses
            # 1000/T; keeping R in J/(mol K) avoids applying the 1000 factor twice.
            / IDEAL_GAS_CONSTANT_J_MOL_K
            * (
                1000.0 / temp_k
                - 1000.0 / float(
                    model.get(
                        "reference_temperature_k",
                        MORAKUL_REFERENCE_TEMPERATURE_K,
                    )
                )
            )
        )
        return float(math.exp(float(np.clip(log_k, -30.0, 5.0))))
    log_k = (
        float(model["logK_ref"])
        + float(model["temp_slope"]) * (float(temp) - 20.0)
        + float(model["ethanol_slope"]) * (float(ethanol) - 50.0)
        + float(model["sugar_slope"]) * (float(sugar) - 100.0)
    )
    return float(math.exp(float(np.clip(log_k, -30.0, 5.0))))


def _co2_gas_density_g_l(temp_c: float, pressure_pa: float = REFERENCE_PRESSURE_PA) -> float:
    """Ideal-gas CO2 density in g/L (numerically equal to kg/m3)."""

    temp_k = float(temp_c) + 273.15
    return float(
        float(pressure_pa)
        * CO2_MOLAR_MASS_KG_MOL
        / (IDEAL_GAS_CONSTANT_PA_M3_MOL_K * temp_k)
    )


def build_forcings(
    tables: pilot_calibration.CalibrationTables,
    calibration_config: dict[str, Any],
    calibration_run: Path,
    species: str,
    analyte_label: str,
    partition: dict[str, float],
    co2_rate_override_by_run: dict[str, tuple[np.ndarray, np.ndarray]] | None = None,
    loss_model: str = LOSS_MODEL_EMPIRICAL,
) -> tuple[dict[str, AromaForcing], pd.DataFrame, pd.DataFrame, dict[str, float]]:
    if loss_model not in {
        LOSS_MODEL_EMPIRICAL,
        LOSS_MODEL_DYNAMIC_TRANSFER,
        LOSS_MODEL_EQUILIBRIUM,
    }:
        raise ValueError(f"Unknown aroma loss model {loss_model!r}")
    parameter_rows = pd.read_csv(calibration_run / "primary_parameter_estimates.csv")
    parameter_values = parameter_rows.set_index("parameter")["estimate"].astype(float).to_dict()
    model = pilot_calibration._model_module()
    theta = dict(model.DEFAULT_THETA)
    for name in calibration_config["primary_fit"]["parameters"]:
        theta[name] = float(parameter_values[name])
    nuisance = {
        "oculyze_biomass_scale_kg_m3_per_million_cells_ml": float(
            parameter_values["oculyze_biomass_scale_kg_m3_per_million_cells_ml"]
        )
    }
    batches = {
        batch.batch: batch
        for batch in pilot_calibration.build_batches(tables, nuisance, calibration_config)
    }
    wine = tables.wine_aroma[tables.wine_aroma["analyte"].eq(analyte_label)].copy()
    condensate = tables.condensate[tables.condensate["analyte"].eq(analyte_label)].copy()
    metadata = tables.metadata.set_index("experiment_id")
    forcings: dict[str, AromaForcing] = {}
    joint_runs = sorted(set(wine["experiment_id"]) & set(condensate["experiment_id"]))
    for run in joint_runs:
        batch = batches[str(run)]
        wine_run = wine[wine["experiment_id"].eq(run)]
        cond_run = condensate[condensate["experiment_id"].eq(run)]
        grid = _dense_grid(batch, wine_run, cond_run)
        core = model.simulate(batch, theta, grid)
        if core is None:
            raise RuntimeError(f"Primary support simulation failed for {run}")
        growth = []
        uptake = []
        loss_basis = []
        co2_rates = []
        partition_bases = []
        gas_turnovers = []
        temperatures = []
        ethanols = []
        total_sugars = []
        for time_h in grid:
            x = max(float(np.interp(time_h, core.index, core["X"])), 0.0)
            n = max(float(np.interp(time_h, core.index, core["N"])), 0.0)
            g = max(float(np.interp(time_h, core.index, core["G"])), 0.0)
            f = max(float(np.interp(time_h, core.index, core["F"])), 0.0)
            e = max(float(np.interp(time_h, core.index, core["E"])), 0.0)
            terms = model.kinetic_terms(theta, model.temperature_at(batch, float(time_h)), x, n, g, f, e)
            total = max(float(terms["sugar_total"]), 1e-8)
            maintenance = float(terms["maintenance"] * terms["maintenance_availability"])
            g_uptake = (
                theta["qXG"] * terms["growth_factor"]
                + theta["qEG"] * terms["glucose_ferm_factor"]
                + maintenance * g / total
            ) * x
            f_uptake = (
                theta["qXF"] * terms["growth_factor"]
                + theta["qEF"] * terms["fructose_ferm_factor"]
                + maintenance * f / total
            ) * x
            ethanol_rate = max(float(terms["beta_g"] + terms["beta_f"]) * x, 0.0)
            co2_rate = (44.01 / (2.0 * 46.07)) * ethanol_rate
            if co2_rate_override_by_run is not None and str(run) in co2_rate_override_by_run:
                override_time, override_rate = co2_rate_override_by_run[str(run)]
                co2_rate = max(
                    float(
                        np.interp(
                            float(time_h),
                            np.asarray(override_time, dtype=float),
                            np.asarray(override_rate, dtype=float),
                        )
                    ),
                    0.0,
                )
            temp_c = float(model.temperature_at(batch, float(time_h)))
            k_part = _partition_basis(partition, temp_c, e, g + f)
            gas_turnover = co2_rate / _co2_gas_density_g_l(temp_c)
            growth.append(n / (n + 0.035))
            uptake.append(max(float(g_uptake + f_uptake), 0.0))
            if loss_model in {
                LOSS_MODEL_DYNAMIC_TRANSFER,
                LOSS_MODEL_EQUILIBRIUM,
            }:
                loss_basis.append(max(k_part * gas_turnover, 0.0))
            else:
                loss_basis.append(max(k_part * co2_rate, 0.0))
            co2_rates.append(co2_rate)
            partition_bases.append(k_part)
            gas_turnovers.append(gas_turnover)
            temperatures.append(temp_c)
            ethanols.append(e)
            total_sugars.append(g + f)
        forcings[str(run)] = AromaForcing(
            experiment_id=str(run),
            time_h=grid,
            growth_fraction=np.asarray(growth),
            sugar_uptake_g_l_h=np.asarray(uptake),
            loss_basis_h_inv=np.asarray(loss_basis),
            initial_concentration_ug_l=_initial_concentration(wine_run),
            volume_l=float(metadata.loc[str(run), "initial_volume_l"]),
            trap_efficiency=float(partition["trap_efficiency"]),
            co2_rate_g_l_h=np.asarray(co2_rates),
            partition_basis_l_g=np.asarray(partition_bases),
            gas_turnover_h_inv=np.asarray(gas_turnovers),
            temperature_c=np.asarray(temperatures),
            ethanol_g_l=np.asarray(ethanols),
            total_sugar_g_l=np.asarray(total_sugars),
            loss_model=loss_model,
        )
    return forcings, wine, condensate, theta


def _parameter_names(forcings: dict[str, AromaForcing]) -> tuple[str, ...]:
    loss_models = {forcing.loss_model for forcing in forcings.values()}
    if len(loss_models) != 1:
        raise ValueError(f"A fit cannot mix aroma loss models: {sorted(loss_models)}")
    loss_model = next(iter(loss_models))
    if loss_model == LOSS_MODEL_DYNAMIC_TRANSFER:
        return DYNAMIC_TRANSFER_PARAMETER_NAMES
    if loss_model == LOSS_MODEL_EQUILIBRIUM:
        return FORMATION_PARAMETER_NAMES
    return PARAMETER_NAMES


def transfer_diagnostics(forcing: AromaForcing, log_values: np.ndarray) -> pd.DataFrame:
    """Return the time-varying physical loss driver used by a simulation."""

    parameter_values = np.exp(np.asarray(log_values, dtype=float))
    if forcing.loss_model == LOSS_MODEL_DYNAMIC_TRANSFER:
        if forcing.gas_turnover_h_inv is None or forcing.partition_basis_l_g is None:
            raise ValueError("Dynamic transfer forcing lacks gas-turnover or partition arrays")
        qgas = np.maximum(np.asarray(forcing.gas_turnover_h_inv, dtype=float), 0.0)
        partition = np.maximum(np.asarray(forcing.partition_basis_l_g, dtype=float), 0.0)
        ethanol = np.asarray(forcing.ethanol_g_l, dtype=float)
        transport_parameter = float(parameter_values[2])
        ethanol_multiplier = float(parameter_values[3])
        log_kla = np.log(transport_parameter) + (
            (ethanol - 50.0) / 10.0
        ) * np.log(ethanol_multiplier)
        kla = np.exp(np.clip(log_kla, -30.0, 30.0))
        transfer_efficiency = np.zeros_like(qgas)
        flowing = qgas > 1e-12
        denominator = partition[flowing] * qgas[flowing]
        transfer_efficiency[flowing] = -np.expm1(
            -kla[flowing] / np.maximum(denominator, 1e-30)
        )
        loss_coefficient = partition * qgas * transfer_efficiency
    elif forcing.loss_model == LOSS_MODEL_EQUILIBRIUM:
        if forcing.gas_turnover_h_inv is None or forcing.partition_basis_l_g is None:
            raise ValueError("Equilibrium forcing lacks gas-turnover or partition arrays")
        qgas = np.maximum(np.asarray(forcing.gas_turnover_h_inv, dtype=float), 0.0)
        partition = np.maximum(np.asarray(forcing.partition_basis_l_g, dtype=float), 0.0)
        transfer_efficiency = np.where(qgas > 1e-12, 1.0, 0.0)
        kla = np.full(len(forcing.time_h), np.nan)
        loss_coefficient = partition * qgas
    else:
        transport_parameter = float(parameter_values[2])
        qgas = (
            np.asarray(forcing.gas_turnover_h_inv, dtype=float)
            if forcing.gas_turnover_h_inv is not None
            else np.full(len(forcing.time_h), np.nan)
        )
        transfer_efficiency = np.full(len(forcing.time_h), np.nan)
        kla = np.full(len(forcing.time_h), np.nan)
        loss_coefficient = transport_parameter * np.asarray(
            forcing.loss_basis_h_inv, dtype=float
        )
    return pd.DataFrame(
        {
            "time_h": forcing.time_h,
            "co2_rate_g_l_h": forcing.co2_rate_g_l_h,
            "gas_turnover_h_inv": qgas,
            "partition_k_gas_over_liquid": forcing.partition_basis_l_g,
            "transfer_efficiency": transfer_efficiency,
            "mass_transfer_kla_h_inv": kla,
            "loss_coefficient_h_inv": loss_coefficient,
            "temperature_c": forcing.temperature_c,
            "ethanol_g_l": forcing.ethanol_g_l,
            "total_sugar_g_l": forcing.total_sugar_g_l,
        }
    )


def simulate_aroma(forcing: AromaForcing, log_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    parameter_values = np.exp(np.asarray(log_values, dtype=float))
    growth_k, stationary_k = parameter_values[:2]
    time = forcing.time_h
    loss_coefficients = transfer_diagnostics(forcing, log_values)[
        "loss_coefficient_h_inv"
    ].to_numpy(dtype=float)
    liquid = np.empty(len(time), dtype=float)
    captured = np.empty(len(time), dtype=float)
    liquid[0] = forcing.initial_concentration_ug_l
    captured[0] = 0.0
    for index in range(1, len(time)):
        dt = float(time[index] - time[index - 1])
        phi = float(forcing.growth_fraction[index - 1])
        production = (
            growth_k * phi + stationary_k * (1.0 - phi)
        ) * float(forcing.sugar_uptake_g_l_h[index - 1])
        loss_coefficient = float(loss_coefficients[index - 1])
        previous = max(float(liquid[index - 1]), 0.0)
        if loss_coefficient > 1e-12:
            decay = math.exp(-loss_coefficient * dt)
            current = previous * decay + production / loss_coefficient * (1.0 - decay)
        else:
            current = previous + production * dt
        current = max(float(current), 0.0)
        volatilized_per_l = max(production * dt - (current - previous), 0.0)
        liquid[index] = current
        captured[index] = captured[index - 1] + (
            volatilized_per_l * forcing.volume_l * forcing.trap_efficiency
        )
    return liquid, captured


def _sigma(observed_or_limit: float, relative: float, floor: float) -> float:
    return max(float(floor), float(relative) * max(abs(float(observed_or_limit)), float(floor)))


def prediction_tables(
    forcings: dict[str, AromaForcing],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    log_values: np.ndarray,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    error = config["error_model"]
    simulations = {run: simulate_aroma(forcing, log_values) for run, forcing in forcings.items()}
    wine_rows = []
    for row in wine.itertuples(index=False):
        run = str(row.experiment_id)
        forcing = forcings[run]
        liquid, _ = simulations[run]
        predicted = float(np.interp(float(row.time_h), forcing.time_h, liquid))
        status = str(row.model_observation_type)
        reference = float(row.observed_value) if status == "observed" else float(row.upper_bound)
        sigma = _sigma(reference, error["wine_relative_sigma"], error["wine_minimum_sigma_ug_l"])
        wine_rows.append(
            {
                "experiment_id": run,
                "sample_id": row.sample_id,
                "time_h": float(row.time_h),
                "status": status,
                "observed_or_upper_bound": reference,
                "predicted_ug_l": predicted,
                "sigma": sigma,
            }
        )
    cond_rows = []
    for row in condensate.itertuples(index=False):
        run = str(row.experiment_id)
        forcing = forcings[run]
        _, captured = simulations[run]
        end = float(row.time_h)
        start = end - float(row.capture_interval_h)
        predicted = float(
            np.interp(end, forcing.time_h, captured)
            - np.interp(start, forcing.time_h, captured)
        )
        status = str(row.mix_model_observation_type)
        reference = float(row.observed_value) if status == "observed" else float(row.upper_bound)
        sigma = _sigma(
            reference,
            error["condensate_relative_sigma"],
            error["condensate_minimum_sigma_ug"],
        )
        cond_rows.append(
            {
                "experiment_id": run,
                "mix_id": row.mix_id,
                "interval_start_h": start,
                "interval_end_h": end,
                "status": status,
                "observed_or_upper_bound": reference,
                "predicted_captured_ug": max(predicted, 0.0),
                "sigma": sigma,
            }
        )
    return pd.DataFrame(wine_rows), pd.DataFrame(cond_rows)


def residual_vector(
    log_values: np.ndarray,
    forcings: dict[str, AromaForcing],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    config: dict[str, Any],
    default_log: np.ndarray,
    *,
    include_prior: bool = True,
) -> np.ndarray:
    wine_pred, cond_pred = prediction_tables(forcings, wine, condensate, log_values, config)
    residuals = []
    for frame, prediction_column in (
        (wine_pred, "predicted_ug_l"),
        (cond_pred, "predicted_captured_ug"),
    ):
        observed = frame["observed_or_upper_bound"].to_numpy(dtype=float)
        predicted = frame[prediction_column].to_numpy(dtype=float)
        sigma = frame["sigma"].to_numpy(dtype=float)
        status = frame["status"].to_numpy(dtype=str)
        direct = status == "observed"
        if direct.any():
            residuals.append((predicted[direct] - observed[direct]) / sigma[direct])
        censored = ~direct
        if censored.any():
            z = (observed[censored] - predicted[censored]) / sigma[censored]
            residuals.append(np.sqrt(np.maximum(-2.0 * log_ndtr(z), 0.0)))
    if include_prior:
        prior_sigma = float(config["optimization"]["weak_log_prior_sigma"])
        residuals.append((np.asarray(log_values) - default_log) / prior_sigma)
    return np.concatenate(residuals)


def _bounds(
    config: dict[str, Any], parameter_names: tuple[str, ...]
) -> tuple[np.ndarray, np.ndarray]:
    payload = config["parameter_bounds"]
    lower = np.log([float(payload[name][0]) for name in parameter_names])
    upper = np.log([float(payload[name][1]) for name in parameter_names])
    return lower, upper


def _fit_with_fixed(
    start: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    fixed_index: int | None,
    fixed_value: float | None,
    fun,
    max_nfev: int,
):
    if fixed_index is None:
        return least_squares(fun, start, bounds=(lower, upper), max_nfev=max_nfev, x_scale="jac")
    free = [idx for idx in range(len(start)) if idx != fixed_index]
    fixed = float(np.clip(fixed_value, lower[fixed_index], upper[fixed_index]))
    start_free = start[free]
    lower_free = lower[free]
    upper_free = upper[free]

    def expand(values):
        full = np.empty(len(start), dtype=float)
        full[fixed_index] = fixed
        full[free] = values
        return full

    result = least_squares(
        lambda values: fun(expand(values)),
        start_free,
        bounds=(lower_free, upper_free),
        max_nfev=max_nfev,
        x_scale="jac",
    )
    result.full_x = expand(result.x)
    return result


def fit_species(
    species: str,
    forcings: dict[str, AromaForcing],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    config: dict[str, Any],
) -> AromaFit:
    parameter_names = _parameter_names(forcings)
    if parameter_names == DYNAMIC_TRANSFER_PARAMETER_NAMES:
        defaults = np.asarray(
            config["parameter_defaults_dynamic_transfer"][species], dtype=float
        )
    elif parameter_names == FORMATION_PARAMETER_NAMES:
        defaults = np.asarray(config["parameter_defaults"][species][:2], dtype=float)
    else:
        defaults = np.asarray(config["parameter_defaults"][species], dtype=float)
    default_log = np.log(defaults)
    lower, upper = _bounds(config, parameter_names)
    starts = int(config["optimization"]["sobol_multistarts"])
    seed = int(config["optimization"]["seed"])
    sampler = qmc.Sobol(d=len(parameter_names), scramble=True, seed=seed)
    sobol = sampler.random_base2(int(math.ceil(math.log2(max(starts, 2)))))[:starts]
    positions = lower + sobol * (upper - lower)
    positions[0] = np.clip(default_log, lower + 1e-8, upper - 1e-8)
    fun = lambda values: residual_vector(
        values, forcings, wine, condensate, config, default_log
    )
    rows = []
    results = []
    for index, start in enumerate(positions):
        initial = fun(start)
        result = least_squares(
            fun,
            start,
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            loss=str(config["optimization"]["loss"]),
            max_nfev=int(config["optimization"]["max_nfev"]),
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
        )
        final = fun(result.x)
        results.append(result)
        rows.append(
            {
                "species": species,
                "start_index": index,
                "success": bool(result.success),
                "nfev": int(result.nfev),
                "initial_objective": float(np.dot(initial, initial)),
                "final_objective": float(np.dot(final, final)),
                **{
                    f"estimate__{name}": float(math.exp(value))
                    for name, value in zip(parameter_names, result.x)
                },
            }
        )
    summary = pd.DataFrame(rows).sort_values("final_objective").reset_index(drop=True)
    best_index = int(summary.iloc[0]["start_index"])
    best = results[best_index]
    residuals = fun(best.x)
    dof = max(len(residuals) - len(parameter_names), 1)
    variance = float(np.dot(residuals, residuals) / dof)
    covariance_array = np.linalg.pinv(best.jac.T @ best.jac, rcond=1e-10) * variance
    covariance = pd.DataFrame(covariance_array, index=parameter_names, columns=parameter_names)
    estimates = np.exp(best.x)
    active = np.isclose(best.x, lower, atol=1e-5) | np.isclose(best.x, upper, atol=1e-5)
    parameter_table = pd.DataFrame(
        {
            "species": species,
            "parameter": parameter_names,
            "estimate": estimates,
            "lower_bound": np.exp(lower),
            "upper_bound": np.exp(upper),
            "active_bound": active,
            "std_log_local": np.sqrt(np.maximum(np.diag(covariance_array), 0.0)),
        }
    )
    wine_pred, cond_pred = prediction_tables(forcings, wine, condensate, best.x, config)
    wine_pred.insert(0, "species", species)
    cond_pred.insert(0, "species", species)

    profile_rows = []
    best_objective = float(np.dot(residuals, residuals))
    offsets = config["optimization"]["profile_grid_log_offsets"]
    for parameter_index, parameter in enumerate(parameter_names):
        for offset in offsets:
            fixed = float(np.clip(best.x[parameter_index] + float(offset), lower[parameter_index], upper[parameter_index]))
            result = _fit_with_fixed(
                best.x,
                lower,
                upper,
                parameter_index,
                fixed,
                fun,
                int(config["optimization"]["max_nfev"]),
            )
            full = result.full_x
            objective = float(np.dot(fun(full), fun(full)))
            profile_rows.append(
                {
                    "species": species,
                    "parameter": parameter,
                    "log_offset_requested": float(offset),
                    "fixed_value": float(math.exp(fixed)),
                    "objective": objective,
                    "delta_objective": objective - best_objective,
                    "success": bool(result.success),
                }
            )
    profiles = pd.DataFrame(profile_rows)
    basin_tolerance = 0.05
    near = summary[
        summary["success"]
        & (summary["final_objective"] <= best_objective * (1.0 + basin_tolerance) + 1e-12)
    ]
    profile_max_delta = profiles.groupby("parameter")["delta_objective"].max().to_dict()
    profile_left_delta = (
        profiles[profiles["log_offset_requested"] < 0.0]
        .groupby("parameter")["delta_objective"]
        .max()
        .to_dict()
    )
    profile_right_delta = (
        profiles[profiles["log_offset_requested"] > 0.0]
        .groupby("parameter")["delta_objective"]
        .max()
        .to_dict()
    )
    std_log = dict(zip(parameter_names, np.sqrt(np.maximum(np.diag(covariance_array), 0.0))))
    parameter_identified = {
        name: bool(
            profile_left_delta.get(name, 0.0) >= 3.84
            and profile_right_delta.get(name, 0.0) >= 3.84
            and std_log[name] <= 1.25
            and not bool(active[index])
        )
        for index, name in enumerate(parameter_names)
    }
    observation_count = len(wine) + len(condensate)
    loss_model = next(iter({forcing.loss_model for forcing in forcings.values()}))
    if loss_model == LOSS_MODEL_EQUILIBRIUM:
        loss_identified: bool | None = None
    else:
        loss_identified = all(
            parameter_identified[name] for name in parameter_names[2:]
        )
    validation = {
        "species": species,
        "converged_multistarts": int(summary["success"].sum()),
        "multistarts_in_best_basin": int(len(near)),
        "objective": best_objective,
        "objective_per_observation": best_objective / max(observation_count, 1),
        "active_bound_fraction": float(np.mean(active)),
        "covariance_finite": bool(np.isfinite(covariance_array).all()),
        "profile_max_delta_objective": {
            key: float(value) for key, value in profile_max_delta.items()
        },
        "profile_left_max_delta_objective": {
            key: float(value) for key, value in profile_left_delta.items()
        },
        "profile_right_max_delta_objective": {
            key: float(value) for key, value in profile_right_delta.items()
        },
        "parameter_identified": parameter_identified,
        "weak_parameters": [
            name for name, identified in parameter_identified.items() if not identified
        ],
        "loss_separately_identified": loss_identified,
        "loss_parameter_fixed_by_physics": loss_model == LOSS_MODEL_EQUILIBRIUM,
        "loss_model": loss_model,
    }
    return AromaFit(
        species,
        parameter_names,
        np.asarray(best.x),
        parameter_table,
        summary,
        wine_pred,
        cond_pred,
        covariance,
        profiles,
        validation,
    )


def evaluate_gate(fits: list[AromaFit]) -> dict[str, Any]:
    checks = {
        "all_species_fitted": len(fits) == 3,
        "all_multistarts_converged": all(
            fit.validation["converged_multistarts"] >= 6 for fit in fits
        ),
        "all_covariances_finite": all(fit.validation["covariance_finite"] for fit in fits),
        "no_fit_fully_bound_dominated": all(
            fit.validation["active_bound_fraction"] < 2.0 / 3.0 for fit in fits
        ),
        "no_extreme_structural_objective": all(
            fit.validation["objective_per_observation"] <= 8.0 for fit in fits
        ),
    }
    loss_status = {
        fit.species: bool(fit.validation["loss_separately_identified"]) for fit in fits
    }
    computational_pass = all(checks.values())
    all_parameter_directions_identified = all(
        all(fit.validation["parameter_identified"].values()) for fit in fits
    )
    full_release = (
        computational_pass and all(loss_status.values()) and all_parameter_directions_identified
    )
    verdict = "PASS" if full_release else "PASS_CONDITIONAL" if computational_pass else "FAIL"
    return {
        "gate": "phase_B_aroma_calibration",
        "verdict": verdict,
        "checks": checks,
        "loss_identifiability": loss_status,
        "all_parameter_directions_identified": all_parameter_directions_identified,
        "weak_parameter_directions": {
            fit.species: fit.validation["weak_parameters"] for fit in fits
        },
        "conditions": [
            "Any non-identified aroma-loss direction must remain broad in the robust MBDoE ensemble.",
            "Independent Ultra verification remains required before physical profile release.",
            "No physical profile is authorized by this gate.",
        ],
        "profiles_for_physical_execution": False,
    }
