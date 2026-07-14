from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import aroma_partition_unifac as unifac_partition
import run_new_must_glycerol_estimability_doe as base
import run_secondary_metabolite_data_review as secondary_review


RESULTS_DIR = SCRIPT_DIR / "results" / "secondary_joint_campaign_doe"
NOTEBOOK_PATH = (
    SCRIPT_DIR
    / "shared"
    / "notebooks"
    / "fermentation_secondary_joint_campaign_doe.ipynb"
)

SECONDARY_STATES = ("Pyr", "AcAld", "Acetate", "O2")
AROMA_SPECIES = ("ethyl_acetate", "isoamyl_acetate", "ethyl_octanoate")
AROMA_SHORT = {"ethyl_acetate": "EA", "isoamyl_acetate": "IAA", "ethyl_octanoate": "EO"}
AROMA_COLUMNS = {
    "ethyl_acetate": "ethyl_acetate_total",
    "isoamyl_acetate": "isoamyl_acetate_total",
    "ethyl_octanoate": "ethyl_octanoate_total",
}

FERMENTATION_TARGETS = (
    "mu0",
    "qN",
    "betaG0",
    "betaF0",
    "qEG",
    "qEF",
    "iG",
    "iE",
    "Kd0",
    "gammaG0",
    "gammaF0",
)
SECONDARY_TARGETS = (
    "kPyrS",
    "kPyrD",
    "kAldPyr",
    "kAldS",
    "kAldRed",
    "kAcAld",
    "kAcStress",
    "kAcAssim",
    "qO2",
    "kLaO2",
    "O2sat",
)
AROMA_TARGETS = (
    "k_EA_growth",
    "k_EA_stationary",
    "k_IAA_growth",
    "k_IAA_stationary",
    "k_EO_growth",
    "k_EO_stationary",
    "alpha_EA_loss",
    "alpha_IAA_loss",
    "alpha_EO_loss",
)
TARGET_PARAMETERS = FERMENTATION_TARGETS + SECONDARY_TARGETS + AROMA_TARGETS

SECONDARY_BOUNDS = {
    "kPyrS": (1e-3, 10.0),
    "kPyrD": (1e-4, 1.0),
    "kAldPyr": (1e-4, 1.0),
    "kAldS": (1e-3, 10.0),
    "kAldRed": (1e-4, 1.0),
    "kAcAld": (1e-5, 1.0),
    "kAcStress": (1e-5, 0.1),
    "kAcAssim": (1e-5, 1.0),
    "qO2": (1e-3, 20.0),
    "kLaO2": (1e-5, 0.2),
    "O2sat": (0.5, 12.0),
}
AROMA_BOUNDS = {
    "k_EA_growth": (1e-4, 2.0),
    "k_EA_stationary": (1e-4, 3.0),
    "k_IAA_growth": (1e-5, 0.2),
    "k_IAA_stationary": (1e-5, 0.3),
    "k_EO_growth": (1e-6, 0.08),
    "k_EO_stationary": (1e-6, 0.12),
    "alpha_EA_loss": (0.2, 5.0),
    "alpha_IAA_loss": (0.2, 5.0),
    "alpha_EO_loss": (0.2, 5.0),
}
PARAMETER_BOUNDS = {**{name: base.PARAMETER_BOUNDS[name] for name in base.FULL17}, **SECONDARY_BOUNDS, **AROMA_BOUNDS}

SECONDARY_DEFAULTS = {
    "kPyrS": 0.8,
    "kPyrD": 0.025,
    "kAldPyr": 0.030,
    "kAldS": 0.35,
    "kAldRed": 0.020,
    "kAcAld": 0.004,
    "kAcStress": 0.0015,
    "kAcAssim": 0.010,
    "qO2": 0.05,
    "kLaO2": 0.003,
    "O2sat": 6.5,
}
AROMA_DEFAULTS = {
    "k_EA_growth": 0.10,
    "k_EA_stationary": 0.28,
    "k_IAA_growth": 0.003,
    "k_IAA_stationary": 0.012,
    "k_EO_growth": 0.0008,
    "k_EO_stationary": 0.0030,
    "alpha_EA_loss": 1.0,
    "alpha_IAA_loss": 1.0,
    "alpha_EO_loss": 1.0,
}

SIGMA = {
    "Pyr": 6.0,
    "AcAld": 8.0,
    "Acetate": 0.035,
    "O2": 0.20,
    "CO2": 2.5,
    "ethyl_acetate_liq": 2.0,
    "isoamyl_acetate_liq": 0.08,
    "ethyl_octanoate_liq": 0.02,
    "ethyl_acetate_cond": 2.0,
    "isoamyl_acetate_cond": 0.08,
    "ethyl_octanoate_cond": 0.02,
}

O2_HALF_MG_L = 0.5
N_PHASE_HALF_KG_M3 = 0.035
CO2_G_PER_G_ETHANOL = 44.01 / (2.0 * 46.07)
PARTITION_MODE = "water_ethanol_total_sugar_as_glucose"


@dataclass(frozen=True)
class JointSimulation:
    core: pd.DataFrame
    secondary: pd.DataFrame
    aromas_liq: pd.DataFrame
    aromas_loss: pd.DataFrame
    aromas_condensate: pd.Series


def parameter_bounds(name: str) -> tuple[float, float]:
    return PARAMETER_BOUNDS[name]


def clip_all(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in PARAMETER_BOUNDS.items():
        if name in out and np.isfinite(out[name]):
            out[name] = float(np.clip(out[name], lb, ub))
    return out


def load_reference_theta() -> dict[str, float]:
    theta = dict(base.DEFAULT_THETA)
    candidates = [
        SCRIPT_DIR / "results" / "new_must_glycerol_overnight_validation" / "theta_multistart_07.csv",
        SCRIPT_DIR / "results" / "new_must_glycerol_overnight_validation" / "theta_reference.csv",
        SCRIPT_DIR / "results" / "new_must_glycerol_estimability_doe" / "theta_mixed_full17_l2.csv",
    ]
    for path in candidates:
        if path.exists():
            loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
            theta.update({name: float(value) for name, value in loaded.items() if name in theta and np.isfinite(value)})
            break
    theta.update(SECONDARY_DEFAULTS)
    theta.update(AROMA_DEFAULTS)
    return clip_all(theta)


def log_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    return np.array([math.log(float(theta[name])) for name in parameters], dtype=float)


def theta_from_log(x: np.ndarray, parameters: tuple[str, ...], base_theta: dict[str, float]) -> dict[str, float]:
    theta = dict(base_theta)
    for value, name in zip(np.asarray(x, dtype=float), parameters):
        theta[name] = float(math.exp(float(value)))
    return clip_all(theta)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    lb = []
    ub = []
    for name in parameters:
        low, high = parameter_bounds(name)
        lb.append(math.log(low))
        ub.append(math.log(high))
    return np.asarray(lb, dtype=float), np.asarray(ub, dtype=float)


def load_current_review_data() -> pd.DataFrame:
    data = secondary_review.load_all_review_data()
    for column in [
        "time_h",
        "pyruvic_acid",
        "acetaldehyde",
        "acetic_acid",
        "DO_mg_l",
        "ethyl_acetate_total",
        "isoamyl_acetate_total",
        "ethyl_octanoate_total",
    ]:
        if column in data.columns:
            data[column] = pd.to_numeric(data[column], errors="coerce")
    return data.sort_values(["medium", "batch", "time_h"]).reset_index(drop=True)


def _dense_time_grid(batch: base.BatchData | base.FutureDesign, sample_times: Iterable[float] | None = None) -> np.ndarray:
    if isinstance(batch, base.BatchData):
        horizon = float(batch.time[-1])
    else:
        horizon = float(batch.horizon_h)
    dense = set(np.arange(0.0, horizon + 1e-9, 2.0).round(8))
    if isinstance(batch, base.BatchData):
        dense.update(float(t) for t in batch.time)
    if sample_times is not None:
        dense.update(float(t) for t in sample_times)
    for schedule in batch.pulses.values():
        for time_h, _amount in schedule:
            dense.add(round(max(0.0, float(time_h) - 2.0), 8))
            dense.add(round(float(time_h), 8))
            dense.add(round(min(horizon, float(time_h) + 2.0), 8))
    return np.asarray(sorted(t for t in dense if 0.0 <= float(t) <= horizon), dtype=float)


def _interp(frame: pd.DataFrame, column: str, t: float) -> float:
    return float(np.interp(float(t), frame.index.to_numpy(dtype=float), frame[column].to_numpy(dtype=float)))


def _core_rates(theta: dict[str, float], batch: base.BatchData | base.FutureDesign, core: pd.DataFrame, t: float) -> dict[str, float]:
    x = max(_interp(core, "X", t), 0.0)
    n = max(_interp(core, "N", t), 0.0)
    g = max(_interp(core, "G", t), 0.0)
    f = max(_interp(core, "F", t), 0.0)
    e = max(_interp(core, "E", t), 0.0)
    terms = base.kinetic_terms(theta, base.temperature_at(batch, float(t)), x, n, g, f, e)
    sugar_total = max(terms["sugar_total"], 1e-8)
    growth_factor = terms["growth_factor"]
    glucose_factor = terms["glucose_ferm_factor"]
    fructose_factor = terms["fructose_ferm_factor"]
    maintenance_flux = terms["maintenance"] * terms["maintenance_availability"]
    g_uptake = (
        theta["qXG"] * growth_factor
        + theta["qEG"] * glucose_factor
        + maintenance_flux * (g / sugar_total)
    ) * x
    f_uptake = (
        theta["qXF"] * growth_factor
        + theta["qEF"] * fructose_factor
        + maintenance_flux * (f / sugar_total)
    ) * x
    ethanol_prod = (terms["beta_g"] + terms["beta_f"]) * x
    return {
        "X": x,
        "N": n,
        "G": g,
        "F": f,
        "E": e,
        "TempC": base.temperature_at(batch, float(t)),
        "sugar_uptake": max(g_uptake + f_uptake, 0.0),
        "ethanol_prod": max(ethanol_prod, 0.0),
    }


@lru_cache(maxsize=20000)
def _partition_k_cached(species: str, temp_c: float, ethanol_g_l: float, sugar_g_l: float) -> float:
    try:
        result = unifac_partition.unifac_partition_K(
            species,
            float(temp_c),
            float(ethanol_g_l),
            glucose_g_l=float(sugar_g_l),
            fructose_g_l=0.0,
            mode=PARTITION_MODE,
        )
        value = float(result["K"])
        return max(value, 1e-12)
    except Exception:
        fallback = {
            "ethyl_acetate": 0.020,
            "isoamyl_acetate": 0.005,
            "ethyl_octanoate": 0.002,
        }[species]
        return fallback * math.exp(0.035 * (float(temp_c) - 20.0) + 0.003 * (float(ethanol_g_l) - 50.0))


def partition_k(species: str, temp_c: float, ethanol_g_l: float, sugar_g_l: float) -> float:
    return _partition_k_cached(
        species,
        round(float(temp_c), 1),
        round(float(ethanol_g_l), 1),
        round(float(sugar_g_l), 1),
    )


def secondary_rhs(t: float, y: np.ndarray, theta: dict[str, float], batch: base.BatchData | base.FutureDesign, core: pd.DataFrame) -> list[float]:
    pyr, ald, acetate, o2, co2 = [max(float(v), 0.0) for v in y]
    r = _core_rates(theta, batch, core, t)
    x = r["X"]
    n_eff = r["N"]
    g_o2 = o2 / (o2 + O2_HALF_MG_L)
    h_ana = O2_HALF_MG_L / (o2 + O2_HALF_MG_L)
    h_n = n_eff / (n_eff + N_PHASE_HALF_KG_M3)
    h_e = r["E"] / (r["E"] + 60.0)
    d_pyr = theta["kPyrS"] * r["sugar_uptake"] + 0.15 * g_o2 * x - theta["kPyrD"] * pyr * x
    d_ald = (
        theta["kAldPyr"] * pyr * x
        + theta["kAldS"] * r["sugar_uptake"]
        - theta["kAldRed"] * ald * x * h_ana
        - theta["kAcAld"] * ald * x * g_o2
    )
    d_ac = theta["kAcAld"] * ald * x * g_o2 / 1000.0 + theta["kAcStress"] * x * h_e - theta["kAcAssim"] * acetate * x * h_n
    d_o2 = theta["kLaO2"] * (theta["O2sat"] - o2) - theta["qO2"] * x * o2 / (o2 + O2_HALF_MG_L)
    d_co2 = CO2_G_PER_G_ETHANOL * r["ethanol_prod"]
    return [d_pyr, d_ald, d_ac, d_o2, d_co2]


def aroma_rhs(t: float, y: np.ndarray, theta: dict[str, float], batch: base.BatchData | base.FutureDesign, core: pd.DataFrame) -> list[float]:
    liquid = {species: max(float(y[idx]), 0.0) for idx, species in enumerate(AROMA_SPECIES)}
    loss = {species: max(float(y[idx + len(AROMA_SPECIES)]), 0.0) for idx, species in enumerate(AROMA_SPECIES)}
    r = _core_rates(theta, batch, core, t)
    phi_growth = r["N"] / (r["N"] + N_PHASE_HALF_KG_M3)
    co2_rate = max(CO2_G_PER_G_ETHANOL * r["ethanol_prod"], 0.0)
    sugar = max(r["G"] + r["F"], 0.0)
    out = []
    loss_rates = []
    for species in AROMA_SPECIES:
        short = AROMA_SHORT[species]
        k_prod = theta[f"k_{short}_growth"] * phi_growth + theta[f"k_{short}_stationary"] * (1.0 - phi_growth)
        prod = k_prod * r["sugar_uptake"]
        alpha = theta[f"alpha_{short}_loss"]
        k_lg = partition_k(species, r["TempC"], r["E"], sugar)
        loss_rate = alpha * k_lg * co2_rate * liquid[species]
        out.append(prod - loss_rate)
        loss_rates.append(loss_rate)
    out.extend(loss_rates)
    return out


def integrate_secondary_euler(batch: base.BatchData | base.FutureDesign, theta: dict[str, float], core: pd.DataFrame) -> pd.DataFrame:
    time = core.index.to_numpy(dtype=float)
    y = np.array(
        [
            float(batch.initials.get("Pyr", 0.0)),
            float(batch.initials.get("AcAld", 0.0)),
            float(batch.initials.get("Acetate", 0.0)),
            float(batch.initials.get("O2", 6.5)),
            float(batch.initials.get("CO2", 0.0)),
        ],
        dtype=float,
    )
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(secondary_rhs(float(time[idx - 1]), y, theta, batch, core), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    secondary = pd.DataFrame(rows, index=time, columns=("Pyr", "AcAld", "Acetate", "O2", "CO2"))
    return secondary


def integrate_aroma_euler(batch: base.BatchData | base.FutureDesign, theta: dict[str, float], core: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    time = core.index.to_numpy(dtype=float)
    y = np.zeros(len(AROMA_SPECIES) * 2, dtype=float)
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(aroma_rhs(float(time[idx - 1]), y, theta, batch, core), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    arr = np.asarray(rows, dtype=float)
    liq = pd.DataFrame(arr[:, : len(AROMA_SPECIES)], index=time, columns=AROMA_SPECIES)
    loss = pd.DataFrame(arr[:, len(AROMA_SPECIES) :], index=time, columns=AROMA_SPECIES)
    trap_eff = {"ethyl_acetate": 0.75, "isoamyl_acetate": 0.80, "ethyl_octanoate": 0.90}
    cond = pd.Series({species: float(loss[species].iloc[-1]) * trap_eff[species] for species in AROMA_SPECIES})
    return liq, loss, cond


def simulate_joint(
    batch: base.BatchData | base.FutureDesign,
    theta: dict[str, float],
    sample_times: Iterable[float] | None = None,
    include_aroma: bool = True,
) -> JointSimulation | None:
    time = _dense_time_grid(batch, sample_times)
    core = base.simulate(batch, theta, time)
    if core is None:
        return None
    core = core.sort_index()
    secondary = integrate_secondary_euler(batch, theta, core)
    if not include_aroma:
        empty = pd.DataFrame(0.0, index=time, columns=AROMA_SPECIES)
        return JointSimulation(
            core=core,
            secondary=secondary,
            aromas_liq=empty,
            aromas_loss=empty.copy(),
            aromas_condensate=pd.Series({species: 0.0 for species in AROMA_SPECIES}),
        )
    try:
        liq, loss, cond = integrate_aroma_euler(batch, theta, core)
    except Exception:
        return None
    return JointSimulation(core=core, secondary=secondary, aromas_liq=liq, aromas_loss=loss, aromas_condensate=cond)


def precompute_core_cache(batches: list[base.BatchData], theta: dict[str, float]) -> dict[str, pd.DataFrame]:
    cache: dict[str, pd.DataFrame] = {}
    for batch in batches:
        time = _dense_time_grid(batch, batch.time)
        core = base.simulate(batch, theta, time)
        if core is not None:
            cache[batch.label] = core.sort_index()
    return cache


def simulate_secondary_from_core(
    batch: base.BatchData,
    theta: dict[str, float],
    core: pd.DataFrame,
) -> pd.DataFrame | None:
    if len(core.index) < 2:
        return None
    return integrate_secondary_euler(batch, theta, core)


def current_residual_vector(
    theta: dict[str, float],
    batches: list[base.BatchData],
    parameter_mode: str = "secondary",
    core_cache: dict[str, pd.DataFrame] | None = None,
) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in batches:
        if core_cache is not None and batch.label in core_cache and parameter_mode == "secondary":
            secondary = simulate_secondary_from_core(batch, theta, core_cache[batch.label])
            if secondary is None:
                return np.ones(1000, dtype=float) * 1e6
            aromas_liq = None
        else:
            sim = simulate_joint(batch, theta, batch.time, include_aroma=parameter_mode == "joint")
            if sim is None:
                return np.ones(1000, dtype=float) * 1e6
            secondary = sim.secondary
            aromas_liq = sim.aromas_liq
        obs_map = {
            "Pyr": batch.observations.get("Pyr", np.full_like(batch.time, np.nan, dtype=float)),
            "AcAld": batch.observations.get("AcAld", np.full_like(batch.time, np.nan, dtype=float)),
            "Acetate": batch.observations.get("Acetate", np.full_like(batch.time, np.nan, dtype=float)),
            "O2": batch.observations.get("O2", np.full_like(batch.time, np.nan, dtype=float)),
        }
        for state, obs in obs_map.items():
            obs = np.asarray(obs, dtype=float)
            mask = np.isfinite(obs)
            if mask.any():
                pred = secondary.loc[batch.time, state].to_numpy(dtype=float)
                residuals.append((pred[mask] - obs[mask]) / SIGMA[state])
        if parameter_mode == "joint":
            for species, column in AROMA_COLUMNS.items():
                obs = batch.observations.get(column, np.full_like(batch.time, np.nan, dtype=float))
                obs = np.asarray(obs, dtype=float)
                mask = np.isfinite(obs)
                if mask.any():
                    if aromas_liq is None:
                        continue
                    pred = aromas_liq.loc[batch.time, species].to_numpy(dtype=float)
                    residuals.append((pred[mask] - obs[mask]) / SIGMA[f"{species}_liq"])
    if not residuals:
        return np.array([], dtype=float)
    return np.concatenate(residuals)


def make_secondary_batches(data: pd.DataFrame) -> list[base.BatchData]:
    batches = base.make_batches(data)
    by_key = {
        (str(batch.medium), str(batch.batch)): batch
        for batch in batches
    }
    out = []
    for (medium, batch_name), group in data.groupby(["medium", "batch"], sort=True):
        key = (str(medium), str(batch_name))
        if key not in by_key:
            continue
        batch = by_key[key]
        group = group.sort_values("time_h").drop_duplicates("time_h")
        group = group.set_index("time_h").reindex(batch.time)
        obs = dict(batch.observations)
        obs["Pyr"] = pd.to_numeric(group.get("pyruvic_acid"), errors="coerce").to_numpy(dtype=float)
        obs["AcAld"] = pd.to_numeric(group.get("acetaldehyde"), errors="coerce").to_numpy(dtype=float)
        obs["Acetate"] = pd.to_numeric(group.get("acetic_acid"), errors="coerce").to_numpy(dtype=float)
        obs["O2"] = pd.to_numeric(group.get("DO_mg_l"), errors="coerce").to_numpy(dtype=float)
        for species, column in AROMA_COLUMNS.items():
            obs[column] = pd.to_numeric(group.get(column), errors="coerce").to_numpy(dtype=float)
        initials = dict(batch.initials)
        for state, column in [
            ("Pyr", "pyruvic_acid"),
            ("AcAld", "acetaldehyde"),
            ("Acetate", "acetic_acid"),
            ("O2", "DO_mg_l"),
        ]:
            series = pd.to_numeric(group.get(column), errors="coerce").dropna()
            initials[state] = float(max(series.iloc[0], 0.0)) if not series.empty else 0.0
        if not np.isfinite(initials.get("O2", np.nan)) or initials.get("O2", 0.0) <= 0.0:
            initials["O2"] = 6.5
        initials["CO2"] = 0.0
        out.append(
            base.BatchData(
                medium=batch.medium,
                batch=batch.batch,
                time=batch.time,
                temperature_c=batch.temperature_c,
                pulses=batch.pulses,
                observations=obs,
                initials=initials,
            )
        )
    return out


def secondary_initials_by_medium(data: pd.DataFrame) -> dict[str, dict[str, float]]:
    defaults = {
        "natural": {"Pyr": 0.0, "AcAld": 0.0, "Acetate": 0.08, "O2": 1.4, "CO2": 0.0},
        "synthetic": {"Pyr": 0.0, "AcAld": 0.0, "Acetate": 0.0, "O2": 6.5, "CO2": 0.0},
    }
    for medium, group in data.groupby("medium"):
        if str(medium) not in defaults:
            continue
        first_rows = []
        for _batch, batch_group in group.sort_values("time_h").groupby("batch"):
            first_rows.append(batch_group.iloc[0])
        if not first_rows:
            continue
        first = pd.DataFrame(first_rows)
        mapping = {"Pyr": "pyruvic_acid", "AcAld": "acetaldehyde", "Acetate": "acetic_acid", "O2": "DO_mg_l"}
        for state, column in mapping.items():
            if column in first.columns:
                values = pd.to_numeric(first[column], errors="coerce").dropna()
                if not values.empty:
                    defaults[str(medium)][state] = float(max(values.median(), 0.0))
    return defaults


def fit_secondary_parameters(theta0: dict[str, float], batches: list[base.BatchData], max_nfev: int) -> tuple[dict[str, float], pd.DataFrame]:
    parameters = SECONDARY_TARGETS
    core_cache = precompute_core_cache(batches, theta0)
    x0 = log_vector(theta0, parameters)
    lb, ub = log_bounds(parameters)
    x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)

    prior_sigma = {
        "qO2": 0.8,
        "kLaO2": 1.0,
        "O2sat": 0.35,
    }

    def fun(x):
        theta = theta_from_log(x, parameters, theta0)
        res = [current_residual_vector(theta, batches, parameter_mode="secondary", core_cache=core_cache)]
        for name, sigma in prior_sigma.items():
            res.append(np.array([math.log(theta[name] / theta0[name]) / sigma], dtype=float))
        return np.concatenate(res)

    start = fun(x0)
    result = least_squares(
        fun,
        x0,
        bounds=(lb, ub),
        method="trf",
        x_scale="jac",
        loss="soft_l1",
        f_scale=2.0,
        max_nfev=int(max_nfev),
        ftol=1e-5,
        xtol=1e-5,
        gtol=1e-5,
    )
    theta_hat = theta_from_log(result.x, parameters, theta0)
    end = fun(result.x)
    rows = [
        {
            "fit": "secondary_current_l2",
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "initial_wsse": float(np.dot(start, start)),
            "final_wsse": float(np.dot(end, end)),
            "n_residuals": int(len(end)),
            "wsse_per_residual": float(np.dot(end, end) / max(len(end), 1)),
        }
    ]
    return theta_hat, pd.DataFrame(rows)


def finite_difference_jacobian(
    theta: dict[str, float],
    parameters: tuple[str, ...],
    residual_fun,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    base_res = residual_fun(theta)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        low, high = parameter_bounds(name)
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), low, high))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), low, high))
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
    if not cols:
        return np.empty((len(base_res), 0)), base_res
    return np.column_stack(cols), base_res


def future_sample_times(design: base.FutureDesign) -> np.ndarray:
    return base.operational_sample_times(design.horizon_h, policy="balanced")


def o2_measurement_times(design: base.FutureDesign) -> np.ndarray:
    operational = base.operational_sample_times(design.horizon_h, policy="balanced")
    early = operational[operational <= min(float(design.horizon_h), 72.0)]
    times = np.unique(np.concatenate([np.array([0.0], dtype=float), early]))
    return times[times <= float(design.horizon_h)]


def co2_measurement_times(design: base.FutureDesign) -> np.ndarray:
    return np.arange(0.0, float(design.horizon_h) + 1e-9, 6.0, dtype=float)


def future_residual_vector(theta: dict[str, float], design: base.FutureDesign, include_core: bool = True) -> np.ndarray:
    sample_times = future_sample_times(design)
    all_times = sorted(set(sample_times).union(o2_measurement_times(design)).union(co2_measurement_times(design)).union({design.horizon_h}))
    sim = simulate_joint(design, theta, all_times)
    if sim is None:
        return np.ones(1000, dtype=float) * 1e6
    residuals: list[np.ndarray] = []
    if include_core:
        core = sim.core.loc[[t for t in sample_times if t in sim.core.index]]
        for state in base.STATE_NAMES:
            center = core[state].to_numpy(dtype=float)
            residuals.append(center / base.sigma_for_state(state, center))
    sec = sim.secondary
    for state in ("Pyr", "AcAld", "Acetate"):
        center = sec.loc[[t for t in sample_times if t in sec.index], state].to_numpy(dtype=float)
        residuals.append(center / SIGMA[state])
    o2_times = [t for t in o2_measurement_times(design) if t in sec.index]
    if o2_times:
        residuals.append(sec.loc[o2_times, "O2"].to_numpy(dtype=float) / SIGMA["O2"])
    co2_times = [t for t in co2_measurement_times(design) if t in sec.index]
    if co2_times:
        residuals.append(sec.loc[co2_times, "CO2"].to_numpy(dtype=float) / SIGMA["CO2"])
    liq = sim.aromas_liq.loc[[t for t in sample_times if t in sim.aromas_liq.index]]
    for species in AROMA_SPECIES:
        residuals.append(liq[species].to_numpy(dtype=float) / SIGMA[f"{species}_liq"])
        residuals.append(np.array([sim.aromas_condensate[species] / SIGMA[f"{species}_cond"]], dtype=float))
    return np.concatenate(residuals)


def fim_from_residuals(theta: dict[str, float], parameters: tuple[str, ...], residual_fun, step: float) -> tuple[np.ndarray, np.ndarray]:
    jac, resid = finite_difference_jacobian(theta, parameters, residual_fun, step)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T), resid


def stable_inverse(fim: np.ndarray, ridge_fraction: float = 1e-9) -> np.ndarray:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    scale = max(float(np.trace(fim)) / max(fim.shape[0], 1), 1.0)
    return np.linalg.pinv(fim + ridge_fraction * scale * np.eye(fim.shape[0]))


def fim_metrics(fim: np.ndarray, parameters: tuple[str, ...], prefix: str = "") -> dict[str, float]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eig = np.linalg.eigvalsh(fim)
    max_eig = float(np.max(eig)) if eig.size else np.nan
    floor = max(max_eig * 1e-12, np.finfo(float).tiny) if np.isfinite(max_eig) and max_eig > 0.0 else np.finfo(float).tiny
    eig_pos = np.clip(eig, floor, None)
    cov = stable_inverse(fim)
    return {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}max_eigenvalue": max_eig,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace_inv": float(np.trace(cov)),
        f"{prefix}rank_1e-8": int(np.sum(eig > max_eig * 1e-8)) if np.isfinite(max_eig) and max_eig > 0.0 else 0,
    }


def parameter_estimability(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...], label: str) -> pd.DataFrame:
    cov = stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[idx, idx], 0.0)))
        value = float(theta[name])
        low, high = parameter_bounds(name)
        active = value <= low * 1.01 or value >= high / 1.01
        if std_log <= 0.35 and not active:
            cls = "well_estimated"
        elif std_log <= 0.75 and not active:
            cls = "moderate"
        elif std_log <= 1.25 and not active:
            cls = "weak_but_actionable"
        else:
            cls = "weak_or_confounded"
        rows.append(
            {
                "analysis": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std_log,
                "approx_95_multiplier": float(math.exp(1.96 * min(std_log, 20.0))),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return pd.DataFrame(rows)


def eigen_diagnostics(fim: np.ndarray, parameters: tuple[str, ...], label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eigvals, eigvecs = np.linalg.eigh(fim)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    max_eig = float(np.max(eigvals)) if eigvals.size else np.nan
    spectrum = pd.DataFrame(
        {
            "analysis": label,
            "direction": np.arange(1, len(eigvals) + 1),
            "eigenvalue": eigvals,
            "relative_eigenvalue": eigvals / max_eig if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        }
    )
    rows = []
    for idx in range(min(8, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        top = np.argsort(np.abs(vec))[::-1][:8]
        rows.append(
            {
                "analysis": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs(vec[i]):.3f}" for i in top),
            }
        )
    return spectrum, pd.DataFrame(rows)


def prior_fim(theta: dict[str, float], batches: list[base.BatchData], step: float) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ferment_batches = [b for b in batches if b.medium in {"natural", "synthetic"}]
    jac_f, _res = base.build_jacobian(theta, FERMENTATION_TARGETS, ferment_batches, step=step)
    fim_f = jac_f.T @ jac_f
    prior = np.eye(len(TARGET_PARAMETERS), dtype=float) * 1e-6
    for i, name_i in enumerate(FERMENTATION_TARGETS):
        ii = TARGET_PARAMETERS.index(name_i)
        for j, name_j in enumerate(FERMENTATION_TARGETS):
            jj = TARGET_PARAMETERS.index(name_j)
            prior[ii, jj] = fim_f[i, j]

    core_cache = precompute_core_cache(batches, theta)
    fim_s, _ = fim_from_residuals(
        theta,
        SECONDARY_TARGETS,
        lambda th: current_residual_vector(th, batches, "secondary", core_cache=core_cache),
        step,
    )
    for i, name_i in enumerate(SECONDARY_TARGETS):
        ii = TARGET_PARAMETERS.index(name_i)
        for j, name_j in enumerate(SECONDARY_TARGETS):
            jj = TARGET_PARAMETERS.index(name_j)
            prior[ii, jj] += fim_s[i, j]

    for name in AROMA_TARGETS:
        idx = TARGET_PARAMETERS.index(name)
        cv = 0.75 if name.startswith("alpha_") else 1.25
        prior[idx, idx] += 1.0 / (cv * cv)

    prior = 0.5 * (prior + prior.T)
    spectrum, loadings = eigen_diagnostics(prior, TARGET_PARAMETERS, "current_prior")
    estim = parameter_estimability(prior, theta, TARGET_PARAMETERS, "current_prior")
    return prior, spectrum, loadings, estim


def extend_candidate_set(data: pd.DataFrame) -> list[base.FutureDesign]:
    sec_initials = secondary_initials_by_medium(data)
    designs = base.make_future_designs(data[data["medium"].isin(["natural", "synthetic"])])
    designs = [
        base.FutureDesign(
            design.name,
            design.family,
            design.medium,
            design.horizon_h,
            {**dict(design.initials), **sec_initials.get(design.medium, {})},
            design.temperature_segments,
            design.pulses,
            design.rationale,
        )
        for design in designs
    ]
    by_name = {design.name: design for design in designs}
    synthetic = by_name["synthetic_high_sugar_reference"].initials
    natural = by_name["natural_control_18C"].initials

    def p(horizon: float, **kwargs) -> dict[str, tuple[tuple[float, float], ...]]:
        schedules = {channel: [] for channel in base.INPUT_CHANNELS}
        for channel, rows in kwargs.items():
            for target, amount in rows:
                schedules[channel].append((base.nearest_operational_time(float(target), horizon), float(amount)))
        return {channel: tuple(sorted(rows)) for channel, rows in schedules.items()}

    extra = [
        base.FutureDesign(
            "synthetic_pyruvate_peak_lowN_warm_shift",
            "secondary_pyruvate",
            "synthetic",
            216.0,
            {**dict(synthetic), "N": 0.045, "G": 115.0, "F": 115.0, "E": 0.0},
            (16.0, 18.0, 24.0, 20.0),
            p(216.0, N=((50.0, 0.035),)),
            "Low initial YAN plus warm shift targets pyruvate peak and acetaldehyde drainage without adding new states.",
        ),
        base.FutureDesign(
            "synthetic_acetaldehyde_redox_proxy_ethanol",
            "secondary_acetaldehyde",
            "synthetic",
            216.0,
            {**dict(synthetic), "N": 0.11, "G": 95.0, "F": 95.0, "E": 18.0},
            (18.0, 23.0, 23.0, 19.0),
            p(216.0, N=((26.0, 0.035),), E=((98.0, 18.0),)),
            "Ethanol stress and temperature perturb acetate/acetaldehyde directions while O2 remains measured, not controlled.",
        ),
        base.FutureDesign(
            "synthetic_fast_CO2_aroma_strip_highN",
            "co2_aroma_loss",
            "synthetic",
            168.0,
            {**dict(synthetic), "N": 0.30, "G": 135.0, "F": 135.0, "E": 0.0},
            (20.0, 24.0, 24.0, 20.0),
            p(168.0),
            "High-rate fermentation excites CO2 stripping and liquid/condensate aroma split.",
        ),
        base.FutureDesign(
            "natural_aroma_matrix_cold_hot",
            "natural_aroma_matrix",
            "natural",
            216.0,
            {**dict(natural), "N": max(float(natural["N"]), 0.16), "E": 0.0},
            (14.0, 18.0, 24.0, 19.0),
            p(216.0, N=((50.0, 0.035),)),
            "Natural-matrix temperature excitation tests whether aroma/secondary parameters transfer to wine must.",
        ),
        base.FutureDesign(
            "natural_CO2_strip_reference_highN",
            "natural_co2_aroma",
            "natural",
            168.0,
            {**dict(natural), "N": max(float(natural["N"]), 0.24), "E": 0.0},
            (20.0, 23.0, 23.0, 19.0),
            p(168.0),
            "Natural high-rate reference for online CO2 and final condensate constraints.",
        ),
    ]
    designs.extend(extra)
    return designs


def candidate_fim(theta: dict[str, float], design: base.FutureDesign, parameters: tuple[str, ...], step: float) -> np.ndarray:
    fim, _ = fim_from_residuals(theta, parameters, lambda th: future_residual_vector(th, design, include_core=True), step)
    return fim


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    post_cov = stable_inverse(combined)
    out = {}
    reductions = []
    for idx, name in enumerate(parameters):
        before = float(prior_cov[idx, idx])
        after = float(post_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        out[f"var_ratio_{name}"] = ratio
        out[f"var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        if name in SECONDARY_TARGETS or name in AROMA_TARGETS:
            reductions.append(out[f"var_reduction_{name}"])
    out["new_param_mean_var_reduction"] = float(np.nanmean(reductions)) if reductions else np.nan
    out["new_param_worst_var_reduction"] = float(np.nanmin(reductions)) if reductions else np.nan
    return out


def score_fim(fim: np.ndarray, parameters: tuple[str, ...], objective: str) -> float:
    metrics = fim_metrics(fim, parameters)
    if objective == "d_opt":
        return float(metrics["logdet"])
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-18)
    return float(metrics["logdet"]) + 2.0 * math.log(min_rel) - 0.05 * math.log(max(float(metrics["trace_inv"]), 1e-18))


def greedy_select(
    fims: dict[str, np.ndarray],
    designs: dict[str, base.FutureDesign],
    prior: np.ndarray,
    parameters: tuple[str, ...],
    campaign_size: int,
    objective: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    selected: list[str] = []
    remaining = set(fims)
    current = prior.copy()
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + fims[name]
            score = score_fim(trial, parameters, objective)
            if score > best_score:
                best = name
                best_score = score
                best_metrics = {
                    **fim_metrics(trial, parameters, prefix="campaign_"),
                    **variance_reduction(prior, trial, parameters),
                }
        if best is None:
            break
        current = current + fims[best]
        remaining.remove(best)
        selected.append(best)
        design = designs[best]
        rows.append(
            {
                "objective": objective,
                "campaign_order": order,
                "candidate": best,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{v:g}" for v in design.temperature_segments),
                "score": float(best_score),
                "rationale": design.rationale,
                **(best_metrics or {}),
            }
        )
    return pd.DataFrame(rows), current


def dopt_exchange(
    seed: list[str],
    fims: dict[str, np.ndarray],
    prior: np.ndarray,
    parameters: tuple[str, ...],
    max_iter: int = 25,
) -> list[str]:
    selected = list(seed)
    pool = sorted(set(fims) - set(selected))
    current_score = score_fim(prior + sum((fims[name] for name in selected), np.zeros_like(prior)), parameters, "d_opt")
    improved = True
    iteration = 0
    while improved and iteration < int(max_iter):
        improved = False
        iteration += 1
        for out_name in list(selected):
            for in_name in list(pool):
                trial = [name for name in selected if name != out_name] + [in_name]
                trial_fim = prior + sum((fims[name] for name in trial), np.zeros_like(prior))
                trial_score = score_fim(trial_fim, parameters, "d_opt")
                if trial_score > current_score + 1e-9:
                    selected = trial
                    pool = sorted(set(fims) - set(selected))
                    current_score = trial_score
                    improved = True
                    break
            if improved:
                break
    return selected


def selection_frame(names: list[str], designs: dict[str, base.FutureDesign], prior: np.ndarray, fims: dict[str, np.ndarray], parameters: tuple[str, ...], label: str) -> tuple[pd.DataFrame, np.ndarray]:
    current = prior.copy()
    rows = []
    for order, name in enumerate(names, start=1):
        current = current + fims[name]
        design = designs[name]
        rows.append(
            {
                "objective": label,
                "campaign_order": order,
                "candidate": name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{v:g}" for v in design.temperature_segments),
                "score": score_fim(current, parameters, "d_opt"),
                "rationale": design.rationale,
                **fim_metrics(current, parameters, prefix="campaign_"),
                **variance_reduction(prior, current, parameters),
            }
        )
    return pd.DataFrame(rows), current


def plot_candidate_inputs(selected: pd.DataFrame, designs: dict[str, base.FutureDesign]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for name in selected["candidate"].drop_duplicates():
        design = designs[str(name)]
        time = np.linspace(0.0, float(design.horizon_h), 300)
        fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
        axes[0].step(time, [base.temperature_at(design, t) for t in time], where="post", color="tab:red")
        axes[0].set_ylabel("T [C]")
        for channel, ax in zip(("N", "G", "F", "E", "X"), [axes[1], axes[1], axes[1], axes[1], axes[1]]):
            schedule = design.pulses.get(channel, tuple())
            if schedule:
                ax.scatter([t for t, _ in schedule], [a for _, a in schedule], label=channel)
        axes[1].set_ylabel("Pulse amount")
        axes[1].legend(loc="best", fontsize=8)
        sample = future_sample_times(design)
        axes[2].vlines(sample, 0.0, 1.0, color="tab:blue", alpha=0.6, label="liquid samples")
        axes[2].vlines(o2_measurement_times(design), 1.05, 1.35, color="tab:green", alpha=0.8, label="O2 points")
        axes[2].vlines(co2_measurement_times(design), 1.40, 1.70, color="tab:purple", alpha=0.25, label="CO2 online grid")
        axes[2].set_ylim(0, 1.85)
        axes[2].set_yticks([])
        axes[2].set_xlabel("time [h]")
        axes[2].legend(loc="upper right", fontsize=8)
        for ax in axes:
            ax.grid(True, alpha=0.25)
        fig.suptitle(name)
        fig.tight_layout()
        fig.savefig(plot_dir / f"secondary_joint_inputs_{name}.png", dpi=170)
        plt.close(fig)


def plot_predictions(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], theta: dict[str, float]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for name in selected["candidate"].drop_duplicates():
        design = designs[str(name)]
        times = sorted(set(future_sample_times(design)).union(co2_measurement_times(design)).union(o2_measurement_times(design)).union({0.0, design.horizon_h}))
        sim = simulate_joint(design, theta, times)
        if sim is None:
            continue
        fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes[:4], ("Pyr", "AcAld", "Acetate", "O2")):
            ax.plot(sim.secondary.index, sim.secondary[state], lw=1.5)
            ax.set_title(state)
            ax.grid(True, alpha=0.25)
        axes[4].plot(sim.secondary.index, sim.secondary["CO2"], lw=1.5, color="tab:purple")
        axes[4].set_title("CO2")
        axes[4].grid(True, alpha=0.25)
        for species in AROMA_SPECIES:
            axes[5].plot(sim.aromas_liq.index, sim.aromas_liq[species], lw=1.3, label=species)
            axes[6].plot(sim.aromas_loss.index, sim.aromas_loss[species], lw=1.3, label=species)
        axes[5].set_title("Liquid aromas")
        axes[6].set_title("Volatilized aromas")
        axes[5].legend(fontsize=7)
        axes[6].legend(fontsize=7)
        axes[7].axis("off")
        axes[7].text(0.02, 0.95, "Final condensate\n" + sim.aromas_condensate.to_string(), va="top", family="monospace")
        for ax in axes[5:7]:
            ax.grid(True, alpha=0.25)
        fig.suptitle(f"Predicted extended outputs: {name}")
        fig.tight_layout()
        fig.savefig(plot_dir / f"secondary_joint_prediction_{name}.png", dpi=170)
        plt.close(fig)


def plot_current_fit(batches: list[base.BatchData], theta: dict[str, float]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        observed_states = [
            state
            for state in ("Pyr", "AcAld", "Acetate", "O2")
            if np.isfinite(np.asarray(batch.observations.get(state, []), dtype=float)).any()
        ]
        if not observed_states:
            continue
        sim = simulate_joint(batch, theta, batch.time, include_aroma=False)
        if sim is None:
            continue
        fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
        axes = axes.ravel()
        for state in ("G", "F", "E"):
            if state in sim.core.columns:
                axes[0].plot(sim.core.index, sim.core[state], lw=1.4, label=f"{state} pred")
        axes[0].set_title("Core substrates/products")
        axes[0].legend(fontsize=8)
        if "X" in sim.core.columns:
            axes[1].plot(sim.core.index, sim.core["X"], lw=1.4, color="tab:green", label="X pred")
            axes[1].set_title("Viable biomass")
            axes[1].legend(fontsize=8)
        for ax, state in zip(axes[2:6], ("Pyr", "AcAld", "Acetate", "O2")):
            ax.plot(sim.secondary.index, sim.secondary[state], lw=1.5, label=f"{state} pred")
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if mask.any():
                ax.scatter(batch.time[mask], obs[mask], s=22, color="black", label=f"{state} obs")
            ax.set_title(state)
            ax.legend(fontsize=8)
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.set_xlabel("time [h]")
        fig.suptitle(f"Current fit: {batch.medium}/{batch.batch}")
        fig.tight_layout()
        fig.savefig(plot_dir / f"secondary_current_fit_{batch.medium}_{batch.batch}.png", dpi=170)
        plt.close(fig)


def run_pyomo_checks(theta: dict[str, float], selected: pd.DataFrame, designs: dict[str, base.FutureDesign], args) -> pd.DataFrame:
    rows = []
    old_results_dir = base.RESULTS_DIR
    base.RESULTS_DIR = RESULTS_DIR
    try:
        for name in selected["candidate"].head(int(args.pyomo_candidates)):
            design = designs[str(name)]
            try:
                row = base.run_pyomo_doe_check(theta, design, base.REDUCED11, "balanced", args.step)
                rows.append({"candidate": name, "status": "ok", **row})
            except Exception as err:
                rows.append({"candidate": name, "status": "failed", "error": f"{type(err).__name__}: {err}"})
    finally:
        base.RESULTS_DIR = old_results_dir
    return pd.DataFrame(rows)


def write_notebook() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Secondary-state joint campaign DOE\n\n"
            "This notebook summarizes the extended fermentation DOE pipeline with pyruvate, acetaldehyde, acetate, dissolved oxygen exposure, CO2, liquid aroma measurements, and final condensate constraints."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "ROOT = Path.cwd()\n"
            "if not (ROOT / 'results').exists() and (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'results/secondary_joint_campaign_doe'\n"
            "fit = pd.read_csv(RESULTS / 'secondary_fit_summary.csv')\n"
            "estim = pd.read_csv(RESULTS / 'current_prior_estimability.csv')\n"
            "ranking = pd.read_csv(RESULTS / 'candidate_ranking.csv')\n"
            "selected = pd.read_csv(RESULTS / 'selected_campaigns.csv')\n"
            "display(fit)\n"
            "display(estim.sort_values(['analysis','classification','std_log_approx']).head(40))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Candidate ranking\n\n"
            "The ranking combines current-data prior information with future outputs. Liquid samples are limited to operational hours; CO2 is treated as online; O2 is measured on the operational early-time grid and is not optimized as a controlled input."
        ),
        nbformat.v4.new_code_cell(
            "display(ranking.sort_values('hybrid_score', ascending=False).head(20))\n"
            "display(selected)"
        ),
        nbformat.v4.new_markdown_cell("## Eigenvalue diagnostics"),
        nbformat.v4.new_code_cell(
            "spectrum = pd.read_csv(RESULTS / 'campaign_eigen_spectrum.csv')\n"
            "loadings = pd.read_csv(RESULTS / 'campaign_weak_loadings.csv')\n"
            "display(spectrum.groupby('analysis').head(6))\n"
            "display(loadings)"
        ),
        nbformat.v4.new_markdown_cell("## Pyomo.DoE checks"),
        nbformat.v4.new_code_cell(
            "pyomo_path = RESULTS / 'pyomo_reduced11_checks.csv'\n"
            "display(pd.read_csv(pyomo_path) if pyomo_path.exists() else 'not run')"
        ),
        nbformat.v4.new_markdown_cell("## Protocol plots"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_joint_inputs_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Current-fit plots"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_current_fit_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Predicted extended outputs"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_joint_prediction_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Report"),
        nbformat.v4.new_code_cell("print((RESULTS / 'secondary_joint_campaign_report.md').read_text(encoding='utf-8'))"),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def write_report(
    fit_summary: pd.DataFrame,
    current_estimability: pd.DataFrame,
    ranking: pd.DataFrame,
    selected: pd.DataFrame,
    dopt_selected: pd.DataFrame,
    pyomo_checks: pd.DataFrame,
) -> None:
    weak_current = current_estimability[current_estimability["classification"].eq("weak_or_confounded")]
    lines = [
        "# Secondary-state joint campaign DOE",
        "",
        "## Scope",
        "",
        "This run extends the current glycerol model with pyruvate, acetaldehyde, acetate, dissolved oxygen exposure, CO2, and liquid/condensate aroma outputs.",
        "",
        "Oxygen is not treated as a manipulated design input. Existing DO measurements are used as direct O2 observations, and future DO points are scheduled on the operational early-time grid because the lab cannot deliver a reliable oxygen pulse without extra experimental work.",
        "",
        "Final condensate is used as an integral volatilization constraint. Liquid aroma trajectories remain net liquid outputs, while alpha-loss parameters are regularized by literature/UNIFAC structure.",
        "",
        "## Secondary fit",
        "",
        fit_summary.to_markdown(index=False),
        "",
        "## Current weak/confounded parameters",
        "",
        weak_current.to_markdown(index=False) if not weak_current.empty else "No weak/confounded parameters by the approximate current FIM criterion.",
        "",
        "## Top candidate ranking",
        "",
        ranking.sort_values("hybrid_score", ascending=False).head(12).to_markdown(index=False),
        "",
        "## Selected hybrid campaign",
        "",
        selected.to_markdown(index=False),
        "",
        "## D-opt exchange benchmark",
        "",
        dopt_selected.to_markdown(index=False),
        "",
    ]
    reduced_path = RESULTS_DIR / "selected_campaigns_reduced_fix_secondary_degenerate_do.csv"
    reduced_label = "Reduced benchmark with secondary degenerate terms fixed"
    reduced_text = "The current DO-aware fit leaves `kAldPyr`, `kAldRed`, and `kAcAssim` weak or boundary-dominated. Keeping them free makes the weakest FIM directions dominated by parameters that the available controls cannot excite directly."
    if not reduced_path.exists():
        reduced_path = RESULTS_DIR / "selected_campaigns_reduced_fix_redox_assim_degenerate.csv"
        reduced_label = "Reduced benchmark with redox/assimilation terms fixed"
        reduced_text = "The current fit pushed `kAldRed` and `kAcAssim` to their lower bounds. Keeping them free makes the weakest FIM directions dominated by parameters that the available controls cannot excite directly."
    if not reduced_path.exists():
        reduced_path = RESULTS_DIR / "selected_campaigns_reduced_fix_acetate_degenerate.csv"
        reduced_label = "Reduced benchmark with degenerate acetate terms fixed"
        reduced_text = "The current fit pushed several acetate/acetaldehyde terms to their lower bounds. Keeping them free makes the weakest FIM directions dominated by parameters that the available controls cannot excite directly."
    if reduced_path.exists():
        reduced = pd.read_csv(reduced_path)
        reduced_hybrid = reduced[reduced["objective"].astype(str).str.contains("hybrid") & ~reduced["objective"].astype(str).str.contains("exchange")].copy()
        reduced_exchange = reduced[reduced["objective"].astype(str).str.contains("exchange")].copy()
        lines.extend(
            [
                f"## {reduced_label}",
                "",
                reduced_text,
                "",
                "Reduced hybrid campaign:",
                "",
                reduced_hybrid.to_markdown(index=False),
                "",
                "Reduced D-opt exchange campaign:",
                "",
                reduced_exchange.to_markdown(index=False),
                "",
            ]
        )
    lines.extend(
        [
        "## Pyomo.DoE reduced checks",
        "",
        pyomo_checks.to_markdown(index=False) if not pyomo_checks.empty else "_Not run._",
        "",
        "## Interpretation",
        "",
        "- If hybrid and D-opt choose similar campaigns, the design is not an artifact of the hybrid score.",
        "- Synthetic experiments dominate kinetic separation because macronutrient composition is controllable.",
        "- Natural experiments remain necessary to test matrix transfer, aroma partition behavior, and model validity in wine-like must.",
        "- Additional fermentations are useful only if they excite a new weak eigen-direction; otherwise replicate effort should go to repeatability and analytical confidence.",
        ]
    )
    (RESULTS_DIR / "secondary_joint_campaign_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Extended secondary-state and aroma DOE campaign.")
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--fit-nfev", type=int, default=80)
    parser.add_argument("--step", type=float, default=1e-2)
    parser.add_argument("--pyomo-candidates", type=int, default=3)
    parser.add_argument("--skip-pyomo", action="store_true")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    theta0 = load_reference_theta()
    data = load_current_review_data()
    data.to_csv(RESULTS_DIR / "secondary_joint_input_data.csv", index=False)
    batches = make_secondary_batches(data)

    theta_hat, fit_summary = fit_secondary_parameters(theta0, batches, max_nfev=args.fit_nfev)
    pd.Series(theta_hat, name="theta_secondary_joint").to_csv(RESULTS_DIR / "theta_secondary_joint.csv")
    fit_summary.to_csv(RESULTS_DIR / "secondary_fit_summary.csv", index=False)
    plot_current_fit(batches, theta_hat)

    prior, prior_spectrum, prior_loadings, prior_estimability = prior_fim(theta_hat, batches, args.step)
    pd.DataFrame(prior, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "current_prior_fim.csv")
    prior_spectrum.to_csv(RESULTS_DIR / "current_prior_eigen_spectrum.csv", index=False)
    prior_loadings.to_csv(RESULTS_DIR / "current_prior_weak_loadings.csv", index=False)
    prior_estimability.to_csv(RESULTS_DIR / "current_prior_estimability.csv", index=False)

    designs = extend_candidate_set(data)
    design_by_name = {design.name: design for design in designs}
    rows = []
    fims: dict[str, np.ndarray] = {}
    for design in designs:
        print(f"[candidate] {design.name}", flush=True)
        try:
            fim = candidate_fim(theta_hat, design, TARGET_PARAMETERS, args.step)
            fims[design.name] = fim
            pd.DataFrame(fim, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / f"candidate_fim_{design.name}.csv")
            combined = prior + fim
            metrics = {
                "candidate": design.name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "hybrid_score": score_fim(combined, TARGET_PARAMETERS, "hybrid"),
                "dopt_score": score_fim(combined, TARGET_PARAMETERS, "d_opt"),
                **fim_metrics(fim, TARGET_PARAMETERS, prefix="new_"),
                **fim_metrics(combined, TARGET_PARAMETERS, prefix="combined_"),
                **variance_reduction(prior, combined, TARGET_PARAMETERS),
            }
            rows.append(metrics)
        except Exception as err:
            rows.append({"candidate": design.name, "family": design.family, "medium": design.medium, "status": "failed", "error": f"{type(err).__name__}: {err}"})
    ranking = pd.DataFrame(rows)
    ranking.to_csv(RESULTS_DIR / "candidate_ranking.csv", index=False)

    selected_hybrid, fim_hybrid = greedy_select(fims, design_by_name, prior, TARGET_PARAMETERS, args.campaign_size, "hybrid")
    selected_dgreedy, fim_dgreedy = greedy_select(fims, design_by_name, prior, TARGET_PARAMETERS, args.campaign_size, "d_opt")
    exchanged_names = dopt_exchange(selected_hybrid["candidate"].astype(str).tolist(), fims, prior, TARGET_PARAMETERS)
    selected_exchange, fim_exchange = selection_frame(exchanged_names, design_by_name, prior, fims, TARGET_PARAMETERS, "dopt_exchange_from_hybrid")

    selected_all = pd.concat([selected_hybrid, selected_dgreedy, selected_exchange], ignore_index=True)
    selected_all.to_csv(RESULTS_DIR / "selected_campaigns.csv", index=False)
    pd.DataFrame(fim_hybrid, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "campaign_fim_hybrid.csv")
    pd.DataFrame(fim_dgreedy, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "campaign_fim_dopt_greedy.csv")
    pd.DataFrame(fim_exchange, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "campaign_fim_dopt_exchange.csv")

    spectra = []
    loadings = []
    estim_rows = []
    for label, fim in [
        ("hybrid", fim_hybrid),
        ("dopt_greedy", fim_dgreedy),
        ("dopt_exchange_from_hybrid", fim_exchange),
    ]:
        spectrum, weak = eigen_diagnostics(fim, TARGET_PARAMETERS, label)
        spectra.append(spectrum)
        loadings.append(weak)
        estim_rows.append(parameter_estimability(fim, theta_hat, TARGET_PARAMETERS, label))
    pd.concat(spectra, ignore_index=True).to_csv(RESULTS_DIR / "campaign_eigen_spectrum.csv", index=False)
    pd.concat(loadings, ignore_index=True).to_csv(RESULTS_DIR / "campaign_weak_loadings.csv", index=False)
    pd.concat(estim_rows, ignore_index=True).to_csv(RESULTS_DIR / "campaign_estimability.csv", index=False)

    plot_candidate_inputs(selected_hybrid, design_by_name)
    plot_predictions(selected_hybrid.head(min(6, len(selected_hybrid))), design_by_name, theta_hat)

    pyomo_checks = pd.DataFrame()
    if not args.skip_pyomo:
        pyomo_checks = run_pyomo_checks(theta_hat, selected_hybrid, design_by_name, args)
        pyomo_checks.to_csv(RESULTS_DIR / "pyomo_reduced11_checks.csv", index=False)

    metadata = {
        "campaign_size": int(args.campaign_size),
        "target_parameters": list(TARGET_PARAMETERS),
        "fermentation_targets": list(FERMENTATION_TARGETS),
        "secondary_targets": list(SECONDARY_TARGETS),
        "aroma_targets": list(AROMA_TARGETS),
        "oxygen_policy": "existing DO is fitted as O2; future O2 is measured on an operational early-time grid; not a manipulated DOE input",
        "sampling_policy": "balanced operational liquid samples; DO on early operational points; online CO2 every 6 h in FIM; final condensate",
    }
    (RESULTS_DIR / "secondary_joint_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_report(fit_summary, prior_estimability, ranking, selected_hybrid, selected_exchange, pyomo_checks)
    write_notebook()

    print("\nSelected hybrid campaign:")
    print(selected_hybrid[["campaign_order", "candidate", "family", "medium", "campaign_logdet", "new_param_worst_var_reduction"]].to_string(index=False))
    print(f"\n[done] results={RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
