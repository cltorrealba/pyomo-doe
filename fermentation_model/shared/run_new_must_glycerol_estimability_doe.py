from __future__ import annotations

import argparse
import math
import os
import sys
import textwrap
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.stats import chi2

try:  # Notebook export is optional for model/calibration reuse.
    import nbformat
except ModuleNotFoundError:  # pragma: no cover - depends on runtime extras
    nbformat = None

try:  # Pyomo-DOE is needed only by the design entry points below.
    import pyomo.environ as pyo
    from pyomo.contrib.doe import DesignOfExperiments
except ModuleNotFoundError:  # pragma: no cover - calibration can use SciPy only
    pyo = None
    DesignOfExperiments = None

SHARED_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SHARED_DIR.parent
SCRIPT_DIR = FERMENTATION_MODEL_DIR
REPO_ROOT = FERMENTATION_MODEL_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from shared.new_must_data_loader import load_new_must_data
from shared.paths import SHARED_RESULTS_DIR

RESULTS_DIR = SHARED_RESULTS_DIR / "new_must_glycerol_estimability_doe"
NORMALIZED_DATA_PATH = SHARED_RESULTS_DIR / "new_must_data_loading" / "new_must_normalized_long.csv"
OLD_FINAL_THETA_PATH = SHARED_RESULTS_DIR / "identifiability_reduction" / "theta_final_identifiable.csv"
NOTEBOOK_PATH = (
    SHARED_DIR
    / "notebooks"
    / "fermentation_new_must_glycerol_estimability_doe.ipynb"
)

FULL15 = (
    "mu0",
    "sN",
    "qN",
    "qXG",
    "qXF",
    "betaG0",
    "sG",
    "betaF0",
    "sF",
    "qEG",
    "qEF",
    "iG",
    "iE",
    "Kd0",
    "m0",
)
GLYCEROL_PARAMETERS = ("gammaG0", "gammaF0")
FULL17 = FULL15 + GLYCEROL_PARAMETERS
REDUCED11 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF", "iG", "iE", "Kd0", "gammaG0", "gammaF0")
CORE_FIT = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF", "iG")
WEAK_FROM_OLD = ("sN", "qXG", "qXF", "sG", "sF", "iE", "Kd0", "m0")
STATE_NAMES = ("X", "Xd", "N", "G", "F", "E", "Gly")
INPUT_CHANNELS = ("N", "G", "F", "E", "X")

PHYSICAL_PARAMETER_BOUNDS = {
    "mu0": (5e-2, 1e0),
    "betaG0": (1e-1, 1e1),
    "betaF0": (1e-1, 1e1),
    "Kn0": (1e-2, 1e-1),
    "Kg0": (1e-2, 1e1),
    "Kf0": (1e-2, 1e1),
    "Kig0": (1e1, 1e3),
    "Kie0": (1e1, 1e3),
    "Kd0": (1e-4, 5e-2),
    "Yxn": (1e0, 1e2),
    "Yxg": (1e-1, 1e1),
    "Yxf": (1e-1, 1e1),
    "Yeg": (1e-1, 1e1),
    "Yef": (1e-1, 1e1),
    "m0": (1e-4, 5e-2),
}

PARAMETER_BOUNDS = {
    "mu0": PHYSICAL_PARAMETER_BOUNDS["mu0"],
    "sN": (
        PHYSICAL_PARAMETER_BOUNDS["mu0"][0] / PHYSICAL_PARAMETER_BOUNDS["Kn0"][1],
        PHYSICAL_PARAMETER_BOUNDS["mu0"][1] / PHYSICAL_PARAMETER_BOUNDS["Kn0"][0],
    ),
    "qN": (
        PHYSICAL_PARAMETER_BOUNDS["mu0"][0] / PHYSICAL_PARAMETER_BOUNDS["Yxn"][1],
        PHYSICAL_PARAMETER_BOUNDS["mu0"][1] / PHYSICAL_PARAMETER_BOUNDS["Yxn"][0],
    ),
    "qXG": (
        PHYSICAL_PARAMETER_BOUNDS["mu0"][0] / PHYSICAL_PARAMETER_BOUNDS["Yxg"][1],
        PHYSICAL_PARAMETER_BOUNDS["mu0"][1] / PHYSICAL_PARAMETER_BOUNDS["Yxg"][0],
    ),
    "qXF": (
        PHYSICAL_PARAMETER_BOUNDS["mu0"][0] / PHYSICAL_PARAMETER_BOUNDS["Yxf"][1],
        PHYSICAL_PARAMETER_BOUNDS["mu0"][1] / PHYSICAL_PARAMETER_BOUNDS["Yxf"][0],
    ),
    "betaG0": PHYSICAL_PARAMETER_BOUNDS["betaG0"],
    "sG": (
        PHYSICAL_PARAMETER_BOUNDS["betaG0"][0] / PHYSICAL_PARAMETER_BOUNDS["Kg0"][1],
        PHYSICAL_PARAMETER_BOUNDS["betaG0"][1] / PHYSICAL_PARAMETER_BOUNDS["Kg0"][0],
    ),
    "betaF0": PHYSICAL_PARAMETER_BOUNDS["betaF0"],
    "sF": (
        PHYSICAL_PARAMETER_BOUNDS["betaF0"][0] / PHYSICAL_PARAMETER_BOUNDS["Kf0"][1],
        PHYSICAL_PARAMETER_BOUNDS["betaF0"][1] / PHYSICAL_PARAMETER_BOUNDS["Kf0"][0],
    ),
    "qEG": (
        PHYSICAL_PARAMETER_BOUNDS["betaG0"][0] / PHYSICAL_PARAMETER_BOUNDS["Yeg"][1],
        PHYSICAL_PARAMETER_BOUNDS["betaG0"][1] / PHYSICAL_PARAMETER_BOUNDS["Yeg"][0],
    ),
    "qEF": (
        PHYSICAL_PARAMETER_BOUNDS["betaF0"][0] / PHYSICAL_PARAMETER_BOUNDS["Yef"][1],
        PHYSICAL_PARAMETER_BOUNDS["betaF0"][1] / PHYSICAL_PARAMETER_BOUNDS["Yef"][0],
    ),
    "iG": (1.0 / PHYSICAL_PARAMETER_BOUNDS["Kig0"][1], 1.0 / PHYSICAL_PARAMETER_BOUNDS["Kig0"][0]),
    "iE": (1.0 / PHYSICAL_PARAMETER_BOUNDS["Kie0"][1], 1.0 / PHYSICAL_PARAMETER_BOUNDS["Kie0"][0]),
    "Kd0": PHYSICAL_PARAMETER_BOUNDS["Kd0"],
    "m0": PHYSICAL_PARAMETER_BOUNDS["m0"],
    "gammaG0": (1e-4, 2.0),
    "gammaF0": (1e-4, 2.0),
}

DEFAULT_THETA = {
    "mu0": 0.18,
    "sN": 18.0,
    "qN": 0.18 / 19.69,
    "qXG": 0.18 / 1.60,
    "qXF": 0.18 / 1.60,
    "betaG0": 0.225,
    "sG": 0.225 / 7.5,
    "betaF0": 0.225,
    "sF": 0.225 / 7.5,
    "qEG": 0.225 / 0.49,
    "qEF": 0.225 / 0.49,
    "iG": 0.02,
    "iE": 1.0 / 40.0,
    "Kd0": 0.00044,
    "m0": 0.01,
    "gammaG0": 0.024,
    "gammaF0": 0.024,
}

FIXED_CONSTANTS = {
    "Cde": 0.0415,
    "Etd": 130000.0,
    "R": 8.314,
    "Eac": 59453.0,
    "Eafe": 11000.0,
    "EaKn": 46055.0,
    "EaKg": 46055.0,
    "EaKf": 46055.0,
    "EaKig": 46055.0,
    "EaKie": 46055.0,
    "Eam": 37681.0,
}

MEASUREMENT_ERROR_FLOOR = {
    "X": 0.06,
    "Xd": 0.06,
    "N": 0.012,
    "G": 2.5,
    "F": 2.5,
    "E": 2.0,
    "Gly": 0.35,
}
MEASUREMENT_ERROR_REL = {
    "X": 0.08,
    "Xd": 0.12,
    "N": 0.08,
    "G": 0.025,
    "F": 0.025,
    "E": 0.025,
    "Gly": 0.05,
}
STATE_COLUMNS = {
    "X": "X_viable_kg_m3",
    "Xd": "X_dead_kg_m3",
    "N": "N_kg_m3",
    "G": "G_g_l",
    "F": "F_g_l",
    "E": "E_g_l",
    "Gly": "glycerol_g_l",
}
STATE_BOUNDS = {
    "X": (0.0, 20.0),
    "Xd": (0.0, 30.0),
    "N": (0.0, 5.0),
    "G": (0.0, 350.0),
    "F": (0.0, 350.0),
    "E": (0.0, 220.0),
    "Gly": (0.0, 30.0),
}
MAINTENANCE_SUGAR_CUTOFF_KG_M3 = 1.0


@dataclass(frozen=True)
class BatchData:
    medium: str
    batch: str
    time: np.ndarray
    temperature_c: np.ndarray
    pulses: dict[str, tuple[tuple[float, float], ...]]
    observations: dict[str, np.ndarray]
    initials: dict[str, float]

    @property
    def label(self) -> str:
        return f"{self.medium}/{self.batch}"


@dataclass(frozen=True)
class FutureDesign:
    name: str
    family: str
    medium: str
    horizon_h: float
    initials: dict[str, float]
    temperature_segments: tuple[float, ...]
    pulses: dict[str, tuple[tuple[float, float], ...]]
    rationale: str


def load_normalized_data() -> pd.DataFrame:
    if NORMALIZED_DATA_PATH.exists():
        data = pd.read_csv(NORMALIZED_DATA_PATH)
    else:
        data = load_new_must_data()
        NORMALIZED_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(NORMALIZED_DATA_PATH, index=False)
    data = data.copy()
    for col in STATE_COLUMNS.values():
        if col in data.columns:
            data[col] = pd.to_numeric(data[col], errors="coerce")
    for col in ["time_h", "temperature_c", "N_pulse_kg_m3"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    for col in ["X_viable_kg_m3", "X_dead_kg_m3", "N_kg_m3", "G_g_l", "F_g_l", "E_g_l", "glycerol_g_l"]:
        if col in data.columns:
            data.loc[data[col] < 0.0, col] = 0.0
    return data.sort_values(["medium", "batch", "time_h"]).reset_index(drop=True)


def _first_finite(values: Iterable[float], default: float) -> float:
    arr = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
    if arr.empty:
        return float(default)
    return float(arr.iloc[0])


def _initials_from_group(group: pd.DataFrame) -> dict[str, float]:
    g0 = _first_finite(group["G_g_l"], 80.0)
    f0 = _first_finite(group["F_g_l"], 80.0)
    return {
        "X": _first_finite(group["X_viable_kg_m3"], 0.45),
        "Xd": _first_finite(group["X_dead_kg_m3"], 0.0),
        "N": _first_finite(group["N_kg_m3"], 0.18),
        "G": g0,
        "F": f0,
        "E": _first_finite(group["E_g_l"], 0.0),
        "Gly": _first_finite(group["glycerol_g_l"], 0.0),
    }


def make_batches(data: pd.DataFrame) -> list[BatchData]:
    batches: list[BatchData] = []
    for (medium, batch), group in data.groupby(["medium", "batch"], sort=True):
        group = group.sort_values("time_h").drop_duplicates("time_h").reset_index(drop=True)
        time = group["time_h"].to_numpy(dtype=float)
        if len(time) < 3 or not np.all(np.diff(time) > 0.0):
            continue
        temperature = group["temperature_c"].interpolate(limit_direction="both").fillna(18.0).to_numpy(dtype=float)
        pulse_rows = group[["time_h", "N_pulse_kg_m3"]].dropna()
        pulses = {
            "N": tuple((float(row.time_h), float(row.N_pulse_kg_m3)) for row in pulse_rows.itertuples() if float(row.N_pulse_kg_m3) > 0.0),
            "G": tuple(),
            "F": tuple(),
            "E": tuple(),
            "X": tuple(),
        }
        observations = {}
        for state, col in STATE_COLUMNS.items():
            values = group[col].to_numpy(dtype=float)
            if state in STATE_BOUNDS:
                lb, _ub = STATE_BOUNDS[state]
                values = np.where(np.isfinite(values), np.maximum(values, lb), values)
            observations[state] = values
        batches.append(
            BatchData(
                medium=str(medium),
                batch=str(batch),
                time=time,
                temperature_c=temperature,
                pulses=pulses,
                observations=observations,
                initials=_initials_from_group(group),
            )
        )
    return batches


def complete_initial_theta(batches: list[BatchData]) -> dict[str, float]:
    theta = dict(DEFAULT_THETA)
    if batches:
        initial_sugars = [b.initials["G"] + b.initials["F"] for b in batches if b.initials["G"] + b.initials["F"] > 0.0]
        if initial_sugars:
            theta["iG"] = 4.0 / float(np.nanmedian(initial_sugars))
    if OLD_FINAL_THETA_PATH.exists():
        previous = pd.read_csv(OLD_FINAL_THETA_PATH, index_col=0).iloc[:, 0].to_dict()
        theta.update({name: float(value) for name, value in previous.items() if name in theta and np.isfinite(value)})
    return clip_theta(theta)


def clip_theta(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in PARAMETER_BOUNDS.items():
        value = float(out.get(name, DEFAULT_THETA[name]))
        out[name] = float(np.clip(value, float(lb), float(ub)))
    return out


def theta_to_log_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    return np.array([math.log(float(theta[name])) for name in parameters], dtype=float)


def log_vector_to_theta(x: np.ndarray, parameters: tuple[str, ...], base_theta: dict[str, float]) -> dict[str, float]:
    theta = dict(base_theta)
    for value, name in zip(np.asarray(x, dtype=float), parameters):
        theta[name] = float(math.exp(float(value)))
    return clip_theta(theta)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    lower = []
    upper = []
    for name in parameters:
        lb, ub = PARAMETER_BOUNDS[name]
        lower.append(math.log(float(lb)))
        upper.append(math.log(float(ub)))
    return np.array(lower, dtype=float), np.array(upper, dtype=float)


def temperature_at(batch: BatchData | FutureDesign, t: float) -> float:
    if isinstance(batch, BatchData):
        return float(np.interp(float(t), batch.time, batch.temperature_c))
    edges = np.linspace(0.0, float(batch.horizon_h), len(batch.temperature_segments) + 1)
    idx = int(np.searchsorted(edges[1:-1], float(t), side="right"))
    return float(batch.temperature_segments[min(max(idx, 0), len(batch.temperature_segments) - 1)])


def pulse_rate(t: float, schedule: tuple[tuple[float, float], ...], width_h: float = 1.5) -> float:
    """Return a causal diagnostic pulse-rate approximation.

    The production simulator below does not smear pulses: it applies exact state
    jumps at their event times.  This one-sided exponential remains only for
    legacy diagnostics that require a rate-valued function.  In particular it
    is identically zero before every event, unlike the former symmetric
    Gaussian approximation.
    """

    width_h = max(float(width_h), 1e-6)
    total = 0.0
    for pulse_time, amount in schedule:
        elapsed = float(t) - float(pulse_time)
        amount = float(amount)
        if amount <= 0.0 or elapsed < 0.0:
            continue
        total += (amount / width_h) * math.exp(-elapsed / width_h)
    return float(total)


def kinetic_terms(theta: dict[str, float], temperature_c: float, x: float, n: float, g: float, f: float, e: float) -> dict[str, float]:
    x = max(float(x), 0.0)
    n = max(float(n), 0.0)
    g = max(float(g), 0.0)
    f = max(float(f), 0.0)
    e = max(float(e), 0.0)
    temp_k = float(temperature_c) + 273.15
    constants = FIXED_CONSTANTS
    eps = 1e-8
    r = constants["R"]
    a_mu = math.exp(constants["Eac"] * (temp_k - 300.0) / (300.0 * r * temp_k))
    a_beta = math.exp(constants["Eafe"] * (temp_k - 296.15) / (296.15 * r * temp_k))
    a_kn = math.exp(constants["EaKn"] * (temp_k - 293.15) / (293.15 * r * temp_k))
    a_kg = math.exp(constants["EaKg"] * (temp_k - 293.15) / (293.15 * r * temp_k))
    a_kf = math.exp(constants["EaKf"] * (temp_k - 293.15) / (293.15 * r * temp_k))
    a_kig = math.exp(constants["EaKig"] * (temp_k - 293.15) / (293.15 * r * temp_k))
    a_kie = math.exp(constants["EaKie"] * (temp_k - 293.15) / (293.15 * r * temp_k))
    kn = (theta["mu0"] / theta["sN"]) * a_kn
    kg = (theta["betaG0"] / theta["sG"]) * a_kg
    kf = (theta["betaF0"] / theta["sF"]) * a_kf
    i_g = theta["iG"] / (a_kig + eps)
    i_e = theta["iE"] / (a_kie + eps)
    n_lim = n / (n + kn + eps)
    g_lim = g / (g + kg + eps)
    f_lim = f / (f + kf + eps)
    e_inhib = 1.0 / (1.0 + i_e * e)
    g_inhib_for_f = 1.0 / (1.0 + i_g * g)
    growth_factor = a_mu * n_lim
    glucose_ferm_factor = a_beta * g_lim * e_inhib
    fructose_ferm_factor = a_beta * f_lim * g_inhib_for_f * e_inhib
    mu = theta["mu0"] * growth_factor
    beta_g = theta["betaG0"] * glucose_ferm_factor
    beta_f = theta["betaF0"] * fructose_ferm_factor
    maintenance = theta["m0"] * math.exp(constants["Eam"] * (temp_k - 293.3) / (293.3 * r * temp_k))
    td = -0.0001 * e**3 + 0.0049 * e**2 - 0.1279 * e + 315.89
    if temp_k >= td:
        kd = theta["Kd0"] * math.exp(
            (constants["Cde"] * e) + (constants["Etd"] * (temp_k - 305.65)) / (305.65 * r * temp_k)
        )
    else:
        kd = 0.0
    sugar_total = g + f + eps
    maintenance_availability = sugar_total / (sugar_total + MAINTENANCE_SUGAR_CUTOFF_KG_M3)
    return {
        "x": x,
        "growth_factor": growth_factor,
        "glucose_ferm_factor": glucose_ferm_factor,
        "fructose_ferm_factor": fructose_ferm_factor,
        "mu": mu,
        "kd": kd,
        "beta_g": beta_g,
        "beta_f": beta_f,
        "maintenance": maintenance,
        "maintenance_availability": maintenance_availability,
        "sugar_total": sugar_total,
    }


def rhs(t: float, y: np.ndarray, theta: dict[str, float], batch: BatchData | FutureDesign) -> list[float]:
    x, xd, n, g, f, e, gly = [max(0.0, float(value)) for value in y]
    terms = kinetic_terms(theta, temperature_at(batch, t), x, n, g, f, e)
    sugar_total = terms["sugar_total"]
    growth_factor = terms["growth_factor"]
    glucose_ferm_factor = terms["glucose_ferm_factor"]
    fructose_ferm_factor = terms["fructose_ferm_factor"]
    maintenance_flux = terms["maintenance"] * terms["maintenance_availability"]
    dx = (terms["mu"] - terms["kd"]) * x
    dxd = terms["kd"] * x
    dn = -theta["qN"] * growth_factor * x
    dg = -(
        theta["qXG"] * growth_factor
        + theta["qEG"] * glucose_ferm_factor
        + maintenance_flux * (g / sugar_total)
    ) * x
    df = -(
        theta["qXF"] * growth_factor
        + theta["qEF"] * fructose_ferm_factor
        + maintenance_flux * (f / sugar_total)
    ) * x
    de = (terms["beta_g"] + terms["beta_f"]) * x
    dgly = (theta["gammaG0"] * glucose_ferm_factor + theta["gammaF0"] * fructose_ferm_factor) * x
    return [dx, dxd, dn, dg, df, de, dgly]


def initial_vector(batch: BatchData | FutureDesign) -> np.ndarray:
    return np.array([float(batch.initials[state]) for state in STATE_NAMES], dtype=float)


def _pulse_events(
    batch: BatchData | FutureDesign,
    t0: float,
    tf: float,
) -> dict[float, np.ndarray]:
    events: dict[float, np.ndarray] = {}
    state_index = {state: index for index, state in enumerate(STATE_NAMES)}
    for channel in INPUT_CHANNELS:
        for event_time, amount in batch.pulses.get(channel, tuple()):
            event_time = float(event_time)
            amount = float(amount)
            if not np.isfinite(event_time) or not np.isfinite(amount) or amount < 0.0:
                raise ValueError("Pulse times and amounts must be finite; amounts must be nonnegative")
            if amount == 0.0 or event_time < t0 - 1e-12 or event_time > tf + 1e-12:
                continue
            delta = events.setdefault(event_time, np.zeros(len(STATE_NAMES), dtype=float))
            delta[state_index[channel]] += amount
    return events


def simulate(
    batch: BatchData | FutureDesign,
    theta: dict[str, float],
    output_time: np.ndarray | None = None,
    *,
    sample_event_order: str = "sample_before_action",
    integration_max_step_h: float = 2.0,
) -> pd.DataFrame | None:
    """Integrate continuous dynamics between events and apply exact pulse jumps.

    ``sample_event_order`` controls the state returned when a requested sample
    shares a timestamp with an action. Pilot 2026 approves
    ``sample_before_action``: the recorded state is pre-jump and the complete
    dose is applied immediately afterward. The alternate convention remains
    available for explicit numerical tests.
    """

    if sample_event_order not in {"sample_before_action", "action_before_sample"}:
        raise ValueError("sample_event_order must be sample_before_action or action_before_sample")
    if not np.isfinite(integration_max_step_h) or integration_max_step_h <= 0.0:
        raise ValueError("integration_max_step_h must be finite and positive")
    if output_time is None:
        if isinstance(batch, BatchData):
            output_time = batch.time
        else:
            output_time = operational_sample_times(batch.horizon_h)
    output_time = np.asarray(sorted(set(float(t) for t in output_time)), dtype=float)
    if len(output_time) == 0:
        return None
    t0 = float(output_time[0])
    tf = float(output_time[-1])
    if isinstance(batch, BatchData):
        y0 = initial_vector(batch)
        if abs(float(batch.time[0]) - t0) > 1e-9:
            t0 = float(batch.time[0])
            output_time = np.asarray(sorted(set([t0] + [float(t) for t in output_time])), dtype=float)
    else:
        y0 = initial_vector(batch)
        if t0 > 0.0:
            output_time = np.asarray(sorted(set([0.0] + [float(t) for t in output_time])), dtype=float)
            t0 = 0.0
    if tf < t0:
        return None
    try:
        events = _pulse_events(batch, t0, tf)
    except ValueError:
        raise
    output_set = set(float(value) for value in output_time)
    timeline = sorted(output_set | set(events))
    current_time = float(t0)
    current_state = np.asarray(y0, dtype=float).copy()
    records: dict[float, np.ndarray] = {}
    try:
        for time_h in timeline:
            if time_h < current_time - 1e-12:
                continue
            if time_h > current_time + 1e-12:
                sol = solve_ivp(
                    lambda t, y: rhs(t, y, theta, batch),
                    (current_time, time_h),
                    current_state,
                    t_eval=[time_h],
                    method="LSODA",
                    rtol=1e-6,
                    atol=1e-8,
                    max_step=float(integration_max_step_h),
                )
                if not sol.success or sol.y.shape[1] != 1:
                    return None
                current_state = np.asarray(sol.y[:, -1], dtype=float)
                current_time = float(time_h)
            delta = events.get(float(time_h))
            if delta is not None and sample_event_order == "sample_before_action":
                if float(time_h) in output_set:
                    records[float(time_h)] = current_state.copy()
                current_state = current_state + delta
            elif delta is not None:
                current_state = current_state + delta
                if float(time_h) in output_set:
                    records[float(time_h)] = current_state.copy()
            elif float(time_h) in output_set:
                records[float(time_h)] = current_state.copy()
    except Exception:
        return None
    if any(float(time_h) not in records for time_h in output_time):
        return None
    values = np.vstack([records[float(time_h)] for time_h in output_time])
    for idx, state in enumerate(STATE_NAMES):
        lb, ub = STATE_BOUNDS[state]
        values[:, idx] = np.clip(values[:, idx], lb, ub)
    return pd.DataFrame(values, index=output_time, columns=STATE_NAMES)


def sigma_for_state(state: str, observed_or_simulated: np.ndarray | float) -> np.ndarray:
    arr = np.asarray(observed_or_simulated, dtype=float)
    floor = float(MEASUREMENT_ERROR_FLOOR[state])
    rel = float(MEASUREMENT_ERROR_REL[state])
    return np.maximum(floor, rel * np.maximum(np.abs(arr), floor))


def residual_vector(
    theta: dict[str, float],
    batches: list[BatchData],
    states: tuple[str, ...] = STATE_NAMES,
    l2_reference: dict[str, float] | None = None,
    l2_parameters: tuple[str, ...] = tuple(),
    l2_lambda: float = 0.0,
) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in batches:
        sim = simulate(batch, theta, batch.time)
        if sim is None:
            return np.ones(1000, dtype=float) * 1e6
        for state in states:
            obs = np.asarray(batch.observations[state], dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sim.loc[batch.time, state].to_numpy(dtype=float)
            sigma = sigma_for_state(state, obs[mask])
            residuals.append((pred[mask] - obs[mask]) / sigma)
    if l2_reference is not None and l2_lambda > 0.0 and l2_parameters:
        l2 = []
        scale = math.sqrt(float(l2_lambda))
        for name in l2_parameters:
            value = max(float(theta[name]), 1e-16)
            ref = max(float(l2_reference[name]), 1e-16)
            l2.append(scale * math.log(value / ref))
        residuals.append(np.asarray(l2, dtype=float))
    if not residuals:
        return np.array([], dtype=float)
    return np.concatenate(residuals)


def fit_parameters(
    label: str,
    batches: list[BatchData],
    base_theta: dict[str, float],
    parameters: tuple[str, ...],
    max_nfev: int,
    l2_reference: dict[str, float] | None = None,
    l2_parameters: tuple[str, ...] = tuple(),
    l2_lambda: float = 0.0,
) -> tuple[dict[str, float], dict[str, float]]:
    x0 = theta_to_log_vector(base_theta, parameters)
    lb, ub = log_bounds(parameters)
    x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)

    def fun(x):
        theta = log_vector_to_theta(x, parameters, base_theta)
        return residual_vector(theta, batches, l2_reference=l2_reference, l2_parameters=l2_parameters, l2_lambda=l2_lambda)

    start_resid = fun(x0)
    result = least_squares(
        fun,
        x0,
        bounds=(lb, ub),
        method="trf",
        x_scale="jac",
        loss="linear",
        max_nfev=int(max_nfev),
        ftol=2e-5,
        xtol=2e-5,
        gtol=2e-5,
        verbose=0,
    )
    theta_hat = log_vector_to_theta(result.x, parameters, base_theta)
    end_resid = fun(result.x)
    summary = {
        "fit": label,
        "n_batches": len(batches),
        "mediums": ",".join(sorted(set(batch.medium for batch in batches))),
        "parameters": ",".join(parameters),
        "n_parameters": len(parameters),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nfev": int(result.nfev),
        "initial_wsse": float(np.dot(start_resid, start_resid)),
        "final_wsse": float(np.dot(end_resid, end_resid)),
        "n_residuals": int(len(end_resid)),
        "dof": int(max(len(end_resid) - len(parameters), 0)),
        "wsse_per_residual": float(np.dot(end_resid, end_resid) / max(len(end_resid), 1)),
        "wsse_per_dof": float(np.dot(end_resid, end_resid) / max(len(end_resid) - len(parameters), 1)),
        "l2_lambda": float(l2_lambda),
        "l2_parameters": ",".join(l2_parameters),
    }
    return theta_hat, summary


def residual_state_summary(theta: dict[str, float], batches: list[BatchData], label: str) -> pd.DataFrame:
    rows = []
    for batch in batches:
        sim = simulate(batch, theta, batch.time)
        if sim is None:
            continue
        for state in STATE_NAMES:
            obs = np.asarray(batch.observations[state], dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sim.loc[batch.time, state].to_numpy(dtype=float)
            sigma = sigma_for_state(state, obs[mask])
            resid = pred[mask] - obs[mask]
            wresid = resid / sigma
            rows.append(
                {
                    "fit": label,
                    "medium": batch.medium,
                    "batch": batch.batch,
                    "state": state,
                    "n": int(mask.sum()),
                    "sse": float(np.dot(resid, resid)),
                    "wsse": float(np.dot(wresid, wresid)),
                    "rmse": float(np.sqrt(np.mean(resid**2))),
                    "weighted_rmse": float(np.sqrt(np.mean(wresid**2))),
                    "mean_bias": float(np.mean(resid)),
                }
            )
    return pd.DataFrame(rows)


def build_jacobian(
    theta: dict[str, float],
    parameters: tuple[str, ...],
    batches: list[BatchData],
    step: float = 1e-2,
) -> tuple[np.ndarray, np.ndarray]:
    base_resid = residual_vector(theta, batches)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), *PARAMETER_BOUNDS[name]))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), *PARAMETER_BOUNDS[name]))
        r_plus = residual_vector(theta_plus, batches)
        r_minus = residual_vector(theta_minus, batches)
        if len(r_plus) != len(base_resid) or len(r_minus) != len(base_resid):
            cols.append(np.zeros_like(base_resid))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
    if not cols:
        return np.empty((len(base_resid), 0)), base_resid
    return np.column_stack(cols), base_resid


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
    out = {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}max_eigenvalue": max_eig,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace": float(np.trace(fim)),
        f"{prefix}trace_inv": float(np.trace(cov)),
        f"{prefix}rank_1e-8": int(np.sum(eig > max_eig * 1e-8)) if np.isfinite(max_eig) and max_eig > 0.0 else 0,
    }
    return out


def fim_diagnostics(fim: np.ndarray, parameters: tuple[str, ...], theta: dict[str, float], label: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eigvals, eigvecs = np.linalg.eigh(fim)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    max_eig = float(np.max(eigvals)) if eigvals.size else np.nan
    eig_rows = []
    load_rows = []
    for idx, value in enumerate(eigvals):
        eig_rows.append(
            {
                "analysis": label,
                "direction": idx + 1,
                "eigenvalue": float(value),
                "relative_eigenvalue": float(value / max_eig) if max_eig > 0.0 else np.nan,
            }
        )
    for idx in range(min(6, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        abs_vec = np.abs(vec)
        top = np.argsort(abs_vec)[::-1][: min(7, len(parameters))]
        load_rows.append(
            {
                "analysis": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs_vec[i]:.3f}" for i in top),
            }
        )
    cov = stable_inverse(fim)
    param_rows = []
    for i, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[i, i], 0.0)))
        value = float(theta[name])
        lb, ub = PARAMETER_BOUNDS[name]
        active = value <= float(lb) * 1.01 or value >= float(ub) / 1.01
        if std_log <= 0.35 and not active:
            classification = "well_estimated"
        elif std_log <= 0.75 and not active:
            classification = "moderate"
        elif name in GLYCEROL_PARAMETERS and std_log <= 1.0 and not active:
            classification = "new_data_limited"
        else:
            classification = "weak_or_confounded"
        param_rows.append(
            {
                "analysis": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std_log,
                "approx_95_multiplier": float(math.exp(1.96 * min(std_log, 20.0))),
                "fim_diag": float(fim[i, i]),
                "active_bound": bool(active),
                "classification": classification,
            }
        )
    return pd.DataFrame(eig_rows), pd.DataFrame(load_rows), pd.DataFrame(param_rows)


def profile_parameters(
    theta_hat: dict[str, float],
    batches: list[BatchData],
    fit_parameters_subset: tuple[str, ...],
    profiled_parameters: tuple[str, ...],
    base_objective: float,
    max_nfev: int,
    grid_points: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    threshold = float(chi2.ppf(0.95, df=1))
    for profiled in profiled_parameters:
        theta_value = float(theta_hat[profiled])
        lb, ub = PARAMETER_BOUNDS[profiled]
        low = max(float(lb), theta_value / 2.0)
        high = min(float(ub), theta_value * 2.0)
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            continue
        grid = np.unique(np.concatenate([np.linspace(low, high, int(grid_points)), np.array([theta_value])]))
        nuisance = tuple(name for name in fit_parameters_subset if name != profiled)
        for fixed_value in sorted(float(v) for v in grid):
            fixed_base = dict(theta_hat)
            fixed_base[profiled] = fixed_value
            if nuisance:
                try:
                    theta_prof, summary = fit_parameters(
                        f"profile_{profiled}_{fixed_value:.5g}",
                        batches,
                        fixed_base,
                        nuisance,
                        max_nfev=max_nfev,
                    )
                except Exception as err:
                    rows.append(
                        {
                            "profiled_parameter": profiled,
                            "theta_value": fixed_value,
                            "success": False,
                            "objective": np.nan,
                            "lr_stat": np.nan,
                            "error": f"{type(err).__name__}: {err}",
                        }
                    )
                    continue
            else:
                theta_prof = fixed_base
                res = residual_vector(theta_prof, batches)
                summary = {"final_wsse": float(np.dot(res, res)), "success": True}
            obj = float(summary["final_wsse"])
            rows.append(
                {
                    "profiled_parameter": profiled,
                    "theta_hat": theta_value,
                    "theta_value": fixed_value,
                    "success": bool(summary.get("success", True)),
                    "objective": obj,
                    "lr_stat": max(0.0, obj - float(base_objective)),
                    "chi2_95_threshold": threshold,
                    "error": "",
                }
            )
    profile_df = pd.DataFrame(rows)
    summary_rows = []
    if not profile_df.empty:
        ok = profile_df[profile_df["success"]].copy()
        for parameter, group in ok.groupby("profiled_parameter"):
            theta0 = float(group["theta_hat"].iloc[0])
            left = group[group["theta_value"] < theta0]
            right = group[group["theta_value"] > theta0]
            summary_rows.append(
                {
                    "parameter": parameter,
                    "n_success": int(len(group)),
                    "max_lr_stat": float(group["lr_stat"].max()),
                    "crosses_left_95": bool((left["lr_stat"] >= threshold).any()),
                    "crosses_right_95": bool((right["lr_stat"] >= threshold).any()),
                    "profile_identifiable_95": bool(
                        (left["lr_stat"] >= threshold).any() and (right["lr_stat"] >= threshold).any()
                    ),
                    "chi2_95_threshold": threshold,
                }
            )
    return profile_df, pd.DataFrame(summary_rows)


def operational_sample_times(
    horizon_h: float,
    start: str = "2026-06-15 08:00",
    policy: str = "front_loaded",
) -> np.ndarray:
    start_dt = datetime.fromisoformat(start)
    horizon = float(horizon_h)
    points: list[float] = []
    process_day = 0
    current_day = start_dt.date()
    while True:
        day_start = datetime.combine(current_day, datetime.min.time())
        if (day_start - start_dt).total_seconds() / 3600.0 > horizon + 24.0:
            break
        if day_start.weekday() < 5:
            if policy == "balanced":
                hours = (10, 12, 14, 16)
            elif policy == "two_per_day":
                hours = (10, 16)
            else:
                hours = (10, 12, 14, 16) if process_day <= 2 else (10, 16)
            for hour in hours:
                dt = datetime.combine(current_day, datetime.min.time()) + timedelta(hours=hour)
                t = (dt - start_dt).total_seconds() / 3600.0
                if 0.0 < t <= horizon:
                    points.append(float(t))
            process_day += 1
        current_day = current_day + timedelta(days=1)
    return np.array(sorted(set(round(t, 8) for t in points)), dtype=float)


def nearest_operational_time(target_h: float, horizon_h: float) -> float:
    times = operational_sample_times(horizon_h, policy="balanced")
    if len(times) == 0:
        return float(np.clip(target_h, 0.0, horizon_h))
    return float(times[np.argmin(np.abs(times - float(target_h)))])


def make_future_designs(data: pd.DataFrame) -> list[FutureDesign]:
    summary = []
    for medium, group in data.groupby("medium"):
        summary.append(
            {
                "medium": medium,
                "X": float(group.groupby("batch")["X_viable_kg_m3"].first().median()),
                "Xd": float(group.groupby("batch")["X_dead_kg_m3"].first().median()),
                "N": float(group.groupby("batch")["N_kg_m3"].first().median()),
                "G": float(group.groupby("batch")["G_g_l"].first().median()),
                "F": float(group.groupby("batch")["F_g_l"].first().median()),
                "E": 0.0,
                "Gly": float(group.groupby("batch")["glycerol_g_l"].first().median()),
            }
        )
    medium_defaults = {row["medium"]: {k: row[k] for k in STATE_NAMES} for row in summary}
    natural = medium_defaults.get("natural", {"X": 0.5, "Xd": 0.02, "N": 0.22, "G": 78, "F": 77, "E": 0, "Gly": 0.8})
    synthetic = medium_defaults.get("synthetic", {"X": 0.45, "Xd": 0.05, "N": 0.22, "G": 124, "F": 118, "E": 0, "Gly": 0.3})

    def p(horizon: float, **kwargs) -> dict[str, tuple[tuple[float, float], ...]]:
        schedules = {channel: [] for channel in INPUT_CHANNELS}
        for channel, rows in kwargs.items():
            for target, amount in rows:
                schedules[channel].append((nearest_operational_time(float(target), horizon), float(amount)))
        return {channel: tuple(sorted(rows)) for channel, rows in schedules.items()}

    designs: list[FutureDesign] = []
    designs.append(
        FutureDesign(
            "natural_control_18C",
            "natural_baseline",
            "natural",
            216.0,
            dict(natural),
            (18.0, 18.0, 18.0, 18.0),
            p(216.0),
            "Natural must baseline to quantify transfer from the current wine matrix.",
        )
    )
    designs.append(
        FutureDesign(
            "natural_cold_hot_switch",
            "natural_temperature",
            "natural",
            216.0,
            dict(natural),
            (15.0, 18.0, 23.0, 20.0),
            p(216.0, N=((26.0, 0.04),)),
            "Natural matrix with temperature excitation and a modest YAN pulse.",
        )
    )
    designs.append(
        FutureDesign(
            "natural_glucose_pulse_after_growth",
            "natural_sugar_pulse",
            "natural",
            216.0,
            dict(natural),
            (18.0, 20.0, 22.0, 19.0),
            p(216.0, G=((50.0, 35.0),), N=((26.0, 0.03),)),
            "Natural must with glucose perturbation to test sugar-transfer validity.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_natural_like_control",
            "synthetic_baseline",
            "synthetic",
            216.0,
            {**dict(synthetic), "G": natural["G"], "F": natural["F"], "N": natural["N"]},
            (18.0, 18.0, 18.0, 18.0),
            p(216.0),
            "Synthetic must matched to natural macronutrients; isolates matrix effect.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_high_sugar_reference",
            "synthetic_baseline",
            "synthetic",
            216.0,
            dict(synthetic),
            (18.0, 20.0, 20.0, 18.0),
            p(216.0, N=((26.0, 0.05),)),
            "Synthetic high-sugar operating point close to the current synthetic dataset.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_glucose_rich_fructose_pulse",
            "sugar_separation",
            "synthetic",
            216.0,
            {**dict(synthetic), "G": 155.0, "F": 35.0, "N": 0.24},
            (17.0, 20.0, 23.0, 20.0),
            p(216.0, F=((74.0, 35.0),), N=((26.0, 0.05),)),
            "High glucose plus fructose pulse separates glucose and fructose uptake/yield directions.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_fructose_rich_glucose_pulse",
            "sugar_separation",
            "synthetic",
            216.0,
            {**dict(synthetic), "G": 35.0, "F": 155.0, "N": 0.24},
            (17.0, 20.0, 23.0, 20.0),
            p(216.0, G=((74.0, 35.0),), N=((26.0, 0.05),)),
            "Fructose-rich run probes iG and fructose kinetic directions.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_low_yan_ladder",
            "nitrogen_saturation",
            "synthetic",
            216.0,
            {**dict(synthetic), "N": 0.035, "G": 95.0, "F": 95.0},
            (15.0, 18.0, 22.0, 22.0),
            p(216.0, N=((26.0, 0.025), (50.0, 0.045), (74.0, 0.065))),
            "Low-YAN ladder targets sN and qN separation.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_low_sugar_saturation_scan",
            "sugar_saturation",
            "synthetic",
            168.0,
            {**dict(synthetic), "N": 0.16, "G": 25.0, "F": 25.0},
            (16.0, 20.0, 24.0, 20.0),
            p(168.0, G=((26.0, 20.0),), F=((50.0, 20.0),), N=((26.0, 0.03),)),
            "Low/intermediate sugar levels excite sG and sF instead of keeping saturation flat.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_ethanol_inhibition_challenge",
            "ethanol_inhibition",
            "synthetic",
            216.0,
            {**dict(synthetic), "E": 30.0, "G": 85.0, "F": 85.0, "N": 0.20},
            (18.0, 22.0, 25.0, 22.0),
            p(216.0, N=((26.0, 0.05),)),
            "Initial ethanol decouples ethanol inhibition from ethanol generated by fermentation.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_late_ethanol_death_probe",
            "death",
            "synthetic",
            216.0,
            {**dict(synthetic), "E": 5.0, "G": 80.0, "F": 80.0, "N": 0.16, "Xd": 0.02},
            (20.0, 25.0, 25.0, 22.0),
            p(216.0, N=((26.0, 0.04),), E=((98.0, 25.0),)),
            "Late ethanol stress plus Xd observation targets Kd0/iE.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_high_biomass_low_N_maintenance",
            "maintenance",
            "synthetic",
            216.0,
            {**dict(synthetic), "X": 2.0, "N": 0.03, "G": 75.0, "F": 75.0, "E": 10.0},
            (16.0, 18.0, 20.0, 18.0),
            p(216.0, G=((50.0, 25.0),), F=((98.0, 25.0),), N=((74.0, 0.015),)),
            "High biomass with low N creates low-growth sugar consumption windows for m0.",
        )
    )
    designs.append(
        FutureDesign(
            "synthetic_viable_biomass_step",
            "biomass_input",
            "synthetic",
            216.0,
            {**dict(synthetic), "X": 0.30, "N": 0.22, "G": 90.0, "F": 90.0},
            (18.0, 22.0, 22.0, 18.0),
            p(216.0, N=((26.0, 0.05),), X=((50.0, 1.0),)),
            "Known viable biomass addition tests rate proportionality to X and supports yield separation.",
        )
    )
    return designs


def future_fim(theta: dict[str, float], design: FutureDesign, parameters: tuple[str, ...], sample_policy: str, step: float) -> np.ndarray:
    sample_times = operational_sample_times(design.horizon_h, policy=sample_policy)
    base = simulate(design, theta, sample_times)
    if base is None:
        raise RuntimeError(f"Simulation failed for candidate {design.name}")
    base = base.loc[[t for t in sample_times if t in base.index]]
    columns = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), *PARAMETER_BOUNDS[name]))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), *PARAMETER_BOUNDS[name]))
        sim_plus = simulate(design, theta_plus, sample_times)
        sim_minus = simulate(design, theta_minus, sample_times)
        if sim_plus is None or sim_minus is None:
            columns.append(np.zeros(len(sample_times) * len(STATE_NAMES), dtype=float))
            continue
        pieces = []
        for state in STATE_NAMES:
            center = base[state].to_numpy(dtype=float)
            sigma = sigma_for_state(state, center)
            pieces.append((sim_plus.loc[base.index, state].to_numpy(dtype=float) - sim_minus.loc[base.index, state].to_numpy(dtype=float)) / (2.0 * step * sigma))
        columns.append(np.concatenate(pieces))
    jac = np.column_stack(columns)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def variance_reduction(prior_fim: np.ndarray, combined_fim: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior_fim)
    combined_cov = stable_inverse(combined_fim)
    rows = {}
    ratios = []
    for name in parameters:
        idx = parameters.index(name)
        before = float(prior_cov[idx, idx])
        after = float(combined_cov[idx, idx])
        ratio = after / before if before > 0.0 else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        if name in WEAK_FROM_OLD or name in GLYCEROL_PARAMETERS:
            ratios.append(ratio)
    if ratios:
        rows["weak_mean_var_reduction"] = float(1.0 - np.nanmean(ratios))
        rows["weak_worst_var_reduction"] = float(1.0 - np.nanmax(ratios))
    return rows


def score_fim(fim: np.ndarray, objective: str) -> float:
    metrics = fim_metrics(fim, FULL17)
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-18)
    if objective == "d_opt":
        return float(metrics["logdet"])
    if objective == "e_opt":
        return math.log(max(float(metrics["min_eigenvalue"]), 1e-18))
    if objective == "a_opt":
        return -math.log(max(float(metrics["trace_inv"]), 1e-18))
    return float(metrics["logdet"]) + 2.0 * math.log(min_rel) - 0.1 * math.log(max(float(metrics["trace_inv"]), 1e-18))


def greedy_campaign(
    candidate_fims: dict[str, np.ndarray],
    designs: dict[str, FutureDesign],
    prior_fim: np.ndarray,
    campaign_size: int,
    objective: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    selected: list[str] = []
    remaining = set(candidate_fims)
    current = prior_fim.copy()
    rows = []
    for idx in range(1, int(campaign_size) + 1):
        best_name = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + candidate_fims[name]
            score = score_fim(trial, objective)
            if score > best_score:
                best_name = name
                best_score = score
                best_metrics = {**fim_metrics(trial, FULL17, prefix="campaign_"), **variance_reduction(prior_fim, trial, FULL17)}
        if best_name is None:
            break
        current = current + candidate_fims[best_name]
        remaining.remove(best_name)
        selected.append(best_name)
        design = designs[best_name]
        rows.append(
            {
                "objective": objective,
                "campaign_order": idx,
                "candidate": best_name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{v:g}" for v in design.temperature_segments),
                "score": float(best_score),
                "rationale": design.rationale,
                **best_metrics,
            }
        )
    return pd.DataFrame(rows), current


class PyomoFutureExperiment:
    def __init__(self, theta: dict[str, float], design: FutureDesign, parameters: tuple[str, ...], sample_policy: str):
        self.theta = dict(theta)
        self.design = design
        self.parameters = tuple(parameters)
        self.sample_policy = sample_policy

    def get_labeled_model(self):
        return build_pyomo_discrete_model(self.theta, self.design, self.parameters, self.sample_policy)


def _pyomo_time_grid(design: FutureDesign, sample_policy: str) -> tuple[np.ndarray, set[float]]:
    sample_times = set(float(t) for t in operational_sample_times(design.horizon_h, policy=sample_policy))
    dense = set(np.arange(0.0, float(design.horizon_h) + 1e-9, 4.0).round(8))
    for schedule in design.pulses.values():
        for time_h, _amount in schedule:
            dense.add(round(max(0.0, float(time_h) - 4.0), 8))
            dense.add(round(float(time_h), 8))
            dense.add(round(min(float(design.horizon_h), float(time_h) + 4.0), 8))
    dense.update(sample_times)
    dense.add(0.0)
    return np.array(sorted(dense), dtype=float), sample_times


def build_pyomo_discrete_model(theta: dict[str, float], design: FutureDesign, parameters: tuple[str, ...], sample_policy: str):
    time, sample_times = _pyomo_time_grid(design, sample_policy)
    base_sim = simulate(design, theta, time)
    if base_sim is None:
        raise RuntimeError("Base simulation failed before building Pyomo DOE model.")
    m = pyo.ConcreteModel(f"new_must_glycerol_doe_{design.name}")
    m.K = pyo.RangeSet(0, len(time) - 1)
    m.Kpos = pyo.RangeSet(1, len(time) - 1)
    m.time = pyo.Param(m.K, initialize={k: float(time[k]) for k in range(len(time))})
    m.dt = pyo.Param(m.Kpos, initialize={k: float(time[k] - time[k - 1]) for k in range(1, len(time))})
    m.TempC = pyo.Param(m.K, initialize={k: temperature_at(design, float(time[k])) for k in range(len(time))})
    for channel in INPUT_CHANNELS:
        setattr(
            m,
            f"{channel}_input",
            pyo.Param(
                m.K,
                initialize={
                    k: (
                        0.0
                        if k == 0
                        else sum(
                            float(amount)
                            for event_time, amount in design.pulses.get(channel, tuple())
                            if float(time[k - 1]) < float(event_time) <= float(time[k])
                            or (k == 1 and float(event_time) == float(time[0]))
                        )
                        / float(time[k] - time[k - 1])
                    )
                    for k in range(len(time))
                },
            ),
        )

    for name in FULL17:
        var = pyo.Var(bounds=PARAMETER_BOUNDS[name], initialize=float(theta[name]))
        setattr(m, name, var)
        var.fix(float(theta[name]))

    for state in STATE_NAMES:
        lb, ub = STATE_BOUNDS[state]
        init = {k: float(base_sim.iloc[k][state]) for k in range(len(time))}
        setattr(m, state, pyo.Var(m.K, bounds=(lb, ub), initialize=init))

    for state in STATE_NAMES:
        state_var = getattr(m, state)
        state_var[0].fix(float(design.initials[state]))

    eps = 1e-8

    def terms_at(_m, k):
        temp_k = _m.TempC[k] + 273.15
        n = _m.N[k]
        g = _m.G[k]
        f = _m.F[k]
        e = _m.E[k]
        a_mu = pyo.exp(FIXED_CONSTANTS["Eac"] * (temp_k - 300.0) / (300.0 * FIXED_CONSTANTS["R"] * temp_k))
        a_beta = pyo.exp(FIXED_CONSTANTS["Eafe"] * (temp_k - 296.15) / (296.15 * FIXED_CONSTANTS["R"] * temp_k))
        a_kn = pyo.exp(FIXED_CONSTANTS["EaKn"] * (temp_k - 293.15) / (293.15 * FIXED_CONSTANTS["R"] * temp_k))
        a_kg = pyo.exp(FIXED_CONSTANTS["EaKg"] * (temp_k - 293.15) / (293.15 * FIXED_CONSTANTS["R"] * temp_k))
        a_kf = pyo.exp(FIXED_CONSTANTS["EaKf"] * (temp_k - 293.15) / (293.15 * FIXED_CONSTANTS["R"] * temp_k))
        a_kig = pyo.exp(FIXED_CONSTANTS["EaKig"] * (temp_k - 293.15) / (293.15 * FIXED_CONSTANTS["R"] * temp_k))
        a_kie = pyo.exp(FIXED_CONSTANTS["EaKie"] * (temp_k - 293.15) / (293.15 * FIXED_CONSTANTS["R"] * temp_k))
        kn = (_m.mu0 / _m.sN) * a_kn
        kg = (_m.betaG0 / _m.sG) * a_kg
        kf = (_m.betaF0 / _m.sF) * a_kf
        i_g = _m.iG / (a_kig + eps)
        i_e = _m.iE / (a_kie + eps)
        n_lim = n / (n + kn + eps)
        g_lim = g / (g + kg + eps)
        f_lim = f / (f + kf + eps)
        e_inhib = 1.0 / (1.0 + i_e * e)
        g_inhib_for_f = 1.0 / (1.0 + i_g * g)
        growth_factor = a_mu * n_lim
        glucose_ferm_factor = a_beta * g_lim * e_inhib
        fructose_ferm_factor = a_beta * f_lim * g_inhib_for_f * e_inhib
        mu = _m.mu0 * growth_factor
        beta_g = _m.betaG0 * glucose_ferm_factor
        beta_f = _m.betaF0 * fructose_ferm_factor
        maintenance = _m.m0 * pyo.exp(FIXED_CONSTANTS["Eam"] * (temp_k - 293.3) / (293.3 * FIXED_CONSTANTS["R"] * temp_k))
        td = -0.0001 * e**3 + 0.0049 * e**2 - 0.1279 * e + 315.89
        kd = pyo.Expr_if(
            IF=(temp_k >= td),
            THEN=_m.Kd0
            * pyo.exp(
                (FIXED_CONSTANTS["Cde"] * e)
                + (FIXED_CONSTANTS["Etd"] * (temp_k - 305.65)) / (305.65 * FIXED_CONSTANTS["R"] * temp_k)
            ),
            ELSE=0.0,
        )
        sugar_total = g + f + eps
        maintenance_availability = sugar_total / (sugar_total + MAINTENANCE_SUGAR_CUTOFF_KG_M3)
        return (
            growth_factor,
            glucose_ferm_factor,
            fructose_ferm_factor,
            mu,
            beta_g,
            beta_f,
            kd,
            maintenance,
            maintenance_availability,
            sugar_total,
        )

    @m.Constraint(m.Kpos)
    def X_balance(_m, k):
        growth_factor, _gf, _ff, mu, _bg, _bf, kd, _maint, _ma, _st = terms_at(_m, k)
        return (_m.X[k] - _m.X[k - 1]) / _m.dt[k] == (mu - kd) * _m.X[k] + _m.X_input[k]

    @m.Constraint(m.Kpos)
    def Xd_balance(_m, k):
        _gf0, _gf, _ff, _mu, _bg, _bf, kd, _maint, _ma, _st = terms_at(_m, k)
        return (_m.Xd[k] - _m.Xd[k - 1]) / _m.dt[k] == kd * _m.X[k]

    @m.Constraint(m.Kpos)
    def N_balance(_m, k):
        growth_factor, _gf, _ff, _mu, _bg, _bf, _kd, _maint, _ma, _st = terms_at(_m, k)
        return (_m.N[k] - _m.N[k - 1]) / _m.dt[k] == -_m.qN * growth_factor * _m.X[k] + _m.N_input[k]

    @m.Constraint(m.Kpos)
    def G_balance(_m, k):
        growth_factor, glucose_factor, _ff, _mu, _bg, _bf, _kd, maintenance, availability, sugar_total = terms_at(_m, k)
        return (_m.G[k] - _m.G[k - 1]) / _m.dt[k] == -(
            _m.qXG * growth_factor + _m.qEG * glucose_factor + maintenance * availability * (_m.G[k] / sugar_total)
        ) * _m.X[k] + _m.G_input[k]

    @m.Constraint(m.Kpos)
    def F_balance(_m, k):
        growth_factor, _gf, fructose_factor, _mu, _bg, _bf, _kd, maintenance, availability, sugar_total = terms_at(_m, k)
        return (_m.F[k] - _m.F[k - 1]) / _m.dt[k] == -(
            _m.qXF * growth_factor + _m.qEF * fructose_factor + maintenance * availability * (_m.F[k] / sugar_total)
        ) * _m.X[k] + _m.F_input[k]

    @m.Constraint(m.Kpos)
    def E_balance(_m, k):
        _growth, _gf, _ff, _mu, beta_g, beta_f, _kd, _maint, _ma, _st = terms_at(_m, k)
        return (_m.E[k] - _m.E[k - 1]) / _m.dt[k] == (beta_g + beta_f) * _m.X[k] + _m.E_input[k]

    @m.Constraint(m.Kpos)
    def Gly_balance(_m, k):
        _growth, glucose_factor, fructose_factor, _mu, _bg, _bf, _kd, _maint, _ma, _st = terms_at(_m, k)
        return (_m.Gly[k] - _m.Gly[k - 1]) / _m.dt[k] == (
            _m.gammaG0 * glucose_factor + _m.gammaF0 * fructose_factor
        ) * _m.X[k]

    m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for name in parameters:
        var = getattr(m, name)
        m.unknown_parameters[var] = float(theta[name])

    m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for k, time_h in enumerate(time):
        if float(time_h) not in sample_times:
            continue
        for state in STATE_NAMES:
            var = getattr(m, state)[k]
            simulated = float(base_sim.iloc[k][state])
            m.experiment_outputs[var] = simulated
            m.measurement_error[var] = float(sigma_for_state(state, simulated))
    return m


def run_pyomo_doe_check(theta: dict[str, float], design: FutureDesign, parameters: tuple[str, ...], sample_policy: str, step: float) -> dict[str, float]:
    solver = pyo.SolverFactory("ipopt")
    solver.options["max_iter"] = 3000
    solver.options["tol"] = 1e-6
    experiment = PyomoFutureExperiment(theta, design, parameters, sample_policy)
    doe = DesignOfExperiments(
        experiment=experiment,
        step=float(step),
        scale_nominal_param_value=True,
        solver=solver,
        tee=False,
    )
    fim = np.asarray(doe.compute_FIM(method="sequential"), dtype=float)
    pd.DataFrame(fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / f"pyomo_doe_fim_{design.name}.csv")
    return {
        "candidate": design.name,
        "status": "ok",
        "parameters": ",".join(parameters),
        **fim_metrics(fim, parameters, prefix="pyomo_"),
    }


def write_fit_plots(theta: dict[str, float], batches: list[BatchData], label: str, max_batches_per_medium: int = 3) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    selected = []
    for medium in sorted(set(batch.medium for batch in batches)):
        medium_batches = [batch for batch in batches if batch.medium == medium]
        selected.extend(medium_batches[:max_batches_per_medium])
    for batch in selected:
        sim = simulate(batch, theta, batch.time)
        if sim is None:
            continue
        fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes, STATE_NAMES):
            obs = batch.observations[state]
            mask = np.isfinite(obs)
            ax.plot(sim.index, sim[state], "-", lw=2, label="model")
            ax.scatter(batch.time[mask], obs[mask], s=22, label="data")
            ax.set_title(state)
            ax.grid(True, alpha=0.25)
        axes[-1].axis("off")
        axes[0].legend(loc="best")
        fig.suptitle(f"{label}: {batch.label}", y=0.995)
        fig.tight_layout()
        fig.savefig(plot_dir / f"fit_{label}_{batch.medium}_{batch.batch}.png", dpi=160)
        plt.close(fig)


def write_eigen_plot(eigen_table: pd.DataFrame) -> None:
    if eigen_table.empty:
        return
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, group in eigen_table.groupby("analysis"):
        ax.semilogy(group["direction"], np.maximum(group["relative_eigenvalue"], 1e-16), marker="o", label=label)
    ax.set_xlabel("FIM eigen-direction, ordered weak to strong")
    ax.set_ylabel("relative eigenvalue")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(plot_dir / "fim_eigen_spectra.png", dpi=170)
    plt.close(fig)


def write_campaign_input_plots(selected: pd.DataFrame, designs: dict[str, FutureDesign], sample_policy: str) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for name in selected["candidate"].drop_duplicates():
        design = designs[str(name)]
        fig, axes = plt.subplots(7, 1, figsize=(10, 12), sharex=True)
        time = np.linspace(0.0, design.horizon_h, 500)
        axes[0].step(time, [temperature_at(design, t) for t in time], where="post", color="tab:red")
        axes[0].set_ylabel("T [C]")
        sample_times = operational_sample_times(design.horizon_h, policy=sample_policy)
        for t in sample_times:
            axes[0].axvline(t, color="0.85", lw=0.5, zorder=0)
        init_states = ["X", "N", "G", "F", "E", "Gly"]
        axes[1].bar(init_states, [design.initials[s] for s in init_states], color="0.35")
        axes[1].set_ylabel("initial")
        for ax, channel in zip(axes[2:], INPUT_CHANNELS):
            rows = design.pulses.get(channel, tuple())
            if rows:
                ax.vlines([t for t, _a in rows], 0.0, [a for _t, a in rows], lw=3)
                ax.scatter([t for t, _a in rows], [a for _t, a in rows], s=35)
            else:
                ax.axhline(0.0, color="0.85", lw=1)
            ax.set_ylabel(channel)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("process time [h]")
        for ax in axes:
            ax.set_xlim(0.0, design.horizon_h)
        fig.suptitle("\n".join(textwrap.wrap(f"{design.name}: {design.rationale}", 100)), fontsize=10, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(plot_dir / f"campaign_inputs_{design.name}.png", dpi=170)
        plt.close(fig)


def write_report(
    fit_summary: pd.DataFrame,
    param_summary: pd.DataFrame,
    profile_summary: pd.DataFrame,
    candidate_ranking: pd.DataFrame,
    selected_hybrid: pd.DataFrame,
    selected_dopt: pd.DataFrame,
    pyomo_row: dict[str, float],
) -> None:
    best_fit = fit_summary.sort_values("wsse_per_residual").iloc[0] if not fit_summary.empty else None
    lines = [
        "# New must glycerol estimability and DOE report",
        "",
        "## Scope",
        "",
        "This run uses the new natural and synthetic must databases, excludes CO2 and aroma states for now, and extends the kinetic state vector with dead biomass and glycerol.",
        "",
        "State vector: `X, Xd, N, G, F, E, Gly`.",
        "",
        "New glycerol equations:",
        "",
        "$$\\frac{dGly}{dt} = \\left(\\gamma_{G,0}\\,\\phi_G + \\gamma_{F,0}\\,\\phi_F\\right)X$$",
        "",
        "where `phi_G` and `phi_F` are the same temperature, substrate saturation, glucose/fructose interaction, and ethanol inhibition factors used by the ethanol-production terms.",
        "",
        "## Calibration summary",
        "",
        fit_summary.to_markdown(index=False) if not fit_summary.empty else "_No fit summary was generated._",
        "",
    ]
    if best_fit is not None:
        lines += [
            f"Lowest normalized weighted objective in this run: `{best_fit['fit']}` with WSSE/residual `{float(best_fit['wsse_per_residual']):.3g}`.",
            "",
        ]
    lines += [
        "## Parameter estimability",
        "",
        param_summary.sort_values(["analysis", "classification", "std_log_approx"]).to_markdown(index=False)
        if not param_summary.empty
        else "_No parameter diagnostics were generated._",
        "",
        "## Profile likelihood check",
        "",
        "When `profile_set` is present, `core` profiles refit the reduced kinetic nuisance set, while `weak` profiles use a more flexible nuisance set around the weak directions. The weak set is the conservative diagnostic for practical identifiability.",
        "",
        profile_summary.to_markdown(index=False) if not profile_summary.empty else "_Profiles were skipped or failed._",
        "",
        "## DOE candidate ranking",
        "",
        candidate_ranking.head(12)[
            [
                "candidate",
                "family",
                "medium",
                "combined_logdet",
                "combined_min_relative_eigenvalue",
                "combined_condition_number",
                "weak_mean_var_reduction",
                "weak_worst_var_reduction",
            ]
        ].to_markdown(index=False)
        if not candidate_ranking.empty
        else "_No candidate ranking was generated._",
        "",
        "## Selected campaign: hybrid D/E/A score",
        "",
        selected_hybrid[
            [
                "campaign_order",
                "candidate",
                "family",
                "medium",
                "campaign_logdet",
                "campaign_min_relative_eigenvalue",
                "weak_mean_var_reduction",
                "weak_worst_var_reduction",
                "rationale",
            ]
        ].to_markdown(index=False)
        if not selected_hybrid.empty
        else "_No hybrid campaign was selected._",
        "",
        "## Pure D-optimal benchmark",
        "",
        selected_dopt[
            [
                "campaign_order",
                "candidate",
                "family",
                "medium",
                "campaign_logdet",
                "campaign_min_relative_eigenvalue",
                "weak_mean_var_reduction",
                "weak_worst_var_reduction",
            ]
        ].to_markdown(index=False)
        if not selected_dopt.empty
        else "_No D-optimal campaign was selected._",
        "",
        "## Pyomo.DoE check",
        "",
        pd.DataFrame([pyomo_row]).to_markdown(index=False) if pyomo_row else "_Pyomo.DoE check was skipped or failed._",
        "",
        "## Operational interpretation",
        "",
        "- Sampling policy used for future designs: Monday-Friday only, 10:00-16:00, four samples per day during the first three process weekdays and two samples per process weekday afterwards.",
        "- Pulse/adition times were snapped to feasible working-window times; no night or weekend manual action is required in the proposed candidate set.",
        "- Temperature setpoints are assumed automatically programmable.",
        "- The selected mixed campaign should be interpreted as model-based screening, not as final protocol approval; the first block of three fermentations should be used to update the prior before committing to the remaining blocks.",
        "",
        "## Main conclusions to verify experimentally",
        "",
        "- Natural and synthetic musts excite different regions: natural contributes matrix-transfer relevance and lower-sugar dynamics; synthetic contributes stronger sugar and pulse perturbations.",
        "- A mixed approach is preferred unless the medium-transfer residuals show that one medium has a systematic model mismatch that the current structure cannot absorb.",
        "- Glycerol observations directly support `gammaG0` and `gammaF0`, but these parameters remain coupled to the glucose/fructose fermentation factors; sugar-composition perturbations are therefore required.",
        "- `Kd0`, `iE`, and `m0` need targeted stress or low-growth windows; ordinary control fermentations are insufficient.",
    ]
    (RESULTS_DIR / "new_must_glycerol_estimability_doe_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_notebook() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# New must glycerol estimability and optimal experimental design\n\n"
            "This notebook documents the calibration-prior, estimability, and model-based design workflow for the new natural and synthetic must datasets. "
            "CO2 and aroma states are intentionally excluded in this iteration; dead biomass and glycerol are included."
        ),
        nbformat.v4.new_markdown_cell(
            "## Model extension\n\n"
            "The state vector is\n\n"
            "$$z(t)=\\left[X, X_d, N, G, F, E, Gly\\right]^T.$$ \n\n"
            "The base kinetic equations retain the effective Zenteno formulation used in the previous notebook. Dead biomass is represented as\n\n"
            "$$\\frac{dX_d}{dt}=K_d(T,E)X,$$\n\n"
            "and viable biomass as\n\n"
            "$$\\frac{dX}{dt}=\\left(\\mu-K_d\\right)X + u_X(t).$$\n\n"
            "Glycerol is added as a fermentation-associated by-product:\n\n"
            "$$\\frac{dGly}{dt}=\\left(\\gamma_{G,0}\\,\\phi_G(T,G,E)+\\gamma_{F,0}\\,\\phi_F(T,F,G,E)\\right)X.$$\n\n"
            "The factors $\\phi_G$ and $\\phi_F$ are the same effective fermentation factors that multiply the ethanol rates. "
            "This makes glycerol informative for the same sugar/ethanol-inhibition subspace, rather than treating it as an independent empirical curve."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "import matplotlib.pyplot as plt\n\n"
            "ROOT = Path.cwd()\n"
            "if (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'shared/results/new_must_glycerol_estimability_doe'\n"
            "fit_summary = pd.read_csv(RESULTS / 'fit_summary.csv')\n"
            "param_summary = pd.read_csv(RESULTS / 'parameter_estimability_summary.csv')\n"
            "candidate_ranking = pd.read_csv(RESULTS / 'candidate_ranking.csv')\n"
            "selected_hybrid = pd.read_csv(RESULTS / 'selected_campaign_hybrid.csv')\n"
            "selected_dopt = pd.read_csv(RESULTS / 'selected_campaign_d_opt.csv')\n"
            "fit_summary"
        ),
        nbformat.v4.new_markdown_cell(
            "## Calibration prior\n\n"
            "Three calibration views are fitted: natural-only, synthetic-only, and mixed. "
            "The practical prior used for DOE is the mixed full-parameter L2 fit, because it preserves the complete parameter covariance while avoiding uncontrolled drift of weak directions."
        ),
        nbformat.v4.new_code_cell(
            "display(fit_summary.sort_values('final_wsse'))\n"
            "display(param_summary.sort_values(['analysis', 'std_log_approx']).head(30))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Medium comparison\n\n"
            "A parameter is considered reliable only if it is locally informed by the FIM and does not sit on an active bound. "
            "Medium-specific differences should be interpreted as structural/matrix-transfer evidence, not only as noise: natural and synthetic musts occupy different sugar and matrix regimes."
        ),
        nbformat.v4.new_code_cell(
            "pivot = param_summary.pivot_table(index='parameter', columns='analysis', values='std_log_approx', aggfunc='first')\n"
            "pivot.sort_index()"
        ),
        nbformat.v4.new_markdown_cell(
            "## Profile likelihood\n\n"
            "The profile likelihood check refits nuisance parameters while fixing selected parameters over a bounded local grid. "
            "A parameter is practically identifiable at 95% confidence when the likelihood-ratio curve crosses the $\\chi^2_1(0.95)$ threshold on both sides of the estimate."
        ),
        nbformat.v4.new_code_cell(
            "profile_path = RESULTS / 'profile_summary_combined.csv'\n"
            "if not profile_path.exists():\n"
            "    profile_path = RESULTS / 'profile_summary.csv'\n"
            "pd.read_csv(profile_path) if profile_path.exists() else 'profile summary not available'"
        ),
        nbformat.v4.new_markdown_cell(
            "## Model-based design of experiments\n\n"
            "Candidate experiments are evaluated by their expected Fisher information matrix. "
            "The current-data FIM is used as prior information, and each new candidate contributes an additive FIM under the local linear approximation.\n\n"
            "The main objective used for the proposed campaign is a hybrid deterministic score:\n\n"
            "$$\\Psi = \\log\\det(F) + 2\\log\\left(\\frac{\\lambda_{min}(F)}{\\lambda_{max}(F)}\\right) - 0.1\\log\\operatorname{tr}(F^{-1}).$$\n\n"
            "This keeps D-optimality as the main information-volume criterion while penalizing designs that leave a very weak eigen-direction. "
            "A pure D-optimal benchmark is also reported."
        ),
        nbformat.v4.new_code_cell(
            "display(candidate_ranking.head(12))\n"
            "display(selected_hybrid)\n"
            "display(selected_dopt)"
        ),
        nbformat.v4.new_markdown_cell(
            "## Input plots\n\n"
            "The generated PNG files in `shared/results/new_must_glycerol_estimability_doe/plots` show temperature setpoints, initial states, pulse timings, and feasible sampling windows for each selected candidate."
        ),
        nbformat.v4.new_code_cell(
            "from IPython.display import Image, display\n"
            "for png in sorted((RESULTS / 'plots').glob('campaign_inputs_*.png'))[:5]:\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Pyomo.DoE usage\n\n"
            "The runner includes a discrete Pyomo model with the same seven-state kinetics and labels `unknown_parameters`, `experiment_outputs`, and `measurement_error` for `pyomo.contrib.doe.DesignOfExperiments`. "
            "The saved `pyomo_doe_check.csv` records the sequential FIM check for the first selected experiment when Ipopt is available."
        ),
        nbformat.v4.new_code_cell(
            "pyomo_path = RESULTS / 'pyomo_doe_check.csv'\n"
            "pd.read_csv(pyomo_path) if pyomo_path.exists() else 'Pyomo.DoE check not available'"
        ),
        nbformat.v4.new_markdown_cell(
            "## How to interpret the output\n\n"
            "Use the calibration tables to decide which parameters are already estimable from current data. "
            "Use the campaign table to decide which new fermentations most reduce weak parameter variances and improve the smallest FIM eigenvalues. "
            "After the first block of three fermentations, rerun this notebook using those data as prior information before executing the next block."
        ),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="New natural/synthetic must glycerol estimability and DOE workflow.")
    parser.add_argument("--max-nfev", type=int, default=120)
    parser.add_argument("--full-max-nfev", type=int, default=90)
    parser.add_argument("--profile-max-nfev", type=int, default=45)
    parser.add_argument("--profile-grid", type=int, default=5)
    parser.add_argument("--max-profile-params", type=int, default=8)
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--sample-policy", choices=["front_loaded", "balanced", "two_per_day"], default="front_loaded")
    parser.add_argument("--sensitivity-step", type=float, default=1e-2)
    parser.add_argument("--skip-profiles", action="store_true")
    parser.add_argument("--skip-pyomo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    data = load_normalized_data()
    batches = make_batches(data)
    natural_batches = [batch for batch in batches if batch.medium == "natural"]
    synthetic_batches = [batch for batch in batches if batch.medium == "synthetic"]
    mixed_batches = list(batches)
    theta0 = complete_initial_theta(batches)
    pd.Series(theta0, name="initial").to_csv(RESULTS_DIR / "theta_initial.csv")

    fit_rows = []
    state_rows = []
    theta_rows = []
    group_fits: dict[str, dict[str, float]] = {}
    for label, group_batches in [
        ("natural_reduced11", natural_batches),
        ("synthetic_reduced11", synthetic_batches),
        ("mixed_reduced11", mixed_batches),
    ]:
        print(f"[fit] {label}: {len(group_batches)} batches", flush=True)
        theta_hat, summary = fit_parameters(label, group_batches, theta0, REDUCED11, max_nfev=args.max_nfev)
        group_fits[label] = theta_hat
        fit_rows.append(summary)
        theta_rows.append({"fit": label, **theta_hat})
        state_rows.append(residual_state_summary(theta_hat, group_batches, label))
        pd.Series(theta_hat, name=label).to_csv(RESULTS_DIR / f"theta_{label}.csv")

    print("[fit] mixed_full17_l2", flush=True)
    theta_full, full_summary = fit_parameters(
        "mixed_full17_l2",
        mixed_batches,
        group_fits["mixed_reduced11"],
        FULL17,
        max_nfev=args.full_max_nfev,
        l2_reference=theta0,
        l2_parameters=tuple(name for name in FULL17 if name not in CORE_FIT),
        l2_lambda=1.0,
    )
    group_fits["mixed_full17_l2"] = theta_full
    fit_rows.append(full_summary)
    theta_rows.append({"fit": "mixed_full17_l2", **theta_full})
    state_rows.append(residual_state_summary(theta_full, mixed_batches, "mixed_full17_l2"))
    pd.Series(theta_full, name="mixed_full17_l2").to_csv(RESULTS_DIR / "theta_mixed_full17_l2.csv")

    fit_summary = pd.DataFrame(fit_rows)
    fit_summary.to_csv(RESULTS_DIR / "fit_summary.csv", index=False)
    pd.DataFrame(theta_rows).to_csv(RESULTS_DIR / "theta_fit_table.csv", index=False)
    residual_by_state = pd.concat(state_rows, ignore_index=True) if state_rows else pd.DataFrame()
    residual_by_state.to_csv(RESULTS_DIR / "residual_by_state.csv", index=False)

    write_fit_plots(theta_full, mixed_batches, "mixed_full17_l2")

    fim_rows = []
    eigen_tables = []
    loading_tables = []
    parameter_tables = []
    group_specs = [
        ("natural_current", natural_batches, group_fits["natural_reduced11"]),
        ("synthetic_current", synthetic_batches, group_fits["synthetic_reduced11"]),
        ("mixed_current_reducedprior", mixed_batches, group_fits["mixed_reduced11"]),
        ("mixed_current_full_l2", mixed_batches, theta_full),
    ]
    prior_fim = None
    for label, group_batches, theta in group_specs:
        print(f"[fim] {label}", flush=True)
        jac, resid = build_jacobian(theta, FULL17, group_batches, step=args.sensitivity_step)
        fim = jac.T @ jac
        fim = 0.5 * (fim + fim.T)
        if label == "mixed_current_full_l2":
            prior_fim = fim
        pd.DataFrame(fim, index=FULL17, columns=FULL17).to_csv(RESULTS_DIR / f"fim_{label}.csv")
        fim_rows.append({"analysis": label, "n_residuals": int(len(resid)), **fim_metrics(fim, FULL17)})
        eig, load, param = fim_diagnostics(fim, FULL17, theta, label)
        eigen_tables.append(eig)
        loading_tables.append(load)
        parameter_tables.append(param)

    fim_summary = pd.DataFrame(fim_rows)
    fim_summary.to_csv(RESULTS_DIR / "fim_summary.csv", index=False)
    eigen_table = pd.concat(eigen_tables, ignore_index=True)
    loading_table = pd.concat(loading_tables, ignore_index=True)
    param_summary = pd.concat(parameter_tables, ignore_index=True)
    eigen_table.to_csv(RESULTS_DIR / "fim_eigenvalues.csv", index=False)
    loading_table.to_csv(RESULTS_DIR / "fim_weak_eigendirections.csv", index=False)
    param_summary.to_csv(RESULTS_DIR / "parameter_estimability_summary.csv", index=False)
    write_eigen_plot(eigen_table)

    profile_df = pd.DataFrame()
    profile_summary = pd.DataFrame()
    if not args.skip_profiles:
        weak_order = (
            param_summary[param_summary["analysis"].eq("mixed_current_full_l2")]
            .sort_values("std_log_approx", ascending=False)["parameter"]
            .tolist()
        )
        profile_params = []
        for name in weak_order:
            if name not in profile_params:
                profile_params.append(name)
        for name in REDUCED11:
            if name not in profile_params:
                profile_params.append(name)
        profile_params = tuple(profile_params[: int(args.max_profile_params)])
        base_resid = residual_vector(theta_full, mixed_batches)
        base_obj = float(np.dot(base_resid, base_resid))
        print(f"[profiles] {profile_params}", flush=True)
        profile_df, profile_summary = profile_parameters(
            theta_full,
            mixed_batches,
            REDUCED11,
            profile_params,
            base_obj,
            max_nfev=args.profile_max_nfev,
            grid_points=args.profile_grid,
        )
    profile_df.to_csv(RESULTS_DIR / "profile_profiles.csv", index=False)
    profile_summary.to_csv(RESULTS_DIR / "profile_summary.csv", index=False)

    if prior_fim is None:
        raise RuntimeError("Prior FIM was not computed.")

    designs = make_future_designs(data)
    design_by_name = {design.name: design for design in designs}
    candidate_rows = []
    candidate_fims: dict[str, np.ndarray] = {}
    for idx, design in enumerate(designs, start=1):
        print(f"[candidate {idx}/{len(designs)}] {design.name}", flush=True)
        try:
            fim = future_fim(theta_full, design, FULL17, args.sample_policy, args.sensitivity_step)
            candidate_fims[design.name] = fim
            combined = prior_fim + fim
            row = {
                "candidate": design.name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "status": "ok",
                "error": "",
                **fim_metrics(fim, FULL17, prefix="new_"),
                **fim_metrics(combined, FULL17, prefix="combined_"),
                **variance_reduction(prior_fim, combined, FULL17),
            }
        except Exception as err:
            row = {
                "candidate": design.name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "status": "failed",
                "error": f"{type(err).__name__}: {err}",
            }
        candidate_rows.append(row)
    candidate_ranking = pd.DataFrame(candidate_rows)
    if "combined_logdet" in candidate_ranking:
        candidate_ranking = candidate_ranking.sort_values("combined_logdet", ascending=False)
    candidate_ranking.to_csv(RESULTS_DIR / "candidate_ranking.csv", index=False)
    for name, fim in candidate_fims.items():
        pd.DataFrame(fim, index=FULL17, columns=FULL17).to_csv(RESULTS_DIR / f"candidate_fim_{name}.csv")

    selected_hybrid, final_hybrid_fim = greedy_campaign(
        candidate_fims,
        design_by_name,
        prior_fim,
        campaign_size=args.campaign_size,
        objective="hybrid",
    )
    selected_dopt, final_dopt_fim = greedy_campaign(
        candidate_fims,
        design_by_name,
        prior_fim,
        campaign_size=args.campaign_size,
        objective="d_opt",
    )
    selected_hybrid.to_csv(RESULTS_DIR / "selected_campaign_hybrid.csv", index=False)
    selected_dopt.to_csv(RESULTS_DIR / "selected_campaign_d_opt.csv", index=False)
    pd.DataFrame(final_hybrid_fim, index=FULL17, columns=FULL17).to_csv(RESULTS_DIR / "final_campaign_fim_hybrid.csv")
    pd.DataFrame(final_dopt_fim, index=FULL17, columns=FULL17).to_csv(RESULTS_DIR / "final_campaign_fim_d_opt.csv")
    final_eig, final_load, final_param = fim_diagnostics(final_hybrid_fim, FULL17, theta_full, "mixed_prior_plus_hybrid_campaign")
    final_eig.to_csv(RESULTS_DIR / "final_campaign_eigenvalues_hybrid.csv", index=False)
    final_load.to_csv(RESULTS_DIR / "final_campaign_weak_eigendirections_hybrid.csv", index=False)
    final_param.to_csv(RESULTS_DIR / "final_campaign_parameter_estimability_hybrid.csv", index=False)
    write_campaign_input_plots(selected_hybrid, design_by_name, args.sample_policy)

    pyomo_row: dict[str, float] = {}
    if not args.skip_pyomo and not selected_hybrid.empty:
        pyomo_candidate = design_by_name[str(selected_hybrid.iloc[0]["candidate"])]
        try:
            print(f"[pyomo doe] {pyomo_candidate.name}", flush=True)
            pyomo_row = run_pyomo_doe_check(theta_full, pyomo_candidate, REDUCED11, args.sample_policy, args.sensitivity_step)
        except Exception as err:
            pyomo_row = {
                "candidate": pyomo_candidate.name,
                "status": "failed",
                "error": f"{type(err).__name__}: {err}",
            }
    pd.DataFrame([pyomo_row]).to_csv(RESULTS_DIR / "pyomo_doe_check.csv", index=False)

    write_report(fit_summary, param_summary, profile_summary, candidate_ranking, selected_hybrid, selected_dopt, pyomo_row)
    write_notebook()
    print(f"[done] Results written to {RESULTS_DIR}", flush=True)
    print(f"[done] Notebook written to {NOTEBOOK_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
