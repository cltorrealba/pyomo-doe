from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PILOT_AROMA_PARTITION_MODE", "water_ethanol_total_sugar_as_glucose")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

try:
    from pyomo.contrib.doe.utils import rescale_FIM
except Exception:  # pragma: no cover
    rescale_FIM = None

SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SCRIPT_DIR.parent
SUPPORT_DIR = SCRIPT_DIR / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from shared import run_new_must_glycerol_estimability_doe as base
import run_pilot_2025_aroma_model_selection_doe as aroma_sel
import run_pilot_2025_calibration_estimability as pilot
import run_pilot_2025_co2_stripping_benchmark as old_co2
import run_pilot_2025_global_model_selection_doe as global_sel
from shared import run_secondary_joint_campaign_doe as joint
from shared import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "co2_solubility_integrated_doe"
NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_co2_solubility_integrated_doe.ipynb"
EXECUTED_NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_co2_solubility_integrated_doe.executed.ipynb"
PLOT_DIR = RESULTS_DIR / "plots"

PARTITION_MODE = "water_ethanol_total_sugar_as_glucose"
CO2_DOWNSAMPLE_H = 2.0
CO2_MW_G_MOL = 44.01
MOLAR_VOLUME_STP_L_MOL = 22.414

CO2_SOLUBILITY_PARAMETERS = ("kCO2_release_h", "CO2sat_scale")
O2_GATED_PARAMETERS = ("O2_qmax_mg_gdw_h",)
CO2_PARAMETER_BOUNDS = {
    "kCO2_release_h": (0.03, 25.0),
    "CO2sat_scale": (0.35, 2.50),
    "O2_qmax_mg_gdw_h": (0.15, 6.0),
}
CO2_DEFAULTS = {
    "kCO2_release_h": 1.2,
    "CO2sat_scale": 1.0,
}
O2_MACRO_DEFAULTS = {
    # Effective macro-oxygen defaults. Units use X in gDW/L, numerically equal
    # to kg/m3 in the existing model.
    "O2_qmax_mg_gdw_h": 1.2,
    "O2_K_mg_l": 0.25,
    "O2_ana_K_mg_l": 0.75,
    "O2_ana_hill": 2.0,
    "O2_crabtree_floor": 0.08,
    "O2sat_scale": 1.0,
    "O2_kLa_h": 0.0,
    "O2_25171_fraction": 0.05,
}
O2_SLOW_TRANSITION_DEFAULTS = {**O2_MACRO_DEFAULTS, "O2_qmax_mg_gdw_h": 0.15}
CO2_G_PER_MG_O2_RESP = (44.01 / 32.00) / 1000.0
CO2_INFORMATION_PARAMETERS = CO2_SOLUBILITY_PARAMETERS + O2_GATED_PARAMETERS

PARAMETER_BOUNDS = {
    **pilot.PARAMETER_BOUNDS,
    **v2.V2_BOUNDS,
    **aroma_sel.PARAMETER_BOUNDS,
    **CO2_PARAMETER_BOUNDS,
}


@dataclass(frozen=True)
class CO2Candidate:
    name: str
    parameters: tuple[str, ...]
    defaults: dict[str, float]
    bounds: dict[str, tuple[float, float]]
    mechanistic_class: str
    description: str


def co2_candidate_library() -> dict[str, CO2Candidate]:
    return {
        "instant": CO2Candidate(
            name="instant",
            parameters=tuple(),
            defaults={},
            bounds={},
            mechanistic_class="rate_proxy",
            description="CO2 sensor flow is proportional to instantaneous ethanol-derived CO2 production.",
        ),
        "old_lag_threshold": CO2Candidate(
            name="old_lag_threshold",
            parameters=("qgas_tau_h", "qgas_threshold_fraction"),
            defaults={"qgas_tau_h": 24.0, "qgas_threshold_fraction": 0.20},
            bounds={"qgas_tau_h": (0.25, 96.0), "qgas_threshold_fraction": (0.0, 0.90)},
            mechanistic_class="empirical_effective",
            description="Previous effective lag-threshold model used as a benchmark.",
        ),
        "solubility_fixed": CO2Candidate(
            name="solubility_fixed",
            parameters=("kCO2_release_h",),
            defaults={"kCO2_release_h": 1.2, "CO2sat_scale": 1.0},
            bounds={"kCO2_release_h": CO2_PARAMETER_BOUNDS["kCO2_release_h"]},
            mechanistic_class="dissolved_co2",
            description="Dissolved CO2 buffer with literature saturation correlation and fitted release rate.",
        ),
        "solubility_scaled": CO2Candidate(
            name="solubility_scaled",
            parameters=CO2_SOLUBILITY_PARAMETERS,
            defaults=dict(CO2_DEFAULTS),
            bounds=dict(CO2_PARAMETER_BOUNDS),
            mechanistic_class="dissolved_co2",
            description="Dissolved CO2 buffer with fitted release rate and saturation scale.",
        ),
        "solubility_o2_literature": CO2Candidate(
            name="solubility_o2_literature",
            parameters=CO2_SOLUBILITY_PARAMETERS,
            defaults={**CO2_DEFAULTS, **O2_MACRO_DEFAULTS},
            bounds={name: CO2_PARAMETER_BOUNDS[name] for name in CO2_SOLUBILITY_PARAMETERS},
            mechanistic_class="o2_gated_dissolved_co2",
            description=(
                "Dissolved CO2 buffer with a literature-informed macro-O2 state, Crabtree floor, "
                "and O2-gated anaerobic ethanol/CO2 production."
            ),
        ),
        "solubility_o2_slow_transition": CO2Candidate(
            name="solubility_o2_slow_transition",
            parameters=CO2_SOLUBILITY_PARAMETERS,
            defaults={**CO2_DEFAULTS, **O2_SLOW_TRANSITION_DEFAULTS},
            bounds={name: CO2_PARAMETER_BOUNDS[name] for name in CO2_SOLUBILITY_PARAMETERS},
            mechanistic_class="o2_gated_dissolved_co2",
            description=(
                "O2-gated dissolved CO2 model with the effective O2 uptake fixed at the slow-transition "
                "value indicated by the qfit benchmark; used to avoid estimating a bound-active O2 parameter."
            ),
        ),
        "solubility_o2_qfit": CO2Candidate(
            name="solubility_o2_qfit",
            parameters=CO2_SOLUBILITY_PARAMETERS + O2_GATED_PARAMETERS,
            defaults={**CO2_DEFAULTS, **O2_MACRO_DEFAULTS},
            bounds={name: CO2_PARAMETER_BOUNDS[name] for name in CO2_SOLUBILITY_PARAMETERS + O2_GATED_PARAMETERS},
            mechanistic_class="o2_gated_dissolved_co2",
            description=(
                "Same O2-gated dissolved CO2 model, but the effective O2 uptake scale is fitted "
                "to test whether the fixed literature value is adequate."
            ),
        ),
    }


def parameter_bounds(name: str) -> tuple[float, float]:
    if name in PARAMETER_BOUNDS:
        return PARAMETER_BOUNDS[name]
    if name == "qgas_tau_h":
        return 0.25, 96.0
    if name == "qgas_threshold_fraction":
        return 0.0, 0.90
    raise KeyError(name)


def clip_theta(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in PARAMETER_BOUNDS.items():
        if name in out and np.isfinite(float(out[name])):
            out[name] = float(np.clip(float(out[name]), lb, ub))
    for name in ("qgas_tau_h", "qgas_threshold_fraction"):
        if name in out and np.isfinite(float(out[name])):
            lb, ub = parameter_bounds(name)
            out[name] = float(np.clip(float(out[name]), lb, ub))
    for name, value in CO2_DEFAULTS.items():
        out.setdefault(name, value)
    for name, value in O2_MACRO_DEFAULTS.items():
        out.setdefault(name, value)
    return out


def load_selected_context() -> tuple[dict[str, float], str, aroma_sel.AromaVariant, str, global_sel.SecondaryVariant]:
    theta = pilot.load_initial_theta()
    theta.update(v2.default_theta_v2(theta))
    for path in [
        SCRIPT_DIR / "results" / "global_sugar" / "theta_selected_global_model.csv",
        SCRIPT_DIR / "results" / "global_state_model_selection_doe" / "theta_selected_global_model.csv",
        aroma_sel.RESULTS_DIR / "theta_selected_aroma_model.csv",
        pilot.RESULTS_DIR / "theta_pilot_extended.csv",
    ]:
        if path.exists():
            loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
            theta.update({str(k): float(v) for k, v in loaded.items() if np.isfinite(float(v))})
            break

    aroma_name = "ea_ethanol_nlimited"
    selected_aroma_path = aroma_sel.RESULTS_DIR / "selected_model.txt"
    if selected_aroma_path.exists():
        aroma_name = selected_aroma_path.read_text(encoding="utf-8").strip()
    aroma_variants = aroma_sel.variant_library()
    aroma_variant = aroma_variants.get(aroma_name, aroma_variants["ea_ethanol_nlimited"])
    for name, value in aroma_sel.EXTRA_DEFAULTS.items():
        theta.setdefault(name, value)

    secondary_name = "secondary_phase_split"
    selected_secondary_path = SCRIPT_DIR / "results" / "global_state_model_selection_doe" / "selected_secondary_model.txt"
    if selected_secondary_path.exists():
        secondary_name = selected_secondary_path.read_text(encoding="utf-8").strip()
    secondary_variants = global_sel.secondary_variant_library()
    secondary_variant = secondary_variants.get(secondary_name, secondary_variants["secondary_phase_split"])
    theta.update(secondary_variant.fixed_overrides)
    for name, value in CO2_DEFAULTS.items():
        theta.setdefault(name, value)
    return clip_theta(theta), aroma_name, aroma_variant, secondary_name, secondary_variant


def candidate_params_from_theta(theta: dict[str, float], candidate: CO2Candidate) -> dict[str, float]:
    params = dict(candidate.defaults)
    for name in candidate.parameters:
        params[name] = float(theta.get(name, candidate.defaults.get(name, CO2_DEFAULTS.get(name, O2_MACRO_DEFAULTS.get(name, 1.0)))))
    if "CO2sat_scale" not in params:
        params["CO2sat_scale"] = float(theta.get("CO2sat_scale", 1.0))
    if "kCO2_release_h" not in params:
        params["kCO2_release_h"] = float(theta.get("kCO2_release_h", 1.2))
    for name, value in O2_MACRO_DEFAULTS.items():
        params.setdefault(name, float(theta.get(name, value)))
    return params


def smooth_positive(value: float, smooth: float = 0.015) -> float:
    value = float(value)
    smooth = max(float(smooth), 1e-9)
    if value / smooth > 50.0:
        return value
    if value / smooth < -50.0:
        return 0.0
    return smooth * math.log1p(math.exp(value / smooth))


def co2_saturation_g_l(temp_c: float, ethanol_g_l: float, glucose_g_l: float, fructose_g_l: float, sat_scale: float = 1.0) -> float:
    """Working CO2 saturation correlation for wine-like water-ethanol-sugar media.

    The base value is CO2 solubility in water near 20 C and 1 atm CO2.
    Ethanol is treated as increasing apparent CO2 solubility, while total sugar
    reduces water activity. The fitted scale absorbs matrix uncertainty.
    """
    temp_factor = math.exp(-0.032 * (float(temp_c) - 20.0))
    ethanol_factor = math.exp(0.0016 * max(float(ethanol_g_l), 0.0))
    sugar_factor = math.exp(-0.0012 * max(float(glucose_g_l) + float(fructose_g_l), 0.0))
    return max(float(sat_scale) * 1.69 * temp_factor * ethanol_factor * sugar_factor, 1e-6)


def candidate_uses_o2_gate(candidate: CO2Candidate | None) -> bool:
    return bool(candidate is not None and candidate.mechanistic_class == "o2_gated_dissolved_co2")


def o2_saturation_mg_l(temp_c: float, ethanol_g_l: float, glucose_g_l: float, fructose_g_l: float, sat_scale: float = 1.0) -> float:
    """Working air-saturated dissolved oxygen concentration in wine-like must."""
    temp_factor = math.exp(-0.024 * (float(temp_c) - 20.0))
    ethanol_factor = math.exp(-0.0020 * max(float(ethanol_g_l), 0.0))
    sugar_factor = math.exp(-0.0010 * max(float(glucose_g_l) + float(fructose_g_l), 0.0))
    return max(float(sat_scale) * 8.6 * temp_factor * ethanol_factor * sugar_factor, 1e-6)


def initial_macro_o2_mg_l(batch: base.BatchData | base.FutureDesign, core: pd.DataFrame, params: dict[str, float]) -> float:
    t0 = float(core.index.min())
    temp = base.temperature_at(batch, t0)
    g = float(core.loc[t0, "G"]) if t0 in core.index else float(batch.initials.get("G", 0.0))
    f = float(core.loc[t0, "F"]) if t0 in core.index else float(batch.initials.get("F", 0.0))
    e = float(core.loc[t0, "E"]) if t0 in core.index else float(batch.initials.get("E", 0.0))
    o2 = o2_saturation_mg_l(temp, e, g, f, params.get("O2sat_scale", 1.0))
    if str(getattr(batch, "batch", "")) == "25171":
        o2 *= float(params.get("O2_25171_fraction", O2_MACRO_DEFAULTS["O2_25171_fraction"]))
    return max(float(o2), 0.0)


def integrate_macro_o2(
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    params: dict[str, float],
) -> pd.DataFrame:
    time = core.index.to_numpy(dtype=float)
    qmax = max(float(params.get("O2_qmax_mg_gdw_h", O2_MACRO_DEFAULTS["O2_qmax_mg_gdw_h"])), 1e-12)
    k_o2 = max(float(params.get("O2_K_mg_l", O2_MACRO_DEFAULTS["O2_K_mg_l"])), 1e-12)
    k_la = max(float(params.get("O2_kLa_h", O2_MACRO_DEFAULTS["O2_kLa_h"])), 0.0)
    c = initial_macro_o2_mg_l(batch, core, params)
    rows = []
    for idx, t in enumerate(time):
        if idx > 0:
            t_prev = float(time[idx - 1])
            dt = max(float(t - t_prev), 1e-9)
            x_prev = max(float(np.interp(t_prev, core.index.to_numpy(dtype=float), core["X"].to_numpy(dtype=float))), 0.0)
            temp_prev = base.temperature_at(batch, t_prev)
            g_prev = float(np.interp(t_prev, core.index.to_numpy(dtype=float), core["G"].to_numpy(dtype=float)))
            f_prev = float(np.interp(t_prev, core.index.to_numpy(dtype=float), core["F"].to_numpy(dtype=float)))
            e_prev = float(np.interp(t_prev, core.index.to_numpy(dtype=float), core["E"].to_numpy(dtype=float)))
            sat_prev = o2_saturation_mg_l(temp_prev, e_prev, g_prev, f_prev, params.get("O2sat_scale", 1.0))
            uptake_prev = qmax * x_prev * c / (k_o2 + c)
            c = max(c + dt * (k_la * (sat_prev - c) - uptake_prev), 0.0)
        x = max(float(np.interp(float(t), core.index.to_numpy(dtype=float), core["X"].to_numpy(dtype=float))), 0.0)
        uptake = qmax * x * c / (k_o2 + c)
        sat = o2_saturation_mg_l(
            base.temperature_at(batch, float(t)),
            float(np.interp(float(t), core.index.to_numpy(dtype=float), core["E"].to_numpy(dtype=float))),
            float(np.interp(float(t), core.index.to_numpy(dtype=float), core["G"].to_numpy(dtype=float))),
            float(np.interp(float(t), core.index.to_numpy(dtype=float), core["F"].to_numpy(dtype=float))),
            params.get("O2sat_scale", 1.0),
        )
        rows.append({"time_h": float(t), "O2_mg_l": float(c), "O2_saturation_mg_l": float(sat), "O2_uptake_mg_l_h": float(uptake)})
    return pd.DataFrame(rows).set_index("time_h")


def anaerobic_gate_from_o2(o2_mg_l: np.ndarray | float, params: dict[str, float]) -> np.ndarray:
    o2 = np.maximum(np.asarray(o2_mg_l, dtype=float), 0.0)
    k_ana = max(float(params.get("O2_ana_K_mg_l", O2_MACRO_DEFAULTS["O2_ana_K_mg_l"])), 1e-12)
    hill = max(float(params.get("O2_ana_hill", O2_MACRO_DEFAULTS["O2_ana_hill"])), 0.25)
    return (k_ana**hill) / (k_ana**hill + o2**hill)


def co2_production_g_l_h(
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    time: np.ndarray,
    candidate: CO2Candidate | None = None,
    params: dict[str, float] | None = None,
) -> np.ndarray:
    values = []
    for t in np.asarray(time, dtype=float):
        rates = joint._core_rates(theta, batch, core, float(t))
        values.append(joint.CO2_G_PER_G_ETHANOL * max(float(rates["ethanol_prod"]), 0.0))
    base_prod = np.asarray(values, dtype=float)
    if candidate is None or candidate.mechanistic_class != "o2_gated_dissolved_co2":
        return base_prod
    params = candidate_params_from_theta(theta, candidate) if params is None else params
    o2 = integrate_macro_o2(theta, batch, core, params)
    o2_at_time = np.interp(np.asarray(time, dtype=float), o2.index.to_numpy(dtype=float), o2["O2_mg_l"].to_numpy(dtype=float))
    uptake_at_time = np.interp(
        np.asarray(time, dtype=float),
        o2.index.to_numpy(dtype=float),
        o2["O2_uptake_mg_l_h"].to_numpy(dtype=float),
    )
    phi_ana = anaerobic_gate_from_o2(o2_at_time, params)
    floor = float(np.clip(params.get("O2_crabtree_floor", O2_MACRO_DEFAULTS["O2_crabtree_floor"]), 0.0, 1.0))
    ferment_fraction = floor + (1.0 - floor) * phi_ana
    respiratory_co2 = CO2_G_PER_MG_O2_RESP * np.maximum(uptake_at_time, 0.0)
    return np.maximum(base_prod * ferment_fraction + respiratory_co2, 0.0)


def integrate_dissolved_co2(
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    params: dict[str, float],
    candidate: CO2Candidate | None = None,
) -> pd.DataFrame:
    time = core.index.to_numpy(dtype=float)
    q_prod = co2_production_g_l_h(theta, batch, core, time, candidate, params)
    k_release = max(float(params.get("kCO2_release_h", CO2_DEFAULTS["kCO2_release_h"])), 1e-9)
    sat_scale = max(float(params.get("CO2sat_scale", CO2_DEFAULTS["CO2sat_scale"])), 1e-6)
    c = 0.0
    rows = []
    for idx, t in enumerate(time):
        rates = joint._core_rates(theta, batch, core, float(t))
        csat = co2_saturation_g_l(rates["TempC"], rates["E"], rates["G"], rates["F"], sat_scale)
        if idx == 0:
            qgas = 0.0
        else:
            dt = max(float(time[idx] - time[idx - 1]), 1e-9)
            excess = smooth_positive(c - csat)
            qgas = min(k_release * excess, c / dt + max(float(q_prod[idx - 1]), 0.0))
            c = max(c + dt * (max(float(q_prod[idx - 1]), 0.0) - qgas), 0.0)
        rows.append(
            {
                "time_h": float(t),
                "qprod_g_l_h": float(q_prod[idx]),
                "CO2_dissolved_g_l": float(c),
                "CO2_saturation_g_l": float(csat),
                "qgas_g_l_h": float(max(qgas, 0.0)),
                "qgas_L_min_per_L": float(max(qgas, 0.0) / CO2_MW_G_MOL * MOLAR_VOLUME_STP_L_MOL / 60.0),
            }
        )
    return pd.DataFrame(rows).set_index("time_h")


def gas_flow_series(
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    candidate: CO2Candidate,
    params: dict[str, float],
) -> pd.Series:
    time = core.index.to_numpy(dtype=float)
    if candidate.name == "instant":
        qgas = co2_production_g_l_h(theta, batch, core, time)
    elif candidate.name == "old_lag_threshold":
        old_mode = old_co2.co2_mode_library()["lag_threshold"]
        old_params = {
            "tau_h": float(params.get("qgas_tau_h", params.get("tau_h", 24.0))),
            "threshold_fraction": float(params.get("qgas_threshold_fraction", params.get("threshold_fraction", 0.20))),
        }
        qprod = co2_production_g_l_h(theta, batch, core, time)
        qgas = old_co2.gas_flow_from_production(time, qprod, old_mode, old_params)
    else:
        qgas = integrate_dissolved_co2(theta, batch, core, params, candidate)["qgas_g_l_h"].to_numpy(dtype=float)
    return pd.Series(np.maximum(qgas, 0.0), index=time)


def co2_prediction_for_batch(
    theta: dict[str, float],
    batch: base.BatchData,
    sensor_times: np.ndarray,
    candidate: CO2Candidate,
    params: dict[str, float],
    dt_grid_h: float = 0.25,
) -> np.ndarray:
    horizon = max(float(np.nanmax(batch.time)), float(np.nanmax(sensor_times)))
    grid = np.arange(0.0, horizon + dt_grid_h, dt_grid_h)
    sim_times = np.asarray(sorted(set(np.round(grid, 8)).union(set(np.round(sensor_times, 8))).union({0.0, horizon})), dtype=float)
    core = base.simulate(batch, theta, sim_times)
    if core is None:
        return np.full_like(sensor_times, np.nan, dtype=float)
    qgas = gas_flow_series(theta, batch, core, candidate, params)
    return np.interp(np.asarray(sensor_times, dtype=float), qgas.index.to_numpy(dtype=float), qgas.to_numpy(dtype=float))


def co2_residual(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    candidate: CO2Candidate,
    params: dict[str, float],
) -> np.ndarray:
    residuals = []
    if co2_down.empty:
        return np.array([], dtype=float)
    for batch_name, group in co2_down.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        obs = group["co2_raw"].to_numpy(dtype=float)
        pred = co2_prediction_for_batch(theta, batch, times, candidate, params)
        mask = np.isfinite(obs) & np.isfinite(pred) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        obs = obs[mask]
        pred = pred[mask]
        denom = float(np.dot(pred, pred))
        scale = float(np.dot(obs, pred) / denom) if denom > 1e-12 else 0.0
        sigma = max(pilot.CO2_RATE_SIGMA, 0.10 * float(np.nanmax(obs)))
        residuals.append((scale * pred - obs) / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def co2_residual_from_cache(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    core_cache: dict[str, pd.DataFrame],
    co2_down: pd.DataFrame,
    candidate: CO2Candidate,
    params: dict[str, float],
) -> np.ndarray:
    residuals = []
    if co2_down.empty:
        return np.array([], dtype=float)
    for batch_name, group in co2_down.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        core = core_cache.get(batch.label)
        if core is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        obs = group["co2_raw"].to_numpy(dtype=float)
        qgas = gas_flow_series(theta, batch, core, candidate, params)
        pred = np.interp(times, qgas.index.to_numpy(dtype=float), qgas.to_numpy(dtype=float))
        mask = np.isfinite(obs) & np.isfinite(pred) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        obs = obs[mask]
        pred = pred[mask]
        denom = float(np.dot(pred, pred))
        scale = float(np.dot(obs, pred) / denom) if denom > 1e-12 else 0.0
        sigma = max(pilot.CO2_RATE_SIGMA, 0.10 * float(np.nanmax(obs)))
        residuals.append((scale * pred - obs) / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def co2_metrics(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    candidate: CO2Candidate,
    params: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    if co2_down.empty:
        return pd.DataFrame(), pd.DataFrame()
    for batch_name, group in co2_down.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        obs = group["co2_raw"].to_numpy(dtype=float)
        raw_pred = co2_prediction_for_batch(theta, batch, times, candidate, params)
        mask = np.isfinite(obs) & np.isfinite(raw_pred) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        obs = obs[mask]
        raw_pred = raw_pred[mask]
        times = times[mask]
        denom = float(np.dot(raw_pred, raw_pred))
        scale = float(np.dot(obs, raw_pred) / denom) if denom > 1e-12 else 0.0
        pred = scale * raw_pred
        err = pred - obs
        corr = np.nan
        if np.nanstd(obs) > 1e-12 and np.nanstd(pred) > 1e-12:
            corr = float(np.corrcoef(obs, pred)[0, 1])
        obs_scale = max(float(np.nanmedian(np.abs(obs))), 0.25 * float(np.nanmax(obs) - np.nanmin(obs)), 1e-6)
        rows.append(
            {
                "mode": candidate.name,
                "batch": str(batch_name),
                "n": int(len(obs)),
                "scale_factor_L_min_per_g_l_h": scale,
                "rmse_L_min": float(np.sqrt(np.mean(err * err))),
                "mae_L_min": float(np.mean(np.abs(err))),
                "bias_L_min": float(np.mean(err)),
                "relative_rmse": float(np.sqrt(np.mean(err * err)) / obs_scale),
                "relative_bias": float(np.mean(err) / obs_scale),
                "corr": corr,
            }
        )
        for t, y, yhat, raw in zip(times, obs, pred, raw_pred):
            pred_rows.append(
                {
                    "mode": candidate.name,
                    "batch": str(batch_name),
                    "time_h": float(t),
                    "state": "CO2_flow_L_min",
                    "obs": float(y),
                    "pred": float(yhat),
                    "qgas_g_l_h_raw": float(raw),
                    "scale_factor": scale,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def fit_co2_candidate(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    candidate: CO2Candidate,
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    core_cache = build_core_cache(theta, list(batches_by_name.values()), extra_times=co2_extra_times(co2_down))
    if not candidate.parameters:
        params = dict(candidate.defaults)
        res = co2_residual_from_cache(theta, batches_by_name, core_cache, co2_down, candidate, params)
        return params, pd.DataFrame(
            [
                {
                    "mode": candidate.name,
                    "start": 0,
                    "success": True,
                    "nfev": 0,
                    "data_wsse": float(np.dot(res, res)),
                    "n_data_residuals": int(len(res)),
                    "n_parameters": 0,
                    "active_bound_count": 0,
                    "params_json": "{}",
                }
            ]
        )

    lb = np.asarray([candidate.bounds[name][0] for name in candidate.parameters], dtype=float)
    ub = np.asarray([candidate.bounds[name][1] for name in candidate.parameters], dtype=float)
    x_ref = np.asarray([candidate.defaults[name] for name in candidate.parameters], dtype=float)
    log_mask = np.asarray([name != "qgas_threshold_fraction" for name in candidate.parameters], dtype=bool)
    x0_ref = x_ref.copy()
    x0_ref[log_mask] = np.log(np.maximum(x0_ref[log_mask], 1e-16))
    lb_vec = lb.copy()
    ub_vec = ub.copy()
    lb_vec[log_mask] = np.log(lb_vec[log_mask])
    ub_vec[log_mask] = np.log(ub_vec[log_mask])
    rng = np.random.default_rng(seed)
    best_score = np.inf
    best_params = dict(candidate.defaults)
    rows = []

    def make_params(x: np.ndarray) -> dict[str, float]:
        out = dict(candidate.defaults)
        for name, value, is_log in zip(candidate.parameters, np.asarray(x, dtype=float), log_mask):
            out[name] = float(math.exp(float(value))) if bool(is_log) else float(value)
        return out

    def objective(x: np.ndarray) -> np.ndarray:
        params = make_params(x)
        res = [co2_residual_from_cache(theta, batches_by_name, core_cache, co2_down, candidate, params)]
        if "CO2sat_scale" in candidate.parameters:
            res.append(np.array([math.log(max(params["CO2sat_scale"], 1e-16)) / 0.45], dtype=float))
        if "O2_qmax_mg_gdw_h" in candidate.parameters:
            ref = max(float(O2_MACRO_DEFAULTS["O2_qmax_mg_gdw_h"]), 1e-16)
            res.append(np.array([math.log(max(params["O2_qmax_mg_gdw_h"], 1e-16) / ref) / 0.70], dtype=float))
        return np.concatenate([r for r in res if len(r)])

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = np.clip(x0_ref, lb_vec + 1e-9, ub_vec - 1e-9)
        else:
            x0 = np.empty_like(x0_ref)
            for pos, name in enumerate(candidate.parameters):
                if log_mask[pos]:
                    x0[pos] = rng.uniform(lb_vec[pos], ub_vec[pos])
                else:
                    x0[pos] = rng.uniform(lb_vec[pos], ub_vec[pos])
        result = least_squares(
            objective,
            x0,
            bounds=(lb_vec, ub_vec),
            method="trf",
            x_scale="jac",
            loss="soft_l1",
            f_scale=2.0,
            max_nfev=int(max_nfev),
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
        )
        params = make_params(result.x)
        data_res = co2_residual_from_cache(theta, batches_by_name, core_cache, co2_down, candidate, params)
        data_wsse = float(np.dot(data_res, data_res))
        active = 0
        for name in candidate.parameters:
            lo, hi = candidate.bounds[name]
            value = float(params[name])
            if value <= lo * 1.01 or value >= hi / 1.01:
                active += 1
        score = data_wsse + 10.0 * active
        row = {
            "mode": candidate.name,
            "start": idx,
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "data_wsse": data_wsse,
            "n_data_residuals": int(len(data_res)),
            "n_parameters": int(len(candidate.parameters)),
            "active_bound_count": int(active),
            "fit_selection_score": score,
            "params_json": json.dumps(params),
        }
        row.update(params)
        rows.append(row)
        if score < best_score:
            best_score = score
            best_params = params
    return best_params, pd.DataFrame(rows)


def information_criteria(data_wsse: float, n: int, k: int) -> tuple[float, float]:
    n = max(int(n), 1)
    k = max(int(k), 1)
    aic = float(data_wsse + 2 * k)
    if n > k + 1:
        aic += float((2 * k * (k + 1)) / (n - k - 1))
    bic = float(data_wsse + k * math.log(n))
    return aic, bic


def build_core_cache(theta: dict[str, float], batches: list[base.BatchData], extra_times: dict[str, np.ndarray] | None = None) -> dict[str, pd.DataFrame]:
    cache = {}
    for batch in batches:
        horizon = float(np.nanmax(batch.time))
        grid = set(np.round(np.arange(0.0, horizon + 0.5, 0.5), 8))
        grid.update(np.round(batch.time, 8))
        if extra_times and batch.batch in extra_times:
            grid.update(np.round(extra_times[batch.batch], 8))
        times = np.asarray(sorted(float(t) for t in grid if float(t) >= 0.0), dtype=float)
        core = base.simulate(batch, theta, times)
        if core is not None:
            cache[batch.label] = core.sort_index()
    return cache


def co2_extra_times(co2_down: pd.DataFrame) -> dict[str, np.ndarray]:
    if co2_down.empty:
        return {}
    return {
        str(batch): group["time_h_effective"].to_numpy(dtype=float)
        for batch, group in co2_down.groupby("batch", sort=True)
    }


def qgas_cache_for_batches(
    theta: dict[str, float],
    batches: list[base.BatchData],
    core_cache: dict[str, pd.DataFrame],
    candidate: CO2Candidate,
    params: dict[str, float],
) -> dict[str, pd.Series]:
    out = {}
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            continue
        out[batch.label] = gas_flow_series(theta, batch, core, candidate, params)
    return out


def aroma_prediction_rows(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    qgas_cache: dict[str, pd.Series],
    candidate_name: str,
) -> pd.DataFrame:
    return old_co2.aroma_prediction_rows_combo(
        theta,
        batches,
        variant,
        core_cache,
        sec_cache,
        qgas_cache,
        PARTITION_MODE,
        candidate_name,
    )


def aroma_residual(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    qgas_cache: dict[str, pd.Series],
) -> np.ndarray:
    return old_co2.aroma_residual_combo(theta, batches, variant, core_cache, sec_cache, qgas_cache, PARTITION_MODE)


def summarize_prediction_rows(pred: pd.DataFrame, group_cols: list[str], group_name: str, floor: float = 1e-9) -> pd.DataFrame:
    rows = []
    if pred.empty:
        return pd.DataFrame()
    for keys, group in pred.groupby(group_cols, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        obs = group["obs"].to_numpy(dtype=float)
        yhat = group["pred"].to_numpy(dtype=float)
        err = yhat - obs
        obs_clean = obs[np.isfinite(obs)]
        if obs_clean.size:
            scale = max(float(np.nanmedian(np.abs(obs_clean))), 0.25 * float(np.nanmax(obs_clean) - np.nanmin(obs_clean)), floor)
        else:
            scale = floor
        corr = np.nan
        if len(group) >= 4 and np.nanstd(obs) > 1e-12 and np.nanstd(yhat) > 1e-12:
            corr = float(np.corrcoef(obs, yhat)[0, 1])
        row = {
            "group": group_name,
            "n": int(len(group)),
            "rmse": float(np.sqrt(np.nanmean(err * err))),
            "mae": float(np.nanmean(np.abs(err))),
            "bias": float(np.nanmean(err)),
            "obs_scale": float(scale),
            "relative_rmse": float(np.sqrt(np.nanmean(err * err)) / scale),
            "relative_bias": float(np.nanmean(err) / scale),
            "corr": corr,
        }
        for col, key in zip(group_cols, keys):
            row[col] = key
        rows.append(row)
    return pd.DataFrame(rows)


def core_prediction_rows(theta: dict[str, float], batches: list[base.BatchData]) -> pd.DataFrame:
    return global_sel.core_prediction_rows(theta, batches)


def secondary_prediction_rows(theta: dict[str, float], batches: list[base.BatchData], core_cache: dict[str, pd.DataFrame]) -> pd.DataFrame:
    return global_sel.secondary_prediction_rows(theta, batches, core_cache)


def primary_residual_from_cache(primary_batches: list[base.BatchData], core_cache: dict[str, pd.DataFrame]) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in primary_batches:
        sim = core_cache.get(batch.label)
        if sim is None:
            return np.ones(1000, dtype=float) * 1e6
        for state in base.STATE_NAMES:
            obs = np.asarray(batch.observations[state], dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sim.loc[batch.time, state].to_numpy(dtype=float)
            sigma = base.sigma_for_state(state, obs[mask])
            residuals.append((pred[mask] - obs[mask]) / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def all_state_predictions(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    variant: aroma_sel.AromaVariant,
    candidate: CO2Candidate,
    params: dict[str, float],
    label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    core_cache = build_core_cache(theta, extended_batches)
    sec_cache = aroma_sel.secondary_cache_for(extended_batches, theta, core_cache)
    qgas_cache = qgas_cache_for_batches(theta, extended_batches, core_cache, candidate, params)
    frames = []
    metrics = []

    core_pred = core_prediction_rows(theta, primary_batches)
    if not core_pred.empty:
        core_pred["model_label"] = label
        frames.append(core_pred.assign(group="core"))
        metrics.append(summarize_prediction_rows(core_pred, ["state"], "core"))

    sec_pred = secondary_prediction_rows(theta, extended_batches, core_cache)
    if not sec_pred.empty:
        sec_pred["model_label"] = label
        frames.append(sec_pred.assign(group="secondary"))
        metrics.append(summarize_prediction_rows(sec_pred, ["state"], "secondary"))

    aroma_pred = aroma_prediction_rows(theta, extended_batches, variant, core_cache, sec_cache, qgas_cache, candidate.name)
    if not aroma_pred.empty:
        aroma_pred = aroma_pred.rename(columns={"pool": "state_pool"})
        aroma_pred["state"] = aroma_pred["species"].astype(str) + ":" + aroma_pred["state_pool"].astype(str)
        aroma_pred["model_label"] = label
        frames.append(aroma_pred.assign(group="aroma"))
        metrics.append(summarize_prediction_rows(aroma_pred, ["species", "state_pool"], "aroma", floor=1e-3))

    co2_met, co2_pred = co2_metrics(theta, {b.batch: b for b in extended_batches}, co2_down, candidate, params)
    if not co2_pred.empty:
        co2_pred["model_label"] = label
        frames.append(co2_pred.assign(group="co2"))
    if not co2_met.empty:
        tmp = co2_met.rename(columns={"rmse_L_min": "rmse", "mae_L_min": "mae", "bias_L_min": "bias"})
        tmp["group"] = "co2"
        tmp["state"] = "CO2_flow_L_min"
        metrics.append(tmp)

    pred_all = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    metric_all = pd.concat(metrics, ignore_index=True, sort=False) if metrics else pd.DataFrame()
    return pred_all, metric_all


def fit_integrated_model(
    theta0: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    variant: aroma_sel.AromaVariant,
    candidate: CO2Candidate,
    start_params: dict[str, float],
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    # CO2 parameters are calibrated in the step-0 gas-flow benchmark and then
    # held fixed during aroma recalibration. Letting aroma residuals move the
    # solubility scale can create non-physical compensation against stripping.
    fit_parameters = tuple(dict.fromkeys(variant.parameters))
    if not fit_parameters:
        return dict(theta0), pd.DataFrame()

    theta_ref = dict(theta0)
    theta_ref.update(start_params)
    lb = []
    ub = []
    x_ref = []
    direct_mask = []
    for name in fit_parameters:
        lo, hi = parameter_bounds(name)
        direct = name == "qgas_threshold_fraction"
        direct_mask.append(direct)
        if direct:
            lb.append(float(lo))
            ub.append(float(hi))
            x_ref.append(float(np.clip(float(theta_ref[name]), lo, hi)))
        else:
            lb.append(math.log(float(lo)))
            ub.append(math.log(float(hi)))
            x_ref.append(math.log(max(float(theta_ref[name]), float(lo) * 1.001)))
    lb = np.asarray(lb, dtype=float)
    ub = np.asarray(ub, dtype=float)
    x_ref = np.clip(np.asarray(x_ref, dtype=float), lb + 1e-9, ub - 1e-9)
    rng = np.random.default_rng(seed)
    rows = []
    best_score = np.inf
    best_theta = dict(theta_ref)
    fixed_core_cache = build_core_cache(theta_ref, extended_batches, extra_times=co2_extra_times(co2_down))
    fixed_sec_cache = aroma_sel.secondary_cache_for(extended_batches, theta_ref, fixed_core_cache)
    batches_by_name = {b.batch: b for b in extended_batches}

    def theta_from_x(x: np.ndarray) -> dict[str, float]:
        theta = dict(theta_ref)
        for name, value, direct in zip(fit_parameters, np.asarray(x, dtype=float), direct_mask):
            theta[name] = float(value) if direct else float(math.exp(float(value)))
        return clip_theta(theta)

    def data_residual(theta: dict[str, float]) -> np.ndarray:
        params = candidate_params_from_theta(theta, candidate)
        qgas_cache = qgas_cache_for_batches(theta, extended_batches, fixed_core_cache, candidate, params)
        return np.concatenate(
            [
                co2_residual_from_cache(theta, batches_by_name, fixed_core_cache, co2_down, candidate, params),
                aroma_residual(theta, extended_batches, variant, fixed_core_cache, fixed_sec_cache, qgas_cache),
            ]
        )

    def objective(x: np.ndarray) -> np.ndarray:
        theta = theta_from_x(x)
        residuals = [data_residual(theta)]
        priors = []
        for name in fit_parameters:
            ref = max(float(theta_ref[name]), 1e-16)
            scale = 0.45 if name == "CO2sat_scale" else 1.2
            if name.startswith("k_EA_") or name == "q10_EA":
                scale = 3.0
            if name.startswith("alpha_"):
                scale = 0.8
            if name == "qgas_threshold_fraction":
                priors.append((float(theta[name]) - float(theta_ref[name])) / 0.30)
            else:
                priors.append(math.log(max(float(theta[name]), 1e-16) / ref) / scale)
        residuals.append(np.asarray(priors, dtype=float))
        return np.concatenate([r for r in residuals if len(r)])

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = x_ref.copy()
        else:
            x0 = np.clip(x_ref + rng.normal(0.0, 0.75, size=len(fit_parameters)), lb + 1e-9, ub - 1e-9)
        start = objective(x0)
        result = least_squares(
            objective,
            x0,
            bounds=(lb, ub),
            method="trf",
            x_scale="jac",
            loss="soft_l1",
            f_scale=2.0,
            max_nfev=int(max_nfev),
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
        )
        theta_hat = theta_from_x(result.x)
        end = objective(result.x)
        data = data_residual(theta_hat)
        active = sum(
            1
            for name in fit_parameters
            if theta_hat[name] <= parameter_bounds(name)[0] * 1.01 or theta_hat[name] >= parameter_bounds(name)[1] / 1.01
        )
        score = float(np.dot(data, data)) + 20.0 * active
        row = {
            "mode": candidate.name,
            "start": idx,
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "initial_objective_wsse": float(np.dot(start, start)),
            "final_objective_wsse": float(np.dot(end, end)),
            "data_wsse": float(np.dot(data, data)),
            "n_data_residuals": int(len(data)),
            "n_parameters": int(len(fit_parameters)),
            "active_bound_count": int(active),
            "fit_selection_score": score,
        }
        row.update({name: float(theta_hat[name]) for name in fit_parameters})
        rows.append(row)
        if score < best_score:
            best_score = score
            best_theta = theta_hat
    return clip_theta(best_theta), pd.DataFrame(rows)


def current_residual(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    variant: aroma_sel.AromaVariant,
    candidate: CO2Candidate,
) -> np.ndarray:
    params = candidate_params_from_theta(theta, candidate)
    core_cache = build_core_cache(theta, extended_batches)
    sec_cache = aroma_sel.secondary_cache_for(extended_batches, theta, core_cache)
    qgas_cache = qgas_cache_for_batches(theta, extended_batches, core_cache, candidate, params)
    residuals = [
        primary_residual_from_cache(primary_batches, core_cache),
        v2.residual_v2(theta, extended_batches, core_cache),
        co2_residual_from_cache(theta, {b.batch: b for b in extended_batches}, core_cache, co2_down, candidate, params),
        aroma_residual(theta, extended_batches, variant, core_cache, sec_cache, qgas_cache),
    ]
    return np.concatenate([r for r in residuals if len(r)])


def finite_difference_jacobian(
    theta: dict[str, float],
    parameters: tuple[str, ...],
    residual_fun: Callable[[dict[str, float]], np.ndarray],
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    base_res = residual_fun(theta)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        lb, ub = parameter_bounds(name)
        if name == "qgas_threshold_fraction":
            delta = step * max(float(ub) - float(lb), 1e-9)
            theta_plus[name] = float(np.clip(float(theta[name]) + delta, lb, ub))
            theta_minus[name] = float(np.clip(float(theta[name]) - delta, lb, ub))
            denom = theta_plus[name] - theta_minus[name]
        else:
            theta_plus[name] = float(np.clip(float(theta[name]) * math.exp(step), lb, ub))
            theta_minus[name] = float(np.clip(float(theta[name]) * math.exp(-step), lb, ub))
            denom = math.log(theta_plus[name] / theta_minus[name]) if theta_plus[name] > 0 and theta_minus[name] > 0 else 2.0 * step
        denom = max(abs(denom), 1e-12)
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / denom)
    return np.column_stack(cols), base_res


def stable_inverse(fim: np.ndarray, ridge_fraction: float = 1e-9) -> np.ndarray:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    scale = max(float(np.trace(fim)) / max(fim.shape[0], 1), 1.0)
    return np.linalg.pinv(fim + ridge_fraction * scale * np.eye(fim.shape[0]))


def fim_metrics(fim: np.ndarray, prefix: str = "") -> dict[str, float]:
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
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace_inv": float(np.trace(cov)),
        f"{prefix}rank_1e-8": int(np.sum(eig > max_eig * 1e-8)) if np.isfinite(max_eig) and max_eig > 0 else 0,
    }


def fim_diagnostics(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...], label: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
            "relative_eigenvalue": eigvals / max_eig if np.isfinite(max_eig) and max_eig > 0 else np.nan,
        }
    )
    weak_rows = []
    for idx in range(min(10, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        top = np.argsort(np.abs(vec))[::-1][:8]
        weak_rows.append(
            {
                "analysis": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "relative_eigenvalue": float(eigvals[idx] / max_eig) if np.isfinite(max_eig) and max_eig > 0 else np.nan,
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs(vec[i]):.3f}" for i in top),
            }
        )
    cov = stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std = float(math.sqrt(max(cov[idx, idx], 0.0)))
        value = float(theta[name])
        lb, ub = parameter_bounds(name)
        active = value <= lb * 1.01 or value >= ub / 1.01
        if std <= 0.35 and not active:
            cls = "well_estimated"
        elif std <= 0.75 and not active:
            cls = "moderate"
        elif std <= 1.25 and not active:
            cls = "weak_but_actionable"
        else:
            cls = "weak_or_confounded"
        rows.append(
            {
                "analysis": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std,
                "approx_95_multiplier": float(math.exp(1.96 * min(std, 20.0))),
                "fim_diag": float(fim[idx, idx]),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return spectrum, pd.DataFrame(weak_rows), pd.DataFrame(rows)


def future_sample_times(design: base.FutureDesign) -> np.ndarray:
    return pilot.future_sample_times(design)


def future_co2_times(design: base.FutureDesign) -> np.ndarray:
    return np.arange(0.0, float(design.horizon_h) + 1e-9, 4.0)


def future_residual(
    theta: dict[str, float],
    design: base.FutureDesign,
    variant: aroma_sel.AromaVariant,
    candidate: CO2Candidate,
) -> np.ndarray:
    sample_times = future_sample_times(design)
    co2_times = future_co2_times(design)
    all_times = np.asarray(sorted(set(sample_times).union(set(co2_times)).union({0.0, float(design.horizon_h)})), dtype=float)
    core = base.simulate(design, theta, all_times)
    if core is None:
        return np.ones(1000, dtype=float) * 1e6
    secondary = v2.integrate_secondary_v2(design, theta, core)
    params = candidate_params_from_theta(theta, candidate)
    qgas = gas_flow_series(theta, design, core, candidate, params)
    try:
        liq, cond = old_co2.integrate_aroma_combo(design, theta, core, secondary, qgas, variant, PARTITION_MODE)
    except Exception:
        return np.ones(1000, dtype=float) * 1e6
    residuals: list[np.ndarray] = []
    valid_samples = [float(t) for t in sample_times if float(t) in core.index]
    if valid_samples:
        core_sample = core.loc[valid_samples]
        for state in base.STATE_NAMES:
            center = core_sample[state].to_numpy(dtype=float)
            residuals.append(center / base.sigma_for_state(state, center))
        sec_sample = secondary.loc[valid_samples]
        for state in ("Pyr", "AcAld", "Acetate", "O2"):
            if state == "O2":
                sigma = np.maximum(0.25, 0.20 * np.maximum(sec_sample[state].to_numpy(dtype=float), 0.25))
            else:
                sigma = joint.SIGMA[state]
            residuals.append(sec_sample[state].to_numpy(dtype=float) / sigma)
        aroma_sample = liq.loc[valid_samples]
        cond_sample = cond.loc[valid_samples]
        for species in joint.AROMA_SPECIES:
            liq_center = aroma_sample[species].to_numpy(dtype=float)
            cond_center = cond_sample[species].to_numpy(dtype=float)
            residuals.append(liq_center / pilot.aroma_sigma(species, "liquid", liq_center))
            residuals.append(cond_center / pilot.aroma_sigma(species, "condensate", cond_center))
    valid_co2 = [float(t) for t in co2_times if float(t) in core.index]
    if valid_co2:
        qvals = np.asarray([old_co2.interp_series(qgas, t) for t in valid_co2], dtype=float)
        sigma = np.maximum(0.03, 0.20 * np.maximum(qvals, 0.03))
        residuals.append(qvals / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def candidate_fim(theta: dict[str, float], design: base.FutureDesign, variant: aroma_sel.AromaVariant, candidate: CO2Candidate, parameters: tuple[str, ...], step: float) -> np.ndarray:
    jac, _ = finite_difference_jacobian(theta, parameters, lambda th: future_residual(th, design, variant, candidate), step)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    post_cov = stable_inverse(combined)
    rows = {}
    reductions = []
    aroma_co2_reductions = []
    for idx, name in enumerate(parameters):
        before = float(prior_cov[idx, idx])
        after = float(post_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        reduction = 1.0 - ratio if np.isfinite(ratio) else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = reduction
        if np.isfinite(reduction):
            reductions.append(reduction)
            if name.startswith("k_") or name.startswith("alpha_") or name in CO2_INFORMATION_PARAMETERS:
                aroma_co2_reductions.append(reduction)
    rows["target_mean_var_reduction"] = float(np.nanmean(reductions)) if reductions else np.nan
    rows["target_worst_var_reduction"] = float(np.nanmin(reductions)) if reductions else np.nan
    rows["aroma_co2_mean_var_reduction"] = float(np.nanmean(aroma_co2_reductions)) if aroma_co2_reductions else np.nan
    rows["aroma_co2_worst_var_reduction"] = float(np.nanmin(aroma_co2_reductions)) if aroma_co2_reductions else np.nan
    return rows


def score_fim(fim: np.ndarray, objective: str = "hybrid") -> float:
    metrics = fim_metrics(fim)
    if objective == "d_opt":
        return float(metrics["logdet"])
    if objective == "e_opt":
        return math.log(max(float(metrics["min_eigenvalue"]), 1e-18))
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-18)
    return float(metrics["logdet"]) + 2.0 * math.log(min_rel) - 0.05 * math.log(max(float(metrics["trace_inv"]), 1e-18))


def rank_candidates(candidate_fims: dict[str, np.ndarray], designs: dict[str, base.FutureDesign], prior: np.ndarray, parameters: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for name, fim in candidate_fims.items():
        combined = prior + fim
        design = designs[name]
        rows.append(
            {
                "candidate": name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                **fim_metrics(fim, prefix="candidate_"),
                **fim_metrics(combined, prefix="combined_"),
                **variance_reduction(prior, combined, parameters),
                "hybrid_score": score_fim(combined, "hybrid"),
                "dopt_score": score_fim(combined, "d_opt"),
                "eopt_score": score_fim(combined, "e_opt"),
                "rationale": design.rationale,
            }
        )
    return pd.DataFrame(rows).sort_values(["hybrid_score", "combined_logdet"], ascending=False).reset_index(drop=True)


def greedy_select(candidate_fims: dict[str, np.ndarray], designs: dict[str, base.FutureDesign], prior: np.ndarray, parameters: tuple[str, ...], campaign_size: int, objective: str) -> tuple[pd.DataFrame, np.ndarray]:
    remaining = set(candidate_fims)
    current = prior.copy()
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + candidate_fims[name]
            score = score_fim(trial, objective)
            if score > best_score:
                best = name
                best_score = score
                best_metrics = {**fim_metrics(trial, prefix="campaign_"), **variance_reduction(prior, trial, parameters)}
        if best is None:
            break
        current = current + candidate_fims[best]
        remaining.remove(best)
        design = designs[best]
        rows.append(
            {
                "objective": objective,
                "campaign_order": order,
                "candidate": best,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                "score": float(best_score),
                "rationale": design.rationale,
                **(best_metrics or {}),
            }
        )
    return pd.DataFrame(rows), current


def add_solubility_designs(model_data: pd.DataFrame, designs: dict[str, base.FutureDesign]) -> dict[str, base.FutureDesign]:
    out = dict(designs)
    initials = pilot.natural_initials(model_data)
    h = pilot.NATURAL_DESIGN_HORIZON_H
    zero = {channel: tuple() for channel in base.INPUT_CHANNELS}
    specs = [
        (
            "natural_pilot_CO2_cold_start_warm_ramp",
            "CO2_solubility_temperature",
            (12.0, 14.0, 22.0, 21.0),
            ((54.0, 0.035),),
            "Cold start increases CO2 capacity, then warm ramp releases dissolved CO2 and excites stripping dynamics.",
        ),
        (
            "natural_pilot_CO2_warm_start_cool_retention",
            "CO2_strip_retention",
            (23.0, 22.0, 16.0, 15.0),
            (),
            "Warm early fermentation maximizes gas release, then cooling separates production from retention.",
        ),
        (
            "natural_pilot_CO2_Npulse_release_probe",
            "CO2_N_temperature",
            (18.0, 18.0, 23.0, 20.0),
            ((42.0, 0.045),),
            "N pulse during active growth perturbs CO2 production while temperature step perturbs dissolved capacity.",
        ),
    ]
    for name, family, temps, n_pulses, rationale in specs:
        pulses = dict(zero)
        pulses["N"] = pilot.snap_pulses(h, n_pulses)
        out[name] = base.FutureDesign(
            name=name,
            family=family,
            medium="natural",
            horizon_h=h,
            initials=dict(initials),
            temperature_segments=tuple(float(t) for t in temps),
            pulses={channel: tuple(rows) for channel, rows in pulses.items()},
            rationale=rationale,
        )
    return out


def plot_data_overview(model_data: pd.DataFrame, co2_down: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cols = [
        ("G_g_l", "Glucose (g/L)"),
        ("F_g_l", "Fructose (g/L)"),
        ("YAN_mg_l", "YAN (mg/L)"),
        ("E_g_l", "Ethanol (g/L)"),
        ("X_viable_kg_m3", "Viable biomass (kg/m3)"),
        ("glycerol_g_l", "Glycerol (g/L)"),
        ("DO_mg_l", "DO (mg/L)"),
        ("ethyl_acetate_total", "Ethyl acetate total (mg/L)"),
        ("ethyl_acetate_condensate", "Ethyl acetate condensate eq. (mg/L)"),
        ("isoamyl_acetate_total", "Isoamyl acetate total (mg/L)"),
        ("isoamyl_acetate_condensate", "Isoamyl acetate condensate eq. (mg/L)"),
    ]
    for batch, group in model_data.groupby("batch", sort=True):
        n = len(cols)
        fig, axes = plt.subplots(math.ceil(n / 2), 2, figsize=(13, 2.3 * math.ceil(n / 2)), sharex=True)
        axes = axes.ravel()
        for ax, (col, title) in zip(axes, cols):
            if col not in group.columns:
                ax.axis("off")
                continue
            sub = group[["time_h", col]].dropna()
            if sub.empty:
                ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            else:
                ax.scatter(sub["time_h"], sub[col], s=18)
            ax.set_title(title)
            ax.grid(True, alpha=0.25)
        for ax in axes[len(cols) :]:
            ax.axis("off")
        fig.suptitle(f"Pilot data read check: {batch}", y=0.995)
        fig.tight_layout()
        fig.savefig(output_dir / f"data_read_check_{batch}.png", dpi=170)
        plt.close(fig)
    if not co2_down.empty:
        for batch, group in co2_down.groupby("batch", sort=True):
            fig, ax = plt.subplots(figsize=(10, 3.2))
            ax.plot(group["time_h_effective"], group["co2_raw"], lw=1.0)
            ax.set_title(f"Curated CO2 sensor signal: {batch}")
            ax.set_xlabel("effective time (h)")
            ax.set_ylabel("CO2 (L/min)")
            ax.grid(True, alpha=0.25)
            fig.tight_layout()
            fig.savefig(output_dir / f"co2_curated_{batch}.png", dpi=170)
            plt.close(fig)


def plot_prediction_by_batch(pred: pd.DataFrame, output_dir: Path, label: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred.empty:
        return
    states = [
        "X",
        "Xd",
        "N",
        "G",
        "F",
        "E",
        "Gly",
        "Pyr",
        "AcAld",
        "Acetate",
        "O2",
        "ethyl_acetate:retained",
        "ethyl_acetate:condensate",
        "ethyl_acetate:total",
        "isoamyl_acetate:retained",
        "isoamyl_acetate:condensate",
        "isoamyl_acetate:total",
        "ethyl_octanoate:retained",
        "ethyl_octanoate:condensate",
        "ethyl_octanoate:total",
        "CO2_flow_L_min",
    ]
    for batch, group in pred.groupby("batch", sort=True):
        present = [s for s in states if group["state"].astype(str).eq(s).any()]
        if not present:
            continue
        ncols = 3
        nrows = int(math.ceil(len(present) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(15, max(3.0, 2.7 * nrows)), sharex=False)
        axes = np.asarray(axes).ravel()
        for ax, state in zip(axes, present):
            sub = group[group["state"].astype(str).eq(state)].sort_values("time_h")
            ax.scatter(sub["time_h"], sub["obs"], s=16, color="black", label="obs")
            ax.plot(sub["time_h"], sub["pred"], color="tab:blue", lw=1.3, label="pred")
            err = sub["pred"].to_numpy(dtype=float) - sub["obs"].to_numpy(dtype=float)
            rmse = float(np.sqrt(np.nanmean(err * err))) if len(err) else np.nan
            ax.set_title(f"{state} | RMSE={rmse:.3g}", fontsize=9)
            ax.grid(True, alpha=0.25)
        for ax in axes[len(present) :]:
            ax.axis("off")
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right")
        fig.suptitle(f"{label}: fitted trajectories for batch {batch}", y=0.997)
        fig.tight_layout()
        fig.savefig(output_dir / f"{label}_fit_{batch}.png", dpi=170)
        plt.close(fig)


def plot_co2_benchmark(pred: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred.empty:
        return
    for batch, group in pred.groupby("batch", sort=True):
        fig, ax = plt.subplots(figsize=(10, 4))
        for mode, sub in group.groupby("mode", sort=True):
            sub = sub.sort_values("time_h")
            if mode == group["mode"].iloc[0]:
                ax.scatter(sub["time_h"], sub["obs"], s=12, color="black", alpha=0.35, label="obs")
            ax.plot(sub["time_h"], sub["pred"], lw=1.2, label=str(mode))
        ax.set_title(f"CO2 benchmark: {batch}")
        ax.set_xlabel("effective time (h)")
        ax.set_ylabel("CO2 (L/min)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output_dir / f"co2_benchmark_{batch}.png", dpi=170)
        plt.close(fig)


def o2_macro_diagnostics(
    theta: dict[str, float],
    batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    candidate: CO2Candidate,
    params: dict[str, float],
) -> pd.DataFrame:
    if candidate.mechanistic_class != "o2_gated_dissolved_co2":
        return pd.DataFrame()
    batches_with_co2 = set(co2_down["batch"].astype(str).unique()) if not co2_down.empty else {batch.batch for batch in batches}
    rows = []
    for batch in batches:
        if str(batch.batch) not in batches_with_co2:
            continue
        horizon = max(float(np.nanmax(batch.time)), float(co2_down.loc[co2_down["batch"].astype(str).eq(str(batch.batch)), "time_h_effective"].max()) if not co2_down.empty and co2_down["batch"].astype(str).eq(str(batch.batch)).any() else 0.0)
        times = np.asarray(sorted(set(np.round(np.arange(0.0, horizon + 0.25, 0.25), 8)).union(set(np.round(batch.time, 8)))), dtype=float)
        core = base.simulate(batch, theta, times)
        if core is None:
            continue
        o2 = integrate_macro_o2(theta, batch, core, params)
        co2 = integrate_dissolved_co2(theta, batch, core, params, candidate)
        qbase = co2_production_g_l_h(theta, batch, core, times)
        qeff = co2_production_g_l_h(theta, batch, core, times, candidate, params)
        o2_vals = np.interp(times, o2.index.to_numpy(dtype=float), o2["O2_mg_l"].to_numpy(dtype=float))
        uptake = np.interp(times, o2.index.to_numpy(dtype=float), o2["O2_uptake_mg_l_h"].to_numpy(dtype=float))
        phi = anaerobic_gate_from_o2(o2_vals, params)
        floor = float(np.clip(params.get("O2_crabtree_floor", O2_MACRO_DEFAULTS["O2_crabtree_floor"]), 0.0, 1.0))
        ferment_fraction = floor + (1.0 - floor) * phi
        for pos, t in enumerate(times):
            rows.append(
                {
                    "batch": str(batch.batch),
                    "time_h": float(t),
                    "O2_mg_l": float(o2_vals[pos]),
                    "O2_uptake_mg_l_h": float(uptake[pos]),
                    "anaerobic_gate": float(phi[pos]),
                    "fermentative_fraction": float(ferment_fraction[pos]),
                    "qCO2_base_g_l_h": float(qbase[pos]),
                    "qCO2_effective_g_l_h": float(qeff[pos]),
                    "CO2_dissolved_g_l": float(co2.loc[float(t), "CO2_dissolved_g_l"]) if float(t) in co2.index else np.nan,
                    "CO2_saturation_g_l": float(co2.loc[float(t), "CO2_saturation_g_l"]) if float(t) in co2.index else np.nan,
                    "qgas_g_l_h": float(co2.loc[float(t), "qgas_g_l_h"]) if float(t) in co2.index else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_o2_macro_diagnostics(diag: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if diag.empty:
        return
    for batch, group in diag.groupby("batch", sort=True):
        group = group.sort_values("time_h")
        fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
        axes[0].plot(group["time_h"], group["O2_mg_l"], color="tab:green", lw=1.5)
        axes[0].set_ylabel("O2 (mg/L)")
        axes[1].plot(group["time_h"], group["fermentative_fraction"], color="tab:orange", lw=1.5)
        axes[1].set_ylabel("ferm. fraction")
        axes[2].plot(group["time_h"], group["qCO2_base_g_l_h"], color="0.45", lw=1.2, label="base")
        axes[2].plot(group["time_h"], group["qCO2_effective_g_l_h"], color="tab:blue", lw=1.4, label="O2-gated")
        axes[2].set_ylabel("qCO2 (g/L/h)")
        axes[2].legend(fontsize=8)
        axes[3].plot(group["time_h"], group["CO2_dissolved_g_l"], color="tab:purple", lw=1.3, label="dissolved")
        axes[3].plot(group["time_h"], group["CO2_saturation_g_l"], color="tab:red", lw=1.1, ls="--", label="saturation")
        axes[3].plot(group["time_h"], group["qgas_g_l_h"], color="tab:brown", lw=1.1, label="gas release")
        axes[3].set_ylabel("CO2 (g/L)")
        axes[3].set_xlabel("effective time (h)")
        axes[3].legend(fontsize=8)
        for ax in axes:
            ax.grid(True, alpha=0.25)
        fig.suptitle(f"Macro-O2 and CO2 internal states: {batch}", y=0.995)
        fig.tight_layout()
        fig.savefig(output_dir / f"o2_macro_{batch}.png", dpi=170)
        plt.close(fig)


def plot_eigen_spectrum(spectrum: pd.DataFrame, output_dir: Path, label: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if spectrum.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.semilogy(spectrum["direction"], np.maximum(spectrum["relative_eigenvalue"], 1e-18), marker="o")
    ax.set_xlabel("FIM eigen-direction")
    ax.set_ylabel("relative eigenvalue")
    ax.set_title(label)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / f"{label}.png", dpi=170)
    plt.close(fig)


def plot_design_inputs(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if selected.empty:
        return
    for _, row in selected.sort_values("campaign_order").iterrows():
        name = str(row["candidate"])
        if name not in designs:
            continue
        design = designs[name]
        time = np.linspace(0.0, float(design.horizon_h), 500)
        fig, axes = plt.subplots(3, 1, figsize=(10, 6.5), sharex=True)
        axes[0].step(time, [base.temperature_at(design, t) for t in time], where="post", color="tab:red")
        axes[0].set_ylabel("Temperature (C)")
        axes[0].grid(True, alpha=0.25)
        n_pulses = design.pulses.get("N", tuple())
        if n_pulses:
            axes[1].vlines([t for t, _a in n_pulses], 0.0, [a * 1000.0 for _t, a in n_pulses], color="tab:green", lw=3)
        axes[1].set_ylabel("N pulse (mg/L)")
        axes[1].grid(True, alpha=0.25)
        samples = future_sample_times(design)
        axes[2].vlines(samples, 0.0, 1.0, color="tab:blue", lw=1)
        axes[2].set_ylabel("Samples")
        axes[2].set_xlabel("time (h)")
        axes[2].grid(True, alpha=0.25)
        fig.suptitle(f"{int(row['campaign_order'])}. {design.name}", fontsize=10, y=0.99)
        fig.tight_layout()
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"design_{int(row['campaign_order']):02d}.png", dpi=170)
        plt.close(fig)


def select_co2_mode(summary: pd.DataFrame) -> str:
    table = summary.copy()
    table["selection_penalty"] = 0.0
    table.loc[table["mechanistic_class"].eq("rate_proxy"), "selection_penalty"] += 25.0
    table.loc[table["mechanistic_class"].eq("empirical_effective"), "selection_penalty"] += 8.0
    table["selection_penalty"] += 35.0 * table["active_bound_count"].fillna(0).astype(float)
    table["selection_score"] = table["bic"] + table["selection_penalty"]
    best_rmse = float(table["mean_relative_rmse"].min()) if "mean_relative_rmse" in table else np.inf
    physical = table[
        table["mechanistic_class"].isin(["dissolved_co2", "o2_gated_dissolved_co2"])
        & table["active_bound_count"].fillna(0).astype(float).le(0)
        & table["mean_relative_rmse"].astype(float).le(best_rmse + 0.18)
    ].copy()
    if not physical.empty:
        return str(physical.sort_values(["selection_score", "bic", "mean_relative_rmse"]).iloc[0]["mode"])
    return str(table.sort_values(["selection_score", "bic", "data_wsse"]).iloc[0]["mode"])


def create_notebook(selected_mode: str, aroma_name: str, secondary_name: str) -> None:
    summary = pd.read_csv(RESULTS_DIR / "co2_model_selection_summary.csv")
    selected = summary[summary["mode"].eq(selected_mode)].iloc[0].to_dict()
    fim_metrics_current = json.loads((RESULTS_DIR / "fim_metrics_current.json").read_text(encoding="utf-8"))
    campaign_path = RESULTS_DIR / "selected_campaign_hybrid.csv"
    try:
        campaign = pd.read_csv(campaign_path) if campaign_path.exists() and campaign_path.stat().st_size > 2 else pd.DataFrame()
    except pd.errors.EmptyDataError:
        campaign = pd.DataFrame()
    top_campaign = ", ".join(campaign["candidate"].astype(str).head(3)) if not campaign.empty else "not computed"

    nb = nbformat.v4.new_notebook()
    cells = [
        nbformat.v4.new_markdown_cell(
            rf"""# Pilot 2025 CO2-solubility/O2 integrated calibration, estimability and DOE

This notebook documents the deterministic workflow used for the pilot-scale natural-must dataset. It starts from the reduced extended fermentation model already used in the previous notebooks, then replaces the purely empirical CO2 gas-flow lag with a dissolved-CO2 buffer model and a parsimonious macroscopic oxygen transition model.

Selected inherited secondary structure: `{secondary_name}`.

Selected inherited aroma structure: `{aroma_name}`.

Selected CO2 structure after the new benchmark: `{selected_mode}`.
"""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import json
import pandas as pd
from IPython.display import Image, display
RESULTS = Path('results/co2_solubility_integrated_doe')
PLOTS = RESULTS / 'plots'
print(RESULTS.resolve())
"""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## 1. Full model used as baseline

The core fermentation model tracks viable biomass \(X\), dead biomass \(X_d\), assimilable nitrogen \(N\), glucose \(G\), fructose \(F\), ethanol \(E\), and glycerol \(Gly\):

$$\frac{dX}{dt}=(\mu-k_d)X+u_X$$

$$\frac{dX_d}{dt}=k_dX$$

$$\frac{dN}{dt}=-q_N f_N(T,N)X+u_N$$

$$\frac{dG}{dt}=-\left(q_{XG}f_N+q_{EG}f_G+m\frac{G}{G+F}\right)X+u_G$$

$$\frac{dF}{dt}=-\left(q_{XF}f_N+q_{EF}f_F+m\frac{F}{G+F}\right)X+u_F$$

$$\frac{dE}{dt}=(\beta_G f_G+\beta_F f_F)X+u_E$$

$$\frac{dGly}{dt}=(\gamma_G f_G+\gamma_F f_F)X$$

The secondary layer keeps pyruvate, acetaldehyde, acetate, and oxygen as the reduced chemical-proxy states. In this pilot notebook, the CO2 block additionally uses a macro-O2 state to delay the effective anaerobic ethanol/CO2 source when the must is initially air-saturated. The aroma layer predicts retained liquid concentration and accumulated condenser-equivalent loss for ethyl acetate, isoamyl acetate, and ethyl octanoate.
"""
        ),
        nbformat.v4.new_markdown_cell("### Data read and curation check"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_curation_decisions.csv')"),
        nbformat.v4.new_code_cell(
            """for path in sorted((PLOTS / 'data').glob('data_read_check_*.png')):
    display(Image(filename=str(path)))
for path in sorted((PLOTS / 'data').glob('co2_curated_*.png')):
    display(Image(filename=str(path)))
"""
        ),
        nbformat.v4.new_markdown_cell(
            r"""**Plain-language interpretation.** The workbook is read batch-by-batch, `25150` and `25151` CO2 files are excluded, and `25171` is restarted at the annotated `Pre reinoculo` sample. This makes the post-reinoculation segment the effective fermentation start instead of treating the non-viable inoculum period as model lag."""
        ),
        nbformat.v4.new_markdown_cell("## 2. Initial simulation before CO2 reformulation"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'initial_fit_metrics.csv').sort_values(['group','relative_rmse']).head(40)"),
        nbformat.v4.new_code_cell(
            """for path in sorted((PLOTS / 'initial_fit').glob('initial_fit_*.png')):
    display(Image(filename=str(path)))
"""
        ),
        nbformat.v4.new_markdown_cell(
            r"""**Plain-language interpretation.** This block is the reference fit before adding the dissolved-CO2 and macro-O2 states. If the model predicts CO2 release far earlier than the sensor while ethanol/sugar curves remain acceptable, the issue is not only gas-liquid solubility: the model is also assuming anaerobic ethanol/CO2 production from the start."""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## 3. CO2 solubility, macro-O2 transition model and benchmark

The dissolved CO2 state is:

$$\frac{dC_{CO2,L}}{dt}=r_{CO2,prod}-r_{CO2,gas}$$

The original gas-source proxy is tied to ethanol production:

$$r_{CO2,prod}=\frac{44.01}{2\cdot46.07}r_E$$

The O2-gated candidates replace that source with:

$$r_{CO2,prod}^{eff}=r_{CO2,ferm}^{base}\left[f_C+(1-f_C)\phi_{ana}(O_2)\right]+r_{CO2,resp}$$

where the Crabtree floor \(f_C\) prevents aerobic conditions from fully shutting down fermentation at high sugar, and

$$\phi_{ana}(O_2)=\frac{K_{ana}^{n}}{K_{ana}^{n}+O_2^{n}}.$$

The macro-O2 state is initialized near air saturation for fresh must:

$$O_2(0)=O_2^*(T,E,G,F),$$

with a low effective fraction for batch 25171 because its new \(t=0\) is post-reinoculation after the original non-viable inoculum period. The state evolves as:

$$\frac{dO_2}{dt}=k_{La,O2}(O_2^*-O_2)-q_{O2}X\frac{O_2}{K_{O2}+O_2}.$$

The respiration contribution is kept stoichiometric and small:

$$r_{CO2,resp}=\frac{44.01}{32.00}\frac{q_{O2}X\,O_2}{1000(K_{O2}+O_2)}.$$

Gas release from the liquid still follows:

$$r_{CO2,gas}=k_{rel}\max(C_{CO2,L}-C^*_{CO2},0).$$

The saturation concentration is represented as a process correlation:

$$C^*_{CO2}=s_{CO2}\,1.69\exp[-0.032(T-20)]\exp(0.0016E)\exp[-0.0012(G+F)].$$

This is not a purely arbitrary lag. It is a two-stage physical/effective model: fresh must can contain oxygen that suppresses the anaerobic ethanol/CO2 source, then generated CO2 fills the dissolved pool before gas flow appears once the pool approaches supersaturation. The sensor reports \(L/min\), so each batch is allowed a linear scale factor from model specific release \((g/L/h)\) to measured gas flow.
"""
        ),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_model_selection_summary.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_metrics_selected.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'o2_macro_diagnostics_selected.csv').head(30)"),
        nbformat.v4.new_code_cell(
            """for path in sorted((PLOTS / 'co2_benchmark').glob('co2_benchmark_*.png')):
    display(Image(filename=str(path)))
for path in sorted((PLOTS / 'o2_macro').glob('o2_macro_*.png')):
    display(Image(filename=str(path)))
"""
        ),
        nbformat.v4.new_markdown_cell(
            rf"""**Plain-language interpretation.** The selected model is `{selected_mode}`. Its benchmark row has BIC `{selected.get('bic', float('nan')):.2f}`, data WSSE `{selected.get('data_wsse', float('nan')):.2f}`, and selection score `{selected.get('selection_score', float('nan')):.2f}`. The O2-gated candidates are accepted only if they improve early CO2 timing without requiring active parameter bounds or a purely empirical lag."""
        ),
        nbformat.v4.new_markdown_cell("## 4. Integrated post-CO2 calibration and validation"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'fit_integrated_selected.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'theta_selected_integrated.csv', index_col=0).head(100)"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'final_fit_metrics.csv').sort_values(['group','relative_rmse']).head(60)"),
        nbformat.v4.new_code_cell(
            """for path in sorted((PLOTS / 'final_fit').glob('final_fit_*.png')):
    display(Image(filename=str(path)))
"""
        ),
        nbformat.v4.new_markdown_cell(
            r"""**Plain-language interpretation.** The integrated fit adjusts the selected aroma parameters while keeping the step-0 CO2 solubility parameters fixed. This prevents aroma residuals from moving the physical CO2 buffer to a non-physical boundary. Good retained-aroma fit but poor condensate fit indicates a partition/stripping issue; poor retained and total fit indicates a synthesis-kinetics issue."""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## 5. CO2/O2 FIM, eigenvalues and estimability

The full coupled FIM, including all aroma states, is computationally expensive because each finite-difference perturbation must reintegrate liquid/condensate aroma partition. For this O2-structure iteration, the notebook reports a focused CO2/O2 Fisher Information Matrix for the selected gas-transfer parameters. This is the correct diagnostic for deciding whether the new O2-gated gas model improves the online CO2 direction.

The current-data Fisher Information Matrix is computed by finite differences in log-parameter coordinates:

$$J_{:,j}\approx \frac{r(\theta_j e^{\Delta})-r(\theta_j e^{-\Delta})}{2\Delta}$$

$$F=J^TJ.$$

Eigenvalues close to zero indicate practically weak directions; the weakest eigenvectors show which parameter combinations are confounded. The full aroma-coupled FIM should be rerun later with an aroma-partition cache if the goal is to redesign the whole aroma campaign.
"""
        ),
        nbformat.v4.new_code_cell("json.load(open(RESULTS / 'co2_o2_target_parameters.json'))"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_eigen_spectrum_current.csv').head(20)"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_weak_directions_current.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_parameter_estimability_current.csv')"),
        nbformat.v4.new_code_cell("display(Image(filename=str(PLOTS / 'fim' / 'co2_o2_current_relative_eigenvalues.png')))"),
        nbformat.v4.new_markdown_cell(
            r"""**Plain-language interpretation.** This focused FIM answers a narrower question than the previous global FIM: can the online CO2 data distinguish the selected gas-release parameters after adding macro-O2 gating? If this block is well-conditioned but aroma fits remain weak, the remaining limitation is not the CO2/O2 gas timing alone."""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## 6. CO2/O2 model-based DOE

Candidate natural-must experiments add an experiment FIM to the current-data FIM:

$$F_{total}=F_{current}+\sum_i F_i.$$

The ranking reports D-optimality through \(\log\det(F)\), E-optimality through the minimum eigenvalue, and a hybrid score:

$$\Phi_{hybrid}=\log\det(F)+2\log(\lambda_{min}/\lambda_{max})-0.05\log(\mathrm{trace}(F^{-1})).$$

This focused DOE uses the same additive multi-experiment logic used in the Dowling/Pyomo DoE examples: the current online CO2 data act as prior information and each candidate design contributes incremental information about the selected CO2/O2 gas-transfer parameters. It is not a replacement for a full aroma-coupled MBDoE.
"""
        ),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_candidate_ranking.csv').head(15)"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_selected_campaign_hybrid.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_o2_eigen_spectrum_current_plus_campaign.csv').head(20)"),
        nbformat.v4.new_code_cell(
            """for path in sorted((PLOTS / 'co2_o2_designs').glob('design_*.png')):
    display(Image(filename=str(path)))
"""
        ),
        nbformat.v4.new_markdown_cell(
            rf"""**Plain-language interpretation.** The selected campaign begins with: {top_campaign}. Designs are chosen because they improve the weakest information directions after accounting for the current pilot data, not because their curves look intuitively different."""
        ),
    ]
    nb["cells"] = cells
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pilot 2025 dissolved-CO2/macro-O2 integrated calibration and DOE.")
    parser.add_argument("--n-starts-co2", type=int, default=6)
    parser.add_argument("--n-starts-integrated", type=int, default=5)
    parser.add_argument("--max-nfev-co2", type=int, default=140)
    parser.add_argument("--max-nfev-integrated", type=int, default=100)
    parser.add_argument("--sensitivity-step", type=float, default=0.015)
    parser.add_argument("--campaign-size", type=int, default=6)
    parser.add_argument("--skip-fim", action="store_true")
    parser.add_argument("--skip-doe", action="store_true")
    parser.add_argument("--modes", default="instant,old_lag_threshold,solubility_fixed,solubility_scaled,solubility_o2_literature,solubility_o2_slow_transition,solubility_o2_qfit")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] pilot workbook, CO2 sensors, inherited model context", flush=True)
    model_data, curated_co2, co2_decisions = pilot.load_pilot_model_data()
    co2_down = pilot.downsample_co2(curated_co2, dt_h=CO2_DOWNSAMPLE_H)
    primary_batches = pilot.make_primary_batches(model_data)
    extended_batches = pilot.make_extended_batches(model_data)
    theta0, aroma_name, aroma_variant, secondary_name, secondary_variant = load_selected_context()

    model_data.to_csv(RESULTS_DIR / "pilot_model_data_used.csv", index=False)
    co2_decisions.to_csv(RESULTS_DIR / "co2_curation_decisions.csv", index=False)
    co2_down.to_csv(RESULTS_DIR / "co2_downsampled.csv", index=False)
    pd.Series(theta0).to_csv(RESULTS_DIR / "theta_start.csv")
    (RESULTS_DIR / "selected_aroma_model_inherited.txt").write_text(aroma_name, encoding="utf-8")
    (RESULTS_DIR / "selected_secondary_model_inherited.txt").write_text(secondary_name, encoding="utf-8")

    print("[plot] data overview", flush=True)
    plot_data_overview(model_data, co2_down, PLOT_DIR / "data")

    candidates = co2_candidate_library()
    modes = [name.strip() for name in str(args.modes).split(",") if name.strip()]
    modes = [name for name in modes if name in candidates]
    if not modes:
        modes = ["solubility_scaled"]

    print("[initial] baseline prediction with instant CO2", flush=True)
    initial_candidate = candidates["instant"]
    initial_pred, initial_metrics = all_state_predictions(
        theta0,
        primary_batches,
        extended_batches,
        co2_down,
        aroma_variant,
        initial_candidate,
        {},
        "initial",
    )
    initial_pred.to_csv(RESULTS_DIR / "initial_predictions.csv", index=False)
    initial_metrics.to_csv(RESULTS_DIR / "initial_fit_metrics.csv", index=False)
    plot_prediction_by_batch(initial_pred, PLOT_DIR / "initial_fit", "initial")

    print("[co2] benchmark candidate structures", flush=True)
    fit_rows = []
    metric_rows = []
    pred_rows = []
    params_by_mode: dict[str, dict[str, float]] = {}
    by_name = {batch.batch: batch for batch in extended_batches}
    for idx, mode_name in enumerate(modes, start=1):
        candidate = candidates[mode_name]
        print(f"[co2] {idx}/{len(modes)} {mode_name}", flush=True)
        params, rows = fit_co2_candidate(
            theta0,
            by_name,
            co2_down,
            candidate,
            n_starts=args.n_starts_co2,
            max_nfev=args.max_nfev_co2,
            seed=314 + idx,
        )
        params_by_mode[mode_name] = params
        rows.to_csv(RESULTS_DIR / f"fit_co2_{mode_name}.csv", index=False)
        fit_rows.append(rows)
        metrics, preds = co2_metrics(theta0, by_name, co2_down, candidate, params)
        if not metrics.empty:
            metrics["params_json"] = json.dumps(params)
            metric_rows.append(metrics)
        if not preds.empty:
            pred_rows.append(preds)

    fit_all = pd.concat(fit_rows, ignore_index=True, sort=False)
    metric_all = pd.concat(metric_rows, ignore_index=True, sort=False) if metric_rows else pd.DataFrame()
    pred_all = pd.concat(pred_rows, ignore_index=True, sort=False) if pred_rows else pd.DataFrame()
    fit_all.to_csv(RESULTS_DIR / "fit_co2_all.csv", index=False)
    metric_all.to_csv(RESULTS_DIR / "co2_metrics_all.csv", index=False)
    pred_all.to_csv(RESULTS_DIR / "co2_predictions_all.csv", index=False)
    plot_co2_benchmark(pred_all, PLOT_DIR / "co2_benchmark")

    summary_rows = []
    for mode_name, group in fit_all.groupby("mode", sort=True):
        best = group.sort_values("fit_selection_score" if "fit_selection_score" in group else "data_wsse").iloc[0]
        candidate = candidates[mode_name]
        data_wsse = float(best["data_wsse"])
        n = int(best["n_data_residuals"])
        k = int(best["n_parameters"])
        aic, bic = information_criteria(data_wsse, n, max(k, 1))
        metrics_mode = metric_all[metric_all["mode"].eq(mode_name)]
        summary_rows.append(
            {
                "mode": mode_name,
                "mechanistic_class": candidate.mechanistic_class,
                "description": candidate.description,
                "data_wsse": data_wsse,
                "n_data_residuals": n,
                "n_parameters": k,
                "aicc": aic,
                "bic": bic,
                "active_bound_count": int(best.get("active_bound_count", 0)),
                "mean_relative_rmse": float(metrics_mode["relative_rmse"].mean()) if not metrics_mode.empty else np.nan,
                "min_corr": float(metrics_mode["corr"].min()) if not metrics_mode.empty else np.nan,
                "params_json": json.dumps(params_by_mode.get(mode_name, {})),
            }
        )
    co2_summary = pd.DataFrame(summary_rows)
    selected_mode = select_co2_mode(co2_summary)
    co2_summary["selected"] = co2_summary["mode"].eq(selected_mode)
    if "selection_score" not in co2_summary.columns:
        tmp = co2_summary.copy()
        tmp["selection_penalty"] = 0.0
        tmp.loc[tmp["mechanistic_class"].eq("rate_proxy"), "selection_penalty"] += 25.0
        tmp.loc[tmp["mechanistic_class"].eq("empirical_effective"), "selection_penalty"] += 8.0
        tmp["selection_penalty"] += 35.0 * tmp["active_bound_count"].fillna(0).astype(float)
        co2_summary["selection_penalty"] = tmp["selection_penalty"]
        co2_summary["selection_score"] = tmp["bic"] + tmp["selection_penalty"]
    co2_summary = co2_summary.sort_values(["selection_score", "bic"]).reset_index(drop=True)
    co2_summary.to_csv(RESULTS_DIR / "co2_model_selection_summary.csv", index=False)
    (RESULTS_DIR / "selected_co2_model.txt").write_text(selected_mode, encoding="utf-8")

    selected_candidate = candidates[selected_mode]
    selected_params = params_by_mode[selected_mode]
    pd.Series(selected_params).to_csv(RESULTS_DIR / "co2_params_selected.csv")
    selected_metrics = metric_all[metric_all["mode"].eq(selected_mode)].copy()
    selected_metrics.to_csv(RESULTS_DIR / "co2_metrics_selected.csv", index=False)
    o2_diag = o2_macro_diagnostics(theta0, extended_batches, co2_down, selected_candidate, selected_params)
    o2_diag.to_csv(RESULTS_DIR / "o2_macro_diagnostics_selected.csv", index=False)
    plot_o2_macro_diagnostics(o2_diag, PLOT_DIR / "o2_macro")

    print(f"[integrated] selected CO2 model: {selected_mode}", flush=True)
    theta_fit, fit_integrated = fit_integrated_model(
        theta0,
        primary_batches,
        extended_batches,
        co2_down,
        aroma_variant,
        selected_candidate,
        selected_params,
        n_starts=args.n_starts_integrated,
        max_nfev=args.max_nfev_integrated,
        seed=771,
    )
    fit_integrated.to_csv(RESULTS_DIR / "fit_integrated_selected.csv", index=False)
    pd.Series(theta_fit).to_csv(RESULTS_DIR / "theta_selected_integrated.csv")

    print("[validate] final prediction plots", flush=True)
    final_pred, final_metrics = all_state_predictions(
        theta_fit,
        primary_batches,
        extended_batches,
        co2_down,
        aroma_variant,
        selected_candidate,
        candidate_params_from_theta(theta_fit, selected_candidate),
        "final",
    )
    final_pred.to_csv(RESULTS_DIR / "final_predictions.csv", index=False)
    final_metrics.to_csv(RESULTS_DIR / "final_fit_metrics.csv", index=False)
    plot_prediction_by_batch(final_pred, PLOT_DIR / "final_fit", "final")

    target_parameters = tuple(
        dict.fromkeys(
            pilot.CORE_TARGETS
            + tuple(secondary_variant.parameters)
            + tuple(aroma_variant.parameters)
            + tuple(selected_candidate.parameters)
        )
    )
    (RESULTS_DIR / "target_parameters.json").write_text(json.dumps(list(target_parameters), indent=2), encoding="utf-8")

    fim = np.eye(len(target_parameters), dtype=float) * 1e-9
    if not args.skip_fim:
        print("[fim] current-data FIM", flush=True)
        residual_fun = lambda th: current_residual(th, primary_batches, extended_batches, co2_down, aroma_variant, selected_candidate)
        jac, base_res = finite_difference_jacobian(theta_fit, target_parameters, residual_fun, args.sensitivity_step)
        fim = jac.T @ jac
        fim = 0.5 * (fim + fim.T)
        pd.DataFrame(jac, columns=target_parameters).to_csv(RESULTS_DIR / "jacobian_current.csv", index=False)
        pd.DataFrame(fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_current.csv")
        fim_current_metrics = fim_metrics(fim)
        (RESULTS_DIR / "fim_metrics_current.json").write_text(json.dumps(fim_current_metrics, indent=2), encoding="utf-8")
        if rescale_FIM is not None:
            try:
                scaled = rescale_FIM(fim, np.asarray([theta_fit[p] for p in target_parameters], dtype=float))
                pd.DataFrame(scaled, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_current_pyomodoe_rescaled.csv")
            except Exception as exc:
                (RESULTS_DIR / "pyomodoe_rescale_warning.txt").write_text(str(exc), encoding="utf-8")
        spectrum, weak, estimability = fim_diagnostics(fim, theta_fit, target_parameters, "current")
        spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_current.csv", index=False)
        weak.to_csv(RESULTS_DIR / "weak_directions_current.csv", index=False)
        estimability.to_csv(RESULTS_DIR / "parameter_estimability_current.csv", index=False)
        plot_eigen_spectrum(spectrum, PLOT_DIR / "fim", "current_relative_eigenvalues")
    else:
        (RESULTS_DIR / "fim_metrics_current.json").write_text(json.dumps(fim_metrics(fim), indent=2), encoding="utf-8")
        spectrum, weak, estimability = fim_diagnostics(fim, theta_fit, target_parameters, "current")
        spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_current.csv", index=False)
        weak.to_csv(RESULTS_DIR / "weak_directions_current.csv", index=False)
        estimability.to_csv(RESULTS_DIR / "parameter_estimability_current.csv", index=False)
        plot_eigen_spectrum(spectrum, PLOT_DIR / "fim", "current_relative_eigenvalues")

    selected_campaign = pd.DataFrame()
    if not args.skip_doe and not args.skip_fim:
        print("[doe] candidate library and FIMs", flush=True)
        designs = add_solubility_designs(model_data, aroma_sel.add_extra_designs(model_data, pilot.natural_candidate_designs(model_data)))
        design_rows = []
        for design in designs.values():
            design_rows.append(
                {
                    "candidate": design.name,
                    "family": design.family,
                    "medium": design.medium,
                    "horizon_h": design.horizon_h,
                    "initials_json": json.dumps(design.initials),
                    "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                    "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                    "rationale": design.rationale,
                }
            )
        pd.DataFrame(design_rows).to_csv(RESULTS_DIR / "candidate_design_library.csv", index=False)
        candidate_fims = {}
        for idx, (name, design) in enumerate(designs.items(), start=1):
            print(f"[doe] {idx}/{len(designs)} {name}", flush=True)
            cfim = candidate_fim(theta_fit, design, aroma_variant, selected_candidate, target_parameters, args.sensitivity_step)
            candidate_fims[name] = cfim
            pd.DataFrame(cfim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / f"candidate_fim_{name}.csv")
        ranking = rank_candidates(candidate_fims, designs, fim, target_parameters)
        ranking.to_csv(RESULTS_DIR / "candidate_ranking.csv", index=False)
        selected_campaign, campaign_fim = greedy_select(candidate_fims, designs, fim, target_parameters, args.campaign_size, objective="hybrid")
        selected_campaign.to_csv(RESULTS_DIR / "selected_campaign_hybrid.csv", index=False)
        pd.DataFrame(campaign_fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_current_plus_campaign.csv")
        spec_after, weak_after, est_after = fim_diagnostics(campaign_fim, theta_fit, target_parameters, "current_plus_campaign")
        spec_after.to_csv(RESULTS_DIR / "eigen_spectrum_current_plus_campaign.csv", index=False)
        weak_after.to_csv(RESULTS_DIR / "weak_directions_current_plus_campaign.csv", index=False)
        est_after.to_csv(RESULTS_DIR / "parameter_estimability_current_plus_campaign.csv", index=False)
        plot_eigen_spectrum(spec_after, PLOT_DIR / "fim", "current_plus_campaign_relative_eigenvalues")
        plot_design_inputs(selected_campaign, designs, PLOT_DIR / "designs")
    else:
        pd.DataFrame().to_csv(RESULTS_DIR / "candidate_ranking.csv", index=False)
        pd.DataFrame().to_csv(RESULTS_DIR / "selected_campaign_hybrid.csv", index=False)

    create_notebook(selected_mode, aroma_name, secondary_name)
    print(f"[done] selected CO2 model: {selected_mode}", flush=True)
    print(f"[done] notebook: {NOTEBOOK_PATH}", flush=True)
    print(f"[done] results: {RESULTS_DIR}", flush=True)


if __name__ == "__main__":
    main()
