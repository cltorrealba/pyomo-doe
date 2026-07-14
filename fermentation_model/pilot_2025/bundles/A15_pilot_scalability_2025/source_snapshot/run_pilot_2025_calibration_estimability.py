from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
PILOT_DIR = SCRIPT_DIR.parent if SCRIPT_DIR.name == "support" else SCRIPT_DIR
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(PILOT_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_DIR))
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

import pilot_2025_data_loader as pilot_loader
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = PILOT_DIR / "results" / "calibration_estimability"
NOTEBOOK_PATH = PILOT_DIR / "pilot_2025_calibration_estimability.ipynb"
EXECUTED_NOTEBOOK_PATH = PILOT_DIR / "pilot_2025_calibration_estimability.executed.ipynb"

CORE_TARGETS = base.REDUCED11
SECONDARY_TARGETS = v2.V2_REDUCED_PARAMETERS
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
TARGET_PARAMETERS = CORE_TARGETS + SECONDARY_TARGETS + AROMA_TARGETS

AROMA_BOUNDS = {name: joint.AROMA_BOUNDS[name] for name in AROMA_TARGETS}
PARAMETER_BOUNDS = {
    **{name: base.PARAMETER_BOUNDS[name] for name in base.FULL17},
    **v2.V2_BOUNDS,
    **joint.AROMA_BOUNDS,
}

CO2_RATE_SIGMA = 0.35
CO2_DOWNSAMPLE_H = 2.0
NATURAL_DESIGN_HORIZON_H = 300.0

AROMA_CONDENSATE_COLUMNS = {
    "ethyl_acetate": "ethyl_acetate_condensate",
    "isoamyl_acetate": "isoamyl_acetate_condensate",
    "ethyl_octanoate": "ethyl_octanoate_condensate",
}
AROMA_LIQUID_FLOOR = {
    "ethyl_acetate": 0.40,
    "isoamyl_acetate": 0.20,
    "ethyl_octanoate": 0.006,
}
AROMA_CONDENSATE_FLOOR = {
    "ethyl_acetate": 0.010,
    "isoamyl_acetate": 0.040,
    "ethyl_octanoate": 0.002,
}


def parameter_bounds(name: str) -> tuple[float, float]:
    return PARAMETER_BOUNDS[name]


def clip_extended(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in PARAMETER_BOUNDS.items():
        if name in out and np.isfinite(float(out[name])):
            out[name] = float(np.clip(float(out[name]), lb, ub))
    return out


def load_initial_theta() -> dict[str, float]:
    theta = joint.load_reference_theta()
    theta.update(v2.default_theta_v2(theta))
    secondary_path = FERMENTATION_MODEL_DIR / "results" / "secondary_v2_model_evaluation" / "theta_secondary_v2_reduced_o2fixed.csv"
    if secondary_path.exists():
        loaded = pd.read_csv(secondary_path, index_col=0).iloc[:, 0].to_dict()
        theta.update({str(k): float(v) for k, v in loaded.items() if str(k) in v2.V2_BOUNDS and np.isfinite(float(v))})
    aroma_path = FERMENTATION_MODEL_DIR / "results" / "secondary_joint_campaign_doe" / "theta_secondary_joint.csv"
    if aroma_path.exists():
        loaded = pd.read_csv(aroma_path, index_col=0).iloc[:, 0].to_dict()
        theta.update({str(k): float(v) for k, v in loaded.items() if str(k) in joint.AROMA_BOUNDS and np.isfinite(float(v))})
    return clip_extended(theta)


def load_pilot_model_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    data = pilot_loader.load_pilot_calibration_data()
    raw_co2 = pilot_loader.load_all_co2_sensor_data(calibration_data=data)
    curated_co2, co2_decisions = pilot_loader.curate_co2_sensor_data(raw_co2, calibration_data=data)
    model = data.copy()
    model["time_h_original"] = model["time_h"]
    for _, decision in co2_decisions.iterrows():
        if not bool(decision.get("used_for_calibration", False)):
            continue
        activation = float(decision.get("activation_time_h", 0.0))
        if not np.isfinite(activation) or activation <= 1e-9:
            continue
        batch = str(decision["batch"])
        batch_rows = model[model["batch"].astype(str).eq(batch)].copy()
        if batch_rows.empty:
            continue
        initial = _interpolate_initial_row(batch_rows, activation)
        shifted = batch_rows[batch_rows["time_h_original"].ge(activation)].copy()
        shifted["time_h"] = shifted["time_h_original"] - activation
        shifted = pd.concat([initial.to_frame().T, shifted], ignore_index=True, sort=False)
        shifted = shifted.sort_values("time_h").drop_duplicates("time_h", keep="first")
        model = pd.concat([model[~model["batch"].astype(str).eq(batch)], shifted], ignore_index=True, sort=False)
    model["medium"] = "natural"
    model["X_dead_kg_m3"] = np.nan
    model["pyruvic_acid"] = pd.to_numeric(model["pyruvic_acid_mg_l"], errors="coerce")
    model["acetaldehyde"] = pd.to_numeric(model["acetaldehyde_mg_l"], errors="coerce")
    model["acetic_acid"] = pd.to_numeric(model.get("acetic_acid_g_l"), errors="coerce")
    model["DO_mg_l"] = pd.to_numeric(model.get("DO_mg_l"), errors="coerce")
    model["glycerol_g_l"] = pd.to_numeric(model["glycerol_g_l"], errors="coerce")
    return model.sort_values(["batch", "time_h"]).reset_index(drop=True), curated_co2, co2_decisions


def _interpolate_initial_row(batch_rows: pd.DataFrame, activation_time_h: float) -> pd.Series:
    rows = batch_rows.sort_values("time_h_original").copy()
    initial = rows.iloc[0].copy()
    after = rows[rows["time_h_original"].ge(float(activation_time_h))]
    if not after.empty:
        initial = after.iloc[0].copy()
    initial["time_h_original"] = float(activation_time_h)
    initial["time_h"] = 0.0
    initial["sample_id"] = f"{initial.get('batch', '')}_effective_initial"
    for col in rows.columns:
        if col in {"batch", "source_file", "sample_id", "system", "medium", "scale", "E_unit_assumption", "fecha_hora"}:
            continue
        values = pd.to_numeric(rows[col], errors="coerce")
        times = pd.to_numeric(rows["time_h_original"], errors="coerce")
        mask = values.notna() & times.notna()
        if mask.sum() >= 2 and float(times[mask].min()) <= activation_time_h <= float(times[mask].max()):
            initial[col] = float(np.interp(float(activation_time_h), times[mask].to_numpy(dtype=float), values[mask].to_numpy(dtype=float)))
    initial["time_h_original"] = float(activation_time_h)
    initial["time_h"] = 0.0
    return initial


def make_primary_batches(model_data: pd.DataFrame) -> list[base.BatchData]:
    return base.make_batches(model_data)


def make_extended_batches(model_data: pd.DataFrame) -> list[base.BatchData]:
    batches = joint.make_secondary_batches(model_data)
    by_key = {(str(batch.medium), str(batch.batch)): batch for batch in batches}
    out: list[base.BatchData] = []
    for (medium, batch_name), group in model_data.groupby(["medium", "batch"], sort=True):
        batch = by_key.get((str(medium), str(batch_name)))
        if batch is None:
            continue
        aligned = group.sort_values("time_h").drop_duplicates("time_h").set_index("time_h").reindex(batch.time)
        obs = dict(batch.observations)
        initials = dict(batch.initials)
        for species, cond_col in AROMA_CONDENSATE_COLUMNS.items():
            total_col = joint.AROMA_COLUMNS[species]
            total = pd.to_numeric(aligned.get(total_col), errors="coerce")
            cond = pd.to_numeric(aligned.get(cond_col), errors="coerce")
            obs[cond_col] = cond.to_numpy(dtype=float)
            first_total = total.dropna()
            if not first_total.empty:
                first_idx = first_total.index[0]
                cond0 = float(cond.loc[first_idx]) if first_idx in cond.index and np.isfinite(cond.loc[first_idx]) else 0.0
                initials[f"{species}_liq"] = max(float(first_total.iloc[0]) - cond0, 0.0)
                initials[f"{species}_cond"] = max(cond0, 0.0)
            else:
                initials[f"{species}_liq"] = 0.0
                initials[f"{species}_cond"] = 0.0
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


def log_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    return np.asarray([math.log(float(theta[name])) for name in parameters], dtype=float)


def theta_from_log(x: np.ndarray, base_theta: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base_theta)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        theta[name] = float(math.exp(float(value)))
    return clip_extended(theta)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    low = []
    high = []
    for name in parameters:
        lb, ub = parameter_bounds(name)
        low.append(math.log(float(lb)))
        high.append(math.log(float(ub)))
    return np.asarray(low, dtype=float), np.asarray(high, dtype=float)


def fit_core_multistart(
    batches: list[base.BatchData],
    theta0: dict[str, float],
    n_starts: int,
    max_nfev: int,
    seed: int = 25,
) -> tuple[dict[str, float], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    summaries = []
    best_theta = dict(theta0)
    best_wsse = np.inf
    for idx in range(int(n_starts)):
        start = dict(theta0)
        if idx > 0:
            for name in CORE_TARGETS:
                lb, ub = base.PARAMETER_BOUNDS[name]
                start[name] = float(np.clip(theta0[name] * math.exp(rng.normal(0.0, 0.45)), lb, ub))
        label = f"pilot_core_l2_multistart_{idx:02d}"
        theta_hat, summary = base.fit_parameters(
            label,
            batches,
            start,
            CORE_TARGETS,
            max_nfev=max_nfev,
            l2_reference=theta0,
            l2_parameters=CORE_TARGETS,
            l2_lambda=0.35,
        )
        summaries.append(summary)
        if float(summary["final_wsse"]) < best_wsse:
            best_wsse = float(summary["final_wsse"])
            best_theta = theta_hat
    return clip_extended(best_theta), pd.DataFrame(summaries)


def fit_secondary_v2(
    batches: list[base.BatchData],
    theta_core: dict[str, float],
    max_nfev: int,
) -> tuple[dict[str, float], pd.DataFrame, dict[str, pd.DataFrame]]:
    core_cache = joint.precompute_core_cache(batches, theta_core)
    theta_sec, summary = v2.fit_v2(
        theta_core,
        batches,
        core_cache,
        max_nfev=max_nfev,
        fit_parameters=SECONDARY_TARGETS,
        model_label="pilot_secondary_v2_reduced_o2fixed",
    )
    theta = dict(theta_core)
    theta.update(theta_sec)
    return clip_extended(theta), summary, core_cache


def aroma_sigma(species: str, pool: str, observed: np.ndarray) -> np.ndarray:
    observed = np.asarray(observed, dtype=float)
    if pool == "condensate":
        floor = AROMA_CONDENSATE_FLOOR[species]
        rel = 0.20
    else:
        floor = AROMA_LIQUID_FLOOR[species]
        rel = 0.12
    return np.maximum(float(floor), rel * np.maximum(np.abs(observed), float(floor)))


def integrate_aroma_pilot_euler(
    batch: base.BatchData | base.FutureDesign,
    theta: dict[str, float],
    core: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    time = core.index.to_numpy(dtype=float)
    liquid0 = [float(batch.initials.get(f"{species}_liq", 0.0)) for species in joint.AROMA_SPECIES]
    loss0 = [float(batch.initials.get(f"{species}_cond", 0.0)) for species in joint.AROMA_SPECIES]
    y = np.asarray(liquid0 + loss0, dtype=float)
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(joint.aroma_rhs(float(time[idx - 1]), y, theta, batch, core), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    arr = np.asarray(rows, dtype=float)
    liquid = pd.DataFrame(arr[:, : len(joint.AROMA_SPECIES)], index=time, columns=joint.AROMA_SPECIES)
    condensate = pd.DataFrame(arr[:, len(joint.AROMA_SPECIES) :], index=time, columns=joint.AROMA_SPECIES)
    return liquid, condensate


def aroma_residual_vector(
    theta: dict[str, float],
    batches: list[base.BatchData],
    core_cache: dict[str, pd.DataFrame] | None = None,
) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in batches:
        if core_cache is not None and batch.label in core_cache:
            core = core_cache[batch.label]
        else:
            time = joint._dense_time_grid(batch, batch.time)
            core = base.simulate(batch, theta, time)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        try:
            aromas_liq, aromas_cond = integrate_aroma_pilot_euler(batch, theta, core)
        except Exception:
            return np.ones(1000, dtype=float) * 1e6
        for species, column in joint.AROMA_COLUMNS.items():
            total_obs = np.asarray(batch.observations.get(column, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            cond_col = AROMA_CONDENSATE_COLUMNS[species]
            cond_obs = np.asarray(batch.observations.get(cond_col, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            pred_liq = aromas_liq.loc[batch.time, species].to_numpy(dtype=float)
            pred_cond = aromas_cond.loc[batch.time, species].to_numpy(dtype=float)

            mask_both = np.isfinite(total_obs) & np.isfinite(cond_obs)
            if mask_both.any():
                retained = np.maximum(total_obs[mask_both] - cond_obs[mask_both], 0.0)
                residuals.append((pred_liq[mask_both] - retained) / aroma_sigma(species, "liquid", retained))
                residuals.append((pred_cond[mask_both] - cond_obs[mask_both]) / aroma_sigma(species, "condensate", cond_obs[mask_both]))

            mask_total_only = np.isfinite(total_obs) & ~np.isfinite(cond_obs)
            if mask_total_only.any():
                pred_total = pred_liq[mask_total_only] + pred_cond[mask_total_only]
                residuals.append((pred_total - total_obs[mask_total_only]) / aroma_sigma(species, "liquid", total_obs[mask_total_only]))

            mask_cond_only = np.isfinite(cond_obs) & ~np.isfinite(total_obs)
            if mask_cond_only.any():
                residuals.append((pred_cond[mask_cond_only] - cond_obs[mask_cond_only]) / aroma_sigma(species, "condensate", cond_obs[mask_cond_only]))
    if not residuals:
        return np.array([], dtype=float)
    return np.concatenate(residuals)


def fit_aroma_parameters(
    batches: list[base.BatchData],
    theta0: dict[str, float],
    core_cache: dict[str, pd.DataFrame],
    max_nfev: int,
    n_starts: int,
    seed: int = 88,
) -> tuple[dict[str, float], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    lb, ub = log_bounds(AROMA_TARGETS)
    summaries = []
    best_theta = dict(theta0)
    best_wsse = np.inf
    for idx in range(int(n_starts)):
        start = dict(theta0)
        if idx > 0:
            for name in AROMA_TARGETS:
                lo, hi = AROMA_BOUNDS[name]
                start[name] = float(np.clip(theta0[name] * math.exp(rng.normal(0.0, 0.65)), lo, hi))
        x0 = np.clip(log_vector(start, AROMA_TARGETS), lb + 1e-9, ub - 1e-9)

        def fun(x: np.ndarray) -> np.ndarray:
            theta = theta_from_log(x, theta0, AROMA_TARGETS)
            res = [aroma_residual_vector(theta, batches, core_cache)]
            prior_scales = [1.1 if name.startswith("alpha_") else 1.5 for name in AROMA_TARGETS]
            prior = np.asarray([math.log(theta[name] / theta0[name]) / scale for name, scale in zip(AROMA_TARGETS, prior_scales)], dtype=float)
            res.append(prior)
            return np.concatenate(res)

        start_res = fun(x0)
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
        theta_hat = theta_from_log(result.x, theta0, AROMA_TARGETS)
        end_res = fun(result.x)
        row = {
            "fit": f"pilot_aroma_multistart_{idx:02d}",
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "initial_wsse": float(np.dot(start_res, start_res)),
            "final_wsse": float(np.dot(end_res, end_res)),
            "n_residuals": int(len(end_res)),
            "wsse_per_residual": float(np.dot(end_res, end_res) / max(len(end_res), 1)),
        }
        summaries.append(row)
        if row["final_wsse"] < best_wsse:
            best_wsse = row["final_wsse"]
            best_theta = theta_hat
    return clip_extended(best_theta), pd.DataFrame(summaries)


def co2_rate_from_core(theta: dict[str, float], batch: base.BatchData | base.FutureDesign, core: pd.DataFrame, times: np.ndarray) -> np.ndarray:
    values = []
    for t in np.asarray(times, dtype=float):
        rates = joint._core_rates(theta, batch, core, float(t))
        values.append(joint.CO2_G_PER_G_ETHANOL * max(float(rates["ethanol_prod"]), 0.0))
    return np.asarray(values, dtype=float)


def downsample_co2(co2: pd.DataFrame, dt_h: float = CO2_DOWNSAMPLE_H) -> pd.DataFrame:
    if co2.empty:
        return co2.copy()
    rows = []
    for batch, group in co2.groupby("batch", sort=True):
        clean = group.dropna(subset=["time_h_effective", "co2_raw"]).copy()
        clean = clean[clean["time_h_effective"].ge(0.0)]
        if clean.empty:
            continue
        clean["bin"] = (clean["time_h_effective"] / float(dt_h)).round().astype(int)
        agg = clean.groupby("bin").agg(time_h_effective=("time_h_effective", "median"), co2_raw=("co2_raw", "median")).reset_index(drop=True)
        agg["batch"] = str(batch)
        rows.append(agg)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["batch", "time_h_effective", "co2_raw"])


def co2_shape_residual(theta: dict[str, float], batches_by_name: dict[str, base.BatchData], co2: pd.DataFrame) -> np.ndarray:
    if co2.empty:
        return np.array([], dtype=float)
    residuals = []
    for batch_name, group in co2.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        horizon = max(float(batch.time[-1]), float(np.nanmax(times)))
        sim_times = np.asarray(sorted(set(np.arange(0.0, horizon + 1e-9, 2.0).round(8)).union(set(times.round(8)))), dtype=float)
        design_like = batch
        core = base.simulate(design_like, theta, sim_times)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        pred = co2_rate_from_core(theta, design_like, core, times)
        obs = group["co2_raw"].to_numpy(dtype=float)
        mask = np.isfinite(pred) & np.isfinite(obs) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        pred = pred[mask]
        obs = obs[mask]
        denom = float(np.dot(pred, pred))
        scale = float(np.dot(obs, pred) / denom) if denom > 1e-12 else 0.0
        resid = (scale * pred - obs) / max(CO2_RATE_SIGMA, 0.1 * float(np.nanmax(obs)))
        residuals.append(resid)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def secondary_residual_with_co2(
    theta: dict[str, float],
    batches: list[base.BatchData],
    co2_downsampled: pd.DataFrame,
    core_cache: dict[str, pd.DataFrame] | None = None,
) -> np.ndarray:
    if core_cache is None:
        core_cache = joint.precompute_core_cache(batches, theta)
    res = [v2.residual_v2(theta, batches, core_cache)]
    by_name = {batch.batch: batch for batch in batches}
    res.append(co2_shape_residual(theta, by_name, co2_downsampled))
    return np.concatenate([r for r in res if len(r)])


def combined_residual_vector(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_downsampled: pd.DataFrame,
) -> np.ndarray:
    residuals = [base.residual_vector(theta, primary_batches)]
    core_cache = joint.precompute_core_cache(extended_batches, theta)
    residuals.append(secondary_residual_with_co2(theta, extended_batches, co2_downsampled, core_cache))
    residuals.append(aroma_residual_vector(theta, extended_batches, core_cache))
    return np.concatenate([r for r in residuals if len(r)])


def finite_difference_jacobian(theta: dict[str, float], parameters: tuple[str, ...], residual_fun, step: float) -> tuple[np.ndarray, np.ndarray]:
    base_res = residual_fun(theta)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        lb, ub = parameter_bounds(name)
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), lb, ub))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), lb, ub))
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
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
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace_inv": float(np.trace(cov)),
        f"{prefix}rank_1e-8": int(np.sum(eig > max_eig * 1e-8)) if np.isfinite(max_eig) and max_eig > 0.0 else 0,
    }


def fim_diagnostics(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...], label: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + fim.T)
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
    weak_rows = []
    for idx in range(min(8, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        top = np.argsort(np.abs(vec))[::-1][:8]
        weak_rows.append(
            {
                "analysis": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs(vec[i]):.3f}" for i in top),
            }
        )
    cov = stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[idx, idx], 0.0)))
        value = float(theta[name])
        lb, ub = parameter_bounds(name)
        active = value <= lb * 1.01 or value >= ub / 1.01
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
                "fim_diag": float(fim[idx, idx]),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return spectrum, pd.DataFrame(weak_rows), pd.DataFrame(rows)


def natural_initials(model_data: pd.DataFrame) -> dict[str, float]:
    first = model_data.sort_values("time_h").groupby("batch").first(numeric_only=True)
    def median_col(col: str, default: float) -> float:
        if col not in first:
            return float(default)
        values = pd.to_numeric(first[col], errors="coerce").dropna()
        return float(values.median()) if not values.empty else float(default)
    initials = {
        "X": median_col("X_viable_kg_m3", 0.45),
        "Xd": 0.0,
        "N": median_col("N_kg_m3", 0.10),
        "G": median_col("G_g_l", 65.0),
        "F": median_col("F_g_l", 73.0),
        "E": median_col("E_g_l", 7.0),
        "Gly": median_col("glycerol_g_l", 2.5),
        "Pyr": median_col("pyruvic_acid", 25.0),
        "AcAld": median_col("acetaldehyde", 170.0),
        "Acetate": 0.08,
        "O2": 6.5,
        "CO2": 0.0,
    }
    for species in joint.AROMA_SPECIES:
        total_col = joint.AROMA_COLUMNS[species]
        cond_col = AROMA_CONDENSATE_COLUMNS[species]
        total_series = first[total_col] if total_col in first.columns else pd.Series(dtype=float)
        cond_series = first[cond_col] if cond_col in first.columns else pd.Series(dtype=float)
        total = pd.to_numeric(total_series, errors="coerce").dropna()
        cond = pd.to_numeric(cond_series, errors="coerce").dropna()
        total0 = float(total.median()) if not total.empty else 0.0
        cond0 = float(cond.median()) if not cond.empty else 0.0
        initials[f"{species}_liq"] = max(total0 - cond0, 0.0)
        initials[f"{species}_cond"] = max(cond0, 0.0)
    return initials


def snap_pulses(horizon: float, rows: Iterable[tuple[float, float]]) -> tuple[tuple[float, float], ...]:
    out = []
    for time_h, amount in rows:
        out.append((base.nearest_operational_time(float(time_h), float(horizon)), float(amount)))
    return tuple(sorted(out))


def natural_candidate_designs(model_data: pd.DataFrame) -> dict[str, base.FutureDesign]:
    initials = natural_initials(model_data)
    h = NATURAL_DESIGN_HORIZON_H
    zero = {channel: tuple() for channel in base.INPUT_CHANNELS}
    specs = [
        ("natural_pilot_reference_18C", "reference", (18.0, 18.0, 18.0, 18.0), (), "Natural must reference at moderate temperature."),
        ("natural_pilot_reference_20C", "reference", (20.0, 20.0, 20.0, 20.0), (), "Natural must warmer reference to increase rate and CO2 information."),
        ("natural_pilot_cold_to_warm_earlyN", "temperature_N", (14.0, 16.0, 22.0, 20.0), ((30.0, 0.045),), "Cold start followed by warm transition and early nitrogen pulse."),
        ("natural_pilot_warm_to_cool_noN", "temperature", (22.0, 22.0, 17.0, 16.0), (), "Warm early phase to excite growth and CO2, then cool aroma-retention phase."),
        ("natural_pilot_midN_temperature_step", "temperature_N", (16.0, 20.0, 23.0, 18.0), ((54.0, 0.040),), "Temperature step with mid-growth nitrogen perturbation."),
        ("natural_pilot_lateN_stationary_probe", "N_timing", (18.0, 20.0, 20.0, 18.0), ((84.0, 0.050),), "Late nitrogen addition to test stationary/growth split in secondary and aroma formation."),
        ("natural_pilot_two_step_N_ladder", "N_timing", (17.0, 19.0, 21.0, 18.0), ((30.0, 0.030), (78.0, 0.035)), "Two smaller nitrogen pulses to separate early growth and later metabolic response."),
        ("natural_pilot_high_rate_strip", "co2_aroma", (21.0, 24.0, 23.0, 19.0), ((30.0, 0.030),), "High-rate natural fermentation to excite CO2 stripping and aroma loss directions."),
        ("natural_pilot_low_temp_aroma_retention", "aroma_retention", (13.0, 15.0, 16.0, 16.0), ((54.0, 0.035),), "Cold profile to contrast aroma retention against high-rate stripping."),
        ("natural_pilot_noN_dynamic_temperature", "temperature", (15.0, 22.0, 18.0, 22.0), (), "Temperature-only perturbation for settings where nutrient action is constrained."),
    ]
    designs = {}
    for name, family, temps, n_pulses, rationale in specs:
        pulses = dict(zero)
        pulses["N"] = snap_pulses(h, n_pulses)
        designs[name] = base.FutureDesign(
            name=name,
            family=family,
            medium="natural",
            horizon_h=h,
            initials=dict(initials),
            temperature_segments=tuple(float(t) for t in temps),
            pulses={channel: tuple(rows) for channel, rows in pulses.items()},
            rationale=rationale,
        )
    return designs


def future_sample_times(design: base.FutureDesign) -> np.ndarray:
    return base.operational_sample_times(design.horizon_h, policy="balanced")


def future_co2_times(design: base.FutureDesign) -> np.ndarray:
    return np.arange(0.0, float(design.horizon_h) + 1e-9, 6.0)


def future_residual_vector(theta: dict[str, float], design: base.FutureDesign) -> np.ndarray:
    sample_times = future_sample_times(design)
    co2_times = future_co2_times(design)
    all_times = np.asarray(sorted(set(sample_times).union(set(co2_times)).union({0.0, float(design.horizon_h)})), dtype=float)
    core = base.simulate(design, theta, all_times)
    if core is None:
        return np.ones(1000, dtype=float) * 1e6
    secondary = v2.integrate_secondary_v2(design, theta, core)
    try:
        aromas_liq, aromas_cond = integrate_aroma_pilot_euler(design, theta, core)
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
        for state in ("Pyr", "AcAld", "Acetate"):
            residuals.append(sec_sample[state].to_numpy(dtype=float) / joint.SIGMA[state])
        aroma_sample = aromas_liq.loc[valid_samples]
        aroma_cond_sample = aromas_cond.loc[valid_samples]
        for species in joint.AROMA_SPECIES:
            liq_center = aroma_sample[species].to_numpy(dtype=float)
            cond_center = aroma_cond_sample[species].to_numpy(dtype=float)
            residuals.append(liq_center / aroma_sigma(species, "liquid", liq_center))
            residuals.append(cond_center / aroma_sigma(species, "condensate", cond_center))
    valid_co2 = [float(t) for t in co2_times if float(t) in core.index]
    if valid_co2:
        rate = co2_rate_from_core(theta, design, core, np.asarray(valid_co2, dtype=float))
        sigma = np.maximum(0.05, 0.20 * np.maximum(rate, 0.05))
        residuals.append(rate / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def candidate_fim(theta: dict[str, float], design: base.FutureDesign, step: float) -> np.ndarray:
    jac, _ = finite_difference_jacobian(theta, TARGET_PARAMETERS, lambda th: future_residual_vector(th, design), step)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def score_fim(fim: np.ndarray, objective: str = "hybrid") -> float:
    metrics = fim_metrics(fim)
    if objective == "d_opt":
        return float(metrics["logdet"])
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-18)
    return float(metrics["logdet"]) + 2.0 * math.log(min_rel) - 0.05 * math.log(max(float(metrics["trace_inv"]), 1e-18))


def variance_reduction(prior: np.ndarray, combined: np.ndarray) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    post_cov = stable_inverse(combined)
    rows = {}
    reductions = []
    for idx, name in enumerate(TARGET_PARAMETERS):
        before = float(prior_cov[idx, idx])
        after = float(post_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        if name in SECONDARY_TARGETS or name in AROMA_TARGETS:
            reductions.append(rows[f"var_reduction_{name}"])
    rows["secondary_aroma_mean_var_reduction"] = float(np.nanmean(reductions)) if reductions else np.nan
    rows["secondary_aroma_worst_var_reduction"] = float(np.nanmin(reductions)) if reductions else np.nan
    return rows


def greedy_select(
    candidate_fims: dict[str, np.ndarray],
    designs: dict[str, base.FutureDesign],
    prior: np.ndarray,
    campaign_size: int,
    objective: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    selected: list[str] = []
    remaining = set(candidate_fims)
    current = prior.copy()
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + candidate_fims[name]
            score = score_fim(trial, objective=objective)
            if score > best_score:
                best = name
                best_score = score
                best_metrics = {**fim_metrics(trial, prefix="campaign_"), **variance_reduction(prior, trial)}
        if best is None:
            break
        current = current + candidate_fims[best]
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
                "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                "score": float(best_score),
                "rationale": design.rationale,
                **(best_metrics or {}),
            }
        )
    return pd.DataFrame(rows), current


def rank_candidates(candidate_fims: dict[str, np.ndarray], designs: dict[str, base.FutureDesign], prior: np.ndarray) -> pd.DataFrame:
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
                **variance_reduction(prior, combined),
                "rationale": design.rationale,
            }
        )
    return pd.DataFrame(rows).sort_values("combined_logdet", ascending=False).reset_index(drop=True)


def plot_fit(theta: dict[str, float], batches: list[base.BatchData], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        time = joint._dense_time_grid(batch, batch.time)
        core = base.simulate(batch, theta, time)
        if core is None:
            continue
        sec = v2.integrate_secondary_v2(batch, theta, core)
        try:
            aromas, aroma_cond = integrate_aroma_pilot_euler(batch, theta, core)
        except Exception:
            aromas = pd.DataFrame(index=core.index)
            aroma_cond = pd.DataFrame(index=core.index)
        fig, axes = plt.subplots(6, 2, figsize=(14, 18), sharex=True)
        panels = [
            ("core", "G", "Glucose"),
            ("core", "F", "Fructose"),
            ("core", "N", "YAN"),
            ("core", "X", "Viable biomass"),
            ("core", "E", "Ethanol"),
            ("core", "Gly", "Glycerol"),
            ("sec", "Pyr", "Pyruvic acid"),
            ("sec", "AcAld", "Acetaldehyde"),
            ("aroma_liq", "ethyl_acetate", "Ethyl acetate retained"),
            ("aroma_cond", "ethyl_acetate", "Ethyl acetate condensate"),
            ("aroma_liq", "isoamyl_acetate", "Isoamyl acetate retained"),
            ("aroma_cond", "isoamyl_acetate", "Isoamyl acetate condensate"),
        ]
        for ax, (kind, state, title) in zip(axes.ravel(), panels):
            if kind == "core":
                ax.plot(core.index, core[state], lw=1.6)
                obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            elif kind == "sec":
                ax.plot(sec.index, sec[state], lw=1.6)
                obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            elif kind == "aroma_liq":
                if state in aromas:
                    ax.plot(aromas.index, aromas[state], lw=1.6)
                column = joint.AROMA_COLUMNS[state]
                total_obs = np.asarray(batch.observations.get(column, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
                cond_obs = np.asarray(batch.observations.get(AROMA_CONDENSATE_COLUMNS[state], np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
                obs = np.where(np.isfinite(total_obs) & np.isfinite(cond_obs), np.maximum(total_obs - cond_obs, 0.0), total_obs)
            else:
                if state in aroma_cond:
                    ax.plot(aroma_cond.index, aroma_cond[state], lw=1.6)
                obs = np.asarray(batch.observations.get(AROMA_CONDENSATE_COLUMNS[state], np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if mask.any():
                ax.scatter(batch.time[mask], obs[mask], s=20, color="black", zorder=3)
            ax.set_title(title)
            ax.grid(True, alpha=0.25)
        axes[-1, 0].set_xlabel("time (h)")
        axes[-1, 1].set_xlabel("time (h)")
        fig.suptitle(f"Pilot calibration fit: {batch.batch}", y=0.995)
        fig.tight_layout()
        fig.savefig(output_dir / f"pilot_fit_{batch.batch}.png", dpi=160)
        plt.close(fig)


def plot_candidate_inputs(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in selected["candidate"].drop_duplicates().astype(str):
        design = designs[name]
        time = np.linspace(0.0, float(design.horizon_h), 400)
        fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
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
        fig.suptitle(f"{design.name}: {design.rationale}", fontsize=10, y=0.99)
        fig.tight_layout()
        fig.savefig(output_dir / f"candidate_inputs_{design.name}.png", dpi=160)
        plt.close(fig)


def write_report(
    fit_summary: pd.DataFrame,
    estimability: pd.DataFrame,
    weak: pd.DataFrame,
    ranking: pd.DataFrame,
    selected: pd.DataFrame,
    co2_decisions: pd.DataFrame,
) -> None:
    path = RESULTS_DIR / "pilot_2025_calibration_estimability_report.md"
    with path.open("w", encoding="utf-8") as f:
        f.write("# Pilot 2025 calibration and estimability\n\n")
        f.write("## Scope\n\n")
        f.write(
            "This run reuses the reduced extended fermentation model for pilot-scale natural-must data: "
            "primary fermentation, glycerol, secondary v2 reduced chemistry, liquid aroma synthesis, and aroma volatilization to condenser. "
            "The future design space is restricted to natural must, temperature setpoints, and nutrient additions.\n\n"
        )
        f.write("## Aroma mass-balance interpretation\n\n")
        f.write(
            "`*_total` is interpreted as retained-in-wine plus accumulated condenser equivalent concentration. "
            "`*_condensate` is interpreted as the accumulated condenser equivalent concentration. "
            "When both are observed, the retained wine observation is reconstructed as `total - condensate` and fitted against the liquid aroma state. "
            "The condenser observation is fitted against the cumulative volatilized/captured state, making `alpha_*_loss` estimable as effective volatilization/capture coefficients.\n\n"
        )
        f.write("## CO2 curation\n\n")
        f.write(co2_decisions.to_markdown(index=False) if not co2_decisions.empty else "_No CO2 decisions available._")
        f.write("\n\n")
        f.write("## Calibration stages\n\n")
        f.write(fit_summary.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Practical estimability from current pilot data\n\n")
        cols = ["parameter", "theta", "std_log_approx", "approx_95_multiplier", "active_bound", "classification"]
        f.write(estimability[cols].to_markdown(index=False))
        f.write("\n\n")
        f.write("## Weak FIM directions\n\n")
        f.write(weak.to_markdown(index=False) if not weak.empty else "_No weak directions table._")
        f.write("\n\n")
        f.write("## Natural-must candidate ranking\n\n")
        show_cols = [
            "candidate",
            "family",
            "combined_logdet",
            "combined_min_relative_eigenvalue",
            "secondary_aroma_mean_var_reduction",
            "secondary_aroma_worst_var_reduction",
            "N_pulses_kg_m3",
            "rationale",
        ]
        f.write(ranking[show_cols].head(10).to_markdown(index=False) if not ranking.empty else "_No ranking generated._")
        f.write("\n\n")
        f.write("## Selected campaign, hybrid criterion\n\n")
        f.write(selected[show_cols[:2] + ["campaign_order", "campaign_logdet", "campaign_min_relative_eigenvalue", "N_pulses_kg_m3", "rationale"]].to_markdown(index=False) if not selected.empty else "_No selected campaign generated._")
        f.write("\n\n")
        f.write("## Interpretation notes\n\n")
        f.write("- `alpha_*_loss` parameters are now estimated because condenser-equivalent aroma observations are available.\n")
        f.write("- CO2 is used as a rate-shape signal with a per-batch scale factor, not as an absolute concentration measurement.\n")
        f.write("- Candidate designs do not include glucose/fructose/ethanol/biomass injections because this pilot setting is treated as natural-must constrained.\n")


def create_notebook() -> None:
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Pilot 2025 Calibration, Estimability, and Natural-Must MBDoE\n\n"
            "This notebook applies the reduced extended fermentation model to the pilot-scale natural-must dataset. "
            "It mirrors the laboratory workflow, but the design space is intentionally narrower: natural must is assumed, and the realistic manipulated inputs are temperature setpoints and nutrient additions.\n\n"
            "The model combines primary fermentation, glycerol production, the reduced secondary v2 model, and empirical aroma synthesis for ethyl acetate, isoamyl acetate, and ethyl octanoate."
        ),
        nbformat.v4.new_markdown_cell(
            "## Mathematical Structure\n\n"
            "The primary state vector is\n\n"
            "$$x_p = [X, X_d, N, G, F, E, Gly]^T$$\n\n"
            "where `X` is viable biomass, `X_d` is dead biomass, `N` is assimilable nitrogen, `G` and `F` are glucose and fructose, `E` is ethanol, and `Gly` is glycerol.\n\n"
            "The secondary model uses\n\n"
            "$$x_s = [Pyr, AcAld, Acetate, O_2, CO_2]^T$$\n\n"
            "with the reduced v2 free set\n\n"
            "$$\\theta_s = \\{k_{PyrS,N}, k_{PyrO2}, k_{PyrDrain}, k_{AldS,N}, k_{AldRed}, k_{AcAld}, k_{AcStress}\\}.$$\n\n"
            "Aroma synthesis is empirical and rate-dependent:\n\n"
            "$$r_i = \\left(k_{i,growth}\\phi_N + k_{i,stationary}(1-\\phi_N)\\right)q_S,$$\n\n"
            "where `i` is ethyl acetate, isoamyl acetate, or ethyl octanoate; `q_S` is the total sugar uptake rate; and `\\phi_N` is the nitrogen-growth phase proxy.\n\n"
            "Volatilization is represented as an effective liquid-to-condenser transfer:\n\n"
            "$$r_{loss,i} = \\alpha_i K_i(T,E,S) q_{CO2} C_{L,i}.$$\n\n"
            "The observation model uses\n\n"
            "$$C_{wine,i}^{obs} = C_{total,i}^{obs} - C_{cond,i}^{obs}, \\qquad C_{cond,i}^{obs} = C_{cond,i}.$$\n\n"
            "Therefore the pilot condenser data make the `\\alpha_i` loss parameters estimable as effective volatilization/capture coefficients."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "cwd = Path.cwd()\n"
            "if (cwd / 'results' / 'calibration_estimability').exists():\n"
            "    RESULTS = cwd / 'results' / 'calibration_estimability'\n"
            "else:\n"
            "    RESULTS = cwd / 'fermentation_model' / 'pilot_2025' / 'results' / 'calibration_estimability'\n"
            "fit_summary = pd.read_csv(RESULTS / 'fit_summary.csv')\n"
            "theta = pd.read_csv(RESULTS / 'theta_pilot_extended.csv', index_col=0)\n"
            "estimability = pd.read_csv(RESULTS / 'parameter_estimability_current.csv')\n"
            "weak = pd.read_csv(RESULTS / 'weak_directions_current.csv')\n"
            "ranking = pd.read_csv(RESULTS / 'candidate_ranking_natural.csv')\n"
            "selected = pd.read_csv(RESULTS / 'selected_campaign_hybrid.csv')\n"
            "fit_summary"
        ),
        nbformat.v4.new_markdown_cell("## Calibrated Parameter Vector"),
        nbformat.v4.new_code_cell("theta"),
        nbformat.v4.new_markdown_cell(
            "## Current-Data Estimability\n\n"
            "The FIM is computed from log-parameter finite differences of the full residual vector. "
            "The approximate standard deviation is therefore in log-parameter space. "
            "Large `approx_95_multiplier` values indicate broad practical uncertainty. "
            "Parameters on active bounds are not considered reliable even when the local curvature appears high."
        ),
        nbformat.v4.new_code_cell("estimability.sort_values('std_log_approx', ascending=False)"),
        nbformat.v4.new_markdown_cell(
            "## Weak Eigen-Directions\n\n"
            "Weak directions identify combinations of parameters that the current pilot data do not separate well. "
            "These directions are more informative than single-parameter rankings when parameters are correlated."
        ),
        nbformat.v4.new_code_cell("weak"),
        nbformat.v4.new_markdown_cell(
            "## Natural-Must Candidate Design Ranking\n\n"
            "Candidate experiments are restricted to natural must. "
            "The manipulated variables are temperature setpoint profiles and nitrogen pulse timing/dose. "
            "No glucose, fructose, ethanol, or biomass injections are included in this pilot-scale design library."
        ),
        nbformat.v4.new_code_cell(
            "cols = ['candidate', 'family', 'combined_logdet', 'combined_min_relative_eigenvalue', "
            "'secondary_aroma_mean_var_reduction', 'secondary_aroma_worst_var_reduction', 'N_pulses_kg_m3', 'rationale']\n"
            "ranking[cols].head(10)"
        ),
        nbformat.v4.new_markdown_cell(
            "## Selected Hybrid Campaign\n\n"
            "The hybrid objective keeps D-optimality as the information-volume term and penalizes designs that leave very weak eigen-directions. "
            "This is the same practical logic used in the laboratory design fork."
        ),
        nbformat.v4.new_code_cell(
            "selected[['campaign_order', 'candidate', 'family', 'campaign_logdet', 'campaign_min_relative_eigenvalue', "
            "'secondary_aroma_mean_var_reduction', 'secondary_aroma_worst_var_reduction', 'N_pulses_kg_m3', 'rationale']]"
        ),
        nbformat.v4.new_markdown_cell(
            "## Generated Plots\n\n"
            "Calibration fit plots are saved in `results/calibration_estimability/plots/fits`. "
            "Selected design input plots are saved in `results/calibration_estimability/plots/designs`."
        ),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-multistarts", type=int, default=4)
    parser.add_argument("--aroma-multistarts", type=int, default=4)
    parser.add_argument("--max-core-nfev", type=int, default=140)
    parser.add_argument("--max-secondary-nfev", type=int, default=140)
    parser.add_argument("--max-aroma-nfev", type=int, default=180)
    parser.add_argument("--sensitivity-step", type=float, default=1e-2)
    parser.add_argument("--campaign-size", type=int, default=6)
    parser.add_argument("--skip-candidates", action="store_true")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = RESULTS_DIR / "plots"
    (plot_dir / "fits").mkdir(parents=True, exist_ok=True)
    (plot_dir / "designs").mkdir(parents=True, exist_ok=True)

    model_data, co2_curated, co2_decisions = load_pilot_model_data()
    co2_down = downsample_co2(co2_curated)
    model_data.to_csv(RESULTS_DIR / "pilot_2025_model_data.csv", index=False)
    co2_down.to_csv(RESULTS_DIR / "pilot_2025_co2_downsampled_for_fit.csv", index=False)
    co2_decisions.to_csv(RESULTS_DIR / "pilot_2025_co2_curation_decisions.csv", index=False)

    primary_batches = make_primary_batches(model_data)
    extended_batches = make_extended_batches(model_data)
    theta0 = load_initial_theta()
    pd.Series(theta0).to_csv(RESULTS_DIR / "theta_initial.csv")

    print("[fit] core", flush=True)
    theta_core, core_summary = fit_core_multistart(primary_batches, theta0, args.core_multistarts, args.max_core_nfev)
    pd.Series(theta_core).to_csv(RESULTS_DIR / "theta_core.csv")
    core_summary.to_csv(RESULTS_DIR / "fit_core_multistart.csv", index=False)

    print("[fit] secondary v2", flush=True)
    theta_sec, sec_summary, core_cache = fit_secondary_v2(extended_batches, theta_core, args.max_secondary_nfev)
    # Rebuild cache at the fitted core/secondary vector. Secondary parameters do not affect core, but this keeps one source of truth.
    core_cache = joint.precompute_core_cache(extended_batches, theta_sec)
    pd.Series(theta_sec).to_csv(RESULTS_DIR / "theta_secondary_v2.csv")

    print("[fit] aroma", flush=True)
    theta_hat, aroma_summary = fit_aroma_parameters(extended_batches, theta_sec, core_cache, args.max_aroma_nfev, args.aroma_multistarts)
    pd.Series(theta_hat).to_csv(RESULTS_DIR / "theta_pilot_extended.csv")

    fit_summary = pd.concat([core_summary, sec_summary.rename(columns={"model": "fit"}), aroma_summary], ignore_index=True, sort=False)
    fit_summary.to_csv(RESULTS_DIR / "fit_summary.csv", index=False)

    print("[fim] current data", flush=True)
    jac, resid = finite_difference_jacobian(
        theta_hat,
        TARGET_PARAMETERS,
        lambda th: combined_residual_vector(th, primary_batches, extended_batches, co2_down),
        args.sensitivity_step,
    )
    fim = 0.5 * ((jac.T @ jac) + (jac.T @ jac).T)
    pd.DataFrame(jac, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "jacobian_current.csv", index=False)
    pd.Series(resid).to_csv(RESULTS_DIR / "residual_current.csv", index=False)
    pd.DataFrame(fim, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "fim_current.csv")
    spectrum, weak, estimability = fim_diagnostics(fim, theta_hat, TARGET_PARAMETERS, "current_pilot")
    spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_current.csv", index=False)
    weak.to_csv(RESULTS_DIR / "weak_directions_current.csv", index=False)
    estimability.to_csv(RESULTS_DIR / "parameter_estimability_current.csv", index=False)

    plot_fit(theta_hat, extended_batches, plot_dir / "fits")

    ranking = pd.DataFrame()
    selected = pd.DataFrame()
    if not args.skip_candidates:
        print("[doe] natural candidates", flush=True)
        designs = natural_candidate_designs(model_data)
        design_rows = []
        for design in designs.values():
            design_rows.append(
                {
                    "candidate": design.name,
                    "family": design.family,
                    "medium": design.medium,
                    "horizon_h": design.horizon_h,
                    "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                    "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                    "initials_json": json.dumps(design.initials, sort_keys=True),
                    "rationale": design.rationale,
                }
            )
        pd.DataFrame(design_rows).to_csv(RESULTS_DIR / "candidate_design_library_natural.csv", index=False)

        candidate_fims = {}
        for idx, (name, design) in enumerate(designs.items(), start=1):
            print(f"[doe] candidate {idx}/{len(designs)} {name}", flush=True)
            cfim = candidate_fim(theta_hat, design, args.sensitivity_step)
            candidate_fims[name] = cfim
            pd.DataFrame(cfim, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / f"candidate_fim_{name}.csv")
        ranking = rank_candidates(candidate_fims, designs, fim)
        ranking.to_csv(RESULTS_DIR / "candidate_ranking_natural.csv", index=False)
        selected, campaign_fim = greedy_select(candidate_fims, designs, fim, args.campaign_size, objective="hybrid")
        selected.to_csv(RESULTS_DIR / "selected_campaign_hybrid.csv", index=False)
        pd.DataFrame(campaign_fim, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "fim_current_plus_selected_hybrid.csv")
        _spec, weak_after, estim_after = fim_diagnostics(campaign_fim, theta_hat, TARGET_PARAMETERS, "current_plus_selected_hybrid")
        weak_after.to_csv(RESULTS_DIR / "weak_directions_current_plus_selected_hybrid.csv", index=False)
        estim_after.to_csv(RESULTS_DIR / "parameter_estimability_current_plus_selected_hybrid.csv", index=False)
        plot_candidate_inputs(selected, designs, plot_dir / "designs")

    write_report(fit_summary, estimability, weak, ranking, selected, co2_decisions)
    create_notebook()
    print(f"[done] report: {RESULTS_DIR / 'pilot_2025_calibration_estimability_report.md'}", flush=True)
    print(f"[done] notebook: {NOTEBOOK_PATH}", flush=True)


if __name__ == "__main__":
    main()
