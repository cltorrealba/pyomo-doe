from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

try:
    from pyomo.contrib.doe.utils import rescale_FIM
except Exception:  # pragma: no cover - optional Pyomo DoE utility
    rescale_FIM = None


SCRIPT_DIR = Path(__file__).resolve().parent
PILOT_DIR = SCRIPT_DIR.parent if SCRIPT_DIR.name == "support" else SCRIPT_DIR
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(PILOT_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_DIR))
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

import run_pilot_2025_aroma_model_selection_doe as aroma_sel
import run_pilot_2025_calibration_estimability as pilot
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = PILOT_DIR / "results" / os.environ.get("PILOT_GLOBAL_RESULTS_DIR", "global_state_model_selection_doe")
NOTEBOOK_PATH = PILOT_DIR / os.environ.get("PILOT_GLOBAL_NOTEBOOK", "pilot_2025_global_model_selection_doe.ipynb")
EXECUTED_NOTEBOOK_PATH = PILOT_DIR / os.environ.get(
    "PILOT_GLOBAL_EXECUTED_NOTEBOOK",
    "pilot_2025_global_model_selection_doe.executed.ipynb",
)

SECONDARY_OPTIONAL_TERMS = (
    "kPyrS_stat",
    "kAldPyr",
    "kAldS_stat",
    "kAldO2",
    "kAcAssim",
)
STATE_GROUPS = {
    "core": tuple(base.STATE_NAMES),
    "secondary": ("Pyr", "AcAld", "Acetate", "O2"),
    "co2": ("CO2_rate_shape",),
}


@dataclass(frozen=True)
class SecondaryVariant:
    name: str
    parameters: tuple[str, ...]
    description: str
    equation_markdown: str
    fixed_overrides: dict[str, float]


def secondary_variant_library() -> dict[str, SecondaryVariant]:
    reduced = tuple(v2.V2_REDUCED_PARAMETERS)

    def fixed_for(parameters: tuple[str, ...]) -> dict[str, float]:
        return {name: v2.V2_BOUNDS[name][0] for name in SECONDARY_OPTIONAL_TERMS if name not in parameters}

    variants = {
        "secondary_reduced_o2fixed": SecondaryVariant(
            name="secondary_reduced_o2fixed",
            parameters=reduced,
            fixed_overrides=fixed_for(reduced),
            description=(
                "Current reduced v2 structure: pyruvate and acetaldehyde production are driven by "
                "nitrogen-associated sugar uptake and oxygen gating; acetate has acetaldehyde oxidation "
                "and ethanol-stress production. Optional stationary, pyruvate-to-acetaldehyde, O2-to-acetaldehyde, "
                "and acetate-assimilation terms are fixed at their lower bounds."
            ),
            equation_markdown=(
                r"$$\frac{dPyr}{dt}=k_{PyrS,N}\phi_N q_S+k_{PyrO2}g_{O2}X"
                r"-k_{PyrDrain}Pyr X(0.25+\phi_{stat}+0.5\phi_E)$$" "\n\n"
                r"$$\frac{dAcAld}{dt}=k_{AldS,N}\phi_N q_S"
                r"-k_{AldRed}AcAld X(\phi_{ana}+0.25\phi_{stat})-k_{AcAld}AcAld Xg_{O2}$$" "\n\n"
                r"$$\frac{dAcetate}{dt}=\frac{k_{AcAld}AcAld Xg_{O2}}{1000}+k_{AcStress}X\phi_E$$"
            ),
        ),
        "secondary_phase_split": SecondaryVariant(
            name="secondary_phase_split",
            parameters=tuple(dict.fromkeys(reduced + ("kPyrS_stat", "kAldS_stat"))),
            fixed_overrides=fixed_for(tuple(dict.fromkeys(reduced + ("kPyrS_stat", "kAldS_stat")))),
            description=(
                "Adds stationary-phase production terms for pyruvate and acetaldehyde. This tests whether "
                "late fermentation metabolite accumulation is a phase effect rather than a pure estimability issue."
            ),
            equation_markdown=(
                r"$$\frac{dPyr}{dt}=(k_{PyrS,N}\phi_N+k_{PyrS,stat}\phi_{stat})q_S"
                r"+k_{PyrO2}g_{O2}X-k_{PyrDrain}Pyr X(0.25+\phi_{stat}+0.5\phi_E)$$" "\n\n"
                r"$$\frac{dAcAld}{dt}=(k_{AldS,N}\phi_N+k_{AldS,stat}\phi_{stat})q_S-\cdots$$"
            ),
        ),
        "secondary_redox_o2": SecondaryVariant(
            name="secondary_redox_o2",
            parameters=tuple(dict.fromkeys(reduced + ("kAldPyr", "kAldO2"))),
            fixed_overrides=fixed_for(tuple(dict.fromkeys(reduced + ("kAldPyr", "kAldO2")))),
            description=(
                "Adds a pyruvate-to-acetaldehyde source and an oxygen-linked acetaldehyde source. This tests "
                "whether acetaldehyde dynamics need an explicit redox/oxygen proxy."
            ),
            equation_markdown=(
                r"$$\frac{dAcAld}{dt}=k_{AldPyr}Pyr X+k_{AldS,N}\phi_Nq_S+k_{AldO2}g_{O2}X"
                r"-k_{AldRed}AcAld X(\phi_{ana}+0.25\phi_{stat})-k_{AcAld}AcAld Xg_{O2}$$"
            ),
        ),
        "secondary_acetate_assimilation": SecondaryVariant(
            name="secondary_acetate_assimilation",
            parameters=tuple(dict.fromkeys(reduced + ("kAcAssim",))),
            fixed_overrides=fixed_for(tuple(dict.fromkeys(reduced + ("kAcAssim",)))),
            description=(
                "Adds acetate assimilation during nitrogen-associated growth. This tests whether acetate bias "
                "is explainable by a sink term rather than by adding new states."
            ),
            equation_markdown=(
                r"$$\frac{dAcetate}{dt}=\frac{k_{AcAld}AcAld Xg_{O2}}{1000}+k_{AcStress}X\phi_E"
                r"-k_{AcAssim}Acetate X\phi_N$$"
            ),
        ),
        "secondary_full_chem_o2fixed": SecondaryVariant(
            name="secondary_full_chem_o2fixed",
            parameters=tuple(v2.V2_CHEM_PARAMETERS),
            fixed_overrides={},
            description=(
                "Full chemical secondary layer with O2 transfer parameters fixed. This is the upper-complexity "
                "check: it is accepted only if its improved residuals justify the extra degrees of freedom."
            ),
            equation_markdown=(
                r"$$\frac{dPyr}{dt}=(k_{PyrS,N}\phi_N+k_{PyrS,stat}\phi_{stat})q_S+k_{PyrO2}g_{O2}X"
                r"-k_{PyrDrain}Pyr X(0.25+\phi_{stat}+0.5\phi_E)$$" "\n\n"
                r"$$\frac{dAcAld}{dt}=k_{AldPyr}PyrX+(k_{AldS,N}\phi_N+k_{AldS,stat}\phi_{stat})q_S"
                r"+k_{AldO2}g_{O2}X-k_{AldRed}AcAldX(\phi_{ana}+0.25\phi_{stat})-k_{AcAld}AcAldXg_{O2}$$" "\n\n"
                r"$$\frac{dAcetate}{dt}=\frac{k_{AcAld}AcAldXg_{O2}}{1000}+k_{AcStress}X\phi_E"
                r"-k_{AcAssim}AcetateX\phi_N$$"
            ),
        ),
    }
    return variants


def load_selected_aroma_context() -> tuple[dict[str, float], str, aroma_sel.AromaVariant]:
    theta = pilot.load_initial_theta()
    theta.update(v2.default_theta_v2(theta))
    for name, value in aroma_sel.EXTRA_DEFAULTS.items():
        theta.setdefault(name, value)

    selected_name = "baseline_phase"
    selected_model_path = aroma_sel.RESULTS_DIR / "selected_model.txt"
    if selected_model_path.exists():
        selected_name = selected_model_path.read_text(encoding="utf-8").strip()

    selected_theta_path = aroma_sel.RESULTS_DIR / "theta_selected_aroma_model.csv"
    if selected_theta_path.exists():
        loaded = pd.read_csv(selected_theta_path, index_col=0).iloc[:, 0].to_dict()
        theta.update({str(k): float(v) for k, v in loaded.items() if np.isfinite(float(v))})

    variants = aroma_sel.variant_library()
    if selected_name not in variants:
        selected_name = "baseline_phase"
    theta = aroma_sel.clip_theta(theta)
    theta.update(v2.clip_v2({**v2.default_theta_v2(theta), **theta}))
    return theta, selected_name, variants[selected_name]


def state_scale(obs: np.ndarray, floor: float = 1e-9) -> float:
    obs = np.asarray(obs, dtype=float)
    obs = obs[np.isfinite(obs)]
    if obs.size == 0:
        return float(floor)
    median_abs = float(np.nanmedian(np.abs(obs)))
    obs_range = float(np.nanmax(obs) - np.nanmin(obs)) if obs.size > 1 else 0.0
    return max(median_abs, 0.25 * obs_range, float(floor))


def summarize_prediction_rows(
    pred: pd.DataFrame,
    group_cols: list[str],
    group_name: str,
    scale_floor: float = 1e-9,
) -> pd.DataFrame:
    rows = []
    if pred.empty:
        return pd.DataFrame()
    for keys, group in pred.groupby(group_cols, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        obs = group["obs"].to_numpy(dtype=float)
        yhat = group["pred"].to_numpy(dtype=float)
        residual = yhat - obs
        scale = state_scale(obs, scale_floor)
        corr = np.nan
        if len(group) >= 4 and np.nanstd(obs) > 1e-12 and np.nanstd(yhat) > 1e-12:
            corr = float(np.corrcoef(obs, yhat)[0, 1])
        row = {
            "group": group_name,
            "n": int(len(group)),
            "rmse": float(np.sqrt(np.nanmean(residual * residual))),
            "mae": float(np.nanmean(np.abs(residual))),
            "bias": float(np.nanmean(residual)),
            "obs_scale": float(scale),
            "median_obs": float(np.nanmedian(obs)),
            "median_pred": float(np.nanmedian(yhat)),
            "relative_rmse": float(np.sqrt(np.nanmean(residual * residual)) / scale),
            "relative_bias": float(np.nanmean(residual) / scale),
            "corr": corr,
        }
        for col, key in zip(group_cols, keys):
            row[col] = key
        rows.append(row)
    return pd.DataFrame(rows)


def classify_adequacy(row: pd.Series) -> tuple[str, str]:
    n = int(row.get("n", 0))
    rel_rmse = float(row.get("relative_rmse", np.nan))
    rel_bias = float(row.get("relative_bias", np.nan))
    corr = float(row.get("corr", np.nan))
    if n < 3:
        return "limited_data", "fewer than three observations"
    if np.isfinite(corr) and corr < 0.25 and rel_rmse > 0.50:
        return "structure_warning", "poor shape correlation and high normalized error"
    if rel_rmse > 1.00:
        return "structure_warning", "normalized RMSE exceeds the observation scale"
    if rel_rmse > 0.65 and abs(rel_bias) > 0.35:
        return "structure_warning", "large normalized error with systematic bias"
    if abs(rel_bias) > 0.70:
        return "bias_warning", "large systematic bias"
    if rel_rmse > 0.65:
        return "weak_fit", "high normalized error but limited systematic bias"
    return "adequate", "error is within the pragmatic adequacy thresholds"


def add_adequacy_flags(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return metrics
    out = metrics.copy()
    flags = [classify_adequacy(row) for _, row in out.iterrows()]
    out["adequacy_flag"] = [flag for flag, _reason in flags]
    out["adequacy_reason"] = [reason for _flag, reason in flags]
    return out


def core_prediction_rows(theta: dict[str, float], batches: list[base.BatchData]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for batch in batches:
        sim = base.simulate(batch, theta, batch.time)
        if sim is None:
            continue
        for state in base.STATE_NAMES:
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan)), dtype=float)
            pred = sim.loc[batch.time, state].to_numpy(dtype=float)
            mask = np.isfinite(obs) & np.isfinite(pred)
            for t, y, yhat in zip(batch.time[mask], obs[mask], pred[mask]):
                rows.append({"batch": batch.batch, "time_h": float(t), "state": state, "obs": float(y), "pred": float(yhat)})
    return pd.DataFrame(rows)


def secondary_prediction_rows(
    theta: dict[str, float],
    batches: list[base.BatchData],
    core_cache: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            continue
        sec = v2.integrate_secondary_v2(batch, theta, core)
        for state in ("Pyr", "AcAld", "Acetate", "O2"):
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan)), dtype=float)
            pred = sec.loc[batch.time, state].to_numpy(dtype=float)
            mask = np.isfinite(obs) & np.isfinite(pred)
            for t, y, yhat in zip(batch.time[mask], obs[mask], pred[mask]):
                rows.append({"batch": batch.batch, "time_h": float(t), "state": state, "obs": float(y), "pred": float(yhat)})
    return pd.DataFrame(rows)


def co2_shape_metrics(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_downsampled: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    if co2_downsampled.empty:
        return pd.DataFrame(), pd.DataFrame()
    for batch_name, group in co2_downsampled.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        horizon = max(float(np.nanmax(batch.time)), float(np.nanmax(times)))
        sim_times = np.asarray(sorted(set(np.arange(0.0, horizon + 1e-9, 2.0).round(8)).union(set(times.round(8)))), dtype=float)
        core = base.simulate(batch, theta, sim_times)
        if core is None:
            continue
        pred = pilot.co2_rate_from_core(theta, batch, core, times)
        obs = group["co2_raw"].to_numpy(dtype=float)
        mask = np.isfinite(pred) & np.isfinite(obs) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        pred = pred[mask]
        obs = obs[mask]
        times_used = times[mask]
        denom = float(np.dot(pred, pred))
        scale = float(np.dot(obs, pred) / denom) if denom > 1e-12 else 0.0
        pred_scaled = scale * pred
        err = pred_scaled - obs
        corr = np.nan
        if np.nanstd(obs) > 1e-12 and np.nanstd(pred_scaled) > 1e-12:
            corr = float(np.corrcoef(obs, pred_scaled)[0, 1])
        rows.append(
            {
                "group": "co2",
                "batch": str(batch_name),
                "state": "CO2_rate_shape",
                "n": int(mask.sum()),
                "scale_factor": scale,
                "rmse": float(np.sqrt(np.mean(err * err))),
                "mae": float(np.mean(np.abs(err))),
                "bias": float(np.mean(err)),
                "obs_scale": state_scale(obs, 1e-6),
                "median_obs": float(np.nanmedian(obs)),
                "median_pred": float(np.nanmedian(pred_scaled)),
                "relative_rmse": float(np.sqrt(np.mean(err * err)) / state_scale(obs, 1e-6)),
                "relative_bias": float(np.mean(err) / state_scale(obs, 1e-6)),
                "corr": corr,
            }
        )
        for t, y, yhat in zip(times_used, obs, pred_scaled):
            pred_rows.append(
                {
                    "batch": str(batch_name),
                    "time_h": float(t),
                    "state": "CO2_rate_shape",
                    "obs": float(y),
                    "pred": float(yhat),
                    "scale_factor": scale,
                }
            )
    metrics = pd.DataFrame(rows)
    if not metrics.empty:
        flagged = []
        reasons = []
        for _, row in metrics.iterrows():
            corr = float(row.get("corr", np.nan))
            rel_rmse = float(row.get("relative_rmse", np.nan))
            if np.isfinite(corr) and corr < 0.50:
                flagged.append("structure_warning")
                reasons.append("CO2 signal shape is poorly correlated with ethanol-rate proxy")
            elif rel_rmse > 0.75:
                flagged.append("weak_fit")
                reasons.append("CO2 shape error is high after per-batch scale fit")
            else:
                flagged.append("adequate")
                reasons.append("CO2 shape is compatible after per-batch scale fit")
        metrics["adequacy_flag"] = flagged
        metrics["adequacy_reason"] = reasons
    return metrics, pd.DataFrame(pred_rows)


def aroma_prediction_metrics(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sec_cache = aroma_sel.secondary_cache_for(batches, theta, core_cache)
    pred = aroma_sel.aroma_prediction_rows(theta, batches, variant, core_cache, sec_cache)
    if pred.empty:
        return pd.DataFrame(), pred
    metrics = summarize_prediction_rows(pred, ["species", "pool"], "aroma", scale_floor=1e-3)
    metrics["state"] = metrics["species"].astype(str) + ":" + metrics["pool"].astype(str)
    return add_adequacy_flags(metrics), pred


def lower_upper_log(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    lb = []
    ub = []
    for name in parameters:
        low, high = v2.V2_BOUNDS[name]
        lb.append(math.log(low))
        ub.append(math.log(high))
    return np.asarray(lb, dtype=float), np.asarray(ub, dtype=float)


def secondary_theta_from_log(x: np.ndarray, base_theta: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base_theta)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        theta[name] = float(math.exp(float(value)))
    return v2.clip_v2(theta)


def active_bound_count(theta: dict[str, float], parameters: tuple[str, ...]) -> int:
    count = 0
    for name in parameters:
        lb, ub = v2.V2_BOUNDS[name]
        value = float(theta[name])
        if value <= lb * 1.01 or value >= ub / 1.01:
            count += 1
    return count


def fit_secondary_variant(
    variant: SecondaryVariant,
    theta_reference: dict[str, float],
    batches: list[base.BatchData],
    core_cache: dict[str, pd.DataFrame],
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    base_theta = dict(theta_reference)
    base_theta.update(v2.default_theta_v2(base_theta))
    base_theta.update(variant.fixed_overrides)
    base_theta = v2.clip_v2(base_theta)
    parameters = variant.parameters
    lb, ub = lower_upper_log(parameters)
    x_ref = np.clip(v2.log_vector(base_theta, parameters), lb + 1e-9, ub - 1e-9)
    rng = np.random.default_rng(seed)

    prior_sigma = {
        "kPyrS_stat": 2.2,
        "kAldPyr": 2.0,
        "kAldS_stat": 2.2,
        "kAldO2": 2.0,
        "kAcAssim": 1.6,
        "kAldRed": 1.8,
        "kAcAld": 1.8,
        "kAcStress": 2.0,
    }
    rows = []
    best_theta = dict(base_theta)
    best_data_wsse = np.inf
    best_fit_score = np.inf
    best_pred = pd.DataFrame()

    def residual_with_priors(x: np.ndarray) -> np.ndarray:
        theta = secondary_theta_from_log(x, base_theta, parameters)
        res = [v2.residual_v2(theta, batches, core_cache)]
        priors = []
        for name in parameters:
            if name not in prior_sigma:
                continue
            ref = max(float(base_theta[name]), 1e-16)
            priors.append(math.log(max(float(theta[name]), 1e-16) / ref) / prior_sigma[name])
        if priors:
            res.append(np.asarray(priors, dtype=float))
        return np.concatenate(res)

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = x_ref.copy()
        else:
            x0 = x_ref + rng.normal(0.0, 1.0, size=len(parameters))
            free_new = [i for i, name in enumerate(parameters) if name in SECONDARY_OPTIONAL_TERMS]
            for pos in free_new:
                if rng.random() < 0.50:
                    x0[pos] = rng.uniform(lb[pos], ub[pos])
            x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)
        start_res = residual_with_priors(x0)
        result = least_squares(
            residual_with_priors,
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
        theta_hat = secondary_theta_from_log(result.x, base_theta, parameters)
        end_res = residual_with_priors(result.x)
        data_res = v2.residual_v2(theta_hat, batches, core_cache)
        data_wsse = float(np.dot(data_res, data_res))
        pred_rows = secondary_prediction_rows(theta_hat, batches, core_cache)
        rows.append(
            {
                "model": variant.name,
                "start": idx,
                "success": bool(result.success),
                "status": int(result.status),
                "message": str(result.message),
                "nfev": int(result.nfev),
                "initial_objective_wsse": float(np.dot(start_res, start_res)),
                "final_objective_wsse": float(np.dot(end_res, end_res)),
                "data_wsse": data_wsse,
                "n_data_residuals": int(len(data_res)),
                "n_parameters": int(len(parameters)),
                "wsse_per_data_residual": data_wsse / max(int(len(data_res)), 1),
                "active_bound_count": active_bound_count(theta_hat, parameters),
            }
        )
        fit_score = data_wsse + 25.0 * float(active_bound_count(theta_hat, parameters))
        rows[-1]["fit_selection_score"] = fit_score
        if fit_score < best_fit_score:
            best_fit_score = fit_score
            best_data_wsse = data_wsse
            best_theta = dict(theta_hat)
            best_pred = pred_rows
    return best_theta, pd.DataFrame(rows), best_pred


def information_criteria(data_wsse: float, n: int, k: int) -> tuple[float, float]:
    n = max(int(n), 1)
    k = max(int(k), 1)
    aic = float(data_wsse + 2 * k)
    if n > k + 1:
        aic += float((2 * k * (k + 1)) / (n - k - 1))
    bic = float(data_wsse + k * math.log(n))
    return aic, bic


def secondary_variant_penalty(metrics: pd.DataFrame, active_bounds: int) -> float:
    if metrics.empty:
        return 1e6
    penalty = 25.0 * float(active_bounds)
    for _, row in metrics.iterrows():
        rel_rmse = float(row.get("relative_rmse", np.nan))
        rel_bias = abs(float(row.get("relative_bias", np.nan)))
        if not np.isfinite(rel_rmse):
            continue
        penalty += 60.0 * max(rel_rmse - 0.85, 0.0)
        if np.isfinite(rel_bias):
            penalty += 70.0 * max(rel_bias - 0.40, 0.0)
    return float(penalty)


def plot_secondary_model_selection(model_summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if model_summary.empty:
        return
    ordered = model_summary.sort_values("selection_score")
    x = np.arange(len(ordered))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].bar(x, ordered["bic"], color="tab:blue")
    axes[0].set_title("BIC")
    axes[1].bar(x, ordered["state_penalty"], color="tab:orange")
    axes[1].set_title("Adequacy penalty")
    axes[2].bar(x, ordered["selection_score"], color="tab:green")
    axes[2].set_title("Selection score")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(ordered["model"], rotation=55, ha="right", fontsize=8)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "secondary_model_selection_scores.png", dpi=180)
    plt.close(fig)


def plot_secondary_selected_fit(pred: pd.DataFrame, selected_model: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred.empty:
        return
    for batch, group in pred.groupby("batch", sort=True):
        states = [state for state in ("Pyr", "AcAld", "Acetate", "O2") if group["state"].eq(state).any()]
        if not states:
            continue
        fig, axes = plt.subplots(len(states), 1, figsize=(10, 2.5 * len(states)), sharex=True)
        axes = np.atleast_1d(axes)
        for ax, state in zip(axes, states):
            sub = group[group["state"].eq(state)].sort_values("time_h")
            ax.scatter(sub["time_h"], sub["obs"], color="black", s=18, label="obs")
            ax.plot(sub["time_h"], sub["pred"], color="tab:blue", lw=1.5, label="pred")
            ax.set_ylabel(state)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("time (h)")
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right")
        fig.suptitle(f"{selected_model}: secondary-state fit for pilot {batch}", y=0.995)
        fig.tight_layout()
        fig.savefig(output_dir / f"selected_secondary_fit_{batch}.png", dpi=170)
        plt.close(fig)


def plot_adequacy_table(table: pd.DataFrame, output_dir: Path, label: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if table.empty:
        return
    plot_table = table.copy()
    plot_table["state_label"] = np.where(
        plot_table["group"].eq("aroma"),
        plot_table.get("state", ""),
        plot_table.get("state", ""),
    )
    plot_table["state_label"] = plot_table["group"].astype(str) + ":" + plot_table["state_label"].astype(str)
    plot_table = plot_table.sort_values("relative_rmse", ascending=False).head(25)
    fig, ax = plt.subplots(figsize=(11, max(4.0, 0.28 * len(plot_table))))
    colors = plot_table["adequacy_flag"].map(
        {
            "adequate": "tab:green",
            "limited_data": "0.6",
            "weak_fit": "tab:orange",
            "bias_warning": "tab:red",
            "structure_warning": "tab:red",
        }
    ).fillna("tab:blue")
    ax.barh(plot_table["state_label"], plot_table["relative_rmse"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("relative RMSE")
    ax.set_title(f"Model adequacy diagnostic: {label}")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / f"adequacy_{label}.png", dpi=180)
    plt.close(fig)


def plot_design_inputs_short(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], output_dir: Path) -> None:
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
        samples = pilot.future_sample_times(design)
        axes[2].vlines(samples, 0.0, 1.0, color="tab:blue", lw=1)
        axes[2].set_ylabel("Samples")
        axes[2].set_xlabel("time (h)")
        axes[2].grid(True, alpha=0.25)
        fig.suptitle(f"{int(row['campaign_order'])}. {design.name}", fontsize=10, y=0.99)
        fig.tight_layout()
        fig.savefig(output_dir / f"design_{int(row['campaign_order']):02d}.png", dpi=170)
        plt.close(fig)


def assemble_adequacy_table(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    aroma_variant: aroma_sel.AromaVariant,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    core_rows = core_prediction_rows(theta, primary_batches)
    core_metrics = add_adequacy_flags(summarize_prediction_rows(core_rows, ["state"], "core"))

    core_cache = joint.precompute_core_cache(extended_batches, theta)
    secondary_rows = secondary_prediction_rows(theta, extended_batches, core_cache)
    secondary_metrics = add_adequacy_flags(summarize_prediction_rows(secondary_rows, ["state"], "secondary"))

    aroma_metrics, aroma_rows = aroma_prediction_metrics(theta, extended_batches, aroma_variant, core_cache)
    co2_metrics, co2_rows = co2_shape_metrics(theta, {batch.batch: batch for batch in extended_batches}, co2_down)

    tables = [tbl for tbl in (core_metrics, secondary_metrics, aroma_metrics, co2_metrics) if not tbl.empty]
    adequacy = pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()
    return adequacy, {
        "core_predictions": core_rows,
        "secondary_predictions": secondary_rows,
        "aroma_predictions": aroma_rows,
        "co2_predictions": co2_rows,
    }


def write_report(
    initial_adequacy: pd.DataFrame,
    final_adequacy: pd.DataFrame,
    secondary_summary: pd.DataFrame,
    selected_secondary: SecondaryVariant,
    selected_aroma_name: str,
    target_parameters: tuple[str, ...],
    estim_current: pd.DataFrame,
    weak_current: pd.DataFrame,
    ranking: pd.DataFrame,
    selected_campaign: pd.DataFrame,
) -> None:
    path = RESULTS_DIR / "pilot_2025_global_model_selection_doe_report.md"
    with path.open("w", encoding="utf-8") as f:
        f.write("# Pilot 2025 global model-selection and DOE report\n\n")
        f.write("## Scope\n\n")
        f.write(
            "This run extends the previous ethyl-acetate model-selection workflow to the whole pilot "
            "fermentation model. The goal is to distinguish poor practical estimability from states whose "
            "current model structure does not adequately explain the observations.\n\n"
        )
        coverage_path = RESULTS_DIR / "observation_coverage.csv"
        if coverage_path.exists():
            coverage = pd.read_csv(coverage_path)
            f.write("## Observation coverage\n\n")
            f.write(coverage.to_markdown(index=False))
            f.write("\n\n")
        f.write("## Initial adequacy flags\n\n")
        if initial_adequacy.empty:
            f.write("_No adequacy metrics available._\n\n")
        else:
            cols = ["group", "state", "species", "pool", "n", "relative_rmse", "relative_bias", "corr", "adequacy_flag", "adequacy_reason"]
            cols = [col for col in cols if col in initial_adequacy.columns]
            f.write(initial_adequacy.sort_values(["adequacy_flag", "relative_rmse"], ascending=[False, False])[cols].to_markdown(index=False))
            f.write("\n\n")
        f.write("## Secondary-state model selection\n\n")
        f.write(f"Selected secondary structure: `{selected_secondary.name}`.\n\n")
        f.write(selected_secondary.description + "\n\n")
        f.write(selected_secondary.equation_markdown + "\n\n")
        if not secondary_summary.empty:
            f.write(secondary_summary.to_markdown(index=False))
            f.write("\n\n")
        f.write("## Final adequacy flags\n\n")
        if final_adequacy.empty:
            f.write("_No final adequacy metrics available._\n\n")
        else:
            cols = ["group", "state", "species", "pool", "n", "relative_rmse", "relative_bias", "corr", "adequacy_flag", "adequacy_reason"]
            cols = [col for col in cols if col in final_adequacy.columns]
            f.write(final_adequacy.sort_values(["adequacy_flag", "relative_rmse"], ascending=[False, False])[cols].to_markdown(index=False))
            f.write("\n\n")
        f.write("## Final target parameter set\n\n")
        f.write(f"Aroma structure inherited from the previous run: `{selected_aroma_name}`.\n\n")
        f.write(", ".join(f"`{p}`" for p in target_parameters))
        f.write("\n\n")
        f.write("## Current-data estimability\n\n")
        if not estim_current.empty:
            f.write(estim_current.to_markdown(index=False))
            f.write("\n\n")
        f.write("## Weak FIM directions\n\n")
        if not weak_current.empty:
            f.write(weak_current.to_markdown(index=False))
            f.write("\n\n")
        f.write("## DOE ranking\n\n")
        if ranking.empty:
            f.write("_DOE ranking skipped._\n\n")
        else:
            cols = ["candidate", "family", "combined_logdet", "combined_min_relative_eigenvalue", "target_worst_var_reduction", "hybrid_score", "rationale"]
            cols = [col for col in cols if col in ranking.columns]
            f.write(ranking.head(12)[cols].to_markdown(index=False))
            f.write("\n\n")
        f.write("## Selected campaign\n\n")
        if selected_campaign.empty:
            f.write("_No campaign selected._\n")
        else:
            f.write(selected_campaign.to_markdown(index=False))
            f.write("\n")


def create_notebook(
    selected_secondary: SecondaryVariant,
    selected_aroma_name: str,
) -> None:
    nb = nbformat.v4.new_notebook()
    cells = []
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""# Pilot 2025 global model selection and DOE

This notebook documents a full deterministic workflow for the pilot-scale fermentation data:

1. Load the curated pilot data and the previous ethyl-acetate aroma model.
2. Diagnose state-level model adequacy before re-fitting any new structure.
3. Compare secondary-metabolite ODE variants for pyruvate, acetaldehyde, and acetate.
4. Select a final ODE structure using weighted residual metrics, information criteria, and boundary diagnostics.
5. Build the current-data Fisher Information Matrix (FIM) with Pyomo DoE-compatible scaling.
6. Rank natural-must candidate experiments using model-based DOE metrics.

The distinction used here is:

$$r(\theta)=\frac{\hat{y}(\theta)-y}{\sigma_y}$$

$$J=\frac{\partial r}{\partial \log\theta}, \qquad F=J^\top J.$$

Poor estimability means that $F$ has weak directions after a model can already reproduce the data. Poor structural adequacy means that the residuals remain biased or shape-incompatible before the FIM interpretation is meaningful.
"""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import json
import pandas as pd
import numpy as np

RESULTS = Path('results/global_state_model_selection_doe')
print(RESULTS.resolve())
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Data and curation

The pilot workbook is treated as natural-must data. Online CO2 files are curated as in the previous pilot workflow: unusable CO2 files for batches 25150 and 25151 are excluded, and batch 25171 is shifted so that the sustained CO2 activation point is the effective fermentation start.

Dissolved oxygen is not used as a pilot observation here because the current pilot workbook does not expose a usable DO column; O2 remains a latent state in the secondary model.
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'co2_curation_decisions.csv')"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Observation coverage

The adequacy diagnostic is only interpreted for states with actual observations. States without observations can still affect the FIM through model coupling and future DOE outputs, but their current-data residual fit is not directly testable.
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'observation_coverage.csv')"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## State-level adequacy diagnostic

For each observed state, species, and pool, the diagnostic computes:

$$RMSE=\sqrt{\frac{1}{n}\sum_i(\hat{y}_i-y_i)^2}$$

$$rRMSE=\frac{RMSE}{\max(\operatorname{median}(|y|),0.25(\max y-\min y),\epsilon)}$$

$$rBias=\frac{\operatorname{mean}(\hat{y}-y)}{\max(\operatorname{median}(|y|),0.25(\max y-\min y),\epsilon)}.$$

The flag is structural when normalized error is high and/or residuals show systematic shape incompatibility. Limited data are not interpreted as structural failure.
"""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """initial = pd.read_csv(RESULTS / 'adequacy_initial_current.csv')
cols = [c for c in ['group','state','species','pool','n','relative_rmse','relative_bias','corr','adequacy_flag','adequacy_reason'] if c in initial.columns]
initial.sort_values(['adequacy_flag','relative_rmse'], ascending=[False, False])[cols]
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Secondary-state structural variants

The selected aroma model from the previous workflow is kept fixed while secondary-metabolite structure is tested. This avoids confusing ethyl-acetate structural error with pyruvate, acetaldehyde, or acetate structural error.

The common gating terms are:

$$\phi_N=\frac{N}{N+K_N}, \qquad \phi_{stat}=1-\phi_N,$$

$$g_{O2}=\frac{O_2}{O_2+K_{O2}}, \qquad \phi_{ana}=\frac{K_{O2}}{O_2+K_{O2}}, \qquad \phi_E=\frac{E}{E+K_E}.$$

Candidate structures are compared using robust log-parameter least squares. The model-selection score is:

$$Score = BIC + P_{bounds}+P_{adequacy}.$$

This penalizes models that fit by driving parameters to bounds or by keeping large state-level residual bias.
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'secondary_model_selection_summary.csv')"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            f"""## Selected secondary structure

Selected structure: `{selected_secondary.name}`

{selected_secondary.description}

{selected_secondary.equation_markdown}
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'theta_selected_global_model.csv', index_col=0).head(60)"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Final adequacy after structural selection

The table below repeats the adequacy diagnostic after substituting the selected secondary structure into the global model. States that remain flagged are not automatically discarded; instead, they indicate where either additional experimental excitation, additional measurements, or stronger literature priors are needed before those parameters should be used as flexible MPCC degrees of freedom.
"""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """final = pd.read_csv(RESULTS / 'adequacy_final_selected.csv')
cols = [c for c in ['group','state','species','pool','n','relative_rmse','relative_bias','corr','adequacy_flag','adequacy_reason'] if c in final.columns]
final.sort_values(['adequacy_flag','relative_rmse'], ascending=[False, False])[cols]
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## FIM and estimability

The FIM is built from current-data residuals using finite differences in log-parameter space:

$$J_{ij}\approx\frac{r_i(\theta_j e^h)-r_i(\theta_j e^{-h})}{2h}.$$

The eigenspectrum of $F$ diagnoses practical identifiability. The weakest eigenvectors identify confounded parameter combinations. The approximate log-standard deviation is read from a stabilized inverse FIM:

$$\Sigma_{\log\theta}\approx F^{-1}.$$

Parameters are classified as well estimated, moderate, weak but actionable, or weak/confounded using the same thresholds as the previous pilot aroma workflow.
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("json.load(open(RESULTS / 'target_parameters_global_selected.json'))"))
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'parameter_estimability_global_current.csv')"))
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'weak_directions_global_current.csv')"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Model-based DOE

Candidate natural-must designs are ranked by adding each candidate FIM to the current-data FIM. The hybrid score keeps D-optimality as the main objective while penalizing weak minimum eigen-directions:

$$\Phi_{hybrid}=\log\det(F)-2|\log(\lambda_{min}/\lambda_{max})|-0.05\log(\operatorname{tr}(F^{-1})).$$

The greedy campaign then adds the experiment that maximizes this score at each step. This is homologous to the previous `doe_multiexperiment` logic: current-data FIM acts as the prior information matrix and candidate experiments add information sequentially.
"""
        )
    )
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'candidate_ranking_global_selected.csv').head(15)"))
    cells.append(nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'selected_campaign_hybrid_global_selected.csv')"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            f"""## Interpretation

The final model combines:

- The primary fermentation and glycerol model used in the previous pilot calibration.
- The selected secondary structure `{selected_secondary.name}`.
- The selected aroma structure `{selected_aroma_name}` with Antoine plus UNIFAC water-ethanol partition and CO2-driven stripping.

The practical decision rule is:

- If a state is structurally adequate and a parameter is weak in the FIM, prioritize DOE excitation.
- If a state remains structurally flagged, do not interpret weak FIM directions as only a sampling/design problem; first improve model structure, measurement interpretation, or priors.
- If a parameter is active at a bound and dominates weak directions, fix or regularize it before using the model as an MPCC constraint.
"""
        )
    )
    nb["cells"] = cells
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pilot 2025 global state model selection and DOE.")
    parser.add_argument("--n-starts-secondary", type=int, default=6)
    parser.add_argument("--max-nfev-secondary", type=int, default=160)
    parser.add_argument("--sensitivity-step", type=float, default=0.015)
    parser.add_argument("--campaign-size", type=int, default=6)
    parser.add_argument("--skip-doe", action="store_true")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("[load] pilot data", flush=True)
    model_data, curated_co2, co2_decisions = pilot.load_pilot_model_data()
    model_data.to_csv(RESULTS_DIR / "pilot_model_data_used.csv", index=False)
    co2_decisions.to_csv(RESULTS_DIR / "co2_curation_decisions.csv", index=False)
    co2_down = pilot.downsample_co2(curated_co2)
    co2_down.to_csv(RESULTS_DIR / "co2_curated_downsampled.csv", index=False)
    primary_batches = pilot.make_primary_batches(model_data)
    extended_batches = pilot.make_extended_batches(model_data)
    coverage_rows = []
    for state in base.STATE_NAMES:
        n_obs = int(sum(pd.Series(batch.observations.get(state, [])).notna().sum() for batch in primary_batches))
        coverage_rows.append({"group": "core", "state": state, "n_observations": n_obs})
    for state in ("Pyr", "AcAld", "Acetate", "O2"):
        n_obs = int(sum(pd.Series(batch.observations.get(state, [])).notna().sum() for batch in extended_batches))
        coverage_rows.append({"group": "secondary", "state": state, "n_observations": n_obs})
    for species in joint.AROMA_SPECIES:
        total_col = joint.AROMA_COLUMNS[species]
        cond_col = pilot.AROMA_CONDENSATE_COLUMNS[species]
        n_total = int(sum(pd.Series(batch.observations.get(total_col, [])).notna().sum() for batch in extended_batches))
        n_cond = int(sum(pd.Series(batch.observations.get(cond_col, [])).gt(aroma_sel.CONDENSATE_POSITIVE_THRESHOLD).sum() for batch in extended_batches))
        coverage_rows.append({"group": "aroma_total", "state": species, "n_observations": n_total})
        coverage_rows.append({"group": "aroma_condensate", "state": species, "n_observations": n_cond})
    pd.DataFrame(coverage_rows).to_csv(RESULTS_DIR / "observation_coverage.csv", index=False)

    theta0, selected_aroma_name, selected_aroma = load_selected_aroma_context()
    pd.Series(theta0).to_csv(RESULTS_DIR / "theta_start_from_selected_aroma.csv")
    (RESULTS_DIR / "selected_aroma_model_inherited.txt").write_text(selected_aroma_name, encoding="utf-8")

    print("[diagnostic] initial state adequacy", flush=True)
    initial_adequacy, initial_predictions = assemble_adequacy_table(theta0, primary_batches, extended_batches, co2_down, selected_aroma)
    initial_adequacy.to_csv(RESULTS_DIR / "adequacy_initial_current.csv", index=False)
    for name, frame in initial_predictions.items():
        frame.to_csv(RESULTS_DIR / f"{name}_initial_current.csv", index=False)
    flagged = initial_adequacy[~initial_adequacy["adequacy_flag"].isin(["adequate", "limited_data"])].copy()
    flagged.to_csv(RESULTS_DIR / "flagged_states_initial_current.csv", index=False)
    plot_adequacy_table(initial_adequacy, plot_dir, "initial_current")

    print("[secondary] model selection", flush=True)
    variants = secondary_variant_library()
    core_cache = joint.precompute_core_cache(extended_batches, theta0)
    all_fit_rows = []
    all_metric_rows = []
    model_rows = []
    theta_by_model: dict[str, dict[str, float]] = {}
    pred_by_model: dict[str, pd.DataFrame] = {}
    for idx, variant in enumerate(variants.values(), start=1):
        print(f"[secondary] {idx}/{len(variants)} {variant.name}", flush=True)
        theta_hat, fit_rows, pred_rows = fit_secondary_variant(
            variant,
            theta0,
            extended_batches,
            core_cache,
            n_starts=args.n_starts_secondary,
            max_nfev=args.max_nfev_secondary,
            seed=20260628 + idx,
        )
        theta_model = dict(theta0)
        theta_model.update(theta_hat)
        theta_by_model[variant.name] = theta_model
        pred_by_model[variant.name] = pred_rows
        pd.Series({p: theta_model[p] for p in variant.parameters}).to_csv(RESULTS_DIR / f"theta_secondary_{variant.name}.csv")
        fit_rows.to_csv(RESULTS_DIR / f"secondary_fit_multistart_{variant.name}.csv", index=False)
        pred_rows.to_csv(RESULTS_DIR / f"secondary_predictions_{variant.name}.csv", index=False)
        all_fit_rows.append(fit_rows)
        metrics = add_adequacy_flags(summarize_prediction_rows(pred_rows, ["state"], "secondary"))
        metrics["model"] = variant.name
        metrics.to_csv(RESULTS_DIR / f"secondary_state_metrics_{variant.name}.csv", index=False)
        all_metric_rows.append(metrics)
        best_row = fit_rows.sort_values("fit_selection_score").iloc[0]
        aic, bic = information_criteria(best_row["data_wsse"], best_row["n_data_residuals"], len(variant.parameters))
        state_penalty = secondary_variant_penalty(metrics, int(best_row["active_bound_count"]))
        model_rows.append(
            {
                "model": variant.name,
                "data_wsse": float(best_row["data_wsse"]),
                "n_data_residuals": int(best_row["n_data_residuals"]),
                "n_parameters": int(len(variant.parameters)),
                "active_bound_count": int(best_row["active_bound_count"]),
                "aicc": aic,
                "bic": bic,
                "state_penalty": state_penalty,
                "selection_score": bic + state_penalty,
                "description": variant.description,
            }
        )

    pd.concat(all_fit_rows, ignore_index=True, sort=False).to_csv(RESULTS_DIR / "secondary_variant_fit_multistart_summary.csv", index=False)
    pd.concat(all_metric_rows, ignore_index=True, sort=False).to_csv(RESULTS_DIR / "secondary_variant_state_metrics.csv", index=False)
    secondary_summary = pd.DataFrame(model_rows).sort_values(["selection_score", "bic", "data_wsse"]).reset_index(drop=True)
    selected_secondary_name = str(secondary_summary.iloc[0]["model"])
    selected_secondary = variants[selected_secondary_name]
    secondary_summary["selected"] = secondary_summary["model"].eq(selected_secondary_name)
    secondary_summary.to_csv(RESULTS_DIR / "secondary_model_selection_summary.csv", index=False)
    (RESULTS_DIR / "selected_secondary_model.txt").write_text(selected_secondary_name, encoding="utf-8")
    plot_secondary_model_selection(secondary_summary, plot_dir)

    theta_selected = dict(theta0)
    theta_selected.update(theta_by_model[selected_secondary_name])
    pd.Series(theta_selected).to_csv(RESULTS_DIR / "theta_selected_global_model.csv")
    pd.Series({p: theta_selected[p] for p in selected_secondary.parameters}).to_csv(RESULTS_DIR / "theta_selected_secondary_model.csv")
    plot_secondary_selected_fit(pred_by_model[selected_secondary_name], selected_secondary_name, plot_dir / "selected_secondary_fit")

    print("[diagnostic] final state adequacy", flush=True)
    final_adequacy, final_predictions = assemble_adequacy_table(theta_selected, primary_batches, extended_batches, co2_down, selected_aroma)
    final_adequacy.to_csv(RESULTS_DIR / "adequacy_final_selected.csv", index=False)
    for name, frame in final_predictions.items():
        frame.to_csv(RESULTS_DIR / f"{name}_final_selected.csv", index=False)
    final_flagged = final_adequacy[~final_adequacy["adequacy_flag"].isin(["adequate", "limited_data"])].copy()
    final_flagged.to_csv(RESULTS_DIR / "flagged_states_final_selected.csv", index=False)
    plot_adequacy_table(final_adequacy, plot_dir, "final_selected")

    target_parameters = tuple(dict.fromkeys(pilot.CORE_TARGETS + selected_secondary.parameters + selected_aroma.parameters))
    (RESULTS_DIR / "target_parameters_global_selected.json").write_text(json.dumps(list(target_parameters), indent=2), encoding="utf-8")

    print("[fim] current-data FIM", flush=True)
    jac, resid = aroma_sel.finite_difference_jacobian_generic(
        theta_selected,
        target_parameters,
        lambda th: aroma_sel.combined_residual_selected(th, primary_batches, extended_batches, co2_down, selected_aroma),
        args.sensitivity_step,
    )
    fim = 0.5 * ((jac.T @ jac) + (jac.T @ jac).T)
    pd.DataFrame(jac, columns=target_parameters).to_csv(RESULTS_DIR / "jacobian_global_current.csv", index=False)
    pd.Series(resid).to_csv(RESULTS_DIR / "residual_global_current.csv", index=False)
    pd.DataFrame(fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_global_current.csv")
    if rescale_FIM is not None:
        try:
            scaled = rescale_FIM(fim, np.asarray([theta_selected[p] for p in target_parameters], dtype=float))
            pd.DataFrame(scaled, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_global_current_pyomo_rescaled.csv")
        except Exception as err:
            (RESULTS_DIR / "fim_global_current_pyomo_rescaled_error.txt").write_text(str(err), encoding="utf-8")
    spectrum, weak, estimability = aroma_sel.fim_diagnostics_generic(fim, theta_selected, target_parameters, "global_current")
    spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_global_current.csv", index=False)
    weak.to_csv(RESULTS_DIR / "weak_directions_global_current.csv", index=False)
    estimability.to_csv(RESULTS_DIR / "parameter_estimability_global_current.csv", index=False)

    ranking = pd.DataFrame()
    selected_campaign = pd.DataFrame()
    if not args.skip_doe:
        print("[doe] candidate FIMs", flush=True)
        designs = aroma_sel.add_extra_designs(model_data, pilot.natural_candidate_designs(model_data))
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
        pd.DataFrame(design_rows).to_csv(RESULTS_DIR / "candidate_design_library_global_selected.csv", index=False)
        candidate_fims = {}
        for idx, (name, design) in enumerate(designs.items(), start=1):
            print(f"[doe] candidate {idx}/{len(designs)} {name}", flush=True)
            cfim = aroma_sel.candidate_fim_selected(theta_selected, design, selected_aroma, target_parameters, args.sensitivity_step)
            candidate_fims[name] = cfim
            pd.DataFrame(cfim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / f"candidate_fim_global_{name}.csv")
        ranking = aroma_sel.rank_candidates(candidate_fims, designs, fim, target_parameters)
        ranking.to_csv(RESULTS_DIR / "candidate_ranking_global_selected.csv", index=False)
        selected_campaign, campaign_fim = aroma_sel.greedy_select(candidate_fims, designs, fim, target_parameters, args.campaign_size, objective="hybrid")
        selected_campaign.to_csv(RESULTS_DIR / "selected_campaign_hybrid_global_selected.csv", index=False)
        pd.DataFrame(campaign_fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_global_current_plus_campaign_hybrid.csv")
        spec_after, weak_after, estim_after = aroma_sel.fim_diagnostics_generic(
            campaign_fim,
            theta_selected,
            target_parameters,
            "global_current_plus_campaign",
        )
        spec_after.to_csv(RESULTS_DIR / "eigen_spectrum_global_current_plus_campaign.csv", index=False)
        weak_after.to_csv(RESULTS_DIR / "weak_directions_global_current_plus_campaign.csv", index=False)
        estim_after.to_csv(RESULTS_DIR / "parameter_estimability_global_current_plus_campaign.csv", index=False)
        plot_design_inputs_short(selected_campaign, designs, plot_dir / "designs")

    write_report(
        initial_adequacy,
        final_adequacy,
        secondary_summary,
        selected_secondary,
        selected_aroma_name,
        target_parameters,
        estimability,
        weak,
        ranking,
        selected_campaign,
    )
    create_notebook(selected_secondary, selected_aroma_name)
    print(f"[done] selected secondary model: {selected_secondary_name}", flush=True)
    print(f"[done] inherited aroma model: {selected_aroma_name}", flush=True)
    print(f"[done] results: {RESULTS_DIR}", flush=True)
    print(f"[done] notebook: {NOTEBOOK_PATH}", flush=True)


if __name__ == "__main__":
    main()
