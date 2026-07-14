from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

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
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

import run_new_must_glycerol_estimability_doe as base
import run_pilot_2025_aroma_model_selection_doe as aroma_sel
import run_pilot_2025_calibration_estimability as pilot
import run_pilot_2025_co2_stripping_benchmark as co2bench
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "integrated_co2_aroma_doe"
NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_integrated_co2_aroma_doe.ipynb"
EXECUTED_NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_integrated_co2_aroma_doe.executed.ipynb"
PLOT_DIR = RESULTS_DIR / "plots"

CO2_MODE_CANDIDATES = ("instant", "threshold", "lag_threshold")
CO2_PARAMETER_NAMES = {
    "tau_h": "qgas_tau_h",
    "threshold_fraction": "qgas_threshold_fraction",
    "power_p": "qgas_power_p",
}
INTEGRATED_BOUNDS = {
    "qgas_tau_h": (0.25, 96.0),
    "qgas_threshold_fraction": (0.0, 0.90),
    "qgas_power_p": (0.35, 2.50),
}
INTEGRATED_DEFAULTS = {
    "qgas_tau_h": 36.0,
    "qgas_threshold_fraction": 0.50,
    "qgas_power_p": 1.0,
}


def mode_parameter_names(mode: co2bench.CO2Mode) -> tuple[str, ...]:
    return tuple(CO2_PARAMETER_NAMES[name] for name in mode.parameters)


def parameter_bounds(name: str) -> tuple[float, float]:
    if name in INTEGRATED_BOUNDS:
        return INTEGRATED_BOUNDS[name]
    return aroma_sel.PARAMETER_BOUNDS[name]


def clip_theta(theta: dict[str, float]) -> dict[str, float]:
    out = aroma_sel.clip_theta(theta)
    for name, (lb, ub) in INTEGRATED_BOUNDS.items():
        if name in out and np.isfinite(float(out[name])):
            out[name] = float(np.clip(float(out[name]), lb, ub))
    return out


def load_start_theta() -> tuple[dict[str, float], str, aroma_sel.AromaVariant]:
    paths = [
        SCRIPT_DIR / "results" / "global_sugar" / "theta_selected_global_model.csv",
        SCRIPT_DIR / "results" / "global_state_model_selection_doe" / "theta_selected_global_model.csv",
        aroma_sel.RESULTS_DIR / "theta_selected_aroma_model.csv",
    ]
    theta = pilot.load_initial_theta()
    theta.update(v2.default_theta_v2(theta))
    for path in paths:
        if path.exists():
            loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
            theta.update({str(k): float(v) for k, v in loaded.items() if np.isfinite(float(v))})
            break
    selected_name = "ea_ethanol_nlimited"
    selected_path = aroma_sel.RESULTS_DIR / "selected_model.txt"
    if selected_path.exists():
        selected_name = selected_path.read_text(encoding="utf-8").strip()
    variants = aroma_sel.variant_library()
    selected_variant = variants.get(selected_name, variants["ea_ethanol_nlimited"])
    for name, value in aroma_sel.EXTRA_DEFAULTS.items():
        theta.setdefault(name, value)
    for name, value in INTEGRATED_DEFAULTS.items():
        theta.setdefault(name, value)
    benchmark_path = SCRIPT_DIR / "results" / "co2_stripping_benchmark" / "co2_model_selection_summary.csv"
    if benchmark_path.exists():
        summary = pd.read_csv(benchmark_path)
        for _, row in summary.iterrows():
            mode_name = str(row["mode"])
            if mode_name not in CO2_MODE_CANDIDATES:
                continue
            try:
                params = eval(str(row["params"]), {"__builtins__": {}})
            except Exception:
                continue
            for source, value in dict(params).items():
                target = CO2_PARAMETER_NAMES.get(source)
                if target:
                    theta[target] = float(value)
            break
    return clip_theta(theta), selected_name, selected_variant


def co2_params_from_theta(theta: dict[str, float], mode: co2bench.CO2Mode) -> dict[str, float]:
    params = dict(mode.defaults)
    for source in mode.parameters:
        target = CO2_PARAMETER_NAMES[source]
        params[source] = float(theta.get(target, INTEGRATED_DEFAULTS[target]))
    return params


def loglike_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    values = []
    for name in parameters:
        value = float(theta[name])
        if name == "qgas_threshold_fraction":
            values.append(value)
        else:
            values.append(math.log(max(value, 1e-16)))
    return np.asarray(values, dtype=float)


def theta_from_vector(x: np.ndarray, base_theta: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base_theta)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        if name == "qgas_threshold_fraction":
            theta[name] = float(value)
        else:
            theta[name] = float(math.exp(float(value)))
    return clip_theta(theta)


def vector_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    lb = []
    ub = []
    for name in parameters:
        low, high = parameter_bounds(name)
        if name == "qgas_threshold_fraction":
            lb.append(low)
            ub.append(high)
        else:
            lb.append(math.log(low))
            ub.append(math.log(high))
    return np.asarray(lb, dtype=float), np.asarray(ub, dtype=float)


def integrated_current_residual(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    mode: co2bench.CO2Mode,
    variant: aroma_sel.AromaVariant,
    include_primary_secondary: bool = True,
) -> np.ndarray:
    residuals: list[np.ndarray] = []
    core_cache = joint.precompute_core_cache(extended_batches, theta)
    sec_cache = aroma_sel.secondary_cache_for(extended_batches, theta, core_cache)
    params = co2_params_from_theta(theta, mode)
    if include_primary_secondary:
        residuals.append(base.residual_vector(theta, primary_batches))
        residuals.append(v2.residual_v2(theta, extended_batches, core_cache))
    residuals.append(co2bench.co2_shape_residual(theta, {batch.batch: batch for batch in extended_batches}, co2_down, mode, params))
    qgas_cache = co2bench.gas_flow_cache_for_batches(theta, extended_batches, core_cache, mode, params)
    residuals.append(
        co2bench.aroma_residual_combo(
            theta,
            extended_batches,
            variant,
            core_cache,
            sec_cache,
            qgas_cache,
            "water_ethanol_total_sugar_as_glucose",
        )
    )
    return np.concatenate([r for r in residuals if len(r)]) if residuals else np.array([], dtype=float)


def fit_integrated_mode(
    theta0: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_down: pd.DataFrame,
    mode: co2bench.CO2Mode,
    variant: aroma_sel.AromaVariant,
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fit_parameters = tuple(dict.fromkeys(variant.parameters + mode_parameter_names(mode)))
    lb, ub = vector_bounds(fit_parameters)
    x_ref = np.clip(loglike_vector(theta0, fit_parameters), lb + 1e-9, ub - 1e-9)
    rng = np.random.default_rng(seed)
    rows = []
    best_score = np.inf
    best_theta = dict(theta0)

    def data_residual(theta: dict[str, float]) -> np.ndarray:
        return integrated_current_residual(theta, primary_batches, extended_batches, co2_down, mode, variant, include_primary_secondary=False)

    def objective(x: np.ndarray) -> np.ndarray:
        theta = theta_from_vector(x, theta0, fit_parameters)
        res = [data_residual(theta)]
        priors = []
        for name in variant.parameters:
            scale = 1.15 if name.startswith("alpha_") else 1.8
            if name.startswith("k_EA_") or name == "q10_EA":
                scale = 4.0
            ref = max(float(theta0[name]), 1e-16)
            priors.append(math.log(max(float(theta[name]), 1e-16) / ref) / scale)
        if priors:
            res.append(np.asarray(priors, dtype=float))
        return np.concatenate([r for r in res if len(r)])

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = x_ref.copy()
        else:
            x0 = x_ref.copy()
            for pos, name in enumerate(fit_parameters):
                low, high = lb[pos], ub[pos]
                if name == "qgas_threshold_fraction":
                    x0[pos] = rng.uniform(low, high)
                elif name.startswith("qgas_"):
                    x0[pos] = rng.uniform(low, high)
                else:
                    x0[pos] = np.clip(x0[pos] + rng.normal(0.0, 0.75), low + 1e-9, high - 1e-9)
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
        theta_hat = theta_from_vector(result.x, theta0, fit_parameters)
        end = objective(result.x)
        data = data_residual(theta_hat)
        active = active_bound_count(theta_hat, fit_parameters)
        score = float(np.dot(data, data)) + 20.0 * active
        rows.append(
            {
                "mode": mode.name,
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
        )
        for name in fit_parameters:
            rows[-1][name] = float(theta_hat[name])
        if score < best_score:
            best_score = score
            best_theta = theta_hat

    core_cache = joint.precompute_core_cache(extended_batches, best_theta)
    sec_cache = aroma_sel.secondary_cache_for(extended_batches, best_theta, core_cache)
    params = co2_params_from_theta(best_theta, mode)
    qgas_cache = co2bench.gas_flow_cache_for_batches(best_theta, extended_batches, core_cache, mode, params)
    aroma_pred, aroma_metrics, aroma_wsse, aroma_n = co2bench.evaluate_aroma_combo(
        best_theta,
        extended_batches,
        variant,
        core_cache,
        sec_cache,
        qgas_cache,
        mode.name,
        "water_ethanol_total_sugar_as_glucose",
    )
    co2_metrics, co2_pred = co2bench.co2_metrics_by_batch(best_theta, {batch.batch: batch for batch in extended_batches}, co2_down, mode, params)
    co2_pred.to_csv(RESULTS_DIR / f"co2_predictions_{mode.name}.csv", index=False)
    aroma_pred.to_csv(RESULTS_DIR / f"aroma_predictions_{mode.name}.csv", index=False)
    if not aroma_metrics.empty:
        aroma_metrics["mode"] = mode.name
        aroma_metrics["aroma_wsse"] = aroma_wsse
        aroma_metrics["aroma_n_residuals"] = aroma_n
    if not co2_metrics.empty:
        co2_metrics["mode"] = mode.name
    return best_theta, pd.DataFrame(rows), aroma_metrics, co2_metrics


def load_completed_mode(mode_name: str) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame, pd.DataFrame] | None:
    theta_path = RESULTS_DIR / f"theta_integrated_{mode_name}.csv"
    fit_path = RESULTS_DIR / f"fit_multistart_integrated_{mode_name}.csv"
    aroma_path = RESULTS_DIR / f"aroma_pool_metrics_{mode_name}.csv"
    co2_path = RESULTS_DIR / f"co2_metrics_{mode_name}.csv"
    if not theta_path.exists() or not fit_path.exists():
        return None
    theta = {str(k): float(v) for k, v in pd.read_csv(theta_path, index_col=0).iloc[:, 0].to_dict().items()}
    fit_rows = pd.read_csv(fit_path)
    aroma_metrics = pd.read_csv(aroma_path) if aroma_path.exists() else pd.DataFrame()
    co2_metrics = pd.read_csv(co2_path) if co2_path.exists() else pd.DataFrame()
    if fit_rows.empty:
        return None
    return theta, fit_rows, aroma_metrics, co2_metrics


def active_bound_count(theta: dict[str, float], parameters: tuple[str, ...]) -> int:
    count = 0
    for name in parameters:
        lb, ub = parameter_bounds(name)
        value = float(theta[name])
        if value <= lb * 1.01 or value >= ub / 1.01:
            count += 1
    return count


def information_criteria(data_wsse: float, n: int, k: int) -> tuple[float, float]:
    n = max(int(n), 1)
    k = max(int(k), 1)
    aic = float(data_wsse + 2 * k)
    if n > k + 1:
        aic += float((2 * k * (k + 1)) / (n - k - 1))
    bic = float(data_wsse + k * math.log(n))
    return aic, bic


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
        lb, ub = parameter_bounds(name)
        if name == "qgas_threshold_fraction":
            delta = step * max(ub - lb, 1e-9)
            theta_plus[name] = float(np.clip(theta[name] + delta, lb, ub))
            theta_minus[name] = float(np.clip(theta[name] - delta, lb, ub))
            denom = max(theta_plus[name] - theta_minus[name], 1e-12)
        else:
            theta_plus[name] = float(np.clip(float(theta[name]) * math.exp(step), lb, ub))
            theta_minus[name] = float(np.clip(float(theta[name]) * math.exp(-step), lb, ub))
            denom = 2.0 * step
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / denom)
    return np.column_stack(cols), base_res


def fim_metrics(fim: np.ndarray, prefix: str = "") -> dict[str, float]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eig = np.linalg.eigvalsh(fim)
    max_eig = float(np.max(eig)) if eig.size else np.nan
    floor = max(max_eig * 1e-12, np.finfo(float).tiny) if np.isfinite(max_eig) and max_eig > 0 else np.finfo(float).tiny
    eig_pos = np.clip(eig, floor, None)
    cov = pilot.stable_inverse(fim)
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
    cov = pilot.stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std = float(math.sqrt(max(cov[idx, idx], 0.0)))
        lb, ub = parameter_bounds(name)
        value = float(theta[name])
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
                "std_approx": std,
                "approx_95_multiplier": float(math.exp(1.96 * min(std, 20.0))),
                "fim_diag": float(fim[idx, idx]),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return spectrum, pd.DataFrame(weak_rows), pd.DataFrame(rows)


def future_residual(
    theta: dict[str, float],
    design: base.FutureDesign,
    mode: co2bench.CO2Mode,
    variant: aroma_sel.AromaVariant,
) -> np.ndarray:
    sample_times = pilot.future_sample_times(design)
    co2_times = pilot.future_co2_times(design)
    all_times = np.asarray(sorted(set(sample_times).union(set(co2_times)).union({0.0, float(design.horizon_h)})), dtype=float)
    core = base.simulate(design, theta, all_times)
    if core is None:
        return np.ones(1000, dtype=float) * 1e6
    secondary = v2.integrate_secondary_v2(design, theta, core)
    params = co2_params_from_theta(theta, mode)
    qprod = co2bench.base_co2_production(theta, design, core, core.index.to_numpy(dtype=float))
    qgas = pd.Series(
        co2bench.gas_flow_from_production(core.index.to_numpy(dtype=float), qprod, mode, params),
        index=core.index,
    )
    try:
        liq, cond = co2bench.integrate_aroma_combo(
            design,
            theta,
            core,
            secondary,
            qgas,
            variant,
            "water_ethanol_total_sugar_as_glucose",
        )
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
        aroma_sample = liq.loc[valid_samples]
        cond_sample = cond.loc[valid_samples]
        for species in joint.AROMA_SPECIES:
            liq_center = aroma_sample[species].to_numpy(dtype=float)
            cond_center = cond_sample[species].to_numpy(dtype=float)
            residuals.append(liq_center / pilot.aroma_sigma(species, "liquid", liq_center))
            residuals.append(cond_center / pilot.aroma_sigma(species, "condensate", cond_center))
    valid_co2 = [float(t) for t in co2_times if float(t) in core.index]
    if valid_co2:
        qgas_vals = np.asarray([co2bench.interp_series(qgas, t) for t in valid_co2], dtype=float)
        sigma = np.maximum(0.05, 0.20 * np.maximum(qgas_vals, 0.05))
        residuals.append(qgas_vals / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def candidate_fim(theta: dict[str, float], design: base.FutureDesign, mode: co2bench.CO2Mode, variant: aroma_sel.AromaVariant, parameters: tuple[str, ...], step: float) -> np.ndarray:
    jac, _ = finite_difference_jacobian(theta, parameters, lambda th: future_residual(th, design, mode, variant), step)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = pilot.stable_inverse(prior)
    post_cov = pilot.stable_inverse(combined)
    rows: dict[str, float] = {}
    reductions = []
    aroma_reductions = []
    for idx, name in enumerate(parameters):
        before = float(prior_cov[idx, idx])
        after = float(post_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        reduction = 1.0 - ratio if np.isfinite(ratio) else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = reduction
        if np.isfinite(reduction):
            reductions.append(reduction)
            if name.startswith("k_") or name.startswith("alpha_") or name.startswith("qgas_"):
                aroma_reductions.append(reduction)
    rows["target_mean_var_reduction"] = float(np.nanmean(reductions)) if reductions else np.nan
    rows["target_worst_var_reduction"] = float(np.nanmin(reductions)) if reductions else np.nan
    rows["aroma_co2_mean_var_reduction"] = float(np.nanmean(aroma_reductions)) if aroma_reductions else np.nan
    rows["aroma_co2_worst_var_reduction"] = float(np.nanmin(aroma_reductions)) if aroma_reductions else np.nan
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
        samples = pilot.future_sample_times(design)
        axes[2].vlines(samples, 0.0, 1.0, color="tab:blue", lw=1)
        axes[2].set_ylabel("Samples")
        axes[2].set_xlabel("time (h)")
        axes[2].grid(True, alpha=0.25)
        fig.suptitle(f"{int(row['campaign_order'])}. {design.name}", fontsize=10, y=0.99)
        fig.tight_layout()
        fig.savefig(output_dir / f"design_{int(row['campaign_order']):02d}.png", dpi=170)
        plt.close(fig)


def write_report(model_summary: pd.DataFrame, estimability: pd.DataFrame, weak: pd.DataFrame, ranking: pd.DataFrame, selected: pd.DataFrame, selected_mode: str, target_parameters: tuple[str, ...]) -> None:
    path = RESULTS_DIR / "integrated_co2_aroma_doe_report.md"
    with path.open("w", encoding="utf-8") as f:
        f.write("# Integrated CO2-aroma DOE report\n\n")
        f.write("Partition mode: `water_ethanol_total_sugar_as_glucose`.\n\n")
        f.write("CO2 gas-flow candidates: `instant`, `threshold`, and `lag_threshold`.\n\n")
        f.write("## Model selection\n\n")
        f.write(model_summary.to_markdown(index=False))
        f.write(f"\n\nSelected mode: `{selected_mode}`.\n\n")
        f.write("## Target parameters\n\n")
        f.write(", ".join(f"`{p}`" for p in target_parameters))
        f.write("\n\n## Current-data estimability\n\n")
        f.write(estimability.to_markdown(index=False))
        f.write("\n\n## Weak directions\n\n")
        f.write(weak.to_markdown(index=False))
        f.write("\n\n## DOE ranking\n\n")
        if not ranking.empty:
            cols = ["candidate", "family", "combined_logdet", "combined_min_relative_eigenvalue", "target_worst_var_reduction", "aroma_co2_worst_var_reduction", "hybrid_score", "rationale"]
            cols = [c for c in cols if c in ranking.columns]
            f.write(ranking.head(12)[cols].to_markdown(index=False))
        f.write("\n\n## Selected campaign\n\n")
        f.write(selected.to_markdown(index=False) if not selected.empty else "_No campaign selected._")
        f.write("\n")


def create_notebook(selected_mode: str) -> None:
    nb = nbformat.v4.new_notebook()
    cells = [
        nbformat.v4.new_markdown_cell(
            r"""# Integrated CO2-aroma model selection and DOE

This notebook evaluates a formal model candidate where the gas-liquid partition model includes sugar:

$$K_i(T,E,G,F)=\frac{\gamma_i^{UNIFAC}(T,x_{water},x_{ethanol},x_{G+F})P_i^{sat,Antoine}(T)}{RTC_{tot,L}}.$$

The aroma loss is driven by an effective gas stripping flow:

$$r_{loss,i}=\alpha_iK_iq_{gas}C_{L,i}.$$

Three alternatives are compared:

- `instant`: \(q_{gas}=q_{CO2,prod}\)
- `threshold`: \(q_{gas}=\max(q_{CO2,prod}-q_0,0)\)
- `lag_threshold`: \(\tau dq_{gas}/dt=\max(q_{CO2,prod}-q_0,0)-q_{gas}\)

The selected model is then used to compute current-data FIM diagnostics and a natural-must DOE campaign.
"""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import json
import pandas as pd
RESULTS = Path('results/integrated_co2_aroma_doe')
print(RESULTS.resolve())
"""
        ),
        nbformat.v4.new_markdown_cell("## Integrated Model Selection"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'integrated_model_selection_summary.csv')"),
        nbformat.v4.new_markdown_cell(f"Selected CO2 mode: `{selected_mode}`."),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'theta_selected_integrated_model.csv', index_col=0).head(80)"),
        nbformat.v4.new_markdown_cell("## Fit Metrics"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'integrated_aroma_pool_metrics.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'integrated_co2_metrics.csv')"),
        nbformat.v4.new_markdown_cell(
            r"""## FIM Diagnostics

The FIM uses finite differences in log-parameter space for positive parameters and direct finite differences for the bounded threshold fraction.
"""
        ),
        nbformat.v4.new_code_cell("json.load(open(RESULTS / 'target_parameters_integrated_selected.json'))"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'parameter_estimability_integrated_current.csv')"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'weak_directions_integrated_current.csv')"),
        nbformat.v4.new_markdown_cell("## DOE"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'candidate_ranking_integrated_selected.csv').head(15)"),
        nbformat.v4.new_code_cell("pd.read_csv(RESULTS / 'selected_campaign_hybrid_integrated_selected.csv')"),
    ]
    nb["cells"] = cells
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrated CO2 effective gas flow and sugar-aware aroma DOE.")
    parser.add_argument("--n-starts", type=int, default=8)
    parser.add_argument("--max-nfev", type=int, default=180)
    parser.add_argument("--sensitivity-step", type=float, default=0.015)
    parser.add_argument("--campaign-size", type=int, default=6)
    parser.add_argument("--skip-doe", action="store_true")
    parser.add_argument(
        "--modes",
        default=",".join(CO2_MODE_CANDIDATES),
        help="Comma-separated CO2 modes to fit/evaluate. Defaults to all candidates.",
    )
    parser.add_argument("--resume", action="store_true", help="Reuse completed per-mode fit files when available.")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] pilot data and selected global sugar theta", flush=True)
    model_data, curated_co2, co2_decisions = pilot.load_pilot_model_data()
    co2_down = pilot.downsample_co2(curated_co2)
    primary_batches = pilot.make_primary_batches(model_data)
    extended_batches = pilot.make_extended_batches(model_data)
    theta0, selected_aroma_name, selected_variant = load_start_theta()
    model_data.to_csv(RESULTS_DIR / "pilot_model_data_used.csv", index=False)
    co2_decisions.to_csv(RESULTS_DIR / "co2_curation_decisions.csv", index=False)
    co2_down.to_csv(RESULTS_DIR / "co2_downsampled.csv", index=False)
    pd.Series(theta0).to_csv(RESULTS_DIR / "theta_start_integrated.csv")
    (RESULTS_DIR / "selected_aroma_model_inherited.txt").write_text(selected_aroma_name, encoding="utf-8")

    modes = co2bench.co2_mode_library()
    all_fit = []
    all_aroma_metrics = []
    all_co2_metrics = []
    theta_by_mode: dict[str, dict[str, float]] = {}
    summary_rows = []
    requested_modes = tuple(name.strip() for name in str(args.modes).split(",") if name.strip())
    invalid_modes = sorted(set(requested_modes).difference(CO2_MODE_CANDIDATES))
    if invalid_modes:
        raise ValueError(f"Unknown CO2 modes: {invalid_modes}. Valid modes: {CO2_MODE_CANDIDATES}")
    if not requested_modes:
        raise ValueError("At least one CO2 mode must be requested.")

    for idx, mode_name in enumerate(requested_modes, start=1):
        mode = modes[mode_name]
        print(f"[fit] {idx}/{len(requested_modes)} {mode.name}", flush=True)
        completed = load_completed_mode(mode.name) if args.resume else None
        if completed is None:
            theta_hat, fit_rows, aroma_metrics, co2_metrics = fit_integrated_mode(
                theta0,
                primary_batches,
                extended_batches,
                co2_down,
                mode,
                selected_variant,
                n_starts=args.n_starts,
                max_nfev=args.max_nfev,
                seed=20260629 + idx,
            )
            fit_rows.to_csv(RESULTS_DIR / f"fit_multistart_integrated_{mode.name}.csv", index=False)
            pd.Series(theta_hat).to_csv(RESULTS_DIR / f"theta_integrated_{mode.name}.csv")
        else:
            print(f"[fit] reusing completed {mode.name}", flush=True)
            theta_hat, fit_rows, aroma_metrics, co2_metrics = completed
        if not aroma_metrics.empty:
            aroma_metrics.to_csv(RESULTS_DIR / f"aroma_pool_metrics_{mode.name}.csv", index=False)
            all_aroma_metrics.append(aroma_metrics)
        if not co2_metrics.empty:
            co2_metrics.to_csv(RESULTS_DIR / f"co2_metrics_{mode.name}.csv", index=False)
            all_co2_metrics.append(co2_metrics)
        all_fit.append(fit_rows)
        theta_by_mode[mode.name] = theta_hat
        best = fit_rows.sort_values("fit_selection_score").iloc[0]
        aic, bic = information_criteria(float(best["data_wsse"]), int(best["n_data_residuals"]), int(best["n_parameters"]))
        summary_rows.append(
            {
                "mode": mode.name,
                "data_wsse": float(best["data_wsse"]),
                "n_data_residuals": int(best["n_data_residuals"]),
                "n_parameters": int(best["n_parameters"]),
                "active_bound_count": int(best["active_bound_count"]),
                "aicc": aic,
                "bic": bic,
                "selection_score": bic + 20.0 * int(best["active_bound_count"]),
                "description": mode.description,
            }
        )

    fit_summary = pd.concat(all_fit, ignore_index=True, sort=False)
    fit_summary.to_csv(RESULTS_DIR / "integrated_fit_multistart_summary.csv", index=False)
    aroma_metrics_all = pd.concat(all_aroma_metrics, ignore_index=True, sort=False) if all_aroma_metrics else pd.DataFrame()
    co2_metrics_all = pd.concat(all_co2_metrics, ignore_index=True, sort=False) if all_co2_metrics else pd.DataFrame()
    aroma_metrics_all.to_csv(RESULTS_DIR / "integrated_aroma_pool_metrics.csv", index=False)
    co2_metrics_all.to_csv(RESULTS_DIR / "integrated_co2_metrics.csv", index=False)
    model_summary = pd.DataFrame(summary_rows).sort_values(["selection_score", "bic", "data_wsse"]).reset_index(drop=True)
    selected_mode_name = str(model_summary.iloc[0]["mode"])
    model_summary["selected"] = model_summary["mode"].eq(selected_mode_name)
    model_summary.to_csv(RESULTS_DIR / "integrated_model_selection_summary.csv", index=False)
    (RESULTS_DIR / "selected_integrated_mode.txt").write_text(selected_mode_name, encoding="utf-8")
    theta_selected = theta_by_mode[selected_mode_name]
    pd.Series(theta_selected).to_csv(RESULTS_DIR / "theta_selected_integrated_model.csv")

    selected_mode = modes[selected_mode_name]
    target_parameters = tuple(
        dict.fromkeys(pilot.CORE_TARGETS + tuple(v2.V2_CHEM_PARAMETERS) + selected_variant.parameters + mode_parameter_names(selected_mode))
    )
    (RESULTS_DIR / "target_parameters_integrated_selected.json").write_text(json.dumps(list(target_parameters), indent=2), encoding="utf-8")

    print("[fim] current integrated model", flush=True)
    jac, resid = finite_difference_jacobian(
        theta_selected,
        target_parameters,
        lambda th: integrated_current_residual(th, primary_batches, extended_batches, co2_down, selected_mode, selected_variant, include_primary_secondary=True),
        args.sensitivity_step,
    )
    fim = 0.5 * ((jac.T @ jac) + (jac.T @ jac).T)
    pd.DataFrame(jac, columns=target_parameters).to_csv(RESULTS_DIR / "jacobian_integrated_current.csv", index=False)
    pd.Series(resid).to_csv(RESULTS_DIR / "residual_integrated_current.csv", index=False)
    pd.DataFrame(fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_integrated_current.csv")
    if rescale_FIM is not None:
        try:
            scaled = rescale_FIM(fim, np.asarray([theta_selected[p] for p in target_parameters], dtype=float))
            pd.DataFrame(scaled, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_integrated_current_pyomo_rescaled.csv")
        except Exception as err:
            (RESULTS_DIR / "fim_integrated_current_pyomo_rescaled_error.txt").write_text(str(err), encoding="utf-8")
    spectrum, weak, estimability = fim_diagnostics(fim, theta_selected, target_parameters, "integrated_current")
    spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_integrated_current.csv", index=False)
    weak.to_csv(RESULTS_DIR / "weak_directions_integrated_current.csv", index=False)
    estimability.to_csv(RESULTS_DIR / "parameter_estimability_integrated_current.csv", index=False)

    ranking = pd.DataFrame()
    selected_campaign = pd.DataFrame()
    if not args.skip_doe:
        print("[doe] integrated candidate FIMs", flush=True)
        designs = aroma_sel.add_extra_designs(model_data, pilot.natural_candidate_designs(model_data))
        pd.DataFrame(
            [
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
                for design in designs.values()
            ]
        ).to_csv(RESULTS_DIR / "candidate_design_library_integrated_selected.csv", index=False)
        candidate_fims = {}
        for idx, (name, design) in enumerate(designs.items(), start=1):
            print(f"[doe] candidate {idx}/{len(designs)} {name}", flush=True)
            cfim = candidate_fim(theta_selected, design, selected_mode, selected_variant, target_parameters, args.sensitivity_step)
            candidate_fims[name] = cfim
            pd.DataFrame(cfim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / f"candidate_fim_integrated_{name}.csv")
        ranking = rank_candidates(candidate_fims, designs, fim, target_parameters)
        ranking.to_csv(RESULTS_DIR / "candidate_ranking_integrated_selected.csv", index=False)
        selected_campaign, campaign_fim = greedy_select(candidate_fims, designs, fim, target_parameters, args.campaign_size, objective="hybrid")
        selected_campaign.to_csv(RESULTS_DIR / "selected_campaign_hybrid_integrated_selected.csv", index=False)
        pd.DataFrame(campaign_fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_integrated_current_plus_campaign_hybrid.csv")
        spec_after, weak_after, estim_after = fim_diagnostics(campaign_fim, theta_selected, target_parameters, "integrated_current_plus_campaign")
        spec_after.to_csv(RESULTS_DIR / "eigen_spectrum_integrated_current_plus_campaign.csv", index=False)
        weak_after.to_csv(RESULTS_DIR / "weak_directions_integrated_current_plus_campaign.csv", index=False)
        estim_after.to_csv(RESULTS_DIR / "parameter_estimability_integrated_current_plus_campaign.csv", index=False)
        plot_design_inputs(selected_campaign, designs, PLOT_DIR / "designs")

    write_report(model_summary, estimability, weak, ranking, selected_campaign, selected_mode_name, target_parameters)
    create_notebook(selected_mode_name)
    print(f"[done] selected integrated CO2 mode: {selected_mode_name}", flush=True)
    print(f"[done] results: {RESULTS_DIR}", flush=True)
    print(f"[done] notebook: {NOTEBOOK_PATH}", flush=True)


if __name__ == "__main__":
    main()
