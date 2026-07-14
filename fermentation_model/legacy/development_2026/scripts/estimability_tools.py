from __future__ import annotations

import math
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, qmc


def summarize_residual_contributions(state_fit_comparison: pd.DataFrame) -> pd.DataFrame:
    """Return residual contribution percentages by fit and state."""
    if state_fit_comparison is None or state_fit_comparison.empty:
        return pd.DataFrame()

    table = state_fit_comparison.copy()
    pieces = []
    for objective_col in ["squared_error", "weighted_squared_error"]:
        if objective_col not in table.columns:
            continue
        contribution = table[objective_col] / table.groupby(level=0)[objective_col].transform("sum")
        pieces.append((100.0 * contribution).rename(f"{objective_col}_pct"))
    if not pieces:
        return pd.DataFrame()
    return pd.concat([table] + pieces, axis=1)


def latin_hypercube_samples(
    parameter_bounds: Mapping[str, tuple[float, float]],
    parameters: Sequence[str],
    n_samples: int,
    seed: int = 532,
    log_scale: bool = True,
    linear_parameters: Sequence[str] = (),
) -> pd.DataFrame:
    """Generate Latin-hypercube samples inside parameter bounds."""
    parameters = list(parameters)
    linear_parameters = set(linear_parameters)
    sampler = qmc.LatinHypercube(d=len(parameters), seed=seed)
    unit_samples = sampler.random(n=int(n_samples))

    rows = []
    for sample in unit_samples:
        row = {}
        for u, name in zip(sample, parameters):
            lb, ub = parameter_bounds[name]
            lb = float(lb)
            ub = float(ub)
            if log_scale and name not in linear_parameters and lb > 0.0 and ub > lb:
                row[name] = 10.0 ** (math.log10(lb) + float(u) * (math.log10(ub) - math.log10(lb)))
            else:
                row[name] = lb + float(u) * (ub - lb)
        rows.append(row)
    return pd.DataFrame(rows, columns=parameters)


def summarize_multistart_results(results: pd.DataFrame, objective_col: str = "final_objective") -> pd.DataFrame:
    """Compact summary by objective label for successful multistart runs."""
    if results is None or results.empty:
        return pd.DataFrame()
    ok = results[results["success"]].copy()
    if ok.empty:
        return pd.DataFrame()
    rows = []
    for label, group in ok.groupby("fit_label"):
        rounded_obj = group[objective_col].round(6)
        best = group.loc[group[objective_col].idxmin()]
        rows.append(
            {
                "fit_label": label,
                "n_success": len(group),
                "n_unique_objectives_rounded": rounded_obj.nunique(),
                "best_run": best["run_id"],
                "best_objective": best[objective_col],
                "median_objective": group[objective_col].median(),
                "worst_success_objective": group[objective_col].max(),
            }
        )
    return pd.DataFrame(rows).set_index("fit_label")


def plot_multistart_results(
    results: pd.DataFrame,
    x_parameter: str,
    y_parameter: str,
    objective_col: str = "final_objective",
):
    """Plot LHS initial points colored by final objective value."""
    if results is None or results.empty:
        return None, None
    x_col = f"initial_{x_parameter}"
    y_col = f"initial_{y_parameter}"
    if x_col not in results.columns or y_col not in results.columns:
        raise KeyError(f"Missing {x_col} or {y_col} in multistart results.")

    labels = list(results["fit_label"].dropna().unique())
    fig, axes = plt.subplots(1, len(labels), figsize=(6 * len(labels), 4.5), squeeze=False)
    axes = axes.ravel()
    for ax, label in zip(axes, labels):
        subset = results[results["fit_label"].eq(label)]
        ok = subset[subset["success"]]
        failed = subset[~subset["success"]]
        if not ok.empty:
            points = ax.scatter(ok[x_col], ok[y_col], c=ok[objective_col], cmap="viridis", s=60, alpha=0.85)
            fig.colorbar(points, ax=ax, label=objective_col)
        if not failed.empty:
            ax.scatter(failed[x_col], failed[y_col], marker="x", color="tab:red", s=60, label="failed")
            ax.legend(loc="best")
        ax.set_title(label)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, axes


def covariance_estimability_summary(covariance_results: Mapping[tuple[str, str], pd.DataFrame]):
    """Compute eigen diagnostics and weak covariance directions from covariance matrices."""
    eigen_rows = []
    loading_rows = []
    for (fit_label, method), cov in covariance_results.items():
        if cov is None or cov.empty:
            continue
        cov = cov.astype(float)
        cov_array = cov.to_numpy(dtype=float)
        cov_symmetric = 0.5 * (cov_array + cov_array.T)
        eigvals, eigvecs = np.linalg.eigh(cov_symmetric)
        order = np.argsort(eigvals)[::-1]
        eigvals = eigvals[order]
        eigvecs = eigvecs[:, order]
        positive = eigvals[eigvals > 0.0]
        eigen_rows.append(
            {
                "fit": fit_label,
                "method": method,
                "n_parameters": len(cov),
                "max_cov_eigenvalue": float(eigvals[0]),
                "min_cov_eigenvalue": float(eigvals[-1]),
                "cov_condition_number": float(positive.max() / positive.min()) if len(positive) else np.nan,
                "is_psd": bool(eigvals[-1] >= -1e-10),
            }
        )
        for direction_idx in range(min(3, eigvecs.shape[1])):
            vector = eigvecs[:, direction_idx]
            abs_vector = np.abs(vector)
            top_order = np.argsort(abs_vector)[::-1][: min(5, len(abs_vector))]
            loading_rows.append(
                {
                    "fit": fit_label,
                    "method": method,
                    "weak_direction": direction_idx + 1,
                    "cov_eigenvalue": float(eigvals[direction_idx]),
                    "dominant_parameters": ", ".join(cov.index[i] for i in top_order),
                    "dominant_abs_loadings": ", ".join(f"{abs_vector[i]:.3f}" for i in top_order),
                }
            )
    eigen_summary = pd.DataFrame(eigen_rows)
    if not eigen_summary.empty:
        eigen_summary = eigen_summary.set_index(["fit", "method"]).sort_index()
    weak_loadings = pd.DataFrame(loading_rows)
    if not weak_loadings.empty:
        weak_loadings = weak_loadings.set_index(["fit", "method", "weak_direction"]).sort_index()
    return eigen_summary, weak_loadings


def build_profile_grid(
    theta_hat: Mapping[str, float],
    parameter_bounds: Mapping[str, tuple[float, float]],
    parameter: str,
    std_dev: float | None = None,
    n_grid: int = 7,
    std_span: float = 3.0,
    relative_span: float = 0.5,
    active_rtol: float = 1e-4,
) -> np.ndarray:
    """Create a bounded profile grid around the current estimate."""
    value = float(theta_hat[parameter])
    lb, ub = parameter_bounds[parameter]
    lb = float(lb)
    ub = float(ub)
    scale = max(abs(value), abs(lb), abs(ub), 1.0)
    near_lb = abs(value - lb) <= active_rtol * scale
    near_ub = abs(value - ub) <= active_rtol * scale

    if std_dev is not None and np.isfinite(std_dev) and std_dev > 0.0:
        span = max(std_span * float(std_dev), relative_span * abs(value))
    else:
        span = max(relative_span * abs(value), 0.05 * (ub - lb), 1e-8)

    if near_lb:
        low = lb
        high = min(ub, max(value + span, value * 10.0 if value > 0 else lb + span))
    elif near_ub:
        low = max(lb, min(value - span, value * 0.5 if value > 0 else ub - span))
        high = ub
    else:
        low = max(lb, value - span)
        high = min(ub, value + span)

    if low == high:
        low = max(lb, value * 0.9)
        high = min(ub, value * 1.1 + 1e-8)
    grid = np.linspace(low, high, int(n_grid))
    if not np.any(np.isclose(grid, value, rtol=1e-10, atol=1e-12)):
        grid = np.unique(np.concatenate([grid, np.array([value])]))
    return np.sort(grid)


def summarize_profile_likelihood(profile_df: pd.DataFrame, alpha: float = 0.95) -> pd.DataFrame:
    """Summarize whether each profile crosses the chi-square threshold."""
    if profile_df is None or profile_df.empty:
        return pd.DataFrame()
    threshold = float(chi2.ppf(alpha, df=1))
    ok = profile_df[profile_df["success"]].copy()
    if ok.empty:
        return pd.DataFrame()
    rows = []
    for parameter, group in ok.groupby("profiled_theta"):
        theta_hat = float(group["theta_hat"].iloc[0])
        left = group[group["theta_value"] < theta_hat]
        right = group[group["theta_value"] > theta_hat]
        rows.append(
            {
                "profiled_theta": parameter,
                "n_success": len(group),
                "theta_hat": theta_hat,
                "min_profile_objective": group["objective"].min(),
                "max_lr_stat": group["lr_stat"].max(),
                "crosses_left": bool((left["lr_stat"] >= threshold).any()),
                "crosses_right": bool((right["lr_stat"] >= threshold).any()),
                "chi2_threshold": threshold,
            }
        )
    return pd.DataFrame(rows).set_index("profiled_theta").sort_index()


def plot_profile_likelihood_profiles(profile_df: pd.DataFrame, alpha: float = 0.95):
    """Plot profile likelihood curves from the custom profile table."""
    ok = profile_df[profile_df["success"]].copy() if profile_df is not None and not profile_df.empty else pd.DataFrame()
    if ok.empty:
        return None, None
    threshold = float(chi2.ppf(alpha, df=1))
    parameters = list(ok["profiled_theta"].drop_duplicates())
    ncols = 2
    nrows = int(math.ceil(len(parameters) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 3.8 * nrows), squeeze=False)
    axes_flat = axes.ravel()
    for ax, parameter in zip(axes_flat, parameters):
        group = ok[ok["profiled_theta"].eq(parameter)].sort_values("theta_value")
        ax.plot(group["theta_value"], group["lr_stat"], marker="o")
        ax.axhline(threshold, color="tab:red", linestyle="--", label=f"{alpha:.0%} threshold")
        ax.axvline(float(group["theta_hat"].iloc[0]), color="black", linestyle=":", label="theta_hat")
        ax.set_title(parameter)
        ax.set_xlabel("fixed parameter value")
        ax.set_ylabel("2 * objective increase")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
    for ax in axes_flat[len(parameters) :]:
        ax.axis("off")
    fig.tight_layout()
    return fig, axes
