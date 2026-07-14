from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_secondary_joint_campaign_doe as v1
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "secondary_fit_capacity"
NOTEBOOK_PATH = SCRIPT_DIR / "fermentation_secondary_fit_capacity.ipynb"

FIT_REDUCED = v2.V2_REDUCED_PARAMETERS
FIT_ALL = v2.V2_REDUCED_PARAMETERS + v2.V2_O2_PARAMETERS


def load_v2_reduced_theta(theta_v1: dict[str, float]) -> dict[str, float]:
    theta = {**theta_v1, **v2.default_theta_v2(theta_v1)}
    path = SCRIPT_DIR / "results" / "secondary_v2_model_evaluation" / "theta_secondary_v2_reduced_o2fixed.csv"
    if path.exists():
        loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
        theta.update({str(k): float(value) for k, value in loaded.items() if np.isfinite(float(value))})
    theta.update(
        {
            "kPyrS_stat": v2.V2_BOUNDS["kPyrS_stat"][0],
            "kAldPyr": v2.V2_BOUNDS["kAldPyr"][0],
            "kAldS_stat": v2.V2_BOUNDS["kAldS_stat"][0],
            "kAldO2": v2.V2_BOUNDS["kAldO2"][0],
            "kAcAssim": v2.V2_BOUNDS["kAcAssim"][0],
        }
    )
    return v2.clip_v2(theta)


def has_secondary_observations(batch) -> bool:
    for state in ("Pyr", "AcAld", "Acetate", "O2"):
        values = np.asarray(batch.observations.get(state, []), dtype=float)
        if np.isfinite(values).any():
            return True
    return False


def residual_for_batches(theta: dict[str, float], batches: list, core_cache: dict[str, pd.DataFrame]) -> np.ndarray:
    return v2.residual_v2(theta, batches, core_cache)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    low = []
    high = []
    for name in parameters:
        lb, ub = v2.V2_BOUNDS[name]
        low.append(math.log(lb))
        high.append(math.log(ub))
    return np.asarray(low, dtype=float), np.asarray(high, dtype=float)


def unpack_theta(x: np.ndarray, base_theta: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base_theta)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        theta[name] = float(math.exp(float(value)))
    return v2.clip_v2(theta)


def fit_capacity(
    label: str,
    base_theta: dict[str, float],
    batches: list,
    core_cache: dict[str, pd.DataFrame],
    parameters: tuple[str, ...],
    max_nfev: int,
    prior_sigma: dict[str, float] | None = None,
) -> tuple[dict[str, float], dict[str, float]]:
    prior_sigma = {} if prior_sigma is None else dict(prior_sigma)
    x0 = np.asarray([math.log(float(base_theta[name])) for name in parameters], dtype=float)
    lb, ub = log_bounds(parameters)
    x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)

    def fun(x: np.ndarray) -> np.ndarray:
        theta = unpack_theta(x, base_theta, parameters)
        pieces = [residual_for_batches(theta, batches, core_cache)]
        for name, sigma in prior_sigma.items():
            if name in parameters:
                pieces.append(np.array([math.log(theta[name] / base_theta[name]) / sigma], dtype=float))
        return np.concatenate(pieces)

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
    theta_hat = unpack_theta(result.x, base_theta, parameters)
    end = fun(result.x)
    summary = {
        "fit_label": label,
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nfev": int(result.nfev),
        "n_parameters_fit": int(len(parameters)),
        "n_batches": int(len(batches)),
        "n_residuals": int(len(end)),
        "initial_wsse": float(np.dot(start, start)),
        "final_wsse": float(np.dot(end, end)),
        "wsse_per_residual": float(np.dot(end, end) / max(len(end), 1)),
    }
    return theta_hat, summary


def state_metrics(label: str, theta_by_batch: dict[str, dict[str, float]], batches: list, core_cache: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for batch in batches:
        theta = theta_by_batch.get(batch.label)
        core = core_cache.get(batch.label)
        if theta is None or core is None:
            continue
        sec = v2.integrate_secondary_v2(batch, theta, core)
        for state in ("Pyr", "AcAld", "Acetate", "O2"):
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sec.loc[batch.time, state].to_numpy(dtype=float)
            err = pred[mask] - obs[mask]
            scaled = err / v2.SIGMA[state]
            rows.append(
                {
                    "fit_label": label,
                    "medium": batch.medium,
                    "batch": batch.batch,
                    "state": state,
                    "n_obs": int(len(err)),
                    "wsse": float(np.dot(scaled, scaled)),
                    "wsse_per_obs": float(np.dot(scaled, scaled) / max(len(err), 1)),
                    "rmse": float(np.sqrt(np.mean(err * err))),
                    "mae": float(np.mean(np.abs(err))),
                    "obs_median": float(np.nanmedian(obs[mask])),
                    "obs_range": float(np.nanmax(obs[mask]) - np.nanmin(obs[mask])),
                }
            )
    return pd.DataFrame(rows)


def aggregate_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (label, state), group in metrics.groupby(["fit_label", "state"], sort=True):
        total_wsse = float(group["wsse"].sum())
        total_n = int(group["n_obs"].sum())
        rows.append(
            {
                "fit_label": label,
                "state": state,
                "n_obs": total_n,
                "wsse": total_wsse,
                "wsse_per_obs": total_wsse / max(total_n, 1),
                "rmse_weighted_mean": float(np.average(group["rmse"], weights=group["n_obs"])),
                "mae_weighted_mean": float(np.average(group["mae"], weights=group["n_obs"])),
            }
        )
    out = pd.DataFrame(rows)
    total_rows = []
    for label, group in metrics.groupby("fit_label", sort=True):
        total_wsse = float(group["wsse"].sum())
        total_n = int(group["n_obs"].sum())
        total_rows.append(
            {
                "fit_label": label,
                "state": "ALL",
                "n_obs": total_n,
                "wsse": total_wsse,
                "wsse_per_obs": total_wsse / max(total_n, 1),
                "rmse_weighted_mean": np.nan,
                "mae_weighted_mean": np.nan,
            }
        )
    return pd.concat([out, pd.DataFrame(total_rows)], ignore_index=True)


def parameter_spread(theta_rows: pd.DataFrame, parameters: tuple[str, ...], label: str) -> pd.DataFrame:
    rows = []
    if theta_rows.empty:
        return pd.DataFrame()
    for name in parameters:
        values = pd.to_numeric(theta_rows[name], errors="coerce").dropna()
        if values.empty:
            continue
        median = float(values.median())
        rows.append(
            {
                "fit_label": label,
                "parameter": name,
                "n": int(len(values)),
                "min": float(values.min()),
                "median": median,
                "max": float(values.max()),
                "fold_range": float(values.max() / max(values.min(), 1e-30)),
                "robust_cv": float((values.quantile(0.75) - values.quantile(0.25)) / max(median, 1e-30)),
            }
        )
    return pd.DataFrame(rows)


def plot_capacity(
    batches: list,
    core_cache: dict[str, pd.DataFrame],
    global_theta: dict[str, float],
    batch_theta: dict[str, dict[str, float]],
) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        if not has_secondary_observations(batch):
            continue
        core = core_cache.get(batch.label)
        theta_batch = batch_theta.get(batch.label)
        if core is None or theta_batch is None:
            continue
        sec_global = v2.integrate_secondary_v2(batch, global_theta, core)
        sec_batch = v2.integrate_secondary_v2(batch, theta_batch, core)
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes, ("Pyr", "AcAld", "Acetate", "O2")):
            ax.plot(sec_global.index, sec_global[state], "--", lw=1.3, label="global v2 reduced")
            ax.plot(sec_batch.index, sec_batch[state], "-", lw=1.7, label="batch-fit flexible")
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if mask.any():
                ax.scatter(batch.time[mask], obs[mask], s=24, color="black", label="obs")
            ax.set_title(state)
            ax.set_xlabel("time [h]")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=8)
        fig.suptitle(f"Fit-capacity diagnostic: {batch.medium}/{batch.batch}")
        fig.tight_layout()
        fig.savefig(plot_dir / f"secondary_fit_capacity_{batch.medium}_{batch.batch}.png", dpi=170)
        plt.close(fig)


def write_report(
    fit_summary: pd.DataFrame,
    aggregate: pd.DataFrame,
    spread: pd.DataFrame,
) -> None:
    lines = [
        "# Secondary fit-capacity diagnostic",
        "",
        "## Purpose",
        "",
        "This diagnostic asks whether the secondary-state data can be fitted at all when parameter estimability and parameter transferability are deliberately relaxed.",
        "",
        "The comparison separates three explanations for poor global curves:",
        "",
        "- If per-batch fits are good but global fits are poor, the equation form has enough local flexibility but the parameters are not transferable across batches/media.",
        "- If per-medium fits improve but per-batch fits improve much more, medium-specific effects are present but do not explain all variability.",
        "- If even per-batch fits are poor, the secondary-state structure, core forcing, or observations are inconsistent with the ODE form.",
        "",
        "## Fit summary",
        "",
        fit_summary.to_markdown(index=False),
        "",
        "## Aggregated state metrics",
        "",
        aggregate.to_markdown(index=False),
        "",
        "## Parameter spread in per-batch all-free fits",
        "",
        spread.to_markdown(index=False),
        "",
        "## Interpretation guide",
        "",
        "Treat per-batch all-free fits as an upper bound on fit capacity, not as a calibratable model. Large parameter fold-ranges mean the data can be matched only by sacrificing parameter transferability and identifiability.",
    ]
    (RESULTS_DIR / "secondary_fit_capacity_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_notebook() -> None:
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Secondary fit-capacity diagnostic\n\n"
            "This notebook tests whether the secondary-state data can be fitted when parameter identifiability and transferability are relaxed."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "ROOT = Path.cwd()\n"
            "if not (ROOT / 'results').exists() and (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'results/secondary_fit_capacity'\n"
            "display(pd.read_csv(RESULTS / 'fit_capacity_summary.csv'))\n"
            "display(pd.read_csv(RESULTS / 'fit_capacity_aggregate_metrics.csv'))\n"
            "display(pd.read_csv(RESULTS / 'batch_all_o2free_parameter_spread.csv'))"
        ),
        nbformat.v4.new_markdown_cell("## Batch-Level Metrics"),
        nbformat.v4.new_code_cell(
            "metrics = pd.read_csv(RESULTS / 'fit_capacity_state_metrics.csv')\n"
            "display(metrics.sort_values(['fit_label','medium','batch','state']).head(80))"
        ),
        nbformat.v4.new_markdown_cell("## Fit-Capacity Plots"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_fit_capacity_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Report"),
        nbformat.v4.new_code_cell("print((RESULTS / 'secondary_fit_capacity_report.md').read_text(encoding='utf-8'))"),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> int:
    parser = argparse.ArgumentParser(description="Secondary-state fit-capacity diagnostic.")
    parser.add_argument("--max-nfev-medium", type=int, default=80)
    parser.add_argument("--max-nfev-batch", type=int, default=80)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)

    data = v1.load_current_review_data()
    batches = [batch for batch in v1.make_secondary_batches(data) if has_secondary_observations(batch)]
    theta_v1 = v2.load_theta_v1()
    theta_global = load_v2_reduced_theta(theta_v1)
    core_cache = v1.precompute_core_cache(batches, theta_v1)

    fit_rows = []
    metric_frames = []

    global_theta_by_batch = {batch.label: theta_global for batch in batches}
    global_res = residual_for_batches(theta_global, batches, core_cache)
    fit_rows.append(
        {
            "fit_label": "global_v2_reduced_o2fixed",
            "success": True,
            "status": np.nan,
            "message": "loaded from v2 reduced evaluation",
            "nfev": np.nan,
            "n_parameters_fit": int(len(FIT_REDUCED)),
            "n_batches": int(len(batches)),
            "n_residuals": int(len(global_res)),
            "initial_wsse": np.nan,
            "final_wsse": float(np.dot(global_res, global_res)),
            "wsse_per_residual": float(np.dot(global_res, global_res) / max(len(global_res), 1)),
        }
    )
    metric_frames.append(state_metrics("global_v2_reduced_o2fixed", global_theta_by_batch, batches, core_cache))

    medium_theta_by_batch: dict[str, dict[str, float]] = {}
    for medium, medium_batches in pd.Series(batches).groupby([batch.medium for batch in batches], sort=True):
        fitted, summary = fit_capacity(
            f"medium_v2_reduced_o2fixed_{medium}",
            theta_global,
            list(medium_batches),
            core_cache,
            FIT_REDUCED,
            args.max_nfev_medium,
            prior_sigma={"kAldRed": 2.0},
        )
        fit_rows.append(summary)
        for batch in medium_batches:
            medium_theta_by_batch[batch.label] = fitted
    metric_frames.append(state_metrics("medium_v2_reduced_o2fixed", medium_theta_by_batch, batches, core_cache))

    batch_reduced_theta: dict[str, dict[str, float]] = {}
    batch_all_theta: dict[str, dict[str, float]] = {}
    batch_theta_rows = []
    for batch in batches:
        fitted_reduced, summary_reduced = fit_capacity(
            f"batch_v2_reduced_o2fixed_{batch.medium}_{batch.batch}",
            theta_global,
            [batch],
            core_cache,
            FIT_REDUCED,
            args.max_nfev_batch,
            prior_sigma={},
        )
        batch_reduced_theta[batch.label] = fitted_reduced
        fit_rows.append(summary_reduced)

        fitted_all, summary_all = fit_capacity(
            f"batch_v2_all_o2free_{batch.medium}_{batch.batch}",
            theta_global,
            [batch],
            core_cache,
            FIT_ALL,
            args.max_nfev_batch,
            prior_sigma={},
        )
        batch_all_theta[batch.label] = fitted_all
        fit_rows.append(summary_all)
        row = {"medium": batch.medium, "batch": batch.batch, "label": batch.label}
        for name in FIT_ALL:
            row[name] = fitted_all.get(name, np.nan)
        batch_theta_rows.append(row)

    metric_frames.append(state_metrics("batch_v2_reduced_o2fixed", batch_reduced_theta, batches, core_cache))
    metric_frames.append(state_metrics("batch_v2_all_o2free", batch_all_theta, batches, core_cache))

    fit_summary = pd.DataFrame(fit_rows)
    state_metric_table = pd.concat(metric_frames, ignore_index=True)
    aggregate = aggregate_metrics(state_metric_table)
    batch_theta_df = pd.DataFrame(batch_theta_rows)
    spread = parameter_spread(batch_theta_df, FIT_ALL, "batch_v2_all_o2free")

    fit_summary.to_csv(RESULTS_DIR / "fit_capacity_summary.csv", index=False)
    state_metric_table.to_csv(RESULTS_DIR / "fit_capacity_state_metrics.csv", index=False)
    aggregate.to_csv(RESULTS_DIR / "fit_capacity_aggregate_metrics.csv", index=False)
    batch_theta_df.to_csv(RESULTS_DIR / "batch_all_o2free_parameters.csv", index=False)
    spread.to_csv(RESULTS_DIR / "batch_all_o2free_parameter_spread.csv", index=False)

    plot_capacity(batches, core_cache, theta_global, batch_all_theta)
    write_report(fit_summary, aggregate, spread)
    write_notebook()

    metadata = {
        "n_batches": len(batches),
        "fit_reduced_parameters": list(FIT_REDUCED),
        "fit_all_parameters": list(FIT_ALL),
        "interpretation": "per-batch all-free fit is an upper bound on fit capacity, not an identifiable model",
    }
    (RESULTS_DIR / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("Aggregate metrics:")
    print(aggregate.to_string(index=False))
    print("\nParameter spread:")
    print(spread.to_string(index=False))
    print(f"\n[done] results={RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
