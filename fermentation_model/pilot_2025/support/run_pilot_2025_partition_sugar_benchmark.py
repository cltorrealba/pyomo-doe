from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
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

from shared import aroma_partition_unifac as unifac_partition
from shared import run_new_must_glycerol_estimability_doe as base
import run_pilot_2025_aroma_model_selection_doe as aroma_sel
import run_pilot_2025_calibration_estimability as pilot
from shared import run_secondary_joint_campaign_doe as joint


RESULTS_DIR = PILOT_DIR / "results" / "partition_sugar_benchmark"
PLOT_DIR = RESULTS_DIR / "plots"
PARTITION_MODES = ("water_ethanol", "water_ethanol_total_sugar_as_glucose")
MODE_LABELS = {
    "water_ethanol": "water_ethanol",
    "water_ethanol_total_sugar_as_glucose": "water_ethanol_GF_as_glucose",
}


@dataclass(frozen=True)
class ModeFit:
    mode: str
    theta: dict[str, float]
    fit_summary: pd.DataFrame
    prediction_rows: pd.DataFrame
    pool_metrics: pd.DataFrame


def load_global_theta() -> tuple[dict[str, float], str, aroma_sel.AromaVariant]:
    global_path = PILOT_DIR / "results" / "global_state_model_selection_doe" / "theta_selected_global_model.csv"
    if global_path.exists():
        theta = pd.read_csv(global_path, index_col=0).iloc[:, 0].to_dict()
        theta = {str(k): float(v) for k, v in theta.items() if np.isfinite(float(v))}
    else:
        theta = aroma_sel.load_reference_theta()
    selected_name = "ea_ethanol_nlimited"
    selected_path = aroma_sel.RESULTS_DIR / "selected_model.txt"
    if selected_path.exists():
        selected_name = selected_path.read_text(encoding="utf-8").strip()
    variant = aroma_sel.variant_library().get(selected_name, aroma_sel.variant_library()["ea_ethanol_nlimited"])
    return aroma_sel.clip_theta(theta), selected_name, variant


@lru_cache(maxsize=120000)
def partition_k_cached(
    species: str,
    temp_c_round: float,
    ethanol_g_l_round: float,
    glucose_g_l_round: float,
    fructose_g_l_round: float,
    mode: str,
) -> float:
    result = unifac_partition.unifac_partition_K_antoine(
        species,
        float(temp_c_round),
        float(ethanol_g_l_round),
        glucose_g_l=float(glucose_g_l_round),
        fructose_g_l=float(fructose_g_l_round),
        mode=mode,
        process_min_temp_c=aroma_sel.PARTITION_MIN_TEMP_C,
        process_max_temp_c=aroma_sel.PARTITION_MAX_TEMP_C,
    )
    return max(float(result["K"]), 1e-12)


def partition_k_mode(species: str, temp_c: float, ethanol_g_l: float, glucose_g_l: float, fructose_g_l: float, mode: str) -> float:
    if mode == "water_ethanol":
        glucose_g_l = 0.0
        fructose_g_l = 0.0
    return partition_k_cached(
        species,
        round(float(temp_c), 1),
        round(max(float(ethanol_g_l), 0.0), 1),
        round(max(float(glucose_g_l), 0.0), 1),
        round(max(float(fructose_g_l), 0.0), 1),
        mode,
    )


def aroma_rhs_mode(
    t: float,
    y: np.ndarray,
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    variant: aroma_sel.AromaVariant,
    mode: str,
) -> list[float]:
    liquid = {species: max(float(y[idx]), 0.0) for idx, species in enumerate(joint.AROMA_SPECIES)}
    rates = joint._core_rates(theta, batch, core, float(t))
    phi_growth = rates["N"] / (rates["N"] + joint.N_PHASE_HALF_KG_M3)
    co2_rate = max(joint.CO2_G_PER_G_ETHANOL * rates["ethanol_prod"], 0.0)
    out: list[float] = []
    loss_rates: list[float] = []
    for species in joint.AROMA_SPECIES:
        short = joint.AROMA_SHORT[species]
        k_prod = theta[f"k_{short}_growth"] * phi_growth + theta[f"k_{short}_stationary"] * (1.0 - phi_growth)
        prod = k_prod * max(float(rates["sugar_uptake"]), 0.0)
        if species == "ethyl_acetate":
            prod += aroma_sel.ea_extra_source(variant, theta, rates, secondary, t)
        k_lg = partition_k_mode(
            species,
            rates["TempC"],
            rates["E"],
            rates.get("G", 0.0),
            rates.get("F", 0.0),
            mode,
        )
        loss_rate = theta[f"alpha_{short}_loss"] * k_lg * co2_rate * liquid[species]
        out.append(float(prod - loss_rate))
        loss_rates.append(float(loss_rate))
    out.extend(loss_rates)
    return out


def integrate_aroma_mode(
    batch: base.BatchData,
    theta: dict[str, float],
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    variant: aroma_sel.AromaVariant,
    mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    time = core.index.to_numpy(dtype=float)
    liquid0 = [float(batch.initials.get(f"{species}_liq", 0.0)) for species in joint.AROMA_SPECIES]
    loss0 = [float(batch.initials.get(f"{species}_cond", 0.0)) for species in joint.AROMA_SPECIES]
    y = np.asarray(liquid0 + loss0, dtype=float)
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(aroma_rhs_mode(float(time[idx - 1]), y, theta, batch, core, secondary, variant, mode), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    arr = np.asarray(rows, dtype=float)
    liquid = pd.DataFrame(arr[:, : len(joint.AROMA_SPECIES)], index=time, columns=joint.AROMA_SPECIES)
    condensate = pd.DataFrame(arr[:, len(joint.AROMA_SPECIES) :], index=time, columns=joint.AROMA_SPECIES)
    return liquid, condensate


def prediction_rows_mode(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for batch in batches:
        core = core_cache[batch.label]
        secondary = sec_cache.get(batch.label)
        liq, cond = integrate_aroma_mode(batch, theta, core, secondary, variant, mode)
        for species, total_col in joint.AROMA_COLUMNS.items():
            cond_col = pilot.AROMA_CONDENSATE_COLUMNS[species]
            total_obs = np.asarray(batch.observations.get(total_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_obs = np.asarray(batch.observations.get(cond_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_valid = np.isfinite(cond_obs) & (cond_obs > aroma_sel.CONDENSATE_POSITIVE_THRESHOLD)
            pred_liq = liq.loc[batch.time, species].to_numpy(dtype=float)
            pred_cond = cond.loc[batch.time, species].to_numpy(dtype=float)
            pred_total = pred_liq + pred_cond
            for obs_idx, t in enumerate(batch.time):
                if np.isfinite(total_obs[obs_idx]) and cond_valid[obs_idx]:
                    retained = max(float(total_obs[obs_idx] - cond_obs[obs_idx]), 0.0)
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "retained", "obs": retained, "pred": float(pred_liq[obs_idx])})
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "condensate", "obs": float(cond_obs[obs_idx]), "pred": float(pred_cond[obs_idx])})
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "total", "obs": float(total_obs[obs_idx]), "pred": float(pred_total[obs_idx])})
                elif np.isfinite(total_obs[obs_idx]):
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "total", "obs": float(total_obs[obs_idx]), "pred": float(pred_total[obs_idx])})
                elif cond_valid[obs_idx]:
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "condensate", "obs": float(cond_obs[obs_idx]), "pred": float(pred_cond[obs_idx])})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["residual"] = out["pred"] - out["obs"]
        out["partition_mode"] = MODE_LABELS.get(mode, mode)
    return out


def aroma_residual_mode(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    mode: str,
) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        secondary = sec_cache.get(batch.label)
        liq, cond = integrate_aroma_mode(batch, theta, core, secondary, variant, mode)
        for species, total_col in joint.AROMA_COLUMNS.items():
            species_weight = float(aroma_sel.AROMA_RESIDUAL_WEIGHTS.get(species, 1.0))
            cond_col = pilot.AROMA_CONDENSATE_COLUMNS[species]
            total_obs = np.asarray(batch.observations.get(total_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_obs = np.asarray(batch.observations.get(cond_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_valid = np.isfinite(cond_obs) & (cond_obs > aroma_sel.CONDENSATE_POSITIVE_THRESHOLD)
            pred_liq = liq.loc[batch.time, species].to_numpy(dtype=float)
            pred_cond = cond.loc[batch.time, species].to_numpy(dtype=float)
            pred_total = pred_liq + pred_cond
            mask_both = np.isfinite(total_obs) & cond_valid
            if mask_both.any():
                retained = np.maximum(total_obs[mask_both] - cond_obs[mask_both], 0.0)
                residuals.append(species_weight * (pred_liq[mask_both] - retained) / pilot.aroma_sigma(species, "liquid", retained))
                residuals.append(species_weight * (pred_cond[mask_both] - cond_obs[mask_both]) / pilot.aroma_sigma(species, "condensate", cond_obs[mask_both]))
                residuals.append(species_weight * (pred_total[mask_both] - total_obs[mask_both]) / pilot.aroma_sigma(species, "total", total_obs[mask_both]))
            mask_total_only = np.isfinite(total_obs) & ~cond_valid
            if mask_total_only.any():
                residuals.append(species_weight * (pred_total[mask_total_only] - total_obs[mask_total_only]) / pilot.aroma_sigma(species, "total", total_obs[mask_total_only]))
            mask_cond_only = ~np.isfinite(total_obs) & cond_valid
            if mask_cond_only.any():
                residuals.append(species_weight * (pred_cond[mask_cond_only] - cond_obs[mask_cond_only]) / pilot.aroma_sigma(species, "condensate", cond_obs[mask_cond_only]))
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def fit_mode(
    mode: str,
    theta0: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> ModeFit:
    parameters = variant.parameters
    lb, ub = aroma_sel.log_bounds(parameters)
    x_ref = np.clip(aroma_sel.log_vector(theta0, parameters), lb + 1e-9, ub - 1e-9)
    rng = np.random.default_rng(seed)
    rows = []
    best_theta = dict(theta0)
    fixed_data_res = aroma_residual_mode(theta0, batches, variant, core_cache, sec_cache, mode)
    fixed_data_wsse = float(np.dot(fixed_data_res, fixed_data_res))
    fixed_active_bounds = aroma_sel.active_bound_count(theta0, parameters)
    best_wsse = fixed_data_wsse + 20.0 * float(fixed_active_bounds)
    best_pred = prediction_rows_mode(theta0, batches, variant, core_cache, sec_cache, mode)
    rows.append(
        {
            "partition_mode": MODE_LABELS.get(mode, mode),
            "start": "fixed_reference",
            "success": True,
            "status": 0,
            "message": "selected global theta before aroma-only refit",
            "nfev": 0,
            "initial_objective_wsse": fixed_data_wsse,
            "final_objective_wsse": fixed_data_wsse,
            "data_wsse": fixed_data_wsse,
            "n_data_residuals": int(len(fixed_data_res)),
            "active_bound_count": fixed_active_bounds,
            "fit_selection_score": best_wsse,
        }
    )

    def fun(x: np.ndarray) -> np.ndarray:
        theta = aroma_sel.theta_from_log(x, theta0, parameters)
        data_res = aroma_residual_mode(theta, batches, variant, core_cache, sec_cache, mode)
        priors = []
        for name in parameters:
            scale = 1.15 if name.startswith("alpha_") else 1.8
            if name.startswith("k_EA_") or name == "q10_EA":
                scale = 4.0
            ref = max(float(theta0[name]), 1e-16)
            priors.append(math.log(max(float(theta[name]), 1e-16) / ref) / scale)
        return np.concatenate([data_res, np.asarray(priors, dtype=float)])

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = x_ref.copy()
        else:
            x0 = x_ref + rng.normal(0.0, 0.75, size=len(parameters))
            if rng.random() < 0.35:
                x0 = rng.uniform(lb, ub)
            x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)
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
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
        )
        theta_hat = aroma_sel.theta_from_log(result.x, theta0, parameters)
        end_res = fun(result.x)
        data_res = aroma_residual_mode(theta_hat, batches, variant, core_cache, sec_cache, mode)
        data_wsse = float(np.dot(data_res, data_res))
        rows.append(
            {
                "partition_mode": MODE_LABELS.get(mode, mode),
                "start": idx,
                "success": bool(result.success),
                "status": int(result.status),
                "message": str(result.message),
                "nfev": int(result.nfev),
                "initial_objective_wsse": float(np.dot(start_res, start_res)),
                "final_objective_wsse": float(np.dot(end_res, end_res)),
                "data_wsse": data_wsse,
                "n_data_residuals": int(len(data_res)),
                "active_bound_count": aroma_sel.active_bound_count(theta_hat, parameters),
            }
        )
        fit_score = data_wsse + 20.0 * rows[-1]["active_bound_count"]
        rows[-1]["fit_selection_score"] = fit_score
        if fit_score < best_wsse:
            best_wsse = fit_score
            best_theta = theta_hat
            best_pred = prediction_rows_mode(theta_hat, batches, variant, core_cache, sec_cache, mode)
    fit_summary = pd.DataFrame(rows)
    pool_metrics = aroma_sel.summarize_prediction_metrics(best_pred)
    if not pool_metrics.empty:
        pool_metrics["partition_mode"] = MODE_LABELS.get(mode, mode)
    return ModeFit(mode=mode, theta=best_theta, fit_summary=fit_summary, prediction_rows=best_pred, pool_metrics=pool_metrics)


def fixed_mode_result(
    mode: str,
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
) -> ModeFit:
    pred = prediction_rows_mode(theta, batches, variant, core_cache, sec_cache, mode)
    res = aroma_residual_mode(theta, batches, variant, core_cache, sec_cache, mode)
    metrics = aroma_sel.summarize_prediction_metrics(pred)
    if not metrics.empty:
        metrics["partition_mode"] = MODE_LABELS.get(mode, mode)
    fit = pd.DataFrame(
        [
            {
                "partition_mode": MODE_LABELS.get(mode, mode),
                "start": "fixed_theta",
                "success": True,
                "data_wsse": float(np.dot(res, res)),
                "n_data_residuals": int(len(res)),
                "active_bound_count": aroma_sel.active_bound_count(theta, variant.parameters),
            }
        ]
    )
    return ModeFit(mode=mode, theta=dict(theta), fit_summary=fit, prediction_rows=pred, pool_metrics=metrics)


def partition_audit_from_process(theta: dict[str, float], batches: list[base.BatchData], core_cache: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            continue
        for t in np.linspace(float(core.index.min()), float(core.index.max()), 16):
            rates = joint._core_rates(theta, batch, core, float(t))
            for species in joint.AROMA_SPECIES:
                k_we = partition_k_mode(species, rates["TempC"], rates["E"], rates.get("G", 0.0), rates.get("F", 0.0), "water_ethanol")
                k_sugar = partition_k_mode(
                    species,
                    rates["TempC"],
                    rates["E"],
                    rates.get("G", 0.0),
                    rates.get("F", 0.0),
                    "water_ethanol_total_sugar_as_glucose",
                )
                rows.append(
                    {
                        "batch": batch.batch,
                        "time_h": float(t),
                        "species": species,
                        "temperature_c": float(rates["TempC"]),
                        "ethanol_g_l": float(rates["E"]),
                        "glucose_g_l": float(rates.get("G", 0.0)),
                        "fructose_g_l": float(rates.get("F", 0.0)),
                        "K_water_ethanol": k_we,
                        "K_GF_as_glucose": k_sugar,
                        "K_ratio_GF_as_glucose_over_water_ethanol": k_sugar / k_we if k_we > 0 else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def plot_metric_comparison(metrics: pd.DataFrame) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    if metrics.empty:
        return
    for pool in sorted(metrics["pool"].dropna().unique()):
        sub = metrics[metrics["pool"].eq(pool)].copy()
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(9, 4.5))
        labels = sorted(sub["species"].unique())
        x = np.arange(len(labels))
        width = 0.35
        for offset, mode in zip((-0.5, 0.5), sorted(sub["partition_mode"].unique())):
            vals = [
                float(sub[sub["species"].eq(species) & sub["partition_mode"].eq(mode)]["relative_rmse_to_median"].iloc[0])
                if not sub[sub["species"].eq(species) & sub["partition_mode"].eq(mode)].empty
                else np.nan
                for species in labels
            ]
            ax.bar(x + offset * width, vals, width=width, label=mode)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        ax.set_ylabel("relative RMSE to median")
        ax.set_title(f"Partition benchmark: {pool}")
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"partition_benchmark_{pool}.png", dpi=170)
        plt.close(fig)


def write_report(selected_model: str, fixed_summary: pd.DataFrame, refit_summary: pd.DataFrame, audit_summary: pd.DataFrame) -> None:
    path = RESULTS_DIR / "partition_sugar_benchmark_report.md"
    with path.open("w", encoding="utf-8") as f:
        f.write("# Water-ethanol vs water-ethanol-sugar partition benchmark\n\n")
        f.write(f"Inherited aroma model: `{selected_model}`.\n\n")
        f.write(
            "The sugar-aware mode uses original UNIFAC with water, ethanol, and total `G+F` treated as glucose-equivalent. "
            "Antoine vapor pressure is unchanged. Only the gas-liquid partition coefficient changes:\n\n"
        )
        f.write(
            "$$K_i=\\frac{\\gamma_i^{UNIFAC}(T,x_{water},x_{ethanol},x_{sugar})P_i^{sat,Antoine}(T)}{RTC_{tot,L}}.$$\n\n"
        )
        f.write("## Partition coefficient effect over pilot trajectories\n\n")
        f.write(audit_summary.to_markdown(index=False))
        f.write("\n\n## Fixed-parameter prediction comparison\n\n")
        f.write(fixed_summary.to_markdown(index=False))
        f.write("\n\n## Refit aroma-parameter comparison\n\n")
        f.write(refit_summary.to_markdown(index=False))
        f.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark sugar-aware UNIFAC partitioning for pilot aroma fits.")
    parser.add_argument("--n-starts", type=int, default=4)
    parser.add_argument("--max-nfev", type=int, default=160)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    print("[load] data and selected global theta", flush=True)
    model_data, _curated_co2, _decisions = pilot.load_pilot_model_data()
    batches = pilot.make_extended_batches(model_data)
    theta0, selected_model, variant = load_global_theta()
    core_cache = joint.precompute_core_cache(batches, theta0)
    sec_cache = aroma_sel.secondary_cache_for(batches, theta0, core_cache)

    audit = partition_audit_from_process(theta0, batches, core_cache)
    audit.to_csv(RESULTS_DIR / "partition_K_process_audit.csv", index=False)
    audit_summary = (
        audit.groupby("species")["K_ratio_GF_as_glucose_over_water_ethanol"]
        .agg(["count", "min", "median", "max"])
        .reset_index()
        .rename(columns={"count": "n_process_points", "median": "median_K_ratio"})
    )
    audit_summary.to_csv(RESULTS_DIR / "partition_K_ratio_summary.csv", index=False)

    print("[fixed] compare modes with selected theta", flush=True)
    fixed_results = []
    for mode in PARTITION_MODES:
        result = fixed_mode_result(mode, theta0, batches, variant, core_cache, sec_cache)
        label = MODE_LABELS[mode]
        result.fit_summary.to_csv(RESULTS_DIR / f"fixed_fit_summary_{label}.csv", index=False)
        result.prediction_rows.to_csv(RESULTS_DIR / f"fixed_predictions_{label}.csv", index=False)
        result.pool_metrics.to_csv(RESULTS_DIR / f"fixed_pool_metrics_{label}.csv", index=False)
        fixed_results.append(result)
    fixed_summary = pd.concat([r.fit_summary for r in fixed_results], ignore_index=True, sort=False)
    fixed_metrics = pd.concat([r.pool_metrics for r in fixed_results], ignore_index=True, sort=False)
    fixed_summary.to_csv(RESULTS_DIR / "fixed_fit_summary.csv", index=False)
    fixed_metrics.to_csv(RESULTS_DIR / "fixed_pool_metrics.csv", index=False)

    print("[refit] compare modes with selected aroma parameters free", flush=True)
    refit_results = []
    for idx, mode in enumerate(PARTITION_MODES, start=1):
        print(f"[refit] {MODE_LABELS[mode]}", flush=True)
        result = fit_mode(
            mode,
            theta0,
            batches,
            variant,
            core_cache,
            sec_cache,
            n_starts=args.n_starts,
            max_nfev=args.max_nfev,
            seed=20260629 + idx,
        )
        label = MODE_LABELS[mode]
        result.fit_summary.to_csv(RESULTS_DIR / f"refit_multistart_{label}.csv", index=False)
        result.prediction_rows.to_csv(RESULTS_DIR / f"refit_predictions_{label}.csv", index=False)
        result.pool_metrics.to_csv(RESULTS_DIR / f"refit_pool_metrics_{label}.csv", index=False)
        pd.Series({name: result.theta[name] for name in variant.parameters}).to_csv(RESULTS_DIR / f"theta_refit_{label}.csv")
        refit_results.append(result)
    refit_summary = pd.concat([r.fit_summary for r in refit_results], ignore_index=True, sort=False)
    best_rows = refit_summary.sort_values("fit_selection_score").groupby("partition_mode", as_index=False).first()
    refit_metrics = pd.concat([r.pool_metrics for r in refit_results], ignore_index=True, sort=False)
    refit_summary.to_csv(RESULTS_DIR / "refit_multistart_summary.csv", index=False)
    best_rows.to_csv(RESULTS_DIR / "refit_best_summary.csv", index=False)
    refit_metrics.to_csv(RESULTS_DIR / "refit_pool_metrics.csv", index=False)

    plot_metric_comparison(fixed_metrics.assign(case="fixed"))
    plot_metric_comparison(refit_metrics.assign(case="refit"))
    write_report(selected_model, fixed_summary, best_rows, audit_summary)
    print("[done]", RESULTS_DIR, flush=True)


if __name__ == "__main__":
    main()
