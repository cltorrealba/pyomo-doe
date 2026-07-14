from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
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

from shared import run_new_must_glycerol_estimability_doe as base
import run_pilot_2025_aroma_model_selection_doe as aroma_sel
import run_pilot_2025_calibration_estimability as pilot
import run_pilot_2025_partition_sugar_benchmark as sugar_bench
from shared import run_secondary_joint_campaign_doe as joint


RESULTS_DIR = PILOT_DIR / "results" / "co2_stripping_benchmark"
PLOT_DIR = RESULTS_DIR / "plots"

PARTITION_MODES = ("water_ethanol", "water_ethanol_total_sugar_as_glucose")
PARTITION_LABELS = {
    "water_ethanol": "water_ethanol",
    "water_ethanol_total_sugar_as_glucose": "water_ethanol_GF_as_glucose",
}


@dataclass(frozen=True)
class CO2Mode:
    name: str
    parameters: tuple[str, ...]
    defaults: dict[str, float]
    bounds: dict[str, tuple[float, float]]
    description: str


def co2_mode_library() -> dict[str, CO2Mode]:
    common_bounds = {
        "tau_h": (0.25, 96.0),
        "power_p": (0.35, 2.50),
        "threshold_fraction": (0.0, 0.90),
    }
    return {
        "instant": CO2Mode(
            name="instant",
            parameters=tuple(),
            defaults={},
            bounds={},
            description="Current model: stripping gas flow is directly proportional to instantaneous ethanol-derived CO2 production.",
        ),
        "lag": CO2Mode(
            name="lag",
            parameters=("tau_h",),
            defaults={"tau_h": 8.0},
            bounds={"tau_h": common_bounds["tau_h"]},
            description="First-order gas-release/headspace lag: tau dq/dt = q_prod - q_out.",
        ),
        "power": CO2Mode(
            name="power",
            parameters=("power_p",),
            defaults={"power_p": 1.0},
            bounds={"power_p": common_bounds["power_p"]},
            description="Nonlinear effective gas flow: q_out = q_ref (q_prod/q_ref)^p.",
        ),
        "threshold": CO2Mode(
            name="threshold",
            parameters=("threshold_fraction",),
            defaults={"threshold_fraction": 0.10},
            bounds={"threshold_fraction": common_bounds["threshold_fraction"]},
            description="Gas flow activates only above a fraction of the median production rate.",
        ),
        "lag_power": CO2Mode(
            name="lag_power",
            parameters=("tau_h", "power_p"),
            defaults={"tau_h": 8.0, "power_p": 1.0},
            bounds={"tau_h": common_bounds["tau_h"], "power_p": common_bounds["power_p"]},
            description="Nonlinear production-rate transform followed by first-order gas-release lag.",
        ),
        "lag_threshold": CO2Mode(
            name="lag_threshold",
            parameters=("tau_h", "threshold_fraction"),
            defaults={"tau_h": 8.0, "threshold_fraction": 0.10},
            bounds={"tau_h": common_bounds["tau_h"], "threshold_fraction": common_bounds["threshold_fraction"]},
            description="Thresholded production-rate transform followed by first-order gas-release lag.",
        ),
    }


def load_global_theta() -> tuple[dict[str, float], str, aroma_sel.AromaVariant]:
    return sugar_bench.load_global_theta()


def positive_reference(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    positive = values[np.isfinite(values) & (values > 1e-12)]
    if positive.size == 0:
        return 1.0
    return max(float(np.nanmedian(positive)), 1e-9)


def apply_lag(time: np.ndarray, input_rate: np.ndarray, tau_h: float) -> np.ndarray:
    time = np.asarray(time, dtype=float)
    input_rate = np.asarray(input_rate, dtype=float)
    tau = max(float(tau_h), 1e-6)
    out = np.zeros_like(input_rate, dtype=float)
    if len(out) == 0:
        return out
    out[0] = max(float(input_rate[0]), 0.0)
    for idx in range(1, len(out)):
        dt = max(float(time[idx] - time[idx - 1]), 0.0)
        decay = math.exp(-dt / tau) if tau > 0 else 0.0
        out[idx] = out[idx - 1] * decay + max(float(input_rate[idx - 1]), 0.0) * (1.0 - decay)
    return np.maximum(out, 0.0)


def base_co2_production(theta: dict[str, float], batch: base.BatchData, core: pd.DataFrame, time: np.ndarray) -> np.ndarray:
    values = []
    for t in np.asarray(time, dtype=float):
        rates = joint._core_rates(theta, batch, core, float(t))
        values.append(joint.CO2_G_PER_G_ETHANOL * max(float(rates["ethanol_prod"]), 0.0))
    return np.asarray(values, dtype=float)


def gas_flow_from_production(
    time: np.ndarray,
    q_prod: np.ndarray,
    mode: CO2Mode,
    params: dict[str, float],
) -> np.ndarray:
    q = np.maximum(np.asarray(q_prod, dtype=float), 0.0)
    if "threshold_fraction" in mode.parameters:
        threshold = float(params.get("threshold_fraction", mode.defaults.get("threshold_fraction", 0.0))) * positive_reference(q)
        q = np.maximum(q - threshold, 0.0)
    if "power_p" in mode.parameters:
        p = float(params.get("power_p", mode.defaults.get("power_p", 1.0)))
        ref = positive_reference(q)
        q = ref * np.power(np.maximum(q, 0.0) / ref, p)
    if "tau_h" in mode.parameters:
        q = apply_lag(time, q, float(params.get("tau_h", mode.defaults.get("tau_h", 8.0))))
    return np.maximum(q, 0.0)


def co2_prediction_for_batch(
    theta: dict[str, float],
    batch: base.BatchData,
    sensor_times: np.ndarray,
    mode: CO2Mode,
    params: dict[str, float],
    dt_grid_h: float = 1.0,
) -> np.ndarray:
    horizon = max(float(np.nanmax(batch.time)), float(np.nanmax(sensor_times)))
    grid = np.arange(0.0, horizon + dt_grid_h, dt_grid_h)
    sim_times = np.asarray(sorted(set(np.round(grid, 8)).union(set(np.round(sensor_times, 8))).union({0.0, horizon})), dtype=float)
    core = base.simulate(batch, theta, sim_times)
    if core is None:
        return np.full_like(sensor_times, np.nan, dtype=float)
    q_prod = base_co2_production(theta, batch, core, sim_times)
    q_gas = gas_flow_from_production(sim_times, q_prod, mode, params)
    return np.interp(np.asarray(sensor_times, dtype=float), sim_times, q_gas)


def co2_shape_residual(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    mode: CO2Mode,
    params: dict[str, float],
) -> np.ndarray:
    residuals = []
    for batch_name, group in co2_down.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        obs = group["co2_raw"].to_numpy(dtype=float)
        pred = co2_prediction_for_batch(theta, batch, times, mode, params)
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


def co2_metrics_by_batch(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    mode: CO2Mode,
    params: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    for batch_name, group in co2_down.groupby("batch", sort=True):
        batch = batches_by_name.get(str(batch_name))
        if batch is None:
            continue
        times = group["time_h_effective"].to_numpy(dtype=float)
        obs = group["co2_raw"].to_numpy(dtype=float)
        pred_raw = co2_prediction_for_batch(theta, batch, times, mode, params)
        mask = np.isfinite(obs) & np.isfinite(pred_raw) & (obs >= 0.0)
        if mask.sum() < 5:
            continue
        obs = obs[mask]
        pred_raw = pred_raw[mask]
        times = times[mask]
        denom = float(np.dot(pred_raw, pred_raw))
        scale = float(np.dot(obs, pred_raw) / denom) if denom > 1e-12 else 0.0
        pred = scale * pred_raw
        err = pred - obs
        corr = np.nan
        if np.nanstd(obs) > 1e-12 and np.nanstd(pred) > 1e-12:
            corr = float(np.corrcoef(obs, pred)[0, 1])
        obs_scale = max(float(np.nanmedian(np.abs(obs))), 0.25 * float(np.nanmax(obs) - np.nanmin(obs)), 1e-6)
        rows.append(
            {
                "mode": mode.name,
                "batch": str(batch_name),
                "n": int(len(obs)),
                "scale_factor": scale,
                "rmse": float(np.sqrt(np.mean(err * err))),
                "mae": float(np.mean(np.abs(err))),
                "bias": float(np.mean(err)),
                "relative_rmse": float(np.sqrt(np.mean(err * err)) / obs_scale),
                "relative_bias": float(np.mean(err) / obs_scale),
                "corr": corr,
            }
        )
        for t, y, yhat in zip(times, obs, pred):
            pred_rows.append({"mode": mode.name, "batch": str(batch_name), "time_h": float(t), "obs": float(y), "pred": float(yhat)})
    return pd.DataFrame(rows), pd.DataFrame(pred_rows)


def fit_co2_mode(
    theta: dict[str, float],
    batches_by_name: dict[str, base.BatchData],
    co2_down: pd.DataFrame,
    mode: CO2Mode,
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    if not mode.parameters:
        params = dict(mode.defaults)
        res = co2_shape_residual(theta, batches_by_name, co2_down, mode, params)
        return params, pd.DataFrame(
            [
                {
                    "mode": mode.name,
                    "start": 0,
                    "success": True,
                    "nfev": 0,
                    "data_wsse": float(np.dot(res, res)),
                    "n_data_residuals": int(len(res)),
                    "active_bound_count": 0,
                    "params_json": "{}",
                }
            ]
        )

    lb = np.asarray([mode.bounds[name][0] for name in mode.parameters], dtype=float)
    ub = np.asarray([mode.bounds[name][1] for name in mode.parameters], dtype=float)
    x_ref = np.asarray([mode.defaults[name] for name in mode.parameters], dtype=float)
    rng = np.random.default_rng(seed)
    rows = []
    best_score = np.inf
    best_params = dict(mode.defaults)

    def make_params(x: np.ndarray) -> dict[str, float]:
        return {name: float(value) for name, value in zip(mode.parameters, np.asarray(x, dtype=float))}

    def fun(x: np.ndarray) -> np.ndarray:
        params = make_params(x)
        return co2_shape_residual(theta, batches_by_name, co2_down, mode, params)

    for idx in range(int(n_starts)):
        if idx == 0:
            x0 = np.clip(x_ref, lb + 1e-9, ub - 1e-9)
        else:
            x0 = np.empty_like(x_ref)
            for j, name in enumerate(mode.parameters):
                lo, hi = mode.bounds[name]
                if name == "tau_h":
                    x0[j] = math.exp(rng.uniform(math.log(lo), math.log(hi)))
                else:
                    x0[j] = rng.uniform(lo, hi)
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
        params = make_params(result.x)
        res = co2_shape_residual(theta, batches_by_name, co2_down, mode, params)
        wsse = float(np.dot(res, res))
        active = 0
        for name, value in params.items():
            lo, hi = mode.bounds[name]
            if value <= lo * 1.01 or value >= hi / 1.01:
                active += 1
        score = wsse + 10.0 * active
        rows.append(
            {
                "mode": mode.name,
                "start": idx,
                "success": bool(result.success),
                "status": int(result.status),
                "message": str(result.message),
                "nfev": int(result.nfev),
                "data_wsse": wsse,
                "n_data_residuals": int(len(res)),
                "active_bound_count": int(active),
                "fit_selection_score": score,
                "params_json": str(params),
                **params,
            }
        )
        if score < best_score:
            best_score = score
            best_params = params
    return best_params, pd.DataFrame(rows)


def gas_flow_cache_for_batches(
    theta: dict[str, float],
    batches: list[base.BatchData],
    core_cache: dict[str, pd.DataFrame],
    mode: CO2Mode,
    params: dict[str, float],
) -> dict[str, pd.Series]:
    out: dict[str, pd.Series] = {}
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            continue
        time = core.index.to_numpy(dtype=float)
        q_prod = base_co2_production(theta, batch, core, time)
        q_gas = gas_flow_from_production(time, q_prod, mode, params)
        out[batch.label] = pd.Series(q_gas, index=time)
    return out


def interp_series(series: pd.Series | None, t: float) -> float:
    if series is None or series.empty:
        return 0.0
    index = series.index.to_numpy(dtype=float)
    values = series.to_numpy(dtype=float)
    return float(np.interp(float(t), index, values))


def aroma_rhs_combo(
    t: float,
    y: np.ndarray,
    theta: dict[str, float],
    batch: base.BatchData,
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    qgas: pd.Series,
    variant: aroma_sel.AromaVariant,
    partition_mode: str,
) -> list[float]:
    liquid = {species: max(float(y[idx]), 0.0) for idx, species in enumerate(joint.AROMA_SPECIES)}
    rates = joint._core_rates(theta, batch, core, float(t))
    phi_growth = rates["N"] / (rates["N"] + joint.N_PHASE_HALF_KG_M3)
    co2_flow = max(interp_series(qgas, t), 0.0)
    out: list[float] = []
    loss_rates: list[float] = []
    for species in joint.AROMA_SPECIES:
        short = joint.AROMA_SHORT[species]
        k_prod = theta[f"k_{short}_growth"] * phi_growth + theta[f"k_{short}_stationary"] * (1.0 - phi_growth)
        prod = k_prod * max(float(rates["sugar_uptake"]), 0.0)
        if species == "ethyl_acetate":
            prod += aroma_sel.ea_extra_source(variant, theta, rates, secondary, t)
        k_lg = sugar_bench.partition_k_mode(
            species,
            rates["TempC"],
            rates["E"],
            rates.get("G", 0.0),
            rates.get("F", 0.0),
            partition_mode,
        )
        loss_rate = theta[f"alpha_{short}_loss"] * k_lg * co2_flow * liquid[species]
        out.append(float(prod - loss_rate))
        loss_rates.append(float(loss_rate))
    out.extend(loss_rates)
    return out


def integrate_aroma_combo(
    batch: base.BatchData,
    theta: dict[str, float],
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    qgas: pd.Series,
    variant: aroma_sel.AromaVariant,
    partition_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    time = core.index.to_numpy(dtype=float)
    liquid0 = [float(batch.initials.get(f"{species}_liq", 0.0)) for species in joint.AROMA_SPECIES]
    loss0 = [float(batch.initials.get(f"{species}_cond", 0.0)) for species in joint.AROMA_SPECIES]
    y = np.asarray(liquid0 + loss0, dtype=float)
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(aroma_rhs_combo(float(time[idx - 1]), y, theta, batch, core, secondary, qgas, variant, partition_mode), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    arr = np.asarray(rows, dtype=float)
    liquid = pd.DataFrame(arr[:, : len(joint.AROMA_SPECIES)], index=time, columns=joint.AROMA_SPECIES)
    condensate = pd.DataFrame(arr[:, len(joint.AROMA_SPECIES) :], index=time, columns=joint.AROMA_SPECIES)
    return liquid, condensate


def aroma_prediction_rows_combo(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    qgas_cache: dict[str, pd.Series],
    partition_mode: str,
    co2_mode_name: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for batch in batches:
        core = core_cache[batch.label]
        secondary = sec_cache.get(batch.label)
        qgas = qgas_cache.get(batch.label)
        if qgas is None:
            continue
        liq, cond = integrate_aroma_combo(batch, theta, core, secondary, qgas, variant, partition_mode)
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
        out["co2_mode"] = co2_mode_name
        out["partition_mode"] = PARTITION_LABELS.get(partition_mode, partition_mode)
    return out


def aroma_residual_combo(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    qgas_cache: dict[str, pd.Series],
    partition_mode: str,
) -> np.ndarray:
    pred = aroma_prediction_rows_combo(theta, batches, variant, core_cache, sec_cache, qgas_cache, partition_mode, "tmp")
    if pred.empty:
        return np.array([], dtype=float)
    residuals = []
    for (species, pool), group in pred.groupby(["species", "pool"], sort=True):
        species_weight = float(aroma_sel.AROMA_RESIDUAL_WEIGHTS.get(species, 1.0))
        obs = group["obs"].to_numpy(dtype=float)
        err = group["pred"].to_numpy(dtype=float) - obs
        if pool == "retained":
            sigma = pilot.aroma_sigma(species, "liquid", obs)
        elif pool == "condensate":
            sigma = pilot.aroma_sigma(species, "condensate", obs)
        else:
            sigma = pilot.aroma_sigma(species, "total", obs)
        residuals.append(species_weight * err / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def evaluate_aroma_combo(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: aroma_sel.AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    qgas_cache: dict[str, pd.Series],
    co2_mode_name: str,
    partition_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, float, int]:
    pred = aroma_prediction_rows_combo(theta, batches, variant, core_cache, sec_cache, qgas_cache, partition_mode, co2_mode_name)
    metrics = aroma_sel.summarize_prediction_metrics(pred)
    if not metrics.empty:
        metrics["co2_mode"] = co2_mode_name
        metrics["partition_mode"] = PARTITION_LABELS.get(partition_mode, partition_mode)
    res = aroma_residual_combo(theta, batches, variant, core_cache, sec_cache, qgas_cache, partition_mode)
    return pred, metrics, float(np.dot(res, res)), int(len(res))


def plot_co2_fits(pred_rows: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if pred_rows.empty:
        return
    for batch, group in pred_rows.groupby("batch", sort=True):
        fig, ax = plt.subplots(figsize=(9, 4))
        for mode, sub in group.groupby("mode", sort=True):
            sub = sub.sort_values("time_h")
            if mode == group["mode"].iloc[0]:
                ax.scatter(sub["time_h"], sub["obs"], s=10, color="black", label="obs")
            ax.plot(sub["time_h"], sub["pred"], lw=1.3, label=mode)
        ax.set_title(f"CO2 shape fit: {batch}")
        ax.set_xlabel("effective time (h)")
        ax.set_ylabel("scaled CO2 signal")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output_dir / f"co2_fit_{batch}.png", dpi=170)
        plt.close(fig)


def write_report(
    co2_summary: pd.DataFrame,
    co2_metrics: pd.DataFrame,
    aroma_summary: pd.DataFrame,
    selected_co2_mode: str,
    selected_partition_mode: str,
) -> None:
    path = RESULTS_DIR / "co2_stripping_benchmark_report.md"
    with path.open("w", encoding="utf-8") as f:
        f.write("# CO2 stripping and sugar-aware partition benchmark\n\n")
        f.write("## CO2 model ranking\n\n")
        f.write(co2_summary.to_markdown(index=False))
        f.write("\n\n## CO2 metrics by batch\n\n")
        f.write(co2_metrics.to_markdown(index=False))
        f.write("\n\n## Aroma fixed-parameter benchmark\n\n")
        f.write(aroma_summary.to_markdown(index=False))
        f.write("\n\n")
        f.write(f"Selected CO2 mode by CO2 BIC: `{selected_co2_mode}`.\n\n")
        f.write(f"Best fixed aroma/partition combination in this run: `{selected_co2_mode}` with `{selected_partition_mode}` if it also ranks best in the aroma table; otherwise inspect the table before adopting it.\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark CO2 gas-flow and sugar-aware partition models.")
    parser.add_argument("--n-starts-co2", type=int, default=8)
    parser.add_argument("--max-nfev-co2", type=int, default=180)
    parser.add_argument("--top-co2-for-aroma", type=int, default=4)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    print("[load] pilot data, selected global theta", flush=True)
    model_data, curated_co2, co2_decisions = pilot.load_pilot_model_data()
    co2_down = pilot.downsample_co2(curated_co2)
    batches = pilot.make_extended_batches(model_data)
    batches_by_name = {batch.batch: batch for batch in batches}
    theta, selected_aroma_name, variant = load_global_theta()
    core_cache = joint.precompute_core_cache(batches, theta)
    sec_cache = aroma_sel.secondary_cache_for(batches, theta, core_cache)
    co2_decisions.to_csv(RESULTS_DIR / "co2_curation_decisions.csv", index=False)
    co2_down.to_csv(RESULTS_DIR / "co2_downsampled.csv", index=False)

    modes = co2_mode_library()
    best_params_by_mode: dict[str, dict[str, float]] = {}
    all_fit_rows = []
    all_metric_rows = []
    all_pred_rows = []
    summary_rows = []
    print("[co2] fitting shape models", flush=True)
    for idx, mode in enumerate(modes.values(), start=1):
        print(f"[co2] {idx}/{len(modes)} {mode.name}", flush=True)
        params, fit_rows = fit_co2_mode(
            theta,
            batches_by_name,
            co2_down,
            mode,
            n_starts=args.n_starts_co2,
            max_nfev=args.max_nfev_co2,
            seed=20260629 + idx,
        )
        best_params_by_mode[mode.name] = params
        fit_rows.to_csv(RESULTS_DIR / f"co2_fit_multistart_{mode.name}.csv", index=False)
        metrics, pred_rows = co2_metrics_by_batch(theta, batches_by_name, co2_down, mode, params)
        metrics.to_csv(RESULTS_DIR / f"co2_metrics_{mode.name}.csv", index=False)
        pred_rows.to_csv(RESULTS_DIR / f"co2_predictions_{mode.name}.csv", index=False)
        all_fit_rows.append(fit_rows)
        all_metric_rows.append(metrics)
        all_pred_rows.append(pred_rows)
        best_fit = fit_rows.sort_values("data_wsse").iloc[0]
        n = int(best_fit["n_data_residuals"])
        k = len(mode.parameters)
        bic = float(best_fit["data_wsse"] + max(k, 1) * math.log(max(n, 1)))
        summary_rows.append(
            {
                "mode": mode.name,
                "data_wsse": float(best_fit["data_wsse"]),
                "n_data_residuals": n,
                "n_parameters": k,
                "active_bound_count": int(best_fit.get("active_bound_count", 0)),
                "bic": bic,
                "params": params,
                "description": mode.description,
            }
        )

    fit_summary = pd.concat(all_fit_rows, ignore_index=True, sort=False)
    co2_metrics = pd.concat(all_metric_rows, ignore_index=True, sort=False)
    co2_preds = pd.concat(all_pred_rows, ignore_index=True, sort=False)
    co2_summary = pd.DataFrame(summary_rows).sort_values(["bic", "data_wsse"]).reset_index(drop=True)
    co2_summary["selected_by_bic"] = co2_summary.index == 0
    fit_summary.to_csv(RESULTS_DIR / "co2_fit_multistart_summary.csv", index=False)
    co2_metrics.to_csv(RESULTS_DIR / "co2_metrics_by_batch.csv", index=False)
    co2_preds.to_csv(RESULTS_DIR / "co2_predictions_all_modes.csv", index=False)
    co2_summary.to_csv(RESULTS_DIR / "co2_model_selection_summary.csv", index=False)
    plot_co2_fits(co2_preds, PLOT_DIR / "co2")

    selected_modes = list(co2_summary.head(int(args.top_co2_for_aroma))["mode"].astype(str))
    if "instant" not in selected_modes:
        selected_modes.append("instant")

    print("[aroma] fixed-parameter benchmark with top CO2 modes and partition modes", flush=True)
    aroma_summary_rows = []
    aroma_metrics_rows = []
    for co2_mode_name in selected_modes:
        mode = modes[co2_mode_name]
        params = best_params_by_mode[co2_mode_name]
        qgas_cache = gas_flow_cache_for_batches(theta, batches, core_cache, mode, params)
        for partition_mode in PARTITION_MODES:
            pred, metrics, aroma_wsse, nres = evaluate_aroma_combo(
                theta,
                batches,
                variant,
                core_cache,
                sec_cache,
                qgas_cache,
                co2_mode_name,
                partition_mode,
            )
            label = f"{co2_mode_name}__{PARTITION_LABELS[partition_mode]}"
            pred.to_csv(RESULTS_DIR / f"aroma_predictions_{label}.csv", index=False)
            metrics.to_csv(RESULTS_DIR / f"aroma_pool_metrics_{label}.csv", index=False)
            aroma_metrics_rows.append(metrics)
            co2_wsse = float(co2_summary.loc[co2_summary["mode"].eq(co2_mode_name), "data_wsse"].iloc[0])
            aroma_summary_rows.append(
                {
                    "co2_mode": co2_mode_name,
                    "partition_mode": PARTITION_LABELS[partition_mode],
                    "aroma_data_wsse": aroma_wsse,
                    "aroma_n_residuals": nres,
                    "co2_data_wsse": co2_wsse,
                    "combined_wsse": aroma_wsse + co2_wsse,
                    "co2_params": params,
                }
            )
    aroma_summary = pd.DataFrame(aroma_summary_rows).sort_values(["combined_wsse", "aroma_data_wsse"]).reset_index(drop=True)
    aroma_metrics = pd.concat(aroma_metrics_rows, ignore_index=True, sort=False) if aroma_metrics_rows else pd.DataFrame()
    aroma_summary.to_csv(RESULTS_DIR / "aroma_co2_partition_fixed_summary.csv", index=False)
    aroma_metrics.to_csv(RESULTS_DIR / "aroma_co2_partition_pool_metrics.csv", index=False)

    selected_co2 = str(co2_summary.iloc[0]["mode"])
    selected_partition = str(aroma_summary.iloc[0]["partition_mode"]) if not aroma_summary.empty else ""
    write_report(co2_summary, co2_metrics, aroma_summary, selected_co2, selected_partition)
    print("[done]", RESULTS_DIR, flush=True)


if __name__ == "__main__":
    main()
