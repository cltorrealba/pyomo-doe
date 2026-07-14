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
from nbclient import NotebookClient
from scipy.stats import chi2

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_final_operational_doe_v2 as final
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "estimability_old_vs_lot1"
LOT1_DIR = SCRIPT_DIR / "results" / "lot1_data_preview" / "processed"
LOT1_ETHANOL_WORKBOOK = SCRIPT_DIR / "data" / "Laboratorio 2026" / "DOE_Lote_1" / "DOE_Lote_1_ethanol_manual_entry.xlsx"
LOT1_ETHANOL_WORKBOOK_FALLBACK = (
    SCRIPT_DIR
    / "results"
    / "lot1_data_preview"
    / "manual_entry_templates"
    / "DOE_Lote_1_ethanol_manual_entry.xlsx"
)
NOTEBOOK_DIR = SCRIPT_DIR / "laboratory_2026" / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "fermentation_estimability_old_vs_lot1.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "fermentation_estimability_old_vs_lot1.executed.ipynb"
)

CORE_PARAMETERS = final.FERMENTATION_TARGETS
PROFILE_DEFAULT = ("Kd0", "qN")
SECONDARY_STATES = ("Pyr", "AcAld", "Acetate", "O2")
CHI2_95 = float(chi2.ppf(0.95, df=1))
PERCENT_VV_TO_G_L = 7.89


def safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def to_numeric_manual(value) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() in {"TBD", "PENDING", "NA", "N/A"}:
            return np.nan
        text = text.replace(",", ".")
        return float(pd.to_numeric(text, errors="coerce"))
    return float(pd.to_numeric(value, errors="coerce"))


def is_pending_marker(value) -> bool:
    return isinstance(value, str) and value.strip().upper() in {"TBD", "PENDING"}


def is_nonempty_entry(value) -> bool:
    if pd.isna(value):
        return False
    return bool(str(value).strip())


def load_lot1_ethanol_manual() -> pd.DataFrame:
    workbook = LOT1_ETHANOL_WORKBOOK if LOT1_ETHANOL_WORKBOOK.exists() else LOT1_ETHANOL_WORKBOOK_FALLBACK
    if not workbook.exists():
        return pd.DataFrame()
    frames = []
    for process in ("F1", "F2", "F3"):
        df = pd.read_excel(workbook, sheet_name=process, dtype=object)
        df["source_sheet"] = process
        frames.append(df)
    ethanol = pd.concat(frames, ignore_index=True)
    for col in ethanol.columns:
        if ethanol[col].dtype == object:
            ethanol[col] = ethanol[col].map(lambda x: x.strip() if isinstance(x, str) else x)

    input_cols = ["ethanol_percent_vv_rep1", "ethanol_percent_vv_rep2", "ethanol_g_l_manual"]
    for col in input_cols:
        if col not in ethanol.columns:
            ethanol[col] = np.nan
    ethanol["ethanol_percent_vv_rep1_num"] = ethanol["ethanol_percent_vv_rep1"].map(to_numeric_manual)
    ethanol["ethanol_percent_vv_rep2_num"] = ethanol["ethanol_percent_vv_rep2"].map(to_numeric_manual)
    ethanol["ethanol_g_l_manual_num"] = ethanol["ethanol_g_l_manual"].map(to_numeric_manual)
    ethanol["ethanol_pending"] = ethanol[input_cols].map(is_pending_marker).any(axis=1)
    ethanol["ethanol_has_entry"] = ethanol[input_cols].map(is_nonempty_entry).any(axis=1)

    ethanol["ethanol_percent_vv_final_calc"] = ethanol[
        ["ethanol_percent_vv_rep1_num", "ethanol_percent_vv_rep2_num"]
    ].mean(axis=1, skipna=True)
    empty_replicates = ethanol[["ethanol_percent_vv_rep1_num", "ethanol_percent_vv_rep2_num"]].notna().sum(axis=1).eq(0)
    ethanol.loc[empty_replicates, "ethanol_percent_vv_final_calc"] = np.nan
    ethanol["ethanol_g_l_from_percent_vv"] = ethanol["ethanol_percent_vv_final_calc"] * PERCENT_VV_TO_G_L
    ethanol["ethanol_g_l_for_model_manual"] = ethanol["ethanol_g_l_manual_num"].combine_first(
        ethanol["ethanol_g_l_from_percent_vv"]
    )
    ethanol["ethanol_status"] = np.select(
        [
            ethanol["ethanol_g_l_for_model_manual"].notna(),
            ethanol["ethanol_pending"],
            ethanol["ethanol_has_entry"],
        ],
        ["measured", "pending_TBD", "non_numeric_entry"],
        default="not_measured",
    )
    ethanol["ethanol_sample_type"] = np.where(
        ethanol["ethanol_status"].isin(["measured", "pending_TBD", "non_numeric_entry"]),
        "ethanol_50mL",
        "small_5mL",
    )
    ethanol["sample_volume_ml"] = np.where(ethanol["ethanol_sample_type"].eq("ethanol_50mL"), 50.0, 5.0)

    keep_cols = [
        "process",
        "sample_id",
        "source_sheet",
        "ethanol_percent_vv_rep1",
        "ethanol_percent_vv_rep2",
        "ethanol_percent_vv_rep1_num",
        "ethanol_percent_vv_rep2_num",
        "ethanol_percent_vv_final_calc",
        "ethanol_g_l_manual_num",
        "ethanol_g_l_from_percent_vv",
        "ethanol_g_l_for_model_manual",
        "ethanol_pending",
        "ethanol_status",
        "ethanol_sample_type",
        "sample_volume_ml",
        "ethanol_measurement_datetime",
        "ethanol_method",
        "operator",
        "qc_flag",
        "notes",
    ]
    ethanol = ethanol[[c for c in keep_cols if c in ethanol.columns]].copy()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ethanol.to_csv(RESULTS_DIR / "lot1_ethanol_manual_processed_from_data.csv", index=False)
    ethanol.to_csv(LOT1_DIR / "lot1_ethanol_manual_processed.csv", index=False)
    return ethanol


def apply_lot1_ethanol_manual(y15: pd.DataFrame) -> pd.DataFrame:
    ethanol = load_lot1_ethanol_manual()
    if ethanol.empty:
        return y15
    merge_cols = [c for c in ethanol.columns if c not in {"source_sheet"}]
    ethanol_cols = [c for c in merge_cols if c not in {"process", "sample_id"}]
    y15 = y15.drop(columns=[c for c in ethanol_cols if c in y15.columns], errors="ignore")
    y15 = y15.merge(ethanol[merge_cols], on=["process", "sample_id"], how="left")
    y15["ethanol_status"] = y15["ethanol_status"].fillna("not_measured")
    y15["ethanol_sample_type"] = y15["ethanol_sample_type"].fillna("small_5mL")
    y15["sample_volume_ml"] = pd.to_numeric(y15["sample_volume_ml"], errors="coerce").fillna(5.0)
    y15.to_csv(RESULTS_DIR / "lot1_y15_wide_with_latest_ethanol.csv", index=False)
    return y15


def load_lot1_processed() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y15 = pd.read_csv(LOT1_DIR / "lot1_y15_wide_processed.csv")
    y15 = apply_lot1_ethanol_manual(y15)
    oculyze = pd.read_csv(LOT1_DIR / "lot1_oculyze_processed.csv")
    inputs = pd.read_csv(LOT1_DIR / "lot1_actual_operational_inputs_loaded.csv")
    temp = pd.read_csv(LOT1_DIR / "lot1_temp_10min_preview.csv")
    co2 = pd.read_csv(LOT1_DIR / "lot1_co2_filt_volume_corrected_10min.csv")
    numeric_cols = {
        "t_h",
        "sample_number",
        "sugar_total_g_l",
        "glucose_g_l",
        "fructose_g_l",
        "yan_for_model_physical_mg_l",
        "ethanol_g_l_for_model_manual",
        "glycerol_g_l",
        "pyruvic_acid_mg_l",
        "acetaldehyde_mg_l",
        "acetic_acid_g_l",
        "x_viable_g_l",
        "x_dead_g_l",
        "x_total_g_l",
        "actual_plot_t_h",
        "actual_amount_for_plot",
        "actual_amount",
        "T",
        "SP",
        "flow_filt_sccm",
        "flow_filt_sccm_v0_corrected",
    }
    for df in (y15, oculyze, inputs, temp, co2):
        for col in set(df.columns).intersection(numeric_cols):
            df[col] = safe_num(df[col])
    return y15, oculyze, inputs, temp, co2


def first_finite(values: pd.Series | np.ndarray, default: float = np.nan) -> float:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if arr.empty:
        return float(default)
    return float(arr.iloc[0])


def values_by_time(df: pd.DataFrame, column: str) -> dict[float, float]:
    if column not in df.columns:
        return {}
    out: dict[float, float] = {}
    for row in df[["t_h", column]].dropna(subset=["t_h"]).itertuples(index=False):
        value = float(getattr(row, column))
        if np.isfinite(value):
            out[round(float(row.t_h), 8)] = value
    return out


def infer_lot1_temperature(process: str, t: float, temp: pd.DataFrame, fallback: float = 18.0) -> float:
    sub = temp[temp["process"].eq(process)].sort_values("t_h")
    if sub.empty or "T" not in sub:
        return float(fallback)
    tt = sub["t_h"].to_numpy(dtype=float)
    yy = sub["T"].to_numpy(dtype=float)
    mask = np.isfinite(tt) & np.isfinite(yy)
    if not mask.any():
        return float(fallback)
    return float(np.interp(float(t), tt[mask], yy[mask]))


def lot1_pulses(process: str, inputs: pd.DataFrame) -> dict[str, tuple[tuple[float, float], ...]]:
    pulses = {channel: [] for channel in base.INPUT_CHANNELS}
    sub = inputs[inputs["process"].eq(process)].copy()
    for row in sub.itertuples(index=False):
        executed = str(getattr(row, "executed_bool", getattr(row, "executed", ""))).lower() in {"true", "1", "yes"}
        if not executed:
            continue
        channel = str(getattr(row, "channel", "")).strip().upper()
        if channel not in pulses:
            continue
        time_h = float(getattr(row, "actual_plot_t_h", np.nan))
        amount = float(getattr(row, "actual_amount_for_plot", getattr(row, "actual_amount", np.nan)))
        if np.isfinite(time_h) and np.isfinite(amount):
            pulses[channel].append((time_h, amount))
    return {channel: tuple(sorted(rows)) for channel, rows in pulses.items()}


def make_lot1_batches(include_f2_tail: bool = True) -> tuple[list[base.BatchData], dict[str, np.ndarray]]:
    y15, oculyze, inputs, temp, co2 = load_lot1_processed()
    batches: list[base.BatchData] = []
    co2_times: dict[str, np.ndarray] = {}
    for process in sorted(y15["process"].dropna().unique()):
        y = y15[y15["process"].eq(process)].copy()
        x = oculyze[oculyze["process"].eq(process)].copy()
        if not include_f2_tail and process == "F2":
            y = y[y["t_h"] <= 262.75]
            x = x[x["t_h"] <= 262.75]
        time = sorted(set(y["t_h"].dropna().round(8)).union(set(x["t_h"].dropna().round(8))))
        if len(time) < 3:
            continue
        frame = pd.DataFrame({"t_h": np.asarray(time, dtype=float)})
        y_maps = {
            "N": values_by_time(y, "yan_for_model_physical_mg_l"),
            "G": values_by_time(y, "glucose_g_l"),
            "F": values_by_time(y, "fructose_g_l"),
            "E": values_by_time(y, "ethanol_g_l_for_model_manual"),
            "Gly": values_by_time(y, "glycerol_g_l"),
            "Pyr": values_by_time(y, "pyruvic_acid_mg_l"),
            "AcAld": values_by_time(y, "acetaldehyde_mg_l"),
            "Acetate": values_by_time(y, "acetic_acid_g_l"),
        }
        x_maps = {
            "X": values_by_time(x, "x_viable_g_l"),
            "Xd": values_by_time(x, "x_dead_g_l"),
        }
        observations: dict[str, np.ndarray] = {}
        for state in base.STATE_NAMES:
            source = x_maps[state] if state in x_maps else y_maps.get(state, {})
            values = []
            for t in frame["t_h"]:
                value = source.get(round(float(t), 8), np.nan)
                if state == "N" and np.isfinite(value):
                    value = value / 1000.0
                values.append(value)
            observations[state] = np.asarray(values, dtype=float)
        for state in ("Pyr", "AcAld", "Acetate"):
            values = [y_maps[state].get(round(float(t), 8), np.nan) for t in frame["t_h"]]
            observations[state] = np.asarray(values, dtype=float)
        observations["O2"] = np.full(len(frame), np.nan, dtype=float)

        temperature = np.asarray([infer_lot1_temperature(process, t, temp) for t in frame["t_h"]], dtype=float)
        initials = {state: first_finite(observations[state], base.DEFAULT_THETA.get(state, 0.0)) for state in base.STATE_NAMES}
        initials.update(
            {
                "Pyr": first_finite(observations.get("Pyr", []), 0.0),
                "AcAld": first_finite(observations.get("AcAld", []), 0.0),
                "Acetate": first_finite(observations.get("Acetate", []), 0.0),
                "O2": 6.5,
                "CO2": 0.0,
            }
        )
        batch = base.BatchData(
            medium="synthetic",
            batch=f"lot1_{process}",
            time=frame["t_h"].to_numpy(dtype=float),
            temperature_c=temperature,
            pulses=lot1_pulses(process, inputs),
            observations=observations,
            initials=initials,
        )
        batches.append(batch)

        sub_co2 = co2[co2["process"].eq(process)].copy()
        if not sub_co2.empty:
            max_t = float(sub_co2["t_h"].max())
            co2_times[batch.label] = np.arange(0.0, max_t + 1e-9, 6.0, dtype=float)
    return batches, co2_times


def summarize_batches(batches: list[base.BatchData], label: str) -> pd.DataFrame:
    rows = []
    for b in batches:
        rows.append(
            {
                "case": label,
                "medium": b.medium,
                "batch": b.batch,
                "n_time": len(b.time),
                "t_min_h": float(np.min(b.time)),
                "t_max_h": float(np.max(b.time)),
                "X0": b.initials.get("X", np.nan),
                "Xd0": b.initials.get("Xd", np.nan),
                "N0": b.initials.get("N", np.nan),
                "G0": b.initials.get("G", np.nan),
                "F0": b.initials.get("F", np.nan),
                "E0": b.initials.get("E", np.nan),
                "Gly0": b.initials.get("Gly", np.nan),
                "n_core_obs": int(
                    sum(np.isfinite(np.asarray(b.observations.get(state, []), dtype=float)).sum() for state in base.STATE_NAMES)
                ),
                "n_secondary_obs": int(
                    sum(np.isfinite(np.asarray(b.observations.get(state, []), dtype=float)).sum() for state in SECONDARY_STATES)
                ),
                "pulses": "; ".join(
                    f"{channel}:{','.join(f'{t:g}@{a:g}' for t, a in rows)}"
                    for channel, rows in b.pulses.items()
                    if rows
                ),
            }
        )
    return pd.DataFrame(rows)


def observed_residual_vector(
    theta: dict[str, float],
    batches: list[base.BatchData],
    co2_times: dict[str, np.ndarray] | None = None,
    include_secondary: bool = True,
    include_co2_information: bool = True,
) -> np.ndarray:
    co2_times = co2_times or {}
    residuals: list[np.ndarray] = []
    for batch in batches:
        extra = co2_times.get(batch.label, np.array([], dtype=float)) if include_co2_information else np.array([], dtype=float)
        output_time = np.asarray(sorted(set(batch.time).union(set(float(t) for t in extra))), dtype=float)
        core = base.simulate(batch, theta, output_time)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        for state in base.STATE_NAMES:
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan)), dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = core.loc[batch.time, state].to_numpy(dtype=float)
            residuals.append((pred[mask] - obs[mask]) / base.sigma_for_state(state, obs[mask]))
        if include_secondary:
            sec = v2.integrate_secondary_v2(batch, theta, core)
            for state in ("Pyr", "AcAld", "Acetate", "O2"):
                obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan)), dtype=float)
                mask = np.isfinite(obs)
                if not mask.any():
                    continue
                pred = sec.loc[batch.time, state].to_numpy(dtype=float)
                residuals.append((pred[mask] - obs[mask]) / joint.SIGMA[state])
            if include_co2_information and len(extra):
                times = [float(t) for t in extra if float(t) in sec.index]
                if times:
                    residuals.append(sec.loc[times, "CO2"].to_numpy(dtype=float) / joint.SIGMA["CO2"])
    if not residuals:
        return np.array([], dtype=float)
    return np.concatenate(residuals)


def estimability_from_fim(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...], label: str) -> pd.DataFrame:
    cov = joint.stable_inverse(fim)
    rows = []
    for i, name in enumerate(parameters):
        std = float(math.sqrt(max(cov[i, i], 0.0)))
        low, high = final.parameter_bounds(name)
        value = float(theta[name])
        active = value <= low * 1.01 or value >= high / 1.01
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
                "case": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std,
                "approx_95_multiplier": float(math.exp(1.96 * min(std, 20.0))),
                "fim_diag": float(fim[i, i]),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return pd.DataFrame(rows)


def eigen_tables(fim: np.ndarray, parameters: tuple[str, ...], label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    vals, vecs = np.linalg.eigh(fim)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    max_val = float(np.max(vals)) if len(vals) else np.nan
    eig = pd.DataFrame(
        {
            "case": label,
            "direction": np.arange(1, len(vals) + 1),
            "eigenvalue": vals,
            "relative_eigenvalue": vals / max_val if np.isfinite(max_val) and max_val > 0 else np.nan,
        }
    )
    load_rows = []
    for idx in range(min(8, vecs.shape[1])):
        abs_vec = np.abs(vecs[:, idx])
        top = np.argsort(abs_vec)[::-1][: min(7, len(parameters))]
        load_rows.append(
            {
                "case": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(vals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs_vec[i]:.3f}" for i in top),
            }
        )
    return eig, pd.DataFrame(load_rows)


def perturb_seed(theta: dict[str, float], multipliers: dict[str, float]) -> dict[str, float]:
    seed = dict(theta)
    for name, multiplier in multipliers.items():
        if name not in seed:
            continue
        low, high = final.parameter_bounds(name)
        value = float(seed[name]) * float(multiplier)
        value = min(max(value, low * 1.001), high / 1.001)
        seed[name] = value
    return seed


def fit_core_case(
    label: str,
    batches: list[base.BatchData],
    theta0: dict[str, float],
    max_nfev: int,
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    candidate_rows = []
    theta_default, summary_default = base.fit_parameters(
        f"{label}__seed_default",
        batches,
        theta0,
        CORE_PARAMETERS,
        max_nfev=max_nfev,
    )
    best_theta = theta_default
    best_summary = dict(summary_default)
    best_seed = "seed_default"
    candidate_rows.append({"case": label, "seed": best_seed, **summary_default})

    seed_specs = {
        "qN_up_25pct": {"qN": 1.25},
        "Kd0_up_2x": {"Kd0": 2.0},
        "qN_up_25pct_Kd0_up_2x": {"qN": 1.25, "Kd0": 2.0},
    }
    for seed_name, multipliers in seed_specs.items():
        seed = perturb_seed(theta_default, multipliers)
        theta_candidate, summary_candidate = base.fit_parameters(
            f"{label}__{seed_name}",
            batches,
            seed,
            CORE_PARAMETERS,
            max_nfev=max_nfev,
        )
        candidate_rows.append({"case": label, "seed": seed_name, **summary_candidate})
        if float(summary_candidate["final_wsse"]) < float(best_summary["final_wsse"]):
            best_theta = theta_candidate
            best_summary = dict(summary_candidate)
            best_seed = seed_name

    best_summary["fit"] = label
    best_summary["best_seed"] = best_seed
    return best_theta, pd.DataFrame([best_summary]), pd.DataFrame(candidate_rows)


def log_progress(message: str) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with (RESULTS_DIR / "run_progress.log").open("a", encoding="utf-8") as fh:
        fh.write(message + "\n")
    print(message, flush=True)


def profile_core_case(
    label: str,
    batches: list[base.BatchData],
    theta_hat: dict[str, float],
    parameters: tuple[str, ...],
    grid_points: int,
    max_nfev: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    resid = base.residual_vector(theta_hat, batches)
    base_obj = float(np.dot(resid, resid))
    prof, summary = base.profile_parameters(
        theta_hat,
        batches,
        CORE_PARAMETERS,
        parameters,
        base_objective=base_obj,
        max_nfev=max_nfev,
        grid_points=grid_points,
    )
    if not prof.empty:
        prof.insert(0, "case", label)
    if not summary.empty:
        summary.insert(0, "case", label)
    return prof, summary


def plot_estimability(estimability: pd.DataFrame) -> None:
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    data = estimability.copy()
    order = (
        data.groupby("parameter")["std_log_approx"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    fig, ax = plt.subplots(figsize=(13, 6))
    width = 0.38
    x = np.arange(len(order))
    for j, case in enumerate(["historical_only", "historical_plus_lot1"]):
        sub = data[data["case"].eq(case)].set_index("parameter").reindex(order)
        ax.bar(x + (j - 0.5) * width, sub["std_log_approx"], width=width, label=case)
    ax.axhline(0.35, color="tab:green", linestyle="--", linewidth=1, label="well/moderate threshold")
    ax.axhline(0.75, color="tab:orange", linestyle="--", linewidth=1, label="moderate/weak threshold")
    ax.set_ylabel("Approx. std in log(parameter)")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=70, ha="right")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "estimability_std_log_comparison.png", dpi=170)
    plt.close(fig)


def plot_profiles(profile: pd.DataFrame) -> None:
    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    if profile.empty:
        return
    for case, group_case in profile.groupby("case"):
        params = list(group_case["profiled_parameter"].dropna().unique())
        ncols = 2
        nrows = int(math.ceil(len(params) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4 * nrows), squeeze=False)
        for ax, parameter in zip(axes.ravel(), params):
            group = group_case[group_case["profiled_parameter"].eq(parameter)].sort_values("theta_value")
            ok = group[group["success"]]
            ax.plot(ok["theta_value"], ok["lr_stat"], marker="o")
            if not ok.empty:
                ax.axvline(float(ok["theta_hat"].iloc[0]), color="black", linestyle=":", linewidth=1)
            ax.axhline(CHI2_95, color="tab:red", linestyle="--", linewidth=1, label="95% chi-square")
            ax.set_title(f"{case}: {parameter}")
            ax.set_xlabel("fixed parameter value")
            ax.set_ylabel("profile LR statistic")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=8)
        for ax in axes.ravel()[len(params) :]:
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(fig_dir / f"profile_likelihood_{case}.png", dpi=170)
        plt.close(fig)


def run_analysis(profile_parameters: tuple[str, ...], grid_points: int, fit_nfev: int, profile_nfev: int, step: float) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "run_progress.log").write_text("", encoding="utf-8")
    data_old = base.load_normalized_data()
    old_batches = base.make_batches(data_old)
    lot1_batches, lot1_co2_times = make_lot1_batches()
    combined_batches = old_batches + lot1_batches
    batch_summary = pd.concat(
        [
            summarize_batches(old_batches, "historical_only"),
            summarize_batches(lot1_batches, "lot1_added"),
            summarize_batches(combined_batches, "historical_plus_lot1"),
        ],
        ignore_index=True,
    )
    batch_summary.to_csv(RESULTS_DIR / "batch_summary.csv", index=False)

    theta0 = final.load_theta_final()
    cases = {
        "historical_only": (old_batches, {}),
        "historical_plus_lot1": (combined_batches, lot1_co2_times),
    }
    fit_summaries = []
    fit_candidate_summaries = []
    theta_rows = []
    fim_metrics_rows = []
    estimability_rows = []
    eig_rows = []
    loading_rows = []
    profiles = []
    profile_summaries = []

    for label, (batches, co2_times) in cases.items():
        log_progress(f"[fit] {label}")
        theta_core, fit_summary, fit_candidates = fit_core_case(label, batches, theta0, fit_nfev)
        fit_summaries.append(fit_summary)
        fit_candidate_summaries.append(fit_candidates)
        theta_all = dict(theta0)
        theta_all.update(theta_core)
        for name, value in theta_all.items():
            if name in final.TARGET_PARAMETERS:
                theta_rows.append({"case": label, "parameter": name, "theta": float(value)})

        log_progress(f"[fim-core] {label}")
        jac_core, _ = base.build_jacobian(theta_core, CORE_PARAMETERS, batches, step=step)
        fim_core = jac_core.T @ jac_core
        core_metrics = joint.fim_metrics(fim_core, CORE_PARAMETERS)
        core_metrics.update({"case": label, "block": "core_profile_parameters", "n_parameters": len(CORE_PARAMETERS)})
        fim_metrics_rows.append(core_metrics)

        log_progress(f"[fim-extended] {label}")
        fim_ext = final.finite_difference_fim(
            theta_all,
            final.TARGET_PARAMETERS,
            lambda th, b=batches, c=co2_times: observed_residual_vector(
                th,
                b,
                c,
                include_secondary=True,
                include_co2_information=bool(c),
            ),
            step,
        )
        ext_metrics = joint.fim_metrics(fim_ext, final.TARGET_PARAMETERS)
        ext_metrics.update({"case": label, "block": "extended_observed_plus_co2_information", "n_parameters": len(final.TARGET_PARAMETERS)})
        fim_metrics_rows.append(ext_metrics)
        estimability_rows.append(estimability_from_fim(fim_ext, theta_all, final.TARGET_PARAMETERS, label))
        eig, loads = eigen_tables(fim_ext, final.TARGET_PARAMETERS, label)
        eig_rows.append(eig)
        loading_rows.append(loads)
        pd.DataFrame(fim_ext, index=final.TARGET_PARAMETERS, columns=final.TARGET_PARAMETERS).to_csv(RESULTS_DIR / f"fim_extended_{label}.csv")
        pd.DataFrame(fim_core, index=CORE_PARAMETERS, columns=CORE_PARAMETERS).to_csv(RESULTS_DIR / f"fim_core_{label}.csv")

        log_progress(f"[profile] {label}: {', '.join(profile_parameters)}")
        prof, prof_summary = profile_core_case(
            label,
            batches,
            theta_core,
            profile_parameters,
            grid_points=grid_points,
            max_nfev=profile_nfev,
        )
        profiles.append(prof)
        profile_summaries.append(prof_summary)

    pd.concat(fit_summaries, ignore_index=True).to_csv(RESULTS_DIR / "core_fit_summary.csv", index=False)
    pd.concat(fit_candidate_summaries, ignore_index=True).to_csv(RESULTS_DIR / "core_fit_candidate_summary.csv", index=False)
    pd.DataFrame(theta_rows).to_csv(RESULTS_DIR / "theta_by_case.csv", index=False)
    pd.DataFrame(fim_metrics_rows).to_csv(RESULTS_DIR / "fim_metrics_by_case.csv", index=False)
    estimability = pd.concat(estimability_rows, ignore_index=True)
    estimability.to_csv(RESULTS_DIR / "estimability_extended_by_case.csv", index=False)
    pd.concat(eig_rows, ignore_index=True).to_csv(RESULTS_DIR / "eigenvalues_extended_by_case.csv", index=False)
    pd.concat(loading_rows, ignore_index=True).to_csv(RESULTS_DIR / "weak_directions_extended_by_case.csv", index=False)
    profile = pd.concat(profiles, ignore_index=True) if profiles else pd.DataFrame()
    profile.to_csv(RESULTS_DIR / "profile_likelihood_core_by_case.csv", index=False)
    profile_summary = pd.concat(profile_summaries, ignore_index=True) if profile_summaries else pd.DataFrame()
    profile_summary.to_csv(RESULTS_DIR / "profile_likelihood_summary_core_by_case.csv", index=False)
    plot_estimability(estimability)
    plot_profiles(profile)

    comparison = estimability.pivot_table(index="parameter", columns="case", values="std_log_approx", aggfunc="first").reset_index()
    if {"historical_only", "historical_plus_lot1"}.issubset(comparison.columns):
        comparison["std_log_delta_plus_minus_old"] = comparison["historical_plus_lot1"] - comparison["historical_only"]
        comparison["std_log_ratio_plus_over_old"] = comparison["historical_plus_lot1"] / comparison["historical_only"].replace(0.0, np.nan)
    comparison.to_csv(RESULTS_DIR / "estimability_change_old_vs_plus_lot1.csv", index=False)

    summary = {
        "profile_parameters": profile_parameters,
        "grid_points": int(grid_points),
        "fit_nfev": int(fit_nfev),
        "profile_nfev": int(profile_nfev),
        "step": float(step),
        "n_historical_batches": len(old_batches),
        "n_lot1_batches": len(lot1_batches),
        "n_combined_batches": len(combined_batches),
        "note": "Extended FIM includes observed core/secondary residuals. Lot 1 online CO2 is included as an information-time block, not as a calibrated flow residual.",
    }
    (RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_markdown_report()
    log_progress("[done] estimability/profile comparison complete")


def write_markdown_report() -> None:
    metrics = pd.read_csv(RESULTS_DIR / "fim_metrics_by_case.csv")
    fit_candidates = pd.read_csv(RESULTS_DIR / "core_fit_candidate_summary.csv")
    estimability = pd.read_csv(RESULTS_DIR / "estimability_change_old_vs_plus_lot1.csv")
    profile_summary = pd.read_csv(RESULTS_DIR / "profile_likelihood_summary_core_by_case.csv")
    weak = pd.read_csv(RESULTS_DIR / "weak_directions_extended_by_case.csv")
    top_change = estimability.sort_values("std_log_delta_plus_minus_old").head(12)
    weak_plus = weak[weak["case"].eq("historical_plus_lot1")].head(6)

    def as_markdown(df: pd.DataFrame) -> str:
        try:
            return df.to_markdown(index=False)
        except Exception:
            return "```\n" + df.to_string(index=False) + "\n```"

    report = f"""# Estimability: historical data vs historical + Lot 1

## What was compared

Two information sets were evaluated with the same model and parameter set:

1. `historical_only`: the normalized historical natural/synthetic must database.
2. `historical_plus_lot1`: historical data plus the executed Lot 1 batches reconstructed from `fermentation_lot1_data_preview.executed.ipynb`.

The extended FIM uses observed core states and secondary states where available. Lot 1 online CO2 is included as an information-time block because the current model state is cumulative CO2, whereas the sensor reports gas flow; this should be interpreted as information availability, not as a final calibrated CO2 likelihood.

## FIM Metrics

{as_markdown(metrics)}

## Core Fit Candidate Audit

{as_markdown(fit_candidates[["case", "seed", "final_wsse", "nfev", "success"]])}

## Largest Approximate Estimability Improvements

Negative `std_log_delta_plus_minus_old` means the parameter became more estimable after adding Lot 1.

{as_markdown(top_change)}

## Profile Likelihood Summary

{as_markdown(profile_summary)}

## Weak Directions After Adding Lot 1

{as_markdown(weak_plus)}

## Generated figures

- `figures/estimability_std_log_comparison.png`
- `figures/profile_likelihood_historical_only.png`
- `figures/profile_likelihood_historical_plus_lot1.png`
"""
    (RESULTS_DIR / "estimability_old_vs_lot1_report.md").write_text(report, encoding="utf-8")


def create_notebook() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            """# Estimability and Profile Likelihood: Historical Data vs Historical + Lot 1

This notebook compares parameter estimability before and after adding the executed Lot 1 data.

The workflow follows the same model-based design logic used in the campaign:

1. Load historical normalized natural/synthetic must data.
2. Reconstruct Lot 1 as executable model batches using actual initial states, sampling times and pulse logs.
3. Fit the core kinetic parameter block for both information sets.
4. Compute FIM, covariance approximations, eigenvalue diagnostics and weak directions.
5. Run profile likelihood for selected weak core kinetic parameters.

The extended FIM includes observed core and secondary states. Online CO2 is included as an information-time block, not as a final flow-calibrated likelihood, because the current model state is cumulative CO2 while the sensor reports gas flow.

The executed notebook uses a screening profile likelihood for `Kd0` and `qN` to keep the old-vs-new comparison operational. To make a more exhaustive profile, rerun the analysis cell with more parameters, for example `("Kd0", "qN", "qEG", "betaG0")`, and increase `grid_points` and `profile_nfev`.

The core fit uses a small deterministic polish/multistart around `qN` and `Kd0`. This is included because the profile likelihood can otherwise reveal a lower objective than the initial least-squares solution.
"""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import pandas as pd
from IPython.display import display, Image

import run_estimability_old_vs_lot1 as analysis

RESULTS = Path("results") / "estimability_old_vs_lot1"
"""
        ),
        nbformat.v4.new_markdown_cell("## Run Analysis"),
        nbformat.v4.new_code_cell(
            """analysis.run_analysis(
    profile_parameters=("Kd0", "qN"),
    grid_points=3,
    fit_nfev=300,
    profile_nfev=45,
    step=0.04,
)"""
        ),
        nbformat.v4.new_markdown_cell("## Batch Inventory"),
        nbformat.v4.new_code_cell(
            """batch_summary = pd.read_csv(RESULTS / "batch_summary.csv")
display(batch_summary)"""
        ),
        nbformat.v4.new_markdown_cell("## FIM Metrics"),
        nbformat.v4.new_code_cell(
            """fit_candidates = pd.read_csv(RESULTS / "core_fit_candidate_summary.csv")
display(fit_candidates[["case", "seed", "final_wsse", "nfev", "success"]])

fim_metrics = pd.read_csv(RESULTS / "fim_metrics_by_case.csv")
display(fim_metrics)"""
        ),
        nbformat.v4.new_markdown_cell("## Extended Estimability Comparison"),
        nbformat.v4.new_code_cell(
            """estimability = pd.read_csv(RESULTS / "estimability_extended_by_case.csv")
change = pd.read_csv(RESULTS / "estimability_change_old_vs_plus_lot1.csv")
display(estimability.sort_values(["case", "std_log_approx"]))
display(change.sort_values("std_log_delta_plus_minus_old"))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=str(RESULTS / "figures" / "estimability_std_log_comparison.png")))"""
        ),
        nbformat.v4.new_markdown_cell("## Weak Eigen Directions"),
        nbformat.v4.new_code_cell(
            """weak = pd.read_csv(RESULTS / "weak_directions_extended_by_case.csv")
display(weak)"""
        ),
        nbformat.v4.new_markdown_cell("## Profile Likelihood"),
        nbformat.v4.new_code_cell(
            """profile_summary = pd.read_csv(RESULTS / "profile_likelihood_summary_core_by_case.csv")
profile = pd.read_csv(RESULTS / "profile_likelihood_core_by_case.csv")
display(profile_summary)
display(profile.head(30))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=str(RESULTS / "figures" / "profile_likelihood_historical_only.png")))
display(Image(filename=str(RESULTS / "figures" / "profile_likelihood_historical_plus_lot1.png")))"""
        ),
        nbformat.v4.new_markdown_cell("## Report"),
        nbformat.v4.new_code_cell(
            """report = (RESULTS / "estimability_old_vs_lot1_report.md").read_text(encoding="utf-8")
print(report)"""
        ),
    ]
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.metadata["language_info"] = {"name": "python", "version": "3.12"}
    nbformat.write(nb, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 7200) -> None:
    create_notebook()
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    client = NotebookClient(nb, timeout=timeout, kernel_name="python3", resources={"metadata": {"path": str(SCRIPT_DIR)}})
    client.execute()
    nbformat.write(nb, EXECUTED_NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare estimability for historical data vs historical + executed Lot 1.")
    parser.add_argument("--create-notebook", action="store_true")
    parser.add_argument("--execute-notebook", action="store_true")
    parser.add_argument("--skip-analysis", action="store_true")
    parser.add_argument("--grid-points", type=int, default=3)
    parser.add_argument("--fit-nfev", type=int, default=300)
    parser.add_argument("--profile-nfev", type=int, default=45)
    parser.add_argument("--step", type=float, default=0.04)
    parser.add_argument("--profile-parameters", nargs="*", default=list(PROFILE_DEFAULT))
    args = parser.parse_args()

    if args.execute_notebook:
        execute_notebook()
        return
    if args.create_notebook:
        create_notebook()
    if not args.skip_analysis:
        run_analysis(
            profile_parameters=tuple(args.profile_parameters),
            grid_points=args.grid_points,
            fit_nfev=args.fit_nfev,
            profile_nfev=args.profile_nfev,
            step=args.step,
        )


if __name__ == "__main__":
    main()
