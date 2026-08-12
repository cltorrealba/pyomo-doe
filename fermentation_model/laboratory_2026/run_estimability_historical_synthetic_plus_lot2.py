from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("JUPYTER_ALLOW_INSECURE_WRITES", "1")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from nbclient import NotebookClient


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SCRIPT_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from laboratory_2026 import run_estimability_historical_by_medium as historical
from laboratory_2026 import run_estimability_old_vs_lot1 as reference
from laboratory_2026 import run_final_operational_doe_v2 as final
from laboratory_2026 import run_lot2_data_preview as lot2_preview
from shared import run_new_must_glycerol_estimability_doe as base
from shared import run_secondary_joint_campaign_doe as joint
from shared import run_secondary_v2_model_evaluation as secondary_v2
from shared.paths import LABORATORY_2026_RESULTS_DIR


RESULTS_DIR = LABORATORY_2026_RESULTS_DIR / "estimability_historical_synthetic_plus_lot2"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "fermentation_estimability_synthetic_prior_plus_lot2.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "fermentation_estimability_synthetic_prior_plus_lot2.executed.ipynb"
)
HISTORICAL_RESULTS_DIR = LABORATORY_2026_RESULTS_DIR / "estimability_historical_synthetic"

CORE_PARAMETERS = final.FERMENTATION_TARGETS
SECONDARY_PARAMETERS = final.SECONDARY_TARGETS
EXTENDED_PARAMETERS = final.TARGET_PARAMETERS
PROFILE_DEFAULT = ("Kd0", "qN")
SECONDARY_STATES = ("Pyr", "AcAld", "Acetate", "O2")
FINITE_DIFFERENCE_STEP = 0.04
TEMPERATURE_GRID_H = 2.0
CO2_INFORMATION_GRID_H = 6.0
REPORT_PATH = RESULTS_DIR / "estimability_lot2_report.md"


def source_fingerprint() -> str:
    paths = [
        lot2_preview.ETHANOL_WORKBOOK,
        lot2_preview.SCHEDULE_PATH,
        Path(__file__),
        Path(lot2_preview.__file__),
    ]
    paths.extend(sorted(lot2_preview.Y15_DIR.glob("*.txt")))
    for config in lot2_preview.PROCESS_CONFIG.values():
        paths.extend(sorted((lot2_preview.LOT2_DIR / config["folder"]).glob("*.csv")))
    digest = hashlib.sha256()
    for path in sorted(set(paths), key=lambda value: str(value).lower()):
        digest.update(str(path.relative_to(FERMENTATION_MODEL_DIR)).encode("utf-8", errors="replace"))
        try:
            digest.update(path.read_bytes())
        except (FileNotFoundError, OSError):
            digest.update(b"<unreadable>")
    return digest.hexdigest()


def completed_outputs_available() -> bool:
    required = (
        "summary.json",
        "fit_summary.csv",
        "secondary_fit_summary.csv",
        "theta_by_case.csv",
        "fim_metrics.csv",
        "estimability_change_historical_vs_plus_lot2.csv",
        "profile_likelihood_summary_combined.csv",
        "weak_directions_combined.csv",
        "lot2_fit_metrics_by_state.csv",
    )
    return all((RESULTS_DIR / name).exists() for name in required)


def classify_std(std_log: float, active_bound: bool) -> str:
    if std_log <= 0.35 and not active_bound:
        return "well_estimated"
    if std_log <= 0.75 and not active_bound:
        return "moderate"
    if std_log <= 1.25 and not active_bound:
        return "weak_but_actionable"
    return "weak_or_confounded"


def augment_derived_diagnostics() -> None:
    comparison_path = RESULTS_DIR / "estimability_change_historical_vs_plus_lot2.csv"
    fit_path = RESULTS_DIR / "fit_summary.csv"
    secondary_path = RESULTS_DIR / "secondary_fit_summary.csv"
    if comparison_path.exists() and fit_path.exists() and secondary_path.exists():
        comparison = pd.read_csv(comparison_path)
        core_fit = pd.read_csv(fit_path)
        secondary_fit = pd.read_csv(secondary_path)
        core_row = core_fit[core_fit["fit"].eq("historical_synthetic_plus_lot2")]
        secondary_row = secondary_fit[
            secondary_fit["case"].eq("historical_synthetic_plus_lot2")
        ]
        core_scale = math.sqrt(float(core_row["wsse_per_dof"].iloc[0]))
        secondary_scale = math.sqrt(float(secondary_row["wsse_per_residual"].iloc[0]))
        comparison["dispersion_scale_combined"] = np.where(
            comparison["parameter"].isin(CORE_PARAMETERS),
            core_scale,
            np.where(
                comparison["parameter"].isin(SECONDARY_PARAMETERS), secondary_scale, np.nan
            ),
        )
        comparison["std_log_dispersion_adjusted"] = (
            comparison["std_log_combined"] * comparison["dispersion_scale_combined"]
        )
        comparison["approx_95_multiplier_dispersion_adjusted"] = np.exp(
            1.96 * comparison["std_log_dispersion_adjusted"].clip(upper=20.0)
        )
        comparison["classification_dispersion_adjusted"] = comparison.apply(
            lambda row: classify_std(
                float(row["std_log_dispersion_adjusted"]),
                str(row["active_bound_combined"]).lower() == "true",
            )
            if np.isfinite(row["std_log_dispersion_adjusted"])
            else "unsupported_by_current_measurements",
            axis=1,
        )
        comparison.to_csv(comparison_path, index=False)

    fim_path = RESULTS_DIR / "fim_combined_observed.csv"
    metrics_path = RESULTS_DIR / "fim_metrics.csv"
    if fim_path.exists() and metrics_path.exists():
        fim = pd.read_csv(fim_path, index_col=0).loc[
            list(CORE_PARAMETERS + SECONDARY_PARAMETERS),
            list(CORE_PARAMETERS + SECONDARY_PARAMETERS),
        ].to_numpy(dtype=float)
        supported_row = metric_row(
            "historical_synthetic_plus_lot2",
            "core_secondary_observed",
            fim,
            CORE_PARAMETERS + SECONDARY_PARAMETERS,
        )
        metrics = pd.read_csv(metrics_path)
        metrics = metrics[
            ~(
                metrics["case"].eq("historical_synthetic_plus_lot2")
                & metrics["block"].eq("core_secondary_observed")
            )
        ]
        metrics = pd.concat([metrics, pd.DataFrame([supported_row])], ignore_index=True)
        metrics.to_csv(metrics_path, index=False)


def log_progress(message: str) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with (RESULTS_DIR / "run_progress.log").open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")
    print(message, flush=True)


def first_finite(values: pd.Series | np.ndarray, default: float) -> float:
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    return float(numeric.iloc[0]) if not numeric.empty else float(default)


def remap_values(
    sample_time: np.ndarray,
    values: np.ndarray,
    model_time: np.ndarray,
) -> np.ndarray:
    positions = {round(float(value), 8): idx for idx, value in enumerate(model_time)}
    out = np.full(len(model_time), np.nan, dtype=float)
    for time_h, value in zip(sample_time, values):
        if np.isfinite(time_h):
            out[positions[round(float(time_h), 8)]] = value
    return out


def make_lot2_batches() -> tuple[list[base.BatchData], dict[str, np.ndarray], pd.DataFrame]:
    lot2_preview.run_pipeline()
    samples = pd.read_csv(lot2_preview.PROCESSED_DIR / "lot2_samples_integrated.csv")
    temperature = pd.read_csv(lot2_preview.PROCESSED_DIR / "lot2_temperature_10min.csv")
    inputs = pd.read_csv(lot2_preview.PROCESSED_DIR / "lot2_operational_inputs.csv")
    co2 = pd.read_csv(lot2_preview.PROCESSED_DIR / "lot2_co2_volume_corrected_10min.csv")
    numeric_columns = {
        "t_h",
        "glucose_for_model_g_l",
        "fructose_for_model_g_l",
        "yan_for_model_physical_mg_l",
        "ethanol_g_l_for_model",
        "glycerol_g_l",
        "x_viable_g_l",
        "x_dead_g_l",
        "pyruvic_acid_mg_l",
        "acetaldehyde_mg_l",
        "acetic_acid_g_l",
        "actual_or_model_time_h",
        "amount",
        "T",
        "SP",
    }
    for frame in (samples, temperature, inputs, co2):
        for column in numeric_columns.intersection(frame.columns):
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

    batches: list[base.BatchData] = []
    co2_times: dict[str, np.ndarray] = {}
    support_rows: list[dict[str, object]] = []
    for process in lot2_preview.PROCESS_CONFIG:
        sample = samples[samples["process"].eq(process)].dropna(subset=["t_h"]).copy()
        sample = sample.sort_values(["t_h", "sample_sequence"]).drop_duplicates("t_h", keep="last")
        sensor = temperature[temperature["process"].eq(process)].dropna(subset=["t_h", "T"]).sort_values("t_h")
        if sample.empty or sensor.empty:
            continue
        horizon = float(sample["t_h"].max())
        sample_time = sample["t_h"].to_numpy(dtype=float)
        model_grid = np.arange(0.0, horizon + 1e-9, TEMPERATURE_GRID_H)
        if "SP" in sensor:
            setpoint_change = sensor.loc[sensor["SP"].diff().abs().gt(0.25), "t_h"].to_numpy(dtype=float)
        else:
            setpoint_change = np.array([], dtype=float)
        pulse_time = inputs.loc[inputs["process"].eq(process), "actual_or_model_time_h"].dropna().to_numpy(dtype=float)
        model_time = np.asarray(
            sorted(
                set(float(value) for value in sample_time)
                .union(float(value) for value in model_grid)
                .union(float(value) for value in setpoint_change)
                .union(float(value) for value in pulse_time)
            ),
            dtype=float,
        )
        sensor_t = sensor["t_h"].to_numpy(dtype=float)
        sensor_y = sensor["T"].to_numpy(dtype=float)
        model_temperature = np.interp(model_time, sensor_t, sensor_y)
        initial_oculyze_temperature = first_finite(
            sample.get("temperature_oculyze_c", pd.Series(dtype=float)), sensor_y[0]
        )
        model_temperature = np.where(
            model_time < sensor_t.min(), initial_oculyze_temperature, model_temperature
        )

        ethanol = sample["ethanol_g_l_for_model"].to_numpy(dtype=float)
        ethanol_use = sample["ethanol_use_for_model"].astype(str).str.lower().eq("true").to_numpy()
        ethanol = np.where(ethanol_use, ethanol, np.nan)
        state_values = {
            "X": sample["x_viable_g_l"].to_numpy(dtype=float),
            "Xd": sample["x_dead_g_l"].to_numpy(dtype=float),
            "N": sample["yan_for_model_physical_mg_l"].to_numpy(dtype=float) / 1000.0,
            "G": sample["glucose_for_model_g_l"].to_numpy(dtype=float),
            "F": sample["fructose_for_model_g_l"].to_numpy(dtype=float),
            "E": ethanol,
            "Gly": sample["glycerol_g_l"].clip(lower=0.0).to_numpy(dtype=float),
            "Pyr": sample["pyruvic_acid_mg_l"].clip(lower=0.0).to_numpy(dtype=float),
            "AcAld": sample["acetaldehyde_mg_l"].clip(lower=0.0).to_numpy(dtype=float),
            "Acetate": sample["acetic_acid_g_l"].clip(lower=0.0).to_numpy(dtype=float),
            "O2": np.full(len(sample), np.nan, dtype=float),
        }
        observations = {
            state: remap_values(sample_time, values, model_time)
            for state, values in state_values.items()
        }
        defaults = {"X": 0.3, "Xd": 0.02, "N": 0.2, "G": 90.0, "F": 90.0, "E": 0.0, "Gly": 0.0}
        initials = {
            state: first_finite(state_values[state], default)
            for state, default in defaults.items()
        }
        initials.update(
            {
                "Pyr": first_finite(state_values["Pyr"], 0.0),
                "AcAld": first_finite(state_values["AcAld"], 0.0),
                "Acetate": first_finite(state_values["Acetate"], 0.0),
                "O2": 6.5,
                "CO2": 0.0,
            }
        )
        pulses = {channel: [] for channel in base.INPUT_CHANNELS}
        for row in inputs[inputs["process"].eq(process)].itertuples(index=False):
            channel = str(row.channel).strip().upper()
            if channel in pulses and np.isfinite(row.actual_or_model_time_h) and np.isfinite(row.amount):
                pulses[channel].append((float(row.actual_or_model_time_h), float(row.amount)))
        pulses = {channel: tuple(sorted(values)) for channel, values in pulses.items()}
        batch = base.BatchData(
            medium="synthetic",
            batch=f"lot2_{process}",
            time=model_time,
            temperature_c=model_temperature,
            pulses=pulses,
            observations=observations,
            initials=initials,
        )
        batches.append(batch)

        gas = co2[co2["process"].eq(process)].dropna(subset=["t_h"])
        if not gas.empty:
            targets = np.arange(0.0, float(gas["t_h"].max()) + 1e-9, CO2_INFORMATION_GRID_H)
            available = gas["t_h"].to_numpy(dtype=float)
            selected = []
            for target in targets:
                nearest = float(available[np.argmin(np.abs(available - target))])
                if abs(nearest - target) <= 0.5:
                    selected.append(nearest)
            co2_times[batch.label] = np.asarray(sorted(set(selected)), dtype=float)

        for state, values in observations.items():
            support_rows.append(
                {
                    "process": process,
                    "state": state,
                    "n_observations": int(np.isfinite(values).sum()),
                    "used_in_core_fit": state in base.STATE_NAMES,
                    "used_in_secondary_fit": state in {"Pyr", "AcAld", "Acetate", "O2"},
                }
            )
    return batches, co2_times, pd.DataFrame(support_rows)


def load_historical_core_theta(theta0: dict[str, float]) -> dict[str, float]:
    path = HISTORICAL_RESULTS_DIR / "theta.csv"
    if not path.exists():
        return dict(theta0)
    table = pd.read_csv(path)
    out = dict(theta0)
    for row in table.itertuples(index=False):
        if row.parameter in out and np.isfinite(float(row.theta)):
            out[str(row.parameter)] = float(row.theta)
    return out


def fit_secondary_block(
    label: str,
    batches: list[base.BatchData],
    theta: dict[str, float],
    max_nfev: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    core_cache = {batch.label: base.simulate(batch, theta, batch.time) for batch in batches}
    if any(value is None for value in core_cache.values()):
        raise RuntimeError(f"Core simulation failed before secondary fit for {label}")
    fitted, summary = secondary_v2.fit_v2(
        theta,
        batches,
        core_cache,
        max_nfev=max_nfev,
        fit_parameters=SECONDARY_PARAMETERS,
        model_label=label,
    )
    out = dict(theta)
    out.update({name: fitted[name] for name in SECONDARY_PARAMETERS})
    summary.insert(0, "case", label)
    return out, summary


def extended_fim(
    theta: dict[str, float],
    batches: list[base.BatchData],
    co2_times: dict[str, np.ndarray],
    include_co2: bool,
    step: float,
) -> np.ndarray:
    return final.finite_difference_fim(
        theta,
        EXTENDED_PARAMETERS,
        lambda candidate: reference.observed_residual_vector(
            candidate,
            batches,
            co2_times,
            include_secondary=True,
            include_co2_information=include_co2,
        ),
        step,
    )


def metric_row(
    case: str,
    block: str,
    fim: np.ndarray,
    parameters: tuple[str, ...],
) -> dict[str, object]:
    row = joint.fim_metrics(fim, parameters)
    row.update({"case": case, "block": block, "n_parameters": len(parameters)})
    return row


def theta_rows(case: str, theta: dict[str, float]) -> list[dict[str, object]]:
    return [
        {"case": case, "parameter": name, "theta": float(theta[name])}
        for name in EXTENDED_PARAMETERS
    ]


def load_theta_checkpoint() -> tuple[dict[str, float], dict[str, float]] | None:
    path = RESULTS_DIR / "theta_by_case.csv"
    if not path.exists():
        return None
    table = pd.read_csv(path)
    required_cases = {"historical_synthetic", "historical_synthetic_plus_lot2"}
    if not required_cases.issubset(set(table["case"])):
        return None
    base_theta = final.load_theta_final()
    output = []
    for case in ("historical_synthetic", "historical_synthetic_plus_lot2"):
        values = table[table["case"].eq(case)].set_index("parameter")["theta"].to_dict()
        if not set(EXTENDED_PARAMETERS).issubset(values):
            return None
        theta = dict(base_theta)
        theta.update({name: float(values[name]) for name in EXTENDED_PARAMETERS})
        output.append(theta)
    return output[0], output[1]


def batch_fit_metrics(
    case: str,
    batches: list[base.BatchData],
    theta: dict[str, float],
) -> pd.DataFrame:
    rows = []
    for batch in batches:
        core = base.simulate(batch, theta, batch.time)
        if core is None:
            continue
        secondary = secondary_v2.integrate_secondary_v2(batch, theta, core)
        for state in base.STATE_NAMES:
            observed = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            mask = np.isfinite(observed)
            if not mask.any():
                continue
            predicted = core.loc[batch.time, state].to_numpy(dtype=float)[mask]
            error = predicted - observed[mask]
            rows.append(
                {
                    "case": case,
                    "batch": batch.batch,
                    "state": state,
                    "block": "core",
                    "n": int(mask.sum()),
                    "rmse": float(np.sqrt(np.mean(error**2))),
                    "mae": float(np.mean(np.abs(error))),
                }
            )
        for state in ("Pyr", "AcAld", "Acetate"):
            observed = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            mask = np.isfinite(observed)
            if not mask.any():
                continue
            predicted = secondary.loc[batch.time, state].to_numpy(dtype=float)[mask]
            error = predicted - observed[mask]
            rows.append(
                {
                    "case": case,
                    "batch": batch.batch,
                    "state": state,
                    "block": "secondary",
                    "n": int(mask.sum()),
                    "rmse": float(np.sqrt(np.mean(error**2))),
                    "mae": float(np.mean(np.abs(error))),
                }
            )
    return pd.DataFrame(rows)


def plot_lot2_fits(batches: list[base.BatchData], theta: dict[str, float]) -> None:
    fit_dir = FIGURE_DIR / "lot2_fit_overlays"
    fit_dir.mkdir(parents=True, exist_ok=True)
    units = {
        "X": "g/L",
        "Xd": "g/L",
        "N": "kg/m3",
        "G": "g/L",
        "F": "g/L",
        "E": "g/L",
        "Gly": "g/L",
        "Pyr": "mg/L",
        "AcAld": "mg/L",
        "Acetate": "g/L",
    }
    for batch in batches:
        core = base.simulate(batch, theta, batch.time)
        if core is None:
            continue
        secondary = secondary_v2.integrate_secondary_v2(batch, theta, core)
        states = tuple(base.STATE_NAMES) + ("Pyr", "AcAld", "Acetate")
        fig, axes = plt.subplots(5, 2, figsize=(13, 17), squeeze=False)
        for ax, state in zip(axes.ravel(), states):
            prediction = core[state] if state in base.STATE_NAMES else secondary[state]
            observed = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            mask = np.isfinite(observed)
            ax.plot(batch.time, prediction, color="#0072B2", linewidth=1.5, label="model")
            ax.scatter(batch.time[mask], observed[mask], color="#D55E00", s=22, label="data", zorder=3)
            for channel, events in batch.pulses.items():
                if channel == state:
                    for time_h, amount in events:
                        ax.axvline(time_h, color="black", linestyle="--", linewidth=0.8)
                        ax.annotate(f"+{amount:g}", (time_h, 0.96), xycoords=("data", "axes fraction"), rotation=90, va="top", fontsize=7)
            ax.set_title(state)
            ax.set_ylabel(units[state])
            ax.set_xlabel("Elapsed time [h]")
            ax.grid(True, alpha=0.2)
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Combined-prior fit: {batch.batch}")
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(fit_dir / f"fit_{batch.batch}.png", dpi=160)
        plt.close(fig)


def plot_estimability(comparison: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    order = comparison.sort_values("std_log_combined", ascending=False)["parameter"].tolist()
    ordered = comparison.set_index("parameter").reindex(order)
    x = np.arange(len(order))
    width = 0.38
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(x - width / 2, ordered["std_log_historical"], width, label="historical synthetic")
    ax.bar(x + width / 2, ordered["std_log_combined"], width, label="historical + Lot 2")
    ax.axhline(0.35, color="#009E73", linestyle="--", linewidth=1)
    ax.axhline(0.75, color="#E69F00", linestyle="--", linewidth=1)
    ax.set_yscale("log")
    ax.set_ylabel("Approximate std in log(parameter)")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=70, ha="right")
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "estimability_historical_vs_plus_lot2.png", dpi=170)
    plt.close(fig)


def write_report() -> None:
    fit = pd.read_csv(RESULTS_DIR / "fit_summary.csv")
    secondary_fit = pd.read_csv(RESULTS_DIR / "secondary_fit_summary.csv")
    metrics = pd.read_csv(RESULTS_DIR / "fim_metrics.csv")
    comparison = pd.read_csv(RESULTS_DIR / "estimability_change_historical_vs_plus_lot2.csv")
    profiles = pd.read_csv(RESULTS_DIR / "profile_likelihood_summary_combined.csv")
    weak = pd.read_csv(RESULTS_DIR / "weak_directions_combined.csv")
    support = pd.read_csv(RESULTS_DIR / "lot2_measurement_support.csv")
    best = comparison.sort_values("std_log_ratio_combined_over_historical").head(12)
    unresolved = comparison.sort_values("std_log_combined", ascending=False).head(12)
    supported = comparison[comparison["parameter"].isin(CORE_PARAMETERS + SECONDARY_PARAMETERS)].copy()
    supported = supported.sort_values("std_log_dispersion_adjusted")
    full_observed = metrics[metrics["case"].eq("combined_observed_only_at_combined_theta")].iloc[0]
    supported_metrics = metrics[metrics["block"].eq("core_secondary_observed")].iloc[0]
    co2_potential = metrics[
        metrics["case"].eq("combined_observed_plus_co2_schedule_at_combined_theta")
    ].iloc[0]
    trace_reduction_percent = 100.0 * (
        float(full_observed["trace_inv"]) - float(co2_potential["trace_inv"])
    ) / float(full_observed["trace_inv"])
    report = f"""# Historical synthetic prior plus MBDoE Lot 2

## Scope

The historical synthetic-must calibration is the prior dataset. The three executed Lot 2 MBDoE fermentations are then added with their reconstructed actual temperature profiles, observed sample times, model-safe chemistry, Oculyze biomass, and supported pulse schedule.

The 11-parameter core block and the reduced six-parameter secondary block are re-estimated. Aroma parameters are carried in the extended FIM but cannot gain observed-data information because Lot 2 does not contain aroma measurements.

## Core fit

{fit.to_markdown(index=False)}

## Secondary fit

{secondary_fit.to_markdown(index=False)}

## Lot 2 measurement support

{support.to_markdown(index=False)}

## FIM metrics

{metrics.to_markdown(index=False)}

The observed core-plus-secondary block is locally full rank: `{int(supported_metrics['rank_1e-8'])}/{int(supported_metrics['n_parameters'])}`. The full matrix is rank `{int(full_observed['rank_1e-8'])}/{int(full_observed['n_parameters'])}` because the nine aroma columns are exactly unsupported. Adding CO2 measurement times without fitting a gas-flow observation equation leaves the rank unchanged and reduces `trace_inv` by only `{trace_reduction_percent:.3f}%`.

## Largest uncertainty reductions

{best[["parameter", "std_log_historical", "std_log_combined", "std_log_ratio_combined_over_historical", "classification_combined"]].to_markdown(index=False)}

## Remaining weak directions

{unresolved[["parameter", "std_log_combined", "approx_95_multiplier_combined", "active_bound_combined", "classification_combined"]].to_markdown(index=False)}

## Dispersion-adjusted screening

The unscaled FIM assumes that the residual noise model is correct. It is not: core WSSE/DOF is `{float(fit.loc[fit['fit'].eq('historical_synthetic_plus_lot2'), 'wsse_per_dof'].iloc[0]):.3f}` and secondary WSSE/residual is `{float(secondary_fit.loc[secondary_fit['case'].eq('historical_synthetic_plus_lot2'), 'wsse_per_residual'].iloc[0]):.3f}`. The table below multiplies local standard errors by the corresponding square-root dispersion. This is the safer practical-identifiability classification until model discrepancy and analytical error are separated.

{supported[["parameter", "std_log_combined", "dispersion_scale_combined", "std_log_dispersion_adjusted", "approx_95_multiplier_dispersion_adjusted", "classification_dispersion_adjusted"]].to_markdown(index=False)}

## Profile likelihood

{profiles.to_markdown(index=False) if not profiles.empty else 'No profile result was produced.'}

## Weak eigen-directions after Lot 2

{weak.head(8).to_markdown(index=False)}

## Interpretation

- `observed_only` is the defensible practical-identifiability result for the measurements currently linked to the model.
- `observed_plus_CO2_schedule` is a conditional sensitivity calculation. CO2 flow values are not fitted as cumulative CO2 until a gas-flow observation equation is calibrated.
- The FIM is also decomposed at the combined parameter point. Therefore, `historical_component_at_combined_theta + lot2_increment_at_combined_theta` reproduces the combined FIM and isolates the information supplied by Lot 2 from changes in the optimum.
- Nitrogen pulse timing in F1/F3 remains protocol-based; the F1 fructose pulse and F3 viable-biomass step have direct trajectory evidence.
- One F2 ethanol point at 168 h is visible in the QC notebook but excluded by the conservative sugar/ethanol screen.
- Numerical convergence and full local rank do not imply that the current equations fit every state adequately. F1 glucose, residual YAN, late dead biomass, pyruvate, acetaldehyde, and acetate retain visible systematic error.
- `kPyrO2` is only indirectly informed because Lot 2 contains no DO observations and the oxygen-transfer parameters remain fixed. Treat its adjusted classification as conditional on that fixed oxygen submodel.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def finalize_existing_results() -> None:
    if not completed_outputs_available():
        raise RuntimeError("The completed Lot 2 result set is incomplete and cannot be finalized")
    augment_derived_diagnostics()
    summary_path = RESULTS_DIR / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["source_fingerprint"] = source_fingerprint()
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report()


def run_analysis(
    fit_nfev: int = 450,
    secondary_nfev: int = 300,
    profile_nfev: int = 50,
    grid_points: int = 3,
    step: float = FINITE_DIFFERENCE_STEP,
    profile_parameters: tuple[str, ...] = PROFILE_DEFAULT,
    reuse_if_current: bool = False,
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "run_progress.log").write_text("", encoding="utf-8")
    current_fingerprint = source_fingerprint()
    summary_path = RESULTS_DIR / "summary.json"
    if reuse_if_current and completed_outputs_available() and summary_path.exists():
        saved_summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if saved_summary.get("source_fingerprint") == current_fingerprint:
            log_progress("[resume] source fingerprint unchanged; using completed result set")
            augment_derived_diagnostics()
            write_report()
            log_progress("[done] historical synthetic plus Lot 2 analysis complete")
            return

    log_progress("[load] historical synthetic batches")
    _data, historical_batches, _historical_co2, _tables = historical.make_medium_batches("synthetic")
    log_progress("[load] processed MBDoE Lot 2 batches")
    lot2_batches, lot2_co2_times, support = make_lot2_batches()
    combined_batches = historical_batches + lot2_batches
    support.to_csv(RESULTS_DIR / "lot2_measurement_support.csv", index=False)
    pd.concat(
        [
            reference.summarize_batches(historical_batches, "historical_synthetic"),
            reference.summarize_batches(lot2_batches, "lot2_added"),
            reference.summarize_batches(combined_batches, "historical_synthetic_plus_lot2"),
        ],
        ignore_index=True,
    ).to_csv(RESULTS_DIR / "batch_summary.csv", index=False)

    checkpoint = load_theta_checkpoint()
    if checkpoint is not None:
        log_progress("[resume] loading completed core and secondary fits")
        theta_historical, theta_combined = checkpoint
        theta_core_combined = theta_combined
    else:
        theta_start = final.load_theta_final()
        theta_historical = load_historical_core_theta(theta_start)
        log_progress("[fit-secondary] historical synthetic baseline")
        theta_historical, secondary_summary_historical = fit_secondary_block(
            "historical_synthetic", historical_batches, theta_historical, secondary_nfev
        )

        log_progress("[fit-core] historical synthetic plus Lot 2")
        theta_core_combined, fit_summary, fit_candidates = reference.fit_core_case(
            "historical_synthetic_plus_lot2", combined_batches, theta_historical, fit_nfev
        )
        theta_combined = dict(theta_historical)
        theta_combined.update(theta_core_combined)
        log_progress("[fit-secondary] historical synthetic plus Lot 2")
        theta_combined, secondary_summary_combined = fit_secondary_block(
            "historical_synthetic_plus_lot2", combined_batches, theta_combined, secondary_nfev
        )

        baseline_fit_path = HISTORICAL_RESULTS_DIR / "core_fit_summary.csv"
        baseline_fit = pd.read_csv(baseline_fit_path) if baseline_fit_path.exists() else pd.DataFrame()
        if not baseline_fit.empty:
            baseline_fit["fit"] = "historical_synthetic"
        pd.concat([baseline_fit, fit_summary], ignore_index=True).to_csv(
            RESULTS_DIR / "fit_summary.csv", index=False
        )
        fit_candidates.to_csv(RESULTS_DIR / "combined_core_fit_candidate_summary.csv", index=False)
        pd.concat([secondary_summary_historical, secondary_summary_combined], ignore_index=True).to_csv(
            RESULTS_DIR / "secondary_fit_summary.csv", index=False
        )
        pd.DataFrame(theta_rows("historical_synthetic", theta_historical) + theta_rows("historical_synthetic_plus_lot2", theta_combined)).to_csv(
            RESULTS_DIR / "theta_by_case.csv", index=False
        )

    log_progress("[fim] historical observed-only at historical optimum")
    fim_historical = extended_fim(theta_historical, historical_batches, {}, False, step)
    log_progress("[fim] historical component at combined optimum")
    fim_hist_component = extended_fim(theta_combined, historical_batches, {}, False, step)
    log_progress("[fim] Lot 2 observed-only increment at combined optimum")
    fim_lot2_observed = extended_fim(theta_combined, lot2_batches, lot2_co2_times, False, step)
    fim_combined_observed = fim_hist_component + fim_lot2_observed
    log_progress("[fim] Lot 2 observed plus CO2-schedule potential")
    fim_lot2_co2 = extended_fim(theta_combined, lot2_batches, lot2_co2_times, True, step)
    fim_combined_co2 = fim_hist_component + fim_lot2_co2

    matrices = {
        "historical_observed_only_at_historical_theta": fim_historical,
        "historical_component_at_combined_theta": fim_hist_component,
        "lot2_increment_observed_only_at_combined_theta": fim_lot2_observed,
        "combined_observed_only_at_combined_theta": fim_combined_observed,
        "lot2_increment_observed_plus_co2_schedule_at_combined_theta": fim_lot2_co2,
        "combined_observed_plus_co2_schedule_at_combined_theta": fim_combined_co2,
    }
    matrix_filenames = {
        "historical_observed_only_at_historical_theta": "fim_hist_baseline.csv",
        "historical_component_at_combined_theta": "fim_hist_at_combined.csv",
        "lot2_increment_observed_only_at_combined_theta": "fim_lot2_observed.csv",
        "combined_observed_only_at_combined_theta": "fim_combined_observed.csv",
        "lot2_increment_observed_plus_co2_schedule_at_combined_theta": "fim_lot2_co2_potential.csv",
        "combined_observed_plus_co2_schedule_at_combined_theta": "fim_combined_co2_potential.csv",
    }
    for label, matrix in matrices.items():
        pd.DataFrame(matrix, index=EXTENDED_PARAMETERS, columns=EXTENDED_PARAMETERS).to_csv(
            RESULTS_DIR / matrix_filenames[label]
        )
    metrics = pd.DataFrame(
        [metric_row(label, "extended", matrix, EXTENDED_PARAMETERS) for label, matrix in matrices.items()]
    )

    log_progress("[fim-core] combined core block")
    jac_core, _ = base.build_jacobian(theta_combined, CORE_PARAMETERS, combined_batches, step=step)
    fim_core = jac_core.T @ jac_core
    pd.DataFrame(fim_core, index=CORE_PARAMETERS, columns=CORE_PARAMETERS).to_csv(
        RESULTS_DIR / "fim_core_combined.csv"
    )
    metrics = pd.concat(
        [metrics, pd.DataFrame([metric_row("historical_synthetic_plus_lot2", "core", fim_core, CORE_PARAMETERS)])],
        ignore_index=True,
    )
    metrics.to_csv(RESULTS_DIR / "fim_metrics.csv", index=False)

    estim_historical = final.parameter_estimability(
        fim_historical, theta_historical, EXTENDED_PARAMETERS, "historical_synthetic"
    ).rename(columns={"analysis": "case"})
    estim_combined = final.parameter_estimability(
        fim_combined_observed, theta_combined, EXTENDED_PARAMETERS, "historical_synthetic_plus_lot2"
    ).rename(columns={"analysis": "case"})
    estim_co2 = final.parameter_estimability(
        fim_combined_co2, theta_combined, EXTENDED_PARAMETERS, "historical_plus_lot2_CO2_schedule_potential"
    ).rename(columns={"analysis": "case"})
    estimability = pd.concat([estim_historical, estim_combined, estim_co2], ignore_index=True)
    estimability.to_csv(RESULTS_DIR / "estimability_by_case.csv", index=False)
    comparison = estim_historical.merge(
        estim_combined,
        on="parameter",
        suffixes=("_historical", "_combined"),
        validate="one_to_one",
    )
    comparison["std_log_ratio_combined_over_historical"] = (
        comparison["std_log_approx_combined"] / comparison["std_log_approx_historical"].replace(0.0, np.nan)
    )
    comparison = comparison.rename(
        columns={
            "std_log_approx_historical": "std_log_historical",
            "std_log_approx_combined": "std_log_combined",
        }
    )
    comparison.to_csv(RESULTS_DIR / "estimability_change_historical_vs_plus_lot2.csv", index=False)

    eig, weak = reference.eigen_tables(
        fim_combined_observed, EXTENDED_PARAMETERS, "historical_synthetic_plus_lot2"
    )
    eig.to_csv(RESULTS_DIR / "eigenvalues_combined.csv", index=False)
    weak.to_csv(RESULTS_DIR / "weak_directions_combined.csv", index=False)

    log_progress("[profile] combined core: " + ", ".join(profile_parameters))
    profile, profile_summary = reference.profile_core_case(
        "historical_synthetic_plus_lot2",
        combined_batches,
        theta_core_combined,
        profile_parameters,
        grid_points,
        profile_nfev,
    )
    profile.to_csv(RESULTS_DIR / "profile_likelihood_combined.csv", index=False)
    profile_summary.to_csv(RESULTS_DIR / "profile_likelihood_summary_combined.csv", index=False)

    fit_metrics = pd.concat(
        [
            batch_fit_metrics("historical_theta_on_lot2", lot2_batches, theta_historical),
            batch_fit_metrics("combined_theta_on_lot2", lot2_batches, theta_combined),
        ],
        ignore_index=True,
    )
    fit_metrics.to_csv(RESULTS_DIR / "lot2_fit_metrics_by_state.csv", index=False)
    plot_lot2_fits(lot2_batches, theta_combined)
    plot_estimability(comparison)

    decomposition_relative_error = float(
        np.linalg.norm(fim_combined_observed - (fim_hist_component + fim_lot2_observed))
        / max(np.linalg.norm(fim_combined_observed), 1.0)
    )
    summary = {
        "n_historical_synthetic_batches": len(historical_batches),
        "n_lot2_batches": len(lot2_batches),
        "n_combined_batches": len(combined_batches),
        "n_core_parameters_reestimated": len(CORE_PARAMETERS),
        "n_secondary_parameters_reestimated": len(SECONDARY_PARAMETERS),
        "n_extended_parameters_screened": len(EXTENDED_PARAMETERS),
        "n_lot2_co2_information_times": int(sum(len(values) for values in lot2_co2_times.values())),
        "co2_values_used_in_fit": False,
        "fim_decomposition_relative_error": decomposition_relative_error,
        "profile_parameters": list(profile_parameters),
        "fit_nfev": fit_nfev,
        "secondary_nfev": secondary_nfev,
        "profile_nfev": profile_nfev,
        "grid_points": grid_points,
        "finite_difference_log_step": step,
        "source_fingerprint": current_fingerprint,
    }
    (RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    augment_derived_diagnostics()
    write_report()
    log_progress("[done] historical synthetic plus Lot 2 analysis complete")


def create_notebook() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    cells = [
        nbformat.v4.new_markdown_cell(
            """# Historical synthetic prior plus MBDoE 2026 Lot 2

## tl;dr

This notebook measures how much the executed synthetic-must Lot 2 campaign changes calibration and practical estimability relative to the existing historical synthetic-must prior. It distinguishes observed-data information from the still-conditional CO2 flow schedule and keeps unsupported aroma directions visible rather than regularizing them into apparent identifiability."""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## Context & Methods

### Key assumptions

The observed-data Fisher information matrix is computed in log-parameter coordinates,

$$F(\theta)=J(\theta)^\mathsf{T}J(\theta),\qquad
J_{ij}=\frac{\partial r_i}{\partial\log\theta_j}.$$

The prior and Lot 2 contributions are evaluated at the same combined optimum to obtain

$$F_{\mathrm{combined}}=F_{\mathrm{historical}}+F_{\mathrm{Lot2}}.$$

This decomposition separates new experimental information from changes caused only by moving the parameter estimate. The core kinetic block and reduced secondary block are re-estimated sequentially. Profile likelihood is used for selected core parameters because local FIM uncertainty alone cannot establish bounded confidence intervals."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import json
import pandas as pd
from IPython.display import Image, Markdown, display

ROOT = Path.cwd()
while ROOT.name != "pyomo-doe" and ROOT.parent != ROOT:
    ROOT = ROOT.parent
if ROOT.name != "pyomo-doe":
    raise RuntimeError("Run this notebook from inside the pyomo-doe repository")

import sys
FERMENTATION_MODEL = ROOT / "fermentation_model"
if str(FERMENTATION_MODEL) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL))

from laboratory_2026 import run_estimability_historical_synthetic_plus_lot2 as analysis

analysis.run_analysis(reuse_if_current=True)
summary = json.loads((analysis.RESULTS_DIR / "summary.json").read_text(encoding="utf-8"))
summary"""
        ),
        nbformat.v4.new_markdown_cell("## Data"),
        nbformat.v4.new_code_cell(
            """support = pd.read_csv(analysis.RESULTS_DIR / "lot2_measurement_support.csv")
batch_summary = pd.read_csv(analysis.RESULTS_DIR / "batch_summary.csv")
display(batch_summary.groupby("case").agg(n_batches=("batch", "nunique"), n_core_obs=("n_core_obs", "sum"), n_secondary_obs=("n_secondary_obs", "sum")))
display(support.pivot_table(index="state", values="n_observations", aggfunc="sum"))"""
        ),
        nbformat.v4.new_markdown_cell("## Results"),
        nbformat.v4.new_code_cell(
            """fit = pd.read_csv(analysis.RESULTS_DIR / "fit_summary.csv")
secondary_fit = pd.read_csv(analysis.RESULTS_DIR / "secondary_fit_summary.csv")
metrics = pd.read_csv(analysis.RESULTS_DIR / "fim_metrics.csv")
comparison = pd.read_csv(analysis.RESULTS_DIR / "estimability_change_historical_vs_plus_lot2.csv")
profiles = pd.read_csv(analysis.RESULTS_DIR / "profile_likelihood_summary_combined.csv")
display(fit)
display(secondary_fit)
display(metrics)
display(comparison.sort_values("std_log_ratio_combined_over_historical").head(15))
display(profiles)"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=str(analysis.FIGURE_DIR / "estimability_historical_vs_plus_lot2.png")))
for process in ("F1", "F2", "F3"):
    display(Image(filename=str(analysis.FIGURE_DIR / "lot2_fit_overlays" / f"fit_lot2_{process}.png")))"""
        ),
        nbformat.v4.new_markdown_cell("## Takeaways"),
        nbformat.v4.new_code_cell(
            """display(Markdown(analysis.REPORT_PATH.read_text(encoding="utf-8")))"""
        ),
    ]
    notebook = nbformat.v4.new_notebook(cells=cells)
    notebook.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook.metadata["language_info"] = {"name": "python", "version": sys.version.split()[0]}
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 14400) -> None:
    notebook = nbformat.read(NOTEBOOK_PATH, as_version=4)
    with tempfile.TemporaryDirectory(prefix="lot2_estimability_runtime_") as runtime:
        old_runtime = os.environ.get("JUPYTER_RUNTIME_DIR")
        os.environ["JUPYTER_RUNTIME_DIR"] = runtime
        try:
            executed = NotebookClient(
                notebook,
                timeout=timeout,
                kernel_name="python3",
                resources={"metadata": {"path": str(SCRIPT_DIR.parent.parent)}},
            ).execute()
        finally:
            if old_runtime is None:
                os.environ.pop("JUPYTER_RUNTIME_DIR", None)
            else:
                os.environ["JUPYTER_RUNTIME_DIR"] = old_runtime
    nbformat.write(executed, EXECUTED_NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-execute", action="store_true")
    parser.add_argument("--timeout", type=int, default=14400)
    parser.add_argument("--fit-nfev", type=int, default=450)
    parser.add_argument("--secondary-nfev", type=int, default=300)
    parser.add_argument("--profile-nfev", type=int, default=50)
    parser.add_argument("--grid-points", type=int, default=3)
    args = parser.parse_args()
    create_notebook()
    if args.no_execute:
        run_analysis(
            fit_nfev=args.fit_nfev,
            secondary_nfev=args.secondary_nfev,
            profile_nfev=args.profile_nfev,
            grid_points=args.grid_points,
        )
    else:
        execute_notebook(timeout=args.timeout)
    print(EXECUTED_NOTEBOOK_PATH if not args.no_execute else NOTEBOOK_PATH)


if __name__ == "__main__":
    main()
