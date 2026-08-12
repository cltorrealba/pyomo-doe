from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import replace
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

from laboratory_2026 import run_estimability_old_vs_lot1 as reference
from laboratory_2026 import run_final_operational_doe_v2 as final
from shared import run_new_must_glycerol_estimability_doe as base
from shared import run_secondary_joint_campaign_doe as joint
from shared import run_secondary_v2_model_evaluation as secondary_v2
from shared.paths import DATA_DIR, LABORATORY_2026_RESULTS_DIR, SHARED_RESULTS_DIR


MEDIA = ("natural", "synthetic")
PROFILE_DEFAULT = ("Kd0", "qN")
SECONDARY_STATES = ("Pyr", "AcAld", "Acetate", "O2")
CORE_PARAMETERS = final.FERMENTATION_TARGETS
EXTENDED_PARAMETERS = final.TARGET_PARAMETERS
CORE_SECONDARY_PARAMETERS = CORE_PARAMETERS + final.SECONDARY_TARGETS
PROFILE_COLUMNS = (
    "case",
    "profiled_parameter",
    "theta_hat",
    "theta_value",
    "success",
    "objective",
    "lr_stat",
    "chi2_95_threshold",
    "error",
)
PROFILE_SUMMARY_COLUMNS = (
    "case",
    "parameter",
    "n_success",
    "max_lr_stat",
    "crosses_left_95",
    "crosses_right_95",
    "profile_identifiable_95",
    "chi2_95_threshold",
)

SECONDARY_DATA_PATH = (
    SHARED_RESULTS_DIR / "secondary_metabolite_data_review" / "secondary_metabolite_long.csv"
)
NATURAL_WORKBOOK_CANDIDATES = (
    DATA_DIR / "Laboratorio 2026" / "Vendimia_2026" / "mosto_natural_xthiol.xlsx",
    DATA_DIR / "Laboratorio 2026" / "mosto_natural_xthiol.xlsx",
)
EXTERNAL_PROCESS_ROOT = Path(
    "C:/Users/ctorrealba/OneDrive - Vi\u00f1a Concha y Toro S.A/Documentos/"
    "Proyectos I+D/PI-4497/Resultados/2026/Lecturas Tablero Laboratorio 2026"
)
LOCAL_PROCESS_FALLBACK = DATA_DIR / "Laboratorio 2026" / "raw_data"

STANDARD_MOLAR_VOLUME_L_MOL = 22.414
CO2_MOLAR_MASS_G_MOL = 44.01
DEFAULT_REACTOR_VOLUME_L = 2.0
MAX_CO2_INTEGRATION_GAP_H = 0.5
CO2_INFORMATION_INTERVAL_H = 6.0
MODEL_TEMPERATURE_GRID_H = 2.0


def results_dir(medium: str) -> Path:
    validate_medium(medium)
    return LABORATORY_2026_RESULTS_DIR / f"estimability_historical_{medium}"


def notebook_path(medium: str, executed: bool = False) -> Path:
    validate_medium(medium)
    suffix = ".executed.ipynb" if executed else ".ipynb"
    return SCRIPT_DIR / "notebooks" / f"fermentation_estimability_historical_{medium}{suffix}"


def validate_medium(medium: str) -> None:
    if medium not in MEDIA:
        raise ValueError(f"medium must be one of {MEDIA}; received {medium!r}")


def log_progress(medium: str, message: str) -> None:
    out = results_dir(medium)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "run_progress.log").open("a", encoding="utf-8") as handle:
        handle.write(message + "\n")
    print(message, flush=True)


def first_existing_path(paths: tuple[Path, ...]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError("None of the expected paths exists: " + "; ".join(map(str, paths)))


def load_historical_medium_data(medium: str) -> pd.DataFrame:
    validate_medium(medium)
    if SECONDARY_DATA_PATH.exists():
        data = pd.read_csv(SECONDARY_DATA_PATH)
        source = SECONDARY_DATA_PATH
    else:
        data = base.load_normalized_data()
        source = base.NORMALIZED_DATA_PATH
    data = data[data["medium"].astype(str).eq(medium)].copy()
    if data.empty:
        raise RuntimeError(f"No historical rows found for medium={medium!r} in {source}")
    numeric = set(base.STATE_COLUMNS.values()).union(
        {
            "time_h",
            "temperature_c",
            "N_pulse_kg_m3",
            "pyruvic_acid",
            "acetaldehyde",
            "acetic_acid",
            "DO_mg_l",
        }
    )
    for column in numeric.intersection(data.columns):
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data.sort_values(["batch", "time_h"]).reset_index(drop=True)


def load_natural_batch_metadata(batches: list[str]) -> pd.DataFrame:
    workbook = first_existing_path(NATURAL_WORKBOOK_CANDIDATES)
    rows: list[dict[str, object]] = []
    for batch in batches:
        frame = pd.read_excel(
            workbook,
            sheet_name=batch,
            usecols=lambda column: column in {"t", "fecha_hora", "Volumen inicial (L)"},
        )
        frame["t"] = pd.to_numeric(frame["t"], errors="coerce")
        zero = frame[frame["t"].eq(0.0)]
        if zero.empty:
            raise RuntimeError(f"No t=0 row with fecha_hora was found in {workbook.name}/{batch}")
        timestamp = pd.to_datetime(zero["fecha_hora"].iloc[0], errors="coerce")
        if pd.isna(timestamp):
            raise RuntimeError(f"Invalid t=0 fecha_hora in {workbook.name}/{batch}")
        volume = pd.to_numeric(zero.get("Volumen inicial (L)"), errors="coerce").dropna()
        rows.append(
            {
                "batch": batch,
                "t0": timestamp,
                "reactor_volume_l": float(volume.iloc[0]) if not volume.empty else DEFAULT_REACTOR_VOLUME_L,
                "reactor_volume_source": "workbook" if not volume.empty else "default_2L",
                "metadata_workbook": str(workbook),
            }
        )
    return pd.DataFrame(rows)


def batch_filename_pattern(batch: str) -> re.Pattern[str]:
    number = int(re.search(r"(\d+)$", batch).group(1))
    return re.compile(rf"lab0*{number}(?!\d)", flags=re.IGNORECASE)


def process_file_kind(path: Path) -> str | None:
    name = path.name.lower()
    if "co2_filt" in name or "co2f" in name:
        return "co2"
    if "temp" in name:
        return "temperature"
    return None


def discover_process_files(batch: str) -> list[dict[str, object]]:
    pattern = batch_filename_pattern(batch)
    roots = (
        ("external_primary", EXTERNAL_PROCESS_ROOT, 0),
        ("repository_fallback", LOCAL_PROCESS_FALLBACK, 1),
    )
    rows: list[dict[str, object]] = []
    for source, root, priority in roots:
        if not root.exists():
            continue
        for path in root.rglob("*.csv"):
            if not pattern.search(path.name):
                continue
            kind = process_file_kind(path)
            if kind is None:
                continue
            rows.append(
                {
                    "batch": batch,
                    "kind": kind,
                    "source": source,
                    "source_priority": priority,
                    "path": path,
                    "relative_path": str(path.relative_to(root)),
                    "part_number": int(re.search(r"pt(\d+)", path.name, flags=re.IGNORECASE).group(1))
                    if re.search(r"pt(\d+)", path.name, flags=re.IGNORECASE)
                    else 1,
                }
            )
    return rows


def read_process_file(file_row: dict[str, object]) -> tuple[pd.DataFrame, dict[str, object]]:
    path = Path(file_row["path"])
    kind = str(file_row["kind"])
    wanted = {"timestamp", "fermentador", "status", "is_spike"}
    wanted.add("flow_filt_sccm" if kind == "co2" else "T")
    if kind == "temperature":
        wanted.add("SP")
    frame = pd.read_csv(path, usecols=lambda column: column in wanted)
    if "timestamp" not in frame:
        raise RuntimeError(f"Missing timestamp column in {path}")
    value_column = "flow_filt_sccm" if kind == "co2" else "T"
    if value_column not in frame:
        raise RuntimeError(f"Missing {value_column} column in {path}")
    frame["timestamp"] = pd.to_datetime(frame["timestamp"], errors="coerce")
    for column in (value_column, "SP", "is_spike"):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    valid = frame.dropna(subset=["timestamp", value_column]).copy()
    valid["source"] = file_row["source"]
    valid["source_priority"] = file_row["source_priority"]
    valid["source_file"] = str(path)
    inventory = {
        **{key: value for key, value in file_row.items() if key != "path"},
        "path": str(path),
        "n_rows_raw": int(len(frame)),
        "n_rows_valid": int(len(valid)),
        "timestamp_min": valid["timestamp"].min() if not valid.empty else pd.NaT,
        "timestamp_max": valid["timestamp"].max() if not valid.empty else pd.NaT,
    }
    return valid, inventory


def co2_rate_from_sccm(flow_sccm: pd.Series, volume_l: float) -> pd.Series:
    gas_l_h = pd.to_numeric(flow_sccm, errors="coerce") * 60.0 / 1000.0
    total_g_h = gas_l_h * CO2_MOLAR_MASS_G_MOL / STANDARD_MOLAR_VOLUME_L_MOL
    return total_g_h / float(volume_l)


def resample_process_signal(
    frame: pd.DataFrame,
    kind: str,
    t0: pd.Timestamp,
    horizon_h: float,
    volume_l: float,
) -> pd.DataFrame:
    frame = frame.sort_values(["timestamp", "source_priority"]).drop_duplicates("timestamp", keep="first")
    frame["t_h"] = (frame["timestamp"] - t0).dt.total_seconds() / 3600.0
    frame = frame[frame["t_h"].between(0.0, horizon_h, inclusive="both")].copy()
    if frame.empty:
        return frame
    frequency = "10min" if kind == "co2" else "30min"
    numeric_columns = ["flow_filt_sccm", "T", "SP", "is_spike"]
    columns = [column for column in numeric_columns if column in frame]
    sampled = frame.set_index("timestamp")[columns].resample(frequency).median()
    sampled["n_raw_rows"] = frame.set_index("timestamp").resample(frequency).size()
    sampled = sampled.reset_index()
    sampled["t_h"] = (sampled["timestamp"] - t0).dt.total_seconds() / 3600.0
    value_column = "flow_filt_sccm" if kind == "co2" else "T"
    sampled = sampled.dropna(subset=[value_column]).copy()
    if kind == "co2":
        sampled["flow_filt_sccm_physical"] = sampled["flow_filt_sccm"].clip(lower=0.0)
        sampled["co2_offgas_equivalent_g_l_h"] = co2_rate_from_sccm(
            sampled["flow_filt_sccm_physical"], volume_l
        )
        dt = sampled["t_h"].diff()
        previous_rate = sampled["co2_offgas_equivalent_g_l_h"].shift()
        increment = 0.5 * (sampled["co2_offgas_equivalent_g_l_h"] + previous_rate) * dt
        valid_interval = dt.gt(0.0) & dt.le(MAX_CO2_INTEGRATION_GAP_H)
        sampled["co2_observed_increment_g_l"] = increment.where(valid_interval, 0.0).fillna(0.0)
        sampled["co2_cumulative_observed_lower_bound_g_l"] = sampled[
            "co2_observed_increment_g_l"
        ].cumsum()
        sampled["co2_segment"] = (~valid_interval).cumsum()
    return sampled.reset_index(drop=True)


def process_qc_row(
    batch: str,
    kind: str,
    raw: pd.DataFrame,
    sampled: pd.DataFrame,
    t0: pd.Timestamp,
    horizon_h: float,
) -> dict[str, object]:
    value_column = "flow_filt_sccm" if kind == "co2" else "T"
    raw = raw.sort_values(["timestamp", "source_priority"]).drop_duplicates(
        "timestamp", keep="first"
    )
    raw["t_h"] = (raw["timestamp"] - t0).dt.total_seconds() / 3600.0
    raw = raw[raw["t_h"].between(0.0, horizon_h, inclusive="both")]
    gaps_min = raw["timestamp"].diff().dt.total_seconds() / 60.0
    observed_h = (gaps_min.where(gaps_min.le(30.0), 0.0).sum() / 60.0) if not raw.empty else 0.0
    row: dict[str, object] = {
        "batch": batch,
        "kind": kind,
        "n_raw_rows_in_horizon": int(len(raw)),
        "n_resampled_rows": int(len(sampled)),
        "t_first_h": float(raw["t_h"].min()) if not raw.empty else np.nan,
        "t_last_h": float(raw["t_h"].max()) if not raw.empty else np.nan,
        "horizon_h": float(horizon_h),
        "observed_hours_gap_capped_30min": float(observed_h),
        "coverage_fraction": float(observed_h / horizon_h) if horizon_h > 0 else np.nan,
        "max_gap_min": float(gaps_min.max()) if gaps_min.notna().any() else np.nan,
        "value_min": float(raw[value_column].min()) if not raw.empty else np.nan,
        "value_median": float(raw[value_column].median()) if not raw.empty else np.nan,
        "value_max": float(raw[value_column].max()) if not raw.empty else np.nan,
        "selected_external_rows": int(raw["source"].eq("external_primary").sum())
        if "source" in raw
        else 0,
        "selected_fallback_rows": int(raw["source"].eq("repository_fallback").sum())
        if "source" in raw
        else 0,
    }
    if kind == "co2" and not raw.empty:
        row.update(
            {
                "negative_fraction": float(raw[value_column].lt(0.0).mean()),
                "zero_fraction": float(raw[value_column].eq(0.0).mean()),
                "spike_flag_fraction": float(raw["is_spike"].fillna(0.0).astype(bool).mean())
                if "is_spike" in raw
                else np.nan,
            }
        )
    return row


def load_natural_process_data(
    historical_data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    batches = sorted(historical_data["batch"].astype(str).unique())
    metadata = load_natural_batch_metadata(batches)
    meta_by_batch = metadata.set_index("batch")
    inventory_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    temperature_frames: list[pd.DataFrame] = []
    co2_frames: list[pd.DataFrame] = []
    co2_times: dict[str, np.ndarray] = {}

    for batch in batches:
        horizon_h = float(historical_data.loc[historical_data["batch"].eq(batch), "time_h"].max())
        t0 = pd.Timestamp(meta_by_batch.loc[batch, "t0"])
        volume_l = float(meta_by_batch.loc[batch, "reactor_volume_l"])
        discovered = discover_process_files(batch)
        for kind in ("temperature", "co2"):
            frames: list[pd.DataFrame] = []
            for file_row in [row for row in discovered if row["kind"] == kind]:
                frame, inventory = read_process_file(file_row)
                frames.append(frame)
                inventory_rows.append(inventory)
            if not frames:
                qc_rows.append(
                    {
                        "batch": batch,
                        "kind": kind,
                        "n_raw_rows_in_horizon": 0,
                        "n_resampled_rows": 0,
                        "horizon_h": horizon_h,
                        "coverage_fraction": 0.0,
                    }
                )
                continue
            raw = pd.concat(frames, ignore_index=True)
            sampled = resample_process_signal(raw, kind, t0, horizon_h, volume_l)
            sampled.insert(0, "batch", batch)
            sampled["t0"] = t0
            sampled["reactor_volume_l"] = volume_l
            qc_rows.append(process_qc_row(batch, kind, raw, sampled, t0, horizon_h))
            if kind == "temperature":
                temperature_frames.append(sampled)
            else:
                co2_frames.append(sampled)
                targets = np.arange(0.0, horizon_h + 1e-9, CO2_INFORMATION_INTERVAL_H)
                selected: list[float] = []
                available = sampled["t_h"].to_numpy(dtype=float)
                for target in targets:
                    if not len(available):
                        continue
                    nearest = float(available[np.argmin(np.abs(available - target))])
                    if abs(nearest - target) <= MAX_CO2_INTEGRATION_GAP_H:
                        selected.append(nearest)
                co2_times[f"natural/{batch}"] = np.asarray(sorted(set(selected)), dtype=float)

    inventory = pd.DataFrame(inventory_rows)
    qc = pd.DataFrame(qc_rows)
    temperature = pd.concat(temperature_frames, ignore_index=True) if temperature_frames else pd.DataFrame()
    co2 = pd.concat(co2_frames, ignore_index=True) if co2_frames else pd.DataFrame()
    return metadata, inventory, qc, temperature, co2, co2_times


def remap_observations(
    old_time: np.ndarray,
    observations: dict[str, np.ndarray],
    new_time: np.ndarray,
) -> dict[str, np.ndarray]:
    positions = {round(float(value), 8): index for index, value in enumerate(new_time)}
    remapped: dict[str, np.ndarray] = {}
    for state, values in observations.items():
        out = np.full(len(new_time), np.nan, dtype=float)
        for time_h, value in zip(old_time, np.asarray(values, dtype=float)):
            out[positions[round(float(time_h), 8)]] = value
        remapped[state] = out
    return remapped


def augment_natural_batches_with_temperature(
    batches: list[base.BatchData],
    temperature: pd.DataFrame,
) -> tuple[list[base.BatchData], pd.DataFrame]:
    out: list[base.BatchData] = []
    alignment_rows: list[dict[str, object]] = []
    for batch in batches:
        sensor = temperature[temperature["batch"].eq(batch.batch)].sort_values("t_h")
        if sensor.empty:
            out.append(batch)
            continue
        sensor_time = sensor["t_h"].to_numpy(dtype=float)
        sensor_temp = sensor["T"].to_numpy(dtype=float)
        model_grid = np.arange(0.0, float(batch.time.max()) + 1e-9, MODEL_TEMPERATURE_GRID_H)
        if "SP" in sensor:
            setpoint = pd.to_numeric(sensor["SP"], errors="coerce")
            changes = sensor.loc[setpoint.diff().abs().gt(0.25), "t_h"].to_numpy(dtype=float)
            model_grid = np.concatenate([model_grid, changes])
        new_time = np.asarray(sorted(set(batch.time).union(float(value) for value in model_grid)), dtype=float)
        actual_temperature = np.interp(new_time, sensor_time, sensor_temp)
        nominal_temperature = np.interp(new_time, batch.time, batch.temperature_c)
        actual_temperature = np.where(new_time < sensor_time.min(), nominal_temperature, actual_temperature)
        actual_temperature = np.where(new_time > sensor_time.max(), nominal_temperature, actual_temperature)

        sample_actual = np.interp(batch.time, sensor_time, sensor_temp)
        sample_actual = np.where(batch.time < sensor_time.min(), batch.temperature_c, sample_actual)
        sample_actual = np.where(batch.time > sensor_time.max(), batch.temperature_c, sample_actual)
        for time_h, nominal, actual in zip(batch.time, batch.temperature_c, sample_actual):
            alignment_rows.append(
                {
                    "batch": batch.batch,
                    "time_h": float(time_h),
                    "temperature_nominal_c": float(nominal),
                    "temperature_sensor_c": float(actual),
                    "sensor_minus_nominal_c": float(actual - nominal),
                }
            )
        out.append(
            replace(
                batch,
                time=new_time,
                temperature_c=actual_temperature,
                observations=remap_observations(batch.time, batch.observations, new_time),
            )
        )
    return out, pd.DataFrame(alignment_rows)


def make_medium_batches(
    medium: str,
) -> tuple[
    pd.DataFrame,
    list[base.BatchData],
    dict[str, np.ndarray],
    dict[str, pd.DataFrame],
]:
    data = load_historical_medium_data(medium)
    batches = joint.make_secondary_batches(data)
    process_tables: dict[str, pd.DataFrame] = {}
    co2_times: dict[str, np.ndarray] = {}
    if medium == "natural":
        metadata, inventory, qc, temperature, co2, co2_times = load_natural_process_data(data)
        batches, alignment = augment_natural_batches_with_temperature(batches, temperature)
        process_tables = {
            "metadata": metadata,
            "inventory": inventory,
            "qc": qc,
            "temperature": temperature,
            "co2": co2,
            "temperature_alignment": alignment,
        }
    return data, batches, co2_times, process_tables


def measurement_support(batches: list[base.BatchData]) -> pd.DataFrame:
    states = tuple(base.STATE_NAMES) + SECONDARY_STATES + tuple(joint.AROMA_COLUMNS.values())
    rows = []
    for state in states:
        counts = []
        for batch in batches:
            values = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            counts.append(int(np.isfinite(values).sum()))
        rows.append(
            {
                "state": state,
                "n_observations": int(sum(counts)),
                "n_batches_with_observations": int(sum(count > 0 for count in counts)),
                "used_in_core_fit": state in base.STATE_NAMES,
                "used_in_extended_fim": state in set(base.STATE_NAMES).union(SECONDARY_STATES),
            }
        )
    return pd.DataFrame(rows)


def fit_core_case_with_profile_polish(
    label: str,
    batches: list[base.BatchData],
    theta0: dict[str, float],
    max_nfev: int,
    polish_parameters: tuple[str, ...],
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    best_theta, fit_summary, candidate_summary = reference.fit_core_case(
        label, batches, theta0, max_nfev
    )
    best_summary = fit_summary.iloc[0].to_dict()
    best_seed = str(best_summary.get("best_seed", "reference_multistart"))
    candidates = candidate_summary.to_dict("records")

    for profiled in polish_parameters:
        if profiled not in CORE_PARAMETERS:
            continue
        nuisance = tuple(name for name in CORE_PARAMETERS if name != profiled)
        nuisance_seed = f"profile_center_{profiled}_nuisance"
        theta_nuisance, nuisance_summary = base.fit_parameters(
            f"{label}__{nuisance_seed}",
            batches,
            best_theta,
            nuisance,
            max_nfev=max_nfev,
        )
        candidates.append({"case": label, "seed": nuisance_seed, **nuisance_summary})
        if float(nuisance_summary["final_wsse"]) < float(best_summary["final_wsse"]):
            best_theta = theta_nuisance
            best_summary = dict(nuisance_summary)
            best_seed = nuisance_seed

        release_seed = f"profile_center_{profiled}_release"
        theta_release, release_summary = base.fit_parameters(
            f"{label}__{release_seed}",
            batches,
            theta_nuisance,
            CORE_PARAMETERS,
            max_nfev=max_nfev,
        )
        candidates.append({"case": label, "seed": release_seed, **release_summary})
        release_wsse = float(release_summary["final_wsse"])
        best_wsse = float(best_summary["final_wsse"])
        release_is_equivalent_full_fit = np.isclose(
            release_wsse, best_wsse, rtol=1e-12, atol=1e-9
        ) and int(release_summary["n_parameters"]) > int(best_summary["n_parameters"])
        if release_wsse < best_wsse or release_is_equivalent_full_fit:
            best_theta = theta_release
            best_summary = dict(release_summary)
            best_seed = release_seed

    best_summary["fit"] = label
    best_summary["best_seed"] = best_seed
    return best_theta, pd.DataFrame([best_summary]), pd.DataFrame(candidates)


def save_process_tables(medium: str, tables: dict[str, pd.DataFrame]) -> None:
    out = results_dir(medium)
    names = {
        "metadata": "process_batch_metadata.csv",
        "inventory": "process_file_inventory.csv",
        "qc": "process_signal_qc.csv",
        "temperature": "process_temperature_30min.csv",
        "co2": "process_co2_10min.csv",
        "temperature_alignment": "temperature_nominal_vs_sensor_at_model_times.csv",
    }
    for key, frame in tables.items():
        if key in names:
            frame.to_csv(out / names[key], index=False)


def estimability_plot(medium: str, estimability: pd.DataFrame) -> None:
    figure_dir = results_dir(medium) / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    data = estimability.sort_values("std_log_approx", ascending=False)
    colors = {
        "well_estimated": "#2a9d8f",
        "moderate": "#457b9d",
        "weak_but_actionable": "#e9c46a",
        "weak_or_confounded": "#e76f51",
    }
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(
        np.arange(len(data)),
        data["std_log_approx"],
        color=[colors.get(value, "#777777") for value in data["classification"]],
    )
    ax.axhline(0.35, color="#2a9d8f", linestyle="--", linewidth=1)
    ax.axhline(0.75, color="#e9c46a", linestyle="--", linewidth=1)
    ax.axhline(1.25, color="#e76f51", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(data)))
    ax.set_xticklabels(data["parameter"], rotation=70, ha="right")
    ax.set_ylabel("Approximate standard deviation in log(parameter)")
    ax.set_title(f"Historical {medium}: local practical estimability")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "estimability_std_log.png", dpi=170)
    plt.close(fig)


def eigen_plot(medium: str, eigenvalues: pd.DataFrame) -> None:
    figure_dir = results_dir(medium) / "figures"
    values = eigenvalues.sort_values("direction")
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.semilogy(values["direction"], np.maximum(values["eigenvalue"], 1e-18), marker="o")
    ax.set_xlabel("FIM eigen-direction (weakest to strongest)")
    ax.set_ylabel("Eigenvalue")
    ax.set_title(f"Historical {medium}: extended FIM spectrum")
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(figure_dir / "extended_fim_eigenvalues.png", dpi=170)
    plt.close(fig)


def profile_plot(medium: str, profile: pd.DataFrame) -> None:
    if profile.empty:
        return
    figure_dir = results_dir(medium) / "figures"
    parameters = list(profile["profiled_parameter"].dropna().unique())
    ncols = 2
    nrows = max(1, int(math.ceil(len(parameters) / ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4 * nrows), squeeze=False)
    for ax, parameter in zip(axes.ravel(), parameters):
        group = profile[profile["profiled_parameter"].eq(parameter)].sort_values("theta_value")
        valid = group[group["success"]]
        ax.plot(valid["theta_value"], valid["lr_stat"], marker="o")
        if not valid.empty:
            ax.axvline(float(valid["theta_hat"].iloc[0]), color="black", linestyle=":", linewidth=1)
        ax.axhline(reference.CHI2_95, color="#d62828", linestyle="--", linewidth=1, label="95% threshold")
        ax.set_title(f"{medium}: {parameter}")
        ax.set_xlabel("Fixed parameter value")
        ax.set_ylabel("Profile likelihood ratio statistic")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    for ax in axes.ravel()[len(parameters) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(figure_dir / "profile_likelihood.png", dpi=170)
    plt.close(fig)


def process_plot(medium: str, tables: dict[str, pd.DataFrame]) -> None:
    if medium != "natural" or not tables:
        return
    temperature = tables["temperature"]
    co2 = tables["co2"]
    batches = sorted(set(temperature.get("batch", [])).union(set(co2.get("batch", []))))
    if not batches:
        return
    figure_dir = results_dir(medium) / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(batches), 2, figsize=(14, 2.6 * len(batches)), squeeze=False)
    for row, batch in enumerate(batches):
        temp = temperature[temperature["batch"].eq(batch)]
        gas = co2[co2["batch"].eq(batch)]
        axes[row, 0].plot(temp["t_h"], temp["T"], color="#d1495b", linewidth=1.0, label="T measured")
        if "SP" in temp and temp["SP"].notna().any():
            axes[row, 0].plot(temp["t_h"], temp["SP"], color="#00798c", linewidth=1.0, alpha=0.8, label="setpoint")
        axes[row, 0].set_ylabel(f"{batch}\nTemperature [C]")
        axes[row, 0].grid(True, alpha=0.2)
        axes[row, 0].legend(fontsize=7, loc="best")
        axes[row, 1].plot(gas["t_h"], gas["flow_filt_sccm"], color="#edae49", linewidth=0.9, label="CO2_FILT raw sign")
        axes[row, 1].plot(gas["t_h"], gas["flow_filt_sccm_physical"], color="#003d5b", linewidth=0.9, alpha=0.8, label="physical clip")
        axes[row, 1].set_ylabel("CO2 flow [sccm]")
        axes[row, 1].grid(True, alpha=0.2)
        axes[row, 1].legend(fontsize=7, loc="best")
    axes[-1, 0].set_xlabel("Process time from t=0 [h]")
    axes[-1, 1].set_xlabel("Process time from t=0 [h]")
    fig.suptitle("Historical natural must: stitched process signals", y=1.0)
    fig.tight_layout()
    fig.savefig(figure_dir / "natural_process_signals.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def fit_diagnostics(
    medium: str,
    batches: list[base.BatchData],
    theta: dict[str, float],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    figure_dir = results_dir(medium) / "figures" / "fit_overlays"
    figure_dir.mkdir(parents=True, exist_ok=True)
    core_rows: list[dict[str, object]] = []
    secondary_rows: list[dict[str, object]] = []
    for batch in batches:
        core = base.simulate(batch, theta, batch.time)
        if core is None:
            continue
        fig, axes = plt.subplots(4, 2, figsize=(12, 13), squeeze=False)
        for ax, state in zip(axes.ravel(), base.STATE_NAMES):
            observed = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            mask = np.isfinite(observed)
            predicted = core.loc[batch.time, state].to_numpy(dtype=float)
            ax.plot(batch.time, predicted, color="#264653", linewidth=1.4, label="model")
            ax.scatter(batch.time[mask], observed[mask], color="#e76f51", s=18, label="data", zorder=3)
            ax.set_title(state)
            ax.grid(True, alpha=0.2)
            if mask.any():
                residual = predicted[mask] - observed[mask]
                core_rows.append(
                    {
                        "batch": batch.batch,
                        "state": state,
                        "n": int(mask.sum()),
                        "rmse": float(np.sqrt(np.mean(residual**2))),
                        "mae": float(np.mean(np.abs(residual))),
                    }
                )
        axes.ravel()[-1].axis("off")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Historical {medium}: core fit for {batch.batch}")
        fig.tight_layout()
        fig.savefig(figure_dir / f"core_fit_{batch.batch}.png", dpi=150)
        plt.close(fig)

        secondary = secondary_v2.integrate_secondary_v2(batch, theta, core)
        if secondary is None:
            continue
        has_secondary = any(
            np.isfinite(np.asarray(batch.observations.get(state, []), dtype=float)).any()
            for state in SECONDARY_STATES
        )
        if not has_secondary:
            continue
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), squeeze=False)
        for ax, state in zip(axes.ravel(), SECONDARY_STATES):
            observed = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            mask = np.isfinite(observed)
            predicted = secondary.loc[batch.time, state].to_numpy(dtype=float)
            ax.plot(batch.time, predicted, color="#457b9d", linewidth=1.4, label="conditional model")
            ax.scatter(batch.time[mask], observed[mask], color="#e63946", s=18, label="data", zorder=3)
            ax.set_title(state)
            ax.grid(True, alpha=0.2)
            if mask.any():
                residual = predicted[mask] - observed[mask]
                secondary_rows.append(
                    {
                        "batch": batch.batch,
                        "state": state,
                        "n": int(mask.sum()),
                        "rmse": float(np.sqrt(np.mean(residual**2))),
                        "mae": float(np.mean(np.abs(residual))),
                    }
                )
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f"Historical {medium}: secondary-state check for {batch.batch}")
        fig.tight_layout()
        fig.savefig(figure_dir / f"secondary_fit_{batch.batch}.png", dpi=150)
        plt.close(fig)
    return pd.DataFrame(core_rows), pd.DataFrame(secondary_rows)


def as_markdown(frame: pd.DataFrame) -> str:
    try:
        return frame.to_markdown(index=False)
    except Exception:
        return "```\n" + frame.to_string(index=False) + "\n```"


def read_csv_allow_empty(path: Path, columns: tuple[str, ...] = ()) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=list(columns))


def write_report(medium: str) -> None:
    out = results_dir(medium)
    fit = pd.read_csv(out / "core_fit_summary.csv")
    metrics = pd.read_csv(out / "fim_metrics.csv")
    estimability = pd.read_csv(out / "estimability_extended.csv")
    profiles = read_csv_allow_empty(
        out / "profile_likelihood_summary_core.csv", PROFILE_SUMMARY_COLUMNS
    )
    weak = pd.read_csv(out / "weak_directions_extended.csv")
    support = pd.read_csv(out / "measurement_support.csv")
    class_counts = (
        estimability.groupby("classification", dropna=False).size().rename("n_parameters").reset_index()
    )
    best = estimability.sort_values("std_log_approx").head(10)
    worst = estimability.sort_values("std_log_approx", ascending=False).head(10)
    process_note = (
        "Actual temperature records replace the nominal temperature input. Valid CO2 sensor times "
        "enter the extended FIM as an online measurement schedule, while CO2 values remain diagnostic."
        if medium == "natural"
        else "No matching continuous process-sensor archive was supplied for these synthetic historical batches."
    )
    report = f"""# Historical {medium} estimability analysis

## Scope

This is an independent calibration and estimability analysis using only the historical `{medium}` medium subset. It is not a row-filtered view of the joint historical fit.

{process_note}

The core kinetic block is re-estimated. Secondary and aroma parameters are evaluated conditionally at the shared prior values, matching the logic of `fermentation_estimability_old_vs_lot1`.

## Core fit

{as_markdown(fit)}

## Measurement support

{as_markdown(support)}

## FIM metrics

{as_markdown(metrics)}

## Estimability classes

{as_markdown(class_counts)}

### Strongest local directions by marginal parameter uncertainty

{as_markdown(best[["parameter", "theta", "std_log_approx", "approx_95_multiplier", "active_bound", "classification"]])}

### Weakest or confounded marginal directions

{as_markdown(worst[["parameter", "theta", "std_log_approx", "approx_95_multiplier", "active_bound", "classification"]])}

## Profile likelihood screening

{as_markdown(profiles)}

## Weak eigen-directions

{as_markdown(weak.head(8))}

## Interpretation limits

- FIM classifications are local to the fitted/prior parameter point and use the same log-parameter finite-difference convention as the reference notebook.
- Aroma parameters without an aroma residual block remain unsupported by these notebooks even if aroma columns exist in the source workbook.
- Online CO2 is a gas-flow measurement. Its values are not treated as direct observations of the current cumulative production state until the gas-liquid observation model is calibrated.
- The secondary-state curves are conditional checks, not a new secondary-parameter calibration.
"""
    (out / f"estimability_historical_{medium}_report.md").write_text(report, encoding="utf-8")


def run_medium_analysis(
    medium: str,
    profile_parameters: tuple[str, ...] = PROFILE_DEFAULT,
    grid_points: int = 3,
    fit_nfev: int = 300,
    profile_nfev: int = 45,
    step: float = 0.04,
) -> None:
    validate_medium(medium)
    out = results_dir(medium)
    out.mkdir(parents=True, exist_ok=True)
    (out / "run_progress.log").write_text("", encoding="utf-8")

    log_progress(medium, f"[load] historical {medium}")
    data, batches, co2_times, process_tables = make_medium_batches(medium)
    if not batches:
        raise RuntimeError(f"No valid model batches were created for medium={medium}")
    reference.summarize_batches(batches, f"historical_{medium}_only").to_csv(
        out / "batch_summary.csv", index=False
    )
    measurement_support(batches).to_csv(out / "measurement_support.csv", index=False)
    if process_tables:
        save_process_tables(medium, process_tables)
        process_plot(medium, process_tables)

    theta0 = final.load_theta_final()
    case = f"historical_{medium}_only"
    log_progress(medium, f"[fit-core] {case}")
    theta_core, fit_summary, fit_candidates = fit_core_case_with_profile_polish(
        case,
        batches,
        theta0,
        fit_nfev,
        polish_parameters=profile_parameters,
    )
    fit_summary.to_csv(out / "core_fit_summary.csv", index=False)
    fit_candidates.to_csv(out / "core_fit_candidate_summary.csv", index=False)
    theta_all = dict(theta0)
    theta_all.update(theta_core)
    pd.DataFrame(
        [
            {"parameter": name, "theta": float(theta_all[name]), "reestimated_in_this_notebook": name in CORE_PARAMETERS}
            for name in EXTENDED_PARAMETERS
        ]
    ).to_csv(out / "theta.csv", index=False)

    log_progress(medium, f"[fim-core] {case}")
    jacobian_core, _ = base.build_jacobian(theta_core, CORE_PARAMETERS, batches, step=step)
    fim_core = jacobian_core.T @ jacobian_core
    core_metrics = joint.fim_metrics(fim_core, CORE_PARAMETERS)
    core_metrics.update(
        {"case": case, "block": "core_observed", "n_parameters": len(CORE_PARAMETERS)}
    )
    pd.DataFrame(fim_core, index=CORE_PARAMETERS, columns=CORE_PARAMETERS).to_csv(
        out / "fim_core.csv"
    )

    log_progress(medium, f"[fim-extended] {case}")
    fim_extended = final.finite_difference_fim(
        theta_all,
        EXTENDED_PARAMETERS,
        lambda theta: reference.observed_residual_vector(
            theta,
            batches,
            co2_times,
            include_secondary=True,
            include_co2_information=bool(co2_times),
        ),
        step,
    )
    extended_metrics = joint.fim_metrics(fim_extended, EXTENDED_PARAMETERS)
    extended_metrics.update(
        {
            "case": case,
            "block": "full_26_candidate_plus_co2_schedule"
            if co2_times
            else "full_26_candidate",
            "n_parameters": len(EXTENDED_PARAMETERS),
        }
    )
    supported_indices = [EXTENDED_PARAMETERS.index(name) for name in CORE_SECONDARY_PARAMETERS]
    fim_core_secondary = fim_extended[np.ix_(supported_indices, supported_indices)]
    core_secondary_metrics = joint.fim_metrics(
        fim_core_secondary, CORE_SECONDARY_PARAMETERS
    )
    core_secondary_metrics.update(
        {
            "case": case,
            "block": "core_secondary_observed_plus_co2_schedule"
            if co2_times
            else "core_secondary_observed",
            "n_parameters": len(CORE_SECONDARY_PARAMETERS),
        }
    )
    pd.DataFrame([core_metrics, core_secondary_metrics, extended_metrics]).to_csv(
        out / "fim_metrics.csv", index=False
    )
    pd.DataFrame(
        fim_core_secondary,
        index=CORE_SECONDARY_PARAMETERS,
        columns=CORE_SECONDARY_PARAMETERS,
    ).to_csv(out / "fim_core_secondary.csv")
    pd.DataFrame(fim_extended, index=EXTENDED_PARAMETERS, columns=EXTENDED_PARAMETERS).to_csv(
        out / "fim_extended.csv"
    )
    estimability = reference.estimability_from_fim(
        fim_extended, theta_all, EXTENDED_PARAMETERS, case
    )
    estimability.to_csv(out / "estimability_extended.csv", index=False)
    eigenvalues, weak = reference.eigen_tables(fim_extended, EXTENDED_PARAMETERS, case)
    eigenvalues.to_csv(out / "eigenvalues_extended.csv", index=False)
    weak.to_csv(out / "weak_directions_extended.csv", index=False)
    estimability_plot(medium, estimability)
    eigen_plot(medium, eigenvalues)

    log_progress(medium, f"[profile] {case}: {', '.join(profile_parameters)}")
    profile, profile_summary = reference.profile_core_case(
        case,
        batches,
        theta_core,
        profile_parameters,
        grid_points=grid_points,
        max_nfev=profile_nfev,
    )
    if profile.empty and not len(profile.columns):
        profile = pd.DataFrame(columns=PROFILE_COLUMNS)
    if profile_summary.empty and not len(profile_summary.columns):
        profile_summary = pd.DataFrame(columns=PROFILE_SUMMARY_COLUMNS)
    profile.to_csv(out / "profile_likelihood_core.csv", index=False)
    profile_summary.to_csv(out / "profile_likelihood_summary_core.csv", index=False)
    profile_plot(medium, profile)

    log_progress(medium, f"[fit-diagnostics] {case}")
    core_fit_metrics, secondary_fit_metrics = fit_diagnostics(medium, batches, theta_all)
    core_fit_metrics.to_csv(out / "core_fit_metrics_by_batch_state.csv", index=False)
    secondary_fit_metrics.to_csv(out / "secondary_fit_metrics_by_batch_state.csv", index=False)

    summary = {
        "medium": medium,
        "case": case,
        "n_batches": len(batches),
        "n_source_rows": int(len(data)),
        "profile_parameters": list(profile_parameters),
        "grid_points": int(grid_points),
        "fit_nfev": int(fit_nfev),
        "profile_nfev": int(profile_nfev),
        "finite_difference_log_step": float(step),
        "uses_actual_process_temperature": medium == "natural" and bool(process_tables),
        "n_co2_information_times": int(sum(len(values) for values in co2_times.values())),
        "co2_values_used_in_parameter_fit": False,
        "external_process_root": str(EXTERNAL_PROCESS_ROOT) if medium == "natural" else None,
        "local_process_fallback": str(LOCAL_PROCESS_FALLBACK) if medium == "natural" else None,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_report(medium)
    log_progress(medium, f"[done] historical {medium} analysis complete")


def create_notebook(medium: str) -> None:
    validate_medium(medium)
    path = notebook_path(medium)
    path.parent.mkdir(parents=True, exist_ok=True)
    process_text = (
        "For natural must, actual temperature is reconstructed from all timestamped fragments. "
        "The external process archive has priority; repository copies fill missing, non-conflicting fragments. "
        "Filtered CO2 flow is stitched and audited, but its values are not fitted as cumulative CO2."
        if medium == "natural"
        else "No continuous process-sensor files are added to the synthetic historical subset in this notebook."
    )
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            f"""# Historical {medium.title()} Must: Calibration and Estimability

This notebook independently reproduces the calibration, Fisher information, eigen-direction and profile-likelihood workflow used in `fermentation_estimability_old_vs_lot1`, restricted to the historical `{medium}` medium.

## tl;dr

The numerical conclusions are generated below from a fresh medium-specific fit. The final report states which parameters are locally supported, which remain confounded, and how well the fitted model reproduces each observed state.

## Context & Methods

The core parameter block is re-estimated in log space using the same bounds, residual scaling and deterministic multistart as the reference notebook. The local Fisher information matrix is

$$
F(\\theta) = J(\\theta)^\\mathsf{{T}} J(\\theta),
$$

where $J$ is the Jacobian of scaled residuals with respect to log-parameters. Profile likelihood is used as a nonlinear check for `Kd0` and `qN`; the three-point grid is a screening profile, not a final confidence interval.

{process_text}

### Key Assumptions

- Secondary and aroma parameters are evaluated conditionally at the shared prior values; only the core fermentation block is re-estimated here.
- CO2 flow and cumulative CO2 production are different observables. Sensor values remain diagnostic until a gas-liquid observation layer is calibrated.
- Missing observations stay missing. They are never replaced with zeros.
"""
        ),
        nbformat.v4.new_code_cell(
            f"""from pathlib import Path
import pandas as pd
from IPython.display import display, Image

from laboratory_2026 import run_estimability_historical_by_medium as analysis

MEDIUM = {medium!r}
RESULTS = analysis.results_dir(MEDIUM)
"""
        ),
        nbformat.v4.new_markdown_cell("## Run Analysis"),
        nbformat.v4.new_code_cell(
            """analysis.run_medium_analysis(
    MEDIUM,
    profile_parameters=("Kd0", "qN"),
    grid_points=3,
    fit_nfev=300,
    profile_nfev=45,
    step=0.04,
)"""
        ),
        nbformat.v4.new_markdown_cell("## Data"),
        nbformat.v4.new_code_cell(
            """batch_summary = pd.read_csv(RESULTS / "batch_summary.csv")
support = pd.read_csv(RESULTS / "measurement_support.csv")
display(batch_summary)
display(support)"""
        ),
    ]
    if medium == "natural":
        nb.cells.extend(
            [
                nbformat.v4.new_markdown_cell("### Process-Signal Audit"),
                nbformat.v4.new_code_cell(
                    """process_qc = pd.read_csv(RESULTS / "process_signal_qc.csv")
file_inventory = pd.read_csv(RESULTS / "process_file_inventory.csv")
temperature_alignment = pd.read_csv(RESULTS / "temperature_nominal_vs_sensor_at_model_times.csv")
display(process_qc)
display(file_inventory[["batch", "kind", "source", "part_number", "relative_path", "n_rows_valid", "timestamp_min", "timestamp_max"]])
display(temperature_alignment.groupby("batch")["sensor_minus_nominal_c"].agg(["mean", "std", "min", "max"]))"""
                ),
                nbformat.v4.new_code_cell(
                    """display(Image(filename=str(RESULTS / "figures" / "natural_process_signals.png")))"""
                ),
            ]
        )
    nb.cells.extend(
        [
            nbformat.v4.new_markdown_cell("## Results"),
            nbformat.v4.new_code_cell(
                """fit_candidates = pd.read_csv(RESULTS / "core_fit_candidate_summary.csv")
fit_summary = pd.read_csv(RESULTS / "core_fit_summary.csv")
fim_metrics = pd.read_csv(RESULTS / "fim_metrics.csv")
display(fit_candidates[["case", "seed", "final_wsse", "n_residuals", "nfev", "success"]])
display(fit_summary)
display(fim_metrics)"""
            ),
            nbformat.v4.new_markdown_cell("### Parameter Estimability"),
            nbformat.v4.new_code_cell(
                """estimability = pd.read_csv(RESULTS / "estimability_extended.csv")
display(estimability.sort_values("std_log_approx"))"""
            ),
            nbformat.v4.new_code_cell(
                """display(Image(filename=str(RESULTS / "figures" / "estimability_std_log.png")))"""
            ),
            nbformat.v4.new_markdown_cell("### Eigenvalue Diagnostics"),
            nbformat.v4.new_code_cell(
                """weak = pd.read_csv(RESULTS / "weak_directions_extended.csv")
display(weak)
display(Image(filename=str(RESULTS / "figures" / "extended_fim_eigenvalues.png")))"""
            ),
            nbformat.v4.new_markdown_cell("### Profile Likelihood"),
            nbformat.v4.new_code_cell(
                """profile_summary = pd.read_csv(RESULTS / "profile_likelihood_summary_core.csv")
profile = pd.read_csv(RESULTS / "profile_likelihood_core.csv")
display(profile_summary)
display(profile)"""
            ),
            nbformat.v4.new_code_cell(
                """display(Image(filename=str(RESULTS / "figures" / "profile_likelihood.png")))"""
            ),
            nbformat.v4.new_markdown_cell("### Curve-Fit Checks"),
            nbformat.v4.new_code_cell(
                """core_curve_metrics = pd.read_csv(RESULTS / "core_fit_metrics_by_batch_state.csv")
secondary_curve_metrics = pd.read_csv(RESULTS / "secondary_fit_metrics_by_batch_state.csv")
display(core_curve_metrics.groupby("state").agg(n=("n", "sum"), median_RMSE=("rmse", "median"), max_RMSE=("rmse", "max")))
display(secondary_curve_metrics.groupby("state").agg(n=("n", "sum"), median_RMSE=("rmse", "median"), max_RMSE=("rmse", "max")))"""
            ),
            nbformat.v4.new_markdown_cell("## Takeaways"),
            nbformat.v4.new_code_cell(
                f"""report = (RESULTS / "estimability_historical_{medium}_report.md").read_text(encoding="utf-8")
print(report)"""
            ),
        ]
    )
    nb.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    nb.metadata["language_info"] = {"name": "python", "version": "3.12"}
    nbformat.write(nb, path)


def execute_notebook(medium: str, timeout: int = 10800) -> None:
    create_notebook(medium)
    path = notebook_path(medium)
    runtime_dir = Path(tempfile.gettempdir()) / "pyomo_doe_jupyter_runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime_dir)
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(FERMENTATION_MODEL_DIR)}},
    )
    client.execute()
    nbformat.write(notebook, notebook_path(medium, executed=True))


def write_cross_medium_comparison() -> None:
    comparison_dir = LABORATORY_2026_RESULTS_DIR / "estimability_historical_by_medium"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    estimability_frames = []
    fit_frames = []
    for medium in MEDIA:
        out = results_dir(medium)
        estimability = pd.read_csv(out / "estimability_extended.csv")
        estimability["medium"] = medium
        estimability_frames.append(estimability)
        fit = pd.read_csv(out / "core_fit_summary.csv")
        fit["medium"] = medium
        fit_frames.append(fit)
    estimability_all = pd.concat(estimability_frames, ignore_index=True)
    estimability_all.to_csv(comparison_dir / "estimability_by_medium_long.csv", index=False)
    fit_all = pd.concat(fit_frames, ignore_index=True)
    fit_all.to_csv(comparison_dir / "core_fit_by_medium.csv", index=False)
    wide = estimability_all.pivot_table(
        index="parameter", columns="medium", values="std_log_approx", aggfunc="first"
    ).reset_index()
    if {"natural", "synthetic"}.issubset(wide.columns):
        wide["natural_over_synthetic_std_ratio"] = wide["natural"] / wide["synthetic"].replace(0.0, np.nan)
    wide.to_csv(comparison_dir / "estimability_std_log_comparison.csv", index=False)

    order = wide.set_index("parameter")[[column for column in MEDIA if column in wide]].max(axis=1).sort_values(ascending=False).index
    plot = wide.set_index("parameter").reindex(order)
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(plot))
    width = 0.38
    for offset, medium in enumerate(MEDIA):
        if medium in plot:
            ax.bar(x + (offset - 0.5) * width, plot[medium], width=width, label=medium)
    ax.set_xticks(x)
    ax.set_xticklabels(plot.index, rotation=70, ha="right")
    ax.set_ylabel("Approximate standard deviation in log(parameter)")
    ax.set_title("Independent historical estimability by medium")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(comparison_dir / "estimability_by_medium.png", dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run independent historical estimability analyses for natural and synthetic must."
    )
    parser.add_argument("--medium", choices=(*MEDIA, "all"), default="all")
    parser.add_argument("--create-notebooks", action="store_true")
    parser.add_argument("--execute-notebooks", action="store_true")
    parser.add_argument("--skip-analysis", action="store_true")
    parser.add_argument("--grid-points", type=int, default=3)
    parser.add_argument("--fit-nfev", type=int, default=300)
    parser.add_argument("--profile-nfev", type=int, default=45)
    parser.add_argument("--step", type=float, default=0.04)
    parser.add_argument("--profile-parameters", nargs="*", default=list(PROFILE_DEFAULT))
    args = parser.parse_args()

    selected = MEDIA if args.medium == "all" else (args.medium,)
    if args.execute_notebooks:
        for medium in selected:
            execute_notebook(medium)
    else:
        for medium in selected:
            if args.create_notebooks:
                create_notebook(medium)
            if not args.skip_analysis:
                run_medium_analysis(
                    medium,
                    profile_parameters=tuple(args.profile_parameters),
                    grid_points=args.grid_points,
                    fit_nfev=args.fit_nfev,
                    profile_nfev=args.profile_nfev,
                    step=args.step,
                )
    if args.medium == "all" and all((results_dir(medium) / "estimability_extended.csv").exists() for medium in MEDIA):
        write_cross_medium_comparison()


if __name__ == "__main__":
    main()
