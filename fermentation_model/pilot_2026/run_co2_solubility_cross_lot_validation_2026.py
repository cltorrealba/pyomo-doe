from __future__ import annotations

"""Pilot-2026 reproduction of the laboratory CO2 transfer-layer analysis.

The workflow calibrates the continuous CO2-release observation layer separately
on Pilot Lot 2 and Lot 3, validates each fit on an internal holdout and transfers
the fitted layer across lots.  It reuses the upstream Pilot-2026 kinetic fit and
the laboratory implementation of the CO2 solubility/O2/nitrogen-transition
equations.  Raw Pilot Lot-1 CO2 remains QC-only by owner decision.

Important: the MassView acquisition starts about 2.4--2.7 h after the first
chemical sample.  Artificial leading zeros are therefore excluded rather than
treated as sensor observations.  No claim of pre-inoculation sensor-zero
identification is made.
"""

import argparse
import json
import math
import os
import sys
from dataclasses import replace
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from nbclient import NotebookClient
from scipy.optimize import least_squares
from scipy.signal import savgol_filter


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SCRIPT_DIR.parent
ROOT_DIR = FERMENTATION_MODEL_DIR.parent
LAB_DIR = FERMENTATION_MODEL_DIR / "laboratory_2026"
for path in (FERMENTATION_MODEL_DIR, SCRIPT_DIR, LAB_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from shared import run_new_must_glycerol_estimability_doe as base
from laboratory_2026 import run_co2_matrix_cross_validation_2026 as lab_co2
from laboratory_2026 import run_estimability_historical_by_medium as natural_loader
from pilot_2026.adaptive_design import pilot_calibration


MODEL_DATASET_DIR = (
    SCRIPT_DIR
    / "results"
    / "adaptive_design_2026"
    / "model_dataset"
    / "20260717T125004Z_b03413"
)
PILOT_CALIBRATION_DIR = (
    SCRIPT_DIR
    / "results"
    / "adaptive_design_2026"
    / "baseline_calibration"
    / "20260717T153612Z_b5b67b"
)
INTEGRATION_DIR = SCRIPT_DIR / "results" / "data_integration_2026"
LAB_RESULTS_DIR = (
    FERMENTATION_MODEL_DIR
    / "laboratory_2026"
    / "results"
    / "co2_matrix_cross_validation_2026"
)
RESULTS_DIR = SCRIPT_DIR / "results" / "co2_solubility_cross_lot_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_co2_solubility_cross_lot_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_co2_solubility_cross_lot_validation.executed.ipynb"
)

MODEL_NAME = lab_co2.MODEL_NAME
THRESHOLD_MODEL_NAME = lab_co2.THRESHOLD_RELEASE_MODEL_NAME
LEGACY_MODEL_NAME = lab_co2.LEGACY_MODEL_NAME
MATRICES = ("pilot_lot2", "pilot_lot3")
HOLDOUTS = {"pilot_lot2": "26158", "pilot_lot3": "26211"}
VALID_RUNS = ("26157", "26158", "26159", "26210", "26211", "26212")
LOT_BY_RUN = {
    "26157": "pilot_lot2",
    "26158": "pilot_lot2",
    "26159": "pilot_lot2",
    "26210": "pilot_lot3",
    "26211": "pilot_lot3",
    "26212": "pilot_lot3",
}

CO2_MOLAR_MASS_G_MOL = 44.01
NORMAL_MOLAR_VOLUME_L_MOL = 22.414
REACTOR_VOLUME_L = 230.0
LN_MIN_TO_G_L_H = (
    60.0 * CO2_MOLAR_MASS_G_MOL / NORMAL_MOLAR_VOLUME_L_MOL / REACTOR_VOLUME_L
)
YAN_MODEL_PULSE_KG_M3 = 0.080
YAN_MASS_AUDIT_KG_M3 = 0.090434783
SAMPLE_WINDOW_H = 1.0
LOCAL_FILTER_WINDOW_POINTS = 13
LOCAL_FILTER_MIN_ABS_G_L_H = 0.08
PULSE_PROTECTION_H = 4.0
COLD_THRESHOLD_C = 16.5
DETECTION_LIMIT_STANDARD_G_L_H = 0.05
DETECTION_LIMIT_COLD_G_L_H = 0.08
RANDOM_SEED = 20260820

COLORS = {
    "ink": "#27313B",
    "grid": "#D8DEE5",
    "blue": "#326FA8",
    "blue_light": "#A7C7E7",
    "gold": "#D6A23A",
    "orange": "#D97732",
    "pink": "#B45A7A",
    "grey": "#7B8794",
    "light_grey": "#D8DEE5",
}


def _read_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _first_finite(values: pd.Series | np.ndarray, default: float) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(finite[0]) if len(finite) else float(default)


def _pilot_upstream_context() -> tuple[dict[str, base.BatchData], dict[str, float], pilot_calibration.CalibrationTables]:
    """Load the fixed Pilot-2026 upstream calibration and construct BatchData."""

    tables = pilot_calibration.load_tables(MODEL_DATASET_DIR)
    calibration_config = json.loads(
        (PILOT_CALIBRATION_DIR / "calibration_config.json").read_text(encoding="utf-8")
    )
    parameter_table = pd.read_csv(PILOT_CALIBRATION_DIR / "primary_parameter_estimates.csv")
    parameter_values = dict(
        zip(parameter_table["parameter"].astype(str), parameter_table["estimate"].astype(float))
    )
    nuisance = {
        "oculyze_biomass_scale_kg_m3_per_million_cells_ml": float(
            parameter_values["oculyze_biomass_scale_kg_m3_per_million_cells_ml"]
        )
    }
    batches = pilot_calibration.build_batches(tables, nuisance, calibration_config)
    theta = dict(base.DEFAULT_THETA)
    for name in calibration_config["primary_fit"]["parameters"]:
        theta[name] = float(parameter_values[name])

    # The existing upstream adapter interpolates temperature only at chemistry
    # times.  Merge the complete Sonda1 grid so programmed transitions are not
    # smoothed away inside temperature_at().
    augmented: dict[str, base.BatchData] = {}
    for batch in batches:
        temperature = (
            tables.temperature[tables.temperature["experiment_id"].astype(str).eq(batch.batch)]
            .dropna(subset=["time_h", "executed_temperature_c"])
            .sort_values("time_h")
            .drop_duplicates("time_h", keep="last")
        )
        sensor_time = temperature["time_h"].to_numpy(dtype=float)
        sensor_temperature = temperature["executed_temperature_c"].to_numpy(dtype=float)
        new_time = np.asarray(
            sorted(set(batch.time.astype(float)).union(sensor_time.astype(float))), dtype=float
        )
        actual_temperature = np.interp(new_time, sensor_time, sensor_temperature)
        remapped = natural_loader.remap_observations(
            np.asarray(batch.time, dtype=float), batch.observations, new_time
        )
        augmented[batch.batch] = replace(
            batch,
            time=new_time,
            temperature_c=actual_temperature,
            observations=remapped,
        )
    return augmented, theta, tables


def _primary_wide(primary: pd.DataFrame) -> pd.DataFrame:
    wide = primary.pivot_table(
        index=["experiment_id", "timestamp", "time_h"],
        columns="state",
        values="observed_value",
        aggfunc="first",
    ).reset_index()
    wide.columns.name = None
    wide["experiment_id"] = wide["experiment_id"].astype(str)
    return wide.sort_values(["experiment_id", "time_h"]).reset_index(drop=True)


def _density_crossing(group: pd.DataFrame, target: float = 1040.0) -> dict[str, float]:
    valid = (
        group[["time_h", "density"]]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .sort_values("time_h")
        .drop_duplicates("time_h", keep="last")
    )
    time_h = valid["time_h"].to_numpy(dtype=float)
    density = valid["density"].to_numpy(dtype=float)
    for index in range(len(valid) - 1):
        left = density[index] - target
        right = density[index + 1] - target
        if abs(left) <= 1e-9:
            return {
                "pulse_time_h": float(time_h[index]),
                "bracket_start_h": float(time_h[index]),
                "bracket_end_h": float(time_h[index]),
                "bracket_width_h": 0.0,
            }
        if left * right <= 0.0 and density[index + 1] != density[index]:
            fraction = (target - density[index]) / (density[index + 1] - density[index])
            crossing = time_h[index] + fraction * (time_h[index + 1] - time_h[index])
            return {
                "pulse_time_h": float(crossing),
                "bracket_start_h": float(time_h[index]),
                "bracket_end_h": float(time_h[index + 1]),
                "bracket_width_h": float(time_h[index + 1] - time_h[index]),
            }
    return {
        "pulse_time_h": np.nan,
        "bracket_start_h": np.nan,
        "bracket_end_h": np.nan,
        "bracket_width_h": np.nan,
    }


def build_pulse_schedule(primary: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
    wide = _primary_wide(primary)
    rows: list[dict[str, object]] = []
    for run in sorted(wide["experiment_id"].unique()):
        chemistry = wide[wide["experiment_id"].eq(run)]
        crossing = _density_crossing(chemistry)
        event = events[
            events["experiment_id"].astype(str).eq(run)
            & pd.to_numeric(events["yan_added_mg_l"], errors="coerce").notna()
            & (pd.to_numeric(events["relative_time_h"], errors="coerce") > 1e-9)
        ]
        proxy_time = float(event["relative_time_h"].iloc[0]) if not event.empty else np.nan
        if np.isfinite(crossing["pulse_time_h"]):
            model_time = float(crossing["pulse_time_h"])
            source = "chemical_density_linear_interpolation_1040"
        else:
            model_time = proxy_time
            source = "master_schedule_proxy_no_density_crossing"
        rows.append(
            {
                "matrix": LOT_BY_RUN.get(run, "pilot_lot1_qc"),
                "batch": run,
                "pulse_time_h": model_time,
                "schedule_proxy_time_h": proxy_time,
                **crossing,
                "density_target_kg_m3": 1040.0,
                "amount_N_kg_m3": YAN_MODEL_PULSE_KG_M3,
                "mass_audit_N_kg_m3": YAN_MASS_AUDIT_KG_M3,
                "timing_source": source,
                "lot3_event_reconstructed": run in {"26210", "26211", "26212"},
            }
        )
    return pd.DataFrame(rows).sort_values("batch").reset_index(drop=True)


def override_pilot_pulses(
    batches: dict[str, base.BatchData], pulse_schedule: pd.DataFrame
) -> dict[str, base.BatchData]:
    lookup = pulse_schedule.set_index("batch")
    output: dict[str, base.BatchData] = {}
    for run, batch in batches.items():
        pulses = dict(batch.pulses)
        if run in lookup.index and np.isfinite(float(lookup.loc[run, "pulse_time_h"])):
            pulses["N"] = (
                (
                    float(lookup.loc[run, "pulse_time_h"]),
                    float(lookup.loc[run, "amount_N_kg_m3"]),
                ),
            )
        output[run] = replace(batch, pulses=pulses)
    return output


def load_sampling_schedule(primary: pd.DataFrame) -> pd.DataFrame:
    schedule = (
        primary[["experiment_id", "sample_id", "timestamp", "time_h"]]
        .drop_duplicates()
        .rename(columns={"experiment_id": "batch", "time_h": "sample_time_h"})
    )
    schedule["batch"] = schedule["batch"].astype(str)
    return schedule.sort_values(["batch", "sample_time_h"]).reset_index(drop=True)


def load_sensor_normalization() -> pd.DataFrame:
    parameters = pd.read_csv(PILOT_CALIBRATION_DIR / "co2_parameter_estimates.csv")
    values = dict(zip(parameters["parameter"].astype(str), parameters["estimate"].astype(float)))
    background = float(values["co2_background_Ln_min"])
    rows = []
    for reactor in ("TK-31", "TK-32", "TK-33"):
        compact = reactor.replace("-", "")
        rows.append(
            {
                "reactor": reactor,
                "sensor_gain": float(values[f"sensor_gain_{compact}"]),
                "shared_background_ln_min": background,
                "normalization": "max((MassView-background)/sensor_gain, 0)",
                "source": str(
                    (PILOT_CALIBRATION_DIR / "co2_parameter_estimates.csv").relative_to(ROOT_DIR)
                ),
                "caveat": "estimated previously from the same campaign; no pre-inoculation baseline",
            }
        )
    return pd.DataFrame(rows)


def _filter_run(
    group: pd.DataFrame,
    sample_times: np.ndarray,
    pulse_times: np.ndarray,
) -> pd.DataFrame:
    group = group.sort_values("t_h").copy()
    time_h = group["t_h"].to_numpy(dtype=float)
    values = group["co2_rate_normalized_g_l_h"].to_numpy(dtype=float)
    series = pd.Series(values)
    local_median = (
        series.rolling(LOCAL_FILTER_WINDOW_POINTS, center=True, min_periods=5).median()
    )
    absolute_deviation = (series - local_median).abs()
    local_mad = absolute_deviation.rolling(
        LOCAL_FILTER_WINDOW_POINTS, center=True, min_periods=5
    ).median()
    threshold = np.maximum(
        4.0 * 1.4826 * local_mad.fillna(0.0).to_numpy(dtype=float),
        LOCAL_FILTER_MIN_ABS_G_L_H,
    )
    short_outlier = (
        np.isfinite(local_median.to_numpy(dtype=float))
        & (absolute_deviation.to_numpy(dtype=float) > threshold)
    )

    if len(sample_times):
        distance = np.min(np.abs(time_h[:, None] - sample_times[None, :]), axis=1)
    else:
        distance = np.full(len(group), np.nan)
    near_sample = np.isfinite(distance) & (distance <= SAMPLE_WINDOW_H)

    sample_drop = np.zeros(len(group), dtype=bool)
    sample_baseline = np.full(len(group), np.nan)
    for sample_time in sample_times:
        before = (time_h >= sample_time - 2.5) & (time_h <= sample_time - 0.75)
        after = (time_h >= sample_time + 0.75) & (time_h <= sample_time + 2.5)
        center = np.abs(time_h - sample_time) <= 0.75
        if before.sum() < 3 or after.sum() < 3:
            continue
        before_level = float(np.nanmedian(values[before]))
        after_level = float(np.nanmedian(values[after]))
        baseline = 0.5 * (before_level + after_level)
        agreement = abs(before_level - after_level) <= max(0.35 * baseline, 0.05)
        event_has_excursion = np.any(
            center
            & (
                (values < 0.65 * baseline)
                | (values > 1.80 * max(baseline, 0.03))
            )
            & (np.abs(values - baseline) > 0.06)
        )
        if agreement and event_has_excursion:
            sample_drop |= center
            sample_baseline[center] = baseline

    pulse_protected = np.zeros(len(group), dtype=bool)
    for pulse_time in pulse_times:
        pulse_protected |= (time_h > pulse_time) & (time_h <= pulse_time + PULSE_PROTECTION_H)

    artifact = ((short_outlier & near_sample) | sample_drop) & ~pulse_protected
    # Preserve large isolated non-sampling spikes as audit flags but do not
    # replace them automatically; at pilot scale they may be real degassing.
    suspicious_non_sample = short_outlier & ~near_sample & ~pulse_protected
    repaired = values.copy()
    valid = ~artifact & np.isfinite(values)
    if artifact.any() and valid.sum() >= 2:
        repaired[artifact] = np.interp(time_h[artifact], time_h[valid], values[valid])
    baseline_mask = artifact & np.isfinite(sample_baseline)
    repaired[baseline_mask] = sample_baseline[baseline_mask]

    group["local_median_g_l_h"] = local_median.to_numpy(dtype=float)
    group["artifact_flag"] = artifact
    group["suspicious_non_sample_outlier"] = suspicious_non_sample
    group["artifact_reason"] = np.select(
        [sample_drop & artifact, short_outlier & near_sample & artifact],
        ["sampling_event_transient", "near_sample_short_outlier"],
        default="kept",
    )
    group["nearest_sample_distance_h"] = distance
    group["pulse_response_protected"] = pulse_protected
    group["co2_rate_repaired_g_l_h"] = repaired
    return group


def load_co2_observations(
    metadata: pd.DataFrame,
    temperature: pd.DataFrame,
    sampling_schedule: pd.DataFrame,
    pulse_schedule: pd.DataFrame,
    sensor_normalization: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(INTEGRATION_DIR / "co2_minute_qc.csv", low_memory=False)
    raw["experiment_id"] = raw["experiment_id"].astype(str)
    raw["minute"] = pd.to_datetime(raw["minute"])
    raw = raw[raw["experiment_id"].isin(VALID_RUNS)].copy()
    raw = raw[_read_bool(raw["co2_model_include"])].copy()
    observed = pd.to_numeric(raw["mean_flow_ln_min_observed"], errors="coerce")
    raw = raw[observed.notna()].copy()
    raw["observed_flow_ln_min"] = observed[observed.notna()].to_numpy(dtype=float)
    raw = raw[
        pd.to_numeric(raw["artificial_initial_zero"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(raw["artificial_final_zero"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(raw["invalid_intermediate_zero"], errors="coerce").fillna(0).eq(0)
        & pd.to_numeric(raw["postprocess_cooling"], errors="coerce").fillna(0).eq(0)
    ].copy()

    metadata = metadata.copy()
    metadata["experiment_id"] = metadata["experiment_id"].astype(str)
    meta_columns = [
        "experiment_id",
        "sampling_start",
        "active_end",
        "lot",
        "reactor",
        "condition",
        "replicate",
        "protocol",
        "initial_volume_l",
    ]
    raw = raw.merge(metadata[meta_columns], on="experiment_id", how="left", validate="many_to_one")
    raw["sampling_start"] = pd.to_datetime(raw["sampling_start"])
    raw["active_end"] = pd.to_datetime(raw["active_end"])
    raw["t_h"] = (raw["minute"] - raw["sampling_start"]).dt.total_seconds() / 3600.0
    raw = raw[(raw["t_h"] >= 0.0) & (raw["minute"] <= raw["active_end"])].copy()
    raw["batch"] = raw["experiment_id"]
    raw["matrix"] = raw["batch"].map(LOT_BY_RUN)

    normalization = sensor_normalization.set_index("reactor")
    raw["sensor_gain"] = raw["reactor"].map(normalization["sensor_gain"])
    raw["sensor_background_ln_min"] = raw["reactor"].map(
        normalization["shared_background_ln_min"]
    )
    raw["co2_rate_uncorrected_g_l_h"] = raw["observed_flow_ln_min"] * LN_MIN_TO_G_L_H
    raw["co2_rate_background_corrected_ln_min"] = np.maximum(
        raw["observed_flow_ln_min"] - raw["sensor_background_ln_min"], 0.0
    )
    raw["co2_rate_normalized_g_l_h"] = (
        raw["co2_rate_background_corrected_ln_min"]
        / raw["sensor_gain"]
        * LN_MIN_TO_G_L_H
    )

    temp = temperature.copy()
    temp["experiment_id"] = temp["experiment_id"].astype(str)
    temp = temp.sort_values(["experiment_id", "time_h"])
    raw["sensor_temperature_c"] = np.nan
    raw["setpoint_c"] = np.nan
    for run, group in raw.groupby("batch"):
        source = temp[temp["experiment_id"].eq(run)].dropna(
            subset=["time_h", "executed_temperature_c"]
        )
        index = group.index
        raw.loc[index, "sensor_temperature_c"] = np.interp(
            group["t_h"].to_numpy(dtype=float),
            source["time_h"].to_numpy(dtype=float),
            source["executed_temperature_c"].to_numpy(dtype=float),
        )
        setpoint_source = source.dropna(subset=["commanded_setpoint_c"])
        raw.loc[index, "setpoint_c"] = np.interp(
            group["t_h"].to_numpy(dtype=float),
            setpoint_source["time_h"].to_numpy(dtype=float),
            setpoint_source["commanded_setpoint_c"].to_numpy(dtype=float),
        )

    filtered_frames = []
    for run, group in raw.groupby("batch", sort=True):
        sample_times = sampling_schedule.loc[
            sampling_schedule["batch"].eq(run), "sample_time_h"
        ].to_numpy(dtype=float)
        pulse_times = pulse_schedule.loc[
            pulse_schedule["batch"].eq(run), "pulse_time_h"
        ].dropna().to_numpy(dtype=float)
        filtered_frames.append(_filter_run(group, sample_times, pulse_times))
    sensor_qc = pd.concat(filtered_frames, ignore_index=True).sort_values(
        ["matrix", "batch", "t_h"]
    )

    hourly_rows: list[pd.DataFrame] = []
    for run, group in sensor_qc.groupby("batch", sort=True):
        group = group.sort_values("t_h").copy()
        group["hour_bin"] = np.floor(group["t_h"]).astype(int)
        hourly = group.groupby("hour_bin", sort=True).agg(
            t_h=("t_h", "mean"),
            co2_rate_hourly_raw_g_l_h=("co2_rate_uncorrected_g_l_h", "median"),
            co2_rate_hourly_normalized_g_l_h=("co2_rate_normalized_g_l_h", "median"),
            co2_rate_hourly_repaired_g_l_h=("co2_rate_repaired_g_l_h", "median"),
            sensor_temperature_c=("sensor_temperature_c", "median"),
            setpoint_c=("setpoint_c", "median"),
            artifact_fraction=("artifact_flag", "mean"),
            suspicious_outlier_fraction=("suspicious_non_sample_outlier", "mean"),
            source_minutes=("observed_flow_ln_min", "size"),
        ).reset_index(drop=True)
        pulse_times = pulse_schedule.loc[
            pulse_schedule["batch"].eq(run), "pulse_time_h"
        ].dropna().to_numpy(dtype=float)
        segment = np.searchsorted(pulse_times, hourly["t_h"].to_numpy(dtype=float), side="right")
        smoothed = np.empty(len(hourly), dtype=float)
        for segment_id in np.unique(segment):
            mask = segment == segment_id
            repaired = hourly.loc[mask, "co2_rate_hourly_repaired_g_l_h"].to_numpy(dtype=float)
            robust = (
                pd.Series(repaired).rolling(3, center=True, min_periods=1).median().to_numpy(dtype=float)
            )
            if len(robust) >= 5:
                robust = savgol_filter(robust, 5, 2, mode="interp")
            smoothed[mask] = np.maximum(robust, 0.0)
        hourly["co2_rate_g_l_h"] = smoothed
        hourly["smoothing_segment"] = segment
        hourly["batch"] = run
        hourly["experiment_code"] = run
        hourly["matrix"] = LOT_BY_RUN[run]
        hourly["lot"] = str(group["lot"].iloc[0])
        hourly["reactor"] = str(group["reactor"].iloc[0])
        hourly["protocol"] = str(group["protocol"].iloc[0])
        hourly["calibratable"] = True
        hourly["chemistry_first_h"] = 0.0
        hourly["chemistry_last_h"] = float(
            (pd.Timestamp(group["active_end"].iloc[0]) - pd.Timestamp(group["sampling_start"].iloc[0])).total_seconds()
            / 3600.0
        )
        hourly["detection_limit_g_l_h"] = np.where(
            hourly["sensor_temperature_c"] <= COLD_THRESHOLD_C,
            DETECTION_LIMIT_COLD_G_L_H,
            DETECTION_LIMIT_STANDARD_G_L_H,
        )
        hourly["left_censored"] = hourly["co2_rate_g_l_h"] <= hourly["detection_limit_g_l_h"]
        hourly_rows.append(hourly)
    observations = pd.concat(hourly_rows, ignore_index=True).sort_values(
        ["matrix", "batch", "t_h"]
    ).reset_index(drop=True)

    unnormalized = observations.copy()
    unnormalized["co2_rate_g_l_h"] = unnormalized["co2_rate_hourly_raw_g_l_h"]
    unnormalized["left_censored"] = (
        unnormalized["co2_rate_g_l_h"] <= unnormalized["detection_limit_g_l_h"]
    )

    audit_rows = []
    for run, group in raw.groupby("batch", sort=True):
        audit_rows.append(
            {
                "matrix": LOT_BY_RUN[run],
                "batch": run,
                "reactor": str(group["reactor"].iloc[0]),
                "first_observed_time_h": float(group["t_h"].min()),
                "last_observed_time_h": float(group["t_h"].max()),
                "n_observed_minutes": int(len(group)),
                "sensor_gain": float(group["sensor_gain"].iloc[0]),
                "shared_background_ln_min": float(group["sensor_background_ln_min"].iloc[0]),
                "background_g_l_h_before_gain": float(
                    group["sensor_background_ln_min"].iloc[0] * LN_MIN_TO_G_L_H
                ),
                "pre_inoculation_baseline_available": False,
                "artificial_leading_zeros_excluded": True,
            }
        )
    return observations, unnormalized, sensor_qc, pd.DataFrame(audit_rows)


def chemistry_support_table(batches: dict[str, base.BatchData]) -> pd.DataFrame:
    rows = []
    for run in VALID_RUNS:
        batch = batches[run]
        observed = np.zeros(len(batch.time), dtype=bool)
        for state in base.STATE_NAMES:
            values = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)))
            observed |= np.isfinite(values.astype(float))
        times = batch.time[observed]
        rows.append(
            {
                "matrix": LOT_BY_RUN[run],
                "batch": run,
                "chemistry_first_h": float(np.min(times)),
                "chemistry_last_h": float(np.max(times)),
                "n_chemical_timepoints": int(len(np.unique(times))),
            }
        )
    return pd.DataFrame(rows)


def ph_do_availability(primary: pd.DataFrame) -> pd.DataFrame:
    wide = _primary_wide(primary)
    columns = ["experiment_id", "time_h", "pH", "dissolved_oxygen"]
    for column in columns:
        if column not in wide:
            wide[column] = np.nan
    result = wide[columns].rename(columns={"experiment_id": "batch"}).copy()
    result["matrix"] = result["batch"].map(LOT_BY_RUN)
    result["pH_valid"] = result["pH"].between(2.0, 5.0) | result["pH"].isna()
    result["do_valid"] = result["dissolved_oxygen"].ge(0.0) | result["dissolved_oxygen"].isna()
    result["model_use"] = False
    result["reason"] = "sparse discrete measurements; first values occur at 16 h and cannot resolve initial onset"
    return result[result["batch"].isin(VALID_RUNS)].reset_index(drop=True)


def fit_explicit_batches(
    label: str,
    batch_names: list[str],
    observations: pd.DataFrame,
    cache: dict[str, lab_co2.DriverCache],
    *,
    n_starts: int,
    max_nfev: int,
    seed: int,
    model_variant: str = MODEL_NAME,
    initial_parameters: dict[str, float] | None = None,
) -> tuple[dict[str, float], pd.DataFrame]:
    names = lab_co2._shape_parameter_names(model_variant)
    lower = np.log([lab_co2.PARAMETER_BOUNDS[name][0] for name in names])
    upper = np.log([lab_co2.PARAMETER_BOUNDS[name][1] for name in names])
    reference_values = {
        "kCO2_release_h": 1.2,
        "CO2sat_scale": 0.50,
        "O2_qmax_mg_gdw_h": 0.60,
        "O2_initial_scale": 0.75,
        "pulse_t_rise_h": 12.0,
        "pulse_activity_gain": 1.25,
        "chem_activation_start_fraction": 0.40,
        "chem_activation_duration_fraction": 0.80,
    }
    if initial_parameters:
        reference_values.update(
            {name: float(initial_parameters[name]) for name in names if name in initial_parameters}
        )
    reference = np.log([reference_values[name] for name in names])
    rng = np.random.default_rng(seed)
    starts = [reference] + [
        rng.uniform(lower + 1e-6, upper - 1e-6) for _ in range(max(n_starts - 1, 0))
    ]
    rows = []
    best_result = None
    best_wsse = np.inf
    for start_index, start in enumerate(starts):
        result = least_squares(
            lambda value: lab_co2._fit_residual(
                value, batch_names, observations, cache, model_variant
            ),
            np.clip(start, lower + 1e-9, upper - 1e-9),
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            max_nfev=max_nfev,
        )
        residual = lab_co2._fit_residual(
            result.x, batch_names, observations, cache, model_variant
        )
        wsse = float(np.dot(residual, residual))
        gain, _ = lab_co2._profile_matrix_gain(
            result.x, batch_names, observations, cache, model_variant
        )
        values = lab_co2._shape_values(result.x, model_variant)
        rows.append(
            {
                "matrix": label,
                "model": model_variant,
                "start": start_index,
                "success": bool(result.success),
                "status": int(result.status),
                "nfev": int(result.nfev),
                "wsse_equal_batch": wsse,
                "matrix_gain": gain,
                "message": str(result.message),
                **values,
            }
        )
        if wsse < best_wsse:
            best_wsse = wsse
            best_result = result
    if best_result is None:
        raise RuntimeError(f"No optimizer result for {label}")
    values = lab_co2._shape_values(best_result.x, model_variant)
    gain, _ = lab_co2._profile_matrix_gain(
        best_result.x, batch_names, observations, cache, model_variant
    )
    parameters = {
        **values,
        "matrix_gain": float(gain),
        "model": model_variant,
        "wsse_equal_batch": best_wsse,
        "n_calibration_batches": len(batch_names),
        "n_calibration_points": int(observations["batch"].isin(batch_names).sum()),
    }
    return parameters, pd.DataFrame(rows)


def score_parameter_set(
    label: str,
    parameters: dict[str, float],
    observations: pd.DataFrame,
    cache: dict[str, lab_co2.DriverCache],
    roles: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    prediction_rows = []
    roles = roles or {}
    for run, group in observations[observations["calibratable"]].groupby("batch", sort=True):
        raw_prediction = lab_co2.raw_qgas_prediction(
            cache[run],
            parameters["kCO2_release_h"],
            parameters["CO2sat_scale"],
            parameters["O2_qmax_mg_gdw_h"],
            parameters["O2_initial_scale"],
            chemistry_aligned=True,
            nitrogen_boost_transition=True,
            pulse_t_rise_h=parameters["pulse_t_rise_h"],
            pulse_activity_gain=parameters["pulse_activity_gain"],
            continuous_release=True,
            bounded_chemical_activation=True,
            chem_activation_start_fraction=parameters["chem_activation_start_fraction"],
            chem_activation_duration_fraction=parameters["chem_activation_duration_fraction"],
        )
        predicted = float(parameters["matrix_gain"]) * raw_prediction
        observed = group["co2_rate_g_l_h"].to_numpy(dtype=float)
        time_h = group["t_h"].to_numpy(dtype=float)
        metric = lab_co2._metric_row(
            label,
            str(group["matrix"].iloc[0]),
            run,
            observed,
            predicted,
            time_h,
            roles.get(run, "external_transfer"),
            group["left_censored"].to_numpy(dtype=bool),
            group["detection_limit_g_l_h"].to_numpy(dtype=float),
            MODEL_NAME,
            cache[run].n_pulse_time_h,
        )
        metric_rows.append(metric)
        for t, obs, pred, censored in zip(
            time_h, observed, predicted, group["left_censored"].to_numpy(dtype=bool)
        ):
            prediction_rows.append(
                {
                    "calibration_source": label,
                    "target_matrix": str(group["matrix"].iloc[0]),
                    "batch": run,
                    "role": roles.get(run, "external_transfer"),
                    "time_h": float(t),
                    "observed_g_l_h": float(obs),
                    "predicted_g_l_h": float(pred),
                    "left_censored": bool(censored),
                    "pulse_time_h": cache[run].n_pulse_time_h,
                }
            )
    return pd.DataFrame(prediction_rows), pd.DataFrame(metric_rows)


def _load_lab_fits() -> dict[str, dict[str, float]]:
    table = pd.read_csv(LAB_RESULTS_DIR / "fit_parameters.csv")
    fits: dict[str, dict[str, float]] = {}
    for matrix in ("natural", "synthetic"):
        selected = table[table["calibration_matrix"].eq(matrix)]
        parameters = dict(
            zip(selected["parameter"].astype(str), selected["estimate"].astype(float))
        )
        parameters["model"] = MODEL_NAME
        fits[matrix] = parameters
    return fits


def cross_scale_transfer(
    pooled_fit: dict[str, float],
    pilot_observations: pd.DataFrame,
    pilot_cache: dict[str, lab_co2.DriverCache],
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    lab_fits = _load_lab_fits()
    pilot_roles = {
        run: ("pilot_holdout" if run in HOLDOUTS.values() else "pilot_context")
        for run in VALID_RUNS
    }
    prediction_frames = []
    metric_frames = []
    for matrix, lab_fit in lab_fits.items():
        source = f"laboratory_{matrix}_fit"
        direction = f"laboratory_{matrix}_to_pilot"
        predictions, metrics = score_parameter_set(
            source, lab_fit, pilot_observations, pilot_cache, pilot_roles
        )
        metrics["transfer_direction"] = direction
        predictions["transfer_direction"] = direction
        prediction_frames.append(predictions)
        metric_frames.append(metrics)

    try:
        lab_batches, lab_theta, _ = lab_co2.load_batches()
        lab_observations = pd.read_csv(LAB_RESULTS_DIR / "co2_observations_hourly.csv")
        lab_observations = lab_observations[
            (
                lab_observations["matrix"].eq("natural")
                & lab_observations["batch"].eq("LAB012")
            )
            | (
                lab_observations["matrix"].eq("synthetic")
                & lab_observations["batch"].eq("lot2_F3")
            )
        ].copy()
        lab_cache, _ = lab_co2.build_driver_cache(
            lab_batches,
            {matrix: lab_theta[matrix] for matrix in ("natural", "synthetic")},
            lab_observations,
        )
        lab_predictions, lab_metrics = score_parameter_set(
            "pilot_pooled_fit",
            pooled_fit,
            lab_observations,
            lab_cache,
            {"LAB012": "laboratory_holdout", "lot2_F3": "laboratory_holdout"},
        )
        lab_metrics["transfer_direction"] = (
            "pilot_to_laboratory_" + lab_metrics["target_matrix"].astype(str)
        )
        lab_predictions["transfer_direction"] = (
            "pilot_to_laboratory_" + lab_predictions["target_matrix"].astype(str)
        )
        prediction_frames.append(lab_predictions)
        metric_frames.append(lab_metrics)
        status = "both_directions_completed"
        return pd.concat(prediction_frames, ignore_index=True), pd.concat(
            metric_frames, ignore_index=True
        ), status
    except Exception as error:  # optional external lab raw-data dependency
        status = f"pilot_to_laboratory_not_run: {type(error).__name__}: {error}"
        return pd.concat(prediction_frames, ignore_index=True), pd.concat(
            metric_frames, ignore_index=True
        ), status


def _fit_table(fits: dict[str, dict[str, float]]) -> pd.DataFrame:
    return lab_co2.fit_parameter_table(fits)


def _figure_path(name: str) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    return FIGURE_DIR / name


def _finalize_figure(fig: plt.Figure, name: str) -> plt.Figure:
    path = _figure_path(name)
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor="white")
    return fig


def plot_data_overview(observations: pd.DataFrame, pulse_schedule: pd.DataFrame) -> plt.Figure:
    runs = list(VALID_RUNS)
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=False, sharey=True)
    for ax, run in zip(axes.flat, runs):
        group = observations[observations["batch"].eq(run)]
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(run), "pulse_time_h"]
        ax.plot(group["t_h"], group["co2_rate_g_l_h"], color=COLORS["blue"], lw=1.8)
        ax.fill_between(
            group["t_h"], 0, group["detection_limit_g_l_h"], color=COLORS["light_grey"], alpha=0.45
        )
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls="--", lw=1.4)
        ax.set_title(f"{run} · {group['reactor'].iloc[0]} · protocolo {group['protocol'].iloc[0]}")
        ax.set_xlabel("Tiempo desde primera muestra química (h)")
        ax.set_ylabel("CO₂ (g L⁻¹ h⁻¹)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    fig.suptitle("Piloto 2026: perfiles CO₂ normalizados y filtrados", fontsize=15)
    fig.tight_layout()
    return _finalize_figure(fig, "co2_data_overview.png")


def plot_sensor_filter_examples(sensor_qc: pd.DataFrame, pulse_schedule: pd.DataFrame) -> plt.Figure:
    examples = ("26158", "26211", "26212")
    fig, axes = plt.subplots(len(examples), 1, figsize=(13, 9), sharex=False)
    for ax, run in zip(axes, examples):
        group = sensor_qc[sensor_qc["batch"].eq(run)].sort_values("t_h")
        ax.plot(group["t_h"], group["co2_rate_uncorrected_g_l_h"], color=COLORS["grey"], alpha=0.45, lw=0.8, label="MassView convertido")
        ax.plot(group["t_h"], group["co2_rate_normalized_g_l_h"], color=COLORS["blue_light"], lw=0.9, label="Normalizado")
        ax.plot(group["t_h"], group["co2_rate_repaired_g_l_h"], color=COLORS["blue"], lw=1.2, label="Artefactos reparados")
        flagged = group[group["artifact_flag"]]
        ax.scatter(flagged["t_h"], flagged["co2_rate_normalized_g_l_h"], s=18, facecolors="none", edgecolors=COLORS["pink"], label="Reemplazado")
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(run), "pulse_time_h"]
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls="--", lw=1.4, label="Pulso N")
        ax.set_title(f"{run}: normalización y filtro a 1 min")
        ax.set_ylabel("CO₂ (g L⁻¹ h⁻¹)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    axes[-1].set_xlabel("Tiempo (h)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return _finalize_figure(fig, "sensor_normalization_and_filter.png")


def plot_temperature_profiles(
    temperature: pd.DataFrame, pulse_schedule: pd.DataFrame
) -> plt.Figure:
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharey=True)
    for ax, run in zip(axes.flat, VALID_RUNS):
        group = temperature[temperature["experiment_id"].astype(str).eq(run)].sort_values("time_h")
        ax.plot(group["time_h"], group["executed_temperature_c"], color=COLORS["orange"], lw=1.5, label="Sonda1")
        ax.step(group["time_h"], group["commanded_setpoint_c"], where="post", color=COLORS["ink"], lw=1.0, ls="--", label="Setpoint")
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(run), "pulse_time_h"]
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls=":", lw=1.5)
        ax.set_title(run)
        ax.set_xlabel("Tiempo (h)")
        ax.set_ylabel("Temperatura (°C)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles.append(Line2D([0], [0], color=COLORS["gold"], ls=":", label="Pulso N"))
    fig.suptitle("Temperatura ejecutada utilizada por el modelo", fontsize=15, y=0.995)
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    return _finalize_figure(fig, "temperature_profiles_model_input.png")


def plot_chemistry_timeline(
    primary: pd.DataFrame, pulse_schedule: pd.DataFrame
) -> plt.Figure:
    wide = _primary_wide(primary)
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharey=False)
    for ax, run in zip(axes.flat, VALID_RUNS):
        group = wide[wide["experiment_id"].eq(run)]
        sugar = group.get("glucose", pd.Series(np.nan, index=group.index)) + group.get("fructose", pd.Series(np.nan, index=group.index))
        ax.plot(group["time_h"], sugar, marker="o", ms=3.5, color=COLORS["blue"], label="G+F")
        ax2 = ax.twinx()
        ax2.plot(group["time_h"], group.get("density"), marker="s", ms=3, color=COLORS["orange"], label="Densidad")
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(run), "pulse_time_h"]
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls="--", lw=1.4)
        ax.set_title(run)
        ax.set_xlabel("Tiempo (h)")
        ax.set_ylabel("Glucosa + fructosa (g L⁻¹)")
        ax2.set_ylabel("Densidad (kg m⁻³)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.6)
    fig.suptitle("Soporte químico y cruce de densidad para el pulso", fontsize=15)
    fig.tight_layout()
    return _finalize_figure(fig, "chemistry_and_pulse_timeline.png")


def plot_calibration_overlays(
    predictions: pd.DataFrame, pulse_schedule: pd.DataFrame, matrix: str
) -> plt.Figure:
    runs = [run for run in VALID_RUNS if LOT_BY_RUN[run] == matrix]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    for ax, run in zip(axes, runs):
        group = predictions[
            predictions["calibration_matrix"].eq(matrix) & predictions["batch"].eq(run)
        ].sort_values("time_h")
        ax.plot(group["time_h"], group["observed_g_l_h"], color=COLORS["ink"], lw=1.4, label="Observado")
        ax.plot(group["time_h"], group["predicted_g_l_h"], color=COLORS["blue"], lw=1.7, label="Modelo")
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(run), "pulse_time_h"]
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls="--", lw=1.3)
        ax.set_title(f"{run}{' · holdout' if run == HOLDOUTS[matrix] else ''}")
        ax.set_xlabel("Tiempo (h)")
        ax.set_ylabel("CO₂ (g L⁻¹ h⁻¹)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle(f"Ajuste y validación interna · {matrix}", fontsize=14, y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    return _finalize_figure(fig, f"calibration_overlays_{matrix}.png")


def plot_cross_lot_validation(
    predictions: pd.DataFrame, validation: pd.DataFrame, pulse_schedule: pd.DataFrame
) -> plt.Figure:
    scenarios = validation.sort_values(["target_matrix", "role"])
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharey=True)
    for ax, row in zip(axes.flat, scenarios.itertuples(index=False)):
        group = predictions[
            predictions["calibration_matrix"].eq(row.calibration_matrix)
            & predictions["batch"].eq(row.batch)
        ].sort_values("time_h")
        ax.plot(group["time_h"], group["observed_g_l_h"], color=COLORS["ink"], lw=1.4, label="Observado")
        ax.plot(group["time_h"], group["predicted_g_l_h"], color=COLORS["blue"], lw=1.7, label="Predicho")
        pulse = pulse_schedule.loc[pulse_schedule["batch"].eq(row.batch), "pulse_time_h"]
        if len(pulse):
            ax.axvline(float(pulse.iloc[0]), color=COLORS["gold"], ls="--", lw=1.3)
        ax.set_title(
            f"{row.calibration_matrix} → {row.batch}\nRMSE={row.rmse_g_l_h:.3f}; Δonset={row.onset_delay_h:.1f} h"
        )
        ax.set_xlabel("Tiempo (h)")
        ax.set_ylabel("CO₂ (g L⁻¹ h⁻¹)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Validación interna y transferencia cruzada entre lotes", fontsize=15, y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=2,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    return _finalize_figure(fig, "heldout_cross_lot_validation.png")


def plot_model_comparison(model_validation: pd.DataFrame) -> plt.Figure:
    table = model_validation.copy()
    model_order = [MODEL_NAME, THRESHOLD_MODEL_NAME, LEGACY_MODEL_NAME]
    model_labels = {
        MODEL_NAME: "Liberación continua",
        THRESHOLD_MODEL_NAME: "Umbral + pulso N",
        LEGACY_MODEL_NAME: "Transición lenta",
    }
    model_colors = {
        MODEL_NAME: COLORS["blue"],
        THRESHOLD_MODEL_NAME: COLORS["gold"],
        LEGACY_MODEL_NAME: COLORS["grey"],
    }
    scenario_order = (
        table[["calibration_matrix", "batch", "role"]]
        .drop_duplicates()
        .sort_values(["batch", "role"])
        .reset_index(drop=True)
    )
    scenario_labels = [
        f"{row.calibration_matrix.replace('pilot_', '')} → {row.batch}\n"
        f"{'interno' if row.role == 'internal_holdout' else 'cruzado'}"
        for row in scenario_order.itertuples(index=False)
    ]
    fig, ax = plt.subplots(figsize=(11.5, 5.5))
    centers = np.arange(len(scenario_order), dtype=float)
    width = 0.24
    for offset, model in enumerate(model_order):
        values = []
        for scenario in scenario_order.itertuples(index=False):
            selected = table[
                table["model"].eq(model)
                & table["calibration_matrix"].eq(scenario.calibration_matrix)
                & table["batch"].eq(scenario.batch)
                & table["role"].eq(scenario.role)
            ]
            values.append(float(selected["rmse_g_l_h"].iloc[0]))
        ax.bar(
            centers + (offset - 1) * width,
            values,
            width,
            color=model_colors[model],
            label=model_labels[model],
        )
    ax.set_xticks(centers, scenario_labels)
    ax.set_ylabel("RMSE (g L⁻¹ h⁻¹)")
    ax.set_title("Comparación de modelos en los holdouts piloto")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.grid(axis="y", color=COLORS["grid"], lw=0.6, alpha=0.7)
    fig.tight_layout()
    return _finalize_figure(fig, "heldout_model_comparison.png")


def plot_parameter_comparison(parameter_table: pd.DataFrame) -> plt.Figure:
    pivot = parameter_table.pivot(index="parameter", columns="calibration_matrix", values="estimate")
    order = list(lab_co2._shape_parameter_names(MODEL_NAME)) + ["matrix_gain"]
    pivot = pivot.reindex(order)
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(pivot))
    width = 0.36
    ax.bar(x - width / 2, pivot[MATRICES[0]], width, color=COLORS["blue"], label=MATRICES[0])
    ax.bar(x + width / 2, pivot[MATRICES[1]], width, color=COLORS["gold"], label=MATRICES[1])
    ax.set_yscale("log")
    ax.set_xticks(x, pivot.index, rotation=45, ha="right")
    ax.set_ylabel("Estimación, escala log")
    ax.set_title("Parámetros del bloque CO₂ por lote")
    ax.legend(frameon=False)
    ax.grid(axis="y", color=COLORS["grid"], lw=0.6, alpha=0.7)
    fig.tight_layout()
    return _finalize_figure(fig, "parameter_comparison.png")


def plot_identifiability_profiles(profiles: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, ((matrix, parameter), group) in zip(
        axes.flat, profiles.groupby(["matrix", "profiled_parameter"], sort=True)
    ):
        ax.plot(group["fixed_value"], group["relative_wsse_increase"], marker="o", color=COLORS["blue"])
        ax.axhline(0.05, color=COLORS["gold"], ls="--", lw=1.2)
        ax.axvline(float(group["full_fit_estimate"].iloc[0]), color=COLORS["ink"], ls=":", lw=1.2)
        ax.set_title(f"{matrix} · {parameter}")
        ax.set_xlabel("Valor fijado")
        ax.set_ylabel("Aumento relativo WSSE")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    fig.suptitle("Perfiles prácticos de identificabilidad", fontsize=15)
    fig.tight_layout()
    return _finalize_figure(fig, "activation_identifiability_profiles.png")


def plot_cross_scale_transfer(predictions: pd.DataFrame, metrics: pd.DataFrame) -> plt.Figure:
    selected = metrics[
        metrics["role"].isin(["pilot_holdout", "laboratory_holdout"])
    ].copy()
    n = max(len(selected), 1)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.4 * n), squeeze=False)
    for ax, row in zip(axes.flat, selected.itertuples(index=False)):
        group = predictions[
            predictions["calibration_source"].eq(row.calibration_matrix)
            & predictions["batch"].eq(row.batch)
        ].sort_values("time_h")
        ax.plot(group["time_h"], group["observed_g_l_h"], color=COLORS["ink"], lw=1.3, label="Observado")
        ax.plot(group["time_h"], group["predicted_g_l_h"], color=COLORS["orange"], lw=1.6, label="Transferido")
        if np.isfinite(float(group["pulse_time_h"].iloc[0])):
            ax.axvline(float(group["pulse_time_h"].iloc[0]), color=COLORS["gold"], ls="--", lw=1.3)
        ax.set_title(f"{row.transfer_direction} · {row.batch} · RMSE={row.rmse_g_l_h:.3f}")
        ax.set_xlabel("Tiempo (h)")
        ax.set_ylabel("CO₂ (g L⁻¹ h⁻¹)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    axes.flat[0].legend(frameon=False)
    fig.suptitle(
        "Transferencia del bloque CO₂ entre laboratorio y piloto",
        fontsize=15,
        y=0.998,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.5)
    return _finalize_figure(fig, "cross_scale_transfer.png")


def plot_ph_do_availability(ph_do: pd.DataFrame) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for run, group in ph_do.groupby("batch"):
        valid_ph = group[group["pH_valid"] & group["pH"].notna()]
        if not valid_ph.empty:
            axes[0].plot(valid_ph["time_h"], valid_ph["pH"], marker="o", ms=3, label=run)
        valid_do = group[group["do_valid"] & group["dissolved_oxygen"].notna()]
        if not valid_do.empty:
            axes[1].plot(
                valid_do["time_h"],
                valid_do["dissolved_oxygen"],
                marker="o",
                ms=3,
                label=run,
            )
    axes[0].set_title("pH discreto disponible")
    axes[1].set_title("O₂ disuelto discreto disponible")
    axes[0].set_ylabel("pH")
    axes[1].set_ylabel("O₂ disuelto (mg L⁻¹)")
    for ax in axes:
        ax.set_xlabel("Tiempo (h)")
        ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
        ax.legend(frameon=False)
    fig.suptitle("Cobertura auxiliar: sólo lote 3 y desde 16 h", fontsize=14)
    fig.tight_layout()
    return _finalize_figure(fig, "ph_do_availability.png")


def run_analysis(
    n_starts: int = 5,
    max_nfev: int = 300,
    seed: int = RANDOM_SEED,
) -> dict[str, object]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    lab_co2.HOLDOUTS = dict(HOLDOUTS)
    lab_co2.SYNTHETIC_CODES = {run: run for run in VALID_RUNS}

    batches, pilot_theta, tables = _pilot_upstream_context()
    pulse_schedule = build_pulse_schedule(tables.primary, tables.events)
    batches = override_pilot_pulses(batches, pulse_schedule)
    sampling_schedule = load_sampling_schedule(tables.primary)
    sensor_normalization = load_sensor_normalization()
    observations, unnormalized_observations, sensor_qc, sensor_audit = load_co2_observations(
        tables.metadata,
        tables.temperature,
        sampling_schedule,
        pulse_schedule,
        sensor_normalization,
    )
    chemistry_support = chemistry_support_table(batches)
    ph_do = ph_do_availability(tables.primary)

    theta_by_matrix = {matrix: dict(pilot_theta) for matrix in MATRICES}
    cache, driver_diagnostics = lab_co2.build_driver_cache(
        batches, theta_by_matrix, observations
    )
    temperature_inputs = lab_co2.temperature_model_input_table(cache)
    temperature_alignment = lab_co2.temperature_alignment_summary(
        observations, temperature_inputs
    )

    fits: dict[str, dict[str, float]] = {}
    fit_starts = []
    for index, matrix in enumerate(MATRICES):
        parameters, starts = lab_co2.fit_matrix(
            matrix,
            observations,
            cache,
            n_starts=n_starts,
            max_nfev=max_nfev,
            seed=seed + index * 100,
            model_variant=MODEL_NAME,
        )
        fits[matrix] = parameters
        fit_starts.append(starts)
    fit_starts_table = pd.concat(fit_starts, ignore_index=True)
    fit_parameters = _fit_table(fits)
    predictions, batch_metrics, validation = lab_co2.predict_and_score(
        fits, observations, cache, MODEL_NAME
    )

    jacobian, singular_values, local_correlations = lab_co2.jacobian_identifiability_diagnostics(
        fits, observations, cache
    )
    activation_profiles = lab_co2.profile_activation_parameters(fits, observations, cache)
    loo_estimates, loo_summary = lab_co2.leave_one_batch_out_stability(
        fits, observations, cache, seed
    )

    comparator_fits: dict[str, dict[str, dict[str, float]]] = {}
    comparator_starts = []
    comparator_predictions = []
    comparator_metrics = []
    comparator_validation = []
    for model_variant in (THRESHOLD_MODEL_NAME, LEGACY_MODEL_NAME):
        model_fits = {}
        for index, matrix in enumerate(MATRICES):
            parameters, starts = lab_co2.fit_matrix(
                matrix,
                observations,
                cache,
                n_starts=max(3, n_starts // 2),
                max_nfev=max(180, max_nfev // 2),
                seed=seed + 500 + index * 100,
                model_variant=model_variant,
            )
            model_fits[matrix] = parameters
            comparator_starts.append(starts)
        comparator_fits[model_variant] = model_fits
        pred, metrics, val = lab_co2.predict_and_score(
            model_fits, observations, cache, model_variant
        )
        comparator_predictions.append(pred)
        comparator_metrics.append(metrics)
        comparator_validation.append(val)
    all_model_validation = pd.concat(
        [validation] + comparator_validation, ignore_index=True
    ).sort_values(["batch", "calibration_matrix", "model"])

    # Sensitivity to the fixed sensor normalization.  The driver cache has the
    # same hourly time grid, so only the observation vector changes.
    unnormalized_fits = {}
    unnormalized_starts = []
    for index, matrix in enumerate(MATRICES):
        parameters, starts = lab_co2.fit_matrix(
            matrix,
            unnormalized_observations,
            cache,
            n_starts=2,
            max_nfev=max(180, max_nfev // 2),
            seed=seed + 900 + index * 100,
            model_variant=MODEL_NAME,
            initial_parameters=fits[matrix],
        )
        unnormalized_fits[matrix] = parameters
        unnormalized_starts.append(starts)
    _, _, unnormalized_validation = lab_co2.predict_and_score(
        unnormalized_fits, unnormalized_observations, cache, MODEL_NAME
    )
    sensor_normalization_sensitivity = validation[
        ["calibration_matrix", "target_matrix", "batch", "role", "rmse_g_l_h", "onset_delay_h"]
    ].merge(
        unnormalized_validation[
            ["calibration_matrix", "target_matrix", "batch", "role", "rmse_g_l_h", "onset_delay_h"]
        ],
        on=["calibration_matrix", "target_matrix", "batch", "role"],
        suffixes=("_normalized", "_raw_massview"),
    )

    pooled_calibration_batches = [
        run for run in VALID_RUNS if run not in set(HOLDOUTS.values())
    ]
    pooled_fit, pooled_starts = fit_explicit_batches(
        "pilot_pooled",
        pooled_calibration_batches,
        observations,
        cache,
        n_starts=n_starts,
        max_nfev=max_nfev,
        seed=seed + 1300,
        model_variant=MODEL_NAME,
    )
    pooled_predictions, pooled_metrics = score_parameter_set(
        "pilot_pooled_fit",
        pooled_fit,
        observations,
        cache,
        {run: ("pooled_holdout" if run in HOLDOUTS.values() else "pooled_calibration") for run in VALID_RUNS},
    )

    cross_scale_predictions, cross_scale_metrics, cross_scale_status = cross_scale_transfer(
        pooled_fit, observations, cache
    )

    inventory = (
        observations.groupby(["matrix", "batch", "lot", "reactor", "protocol"], as_index=False)
        .agg(
            co2_hourly_rows=("t_h", "size"),
            co2_first_h=("t_h", "min"),
            co2_last_h=("t_h", "max"),
            observed_peak_g_l_h=("co2_rate_g_l_h", "max"),
            artifact_fraction=("artifact_fraction", "mean"),
            censored_fraction=("left_censored", "mean"),
        )
        .merge(chemistry_support, on=["matrix", "batch"], how="left")
    )
    inventory["role"] = np.where(
        inventory.apply(lambda row: HOLDOUTS[row["matrix"]] == row["batch"], axis=1),
        "holdout",
        "calibration",
    )

    fit_parameters.to_csv(RESULTS_DIR / "fit_parameters.csv", index=False)
    fit_starts_table.to_csv(RESULTS_DIR / "fit_start_diagnostics.csv", index=False)
    observations.to_csv(RESULTS_DIR / "co2_observations_hourly.csv", index=False)
    sensor_qc.to_csv(RESULTS_DIR / "co2_sensor_qc_minute.csv", index=False)
    sensor_audit.to_csv(RESULTS_DIR / "sensor_observation_audit.csv", index=False)
    sensor_normalization.to_csv(RESULTS_DIR / "sensor_normalization_prior.csv", index=False)
    sensor_normalization_sensitivity.to_csv(
        RESULTS_DIR / "sensor_normalization_sensitivity.csv", index=False
    )
    sampling_schedule.to_csv(RESULTS_DIR / "sampling_schedule_used.csv", index=False)
    pulse_schedule.to_csv(RESULTS_DIR / "effective_nutrient_pulses.csv", index=False)
    chemistry_support.to_csv(RESULTS_DIR / "chemical_process_support.csv", index=False)
    ph_do.to_csv(RESULTS_DIR / "ph_do_availability.csv", index=False)
    temperature_inputs.to_csv(RESULTS_DIR / "temperature_model_inputs.csv", index=False)
    temperature_alignment.to_csv(RESULTS_DIR / "temperature_alignment_summary.csv", index=False)
    driver_diagnostics.to_csv(RESULTS_DIR / "driver_diagnostics.csv", index=False)
    inventory.to_csv(RESULTS_DIR / "experiment_inventory.csv", index=False)
    predictions.to_csv(RESULTS_DIR / "prediction_rows.csv", index=False)
    batch_metrics.to_csv(RESULTS_DIR / "batch_metrics.csv", index=False)
    validation.to_csv(RESULTS_DIR / "validation_summary.csv", index=False)
    jacobian.to_csv(RESULTS_DIR / "activation_jacobian_identifiability.csv", index=False)
    singular_values.to_csv(RESULTS_DIR / "activation_jacobian_singular_values.csv", index=False)
    local_correlations.to_csv(RESULTS_DIR / "activation_local_parameter_correlations.csv", index=False)
    activation_profiles.to_csv(RESULTS_DIR / "activation_objective_profiles.csv", index=False)
    loo_estimates.to_csv(RESULTS_DIR / "activation_leave_one_batch_out_estimates.csv", index=False)
    loo_summary.to_csv(RESULTS_DIR / "activation_leave_one_batch_out_summary.csv", index=False)
    pd.concat(comparator_starts, ignore_index=True).to_csv(
        RESULTS_DIR / "comparator_fit_start_diagnostics.csv", index=False
    )
    pd.concat(comparator_predictions, ignore_index=True).to_csv(
        RESULTS_DIR / "comparator_prediction_rows.csv", index=False
    )
    pd.concat(comparator_metrics, ignore_index=True).to_csv(
        RESULTS_DIR / "comparator_batch_metrics.csv", index=False
    )
    all_model_validation.to_csv(RESULTS_DIR / "model_comparison_validation.csv", index=False)
    pd.DataFrame(
        [
            {"parameter": key, "estimate": value}
            for key, value in pooled_fit.items()
            if isinstance(value, (int, float))
        ]
    ).to_csv(RESULTS_DIR / "pooled_pilot_fit_parameters.csv", index=False)
    pooled_starts.to_csv(RESULTS_DIR / "pooled_pilot_fit_starts.csv", index=False)
    pooled_predictions.to_csv(RESULTS_DIR / "pooled_pilot_prediction_rows.csv", index=False)
    pooled_metrics.to_csv(RESULTS_DIR / "pooled_pilot_batch_metrics.csv", index=False)
    cross_scale_predictions.to_csv(RESULTS_DIR / "cross_scale_prediction_rows.csv", index=False)
    cross_scale_metrics.to_csv(RESULTS_DIR / "cross_scale_validation.csv", index=False)

    figures = [
        plot_data_overview(observations, pulse_schedule),
        plot_sensor_filter_examples(sensor_qc, pulse_schedule),
        plot_temperature_profiles(tables.temperature, pulse_schedule),
        plot_chemistry_timeline(tables.primary, pulse_schedule),
        plot_calibration_overlays(predictions, pulse_schedule, MATRICES[0]),
        plot_calibration_overlays(predictions, pulse_schedule, MATRICES[1]),
        plot_cross_lot_validation(predictions, validation, pulse_schedule),
        plot_model_comparison(all_model_validation),
        plot_parameter_comparison(fit_parameters),
        plot_identifiability_profiles(activation_profiles),
        plot_cross_scale_transfer(cross_scale_predictions, cross_scale_metrics),
        plot_ph_do_availability(ph_do),
    ]
    for figure in figures:
        plt.close(figure)

    chart_map = pd.DataFrame(
        [
            ("co2_data_overview.png", "Forma y cobertura de los seis perfiles utilizables", "line_small_multiples"),
            ("sensor_normalization_and_filter.png", "Efecto de normalización y filtro", "multi_series_line"),
            ("temperature_profiles_model_input.png", "Temperatura realmente ejecutada", "line_small_multiples"),
            ("chemistry_and_pulse_timeline.png", "Soporte químico del pulso", "dual_axis_small_multiples"),
            ("calibration_overlays_pilot_lot2.png", "Ajuste lote 2", "observed_vs_predicted_line"),
            ("calibration_overlays_pilot_lot3.png", "Ajuste lote 3", "observed_vs_predicted_line"),
            ("heldout_cross_lot_validation.png", "Transferencia lote 2 y lote 3", "observed_vs_predicted_line"),
            ("heldout_model_comparison.png", "Comparación de modelos", "bar"),
            ("parameter_comparison.png", "Parámetros por lote", "grouped_log_bar"),
            ("activation_identifiability_profiles.png", "Identificabilidad práctica", "profile_line"),
            ("cross_scale_transfer.png", "Transferencia laboratorio-piloto", "observed_vs_predicted_line"),
            ("ph_do_availability.png", "Cobertura pH y O2", "line_with_markers"),
        ],
        columns=["figure", "analytical_question", "chart_variant"],
    )
    chart_map.to_csv(RESULTS_DIR / "chart_map.csv", index=False)

    manifest = {
        "model": MODEL_NAME,
        "scope": "Pilot-2026 CO2 observation layer conditional on fixed Pilot-2026 upstream kinetics",
        "valid_co2_runs": list(VALID_RUNS),
        "excluded_co2_runs": ["26134", "26135", "26136"],
        "holdouts": HOLDOUTS,
        "conversion_ln_min_to_g_l_h": LN_MIN_TO_G_L_H,
        "reactor_volume_l": REACTOR_VOLUME_L,
        "co2_preprocessing": {
            "artificial_leading_zeros": "excluded",
            "observed_first_acquisition_delay_h": sensor_audit.set_index("batch")["first_observed_time_h"].to_dict(),
            "sensor_normalization": "fixed prior background and reactor gains from baseline calibration",
            "sensor_normalization_caveat": "same-campaign prior; no pre-inoculation zero available",
            "sampling_artifact_filter": "local MAD plus sample-window excursion; first 4 h after nutrient pulse protected",
            "hourly_smoother": "3-point median followed by 5-point quadratic Savitzky-Golay within pulse segments",
        },
        "time_alignment": {
            "origin": "first primary chemical sample",
            "end": "active_end from first cooling event or last primary sample",
            "temperature": "Sonda1 full active-process grid",
        },
        "nutrient_pulse": {
            "timing": "linear density crossing at 1040 kg/m3",
            "model_amount_N_kg_m3": YAN_MODEL_PULSE_KG_M3,
            "mass_audit_amount_N_kg_m3": YAN_MASS_AUDIT_KG_M3,
            "lot3_events": "owner-confirmed reconstruction from matching Lot-1 schedule",
        },
        "identifiability_limitations": [
            "no pre-inoculation CO2 baseline",
            "no manual dissolved CO2",
            "pH and dissolved O2 only in Lot 3 at sparse times starting at 16 h",
            "one invalid pH=30.3 row excluded visually",
            "only two calibration fermentations per lot after the holdout",
        ],
        "cross_scale_status": cross_scale_status,
        "cross_scale_interpretation": {
            "laboratory_natural_to_pilot": "primarily a scale-transfer test",
            "laboratory_synthetic_to_pilot": "matrix and scale change together; effects are confounded",
            "pilot_to_laboratory_natural": "reverse scale-transfer test",
            "pilot_to_laboratory_synthetic": "matrix and scale change together; effects are confounded",
        },
        "source_files": {
            "model_dataset": str(MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
            "pilot_upstream_calibration": str(PILOT_CALIBRATION_DIR.relative_to(ROOT_DIR)),
            "co2_minute_qc": str((INTEGRATION_DIR / "co2_minute_qc.csv").relative_to(ROOT_DIR)),
            "laboratory_co2_fit": str((LAB_RESULTS_DIR / "fit_parameters.csv").relative_to(ROOT_DIR)),
        },
    }
    (RESULTS_DIR / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    return {
        "observations": observations,
        "sensor_qc": sensor_qc,
        "sensor_audit": sensor_audit,
        "sensor_normalization": sensor_normalization,
        "sensor_normalization_sensitivity": sensor_normalization_sensitivity,
        "pulse_schedule": pulse_schedule,
        "sampling_schedule": sampling_schedule,
        "chemistry_support": chemistry_support,
        "ph_do": ph_do,
        "inventory": inventory,
        "driver_diagnostics": driver_diagnostics,
        "temperature_inputs": temperature_inputs,
        "temperature_alignment": temperature_alignment,
        "fits": fits,
        "fit_parameters": fit_parameters,
        "fit_starts": fit_starts_table,
        "predictions": predictions,
        "batch_metrics": batch_metrics,
        "validation": validation,
        "jacobian": jacobian,
        "singular_values": singular_values,
        "local_correlations": local_correlations,
        "activation_profiles": activation_profiles,
        "loo_estimates": loo_estimates,
        "loo_summary": loo_summary,
        "model_validation": all_model_validation,
        "pooled_fit": pooled_fit,
        "pooled_metrics": pooled_metrics,
        "cross_scale_predictions": cross_scale_predictions,
        "cross_scale_metrics": cross_scale_metrics,
        "manifest": manifest,
    }


def create_notebook() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    notebook = nbformat.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3"}
    cells = [
        nbformat.v4.new_markdown_cell(
            """# Piloto 2026 · calibración y transferencia del bloque CO₂

## tl;dr

Este notebook reproduce el análisis `solubility_o2_nitrogen_boost_continuous_release` a escala piloto. Calibra por separado los lotes 2 y 3, reserva un ensayo A de cada lote como holdout, transfiere parámetros entre lotes y evalúa transferencia laboratorio↔piloto.

La interpretación del inicio sigue siendo condicional: no existe CO₂ preinoculación observado, no hay CO₂ disuelto manual y pH/O₂ sólo aparecen de forma discreta en el lote 3 desde 16 h. Los ceros artificiales anteriores al inicio de adquisición se excluyen."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contexto y métodos

### Supuestos clave

- Lote 1 (`26134–26136`) permanece fuera de la calibración por perfiles CO₂ no confiables.
- El tiempo cero es la primera muestra química y el final es `active_end`.
- MassView se convierte de Ln/min a g CO₂ L⁻¹ h⁻¹ usando 230 L y 22,414 L mol⁻¹.
- La normalización primaria fija el background y las ganancias por TK de la calibración piloto previa; se incluye sensibilidad sin normalizar.
- Los pulsos se ubican por interpolación del cruce de densidad 1040 kg m⁻³ y se marcan en todos los gráficos CO₂.
- Los parámetros upstream de fermentación quedan fijos en su calibración piloto; sólo se ajusta la capa CO₂."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import sys
from IPython.display import display, Image

ROOT = Path.cwd()
while ROOT != ROOT.parent and not (ROOT / "fermentation_model").exists():
    ROOT = ROOT.parent
if not (ROOT / "fermentation_model").exists():
    raise RuntimeError("Execute from the repository or a descendant directory")
sys.path.insert(0, str(ROOT / "fermentation_model"))
sys.path.insert(0, str(ROOT / "fermentation_model" / "pilot_2026"))

from pilot_2026 import run_co2_solubility_cross_lot_validation_2026 as analysis
result = analysis.run_analysis()
print("Resultados:", analysis.RESULTS_DIR.relative_to(ROOT))"""
        ),
        nbformat.v4.new_markdown_cell("## Datos"),
        nbformat.v4.new_code_cell(
            """display(result["inventory"].round(4))
display(result["sensor_audit"].round(4))
display(result["pulse_schedule"][[
    "matrix", "batch", "pulse_time_h", "schedule_proxy_time_h",
    "bracket_start_h", "bracket_end_h", "timing_source", "lot3_event_reconstructed"
]].round(3))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / "co2_data_overview.png"))
display(Image(filename=analysis.FIGURE_DIR / "sensor_normalization_and_filter.png"))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / "temperature_profiles_model_input.png"))
display(Image(filename=analysis.FIGURE_DIR / "chemistry_and_pulse_timeline.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Resultados de calibración y validación"),
        nbformat.v4.new_code_cell(
            """display(result["fit_parameters"].round(5))
display(result["validation"][[
    "calibration_matrix", "target_matrix", "batch", "role", "rmse_g_l_h",
    "r2", "observed_onset_h", "predicted_onset_h", "onset_delay_h",
    "integral_ratio_pred_over_obs"
]].round(4))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / "calibration_overlays_pilot_lot2.png"))
display(Image(filename=analysis.FIGURE_DIR / "calibration_overlays_pilot_lot3.png"))
display(Image(filename=analysis.FIGURE_DIR / "heldout_cross_lot_validation.png"))"""
        ),
        nbformat.v4.new_code_cell(
            """display(result["model_validation"][[
    "model", "calibration_matrix", "target_matrix", "batch", "role", "rmse_g_l_h", "onset_delay_h"
]].round(4))
display(Image(filename=analysis.FIGURE_DIR / "heldout_model_comparison.png"))
display(Image(filename=analysis.FIGURE_DIR / "parameter_comparison.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Identificabilidad práctica"),
        nbformat.v4.new_code_cell(
            """display(result["jacobian"].round(5))
display(result["loo_summary"].round(5))
activation_corr = result["local_correlations"].query(
    "parameter_1 == 'chem_activation_start_fraction' and parameter_2 == 'chem_activation_duration_fraction'"
)
display(activation_corr.round(5))
display(Image(filename=analysis.FIGURE_DIR / "activation_identifiability_profiles.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Transferencia entre escalas y cobertura auxiliar"),
        nbformat.v4.new_code_cell(
            """display(result["cross_scale_metrics"][[
    "transfer_direction", "calibration_matrix", "target_matrix", "batch", "role", "rmse_g_l_h", "onset_delay_h"
]].round(4))
display(Image(filename=analysis.FIGURE_DIR / "cross_scale_transfer.png"))
display(Image(filename=analysis.FIGURE_DIR / "ph_do_availability.png"))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

La validación cruzada entre lotes prueba si el bloque de liberación/observación CO₂ es estable frente a cambios de lote y programa térmico. La transferencia laboratorio↔piloto evalúa el mismo bloque condicionalmente a los drivers upstream específicos de cada escala.

Los diagnósticos de Jacobiano, perfiles y leave-one-batch-out deben interpretarse con cautela: cada calibración por lote utiliza sólo dos fermentaciones después de reservar el holdout. Esta reproducción no transforma los ceros artificiales en evidencia de ausencia de fermentación y no usa el pH/O₂ tardío para localizar el inicio."""
        ),
        nbformat.v4.new_code_cell(
            """print("Notebook ejecutado sin errores.")
print("Figuras:", len(list(analysis.FIGURE_DIR.glob("*.png"))))
print("Cross-scale:", result["manifest"]["cross_scale_status"])
print("Artefactos:", analysis.RESULTS_DIR.relative_to(ROOT))"""
        ),
    ]
    notebook["cells"] = cells
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    import tempfile

    # Keep this path short: the project lives under a long OneDrive path and
    # IPython otherwise exceeds the legacy Windows MAX_PATH limit.
    runtime_dir = Path(tempfile.gettempdir()) / "pilot26_jupyter"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    # Override (rather than setdefault) because managed Windows sessions can
    # expose a global Jupyter runtime whose ACL cannot be modified here.
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime_dir)
    os.environ["IPYTHONDIR"] = str(runtime_dir / "ipython")
    os.environ["JUPYTER_ALLOW_INSECURE_WRITES"] = "1"
    # jupyter_core reads this flag at import time; nbclient is imported above,
    # so update the module flag explicitly for the OneDrive-mounted workspace.
    from jupyter_core import paths as jupyter_paths

    jupyter_paths.allow_insecure_writes = True
    notebook = nbformat.read(NOTEBOOK_PATH, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT_DIR)}},
    )
    client.execute()
    nbformat.write(notebook, EXECUTED_NOTEBOOK_PATH)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--create-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    create_notebook()
    if not args.create_only:
        execute_notebook(timeout=args.timeout)
    print(NOTEBOOK_PATH)
    if not args.create_only:
        print(EXECUTED_NOTEBOOK_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
