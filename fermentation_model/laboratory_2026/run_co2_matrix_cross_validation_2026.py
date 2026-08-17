from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
_BOOTSTRAP_RUNTIME_DIR = Path(__file__).resolve().parent / "results" / "co2_matrix_cross_validation_2026" / ".jupyter_runtime"
_BOOTSTRAP_RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("JUPYTER_RUNTIME_DIR", str(_BOOTSTRAP_RUNTIME_DIR))
os.environ.setdefault("IPYTHONDIR", str(_BOOTSTRAP_RUNTIME_DIR / "ipython"))
# OneDrive-mounted workspaces can reject the Windows ACL mutation used by
# jupyter_core even though the directory itself is private to this analysis.
os.environ.setdefault("JUPYTER_ALLOW_INSECURE_WRITES", "1")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from nbclient import NotebookClient
from scipy.optimize import least_squares
from scipy.signal import savgol_filter


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SCRIPT_DIR.parent
ROOT_DIR = FERMENTATION_MODEL_DIR.parent
PILOT_2025_DIR = FERMENTATION_MODEL_DIR / "pilot_2025"
for path in (FERMENTATION_MODEL_DIR, SCRIPT_DIR, PILOT_2025_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from shared import run_new_must_glycerol_estimability_doe as base
from laboratory_2026 import run_estimability_historical_by_medium as natural_loader
from laboratory_2026 import run_estimability_historical_synthetic_plus_lot2 as lot2_loader
from laboratory_2026 import run_estimability_old_vs_lot1 as lot1_loader
import run_pilot_2025_co2_solubility_integrated_doe as co2_model


RESULTS_DIR = SCRIPT_DIR / "results" / "co2_matrix_cross_validation_2026"
PLOT_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "co2_solubility_o2_cross_matrix_2026.ipynb"
EXECUTED_NOTEBOOK_PATH = NOTEBOOK_DIR / "co2_solubility_o2_cross_matrix_2026.executed.ipynb"

LOT1_CO2_PATH = (
    SCRIPT_DIR / "results" / "lot1_data_preview" / "processed" / "lot1_co2_filt_volume_corrected_10min.csv"
)
LOT2_CO2_PATH = (
    SCRIPT_DIR / "results" / "lot2_data_preview" / "processed" / "lot2_co2_volume_corrected_10min.csv"
)
LOT1_TEMP_PATH = SCRIPT_DIR / "results" / "lot1_data_preview" / "processed" / "lot1_temp_10min_preview.csv"
LOT2_TEMP_PATH = SCRIPT_DIR / "results" / "lot2_data_preview" / "processed" / "lot2_temperature_10min.csv"
LOT1_SAMPLE_PATH = SCRIPT_DIR / "results" / "lot1_data_preview" / "processed" / "lot1_y15_wide_processed.csv"
LOT2_SAMPLE_PATH = SCRIPT_DIR / "results" / "lot2_data_preview" / "processed" / "lot2_samples_integrated.csv"
NATURAL_THETA_PATH = SCRIPT_DIR / "results" / "estimability_historical_natural" / "theta.csv"
SYNTHETIC_THETA_PATH = (
    SCRIPT_DIR / "results" / "estimability_historical_synthetic_plus_lot2" / "theta_by_case.csv"
)
RAW_NATURAL_DIR = FERMENTATION_MODEL_DIR / "data" / "Laboratorio 2026" / "raw_data"
NUTRIENT_CALENDAR_PATH = Path(r"C:\Users\ctorrealba\Downloads\Fernanda Folch.ics")

MODEL_NAME = "solubility_o2_nitrogen_boost_continuous_release"
THRESHOLD_RELEASE_MODEL_NAME = "solubility_o2_nitrogen_boost_threshold_release"
CHEMISTRY_ALIGNED_MODEL_NAME = "solubility_o2_chemistry_aligned_transition"
LEGACY_MODEL_NAME = "solubility_o2_slow_transition"
HOLDOUTS = {"synthetic": "lot2_F3", "natural": "LAB012"}
SYNTHETIC_CODES = {
    "lot1_F1": "DOE-F01",
    "lot1_F2": "DOE-F02",
    "lot1_F3": "DOE-F03",
    "lot2_F1": "DOE-F04",
    "lot2_F2": "DOE-F05",
    "lot2_F3": "DOE-F06",
}
EXCLUDED_BATCHES = {
    "LAB001": "discarded: unreliable CO2 implementation/profile",
    "LAB002": "discarded: unreliable CO2 implementation/profile",
    "LAB003": "discarded: unreliable CO2 implementation/profile",
    "LAB009": "discarded: unreliable CO2 implementation/profile",
    "lot1_F2": "discarded: unreliable CO2 implementation/profile (DOE-F02)",
}
OBSERVATION_BIN_H = 1.0
DRIVER_GRID_H = 0.25
SAMPLE_ARTIFACT_WINDOW_H = 3.0
SAMPLE_EVENT_HALF_WIDTH_H = 1.25
SAMPLE_DROP_RATIO = 0.65
TRANSIENT_CONTEXT_H = 1.0
TRANSIENT_LOW_RATIO = 0.55
TRANSIENT_HIGH_RATIO = 1.80
TRANSIENT_MIN_DROP_G_L_H = 0.05
TRANSIENT_MIN_SPIKE_G_L_H = 0.08
CO2_DETECTION_LIMIT_G_L_H = 0.05
COLD_CO2_DETECTION_LIMIT_G_L_H = 0.10
COLD_SETPOINT_THRESHOLD_C = 15.5
SMOOTHING_MEDIAN_WINDOW_POINTS = 3
SMOOTHING_SAVGOL_WINDOW_POINTS = 5
SMOOTHING_SAVGOL_POLYORDER = 2
NUTRIENT_DENSITY_TARGET_G_L = 1040.0
MAX_DENSITY_INTERPOLATION_BRACKET_H = 24.0
SPRINGFERM_XTREM_G = 1.00
FDA_G = 0.40
YAN_MG_PER_MG_PRODUCT = 0.20
PULSE_RESPONSE_PROTECTION_H = 4.0
PULSE_PEAK_WINDOW_H = 72.0
PULSE_PEAK_SIGMA_H = 12.0
CHEMISTRY_SUGAR_DROP_G_L = 5.0
CHEMISTRY_ETHANOL_RISE_G_L = 2.0
CHEMISTRY_ACTIVATION_MIN_SCALE_H = 1.0
ONSET_SIGMA_H = 12.0
ONSET_BASELINE_WINDOW_H = 12.0
ONSET_DYNAMIC_RANGE_FRACTION = 0.10
INITIAL_RISE_LOW_FRACTION = 0.02
INITIAL_RISE_HIGH_FRACTION = 0.10
O2_K_MG_L = 0.25
O2_ANA_K_MG_L = 0.75
O2_ANA_HILL = 2.0
O2_CRABTREE_FLOOR = 0.08
CO2_CONTINUOUS_RELEASE_FLOOR = 0.05
EARLY_EMISSION_THRESHOLD_G_L_H = 0.005
SENSOR_ZERO_WINDOW_H = 12.0
SENSOR_ZERO_QUANTILE = 0.10
SENSOR_ZERO_MIN_POINTS = 12
IDENTIFIABILITY_PROFILE_POINTS = 7
IDENTIFIABILITY_PROFILE_MAX_NFEV = 140
LOO_MAX_NFEV = 180
CO2_MOLAR_MASS_G_MOL = 44.01
STANDARD_MOLAR_VOLUME_L_MOL = 22.414
PARAMETER_BOUNDS = {
    "kCO2_release_h": (0.03, 25.0),
    "CO2sat_scale": (0.35, 2.50),
    "O2_qmax_mg_gdw_h": (0.15, 6.0),
    "O2_initial_scale": (0.05, 1.50),
    "pulse_t_rise_h": (1.0, 72.0),
    "pulse_activity_gain": (0.25, 4.0),
    "chem_activation_start_fraction": (0.01, 0.95),
    "chem_activation_duration_fraction": (0.35, 2.0),
    "matrix_gain": (0.10, 20.0),
}

# The physical acquisition files expose channels F1--F3.  DOE-F06 is retained
# with the user's explicit sensor label 6, while its acquisition channel remains
# F3.  Offsets are estimated per experimental run, so the label ambiguity does
# not pool or alter the numerical correction.
SYNTHETIC_SENSOR_LABELS = {
    "lot1_F1": 1,
    "lot1_F2": 2,
    "lot1_F3": 3,
    "lot2_F1": 1,
    "lot2_F2": 2,
    "lot2_F3": 6,
}


@dataclass
class DriverCache:
    batch: str
    matrix: str
    time_h: np.ndarray
    base_qprod_g_l_h: np.ndarray
    base_qprod_no_n_pulse_g_l_h: np.ndarray
    n_pulse_qprod_increment_g_l_h: np.ndarray
    biomass_g_l: np.ndarray
    biomass_no_n_pulse_g_l: np.ndarray
    n_pulse_biomass_increment_g_l: np.ndarray
    n_pulse_time_h: float
    n_pulse_amount_kg_m3: float
    o2_saturation_base_mg_l: np.ndarray
    chemical_activation: np.ndarray
    chemical_activation_previous: np.ndarray
    chemical_activity_lower_h: float
    chemical_activity_upper_h: float
    chemical_activity_center_h: float
    chemical_activity_scale_h: float
    csat_base_g_l: np.ndarray
    temperature_used_c: np.ndarray
    observation_time_h: np.ndarray
    theta_source: str


def _load_theta(path: Path, case: str | None = None) -> dict[str, float]:
    frame = pd.read_csv(path)
    if case is not None:
        frame = frame[frame["case"].astype(str).eq(case)].copy()
    if frame.empty:
        raise RuntimeError(f"No parameter rows found in {path} for case={case!r}")
    theta = dict(base.DEFAULT_THETA)
    theta.update(dict(zip(frame["parameter"].astype(str), pd.to_numeric(frame["theta"], errors="raise"))))
    return {name: float(value) for name, value in theta.items()}


def augment_lot1_batches_with_full_temperature(
    batches: list[base.BatchData],
) -> list[base.BatchData]:
    """Put the native Lot-1 temperature series on the model time grid.

    The original Lot-1 loader sampled temperature only at chemistry times.
    Interpolation between those sparse points smoothed programmed temperature
    steps. Here the full sensor series is merged with state-observation times.
    """

    temperature = pd.read_csv(LOT1_TEMP_PATH)
    for column in ("t_h", "T"):
        temperature[column] = pd.to_numeric(temperature[column], errors="coerce")
    output = []
    for batch in batches:
        process = batch.batch.removeprefix("lot1_")
        sensor = (
            temperature[temperature["process"].astype(str).eq(process)]
            .dropna(subset=["t_h", "T"])
            .sort_values("t_h")
            .drop_duplicates("t_h", keep="last")
        )
        if sensor.empty:
            output.append(batch)
            continue
        sensor_time = sensor["t_h"].to_numpy(dtype=float)
        sensor_temperature = sensor["T"].to_numpy(dtype=float)
        new_time = np.asarray(
            sorted(set(float(value) for value in batch.time).union(float(value) for value in sensor_time)),
            dtype=float,
        )
        actual_temperature = np.interp(new_time, sensor_time, sensor_temperature)
        nominal_temperature = np.interp(new_time, batch.time, batch.temperature_c)
        actual_temperature = np.where(
            (new_time < sensor_time.min()) | (new_time > sensor_time.max()),
            nominal_temperature,
            actual_temperature,
        )
        output.append(
            replace(
                batch,
                time=new_time,
                temperature_c=actual_temperature,
                observations=natural_loader.remap_observations(
                    batch.time, batch.observations, new_time
                ),
            )
        )
    return output


def _unfold_ics_lines(text: str) -> list[str]:
    """Undo RFC 5545 line folding without requiring a calendar dependency."""

    unfolded: list[str] = []
    for line in text.splitlines():
        if line.startswith((" ", "\t")) and unfolded:
            unfolded[-1] += line[1:]
        else:
            unfolded.append(line)
    return unfolded


def _calendar_nutrient_pulse_events(path: Path) -> pd.DataFrame:
    """Extract the in-process LAB nutrient events (pulse 2) from the ICS file."""

    if not path.exists():
        raise FileNotFoundError(f"Nutrient calendar not found: {path}")
    events: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    for line in _unfold_ics_lines(path.read_text(encoding="utf-8-sig")):
        if line == "BEGIN:VEVENT":
            current = {}
            continue
        if line == "END:VEVENT":
            if current is not None:
                events.append(current)
            current = None
            continue
        if current is None or ":" not in line:
            continue
        key, value = line.split(":", 1)
        current[key] = value.replace(r"\n", "\n").replace(r"\,", ",").replace(r"\;", ";")

    rows: list[dict[str, object]] = []
    for event in events:
        text = "\n".join([event.get("SUMMARY", ""), event.get("DESCRIPTION", "")])
        if "Pulso nutricional 2" not in text:
            continue
        match = re.search(r"(?:ING26-)?(LAB\d{3})", text, flags=re.IGNORECASE)
        if match is None:
            continue
        batch = match.group(1).upper()
        number = int(batch[-3:])
        if number < 4:
            continue
        dt_key = next((key for key in event if key.startswith("DTSTART")), None)
        if dt_key is None:
            raise RuntimeError(f"No DTSTART found for nutrient event {batch}")
        raw_dt = event[dt_key]
        if raw_dt.endswith("Z"):
            timestamp_utc = pd.to_datetime(raw_dt, format="%Y%m%dT%H%M%SZ", utc=True)
            timestamp_local = timestamp_utc.tz_convert("America/Santiago")
        else:
            timestamp_local = pd.to_datetime(raw_dt, format="%Y%m%dT%H%M%S").tz_localize(
                "America/Santiago"
            )
            timestamp_utc = timestamp_local.tz_convert("UTC")
        product_match = re.search(r"Producto/Nutriente:\s*([^\n]+)", text)
        dose_match = re.search(r"Dosis:\s*([^\n]+)", text)
        state_match = re.search(r"Estado fermentación:\s*([^\n]+)", text)
        rows.append(
            {
                "batch": batch,
                "calendar_timestamp_utc": timestamp_utc.isoformat(),
                "calendar_timestamp_local": timestamp_local.isoformat(),
                "calendar_product": product_match.group(1).strip() if product_match else "",
                "calendar_dose": dose_match.group(1).strip() if dose_match else "",
                "calendar_process_state": state_match.group(1).strip() if state_match else "",
                "calendar_summary": event.get("SUMMARY", ""),
            }
        )
    frame = pd.DataFrame(rows).sort_values("batch").reset_index(drop=True)
    expected = {f"LAB{number:03d}" for number in range(4, 13)}
    missing = sorted(expected.difference(frame["batch"].astype(str)))
    if missing:
        raise RuntimeError(f"Missing in-process nutrient events in {path.name}: {missing}")
    if frame["batch"].duplicated().any():
        duplicates = frame.loc[frame["batch"].duplicated(False), "batch"].tolist()
        raise RuntimeError(f"Duplicate in-process nutrient events in {path.name}: {duplicates}")
    return frame


def _density_crossing(
    group: pd.DataFrame,
    target_g_l: float,
    reference_time_h: float,
) -> dict[str, float]:
    """Locate a density crossing and retain the measurement bracket used."""

    valid = (
        group[["time_h", "density"]]
        .apply(pd.to_numeric, errors="coerce")
        .dropna()
        .sort_values("time_h")
        .drop_duplicates("time_h", keep="last")
    )
    candidates: list[dict[str, float]] = []
    time_h = valid["time_h"].to_numpy(dtype=float)
    density = valid["density"].to_numpy(dtype=float)
    for index in range(len(valid)):
        if abs(density[index] - target_g_l) <= 1e-9:
            candidates.append(
                {
                    "density_crossing_t_h": float(time_h[index]),
                    "density_bracket_start_h": float(time_h[index]),
                    "density_bracket_end_h": float(time_h[index]),
                    "density_bracket_width_h": 0.0,
                }
            )
        if index == len(valid) - 1:
            continue
        delta_left = density[index] - target_g_l
        delta_right = density[index + 1] - target_g_l
        if delta_left * delta_right < 0.0:
            fraction = (target_g_l - density[index]) / (density[index + 1] - density[index])
            crossing = time_h[index] + fraction * (time_h[index + 1] - time_h[index])
            candidates.append(
                {
                    "density_crossing_t_h": float(crossing),
                    "density_bracket_start_h": float(time_h[index]),
                    "density_bracket_end_h": float(time_h[index + 1]),
                    "density_bracket_width_h": float(time_h[index + 1] - time_h[index]),
                }
            )
    if not candidates:
        return {
            "density_crossing_t_h": np.nan,
            "density_bracket_start_h": np.nan,
            "density_bracket_end_h": np.nan,
            "density_bracket_width_h": np.nan,
        }
    return min(candidates, key=lambda row: abs(row["density_crossing_t_h"] - reference_time_h))


def natural_nutrient_pulse_schedule(
    chemistry: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Build one effective in-process N pulse per LAB004--LAB012.

    Chemical density controls timing for every LAB experiment whenever a 1040
    g/L crossing can be interpolated.  The interpolation bracket is retained as
    timing uncertainty even when it is wide.  Calendar time is only a fallback
    if the chemistry never brackets 1040 g/L.  The inoculation pulse is omitted
    because its YAN is already present in the t=0 chemistry.
    """

    calendar = _calendar_nutrient_pulse_events(NUTRIENT_CALENDAR_PATH)
    meta = metadata.set_index("batch")
    total_yan_mg = (
        SPRINGFERM_XTREM_G * 1000.0 * YAN_MG_PER_MG_PRODUCT
        + FDA_G * 1000.0 * YAN_MG_PER_MG_PRODUCT
    )
    rows: list[dict[str, object]] = []
    for event in calendar.itertuples(index=False):
        batch = str(event.batch)
        if batch not in meta.index:
            raise RuntimeError(f"No natural-batch metadata found for {batch}")
        t0 = pd.Timestamp(meta.loc[batch, "t0"]).tz_localize("America/Santiago")
        calendar_local = pd.Timestamp(event.calendar_timestamp_local)
        calendar_t_h = float((calendar_local - t0).total_seconds() / 3600.0)
        crossing = _density_crossing(
            chemistry[chemistry["batch"].astype(str).eq(batch)],
            NUTRIENT_DENSITY_TARGET_G_L,
            calendar_t_h,
        )
        use_density = np.isfinite(crossing["density_crossing_t_h"])
        model_time_h = (
            crossing["density_crossing_t_h"] if use_density else calendar_t_h
        )
        if use_density:
            timing_source = (
                "chemical_density_linear_interpolation"
                if crossing["density_bracket_width_h"] <= MAX_DENSITY_INTERPOLATION_BRACKET_H
                else "chemical_density_linear_interpolation_sparse_bracket"
            )
        else:
            timing_source = "calendar_timestamp_missing_density_crossing"
        volume_l = float(meta.loc[batch, "reactor_volume_l"])
        rows.append(
            {
                **event._asdict(),
                "matrix": "natural",
                "calendar_t_h": calendar_t_h,
                **crossing,
                "density_target_g_l": NUTRIENT_DENSITY_TARGET_G_L,
                "model_pulse_time_h": float(model_time_h),
                "timing_source": timing_source,
                "calendar_minus_model_h": float(calendar_t_h - model_time_h),
                "reactor_volume_l": volume_l,
                "springferm_xtrem_g": SPRINGFERM_XTREM_G,
                "fda_g": FDA_G,
                "yan_mg_per_mg_product": YAN_MG_PER_MG_PRODUCT,
                "total_yan_mg": total_yan_mg,
                "amount_N_kg_m3": total_yan_mg / volume_l / 1000.0,
                "initial_pulse_at_t0_in_model": False,
                "composition_source": "user-specified common LAB protocol",
                "calendar_source": str(NUTRIENT_CALENDAR_PATH),
                "excluded_from_co2_analysis": batch in EXCLUDED_BATCHES,
            }
        )
    return pd.DataFrame(rows).sort_values("batch").reset_index(drop=True)


def override_natural_nutrient_pulses(
    batches: list[base.BatchData],
    schedule: pd.DataFrame,
) -> list[base.BatchData]:
    pulse_by_batch = schedule.set_index("batch")
    output = []
    for batch in batches:
        if batch.batch not in pulse_by_batch.index:
            output.append(batch)
            continue
        row = pulse_by_batch.loc[batch.batch]
        pulses = dict(batch.pulses)
        pulses["N"] = ((float(row["model_pulse_time_h"]), float(row["amount_N_kg_m3"])),)
        output.append(replace(batch, pulses=pulses))
    return output


def chemical_support_table(batches: dict[str, base.BatchData]) -> pd.DataFrame:
    """First/last analytical state observation for each modeled fermentation."""

    rows = []
    for batch_name, batch in sorted(batches.items()):
        observed = np.zeros(len(batch.time), dtype=bool)
        for state in base.STATE_NAMES:
            values = np.asarray(batch.observations.get(state, np.full(len(batch.time), np.nan)), dtype=float)
            observed |= np.isfinite(values)
        times = np.asarray(batch.time, dtype=float)[observed]
        if not len(times):
            continue
        rows.append(
            {
                "batch": batch_name,
                "matrix": "synthetic" if batch_name.startswith("lot") else "natural",
                "chemistry_first_h": float(np.min(times)),
                "chemistry_last_h": float(np.max(times)),
                "n_chemical_timepoints": int(len(np.unique(times))),
            }
        )
    return pd.DataFrame(rows).sort_values(["matrix", "batch"]).reset_index(drop=True)


def effective_nutrient_pulse_table(
    batches: dict[str, base.BatchData],
    natural_schedule: pd.DataFrame,
) -> pd.DataFrame:
    """Combine overridden LAB pulses and existing synthetic N-pulse schedules."""

    natural_lookup = natural_schedule.set_index("batch")
    rows: list[dict[str, object]] = []
    for batch_name, batch in sorted(batches.items()):
        for pulse_time_h, amount in batch.pulses.get("N", tuple()):
            row = {
                "matrix": "synthetic" if batch_name.startswith("lot") else "natural",
                "batch": batch_name,
                "experiment_code": SYNTHETIC_CODES.get(batch_name, batch_name),
                "pulse_time_h": float(pulse_time_h),
                "amount_N_kg_m3": float(amount),
                "timing_source": "process_time_existing_batch_schedule",
                "density_target_g_l": np.nan,
                "density_bracket_width_h": np.nan,
                "calendar_t_h": np.nan,
                "excluded_from_co2_analysis": batch_name in EXCLUDED_BATCHES,
            }
            if batch_name in natural_lookup.index:
                detail = natural_lookup.loc[batch_name]
                row.update(
                    {
                        "timing_source": str(detail["timing_source"]),
                        "density_target_g_l": float(detail["density_target_g_l"]),
                        "density_bracket_width_h": float(detail["density_bracket_width_h"]),
                        "calendar_t_h": float(detail["calendar_t_h"]),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["matrix", "batch", "pulse_time_h"]).reset_index(drop=True)


def load_batches() -> tuple[dict[str, base.BatchData], dict[str, dict[str, float]], dict[str, pd.DataFrame]]:
    """Load the 2026 batch objects and fixed upstream parameter sets.

    Only the CO2 observation layer is calibrated below. The core fermentation
    drivers are retained from the existing medium-specific 2026 calibrations.
    """

    natural_chemistry, natural_batches, _, natural_tables = natural_loader.make_medium_batches("natural")
    natural_schedule = natural_nutrient_pulse_schedule(
        natural_chemistry, natural_tables["metadata"]
    )
    natural_batches = override_natural_nutrient_pulses(natural_batches, natural_schedule)
    natural_tables["chemistry"] = natural_chemistry
    natural_tables["nutrient_pulses"] = natural_schedule
    lot1_batches, _ = lot1_loader.make_lot1_batches()
    lot1_batches = augment_lot1_batches_with_full_temperature(lot1_batches)
    lot2_batches, _, _ = lot2_loader.make_lot2_batches()
    batches = {batch.batch: batch for batch in natural_batches + lot1_batches + lot2_batches}

    theta_by_matrix = {
        "natural": _load_theta(NATURAL_THETA_PATH),
        "synthetic": _load_theta(SYNTHETIC_THETA_PATH, "historical_synthetic_plus_lot2"),
    }
    return batches, theta_by_matrix, natural_tables


def _sccm_to_g_l_h(flow_sccm: pd.Series, volume_l: float) -> pd.Series:
    gas_l_h = pd.to_numeric(flow_sccm, errors="coerce") * 60.0 / 1000.0
    return gas_l_h * CO2_MOLAR_MASS_G_MOL / STANDARD_MOLAR_VOLUME_L_MOL / float(volume_l)


def _sensor_assignment(batch: str) -> tuple[str, int, str]:
    """Return acquisition channel, user-facing sensor label and provenance."""

    batch = str(batch)
    synthetic_match = re.search(r"_(F[123])$", batch)
    if synthetic_match:
        channel = synthetic_match.group(1)
        return (
            channel,
            int(SYNTHETIC_SENSOR_LABELS.get(batch, int(channel[1:]))),
            "explicit DOE mapping; DOE-F06 label 6 is acquired on physical channel F3",
        )
    natural_match = re.fullmatch(r"LAB(\d{3})", batch, flags=re.IGNORECASE)
    if natural_match:
        channel_number = (int(natural_match.group(1)) - 1) % 3 + 1
        return (
            f"F{channel_number}",
            channel_number,
            "cyclic LAB assignment inferred from repository F1/F2/F3 source filenames",
        )
    return "unknown", -1, "unresolved"


def apply_sensor_zero_correction(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Estimate a run-specific sensor zero from the initial low-flow plateau.

    The offset is the non-negative 10th percentile within the first 12 h of
    available process data.  A low quantile is deliberately used instead of the
    initial median because warm fermentations may already be increasing, while
    isolated start-up spikes should not define zero.  Correction happens in
    native sccm before physical clipping and volume conversion.
    """

    corrected: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for batch, group in raw.groupby("batch", sort=False):
        group = group.sort_values("t_h").copy()
        channel, sensor_id, assignment_source = _sensor_assignment(str(batch))
        time_h = pd.to_numeric(group["t_h"], errors="coerce")
        native = pd.to_numeric(group["sensor_flow_sccm"], errors="coerce")
        volume_l = pd.to_numeric(group["reactor_volume_l"], errors="coerce")
        first_h = float(time_h.min())
        baseline_mask = time_h.le(first_h + SENSOR_ZERO_WINDOW_H) & native.notna()
        candidates = native[baseline_mask]
        if len(candidates) < SENSOR_ZERO_MIN_POINTS:
            candidates = native.dropna().iloc[:SENSOR_ZERO_MIN_POINTS]
            method = "first_available_points_q10"
        else:
            method = "initial_12h_q10"
        raw_quantile = float(candidates.quantile(SENSOR_ZERO_QUANTILE)) if len(candidates) else 0.0
        offset_sccm = max(0.0, raw_quantile)
        corrected_native = native - offset_sccm
        conversion = (
            60.0 / 1000.0 * CO2_MOLAR_MASS_G_MOL / STANDARD_MOLAR_VOLUME_L_MOL
            / volume_l
        )
        group["acquisition_channel"] = channel
        group["sensor_id"] = sensor_id
        group["sensor_assignment_source"] = assignment_source
        group["sensor_zero_method"] = method
        group["sensor_zero_offset_sccm"] = offset_sccm
        group["sensor_flow_zero_corrected_sccm"] = corrected_native
        group["co2_rate_signed_g_l_h"] = native * conversion
        group["co2_rate_uncorrected_physical_g_l_h"] = np.maximum(
            group["co2_rate_signed_g_l_h"], 0.0
        )
        group["co2_rate_zero_corrected_signed_g_l_h"] = corrected_native * conversion
        group["co2_rate_raw_g_l_h"] = np.maximum(
            group["co2_rate_zero_corrected_signed_g_l_h"], 0.0
        )
        corrected.append(group)
        initial_corrected = corrected_native[baseline_mask]
        audit_rows.append(
            {
                "matrix": str(group["matrix"].iloc[0]),
                "batch": str(batch),
                "experiment_code": str(group["experiment_code"].iloc[0]),
                "lot": str(group["lot"].iloc[0]),
                "acquisition_channel": channel,
                "sensor_id": sensor_id,
                "sensor_assignment_source": assignment_source,
                "zero_method": method,
                "zero_window_h": SENSOR_ZERO_WINDOW_H,
                "zero_quantile": SENSOR_ZERO_QUANTILE,
                "n_zero_candidates": int(len(candidates)),
                "initial_native_q10_sccm": raw_quantile,
                "estimated_zero_offset_sccm": offset_sccm,
                "estimated_zero_offset_g_l_h": float(
                    np.nanmedian(offset_sccm * conversion[baseline_mask])
                ),
                "initial_native_median_sccm": float(candidates.median()) if len(candidates) else np.nan,
                "initial_corrected_q10_sccm": float(
                    pd.Series(initial_corrected).quantile(SENSOR_ZERO_QUANTILE)
                ) if len(initial_corrected) else np.nan,
                "initial_corrected_median_sccm": float(np.nanmedian(initial_corrected))
                if len(initial_corrected) else np.nan,
            }
        )
    return (
        pd.concat(corrected, ignore_index=True),
        pd.DataFrame(audit_rows).sort_values(["matrix", "batch"]).reset_index(drop=True),
    )


def load_sampling_schedule() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    lot1 = pd.read_csv(LOT1_SAMPLE_PATH)
    for process, group in lot1.groupby("process"):
        for time_h in pd.to_numeric(group["t_h"], errors="coerce").dropna().unique():
            rows.append({"batch": f"lot1_{process}", "sample_time_h": float(time_h), "source": "lot1_y15"})

    lot2 = pd.read_csv(LOT2_SAMPLE_PATH)
    for process, group in lot2.groupby("process"):
        for time_h in pd.to_numeric(group["t_h"], errors="coerce").dropna().unique():
            rows.append({"batch": f"lot2_{process}", "sample_time_h": float(time_h), "source": "lot2_integrated_samples"})

    natural = natural_loader.load_historical_medium_data("natural")
    for batch, group in natural.groupby("batch"):
        for time_h in pd.to_numeric(group["time_h"], errors="coerce").dropna().unique():
            rows.append({"batch": str(batch), "sample_time_h": float(time_h), "source": "natural_homologated_samples"})
    return pd.DataFrame(rows).drop_duplicates(["batch", "sample_time_h"]).sort_values(["batch", "sample_time_h"])


def excluded_experiment_table() -> pd.DataFrame:
    rows = []
    for batch, reason in EXCLUDED_BATCHES.items():
        matrix = "synthetic" if batch.startswith("lot") else "natural"
        rows.append(
            {
                "matrix": matrix,
                "batch": batch,
                "experiment_code": SYNTHETIC_CODES.get(batch, batch),
                "reason": reason,
                "used_for_calibration": False,
                "used_for_validation": False,
                "shown_in_profile_figures": False,
            }
        )
    return pd.DataFrame(rows).sort_values(["matrix", "batch"]).reset_index(drop=True)


def _attach_temperature_context(raw: pd.DataFrame, natural_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    temperature_frames = []
    lot1 = pd.read_csv(LOT1_TEMP_PATH)
    lot1["batch"] = "lot1_" + lot1["process"].astype(str)
    temperature_frames.append(lot1[["batch", "t_h", "T", "SP"]])

    lot2 = pd.read_csv(LOT2_TEMP_PATH)
    lot2["batch"] = "lot2_" + lot2["process"].astype(str)
    temperature_frames.append(lot2[["batch", "t_h", "T", "SP"]])

    natural = natural_tables["temperature"].copy()
    temperature_frames.append(natural[["batch", "t_h", "T", "SP"]])
    temperature = pd.concat(temperature_frames, ignore_index=True)
    for column in ("t_h", "T", "SP"):
        temperature[column] = pd.to_numeric(temperature[column], errors="coerce")

    out = []
    for batch, group in raw.groupby("batch", sort=False):
        group = group.copy()
        sensor = temperature[temperature["batch"].eq(batch)].dropna(subset=["t_h"]).sort_values("t_h")
        for source, target in (("T", "sensor_temperature_c"), ("SP", "setpoint_c")):
            valid = sensor.dropna(subset=[source])
            if valid.empty:
                group[target] = np.nan
            else:
                group[target] = np.interp(
                    group["t_h"].to_numpy(dtype=float),
                    valid["t_h"].to_numpy(dtype=float),
                    valid[source].to_numpy(dtype=float),
                    left=np.nan,
                    right=np.nan,
                )
        out.append(group)
    return pd.concat(out, ignore_index=True)


def filter_co2_sensor_artifacts(
    raw: pd.DataFrame,
    sampling_schedule: pd.DataFrame,
    nutrient_pulses: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Replace only short, context-reversing CO2 excursions.

    A point is eligible when the medians immediately before and after it agree,
    while the point/short run is much lower or higher than that local baseline.
    Around a recorded sampling time, a second event-level rule captures the
    complete downward pulse when the pre/post-event medians agree. This protects
    sustained biological transitions while avoiding partial replacement of a
    sampling fall/rebound event.
    """

    output = []
    schedule_by_batch = {
        str(batch): group["sample_time_h"].to_numpy(dtype=float)
        for batch, group in sampling_schedule.groupby("batch")
    }
    pulse_by_batch = (
        {
            str(batch): group["pulse_time_h"].to_numpy(dtype=float)
            for batch, group in nutrient_pulses.groupby("batch")
        }
        if nutrient_pulses is not None and not nutrient_pulses.empty
        else {}
    )
    for batch, group in raw.groupby("batch", sort=False):
        group = group.sort_values("t_h").copy()
        time_h = group["t_h"].to_numpy(dtype=float)
        physical = group["co2_rate_raw_g_l_h"].to_numpy(dtype=float)
        if len(group) > 1:
            median_dt = float(np.nanmedian(np.diff(time_h)))
        else:
            median_dt = OBSERVATION_BIN_H
        context_points = max(3, int(round(TRANSIENT_CONTEXT_H / max(median_dt, 1e-6))))
        series = pd.Series(physical, index=group.index)
        before = series.rolling(context_points, min_periods=max(2, context_points // 2)).median().shift(1)
        after = (
            series.iloc[::-1]
            .rolling(context_points, min_periods=max(2, context_points // 2))
            .median()
            .shift(1)
            .iloc[::-1]
        )
        baseline = 0.5 * (before + after)
        context_agreement = (before - after).abs() <= np.maximum(0.35 * baseline.abs(), 0.05)
        low_excursion = (
            baseline.gt(0.04)
            & series.lt(TRANSIENT_LOW_RATIO * baseline)
            & (baseline - series).gt(TRANSIENT_MIN_DROP_G_L_H)
        )
        high_excursion = (
            baseline.gt(0.04)
            & series.gt(TRANSIENT_HIGH_RATIO * baseline)
            & (series - baseline).gt(TRANSIENT_MIN_SPIKE_G_L_H)
        )
        context_artifact = (
            context_agreement & (low_excursion | high_excursion)
        ).fillna(False).to_numpy(dtype=bool)

        pulse_times = pulse_by_batch.get(str(batch), np.array([], dtype=float))
        pulse_response_protected = np.zeros(len(group), dtype=bool)
        for pulse_time in pulse_times:
            pulse_response_protected |= (
                (time_h > pulse_time) & (time_h <= pulse_time + PULSE_RESPONSE_PROTECTION_H)
            )
        context_artifact &= ~pulse_response_protected

        sample_times = schedule_by_batch.get(str(batch), np.array([], dtype=float))
        if len(sample_times):
            distance = np.min(np.abs(time_h[:, None] - sample_times[None, :]), axis=1)
        else:
            distance = np.full(len(group), np.nan)
        near_sample = np.isfinite(distance) & (distance <= SAMPLE_ARTIFACT_WINDOW_H)

        sample_drop = np.zeros(len(group), dtype=bool)
        sample_event_baseline = np.full(len(group), np.nan, dtype=float)
        sample_event_time = np.full(len(group), np.nan, dtype=float)
        for sample_time in sample_times:
            before_event = (
                (time_h >= sample_time - SAMPLE_ARTIFACT_WINDOW_H)
                & (time_h <= sample_time - SAMPLE_EVENT_HALF_WIDTH_H)
            )
            after_event = (
                (time_h >= sample_time + SAMPLE_EVENT_HALF_WIDTH_H)
                & (time_h <= sample_time + SAMPLE_ARTIFACT_WINDOW_H)
            )
            event_window = np.abs(time_h - sample_time) <= SAMPLE_EVENT_HALF_WIDTH_H
            replaceable_event_window = event_window & ~pulse_response_protected
            if before_event.sum() < 3 or after_event.sum() < 3:
                continue
            before_level = float(np.nanmedian(physical[before_event]))
            after_level = float(np.nanmedian(physical[after_event]))
            event_baseline = 0.5 * (before_level + after_level)
            levels_agree = abs(before_level - after_level) <= max(
                0.35 * abs(event_baseline), 0.05
            )
            event_drop = (
                event_window
                & levels_agree
                & (event_baseline > 0.04)
                & (physical < SAMPLE_DROP_RATIO * event_baseline)
                & ((event_baseline - physical) > TRANSIENT_MIN_DROP_G_L_H)
            )
            if event_drop.any():
                # Once a sampling disturbance is confirmed, replace the whole
                # central event. This removes the fall/rebound pair instead of
                # leaving an artificial overshoot at the edge of a partial mask.
                assign = replaceable_event_window & (
                    ~np.isfinite(sample_event_time)
                    | (np.abs(time_h - sample_time) < np.abs(time_h - sample_event_time))
                )
                sample_drop |= replaceable_event_window
                sample_event_baseline[assign] = event_baseline
                sample_event_time[assign] = sample_time

        artifact = context_artifact | sample_drop

        filtered = physical.copy()
        valid = ~artifact & np.isfinite(physical)
        if valid.sum() >= 2 and artifact.any():
            filtered[artifact] = np.interp(time_h[artifact], time_h[valid], physical[valid])
        use_event_baseline = sample_drop & np.isfinite(sample_event_baseline)
        filtered[use_event_baseline] = sample_event_baseline[use_event_baseline]
        group["artifact_flag"] = artifact
        group["artifact_reason"] = np.select(
            [artifact & near_sample, artifact & ~near_sample],
            ["sampling_window_transient", "short_transient_outlier"],
            default="kept",
        )
        group["nearest_sample_distance_h"] = distance
        if len(sample_times):
            nearest_index = np.argmin(np.abs(time_h[:, None] - sample_times[None, :]), axis=1)
            nearest_sample_time = sample_times[nearest_index]
            fallback = context_artifact & near_sample & ~np.isfinite(sample_event_time)
            sample_event_time[fallback] = nearest_sample_time[fallback]
        group["artifact_sample_time_h"] = sample_event_time
        audit_baseline = baseline.to_numpy(dtype=float)
        audit_baseline[sample_drop] = sample_event_baseline[sample_drop]
        group["local_context_baseline_g_l_h"] = audit_baseline
        group["co2_rate_filtered_g_l_h"] = filtered
        group["pulse_response_protected"] = pulse_response_protected
        group["cold_operation"] = pd.to_numeric(group["setpoint_c"], errors="coerce").le(
            COLD_SETPOINT_THRESHOLD_C
        )
        group["detection_limit_g_l_h"] = np.where(
            group["cold_operation"],
            COLD_CO2_DETECTION_LIMIT_G_L_H,
            CO2_DETECTION_LIMIT_G_L_H,
        )
        group["below_detection_limit"] = group["co2_rate_filtered_g_l_h"].le(
            group["detection_limit_g_l_h"]
        )
        group["cold_low_sensitivity"] = group["cold_operation"] & group["below_detection_limit"]
        output.append(group)
    return pd.concat(output, ignore_index=True)


def smooth_hourly_co2_profiles(
    hourly: pd.DataFrame,
    nutrient_pulses: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Apply a short robust smoother after transient reconstruction.

    A centered 3-point median suppresses remaining one-hour impulses and a
    5-point quadratic Savitzky--Golay pass reduces curvature noise. The raw and
    transient-repaired signals remain in separate columns for audit.
    """

    output = []
    pulse_by_batch = (
        {
            str(batch): np.sort(group["pulse_time_h"].to_numpy(dtype=float))
            for batch, group in nutrient_pulses.groupby("batch")
        }
        if nutrient_pulses is not None and not nutrient_pulses.empty
        else {}
    )
    for batch, group in hourly.groupby("batch", sort=False):
        group = group.sort_values("t_h").copy()
        repaired = group["co2_rate_filtered_g_l_h"].to_numpy(dtype=float)
        time_h = group["t_h"].to_numpy(dtype=float)
        segment_id = np.searchsorted(
            pulse_by_batch.get(str(batch), np.array([], dtype=float)), time_h, side="right"
        )
        robust = np.empty_like(repaired)
        smoothed = np.empty_like(repaired)
        for segment in np.unique(segment_id):
            mask = segment_id == segment
            segment_repaired = repaired[mask]
            segment_robust = (
                pd.Series(segment_repaired)
                .rolling(SMOOTHING_MEDIAN_WINDOW_POINTS, center=True, min_periods=1)
                .median()
                .to_numpy(dtype=float)
            )
            if len(segment_robust) >= SMOOTHING_SAVGOL_WINDOW_POINTS:
                segment_smoothed = savgol_filter(
                    segment_robust,
                    window_length=SMOOTHING_SAVGOL_WINDOW_POINTS,
                    polyorder=SMOOTHING_SAVGOL_POLYORDER,
                    mode="interp",
                )
            else:
                segment_smoothed = segment_robust
            robust[mask] = segment_robust
            smoothed[mask] = segment_smoothed
        group["co2_rate_robust_median_g_l_h"] = robust
        group["co2_rate_smoothed_g_l_h"] = np.maximum(smoothed, 0.0)
        group["smoothing_segment"] = segment_id
        output.append(group)
    return pd.concat(output, ignore_index=True).sort_values(["matrix", "batch", "t_h"]).reset_index(drop=True)


def _natural_early_sensor_rows() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the CO2-only LAB001-LAB003 records for coverage/QC.

    These runs have process sensor files in the repository but no matching
    chemistry/state sheets in mosto_natural_xthiol.xlsx. They therefore remain
    descriptive and are not silently used as model-calibration batches.
    """

    rows: list[pd.DataFrame] = []
    inventory: list[dict[str, object]] = []
    top_level_files = list(RAW_NATURAL_DIR.glob("*.csv"))
    for number in range(1, 4):
        batch = f"LAB{number:03d}"
        pattern = re.compile(rf"lab0*{number}(?!\d)", flags=re.IGNORECASE)
        files = [
            path
            for path in top_level_files
            if pattern.search(path.name)
            and ("co2_filt" in path.name.lower() or "co2f" in path.name.lower())
        ]
        input_path = RAW_NATURAL_DIR / f"ing26_lab{number:03d}.csv"
        t0 = pd.NaT
        if input_path.exists():
            input_frame = pd.read_csv(input_path, usecols=lambda column: column == "timestamp")
            t0 = pd.to_datetime(input_frame["timestamp"], errors="coerce").min()

        raw_frames: list[pd.DataFrame] = []
        for path in sorted(files):
            frame = pd.read_csv(path, usecols=lambda column: column in {"timestamp", "flow_filt_sccm"})
            if {"timestamp", "flow_filt_sccm"}.issubset(frame.columns):
                frame["timestamp"] = pd.to_datetime(frame["timestamp"], errors="coerce")
                frame["flow_filt_sccm"] = pd.to_numeric(frame["flow_filt_sccm"], errors="coerce")
                frame["source_file"] = str(path.relative_to(ROOT_DIR))
                raw_frames.append(frame.dropna(subset=["timestamp", "flow_filt_sccm"]))

        if raw_frames and pd.notna(t0):
            raw = (
                pd.concat(raw_frames, ignore_index=True)
                .sort_values(["timestamp", "source_file"])
                .drop_duplicates("timestamp", keep="first")
            )
            raw["t_h"] = (raw["timestamp"] - t0).dt.total_seconds() / 3600.0
            raw = raw[raw["t_h"].ge(0.0)].copy()
            raw["co2_rate_signed_g_l_h"] = _sccm_to_g_l_h(raw["flow_filt_sccm"], 2.0)
            raw["co2_rate_raw_g_l_h"] = raw["co2_rate_signed_g_l_h"].clip(lower=0.0)
            raw["matrix"] = "natural"
            raw["batch"] = batch
            raw["experiment_code"] = batch
            raw["lot"] = "natural_lot1"
            raw["calibratable"] = False
            raw["chemistry_available"] = False
            raw["source"] = "repository raw sensor; assumed 2 L"
            rows.append(raw)
            first_h = float(raw["t_h"].min()) if not raw.empty else np.nan
            last_h = float(raw["t_h"].max()) if not raw.empty else np.nan
            n_raw = int(len(raw))
        else:
            first_h = np.nan
            last_h = np.nan
            n_raw = 0

        inventory.append(
            {
                "matrix": "natural",
                "batch": batch,
                "experiment_code": batch,
                "lot": "natural_lot1",
                "chemistry_available": False,
                "co2_available": bool(n_raw),
                "calibratable": False,
                "split": "descriptive_only",
                "n_co2_raw": n_raw,
                "co2_first_h": first_h,
                "co2_last_h": last_h,
                "volume_basis": "assumed_2L",
                "note": "CO2/temperature files exist; no homologated chemistry sheet is available for the upstream model.",
            }
        )
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(), pd.DataFrame(inventory))


def load_co2_data(
    natural_tables: dict[str, pd.DataFrame],
    chemical_support: pd.DataFrame,
    nutrient_pulses: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    frames: list[pd.DataFrame] = []

    lot1 = pd.read_csv(LOT1_CO2_PATH)
    lot1["batch"] = "lot1_" + lot1["process"].astype(str)
    lot1["matrix"] = "synthetic"
    lot1["experiment_code"] = lot1["batch"].map(SYNTHETIC_CODES)
    lot1["lot"] = "MBDOE_lot1"
    lot1["sensor_flow_sccm"] = pd.to_numeric(lot1["flow_filt_sccm"], errors="coerce")
    lot1["reactor_volume_l"] = pd.to_numeric(lot1["reactor_volume_l"], errors="coerce")
    lot1["calibratable"] = True
    lot1["chemistry_available"] = True
    lot1["source"] = str(LOT1_CO2_PATH.relative_to(ROOT_DIR))
    frames.append(lot1)

    lot2 = pd.read_csv(LOT2_CO2_PATH)
    lot2["batch"] = "lot2_" + lot2["process"].astype(str)
    lot2["matrix"] = "synthetic"
    lot2["experiment_code"] = lot2["batch"].map(SYNTHETIC_CODES)
    lot2["lot"] = "MBDOE_lot2"
    lot2_volume_l = pd.to_numeric(lot2["reactor_volume_ml_sampling_only"], errors="coerce") / 1000.0
    lot2["sensor_flow_sccm"] = pd.to_numeric(lot2["flow_filt_sccm"], errors="coerce")
    lot2["reactor_volume_l"] = lot2_volume_l
    lot2["calibratable"] = True
    lot2["chemistry_available"] = True
    lot2["source"] = str(LOT2_CO2_PATH.relative_to(ROOT_DIR))
    frames.append(lot2)

    natural = natural_tables["co2"].copy()
    natural["matrix"] = "natural"
    natural["experiment_code"] = natural["batch"].astype(str)
    natural["lot"] = "natural_2026"
    natural_volume_l = pd.to_numeric(natural["reactor_volume_l"], errors="coerce")
    natural["sensor_flow_sccm"] = pd.to_numeric(natural["flow_filt_sccm"], errors="coerce")
    natural["reactor_volume_l"] = natural_volume_l
    natural["calibratable"] = True
    natural["chemistry_available"] = True
    natural["source"] = "homologated natural workbook + repository process sensors"
    frames.append(natural)

    keep = [
        "matrix",
        "batch",
        "experiment_code",
        "lot",
        "t_h",
        "sensor_flow_sccm",
        "reactor_volume_l",
        "calibratable",
        "chemistry_available",
        "source",
    ]
    raw = pd.concat([frame[keep] for frame in frames], ignore_index=True)
    raw["t_h"] = pd.to_numeric(raw["t_h"], errors="coerce")
    raw["sensor_flow_sccm"] = pd.to_numeric(raw["sensor_flow_sccm"], errors="coerce")
    raw["reactor_volume_l"] = pd.to_numeric(raw["reactor_volume_l"], errors="coerce")
    raw = raw.dropna(subset=["t_h", "sensor_flow_sccm", "reactor_volume_l"])
    raw = raw[raw["t_h"].ge(0.0)].copy()
    raw = raw[~raw["batch"].isin(EXCLUDED_BATCHES)].copy()
    source_coverage = (
        raw.groupby("batch", as_index=False)
        .agg(
            co2_source_first_h=("t_h", "min"),
            co2_source_last_h=("t_h", "max"),
            n_co2_source_raw=("t_h", "size"),
        )
    )
    raw = raw.merge(
        chemical_support[
            ["batch", "chemistry_first_h", "chemistry_last_h", "n_chemical_timepoints"]
        ],
        on="batch",
        how="inner",
        validate="many_to_one",
    ).merge(source_coverage, on="batch", how="left", validate="many_to_one")
    raw = raw[
        raw["t_h"].ge(raw["chemistry_first_h"] - 1e-9)
        & raw["t_h"].le(raw["chemistry_last_h"] + 1e-9)
    ].copy()
    raw = _attach_temperature_context(raw, natural_tables)
    raw, sensor_zero_offsets = apply_sensor_zero_correction(raw)
    sampling_schedule = load_sampling_schedule()
    sampling_schedule = sampling_schedule[sampling_schedule["batch"].isin(raw["batch"].unique())].copy()
    sampling_schedule = sampling_schedule.merge(
        chemical_support[["batch", "chemistry_first_h", "chemistry_last_h"]],
        on="batch",
        how="left",
        validate="many_to_one",
    )
    sampling_schedule = sampling_schedule[
        sampling_schedule["sample_time_h"].ge(sampling_schedule["chemistry_first_h"] - 1e-9)
        & sampling_schedule["sample_time_h"].le(sampling_schedule["chemistry_last_h"] + 1e-9)
    ].drop(columns=["chemistry_first_h", "chemistry_last_h"])
    active_pulses = nutrient_pulses[
        nutrient_pulses["batch"].isin(raw["batch"].unique())
    ].copy()
    raw = filter_co2_sensor_artifacts(raw, sampling_schedule, active_pulses)
    raw["time_bin_h"] = np.round(raw["t_h"] / OBSERVATION_BIN_H) * OBSERVATION_BIN_H
    hourly = (
        raw.groupby(
            [
                "matrix",
                "batch",
                "experiment_code",
                "lot",
                "calibratable",
                "chemistry_available",
                "source",
                "chemistry_first_h",
                "chemistry_last_h",
                "n_chemical_timepoints",
                "co2_source_first_h",
                "co2_source_last_h",
                "n_co2_source_raw",
                "time_bin_h",
            ],
            as_index=False,
            dropna=False,
        )
        .agg(
            t_h=("t_h", "median"),
            sensor_flow_sccm=("sensor_flow_sccm", "median"),
            sensor_flow_zero_corrected_sccm=("sensor_flow_zero_corrected_sccm", "median"),
            sensor_zero_offset_sccm=("sensor_zero_offset_sccm", "first"),
            sensor_id=("sensor_id", "first"),
            acquisition_channel=("acquisition_channel", "first"),
            co2_rate_signed_g_l_h=("co2_rate_signed_g_l_h", "median"),
            co2_rate_uncorrected_physical_g_l_h=(
                "co2_rate_uncorrected_physical_g_l_h",
                "median",
            ),
            co2_rate_zero_corrected_signed_g_l_h=(
                "co2_rate_zero_corrected_signed_g_l_h",
                "median",
            ),
            co2_rate_raw_g_l_h=("co2_rate_raw_g_l_h", "median"),
            co2_rate_filtered_g_l_h=("co2_rate_filtered_g_l_h", "median"),
            setpoint_c=("setpoint_c", "median"),
            sensor_temperature_c=("sensor_temperature_c", "median"),
            n_raw=("t_h", "size"),
            n_artifact=("artifact_flag", "sum"),
            artifact_fraction=("artifact_flag", "mean"),
            sampling_artifact_fraction=(
                "artifact_reason",
                lambda values: float(pd.Series(values).eq("sampling_window_transient").mean()),
            ),
            cold_operation_fraction=("cold_operation", "mean"),
            cold_low_sensitivity_fraction=("cold_low_sensitivity", "mean"),
            pulse_response_protected_fraction=("pulse_response_protected", "mean"),
        )
        .sort_values(["matrix", "batch", "t_h"])
        .reset_index(drop=True)
    )
    hourly = smooth_hourly_co2_profiles(hourly, active_pulses)
    hourly["co2_rate_g_l_h"] = hourly["co2_rate_smoothed_g_l_h"]
    hourly["cold_operation"] = hourly["setpoint_c"].le(COLD_SETPOINT_THRESHOLD_C)
    hourly["detection_limit_g_l_h"] = np.where(
        hourly["cold_operation"],
        COLD_CO2_DETECTION_LIMIT_G_L_H,
        CO2_DETECTION_LIMIT_G_L_H,
    )
    hourly["left_censored"] = hourly["co2_rate_smoothed_g_l_h"].le(
        hourly["detection_limit_g_l_h"]
    )
    hourly["cold_low_sensitivity"] = hourly["cold_operation"] & hourly["left_censored"]

    qc_rows = []
    for (matrix, batch), group in raw.groupby(["matrix", "batch"], sort=True):
        hgroup = hourly[hourly["batch"].eq(batch)]
        qc_rows.append(
            {
                "matrix": str(matrix),
                "batch": str(batch),
                "experiment_code": SYNTHETIC_CODES.get(str(batch), str(batch)),
                "acquisition_channel": str(group["acquisition_channel"].iloc[0]),
                "sensor_id": int(group["sensor_id"].iloc[0]),
                "sensor_zero_offset_sccm": float(group["sensor_zero_offset_sccm"].iloc[0]),
                "sensor_zero_offset_g_l_h": float(
                    group["co2_rate_uncorrected_physical_g_l_h"].sub(
                        group["co2_rate_raw_g_l_h"]
                    ).median()
                ),
                "n_raw": int(len(group)),
                "negative_signed_fraction": float(group["co2_rate_signed_g_l_h"].lt(0.0).mean()),
                "n_artifacts_replaced": int(group["artifact_flag"].sum()),
                "artifact_fraction": float(group["artifact_flag"].mean()),
                "n_sampling_window_transients": int(group["artifact_reason"].eq("sampling_window_transient").sum()),
                "n_other_short_transients": int(group["artifact_reason"].eq("short_transient_outlier").sum()),
                "median_setpoint_c": float(group["setpoint_c"].median()) if group["setpoint_c"].notna().any() else np.nan,
                "cold_operation_fraction": float(group["cold_operation"].mean()),
                "n_hourly": int(len(hgroup)),
                "n_left_censored_hourly": int(hgroup["left_censored"].sum()),
                "left_censored_fraction_hourly": float(hgroup["left_censored"].mean()),
                "n_cold_low_sensitivity_hourly": int(hgroup["cold_low_sensitivity"].sum()),
                "smoothing_integral_ratio": float(
                    np.trapz(hgroup["co2_rate_smoothed_g_l_h"], hgroup["t_h"])
                    / max(np.trapz(hgroup["co2_rate_filtered_g_l_h"], hgroup["t_h"]), 1e-12)
                ),
                "smoothing_peak_ratio": float(
                    hgroup["co2_rate_smoothed_g_l_h"].max()
                    / max(hgroup["co2_rate_filtered_g_l_h"].max(), 1e-12)
                ),
                "smoothing_roughness_ratio": float(
                    np.mean(np.abs(np.diff(hgroup["co2_rate_smoothed_g_l_h"], n=2)))
                    / max(np.mean(np.abs(np.diff(hgroup["co2_rate_filtered_g_l_h"], n=2))), 1e-12)
                ) if len(hgroup) >= 3 else np.nan,
            }
        )
    qc_summary = pd.DataFrame(qc_rows)

    inventory_rows = []
    for (matrix, batch), group in raw.groupby(["matrix", "batch"], sort=True):
        calibratable = bool(group["calibratable"].iloc[0])
        split = "validation" if batch == HOLDOUTS.get(str(matrix)) else ("calibration" if calibratable else "descriptive_only")
        inventory_rows.append(
            {
                "matrix": str(matrix),
                "batch": str(batch),
                "experiment_code": str(group["experiment_code"].iloc[0]),
                "lot": str(group["lot"].iloc[0]),
                "acquisition_channel": str(group["acquisition_channel"].iloc[0]),
                "sensor_id": int(group["sensor_id"].iloc[0]),
                "sensor_zero_offset_sccm": float(group["sensor_zero_offset_sccm"].iloc[0]),
                "chemistry_available": bool(group["chemistry_available"].iloc[0]),
                "co2_available": True,
                "calibratable": calibratable,
                "split": split,
                "n_co2_raw": int(len(group)),
                "n_co2_source_raw": int(group["n_co2_source_raw"].iloc[0]),
                "n_co2_hourly": int(hourly["batch"].eq(batch).sum()),
                "co2_first_h": float(group["t_h"].min()),
                "co2_last_h": float(group["t_h"].max()),
                "co2_source_first_h": float(group["co2_source_first_h"].iloc[0]),
                "co2_source_last_h": float(group["co2_source_last_h"].iloc[0]),
                "chemistry_first_h": float(group["chemistry_first_h"].iloc[0]),
                "chemistry_last_h": float(group["chemistry_last_h"].iloc[0]),
                "n_chemical_timepoints": int(group["n_chemical_timepoints"].iloc[0]),
                "co2_start_minus_chemistry_h": float(
                    group["t_h"].min() - group["chemistry_first_h"].iloc[0]
                ),
                "co2_end_minus_chemistry_h": float(
                    group["t_h"].max() - group["chemistry_last_h"].iloc[0]
                ),
                "co2_peak_g_l_h": float(
                    hourly.loc[hourly["batch"].eq(batch), "co2_rate_smoothed_g_l_h"].max()
                ),
                "median_setpoint_c": float(group["setpoint_c"].median()) if group["setpoint_c"].notna().any() else np.nan,
                "n_artifacts_replaced": int(group["artifact_flag"].sum()),
                "n_left_censored_hourly": int(hourly.loc[hourly["batch"].eq(batch), "left_censored"].sum()),
                "volume_basis": "sampling-corrected" if matrix == "synthetic" else "workbook/default reactor volume",
                "note": "",
            }
        )
    inventory = pd.DataFrame(inventory_rows)
    return hourly, inventory, raw, qc_summary, sampling_schedule, sensor_zero_offsets


def _chemical_activity_bracket(batch: base.BatchData) -> dict[str, float]:
    """Bracket fermentation activation using chemistry only, never CO2.

    Activity is first supported when glucose+fructose has fallen by 5 g/L or
    ethanol has risen by 2 g/L relative to the first finite chemical sample.
    The previous chemical time is the lower bracket. A logistic transition is
    fixed to 10% at the lower bound and 90% at the upper bound.
    """

    time_h = np.asarray(batch.time, dtype=float)
    glucose = np.asarray(batch.observations.get("G", np.full(len(time_h), np.nan)), dtype=float)
    fructose = np.asarray(batch.observations.get("F", np.full(len(time_h), np.nan)), dtype=float)
    ethanol = np.asarray(batch.observations.get("E", np.full(len(time_h), np.nan)), dtype=float)
    sugar = glucose + fructose
    chemical = np.isfinite(sugar) | np.isfinite(ethanol)
    chemical_indices = np.flatnonzero(chemical)
    if not len(chemical_indices):
        return {"lower_h": 0.0, "upper_h": 0.0, "center_h": 0.0, "scale_h": 1.0}
    sugar_values = sugar[np.isfinite(sugar)]
    ethanol_values = ethanol[np.isfinite(ethanol)]
    sugar_initial = float(sugar_values[0]) if len(sugar_values) else np.nan
    ethanol_initial = float(ethanol_values[0]) if len(ethanol_values) else np.nan
    active = (
        (np.isfinite(sugar) & np.isfinite(sugar_initial) & (sugar <= sugar_initial - CHEMISTRY_SUGAR_DROP_G_L))
        | (
            np.isfinite(ethanol)
            & np.isfinite(ethanol_initial)
            & (ethanol >= ethanol_initial + CHEMISTRY_ETHANOL_RISE_G_L)
        )
    )
    active_indices = np.flatnonzero(active)
    if not len(active_indices):
        lower_h = float(time_h[chemical_indices[0]])
        upper_h = lower_h
    else:
        upper_index = int(active_indices[0])
        earlier = chemical_indices[chemical_indices < upper_index]
        lower_index = int(earlier[-1]) if len(earlier) else upper_index
        lower_h = float(time_h[lower_index])
        upper_h = float(time_h[upper_index])
    center_h = 0.5 * (lower_h + upper_h)
    scale_h = max(
        CHEMISTRY_ACTIVATION_MIN_SCALE_H,
        (upper_h - lower_h) / (2.0 * math.log(9.0)) if upper_h > lower_h else 1.0,
    )
    return {
        "lower_h": lower_h,
        "upper_h": upper_h,
        "center_h": center_h,
        "scale_h": scale_h,
    }


def _logistic_activation(time_h: np.ndarray, center_h: float, scale_h: float) -> np.ndarray:
    z = np.clip((np.asarray(time_h, dtype=float) - float(center_h)) / max(float(scale_h), 1e-9), -50, 50)
    return 1.0 / (1.0 + np.exp(-z))


def _bounded_smoothstep_activation(
    time_h: np.ndarray,
    lower_h: float,
    upper_h: float,
    start_fraction: float = 0.0,
    duration_fraction: float = 1.0,
) -> np.ndarray:
    """Chemistry-anchored gradual activation shared at matrix level."""

    time_h = np.asarray(time_h, dtype=float)
    if upper_h <= lower_h + 1e-9:
        return (time_h >= float(upper_h)).astype(float)
    bracket_h = float(upper_h) - float(lower_h)
    start_h = float(lower_h) + float(start_fraction) * bracket_h
    duration_h = max(float(duration_fraction) * bracket_h, DRIVER_GRID_H)
    fraction = np.clip(
        (time_h - start_h) / duration_h,
        0.0,
        1.0,
    )
    return fraction * fraction * (3.0 - 2.0 * fraction)


def build_driver_cache(
    batches: dict[str, base.BatchData],
    theta_by_matrix: dict[str, dict[str, float]],
    observations: pd.DataFrame,
) -> tuple[dict[str, DriverCache], pd.DataFrame]:
    cache: dict[str, DriverCache] = {}
    diagnostic_rows: list[dict[str, object]] = []
    usable = observations[observations["calibratable"].astype(bool)].copy()

    for batch_name, group in usable.groupby("batch", sort=True):
        if batch_name not in batches:
            raise RuntimeError(f"No BatchData object is available for {batch_name}")
        batch = batches[batch_name]
        matrix = str(group["matrix"].iloc[0])
        theta = theta_by_matrix[matrix]
        horizon = float(group["chemistry_last_h"].iloc[0])
        time_h = np.arange(0.0, horizon + DRIVER_GRID_H + 1e-9, DRIVER_GRID_H)
        time_h = time_h[time_h <= horizon + 1e-9]
        if not np.isclose(time_h[-1], horizon):
            time_h = np.append(time_h, horizon)
        core = base.simulate(batch, theta, time_h)
        if core is None:
            raise RuntimeError(f"The fixed upstream model failed for {batch_name}")

        base_qprod = co2_model.co2_production_g_l_h(theta, batch, core, time_h)
        n_pulses = tuple(batch.pulses.get("N", tuple()))
        if len(n_pulses) > 1:
            raise RuntimeError(
                f"The nitrogen-boost transition currently supports one in-process N pulse; "
                f"{batch_name} has {len(n_pulses)}"
            )
        pulses_without_n = dict(batch.pulses)
        pulses_without_n["N"] = tuple()
        batch_without_n_pulse = replace(batch, pulses=pulses_without_n)
        core_without_n_pulse = base.simulate(batch_without_n_pulse, theta, time_h)
        if core_without_n_pulse is None:
            raise RuntimeError(f"The no-N-pulse counterfactual failed for {batch_name}")
        qprod_without_n_pulse = co2_model.co2_production_g_l_h(
            theta, batch_without_n_pulse, core_without_n_pulse, time_h
        )
        n_pulse_time_h = float(n_pulses[0][0]) if n_pulses else np.nan
        n_pulse_amount_kg_m3 = float(n_pulses[0][1]) if n_pulses else 0.0
        biomass = np.maximum(core["X"].to_numpy(dtype=float), 0.0)
        biomass_without_n_pulse = np.maximum(
            core_without_n_pulse["X"].to_numpy(dtype=float), 0.0
        )
        temperature = np.asarray([base.temperature_at(batch, float(t)) for t in time_h], dtype=float)
        csat_base = (
            1.69
            * np.exp(-0.032 * (temperature - 20.0))
            * np.exp(0.0016 * np.maximum(core["E"].to_numpy(dtype=float), 0.0))
            * np.exp(
                -0.0012
                * np.maximum(core["G"].to_numpy(dtype=float) + core["F"].to_numpy(dtype=float), 0.0)
            )
        )
        o2_saturation_base = np.asarray(
            [
                co2_model.o2_saturation_mg_l(
                    float(temp), float(ethanol), float(glucose), float(fructose), 1.0
                )
                for temp, ethanol, glucose, fructose in zip(
                    temperature,
                    core["E"].to_numpy(dtype=float),
                    core["G"].to_numpy(dtype=float),
                    core["F"].to_numpy(dtype=float),
                )
            ],
            dtype=float,
        )
        activity = _chemical_activity_bracket(batch)
        previous_activation = _logistic_activation(
            time_h, activity["center_h"], activity["scale_h"]
        )
        activation = _bounded_smoothstep_activation(
            time_h, activity["lower_h"], activity["upper_h"]
        )
        theta_source = (
            "estimability_historical_natural/theta.csv"
            if matrix == "natural"
            else "estimability_historical_synthetic_plus_lot2/theta_by_case.csv"
        )
        cache[batch_name] = DriverCache(
            batch=batch_name,
            matrix=matrix,
            time_h=time_h,
            base_qprod_g_l_h=np.asarray(base_qprod, dtype=float),
            base_qprod_no_n_pulse_g_l_h=np.asarray(qprod_without_n_pulse, dtype=float),
            n_pulse_qprod_increment_g_l_h=np.asarray(
                base_qprod - qprod_without_n_pulse, dtype=float
            ),
            biomass_g_l=biomass,
            biomass_no_n_pulse_g_l=biomass_without_n_pulse,
            n_pulse_biomass_increment_g_l=biomass - biomass_without_n_pulse,
            n_pulse_time_h=n_pulse_time_h,
            n_pulse_amount_kg_m3=n_pulse_amount_kg_m3,
            o2_saturation_base_mg_l=o2_saturation_base,
            chemical_activation=activation,
            chemical_activation_previous=previous_activation,
            chemical_activity_lower_h=float(activity["lower_h"]),
            chemical_activity_upper_h=float(activity["upper_h"]),
            chemical_activity_center_h=float(activity["center_h"]),
            chemical_activity_scale_h=float(activity["scale_h"]),
            csat_base_g_l=np.asarray(csat_base, dtype=float),
            temperature_used_c=np.asarray(temperature, dtype=float),
            observation_time_h=group["t_h"].to_numpy(dtype=float),
            theta_source=theta_source,
        )
        diagnostic_rows.append(
            {
                "matrix": matrix,
                "batch": batch_name,
                "theta_source": theta_source,
                "driver_grid_h": DRIVER_GRID_H,
                "driver_horizon_h": horizon,
                "base_qprod_peak_g_l_h": float(np.max(base_qprod)),
                "base_qprod_no_n_pulse_peak_g_l_h": float(np.max(qprod_without_n_pulse)),
                "n_pulse_qprod_increment_peak_g_l_h": float(
                    np.max(base_qprod - qprod_without_n_pulse)
                ),
                "n_pulse_time_h": n_pulse_time_h,
                "n_pulse_amount_kg_m3": n_pulse_amount_kg_m3,
                "chemical_activity_lower_h": float(activity["lower_h"]),
                "chemical_activity_upper_h": float(activity["upper_h"]),
                "chemical_activity_center_h": float(activity["center_h"]),
                "chemical_activity_scale_h": float(activity["scale_h"]),
                "initial_o2_saturation_base_mg_l": float(o2_saturation_base[0]),
                "csat_base_min_g_l": float(np.min(csat_base)),
                "csat_base_max_g_l": float(np.max(csat_base)),
                "temperature_used_min_c": float(np.min(temperature)),
                "temperature_used_max_c": float(np.max(temperature)),
                "temperature_used_mean_c": float(np.mean(temperature)),
                "core_E_end_g_l": float(core["E"].iloc[-1]),
                "core_G_end_g_l": float(core["G"].iloc[-1]),
                "core_F_end_g_l": float(core["F"].iloc[-1]),
            }
        )
    return cache, pd.DataFrame(diagnostic_rows)


def temperature_model_input_table(cache: dict[str, DriverCache]) -> pd.DataFrame:
    rows = []
    for batch_name, item in cache.items():
        role = "validation" if batch_name == HOLDOUTS[item.matrix] else "calibration"
        for time_h, temperature_c, csat_base in zip(
            item.time_h, item.temperature_used_c, item.csat_base_g_l
        ):
            rows.append(
                {
                    "matrix": item.matrix,
                    "batch": batch_name,
                    "experiment_code": SYNTHETIC_CODES.get(batch_name, batch_name),
                    "role": role,
                    "time_h": float(time_h),
                    "temperature_used_c": float(temperature_c),
                    "csat_base_g_l": float(csat_base),
                }
            )
    return pd.DataFrame(rows).sort_values(["matrix", "batch", "time_h"]).reset_index(drop=True)


def temperature_alignment_summary(
    observations: pd.DataFrame,
    temperature_inputs: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for batch, group in observations.groupby("batch", sort=True):
        model = temperature_inputs[temperature_inputs["batch"].eq(batch)].sort_values("time_h")
        valid = group["sensor_temperature_c"].notna().to_numpy()
        used = np.interp(
            group["t_h"].to_numpy(dtype=float),
            model["time_h"].to_numpy(dtype=float),
            model["temperature_used_c"].to_numpy(dtype=float),
        )
        error = used[valid] - group.loc[valid, "sensor_temperature_c"].to_numpy(dtype=float)
        rows.append(
            {
                "matrix": str(group["matrix"].iloc[0]),
                "batch": str(batch),
                "experiment_code": SYNTHETIC_CODES.get(str(batch), str(batch)),
                "n_temperature_comparisons": int(len(error)),
                "model_minus_sensor_bias_c": float(np.mean(error)),
                "model_sensor_mae_c": float(np.mean(np.abs(error))),
                "model_sensor_rmse_c": float(np.sqrt(np.mean(error * error))),
                "model_sensor_max_abs_c": float(np.max(np.abs(error))),
            }
        )
    return pd.DataFrame(rows)


def _smooth_positive(value: float, smooth: float = 0.015) -> float:
    ratio = float(value) / float(smooth)
    if ratio > 50.0:
        return float(value)
    if ratio < -50.0:
        return 0.0
    return float(smooth) * math.log1p(math.exp(ratio))


def _shape_parameter_names(model_variant: str) -> tuple[str, ...]:
    if model_variant == LEGACY_MODEL_NAME:
        return ("kCO2_release_h", "CO2sat_scale")
    if model_variant == CHEMISTRY_ALIGNED_MODEL_NAME:
        return (
            "kCO2_release_h",
            "CO2sat_scale",
            "O2_qmax_mg_gdw_h",
            "O2_initial_scale",
        )
    if model_variant == MODEL_NAME:
        return (
            "kCO2_release_h",
            "CO2sat_scale",
            "O2_qmax_mg_gdw_h",
            "O2_initial_scale",
            "pulse_t_rise_h",
            "pulse_activity_gain",
            "chem_activation_start_fraction",
            "chem_activation_duration_fraction",
        )
    if model_variant == THRESHOLD_RELEASE_MODEL_NAME:
        return (
            "kCO2_release_h",
            "CO2sat_scale",
            "O2_qmax_mg_gdw_h",
            "O2_initial_scale",
            "pulse_t_rise_h",
            "pulse_activity_gain",
        )
    raise ValueError(f"Unknown model variant: {model_variant}")


def _shape_values(log_shape: np.ndarray, model_variant: str) -> dict[str, float]:
    names = _shape_parameter_names(model_variant)
    values = np.exp(np.asarray(log_shape, dtype=float))
    parameters = dict(zip(names, (float(value) for value in values)))
    if model_variant == LEGACY_MODEL_NAME:
        parameters.update({"O2_qmax_mg_gdw_h": 0.15, "O2_initial_scale": 1.0})
    if model_variant not in {MODEL_NAME, THRESHOLD_RELEASE_MODEL_NAME}:
        parameters.update({"pulse_t_rise_h": 1.0, "pulse_activity_gain": 1.0})
    if model_variant != MODEL_NAME:
        parameters.update(
            {
                "chem_activation_start_fraction": 0.0,
                "chem_activation_duration_fraction": 1.0,
            }
        )
    return parameters


def _uses_chemistry_alignment(model_variant: str) -> bool:
    return model_variant in {
        MODEL_NAME,
        THRESHOLD_RELEASE_MODEL_NAME,
        CHEMISTRY_ALIGNED_MODEL_NAME,
    }


def _uses_nitrogen_boost_transition(model_variant: str) -> bool:
    return model_variant in {MODEL_NAME, THRESHOLD_RELEASE_MODEL_NAME}


def _uses_continuous_release(model_variant: str) -> bool:
    return model_variant == MODEL_NAME


def _uses_bounded_chemical_activation(model_variant: str) -> bool:
    return model_variant == MODEL_NAME


def _pulse_utilization_fraction(cache: DriverCache, pulse_t_rise_h: float) -> np.ndarray:
    """Fraction of the added nitrogen physiologically available after the pulse.

    David et al. distribute the nitrogen input uniformly from ``t_add`` to
    ``t_add + t_rise``.  The resulting cumulative availability is therefore a
    causal linear ramp.  Batches without an N pulse return zeros because their
    counterfactual increment is also zero.
    """

    if not np.isfinite(cache.n_pulse_time_h) or cache.n_pulse_amount_kg_m3 <= 0.0:
        return np.zeros(len(cache.time_h), dtype=float)
    elapsed = np.maximum(cache.time_h - float(cache.n_pulse_time_h), 0.0)
    return np.clip(elapsed / max(float(pulse_t_rise_h), 1e-9), 0.0, 1.0)


def effective_qprod_grid(
    cache: DriverCache,
    o2_qmax_mg_gdw_h: float,
    o2_initial_scale: float,
    chemistry_aligned: bool,
    nitrogen_boost_transition: bool = False,
    pulse_t_rise_h: float = 1.0,
    pulse_activity_gain: float = 1.0,
    bounded_chemical_activation: bool = False,
    chem_activation_start_fraction: float = 0.0,
    chem_activation_duration_fraction: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute O2-gated production on the cached upstream trajectory.

    The nitrogen-boost variant replaces the upstream model's instantaneous
    N-pulse effect by a causal finite utilization ramp. The biomass increment
    follows the unscaled ramp. ``pulse_activity_gain`` multiplies the metabolic
    activity already present in the no-pulse counterfactual after the nitrogen
    becomes available; this permits a transport/activity boost without requiring
    a biomass increment, as reported for stationary-phase additions.
    """

    time_h = cache.time_h
    if nitrogen_boost_transition:
        pulse_fraction = _pulse_utilization_fraction(cache, pulse_t_rise_h)
        biomass = np.maximum(
            cache.biomass_no_n_pulse_g_l
            + pulse_fraction * cache.n_pulse_biomass_increment_g_l,
            0.0,
        )
        activity_multiplier = 1.0 + (
            float(pulse_activity_gain) - 1.0
        ) * pulse_fraction
        biological_qprod = np.maximum(
            activity_multiplier * cache.base_qprod_no_n_pulse_g_l_h
            + pulse_fraction * cache.n_pulse_qprod_increment_g_l_h,
            0.0,
        )
    else:
        biomass = cache.biomass_g_l
        biological_qprod = np.maximum(cache.base_qprod_g_l_h, 0.0)
    o2 = np.zeros(len(time_h), dtype=float)
    uptake = np.zeros(len(time_h), dtype=float)
    o2[0] = max(float(o2_initial_scale), 0.0) * float(cache.o2_saturation_base_mg_l[0])
    qmax = max(float(o2_qmax_mg_gdw_h), 1e-12)
    for idx in range(len(time_h)):
        uptake[idx] = (
            qmax
            * max(float(biomass[idx]), 0.0)
            * o2[idx]
            / (O2_K_MG_L + o2[idx])
        )
        if idx < len(time_h) - 1:
            dt = max(float(time_h[idx + 1] - time_h[idx]), 1e-9)
            o2[idx + 1] = max(o2[idx] - dt * uptake[idx], 0.0)
    phi_ana = (O2_ANA_K_MG_L**O2_ANA_HILL) / (
        O2_ANA_K_MG_L**O2_ANA_HILL + np.maximum(o2, 0.0) ** O2_ANA_HILL
    )
    ferment_fraction = O2_CRABTREE_FLOOR + (1.0 - O2_CRABTREE_FLOOR) * phi_ana
    respiratory_co2 = co2_model.CO2_G_PER_MG_O2_RESP * np.maximum(uptake, 0.0)
    if chemistry_aligned:
        activation = (
            _bounded_smoothstep_activation(
                time_h,
                cache.chemical_activity_lower_h,
                cache.chemical_activity_upper_h,
                chem_activation_start_fraction,
                chem_activation_duration_fraction,
            )
            if bounded_chemical_activation
            else cache.chemical_activation_previous
        )
    else:
        activation = np.ones(len(time_h), dtype=float)
    qprod = activation * (
        biological_qprod * ferment_fraction + respiratory_co2
    )
    return np.maximum(qprod, 0.0), o2, phi_ana


def raw_qgas_grid_prediction(
    cache: DriverCache,
    k_release_h: float,
    sat_scale: float,
    o2_qmax_mg_gdw_h: float,
    o2_initial_scale: float,
    chemistry_aligned: bool,
    nitrogen_boost_transition: bool = False,
    pulse_t_rise_h: float = 1.0,
    pulse_activity_gain: float = 1.0,
    continuous_release: bool = False,
    bounded_chemical_activation: bool = False,
    chem_activation_start_fraction: float = 0.0,
    chem_activation_duration_fraction: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    time_h = cache.time_h
    qprod, o2, phi_ana = effective_qprod_grid(
        cache,
        o2_qmax_mg_gdw_h,
        o2_initial_scale,
        chemistry_aligned,
        nitrogen_boost_transition,
        pulse_t_rise_h,
        pulse_activity_gain,
        bounded_chemical_activation,
        chem_activation_start_fraction,
        chem_activation_duration_fraction,
    )
    csat = float(sat_scale) * cache.csat_base_g_l
    dissolved = 0.0
    qgas = np.zeros(len(time_h), dtype=float)
    for idx in range(1, len(time_h)):
        dt = max(float(time_h[idx] - time_h[idx - 1]), 1e-9)
        available = dissolved / dt + max(float(qprod[idx - 1]), 0.0)
        if continuous_release:
            # Open/off-gas operation does not require the liquid to reach the
            # pure-CO2 saturation concentration before any CO2 can leave.  A
            # small linear path is therefore present from the first dissolved
            # CO2, while the transfer fraction increases smoothly as the pool
            # approaches its solubility scale. This removes the artificial
            # zero plateau without introducing a direct-feed shortcut.
            dissolved_fraction = dissolved / max(dissolved + float(csat[idx]), 1e-12)
            release_fraction = CO2_CONTINUOUS_RELEASE_FLOOR + (
                1.0 - CO2_CONTINUOUS_RELEASE_FLOOR
            ) * dissolved_fraction
            potential_release = (
                max(float(k_release_h), 1e-9) * dissolved * release_fraction
            )
        else:
            excess = _smooth_positive(dissolved - float(csat[idx]))
            potential_release = max(float(k_release_h), 1e-9) * excess
        qgas[idx] = min(potential_release, available)
        dissolved = max(dissolved + dt * (max(float(qprod[idx - 1]), 0.0) - qgas[idx]), 0.0)
    return qgas, qprod, o2, phi_ana


def raw_qgas_prediction(
    cache: DriverCache,
    k_release_h: float,
    sat_scale: float,
    o2_qmax_mg_gdw_h: float = 0.15,
    o2_initial_scale: float = 1.0,
    chemistry_aligned: bool = True,
    nitrogen_boost_transition: bool = False,
    pulse_t_rise_h: float = 1.0,
    pulse_activity_gain: float = 1.0,
    continuous_release: bool = False,
    bounded_chemical_activation: bool = False,
    chem_activation_start_fraction: float = 0.0,
    chem_activation_duration_fraction: float = 1.0,
) -> np.ndarray:
    qgas, _, _, _ = raw_qgas_grid_prediction(
        cache,
        k_release_h,
        sat_scale,
        o2_qmax_mg_gdw_h,
        o2_initial_scale,
        chemistry_aligned,
        nitrogen_boost_transition,
        pulse_t_rise_h,
        pulse_activity_gain,
        continuous_release,
        bounded_chemical_activation,
        chem_activation_start_fraction,
        chem_activation_duration_fraction,
    )
    return np.interp(cache.observation_time_h, cache.time_h, qgas)


def _onset_threshold(group: pd.DataFrame) -> float:
    early = group[group["t_h"].le(ONSET_BASELINE_WINDOW_H)]
    baseline = float(early["co2_rate_g_l_h"].median()) if not early.empty else 0.0
    peak = float(group["co2_rate_g_l_h"].max())
    early_limit = (
        float(early["detection_limit_g_l_h"].median())
        if not early.empty
        else CO2_DETECTION_LIMIT_G_L_H
    )
    return max(early_limit, baseline + ONSET_DYNAMIC_RANGE_FRACTION * max(peak - baseline, 0.0))


def _sustained_onset_h(
    time_h: np.ndarray,
    values: np.ndarray,
    threshold: float,
    consecutive: int = 3,
) -> float:
    time_h = np.asarray(time_h, dtype=float)
    values = np.asarray(values, dtype=float)
    above = np.isfinite(values) & (values >= float(threshold))
    for idx in range(max(len(above) - consecutive + 1, 0)):
        if not above[idx : idx + consecutive].all():
            continue
        if idx == 0 or not np.isfinite(values[idx - 1]) or values[idx] <= values[idx - 1]:
            return float(time_h[idx])
        fraction = (float(threshold) - values[idx - 1]) / max(values[idx] - values[idx - 1], 1e-12)
        return float(time_h[idx - 1] + np.clip(fraction, 0.0, 1.0) * (time_h[idx] - time_h[idx - 1]))
    return np.nan


def _postpulse_peak_time_h(
    time_h: np.ndarray,
    values: np.ndarray,
    pulse_time_h: float,
    window_h: float = PULSE_PEAK_WINDOW_H,
) -> float:
    """Return the post-pulse peak time on a fixed causal comparison window."""

    if not np.isfinite(pulse_time_h):
        return np.nan
    time_h = np.asarray(time_h, dtype=float)
    values = np.asarray(values, dtype=float)
    mask = (
        np.isfinite(time_h)
        & np.isfinite(values)
        & (time_h >= float(pulse_time_h))
        & (time_h <= float(pulse_time_h) + float(window_h))
    )
    if int(mask.sum()) < 3:
        return np.nan
    selected_time = time_h[mask]
    selected_values = values[mask]
    return float(selected_time[int(np.argmax(selected_values))])


def _initial_rise_times_h(
    time_h: np.ndarray,
    values: np.ndarray,
    low_threshold: float,
    high_threshold: float,
) -> tuple[float, float, float]:
    """Return sustained low/high crossings and their duration for the early rise."""

    low_h = _sustained_onset_h(time_h, values, low_threshold, consecutive=2)
    high_h = _sustained_onset_h(time_h, values, high_threshold, consecutive=3)
    duration_h = (
        high_h - low_h if np.isfinite(low_h) and np.isfinite(high_h) else np.nan
    )
    return low_h, high_h, duration_h


def _profile_matrix_gain(
    log_shape: np.ndarray,
    batch_names: list[str],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    model_variant: str = MODEL_NAME,
) -> tuple[float, dict[str, np.ndarray]]:
    parameters = _shape_values(log_shape, model_variant)
    numerator = 0.0
    denominator = 0.0
    raw_by_batch: dict[str, np.ndarray] = {}
    for batch_name in batch_names:
        group = observations[observations["batch"].eq(batch_name)]
        prediction = raw_qgas_prediction(
            cache[batch_name],
            parameters["kCO2_release_h"],
            parameters["CO2sat_scale"],
            parameters["O2_qmax_mg_gdw_h"],
            parameters["O2_initial_scale"],
            chemistry_aligned=_uses_chemistry_alignment(model_variant),
            nitrogen_boost_transition=_uses_nitrogen_boost_transition(model_variant),
            pulse_t_rise_h=parameters["pulse_t_rise_h"],
            pulse_activity_gain=parameters["pulse_activity_gain"],
            continuous_release=_uses_continuous_release(model_variant),
            bounded_chemical_activation=_uses_bounded_chemical_activation(model_variant),
            chem_activation_start_fraction=parameters["chem_activation_start_fraction"],
            chem_activation_duration_fraction=parameters["chem_activation_duration_fraction"],
        )
        raw_by_batch[batch_name] = prediction
        quantifiable = ~group["left_censored"].astype(bool).to_numpy()
        if quantifiable.any():
            sigma = max(0.04, 0.10 * float(group.loc[~group["left_censored"], "co2_rate_g_l_h"].max()))
            weight = 1.0 / (sigma * math.sqrt(len(group)))
            observed = group["co2_rate_g_l_h"].to_numpy(dtype=float)[quantifiable]
            prediction_quantifiable = prediction[quantifiable]
            numerator += float(np.dot(weight * prediction_quantifiable, weight * observed))
            denominator += float(np.dot(weight * prediction_quantifiable, weight * prediction_quantifiable))
    gain = numerator / denominator if denominator > 1e-16 else PARAMETER_BOUNDS["matrix_gain"][0]
    gain = float(np.clip(gain, *PARAMETER_BOUNDS["matrix_gain"]))
    return gain, raw_by_batch


def _fit_residual(
    log_shape: np.ndarray,
    batch_names: list[str],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    model_variant: str = MODEL_NAME,
) -> np.ndarray:
    gain, raw_by_batch = _profile_matrix_gain(
        log_shape, batch_names, observations, cache, model_variant
    )
    residuals = []
    for batch_name in batch_names:
        group = observations[observations["batch"].eq(batch_name)]
        observed = group["co2_rate_g_l_h"].to_numpy(dtype=float)
        sigma = max(0.04, 0.10 * float(group["co2_rate_g_l_h"].max()))
        prediction = gain * raw_by_batch[batch_name]
        censored = group["left_censored"].astype(bool).to_numpy()
        batch_residuals = []
        if (~censored).any():
            batch_residuals.append((prediction[~censored] - observed[~censored]) / sigma)
        if censored.any():
            limit = group["detection_limit_g_l_h"].to_numpy(dtype=float)[censored]
            batch_residuals.append(np.maximum(prediction[censored] - limit, 0.0) / sigma)
        residuals.append(np.concatenate(batch_residuals) / math.sqrt(len(group)))
        if _uses_chemistry_alignment(model_variant):
            threshold = _onset_threshold(group)
            observed_onset = _sustained_onset_h(
                group["t_h"].to_numpy(dtype=float), observed, threshold
            )
            predicted_onset = _sustained_onset_h(
                group["t_h"].to_numpy(dtype=float), prediction, threshold
            )
            if np.isfinite(observed_onset):
                if not np.isfinite(predicted_onset):
                    predicted_onset = float(group["t_h"].max()) + ONSET_SIGMA_H
                residuals.append(
                    np.asarray([(predicted_onset - observed_onset) / ONSET_SIGMA_H], dtype=float)
                )
        if _uses_nitrogen_boost_transition(model_variant):
            pulse_time_h = cache[batch_name].n_pulse_time_h
            observed_peak_h = _postpulse_peak_time_h(
                group["t_h"].to_numpy(dtype=float), observed, pulse_time_h
            )
            predicted_peak_h = _postpulse_peak_time_h(
                group["t_h"].to_numpy(dtype=float), prediction, pulse_time_h
            )
            if np.isfinite(observed_peak_h):
                if not np.isfinite(predicted_peak_h):
                    predicted_peak_h = min(
                        float(group["t_h"].max()), pulse_time_h + PULSE_PEAK_WINDOW_H
                    )
                residuals.append(
                    np.asarray(
                        [(predicted_peak_h - observed_peak_h) / PULSE_PEAK_SIGMA_H],
                        dtype=float,
                    )
                )
    return np.concatenate(residuals)


def fit_matrix(
    matrix: str,
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    n_starts: int = 5,
    max_nfev: int = 300,
    seed: int = 20260812,
    model_variant: str = MODEL_NAME,
    initial_parameters: dict[str, float] | None = None,
) -> tuple[dict[str, float], pd.DataFrame]:
    batch_names = sorted(
        observations.loc[
            observations["matrix"].eq(matrix)
            & observations["calibratable"].astype(bool)
            & ~observations["batch"].eq(HOLDOUTS[matrix]),
            "batch",
        ].unique()
    )
    if not batch_names:
        raise RuntimeError(f"No calibration batches are available for {matrix}")

    names = _shape_parameter_names(model_variant)
    lower = np.log([PARAMETER_BOUNDS[name][0] for name in names])
    upper = np.log([PARAMETER_BOUNDS[name][1] for name in names])
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
    if initial_parameters is not None:
        reference_values.update(
            {
                name: float(initial_parameters[name])
                for name in names
                if name in initial_parameters
            }
        )
    reference = np.log([reference_values[name] for name in names])
    rng = np.random.default_rng(seed + (0 if matrix == "synthetic" else 1000))
    starts = [reference] + [rng.uniform(lower + 1e-6, upper - 1e-6) for _ in range(max(n_starts - 1, 0))]
    rows = []
    best = None
    best_wsse = np.inf
    for start_idx, start in enumerate(starts):
        result = least_squares(
            lambda value: _fit_residual(
                value, batch_names, observations, cache, model_variant
            ),
            np.clip(start, lower + 1e-9, upper - 1e-9),
            bounds=(lower, upper),
            method="trf",
            x_scale="jac",
            max_nfev=max_nfev,
        )
        residual = _fit_residual(
            result.x, batch_names, observations, cache, model_variant
        )
        wsse = float(np.dot(residual, residual))
        gain, _ = _profile_matrix_gain(
            result.x, batch_names, observations, cache, model_variant
        )
        values = _shape_values(result.x, model_variant)
        row = {
                "matrix": matrix,
                "model": model_variant,
                "start": start_idx,
                "success": bool(result.success),
                "status": int(result.status),
                "nfev": int(result.nfev),
                "wsse_equal_batch": wsse,
                "matrix_gain": gain,
                "message": str(result.message),
            }
        row.update(values)
        rows.append(row)
        if wsse < best_wsse:
            best_wsse = wsse
            best = result

    if best is None:
        raise RuntimeError(f"Optimization did not return a result for {matrix}")
    values = _shape_values(best.x, model_variant)
    gain, _ = _profile_matrix_gain(
        best.x, batch_names, observations, cache, model_variant
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


def _calibration_batch_names(matrix: str, observations: pd.DataFrame) -> list[str]:
    return sorted(
        observations.loc[
            observations["matrix"].eq(matrix)
            & observations["calibratable"].astype(bool)
            & ~observations["batch"].eq(HOLDOUTS[matrix]),
            "batch",
        ].unique()
    )


def jacobian_identifiability_diagnostics(
    fits: dict[str, dict[str, float]],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Local practical-identifiability audit on log-parameter sensitivities."""

    summary_rows: list[dict[str, object]] = []
    singular_rows: list[dict[str, object]] = []
    correlation_rows: list[dict[str, object]] = []
    for matrix, fitted in fits.items():
        names = list(_shape_parameter_names(MODEL_NAME))
        batches = _calibration_batch_names(matrix, observations)
        center = np.log([float(fitted[name]) for name in names])
        base_residual = _fit_residual(center, batches, observations, cache, MODEL_NAME)
        jacobian = np.empty((len(base_residual), len(names)), dtype=float)
        step = 1e-4
        for column in range(len(names)):
            plus = center.copy()
            minus = center.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (
                _fit_residual(plus, batches, observations, cache, MODEL_NAME)
                - _fit_residual(minus, batches, observations, cache, MODEL_NAME)
            ) / (2.0 * step)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        tolerance = max(jacobian.shape) * np.finfo(float).eps * max(float(singular[0]), 1e-15)
        rank = int(np.sum(singular > tolerance))
        condition = float(singular[0] / singular[-1]) if singular[-1] > tolerance else np.inf
        dof = max(len(base_residual) - len(names), 1)
        residual_variance = float(np.dot(base_residual, base_residual) / dof)
        covariance = residual_variance * np.linalg.pinv(jacobian.T @ jacobian, rcond=1e-12)
        standard_error = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        denominator = np.outer(standard_error, standard_error)
        correlation = np.divide(
            covariance,
            denominator,
            out=np.full_like(covariance, np.nan),
            where=denominator > 0,
        )
        summary_rows.append(
            {
                "matrix": matrix,
                "n_residuals": int(len(base_residual)),
                "n_shape_parameters": int(len(names)),
                "effective_rank": rank,
                "condition_number_log_parameter_jacobian": condition,
                "smallest_singular_value": float(singular[-1]),
                "largest_singular_value": float(singular[0]),
                "residual_variance_scale": residual_variance,
                "interpretation": "ill_conditioned" if condition > 1e4 else "moderate_or_better",
            }
        )
        for index, value in enumerate(singular, start=1):
            singular_rows.append(
                {"matrix": matrix, "singular_index": index, "singular_value": float(value)}
            )
        for i, first in enumerate(names):
            for j, second in enumerate(names):
                correlation_rows.append(
                    {
                        "matrix": matrix,
                        "parameter_1": first,
                        "parameter_2": second,
                        "local_log_parameter_correlation": float(correlation[i, j]),
                        "log_parameter_standard_error_1": float(standard_error[i]),
                    }
                )
    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(singular_rows),
        pd.DataFrame(correlation_rows),
    )


def profile_activation_parameters(
    fits: dict[str, dict[str, float]],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
) -> pd.DataFrame:
    """One-dimensional objective profiles with all other parameters reoptimized."""

    rows: list[dict[str, object]] = []
    profiled_parameters = (
        "chem_activation_start_fraction",
        "chem_activation_duration_fraction",
    )
    for matrix, fitted in fits.items():
        names = list(_shape_parameter_names(MODEL_NAME))
        batches = _calibration_batch_names(matrix, observations)
        optimum = np.log([float(fitted[name]) for name in names])
        for profiled_name in profiled_parameters:
            fixed_index = names.index(profiled_name)
            lower_value, upper_value = PARAMETER_BOUNDS[profiled_name]
            if profiled_name.endswith("duration_fraction"):
                grid = np.geomspace(lower_value, upper_value, IDENTIFIABILITY_PROFILE_POINTS)
            else:
                grid = np.linspace(lower_value, upper_value, IDENTIFIABILITY_PROFILE_POINTS)
            grid = np.unique(np.append(grid, float(fitted[profiled_name])))
            free_indices = [index for index in range(len(names)) if index != fixed_index]
            free_lower = np.log([PARAMETER_BOUNDS[names[index]][0] for index in free_indices])
            free_upper = np.log([PARAMETER_BOUNDS[names[index]][1] for index in free_indices])
            for fixed_value in grid:
                fixed_log = math.log(float(fixed_value))

                def residual_free(free_log: np.ndarray) -> np.ndarray:
                    full = optimum.copy()
                    full[fixed_index] = fixed_log
                    full[free_indices] = free_log
                    return _fit_residual(full, batches, observations, cache, MODEL_NAME)

                result = least_squares(
                    residual_free,
                    np.clip(optimum[free_indices], free_lower + 1e-9, free_upper - 1e-9),
                    bounds=(free_lower, free_upper),
                    method="trf",
                    x_scale="jac",
                    max_nfev=IDENTIFIABILITY_PROFILE_MAX_NFEV,
                )
                residual = residual_free(result.x)
                rows.append(
                    {
                        "matrix": matrix,
                        "profiled_parameter": profiled_name,
                        "fixed_value": float(fixed_value),
                        "wsse_equal_batch": float(np.dot(residual, residual)),
                        "success": bool(result.success),
                        "nfev": int(result.nfev),
                        "full_fit_estimate": float(fitted[profiled_name]),
                    }
                )
    profiles = pd.DataFrame(rows)
    profiles["delta_wsse"] = profiles["wsse_equal_batch"] - profiles.groupby(
        ["matrix", "profiled_parameter"]
    )["wsse_equal_batch"].transform("min")
    profiles["relative_wsse_increase"] = profiles["delta_wsse"] / profiles.groupby(
        ["matrix", "profiled_parameter"]
    )["wsse_equal_batch"].transform("min").clip(lower=1e-12)
    return profiles.sort_values(["matrix", "profiled_parameter", "fixed_value"]).reset_index(drop=True)


def leave_one_batch_out_stability(
    fits: dict[str, dict[str, float]],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refit after omitting each calibration fermentation in turn."""

    rows: list[dict[str, object]] = []
    for matrix, full_fit in fits.items():
        for omitted_batch in _calibration_batch_names(matrix, observations):
            reduced = observations.copy()
            reduced.loc[reduced["batch"].eq(omitted_batch), "calibratable"] = False
            parameters, starts = fit_matrix(
                matrix,
                reduced,
                cache,
                n_starts=1,
                max_nfev=LOO_MAX_NFEV,
                seed=seed,
                model_variant=MODEL_NAME,
                initial_parameters=full_fit,
            )
            best_start = starts.sort_values("wsse_equal_batch").iloc[0]
            row = {
                "matrix": matrix,
                "omitted_batch": omitted_batch,
                "success": bool(best_start["success"]),
                "nfev": int(best_start["nfev"]),
                "wsse_equal_batch": float(parameters["wsse_equal_batch"]),
            }
            for name in _shape_parameter_names(MODEL_NAME):
                row[name] = float(parameters[name])
                row[f"full_fit_{name}"] = float(full_fit[name])
            rows.append(row)
    estimates = pd.DataFrame(rows)
    summaries: list[dict[str, object]] = []
    for matrix, group in estimates.groupby("matrix"):
        for parameter in (
            "chem_activation_start_fraction",
            "chem_activation_duration_fraction",
        ):
            values = group[parameter].to_numpy(dtype=float)
            full_value = float(group[f"full_fit_{parameter}"].iloc[0])
            summaries.append(
                {
                    "matrix": matrix,
                    "parameter": parameter,
                    "full_fit_estimate": full_value,
                    "loo_min": float(np.min(values)),
                    "loo_median": float(np.median(values)),
                    "loo_max": float(np.max(values)),
                    "loo_range": float(np.ptp(values)),
                    "loo_relative_range_to_full": float(np.ptp(values) / max(abs(full_value), 1e-12)),
                    "n_omissions": int(len(values)),
                    "n_at_lower_bound": int(
                        np.isclose(values, PARAMETER_BOUNDS[parameter][0], rtol=0.0, atol=1e-4).sum()
                    ),
                    "n_at_upper_bound": int(
                        np.isclose(values, PARAMETER_BOUNDS[parameter][1], rtol=0.0, atol=1e-4).sum()
                    ),
                }
            )
    return estimates, pd.DataFrame(summaries)


def _metric_row(
    calibration_matrix: str,
    target_matrix: str,
    batch: str,
    observed: np.ndarray,
    predicted: np.ndarray,
    time_h: np.ndarray,
    role: str,
    left_censored: np.ndarray,
    detection_limit: np.ndarray,
    model_variant: str = MODEL_NAME,
    pulse_time_h: float = np.nan,
) -> dict[str, object]:
    left_censored = np.asarray(left_censored, dtype=bool)
    quantifiable = ~left_censored
    observed_quantifiable = observed[quantifiable]
    predicted_quantifiable = predicted[quantifiable]
    time_quantifiable = time_h[quantifiable]
    error = predicted_quantifiable - observed_quantifiable
    rmse = float(np.sqrt(np.mean(error * error)))
    observed_peak = float(np.max(observed_quantifiable))
    corr = np.nan
    if np.std(observed_quantifiable) > 1e-12 and np.std(predicted_quantifiable) > 1e-12:
        corr = float(np.corrcoef(observed_quantifiable, predicted_quantifiable)[0, 1])
    ss_tot = float(np.sum((observed_quantifiable - np.mean(observed_quantifiable)) ** 2))
    r2 = float(1.0 - np.sum(error * error) / ss_tot) if ss_tot > 1e-16 else np.nan
    observed_integral = float(np.trapz(observed, time_h))
    predicted_integral = float(np.trapz(predicted, time_h))
    onset_frame = pd.DataFrame(
        {
            "t_h": time_h,
            "co2_rate_g_l_h": observed,
            "detection_limit_g_l_h": detection_limit,
        }
    )
    onset_threshold = _onset_threshold(onset_frame)
    observed_onset = _sustained_onset_h(time_h, observed, onset_threshold)
    predicted_onset = _sustained_onset_h(time_h, predicted, onset_threshold)
    early_mask = time_h <= ONSET_BASELINE_WINDOW_H
    observed_baseline = float(np.median(observed[early_mask])) if early_mask.any() else float(observed[0])
    dynamic_range = max(float(np.max(observed)) - observed_baseline, 0.0)
    visual_low_threshold = max(
        EARLY_EMISSION_THRESHOLD_G_L_H,
        observed_baseline + INITIAL_RISE_LOW_FRACTION * dynamic_range,
    )
    visual_high_threshold = max(
        visual_low_threshold,
        observed_baseline + INITIAL_RISE_HIGH_FRACTION * dynamic_range,
    )
    observed_rise_low_h, observed_rise_high_h, observed_rise_duration_h = _initial_rise_times_h(
        time_h, observed, visual_low_threshold, visual_high_threshold
    )
    predicted_rise_low_h, predicted_rise_high_h, predicted_rise_duration_h = _initial_rise_times_h(
        time_h, predicted, visual_low_threshold, visual_high_threshold
    )
    predicted_first_emission_h = _sustained_onset_h(
        time_h, predicted, EARLY_EMISSION_THRESHOLD_G_L_H, consecutive=2
    )
    observed_postpulse_peak = _postpulse_peak_time_h(
        time_h, observed, pulse_time_h
    )
    predicted_postpulse_peak = _postpulse_peak_time_h(
        time_h, predicted, pulse_time_h
    )
    return {
        "model": model_variant,
        "calibration_matrix": calibration_matrix,
        "target_matrix": target_matrix,
        "batch": batch,
        "experiment_code": SYNTHETIC_CODES.get(batch, batch),
        "role": role,
        "n": int(quantifiable.sum()),
        "n_total": int(len(observed)),
        "n_left_censored": int(left_censored.sum()),
        "rmse_g_l_h": rmse,
        "mae_g_l_h": float(np.mean(np.abs(error))),
        "bias_g_l_h": float(np.mean(error)),
        "nrmse_peak": rmse / max(observed_peak, 1e-12),
        "correlation": corr,
        "r2": r2,
        "observed_peak_g_l_h": observed_peak,
        "predicted_peak_g_l_h": float(np.max(predicted_quantifiable)),
        "observed_peak_time_h": float(time_quantifiable[int(np.argmax(observed_quantifiable))]),
        "predicted_peak_time_h": float(time_quantifiable[int(np.argmax(predicted_quantifiable))]),
        "integral_ratio_pred_over_obs": predicted_integral / max(observed_integral, 1e-12),
        "integral_ratio_pred_over_observed_lower_bound": predicted_integral / max(observed_integral, 1e-12),
        "onset_threshold_g_l_h": onset_threshold,
        "observed_onset_h": observed_onset,
        "predicted_onset_h": predicted_onset,
        "onset_delay_h": predicted_onset - observed_onset
        if np.isfinite(observed_onset) and np.isfinite(predicted_onset)
        else np.nan,
        "visual_rise_low_threshold_g_l_h": visual_low_threshold,
        "visual_rise_high_threshold_g_l_h": visual_high_threshold,
        "observed_visual_rise_low_h": observed_rise_low_h,
        "observed_visual_rise_high_h": observed_rise_high_h,
        "observed_visual_rise_duration_h": observed_rise_duration_h,
        "predicted_first_emission_0p005_h": predicted_first_emission_h,
        "predicted_visual_rise_low_h": predicted_rise_low_h,
        "predicted_visual_rise_high_h": predicted_rise_high_h,
        "predicted_visual_rise_duration_h": predicted_rise_duration_h,
        "visual_rise_duration_error_h": predicted_rise_duration_h - observed_rise_duration_h
        if np.isfinite(observed_rise_duration_h) and np.isfinite(predicted_rise_duration_h)
        else np.nan,
        "pulse_time_h": pulse_time_h,
        "observed_postpulse_peak_h": observed_postpulse_peak,
        "predicted_postpulse_peak_h": predicted_postpulse_peak,
        "observed_pulse_t_rise_h": observed_postpulse_peak - pulse_time_h
        if np.isfinite(observed_postpulse_peak) and np.isfinite(pulse_time_h)
        else np.nan,
        "predicted_pulse_t_rise_h": predicted_postpulse_peak - pulse_time_h
        if np.isfinite(predicted_postpulse_peak) and np.isfinite(pulse_time_h)
        else np.nan,
        "postpulse_peak_delay_h": predicted_postpulse_peak - observed_postpulse_peak
        if np.isfinite(observed_postpulse_peak) and np.isfinite(predicted_postpulse_peak)
        else np.nan,
    }


def predict_and_score(
    fits: dict[str, dict[str, float]],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    model_variant: str = MODEL_NAME,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prediction_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    usable = observations[observations["calibratable"].astype(bool)].copy()
    for calibration_matrix, parameters in fits.items():
        for batch_name, group in usable.groupby("batch", sort=True):
            target_matrix = str(group["matrix"].iloc[0])
            if target_matrix == calibration_matrix:
                role = "internal_holdout" if batch_name == HOLDOUTS[target_matrix] else "calibration"
            else:
                role = "cross_holdout" if batch_name == HOLDOUTS[target_matrix] else "cross_context"
            raw = raw_qgas_prediction(
                cache[batch_name],
                parameters["kCO2_release_h"],
                parameters["CO2sat_scale"],
                parameters["O2_qmax_mg_gdw_h"],
                parameters["O2_initial_scale"],
                chemistry_aligned=_uses_chemistry_alignment(model_variant),
                nitrogen_boost_transition=_uses_nitrogen_boost_transition(model_variant),
                pulse_t_rise_h=parameters["pulse_t_rise_h"],
                pulse_activity_gain=parameters["pulse_activity_gain"],
                continuous_release=_uses_continuous_release(model_variant),
                bounded_chemical_activation=_uses_bounded_chemical_activation(model_variant),
                chem_activation_start_fraction=parameters["chem_activation_start_fraction"],
                chem_activation_duration_fraction=parameters["chem_activation_duration_fraction"],
            )
            predicted = parameters["matrix_gain"] * raw
            observed = group["co2_rate_g_l_h"].to_numpy(dtype=float)
            time_h = group["t_h"].to_numpy(dtype=float)
            left_censored = group["left_censored"].astype(bool).to_numpy()
            detection_limit = group["detection_limit_g_l_h"].to_numpy(dtype=float)
            metric_rows.append(
                _metric_row(
                    calibration_matrix,
                    target_matrix,
                    batch_name,
                    observed,
                    predicted,
                    time_h,
                    role,
                    left_censored,
                    detection_limit,
                    model_variant,
                    cache[batch_name].n_pulse_time_h,
                )
            )
            for time_value, observed_value, predicted_value, raw_value, censored, artifact_fraction, detection_limit in zip(
                time_h,
                observed,
                predicted,
                raw,
                left_censored,
                group["artifact_fraction"].to_numpy(dtype=float),
                group["detection_limit_g_l_h"].to_numpy(dtype=float),
            ):
                prediction_rows.append(
                    {
                        "calibration_matrix": calibration_matrix,
                        "model": model_variant,
                        "target_matrix": target_matrix,
                        "batch": batch_name,
                        "experiment_code": SYNTHETIC_CODES.get(batch_name, batch_name),
                        "role": role,
                        "time_h": float(time_value),
                        "observed_g_l_h": float(observed_value),
                        "predicted_g_l_h": float(predicted_value),
                        "raw_model_qgas_g_l_h": float(raw_value),
                        "left_censored": bool(censored),
                        "detection_limit_g_l_h": float(detection_limit),
                        "artifact_fraction": float(artifact_fraction),
                        "pulse_time_h": cache[batch_name].n_pulse_time_h,
                        "pulse_amount_N_kg_m3": cache[batch_name].n_pulse_amount_kg_m3,
                    }
                )
    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    validation = metrics[metrics["role"].isin(["internal_holdout", "cross_holdout"])].copy()
    validation["scenario"] = np.where(
        validation["role"].eq("internal_holdout"), "within_matrix_holdout", "cross_matrix_holdout"
    )
    validation = validation.sort_values(["target_matrix", "scenario", "calibration_matrix"]).reset_index(drop=True)
    return predictions, metrics, validation


def fit_parameter_table(fits: dict[str, dict[str, float]]) -> pd.DataFrame:
    rows = []
    for matrix, parameters in fits.items():
        model_variant = str(parameters["model"])
        for parameter in (*_shape_parameter_names(model_variant), "matrix_gain"):
            lower, upper = PARAMETER_BOUNDS[parameter]
            value = float(parameters[parameter])
            relative_bound_distance = min(
                math.log(value / lower), math.log(upper / value)
            ) / math.log(upper / lower)
            active = bool(relative_bound_distance <= 0.01)
            rows.append(
                {
                    "calibration_matrix": matrix,
                    "parameter": parameter,
                    "estimate": value,
                    "lower_bound": lower,
                    "upper_bound": upper,
                    "active_bound": active,
                    "relative_bound_distance": relative_bound_distance,
                }
            )
    return pd.DataFrame(rows)


def plot_data_overview(
    observations: pd.DataFrame,
    nutrient_pulses: pd.DataFrame | None = None,
    save: bool = True,
) -> plt.Figure:
    batches = list(observations.sort_values(["matrix", "batch"])["batch"].unique())
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 2.45 * nrows), sharey=False)
    axes = np.asarray(axes).reshape(-1)
    for axis, batch in zip(axes, batches):
        group = observations[observations["batch"].eq(batch)]
        color = "#0072B2" if group["matrix"].iloc[0] == "synthetic" else "#009E73"
        axis.plot(
            group["t_h"],
            group["co2_rate_raw_g_l_h"],
            color="#9e9e9e",
            linewidth=0.8,
            alpha=0.8,
            label="raw hourly",
        )
        axis.plot(
            group["t_h"], group["co2_rate_g_l_h"],
            color=color, linewidth=1.6, label="robust smoothed",
        )
        artifact = group["artifact_fraction"].gt(0.0)
        if artifact.any():
            axis.scatter(
                group.loc[artifact, "t_h"],
                group.loc[artifact, "co2_rate_raw_g_l_h"],
                marker="x",
                color="#D55E00",
                s=20,
                linewidth=1.0,
                label="transient replaced",
            )
        censored = group["left_censored"].astype(bool)
        if censored.any():
            axis.scatter(
                group.loc[censored, "t_h"],
                group.loc[censored, "detection_limit_g_l_h"],
                facecolors="none",
                edgecolors="#CC79A7",
                s=16,
                linewidth=0.8,
                label="below detection",
            )
        if nutrient_pulses is not None:
            pulse_times = nutrient_pulses.loc[
                nutrient_pulses["batch"].eq(batch), "pulse_time_h"
            ]
            for pulse_time in pulse_times:
                axis.axvline(
                    float(pulse_time), color="#CC79A7", linestyle="--", linewidth=1.0,
                    alpha=0.85, label="N pulse" if pulse_time == pulse_times.iloc[0] else None,
                )
        if not bool(group["calibratable"].iloc[0]):
            axis.set_facecolor("#f3f3f3")
            axis.text(0.98, 0.93, "CO2-only", transform=axis.transAxes, ha="right", va="top", fontsize=8)
        if batch in HOLDOUTS.values():
            axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch}) — HOLDOUT", fontsize=9)
        else:
            axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch})", fontsize=9)
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
    for axis in axes[len(batches) :]:
        axis.axis("off")
    if len(batches):
        axes[0].legend(loc="upper left", frameon=False, fontsize=7)
    fig.suptitle("2026 CO2 sensor QC: raw, robust-smoothed and below-detection intervals", y=1.005)
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "co2_data_overview.png", dpi=170, bbox_inches="tight")
    return fig


def plot_sensor_filter_examples(
    sensor_qc: pd.DataFrame,
    sampling_schedule: pd.DataFrame,
    nutrient_pulses: pd.DataFrame | None = None,
    save: bool = True,
    n_examples: int = 4,
) -> plt.Figure:
    """Zoom the strongest audited sampling disturbances at the native grain."""

    candidates = []
    audited = sensor_qc[
        sensor_qc["artifact_reason"].eq("sampling_window_transient")
        & sensor_qc["artifact_sample_time_h"].notna()
    ].copy()
    for (batch, sample_time), window in audited.groupby(
        ["batch", "artifact_sample_time_h"], sort=False
    ):
        score = float(
            (window["co2_rate_raw_g_l_h"] - window["co2_rate_filtered_g_l_h"])
            .abs()
            .max()
        )
        candidates.append((score, str(batch), float(sample_time)))
    candidates = sorted(candidates, reverse=True)[:n_examples]
    if not candidates:
        raise ValueError("No sampling-window artifacts available for the zoom figure")

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), squeeze=False)
    axes = axes.reshape(-1)
    for axis, (_, batch, sample_time) in zip(axes, candidates):
        group = sensor_qc[
            sensor_qc["batch"].eq(batch)
            & sensor_qc["t_h"].between(sample_time - 4.0, sample_time + 4.0)
        ].sort_values("t_h")
        color = "#0072B2" if group["matrix"].iloc[0] == "synthetic" else "#009E73"
        axis.plot(
            group["t_h"], group["co2_rate_raw_g_l_h"],
            color="#777777", linewidth=1.4, marker=".", markersize=3, label="raw 10 min",
        )
        axis.plot(
            group["t_h"], group["co2_rate_filtered_g_l_h"],
            color=color, linewidth=2.0, label="filtered",
        )
        artifact = group["artifact_sample_time_h"].sub(sample_time).abs().lt(1e-8)
        axis.scatter(
            group.loc[artifact, "t_h"], group.loc[artifact, "co2_rate_raw_g_l_h"],
            marker="x", s=42, color="#D55E00", zorder=4, label="replaced",
        )
        axis.axvspan(
            sample_time - SAMPLE_EVENT_HALF_WIDTH_H,
            sample_time + SAMPLE_EVENT_HALF_WIDTH_H,
            color="#E69F00", alpha=0.10, label="event window",
        )
        axis.axvline(sample_time, color="#D55E00", linestyle="--", linewidth=1.2, label="sample")
        if nutrient_pulses is not None:
            pulse_times = nutrient_pulses.loc[
                nutrient_pulses["batch"].eq(batch), "pulse_time_h"
            ]
            for pulse_time in pulse_times:
                if float(group["t_h"].min()) <= float(pulse_time) <= float(group["t_h"].max()):
                    axis.axvline(
                        float(pulse_time), color="#CC79A7", linestyle="-.", linewidth=1.2,
                        label="N pulse",
                    )
        axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch}) · sample {sample_time:.1f} h")
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
    for axis in axes[len(candidates):]:
        axis.set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.005),
        ncol=5, frameon=False,
    )
    fig.suptitle("Audited CO2 sampling disturbances: native signal and reconstruction", y=0.985)
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.94))
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "sensor_artifact_filter_examples.png", dpi=180, bbox_inches="tight")
    return fig


def plot_temperature_profiles(
    observations: pd.DataFrame,
    temperature_inputs: pd.DataFrame,
    nutrient_pulses: pd.DataFrame | None = None,
    save: bool = True,
) -> plt.Figure:
    """Plot measured temperature, controller setpoint and model input."""

    batches = list(temperature_inputs.sort_values(["matrix", "batch"])["batch"].unique())
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 2.55 * nrows), sharey=True)
    axes = np.asarray(axes).reshape(-1)
    for axis, batch in zip(axes, batches):
        model = temperature_inputs[temperature_inputs["batch"].eq(batch)]
        sensor = observations[observations["batch"].eq(batch)]
        axis.plot(
            model["time_h"], model["temperature_used_c"],
            color="#0072B2", linewidth=2.0, label="temperature used by model",
        )
        axis.plot(
            sensor["t_h"], sensor["sensor_temperature_c"],
            color="#333333", linewidth=1.0, linestyle="--", alpha=0.85,
            label="measured temperature",
        )
        if sensor["setpoint_c"].notna().any():
            axis.step(
                sensor["t_h"], sensor["setpoint_c"], where="post",
                color="#E69F00", linewidth=1.2, alpha=0.9, label="setpoint",
            )
        if nutrient_pulses is not None:
            pulse_times = nutrient_pulses.loc[
                nutrient_pulses["batch"].eq(batch), "pulse_time_h"
            ]
            for pulse_time in pulse_times:
                axis.axvline(
                    float(pulse_time), color="#CC79A7", linestyle="--", linewidth=1.0,
                    alpha=0.85, label="N pulse" if pulse_time == pulse_times.iloc[0] else None,
                )
        suffix = " — HOLDOUT" if batch in HOLDOUTS.values() else ""
        axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch}){suffix}", fontsize=9)
        axis.set_xlabel("time [h]")
        axis.set_ylabel("temperature [°C]")
        axis.grid(alpha=0.2)
    for axis in axes[len(batches):]:
        axis.set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3, frameon=False)
    fig.suptitle("Temperature profiles used by the CO2 model", y=0.995)
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.975))
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "temperature_profiles_model_input.png", dpi=180, bbox_inches="tight")
    return fig


def plot_process_timeline_alignment(
    inventory: pd.DataFrame,
    nutrient_pulses: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Show analytical horizon, available/used CO2 support and nutrient pulses."""

    frame = inventory.sort_values(["matrix", "batch"]).reset_index(drop=True)
    fig, axis = plt.subplots(figsize=(12, max(5.0, 0.48 * len(frame) + 1.5)))
    for position, row in enumerate(frame.itertuples(index=False)):
        axis.plot(
            [row.chemistry_first_h, row.chemistry_last_h], [position, position],
            color="#4D4D4D", linewidth=5.0, alpha=0.32,
            label="chemical process interval" if position == 0 else None,
        )
        axis.plot(
            [row.co2_source_first_h, row.co2_source_last_h],
            [position + 0.13, position + 0.13],
            color="#E69F00", linewidth=2.2,
            label="CO2 source coverage" if position == 0 else None,
        )
        axis.plot(
            [row.co2_first_h, row.co2_last_h], [position - 0.13, position - 0.13],
            color="#0072B2", linewidth=2.8,
            label="CO2 used in fit/validation" if position == 0 else None,
        )
        for pulse_time in nutrient_pulses.loc[
            nutrient_pulses["batch"].eq(row.batch), "pulse_time_h"
        ]:
            axis.scatter(
                float(pulse_time), position, marker="D", s=34, color="#CC79A7", zorder=4,
                label="N pulse" if position == 0 else None,
            )
    axis.set_yticks(np.arange(len(frame)))
    axis.set_yticklabels(
        [f"{SYNTHETIC_CODES.get(row.batch, row.batch)} ({row.batch})" for row in frame.itertuples()]
    )
    axis.set_xlabel("time from first chemical sample [h]")
    axis.set_title("Temporal alignment: chemistry defines the process interval")
    axis.grid(axis="x", alpha=0.2)
    handles, labels = axis.get_legend_handles_labels()
    deduplicated = dict(zip(labels, handles))
    axis.legend(deduplicated.values(), deduplicated.keys(), loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=4, frameon=False)
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "process_timeline_alignment.png", dpi=180, bbox_inches="tight")
    return fig


def plot_parameter_comparison(parameter_table: pd.DataFrame, save: bool = True) -> plt.Figure:
    names = parameter_table["parameter"].drop_duplicates().tolist()
    ncols = 3
    nrows = int(math.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.1 * nrows))
    axes = axes.reshape(-1)
    colors = {"synthetic": "#0072B2", "natural": "#009E73"}
    for axis, parameter in zip(axes, names):
        group = parameter_table[parameter_table["parameter"].eq(parameter)]
        axis.bar(group["calibration_matrix"], group["estimate"], color=[colors[value] for value in group["calibration_matrix"]])
        lower, _ = PARAMETER_BOUNDS[parameter]
        axis.axhline(lower, color="black", linestyle=":", linewidth=0.8)
        axis.set_ylim(0.0, max(float(group["estimate"].max()) * 1.25, lower * 1.5))
        axis.set_title(parameter)
        axis.grid(axis="y", alpha=0.2)
    for axis in axes[len(names) :]:
        axis.set_visible(False)
    fig.suptitle("Medium-specific CO2-layer and nitrogen-boost calibration")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "parameter_comparison.png", dpi=170, bbox_inches="tight")
    return fig


def plot_calibration_overlays(
    predictions: pd.DataFrame,
    matrix: str,
    save: bool = True,
) -> plt.Figure:
    selected = predictions[
        predictions["calibration_matrix"].eq(matrix)
        & predictions["target_matrix"].eq(matrix)
        & predictions["role"].eq("calibration")
    ].copy()
    batches = sorted(selected["batch"].unique())
    if not batches:
        raise RuntimeError(f"No calibration prediction rows are available for {matrix}")
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 2.8 * nrows))
    axes = np.asarray(axes).reshape(-1)
    color = "#0072B2" if matrix == "synthetic" else "#009E73"
    for axis, batch in zip(axes, batches):
        group = selected[selected["batch"].eq(batch)]
        observed = group["observed_g_l_h"].to_numpy(dtype=float)
        predicted = group["predicted_g_l_h"].to_numpy(dtype=float)
        quantifiable = ~group["left_censored"].astype(bool).to_numpy()
        rmse = float(np.sqrt(np.mean((predicted[quantifiable] - observed[quantifiable]) ** 2)))
        axis.plot(group["time_h"], observed, color="black", linewidth=1.25, label="observed")
        axis.plot(group["time_h"], predicted, color=color, linewidth=1.4, label="calibrated model")
        pulse_time_h = float(group["pulse_time_h"].iloc[0])
        if np.isfinite(pulse_time_h):
            axis.axvline(
                pulse_time_h,
                color="#CC79A7",
                linestyle="--",
                linewidth=1.2,
                label="N pulse",
            )
        censored = group["left_censored"].astype(bool)
        if censored.any():
            axis.scatter(
                group.loc[censored, "time_h"],
                group.loc[censored, "detection_limit_g_l_h"],
                facecolors="none",
                edgecolors="#CC79A7",
                s=18,
                linewidth=0.8,
                label="left-censored",
            )
        axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch})", fontsize=9)
        axis.text(
            0.98,
            0.94,
            f"RMSE={rmse:.3f}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False, fontsize=8)
    for axis in axes[len(batches) :]:
        axis.axis("off")
    fig.suptitle(f"CO2 calibration overlays — {matrix} matrix ({len(batches)} complete fermentations)")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / f"calibration_overlays_{matrix}.png", dpi=180, bbox_inches="tight")
    return fig


def plot_validation(predictions: pd.DataFrame, validation: pd.DataFrame, save: bool = True) -> plt.Figure:
    scenarios = [
        ("synthetic", "synthetic", HOLDOUTS["synthetic"], "Within synthetic"),
        ("natural", "synthetic", HOLDOUTS["synthetic"], "Natural → synthetic"),
        ("natural", "natural", HOLDOUTS["natural"], "Within natural"),
        ("synthetic", "natural", HOLDOUTS["natural"], "Synthetic → natural"),
    ]
    colors = {"synthetic": "#0072B2", "natural": "#009E73"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for axis, (calibration_matrix, target_matrix, batch, title) in zip(axes.reshape(-1), scenarios):
        group = predictions[
            predictions["calibration_matrix"].eq(calibration_matrix) & predictions["batch"].eq(batch)
        ].copy()
        metric = validation[
            validation["calibration_matrix"].eq(calibration_matrix) & validation["batch"].eq(batch)
        ].iloc[0]
        axis.plot(group["time_h"], group["observed_g_l_h"], color="black", linewidth=1.4, label="observed")
        axis.plot(
            group["time_h"],
            group["predicted_g_l_h"],
            color=colors[calibration_matrix],
            linewidth=1.5,
            label=f"fit on {calibration_matrix}",
        )
        pulse_time_h = float(group["pulse_time_h"].iloc[0])
        if np.isfinite(pulse_time_h):
            axis.axvline(
                pulse_time_h,
                color="#CC79A7",
                linestyle="--",
                linewidth=1.2,
                label="N pulse",
            )
        censored = group["left_censored"].astype(bool)
        if censored.any():
            axis.scatter(
                group.loc[censored, "time_h"],
                group.loc[censored, "detection_limit_g_l_h"],
                facecolors="none",
                edgecolors="#CC79A7",
                s=20,
                linewidth=0.8,
                label="left-censored",
            )
        axis.set_title(f"{title}: {SYNTHETIC_CODES.get(batch, batch)}")
        axis.text(
            0.98,
            0.95,
            f"RMSE={metric.rmse_g_l_h:.3f}\nNRMSEpeak={metric.nrmse_peak:.2f}\nr={metric.correlation:.2f}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False)
    fig.suptitle("Strict held-out validation: no CO2 refit on DOE-F06 or LAB012")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "heldout_cross_matrix_validation.png", dpi=180, bbox_inches="tight")
    return fig


def transition_diagnostics(
    fits: dict[str, dict[str, float]],
    observations: pd.DataFrame,
    cache: dict[str, DriverCache],
    model_variant: str = MODEL_NAME,
) -> pd.DataFrame:
    """Audit the modeled source transition without using CO2 to locate chemistry activity."""

    rows: list[dict[str, object]] = []
    usable = observations[observations["calibratable"].astype(bool)].copy()
    for calibration_matrix, parameters in fits.items():
        for batch_name, group in usable.groupby("batch", sort=True):
            item = cache[batch_name]
            qgas_grid, qprod_grid, o2_grid, phi_ana_grid = raw_qgas_grid_prediction(
                item,
                parameters["kCO2_release_h"],
                parameters["CO2sat_scale"],
                parameters["O2_qmax_mg_gdw_h"],
                parameters["O2_initial_scale"],
                chemistry_aligned=_uses_chemistry_alignment(model_variant),
                nitrogen_boost_transition=_uses_nitrogen_boost_transition(model_variant),
                pulse_t_rise_h=parameters["pulse_t_rise_h"],
                pulse_activity_gain=parameters["pulse_activity_gain"],
                continuous_release=_uses_continuous_release(model_variant),
                bounded_chemical_activation=_uses_bounded_chemical_activation(model_variant),
                chem_activation_start_fraction=parameters["chem_activation_start_fraction"],
                chem_activation_duration_fraction=parameters["chem_activation_duration_fraction"],
            )
            predicted = parameters["matrix_gain"] * np.interp(
                group["t_h"].to_numpy(dtype=float), item.time_h, qgas_grid
            )
            threshold = _onset_threshold(group)
            observed_onset = _sustained_onset_h(
                group["t_h"].to_numpy(dtype=float),
                group["co2_rate_g_l_h"].to_numpy(dtype=float),
                threshold,
            )
            predicted_onset = _sustained_onset_h(
                group["t_h"].to_numpy(dtype=float), predicted, threshold
            )
            observed_postpulse_peak = _postpulse_peak_time_h(
                group["t_h"].to_numpy(dtype=float),
                group["co2_rate_g_l_h"].to_numpy(dtype=float),
                item.n_pulse_time_h,
            )
            predicted_postpulse_peak = _postpulse_peak_time_h(
                group["t_h"].to_numpy(dtype=float), predicted, item.n_pulse_time_h
            )
            o2_crossings = np.flatnonzero(o2_grid <= O2_ANA_K_MG_L)
            source_threshold = max(0.02, 0.10 * float(np.max(qprod_grid)))
            source_crossings = np.flatnonzero(qprod_grid >= source_threshold)
            rows.append(
                {
                    "model": model_variant,
                    "calibration_matrix": calibration_matrix,
                    "target_matrix": item.matrix,
                    "native_matrix_fit": bool(calibration_matrix == item.matrix),
                    "batch": batch_name,
                    "experiment_code": SYNTHETIC_CODES.get(batch_name, batch_name),
                    "median_setpoint_c": float(group["setpoint_c"].median())
                    if group["setpoint_c"].notna().any()
                    else np.nan,
                    "cold_operation": bool(
                        group["setpoint_c"].notna().any()
                        and float(group["setpoint_c"].median()) <= COLD_SETPOINT_THRESHOLD_C
                    ),
                    "chemical_activity_lower_h": item.chemical_activity_lower_h,
                    "chemical_activity_upper_h": item.chemical_activity_upper_h,
                    "chemical_activity_center_h": item.chemical_activity_center_h,
                    "source_10pct_peak_onset_h": float(item.time_h[source_crossings[0]])
                    if len(source_crossings)
                    else np.nan,
                    "o2_below_anaerobic_halfpoint_h": float(item.time_h[o2_crossings[0]])
                    if len(o2_crossings)
                    else np.nan,
                    "observed_onset_h": observed_onset,
                    "predicted_onset_h": predicted_onset,
                    "onset_delay_h": predicted_onset - observed_onset
                    if np.isfinite(observed_onset) and np.isfinite(predicted_onset)
                    else np.nan,
                    "onset_threshold_g_l_h": threshold,
                    "O2_qmax_mg_gdw_h": parameters["O2_qmax_mg_gdw_h"],
                    "O2_initial_scale": parameters["O2_initial_scale"],
                    "n_pulse_time_h": item.n_pulse_time_h,
                    "n_pulse_amount_kg_m3": item.n_pulse_amount_kg_m3,
                    "pulse_t_rise_h": parameters["pulse_t_rise_h"],
                    "pulse_activity_gain": parameters["pulse_activity_gain"],
                    "chem_activation_start_fraction": parameters["chem_activation_start_fraction"],
                    "chem_activation_duration_fraction": parameters["chem_activation_duration_fraction"],
                    "pulse_utilization_complete_h": item.n_pulse_time_h
                    + parameters["pulse_t_rise_h"]
                    if np.isfinite(item.n_pulse_time_h)
                    else np.nan,
                    "observed_postpulse_peak_h": observed_postpulse_peak,
                    "predicted_postpulse_peak_h": predicted_postpulse_peak,
                    "postpulse_peak_delay_h": predicted_postpulse_peak
                    - observed_postpulse_peak
                    if np.isfinite(observed_postpulse_peak)
                    and np.isfinite(predicted_postpulse_peak)
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_onset_model_comparison(
    batch_metrics: pd.DataFrame,
    previous_batch_metrics: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Compare native-fit onset error for threshold and continuous release."""

    keys = ["target_matrix", "batch", "experiment_code"]
    columns = keys + ["onset_delay_h"]
    new = batch_metrics[batch_metrics["calibration_matrix"].eq(batch_metrics["target_matrix"])][columns]
    old = previous_batch_metrics[
        previous_batch_metrics["calibration_matrix"].eq(previous_batch_metrics["target_matrix"])
    ][columns]
    frame = old.merge(new, on=keys, suffixes=("_previous", "_new"), validate="one_to_one")
    frame = frame.sort_values(["target_matrix", "batch"]).reset_index(drop=True)
    positions = np.arange(len(frame))
    fig, axis = plt.subplots(figsize=(10.5, max(5.2, 0.42 * len(frame) + 1.5)))
    for idx, row in frame.iterrows():
        axis.plot(
            [row["onset_delay_h_previous"], row["onset_delay_h_new"]],
            [idx, idx],
            color="#B0B0B0",
            linewidth=1.0,
            zorder=1,
        )
    axis.scatter(
        frame["onset_delay_h_previous"], positions, marker="x", s=42,
        color="#D55E00", label="threshold release (previous)", zorder=3,
    )
    axis.scatter(
        frame["onset_delay_h_new"], positions, marker="o", s=34,
        color="#0072B2", label="continuous liquid-gas release", zorder=3,
    )
    axis.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    axis.set_yticks(positions)
    axis.set_yticklabels(
        [f"{row.experiment_code} ({row.batch})" for row in frame.itertuples(index=False)]
    )
    axis.set_xlabel("predicted onset minus observed onset [h]")
    axis.set_title("CO2 onset error: threshold versus continuous release")
    axis.grid(axis="x", alpha=0.2)
    axis.legend(frameon=False, loc="best")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "onset_model_comparison.png", dpi=180, bbox_inches="tight")
    return fig


def plot_validation_model_comparison(
    predictions: pd.DataFrame,
    previous_predictions: pd.DataFrame,
    validation: pd.DataFrame,
    previous_validation: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Overlay both model variants on the four strict validation scenarios."""

    scenarios = [
        ("synthetic", "synthetic", HOLDOUTS["synthetic"], "Within synthetic"),
        ("natural", "synthetic", HOLDOUTS["synthetic"], "Natural to synthetic"),
        ("natural", "natural", HOLDOUTS["natural"], "Within natural"),
        ("synthetic", "natural", HOLDOUTS["natural"], "Synthetic to natural"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5))
    for axis, (calibration_matrix, _, batch, title) in zip(axes.reshape(-1), scenarios):
        new = predictions[
            predictions["calibration_matrix"].eq(calibration_matrix)
            & predictions["batch"].eq(batch)
        ].copy()
        old = previous_predictions[
            previous_predictions["calibration_matrix"].eq(calibration_matrix)
            & previous_predictions["batch"].eq(batch)
        ].copy()
        new_metric = validation[
            validation["calibration_matrix"].eq(calibration_matrix)
            & validation["batch"].eq(batch)
        ].iloc[0]
        old_metric = previous_validation[
            previous_validation["calibration_matrix"].eq(calibration_matrix)
            & previous_validation["batch"].eq(batch)
        ].iloc[0]
        axis.plot(new["time_h"], new["observed_g_l_h"], color="black", linewidth=1.4, label="observed")
        axis.plot(
            old["time_h"], old["predicted_g_l_h"], color="#D55E00", linewidth=1.2,
            linestyle=":", label="threshold release",
        )
        axis.plot(
            new["time_h"], new["predicted_g_l_h"], color="#0072B2", linewidth=1.6,
            label="continuous release",
        )
        pulse_time_h = float(new["pulse_time_h"].iloc[0])
        if np.isfinite(pulse_time_h):
            axis.axvline(
                pulse_time_h,
                color="#CC79A7",
                linestyle="--",
                linewidth=1.2,
                label="N pulse",
            )
        axis.set_title(f"{title}: {SYNTHETIC_CODES.get(batch, batch)}")
        axis.text(
            0.98, 0.95,
            f"onset error previous/new: {old_metric.onset_delay_h:.1f}/{new_metric.onset_delay_h:.1f} h\n"
            f"RMSE previous/new: {old_metric.rmse_g_l_h:.3f}/{new_metric.rmse_g_l_h:.3f}",
            transform=axis.transAxes, ha="right", va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False, fontsize=8)
    fig.suptitle("Strict holdouts: threshold and continuous CO2 release")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "heldout_model_comparison.png", dpi=180, bbox_inches="tight")
    return fig


def plot_pulse_response_comparison(
    predictions: pd.DataFrame,
    previous_predictions: pd.DataFrame,
    batch_metrics: pd.DataFrame,
    previous_batch_metrics: pd.DataFrame,
    source_transition_diagnostics: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Show the finite N-pulse response against the instantaneous prior variant."""

    native = predictions[
        predictions["calibration_matrix"].eq(predictions["target_matrix"])
        & predictions["pulse_time_h"].notna()
    ].copy()
    batches = sorted(native["batch"].unique())
    if not batches:
        raise RuntimeError("No pulsed batch is available for the pulse-response figure")
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.0 * nrows))
    axes = np.asarray(axes).reshape(-1)
    for axis, batch in zip(axes, batches):
        new = native[native["batch"].eq(batch)].sort_values("time_h")
        calibration_matrix = str(new["calibration_matrix"].iloc[0])
        old = previous_predictions[
            previous_predictions["calibration_matrix"].eq(calibration_matrix)
            & previous_predictions["batch"].eq(batch)
        ].sort_values("time_h")
        new_metric = batch_metrics[
            batch_metrics["calibration_matrix"].eq(calibration_matrix)
            & batch_metrics["batch"].eq(batch)
        ].iloc[0]
        old_metric = previous_batch_metrics[
            previous_batch_metrics["calibration_matrix"].eq(calibration_matrix)
            & previous_batch_metrics["batch"].eq(batch)
        ].iloc[0]
        diagnostic = source_transition_diagnostics[
            source_transition_diagnostics["calibration_matrix"].eq(calibration_matrix)
            & source_transition_diagnostics["batch"].eq(batch)
        ].iloc[0]
        pulse_time_h = float(new["pulse_time_h"].iloc[0])
        axis.plot(new["time_h"], new["observed_g_l_h"], color="black", linewidth=1.3, label="observed")
        axis.plot(
            old["time_h"], old["predicted_g_l_h"], color="#D55E00", linestyle=":",
            linewidth=1.2, label="threshold release",
        )
        axis.plot(
            new["time_h"], new["predicted_g_l_h"], color="#0072B2", linewidth=1.5,
            label="continuous release",
        )
        axis.axvline(
            pulse_time_h, color="#CC79A7", linestyle="--", linewidth=1.3, label="N pulse"
        )
        completion_h = float(diagnostic["pulse_utilization_complete_h"])
        if np.isfinite(completion_h):
            axis.axvline(
                completion_h, color="#7B3294", linestyle="-.", linewidth=1.0,
                label="fitted response complete",
            )
        axis.set_xlim(
            max(float(new["time_h"].min()), pulse_time_h - 24.0),
            min(float(new["time_h"].max()), pulse_time_h + PULSE_PEAK_WINDOW_H),
        )
        axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch})", fontsize=9)
        axis.text(
            0.98,
            0.95,
            f"peak delay previous/new: "
            f"{old_metric.postpulse_peak_delay_h:.1f}/{new_metric.postpulse_peak_delay_h:.1f} h\n"
            f"t_rise={diagnostic.pulse_t_rise_h:.1f} h; gain={diagnostic.pulse_activity_gain:.2f}",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False, fontsize=7)
    for axis in axes[len(batches) :]:
        axis.axis("off")
    fig.suptitle("CO2 response around the nutritional pulse (native matrix fits)")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "nutrient_pulse_response_comparison.png", dpi=180, bbox_inches="tight")
    return fig


def plot_initial_release_comparison(
    predictions: pd.DataFrame,
    previous_predictions: pd.DataFrame,
    batch_metrics: pd.DataFrame,
    previous_batch_metrics: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Zoom the initial CO2 rise for threshold and continuous-release variants."""

    native = predictions[
        predictions["calibration_matrix"].eq(predictions["target_matrix"])
    ].copy()
    batches = sorted(native["batch"].unique())
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 2.9 * nrows))
    axes = np.asarray(axes).reshape(-1)
    for axis, batch in zip(axes, batches):
        new = native[native["batch"].eq(batch)].sort_values("time_h")
        calibration_matrix = str(new["calibration_matrix"].iloc[0])
        old = previous_predictions[
            previous_predictions["calibration_matrix"].eq(calibration_matrix)
            & previous_predictions["batch"].eq(batch)
        ].sort_values("time_h")
        new_metric = batch_metrics[
            batch_metrics["calibration_matrix"].eq(calibration_matrix)
            & batch_metrics["batch"].eq(batch)
        ].iloc[0]
        old_metric = previous_batch_metrics[
            previous_batch_metrics["calibration_matrix"].eq(calibration_matrix)
            & previous_batch_metrics["batch"].eq(batch)
        ].iloc[0]
        observed_onset_h = float(new_metric["observed_onset_h"])
        end_h = min(
            float(new["time_h"].max()),
            observed_onset_h + 30.0 if np.isfinite(observed_onset_h) else 96.0,
        )
        start_h = float(new["time_h"].min())
        window = new["time_h"].between(start_h, end_h)
        axis.plot(
            new["time_h"], new["observed_g_l_h"], color="black", linewidth=1.25,
            label="observed",
        )
        axis.plot(
            old["time_h"], old["predicted_g_l_h"], color="#D55E00", linestyle=":",
            linewidth=1.25, label="threshold release",
        )
        axis.plot(
            new["time_h"], new["predicted_g_l_h"], color="#0072B2", linewidth=1.55,
            label="continuous release",
        )
        detection_limit = float(new.loc[window, "detection_limit_g_l_h"].median())
        axis.axhline(
            detection_limit, color="#CC79A7", linestyle="--", linewidth=0.9,
            label="operational detection limit",
        )
        pulse_time_h = float(new["pulse_time_h"].iloc[0])
        if np.isfinite(pulse_time_h) and start_h <= pulse_time_h <= end_h:
            axis.axvline(
                pulse_time_h, color="#CC79A7", linestyle="-.", linewidth=1.0,
                label="N pulse",
            )
        axis.set_xlim(start_h, end_h)
        y_values = np.concatenate(
            [
                new.loc[window, "observed_g_l_h"].to_numpy(dtype=float),
                new.loc[window, "predicted_g_l_h"].to_numpy(dtype=float),
                old.loc[old["time_h"].between(start_h, end_h), "predicted_g_l_h"].to_numpy(dtype=float),
            ]
        )
        axis.set_ylim(0.0, max(float(np.nanmax(y_values)) * 1.12, detection_limit * 1.5))
        axis.set_title(f"{SYNTHETIC_CODES.get(batch, batch)} ({batch})", fontsize=9)
        axis.text(
            0.98,
            0.95,
            f"first >0.005 old/new: "
            f"{old_metric.predicted_first_emission_0p005_h:.1f}/"
            f"{new_metric.predicted_first_emission_0p005_h:.1f} h\n"
            f"2-10% rise obs/old/new: {new_metric.observed_visual_rise_duration_h:.1f}/"
            f"{old_metric.predicted_visual_rise_duration_h:.1f}/"
            f"{new_metric.predicted_visual_rise_duration_h:.1f} h",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=7.2,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False, fontsize=6.8)
    for axis in axes[len(batches) :]:
        axis.axis("off")
    fig.suptitle("Initial CO2 rise: threshold versus continuous liquid-gas release")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "initial_gradual_release_comparison.png", dpi=180, bbox_inches="tight")
    return fig


def plot_sensor_zero_correction(
    sensor_qc: pd.DataFrame,
    sensor_zero_offsets: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Audit the run-specific sensor zero before artifact filtering."""

    batches = sorted(sensor_qc["batch"].unique())
    ncols = 3
    nrows = int(math.ceil(len(batches) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 2.75 * nrows), sharey=False)
    axes = np.asarray(axes).reshape(-1)
    for axis, batch in zip(axes, batches):
        group = sensor_qc[sensor_qc["batch"].eq(batch)].sort_values("t_h")
        first_h = float(group["t_h"].min())
        view = group[group["t_h"].le(first_h + 24.0)]
        audit = sensor_zero_offsets[sensor_zero_offsets["batch"].eq(batch)].iloc[0]
        axis.plot(
            view["t_h"],
            view["co2_rate_uncorrected_physical_g_l_h"],
            color="#D55E00",
            linewidth=1.0,
            alpha=0.85,
            label="before zero correction",
        )
        axis.plot(
            view["t_h"],
            view["co2_rate_raw_g_l_h"],
            color="#0072B2",
            linewidth=1.15,
            label="after zero correction",
        )
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.set_title(
            f"{SYNTHETIC_CODES.get(batch, batch)} | {audit.acquisition_channel} / sensor {int(audit.sensor_id)}\n"
            f"offset={audit.estimated_zero_offset_sccm:.3f} sccm",
            fontsize=8.5,
        )
        axis.set_xlabel("time [h]")
        axis.set_ylabel("CO2 [g L$^{-1}$ h$^{-1}$]")
        axis.grid(alpha=0.2)
        axis.legend(loc="upper left", frameon=False, fontsize=6.5)
    for axis in axes[len(batches) :]:
        axis.axis("off")
    fig.suptitle("Run-specific CO2 sensor zero correction — first 24 process hours", y=0.997)
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "sensor_zero_offset_correction.png", dpi=180, bbox_inches="tight")
    return fig


def plot_activation_identifiability(
    profiles: pd.DataFrame,
    loo_summary: pd.DataFrame,
    save: bool = True,
) -> plt.Figure:
    """Objective profiles with leave-one-batch-out ranges overlaid."""

    parameters = (
        "chem_activation_start_fraction",
        "chem_activation_duration_fraction",
    )
    labels = {
        "chem_activation_start_fraction": "activation start fraction",
        "chem_activation_duration_fraction": "activation duration fraction",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
    for row, matrix in enumerate(("synthetic", "natural")):
        for column, parameter in enumerate(parameters):
            axis = axes[row, column]
            group = profiles[
                profiles["matrix"].eq(matrix)
                & profiles["profiled_parameter"].eq(parameter)
            ].sort_values("fixed_value")
            summary = loo_summary[
                loo_summary["matrix"].eq(matrix)
                & loo_summary["parameter"].eq(parameter)
            ].iloc[0]
            axis.plot(
                group["fixed_value"], group["relative_wsse_increase"],
                color="#0072B2", marker="o", linewidth=1.4,
                label="reoptimized objective profile",
            )
            axis.axvspan(
                float(summary.loo_min), float(summary.loo_max),
                color="#E69F00", alpha=0.20, label="leave-one-batch-out range",
            )
            axis.axvline(
                float(summary.full_fit_estimate), color="black", linestyle="--",
                linewidth=1.0, label="full-data estimate",
            )
            axis.axhline(
                0.05, color="#CC79A7", linestyle=":", linewidth=0.9,
                label="5% objective increase (heuristic)",
            )
            axis.set_title(f"{matrix}: {labels[parameter]}")
            axis.set_xlabel(labels[parameter])
            axis.set_ylabel("relative WSSE increase")
            axis.set_ylim(bottom=-0.005)
            axis.grid(alpha=0.2)
            axis.legend(frameon=False, fontsize=7)
    fig.suptitle("Practical identifiability of gradual activation parameters")
    fig.tight_layout()
    if save:
        PLOT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(PLOT_DIR / "activation_identifiability_profiles.png", dpi=180, bbox_inches="tight")
    return fig


def run_analysis(n_starts: int = 5, max_nfev: int = 300, seed: int = 20260812) -> dict[str, object]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    batches, theta_by_matrix, natural_tables = load_batches()
    chemical_support = chemical_support_table(batches)
    natural_nutrient_pulses = natural_tables["nutrient_pulses"].copy()
    nutrient_pulses = effective_nutrient_pulse_table(batches, natural_nutrient_pulses)
    (
        observations,
        inventory,
        sensor_qc,
        qc_summary,
        sampling_schedule,
        sensor_zero_offsets,
    ) = load_co2_data(natural_tables, chemical_support, nutrient_pulses)
    cache, driver_diagnostics = build_driver_cache(batches, theta_by_matrix, observations)
    temperature_inputs = temperature_model_input_table(cache)
    temperature_alignment = temperature_alignment_summary(observations, temperature_inputs)
    exclusions = excluded_experiment_table()

    fits: dict[str, dict[str, float]] = {}
    start_frames = []
    for matrix in ("synthetic", "natural"):
        parameters, starts = fit_matrix(
            matrix,
            observations,
            cache,
            n_starts=n_starts,
            max_nfev=max_nfev,
            seed=seed,
        )
        fits[matrix] = parameters
        start_frames.append(starts)

    fit_parameters = fit_parameter_table(fits)
    fit_starts = pd.concat(start_frames, ignore_index=True)
    predictions, batch_metrics, validation = predict_and_score(fits, observations, cache)
    (
        jacobian_identifiability,
        jacobian_singular_values,
        local_parameter_correlations,
    ) = jacobian_identifiability_diagnostics(fits, observations, cache)
    activation_profiles = profile_activation_parameters(fits, observations, cache)
    loo_activation_estimates, loo_activation_summary = leave_one_batch_out_stability(
        fits, observations, cache, seed
    )

    previous_fits: dict[str, dict[str, float]] = {}
    previous_start_frames = []
    for matrix in ("synthetic", "natural"):
        parameters, starts = fit_matrix(
            matrix,
            observations,
            cache,
            n_starts=max(3, min(n_starts, 5)),
            max_nfev=max_nfev,
            seed=seed,
            model_variant=THRESHOLD_RELEASE_MODEL_NAME,
        )
        previous_fits[matrix] = parameters
        previous_start_frames.append(starts)
    previous_fit_parameters = fit_parameter_table(previous_fits)
    previous_fit_starts = pd.concat(previous_start_frames, ignore_index=True)
    previous_predictions, previous_batch_metrics, previous_validation = predict_and_score(
        previous_fits, observations, cache, model_variant=THRESHOLD_RELEASE_MODEL_NAME
    )

    legacy_fits: dict[str, dict[str, float]] = {}
    legacy_start_frames = []
    for matrix in ("synthetic", "natural"):
        parameters, starts = fit_matrix(
            matrix,
            observations,
            cache,
            n_starts=max(3, min(n_starts, 5)),
            max_nfev=max_nfev,
            seed=seed,
            model_variant=LEGACY_MODEL_NAME,
        )
        legacy_fits[matrix] = parameters
        legacy_start_frames.append(starts)
    legacy_fit_parameters = fit_parameter_table(legacy_fits)
    legacy_fit_starts = pd.concat(legacy_start_frames, ignore_index=True)
    legacy_predictions, legacy_batch_metrics, legacy_validation = predict_and_score(
        legacy_fits, observations, cache, model_variant=LEGACY_MODEL_NAME
    )
    source_transition_diagnostics = transition_diagnostics(
        fits, observations, cache, model_variant=MODEL_NAME
    )
    previous_source_transition_diagnostics = transition_diagnostics(
        previous_fits, observations, cache, model_variant=THRESHOLD_RELEASE_MODEL_NAME
    )
    legacy_source_transition_diagnostics = transition_diagnostics(
        legacy_fits, observations, cache, model_variant=LEGACY_MODEL_NAME
    )
    comparison_keys = ["scenario", "calibration_matrix", "target_matrix", "batch", "experiment_code"]
    comparison_metrics = [
        "rmse_g_l_h",
        "nrmse_peak",
        "correlation",
        "r2",
        "integral_ratio_pred_over_observed_lower_bound",
        "observed_onset_h",
        "predicted_onset_h",
        "onset_delay_h",
        "predicted_first_emission_0p005_h",
        "observed_visual_rise_duration_h",
        "predicted_visual_rise_duration_h",
        "visual_rise_duration_error_h",
        "observed_postpulse_peak_h",
        "predicted_postpulse_peak_h",
        "observed_pulse_t_rise_h",
        "predicted_pulse_t_rise_h",
        "postpulse_peak_delay_h",
    ]
    model_comparison = previous_validation[comparison_keys + comparison_metrics].merge(
        validation[comparison_keys + comparison_metrics],
        on=comparison_keys,
        suffixes=("_threshold_release", "_continuous_release"),
        validate="one_to_one",
    )
    model_comparison["rmse_change_continuous_minus_threshold_release"] = (
        model_comparison["rmse_g_l_h_continuous_release"]
        - model_comparison["rmse_g_l_h_threshold_release"]
    )
    model_comparison["abs_onset_error_change_continuous_minus_threshold_release"] = (
        model_comparison["onset_delay_h_continuous_release"].abs()
        - model_comparison["onset_delay_h_threshold_release"].abs()
    )
    model_comparison["abs_postpulse_peak_error_change_continuous_minus_threshold_release"] = (
        model_comparison["postpulse_peak_delay_h_continuous_release"].abs()
        - model_comparison["postpulse_peak_delay_h_threshold_release"].abs()
    )
    model_comparison["abs_visual_rise_duration_error_change_continuous_minus_threshold_release"] = (
        model_comparison["visual_rise_duration_error_h_continuous_release"].abs()
        - model_comparison["visual_rise_duration_error_h_threshold_release"].abs()
    )
    legacy_model_comparison = legacy_validation[comparison_keys + comparison_metrics].merge(
        validation[comparison_keys + comparison_metrics],
        on=comparison_keys,
        suffixes=("_legacy_slow_transition", "_continuous_release"),
        validate="one_to_one",
    )

    baseline_observations = observations.copy()
    baseline_observations["co2_rate_g_l_h"] = baseline_observations["co2_rate_raw_g_l_h"]
    baseline_observations["left_censored"] = False
    baseline_fits: dict[str, dict[str, float]] = {}
    for matrix in ("synthetic", "natural"):
        parameters, _ = fit_matrix(
            matrix,
            baseline_observations,
            cache,
            n_starts=max(3, min(n_starts, 5)),
            max_nfev=max_nfev,
            seed=seed,
            model_variant=MODEL_NAME,
        )
        baseline_fits[matrix] = parameters
    baseline_fit_parameters = fit_parameter_table(baseline_fits)
    baseline_predictions, _, baseline_validation = predict_and_score(
        baseline_fits, baseline_observations, cache
    )
    common_support = baseline_predictions.drop(
        columns=["observed_g_l_h", "left_censored", "detection_limit_g_l_h"]
    ).merge(
        observations[
            ["batch", "t_h", "co2_rate_g_l_h", "left_censored", "detection_limit_g_l_h"]
        ],
        left_on=["batch", "time_h"],
        right_on=["batch", "t_h"],
        how="left",
        validate="many_to_one",
    )
    common_support_rows = []
    for row in baseline_validation.itertuples(index=False):
        group = common_support[
            common_support["calibration_matrix"].eq(row.calibration_matrix)
            & common_support["batch"].eq(row.batch)
        ].sort_values("time_h")
        common_metric = _metric_row(
            row.calibration_matrix,
            row.target_matrix,
            row.batch,
            group["co2_rate_g_l_h"].to_numpy(dtype=float),
            group["predicted_g_l_h"].to_numpy(dtype=float),
            group["time_h"].to_numpy(dtype=float),
            row.role,
            group["left_censored"].astype(bool).to_numpy(),
            group["detection_limit_g_l_h"].to_numpy(dtype=float),
            MODEL_NAME,
            cache[row.batch].n_pulse_time_h,
        )
        common_metric["scenario"] = row.scenario
        common_support_rows.append(common_metric)
    baseline_common_support_validation = pd.DataFrame(common_support_rows)
    comparison_keys = ["calibration_matrix", "target_matrix", "batch", "experiment_code"]
    filter_impact = baseline_common_support_validation[
        comparison_keys + ["rmse_g_l_h", "nrmse_peak", "bias_g_l_h", "correlation", "r2"]
    ].merge(
        validation[comparison_keys + ["rmse_g_l_h", "nrmse_peak", "bias_g_l_h", "correlation", "r2"]],
        on=comparison_keys,
        suffixes=("_unfiltered_fit_common_support", "_filtered_censored"),
    )
    filter_impact["rmse_change"] = (
        filter_impact["rmse_g_l_h_filtered_censored"]
        - filter_impact["rmse_g_l_h_unfiltered_fit_common_support"]
    )

    observations.to_csv(RESULTS_DIR / "co2_observations_hourly.csv", index=False)
    sensor_qc.to_csv(RESULTS_DIR / "co2_sensor_qc_10min.csv", index=False)
    qc_summary.to_csv(RESULTS_DIR / "co2_sensor_qc_summary.csv", index=False)
    sensor_zero_offsets.to_csv(RESULTS_DIR / "sensor_zero_offsets.csv", index=False)
    sampling_schedule.to_csv(RESULTS_DIR / "sampling_schedule_used.csv", index=False)
    chemical_support.to_csv(RESULTS_DIR / "chemical_process_support.csv", index=False)
    natural_nutrient_pulses.to_csv(RESULTS_DIR / "natural_nutrient_pulse_details.csv", index=False)
    nutrient_pulses.to_csv(RESULTS_DIR / "effective_nutrient_pulses.csv", index=False)
    exclusions.to_csv(RESULTS_DIR / "excluded_experiments.csv", index=False)
    temperature_inputs.to_csv(RESULTS_DIR / "temperature_model_inputs.csv", index=False)
    temperature_alignment.to_csv(RESULTS_DIR / "temperature_alignment_summary.csv", index=False)
    inventory.to_csv(RESULTS_DIR / "experiment_inventory.csv", index=False)
    driver_diagnostics.to_csv(RESULTS_DIR / "driver_diagnostics.csv", index=False)
    fit_parameters.to_csv(RESULTS_DIR / "fit_parameters.csv", index=False)
    fit_starts.to_csv(RESULTS_DIR / "fit_start_diagnostics.csv", index=False)
    jacobian_identifiability.to_csv(
        RESULTS_DIR / "activation_jacobian_identifiability.csv", index=False
    )
    jacobian_singular_values.to_csv(
        RESULTS_DIR / "activation_jacobian_singular_values.csv", index=False
    )
    local_parameter_correlations.to_csv(
        RESULTS_DIR / "activation_local_parameter_correlations.csv", index=False
    )
    activation_profiles.to_csv(
        RESULTS_DIR / "activation_objective_profiles.csv", index=False
    )
    loo_activation_estimates.to_csv(
        RESULTS_DIR / "activation_leave_one_batch_out_estimates.csv", index=False
    )
    loo_activation_summary.to_csv(
        RESULTS_DIR / "activation_leave_one_batch_out_summary.csv", index=False
    )
    previous_fit_parameters.to_csv(
        RESULTS_DIR / "threshold_release_fit_parameters.csv", index=False
    )
    previous_fit_starts.to_csv(
        RESULTS_DIR / "threshold_release_fit_starts.csv", index=False
    )
    previous_predictions.to_csv(
        RESULTS_DIR / "threshold_release_prediction_rows.csv", index=False
    )
    previous_batch_metrics.to_csv(
        RESULTS_DIR / "threshold_release_batch_metrics.csv", index=False
    )
    previous_validation.to_csv(
        RESULTS_DIR / "threshold_release_validation_summary.csv", index=False
    )
    legacy_fit_parameters.to_csv(RESULTS_DIR / "legacy_slow_transition_fit_parameters.csv", index=False)
    legacy_fit_starts.to_csv(RESULTS_DIR / "legacy_slow_transition_fit_starts.csv", index=False)
    legacy_predictions.to_csv(RESULTS_DIR / "legacy_slow_transition_prediction_rows.csv", index=False)
    legacy_batch_metrics.to_csv(RESULTS_DIR / "legacy_slow_transition_batch_metrics.csv", index=False)
    legacy_validation.to_csv(RESULTS_DIR / "legacy_slow_transition_validation_summary.csv", index=False)
    source_transition_diagnostics.to_csv(
        RESULTS_DIR / "source_transition_diagnostics.csv", index=False
    )
    previous_source_transition_diagnostics.to_csv(
        RESULTS_DIR / "threshold_release_source_diagnostics.csv", index=False
    )
    legacy_source_transition_diagnostics.to_csv(
        RESULTS_DIR / "legacy_slow_transition_source_diagnostics.csv", index=False
    )
    model_comparison.to_csv(RESULTS_DIR / "model_comparison_validation.csv", index=False)
    legacy_model_comparison.to_csv(
        RESULTS_DIR / "legacy_vs_nitrogen_boost_validation.csv", index=False
    )
    predictions.to_csv(RESULTS_DIR / "prediction_rows.csv", index=False)
    batch_metrics.to_csv(RESULTS_DIR / "batch_metrics.csv", index=False)
    validation.to_csv(RESULTS_DIR / "validation_summary.csv", index=False)
    baseline_fit_parameters.to_csv(RESULTS_DIR / "baseline_unfiltered_fit_parameters.csv", index=False)
    baseline_validation.to_csv(RESULTS_DIR / "baseline_unfiltered_validation_summary.csv", index=False)
    baseline_common_support_validation.to_csv(
        RESULTS_DIR / "baseline_unfiltered_fit_common_support_validation.csv", index=False
    )
    filter_impact.to_csv(RESULTS_DIR / "filter_impact_validation.csv", index=False)

    plot_data_overview(observations, nutrient_pulses)
    plt.close()
    plot_sensor_zero_correction(sensor_qc, sensor_zero_offsets)
    plt.close()
    plot_sensor_filter_examples(sensor_qc, sampling_schedule, nutrient_pulses)
    plt.close()
    plot_temperature_profiles(observations, temperature_inputs, nutrient_pulses)
    plt.close()
    plot_process_timeline_alignment(inventory, nutrient_pulses)
    plt.close()
    plot_parameter_comparison(fit_parameters)
    plt.close()
    plot_activation_identifiability(activation_profiles, loo_activation_summary)
    plt.close()
    for matrix in ("synthetic", "natural"):
        plot_calibration_overlays(predictions, matrix)
        plt.close()
    plot_validation(predictions, validation)
    plt.close()
    plot_onset_model_comparison(batch_metrics, previous_batch_metrics)
    plt.close()
    plot_validation_model_comparison(
        predictions, previous_predictions, validation, previous_validation
    )
    plt.close()
    plot_initial_release_comparison(
        predictions,
        previous_predictions,
        batch_metrics,
        previous_batch_metrics,
    )
    plt.close()
    plot_pulse_response_comparison(
        predictions,
        previous_predictions,
        batch_metrics,
        previous_batch_metrics,
        source_transition_diagnostics,
    )
    plt.close()

    manifest = {
        "model": MODEL_NAME,
        "scope": "CO2 observation/transfer layer conditional on fixed medium-specific upstream fermentation drivers",
        "holdouts": HOLDOUTS,
        "observation_bin_h": OBSERVATION_BIN_H,
        "driver_grid_h": DRIVER_GRID_H,
        "sensor_preprocessing": {
            "sensor_zero_correction": {
                "native_unit": "sccm",
                "scope": "estimated independently for every experimental run before clipping, artifact filtering and volume conversion",
                "window_h": SENSOR_ZERO_WINDOW_H,
                "quantile": SENSOR_ZERO_QUANTILE,
                "nonnegative_offset": True,
                "assignment": "physical F1/F2/F3 acquisition channels; explicit DOE-F06 user label sensor 6 retained on channel F3",
                "caveat": "the initial low quantile is assumed to represent instrumental zero; true biological flow already present in the first 12 h may be partially removed",
            },
            "sample_artifact_window_h": SAMPLE_ARTIFACT_WINDOW_H,
            "sample_event_half_width_h": SAMPLE_EVENT_HALF_WIDTH_H,
            "sample_drop_ratio": SAMPLE_DROP_RATIO,
            "transient_context_h": TRANSIENT_CONTEXT_H,
            "transient_low_ratio": TRANSIENT_LOW_RATIO,
            "transient_high_ratio": TRANSIENT_HIGH_RATIO,
            "co2_detection_limit_standard_g_l_h": CO2_DETECTION_LIMIT_G_L_H,
            "co2_detection_limit_cold_g_l_h": COLD_CO2_DETECTION_LIMIT_G_L_H,
            "cold_setpoint_threshold_c": COLD_SETPOINT_THRESHOLD_C,
            "hourly_robust_median_window_points": SMOOTHING_MEDIAN_WINDOW_POINTS,
            "hourly_savgol_window_points": SMOOTHING_SAVGOL_WINDOW_POINTS,
            "hourly_savgol_polyorder": SMOOTHING_SAVGOL_POLYORDER,
            "censored_residual": "one-sided hinge above the detection limit; excluded from profiled gain and point metrics",
        },
        "excluded_batches": EXCLUDED_BATCHES,
        "temperature_use": {
            "core_dynamics": "BatchData.temperature_c via temperature_at(batch, t)",
            "co2_solubility": "Cstar base includes exp[-0.032*(temperature_c-20)]",
            "setpoint": "context/QC only; measured/reconstructed temperature is the model input",
            "lot1_correction": "full native sensor temperature grid merged into BatchData; not sparse chemistry-time interpolation",
        },
        "time_alignment": {
            "origin": "first chemical sample (t=0)",
            "end": "last chemical sample",
            "co2_policy": "CO2 rows outside the chemical interval are excluded; absent late CO2 is not imputed",
            "temperature_policy": "temperature is interpolated on the dynamic grid through the last chemical sample",
        },
        "natural_nutrient_protocol": {
            "calendar_source": str(NUTRIENT_CALENDAR_PATH),
            "initial_pulse": "ignored because initial YAN is already represented at t=0",
            "in_process_springferm_xtrem_g": SPRINGFERM_XTREM_G,
            "in_process_fda_g": FDA_G,
            "yan_mg_per_mg_product": YAN_MG_PER_MG_PRODUCT,
            "total_yan_mg_per_2L_reactor": float(
                (SPRINGFERM_XTREM_G + FDA_G) * 1000.0 * YAN_MG_PER_MG_PRODUCT
            ),
            "modeled_N_increment_kg_m3": 0.14,
            "density_target_g_l": NUTRIENT_DENSITY_TARGET_G_L,
            "timing_rule": "linear chemical-density crossing at 1040 g/L whenever a crossing is available",
            "timing_uncertainty": "crossings bracketed by more than 24 h are retained and labeled sparse_bracket; calendar is audit/fallback only",
        },
        "synthetic_nutrient_protocol": {
            "timing_rule": "existing DOE batch schedule expressed as process time from the first chemical sample",
            "calendar_or_density_retiming": False,
        },
        "n_starts": n_starts,
        "max_nfev": max_nfev,
        "random_seed": seed,
        "chemistry_aligned_transition": {
            "activation_evidence": "first chemical sample with glucose+fructose loss >=5 g/L or ethanol increase >=2 g/L",
            "activation_shape": "smoothstep; exactly zero before a fitted fractional start inside the chemistry bracket and gradual over a fitted multiple of bracket width",
            "estimated_by_matrix": [
                "chem_activation_start_fraction",
                "chem_activation_duration_fraction",
            ],
            "co2_used_to_locate_activation": "only within the chemistry-defined bracket through shared matrix parameters; no batch-specific lag",
            "onset_definition": "first three consecutive hourly values above early baseline plus 10% dynamic range, never below the row-specific detection limit",
            "onset_residual_scale_h": ONSET_SIGMA_H,
            "direct_comparator": THRESHOLD_RELEASE_MODEL_NAME,
            "legacy_comparator": LEGACY_MODEL_NAME,
        },
        "nitrogen_boost_transition": {
            "counterfactual": "simulate each batch with and without its in-process N pulse",
            "availability_shape": "causal linear ramp from pulse time to pulse time + fitted pulse_t_rise_h",
            "biomass_effect": "unscaled ramped difference between with-pulse and no-pulse upstream trajectories",
            "activity_effect": "pulse_activity_gain is the final multiplier on no-pulse metabolic activity; the term is independent of pulse-induced biomass growth",
            "estimated_by_matrix": ["pulse_t_rise_h", "pulse_activity_gain"],
            "postpulse_peak_window_h": PULSE_PEAK_WINDOW_H,
            "postpulse_peak_residual_scale_h": PULSE_PEAK_SIGMA_H,
            "literature_basis": [
                "David et al.: finite nitrogen uptake interval and transporter/activity state after addition",
                "Seguinot et al.: immediate activity boost after stationary-phase nitrogen addition without requiring biomass growth",
            ],
        },
        "continuous_co2_release": {
            "previous_problem": "qgas was driven only by dissolved CO2 above the full pure-CO2 saturation scale, producing a zero plateau and abrupt release",
            "new_flux": "k_release * C_dissolved * [floor + (1-floor)*C_dissolved/(C_dissolved+Cstar)]",
            "fixed_early_release_floor": CO2_CONTINUOUS_RELEASE_FLOOR,
            "first_emission_diagnostic_g_l_h": EARLY_EMISSION_THRESHOLD_G_L_H,
            "direct_comparator": THRESHOLD_RELEASE_MODEL_NAME,
        },
        "practical_identifiability": {
            "local_test": "SVD and covariance of the equal-batch residual Jacobian with respect to log shape parameters",
            "global_direction_test": "one-dimensional objective profiles with all other shape parameters reoptimized",
            "stability_test": "leave-one-calibration-batch-out refits initialized from the full-data optimum",
            "profile_points_per_parameter_plus_optimum": IDENTIFIABILITY_PROFILE_POINTS,
            "profile_relative_wsse_heuristic": 0.05,
            "warning": "the 5% objective line is a sensitivity heuristic, not a likelihood confidence interval",
        },
        "natural_requested": [f"LAB{value:03d}" for value in range(1, 13)],
        "natural_model_ready": sorted(inventory.loc[inventory["matrix"].eq("natural") & inventory["calibratable"], "batch"].tolist()),
        "natural_descriptive_only": sorted(inventory.loc[inventory["matrix"].eq("natural") & ~inventory["calibratable"], "batch"].tolist()),
        "source_files": {
            "synthetic_lot1_co2": str(LOT1_CO2_PATH.relative_to(ROOT_DIR)),
            "synthetic_lot2_co2": str(LOT2_CO2_PATH.relative_to(ROOT_DIR)),
            "natural_workbook": str(natural_loader.first_existing_path(natural_loader.NATURAL_WORKBOOK_CANDIDATES).relative_to(ROOT_DIR)),
            "natural_process_dir": str(RAW_NATURAL_DIR.relative_to(ROOT_DIR)),
            "natural_theta": str(NATURAL_THETA_PATH.relative_to(ROOT_DIR)),
            "synthetic_theta": str(SYNTHETIC_THETA_PATH.relative_to(ROOT_DIR)),
            "natural_nutrient_calendar": str(NUTRIENT_CALENDAR_PATH),
        },
        "oxygen_transition": {
            "estimated_by_matrix": ["O2_qmax_mg_gdw_h", "O2_initial_scale"],
            "fixed_O2_K_mg_l": O2_K_MG_L,
            "fixed_O2_ana_K_mg_l": O2_ANA_K_MG_L,
            "fixed_O2_ana_hill": O2_ANA_HILL,
            "fixed_O2_crabtree_floor": O2_CRABTREE_FLOOR,
            "O2_kLa_h": 0.0,
        },
    }
    (RESULTS_DIR / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {
        "observations": observations,
        "sensor_qc": sensor_qc,
        "sensor_zero_offsets": sensor_zero_offsets,
        "qc_summary": qc_summary,
        "sampling_schedule": sampling_schedule,
        "chemical_support": chemical_support,
        "natural_nutrient_pulses": natural_nutrient_pulses,
        "nutrient_pulses": nutrient_pulses,
        "exclusions": exclusions,
        "temperature_inputs": temperature_inputs,
        "temperature_alignment": temperature_alignment,
        "inventory": inventory,
        "driver_diagnostics": driver_diagnostics,
        "fits": fits,
        "fit_parameters": fit_parameters,
        "fit_starts": fit_starts,
        "jacobian_identifiability": jacobian_identifiability,
        "jacobian_singular_values": jacobian_singular_values,
        "local_parameter_correlations": local_parameter_correlations,
        "activation_profiles": activation_profiles,
        "loo_activation_estimates": loo_activation_estimates,
        "loo_activation_summary": loo_activation_summary,
        "previous_fits": previous_fits,
        "previous_fit_parameters": previous_fit_parameters,
        "previous_fit_starts": previous_fit_starts,
        "previous_predictions": previous_predictions,
        "previous_batch_metrics": previous_batch_metrics,
        "previous_validation": previous_validation,
        "legacy_fits": legacy_fits,
        "legacy_fit_parameters": legacy_fit_parameters,
        "legacy_fit_starts": legacy_fit_starts,
        "legacy_predictions": legacy_predictions,
        "legacy_batch_metrics": legacy_batch_metrics,
        "legacy_validation": legacy_validation,
        "source_transition_diagnostics": source_transition_diagnostics,
        "previous_source_transition_diagnostics": previous_source_transition_diagnostics,
        "legacy_source_transition_diagnostics": legacy_source_transition_diagnostics,
        "model_comparison": model_comparison,
        "legacy_model_comparison": legacy_model_comparison,
        "predictions": predictions,
        "batch_metrics": batch_metrics,
        "validation": validation,
        "baseline_fit_parameters": baseline_fit_parameters,
        "baseline_validation": baseline_validation,
        "baseline_common_support_validation": baseline_common_support_validation,
        "filter_impact": filter_impact,
        "manifest": manifest,
    }


def _create_notebook_legacy() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb["metadata"]["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb["metadata"]["language_info"] = {"name": "python", "version": f"{sys.version_info.major}.{sys.version_info.minor}"}
    cells = []
    cells.append(
        nbformat.v4.new_markdown_cell(
            """# CO₂ 2026: calibración por matriz y validación cruzada

## TL;DR

Este notebook calibra la capa **`solubility_o2_slow_transition`** por separado con mosto sintético (MBDOE lotes 1–2) y mosto natural, y luego transfiere los parámetros sin reajuste entre matrices.

- Holdout sintético predefinido: **DOE-F06 (`lot2_F3`)**.
- Holdout natural predefinido: **LAB012**.
- Se excluyen completamente **LAB001, LAB002, LAB003, LAB009 y DOE-F02 (`lot1_F2`)** por perfiles no confiables/problemas de implementación. No participan en ajuste, validación ni figuras de perfiles.
- En torno a una toma de muestra se reconstruye el evento completo de ±1.25 h sólo si existe una caída y los niveles pre/post concuerdan. Luego se aplica mediana robusta de 3 h y Savitzky–Golay cuadrático de 5 h; la señal cruda y la reconstruida quedan guardadas para auditoría.
- Las lecturas bajo 0.05 g CO₂ L⁻¹ h⁻¹ se tratan como censura izquierda. En operación fría (setpoint ≤15.5 °C) se reportan además como intervalos de baja sensibilidad.
- La temperatura **sí entra al modelo**: alimenta la cinética upstream y la solubilidad de CO₂. Se grafica el perfil efectivamente usado junto con temperatura medida y setpoint.
- La primera y la última muestra química definen el intervalo del proceso. CO₂ fuera de ese intervalo se descarta; si el sensor termina antes, el tramo faltante no se inventa.
- LAB004–LAB012 reciben un pulso efectivo de **0.14 kg N m⁻³** (280 mg YAN en 2 L). El pulso inicial se omite porque ya está representado en el punto químico t=0.
- La validación es de la **capa CO₂ condicionada por los drivers de fermentación existentes**. No constituye todavía una validación end-to-end de todos los parámetros cinéticos.
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Contexto y método

El modelo mantiene la estructura seleccionada previamente:

\[
C^*_{CO_2}=s_{sat}\,1.69\,e^{-0.032(T-20)}e^{0.0016E}e^{-0.0012(G+F)}
\]

\[
\frac{dO_2}{dt}=k_La(O_2^*-O_2)-q_{O_2,max}X\frac{O_2}{K_{O_2}+O_2},\qquad
\phi_{ana}=\frac{K_{ana}^{h}}{K_{ana}^{h}+O_2^{h}}
\]

\[
\frac{dC_{CO_2}}{dt}=q_{prod}-q_{gas},\qquad
q_{gas}=k_{release}\,\operatorname{softplus}(C_{CO_2}-C^*_{CO_2})
\]

Se estiman por matriz `kCO2_release_h`, `CO2sat_scale` y una única `matrix_gain`. Esta última reemplaza la escala libre por fermentación del benchmark piloto y queda congelada durante todas las validaciones. Las señales se convierten a `g CO2 L^-1 h^-1`, se agregan por mediana cada 1 h y cada fermentación recibe igual peso en la función objetivo.

Antes de agregar se reemplazan excursiones cortas cuyo nivel previo y posterior es consistente. En cada muestreo registrado se evalúa una ventana central de ±1.25 h contra medianas pre/post separadas: si se confirma una caída, se reconstruye toda la ventana para eliminar conjuntamente caída y rebote. Sobre la mediana horaria se aplica una mediana móvil centrada de 3 puntos y un Savitzky–Golay de 5 puntos, orden 2. Se reportan por fermentación los cambios de integral, peak y rugosidad. Los puntos bajo el límite de detección no se ajustan como ceros: aportan sólo penalización unilateral cuando el modelo supera dicho límite y se excluyen de la estimación de `matrix_gain` y de las métricas punto a punto.

Los estados upstream (`X`, `G`, `F`, `E`, temperatura y producción biológica) se calculan con las calibraciones 2026 ya existentes para cada matriz. `BatchData.temperature_c` se interpola con `temperature_at(batch,t)` y además entra explícitamente en `C*CO₂` mediante `exp[-0.032(T-20)]`. Se usa temperatura medida/reconstruida; el setpoint se conserva como referencia/QC. Por tanto, el holdout evita fuga de **observaciones CO₂**, pero no pretende reestimar/validar nuevamente el bloque cinético completo.

La escala temporal se ancla en la primera muestra química (t=0) y termina en la última. Los timestamps del calendario se convierten de UTC a `America/Santiago` antes de restar t=0. Para LAB004–LAB009 el tiempo del pulso se obtiene por interpolación lineal del cruce de densidad 1040 g/L cuando el bracket químico es ≤24 h; con brackets más amplios (LAB010–LAB012) se conserva el timestamp del calendario. Ambas referencias quedan en la tabla de auditoría.
"""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import sys
import pandas as pd
from IPython import get_ipython
from IPython.display import display

get_ipython().run_line_magic("matplotlib", "inline")
import matplotlib.pyplot as plt

ROOT = Path.cwd()
while ROOT.name != "pyomo-doe" and ROOT.parent != ROOT:
    ROOT = ROOT.parent
if ROOT.name != "pyomo-doe":
    raise RuntimeError("Execute this notebook from the repository or a descendant directory")

FM = ROOT / "fermentation_model"
for path in (FM, FM / "laboratory_2026", FM / "pilot_2025"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from laboratory_2026 import run_co2_matrix_cross_validation_2026 as analysis

pd.set_option("display.max_columns", 30)
pd.set_option("display.width", 160)
print("Repository:", ROOT)
print("Model:", analysis.MODEL_NAME)
print("Holdouts:", analysis.HOLDOUTS)
print("Excluded:", analysis.EXCLUDED_BATCHES)
"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Datos y partición experimental"))
    cells.append(
        nbformat.v4.new_code_cell(
            """result = analysis.run_analysis(n_starts=5, max_nfev=300, seed=20260812)
inventory = result["inventory"]
display(result["exclusions"])
display(inventory[[
    "matrix", "experiment_code", "batch", "lot", "chemistry_available",
    "co2_available", "calibratable", "split", "n_co2_hourly",
    "chemistry_first_h", "chemistry_last_h", "co2_source_last_h",
    "co2_first_h", "co2_last_h", "co2_end_minus_chemistry_h",
    "co2_peak_g_l_h", "median_setpoint_c",
    "n_artifacts_replaced", "n_left_censored_hourly", "note"
]])
display(result["natural_nutrient_pulses"][[
    "batch", "calendar_t_h", "density_crossing_t_h", "density_bracket_width_h",
    "model_pulse_time_h", "timing_source", "amount_N_kg_m3",
    "calendar_product", "calendar_dose", "excluded_from_co2_analysis"
]])
display(result["nutrient_pulses"][[
    "matrix", "experiment_code", "batch", "pulse_time_h", "amount_N_kg_m3",
    "timing_source"
]])
fig = analysis.plot_process_timeline_alignment(
    result["inventory"], result["nutrient_pulses"], save=False
)
plt.show()
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """El split se hace por fermentación completa. No se mezclan puntos de una misma curva entre calibración y validación. Los experimentos descartados quedan únicamente en `excluded_experiments.csv`; no ingresan a observaciones, ajuste, métricas ni figuras."""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["qc_summary"][[
    "matrix", "experiment_code", "batch", "median_setpoint_c",
    "n_artifacts_replaced", "n_sampling_window_transients",
    "n_other_short_transients", "negative_signed_fraction",
    "n_left_censored_hourly", "n_cold_low_sensitivity_hourly",
    "smoothing_roughness_ratio", "smoothing_integral_ratio", "smoothing_peak_ratio"
]])
fig = analysis.plot_data_overview(
    result["observations"], result["nutrient_pulses"], save=False
)
plt.show()
fig = analysis.plot_sensor_filter_examples(
    result["sensor_qc"], result["sampling_schedule"], save=False
)
plt.show()
"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Temperatura usada por el modelo"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["driver_diagnostics"][[
    "matrix", "batch", "temperature_used_min_c", "temperature_used_mean_c",
    "temperature_used_max_c", "csat_base_min_g_l", "csat_base_max_g_l"
]])
display(result["temperature_alignment"])
fig = analysis.plot_temperature_profiles(
    result["observations"], result["temperature_inputs"],
    result["nutrient_pulses"], save=False
)
plt.show()
"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Resultados de calibración"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["fit_parameters"])
best_starts = (
    result["fit_starts"].sort_values("wsse_equal_batch")
    .groupby("matrix", as_index=False).first()
)
display(best_starts[["matrix", "success", "nfev", "wsse_equal_batch", "kCO2_release_h", "CO2sat_scale", "matrix_gain"]])
fig = analysis.plot_parameter_comparison(result["fit_parameters"], save=False)
plt.show()
for matrix in ("synthetic", "natural"):
    fig = analysis.plot_calibration_overlays(result["predictions"], matrix, save=False)
    plt.show()
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """Un parámetro marcado `active_bound=True` llegó al límite permitido por la definición original del candidato. Eso debe interpretarse como señal de tensión estructural/identificabilidad, no como una estimación interior bien resuelta."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Validación interna y transferencia cruzada"))
    cells.append(
        nbformat.v4.new_markdown_cell(
            """`filter_impact` compara ambos ajustes sobre exactamente los mismos puntos cuantificables. El ajuste base se estima con la señal cruda, pero para esta comparación se evalúa contra la señal filtrada y se excluyen los mismos puntos censurados; así el cambio de RMSE no proviene de alterar el denominador."""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """columns = [
    "scenario", "calibration_matrix", "target_matrix", "experiment_code", "batch",
    "n", "n_total", "n_left_censored", "rmse_g_l_h", "nrmse_peak",
    "bias_g_l_h", "correlation", "r2", "integral_ratio_pred_over_observed_lower_bound",
    "observed_peak_g_l_h", "predicted_peak_g_l_h"
]
display(result["validation"][columns])
display(result["filter_impact"])
fig = analysis.plot_validation(result["predictions"], result["validation"], save=False)
plt.show()
"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Diagnóstico por fermentación"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(
    result["batch_metrics"][[
        "calibration_matrix", "target_matrix", "experiment_code", "batch", "role",
        "n", "n_left_censored", "rmse_g_l_h", "nrmse_peak", "bias_g_l_h",
        "correlation", "integral_ratio_pred_over_observed_lower_bound"
    ]].sort_values(["calibration_matrix", "target_matrix", "batch"])
)
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """## Sensibilidad a baja temperatura y experimento de verificación

El umbral de **0.05 g CO₂ L⁻¹ h⁻¹** es operacional y queda explícito en el manifiesto; no debe confundirse con un límite metrológico certificado. Hasta medirlo, los tramos bajo ese nivel se consideran **censura izquierda**, no producción nula. Esto afecta especialmente a LAB010 (setpoint mediano 15 °C).

Para reemplazar el supuesto por evidencia, ejecutar un ensayo de sensor a **15, 18 y 21 °C**, usando el mismo medio/volumen y una referencia independiente (pérdida gravimétrica o analizador de gas calibrado). En cada temperatura aplicar escalones de flujo que cubran 0.01–0.30 g CO₂ L⁻¹ h⁻¹, con blancos y réplicas, y estimar: límite de detección, sesgo, repetibilidad y tiempo de respuesta. El criterio de inclusión posterior debe definirse por temperatura: dato cuantitativo sobre el LOQ, censurado entre LOD–LOQ y no informativo si el sensor no responde al escalón. El notebook se reejecuta cambiando `CO2_DETECTION_LIMIT_G_L_H` y, si corresponde, incorporando un LOQ dependiente de temperatura.
"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Conclusiones reproducibles"))
    cells.append(
        nbformat.v4.new_code_cell(
            """validation = result["validation"].copy()
for row in validation.itertuples(index=False):
    direction = f"{row.calibration_matrix} → {row.target_matrix}"
    print(
        f"{direction:22s} | {row.experiment_code:7s} | "
        f"RMSE={row.rmse_g_l_h:.3f} g/L/h | NRMSEpeak={row.nrmse_peak:.2f} | "
        f"r={row.correlation:.2f} | n_censored={row.n_left_censored} | "
        f"integral pred/lower-bound={row.integral_ratio_pred_over_observed_lower_bound:.2f}"
    )

active = result["fit_parameters"].query("active_bound")
if not active.empty:
    print("\\nAdvertencia: parámetros en límite:")
    display(active)

print("\\nLectura recomendada:")
print("- La transferencia se considera apoyada sólo si mantiene simultáneamente error, forma e integral en el holdout.")
print("- Las métricas punto a punto excluyen censura izquierda; la integral observada sigue siendo sólo una cota inferior.")
print("- Un buen ajuste de calibración no sustituye el resultado retenido.")
print("- LAB001-LAB003, LAB009 y DOE-F02 están totalmente excluidos de esta corrida.")
print("- La temperatura usada por el modelo y el setpoint quedan auditados en CSV y figura.")
print("- La analítica química fija t=0 y el final; la cobertura CO2 efectiva queda auditada por lote.")
print("- Los pulsos LAB usan 0.14 kg N/m3 y el filtro/suavizado no cruza la discontinuidad del pulso.")
print("- Todos los CSV, figuras y el manifiesto quedaron en:", analysis.RESULTS_DIR.relative_to(ROOT))
"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """## Artefactos

La ejecución guarda observaciones horarias, soporte químico, pulsos efectivos y su detalle de calendario/densidad, inventario, parámetros, diagnósticos de multistart, predicciones, métricas por lote, resumen de validación, figuras y un manifiesto JSON con fuentes y supuestos. Esto permite auditar el resultado sin depender únicamente de las salidas embebidas del notebook.
"""
        )
    )
    nb["cells"] = cells
    nbformat.write(nb, NOTEBOOK_PATH)


def create_notebook() -> None:
    """Create the reproducible notebook for continuous CO2 release and N response."""

    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    nb["metadata"]["language_info"] = {
        "name": "python",
        "version": f"{sys.version_info.major}.{sys.version_info.minor}",
    }
    cells: list[object] = []
    cells.append(
        nbformat.v4.new_markdown_cell(
            """# CO₂ 2026: cero de sensor, inicio gradual y validación cruzada

## Resultado que se evalúa

Esta corrida corrige dos fuentes de la cola inicial artificial: la activación casi escalonada de la fuente biológica y la liberación gaseosa nula hasta superar la saturación. La comparación primaria conserva la respuesta de N pero usa la liberación por umbral anterior; el modelo lento legado se conserva como referencia histórica.

- Antes de filtrar, cada corrida recibe una corrección de cero en sccm estimada con el percentil 10 de sus primeras 12 h. El cero no se comparte entre campañas.
- Los archivos exponen canales físicos F1–F3. DOE-F06 conserva la etiqueta solicitada `sensor 6`, aunque se adquiere por F3; el offset se estima por corrida y no depende de agrupar ambas etiquetas.
- La química define el intervalo admisible de activación. Dos parámetros compartidos por matriz sitúan el comienzo y la duración de una rampa suave dentro de esa evidencia; no se estima un retardo libre por lote.
- `O2_qmax_mg_gdw_h` y `O2_initial_scale` se estiman por matriz junto con transferencia, solubilidad y ganancia.
- La fase gaseosa puede recibir un flujo pequeño desde que existe CO₂ disuelto; ya no se exige superar (C^*_{CO_2}) para abandonar exactamente cero.
- El efecto del pulso se construye como la diferencia causal entre simulaciones upstream con y sin adición de N. `pulse_t_rise_h` distribuye su disponibilidad en el tiempo y `pulse_activity_gain` representa capacidad metabólica adicional sin convertirla automáticamente en biomasa.
- En LAB el pulso se ubica por cruce de densidad 1040 g/L; en DOE-F0X se respeta el tiempo de proceso ya registrado en cada fermentación.
- Todos los gráficos de perfil de CO₂ marcan el instante del pulso nutricional.
- A setpoint ≤15.5 °C se usa un umbral operacional más conservador de 0.10 g CO₂ L⁻¹ h⁻¹; en el resto, 0.05.
- Se calibran por separado mosto sintético y natural. DOE-F06 y LAB012 son holdouts completos; después se cruzan los modelos sin reajuste.
- Se excluyen LAB001–LAB003, LAB009 y DOE-F02 por perfiles no confiables.

La analítica química sigue definiendo t=0 y el término del proceso. Esta es una validación de la capa CO₂ condicionada por los drivers upstream disponibles, no una recalibración end-to-end del modelo de fermentación."""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            r"""## Estructura del modelo

Para la corrida (r), primero se corrige el cero instrumental en la señal nativa:

\[
b_r=\max\left[0,Q_{0.10}\{y_{sccm}(t):t\leq t_0+12\,h\}\right],\qquad
y_{corr}(t)=\max[y_{sccm}(t)-b_r,0]
\]

La conversión a g L⁻¹ h⁻¹, el filtrado de transientes y el suavizado ocurren después. Esta regla supone que el extremo inferior de las primeras 12 h representa cero instrumental; si ya existe flujo biológico durante toda esa ventana, puede sustraer parte de la señal real.

La primera muestra química que cumple pérdida de glucosa+fructosa ≥5 g/L o aumento de etanol ≥2 g/L marca el extremo superior del intervalo de activación. La muestra química anterior fija el extremo inferior. La fuente usa una rampa *smoothstep* continua:

\[
t_s=t_L+f_s(t_U-t_L),\quad
\Delta t=f_d(t_U-t_L),\quad
a_{chem}=3z^2-2z^3,\quad
z=\operatorname{clip}\!\left(\frac{t-t_s}{\Delta t},0,1\right)
\]

Los parámetros (f_s) y (f_d) son compartidos por matriz. Esto permite un aumento pequeño y gradual desde (t_s), pero mantiene el comienzo ligado al bracket químico en vez de introducir una latencia independiente por fermentación.

\[
\frac{dO_2}{dt}=-q_{O_2,max}X\frac{O_2}{K_{O_2}+O_2},\qquad
\phi_{ana}=\frac{K_{ana}^{h}}{K_{ana}^{h}+O_2^{h}}
\]

\[
q_{prod}=a_{chem}(t)\left[q_{bio}\left(f_{Crabtree}+(1-f_{Crabtree})\phi_{ana}\right)+q_{resp}\right]
\]

Para una adición en \(t_N\), la fracción utilizada es causal:

\[
f_N(t)=\operatorname{clip}\left(\frac{t-t_N}{t_{rise}},0,1\right),\qquad
q_{bio,N}=\left[1+(g_N-1)f_N(t)\right]q_{bio,0}
+f_N(t)\left(q_{bio,+N}-q_{bio,0}\right)
\]

La biomasa usa la misma diferencia con/sin pulso, pero sin \(g_N\). Así, el término de actividad puede capturar el aumento de capacidad de transporte descrito por el grupo de Sablayrolles sin imponer crecimiento ficticio.

\[
C^*_{CO_2}=s_{sat}\,1.69\,e^{-0.032(T-20)}e^{0.0016E}e^{-0.0012(G+F)}
\]

\[
\frac{dC_{CO_2}}{dt}=q_{prod}-q_{gas},\qquad
q_{gas}=\min\!\left[k_{release}C_{CO_2}
\left(f_0+(1-f_0)\frac{C_{CO_2}}{C_{CO_2}+C^*_{CO_2}}\right),
\frac{C_{CO_2}}{\Delta t}+q_{prod}\right],\quad f_0=0.05
\]

La formulación anterior, usada como comparador directo, tenía (q_{gas}=k_{release}\operatorname{softplus}(C_{CO_2}-C^*_{CO_2})): eso mantenía la predicción pegada a cero hasta acumular suficiente CO₂ disuelto. El nuevo piso (f_0) representa liberación sub-saturada pequeña y evita ese umbral duro.

El inicio observado se define como el primero de tres puntos horarios consecutivos sobre el mayor valor entre: límite de detección local y línea base inicial +10 % del rango dinámico. Además se reporta la primera emisión sostenida sobre 0.005 g L⁻¹ h⁻¹ y la duración visual 2–10 % del ascenso. El ajuste conserva el residuo de tiempo de inicio y, en lotes pulsados, el desfase del máximo durante las 72 h posteriores a la adición, además de los residuos del perfil completo."""
        )
    )
    cells.append(
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import sys
import pandas as pd
from IPython import get_ipython
from IPython.display import display

get_ipython().run_line_magic("matplotlib", "inline")
import matplotlib.pyplot as plt

ROOT = Path.cwd()
while ROOT.name != "pyomo-doe" and ROOT.parent != ROOT:
    ROOT = ROOT.parent
if ROOT.name != "pyomo-doe":
    raise RuntimeError("Execute this notebook from the repository or a descendant")

FM = ROOT / "fermentation_model"
for path in (FM, FM / "laboratory_2026", FM / "pilot_2025"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from laboratory_2026 import run_co2_matrix_cross_validation_2026 as analysis

pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 180)
print("Repository:", ROOT)
print("Current model:", analysis.MODEL_NAME)
print("Direct comparator:", analysis.THRESHOLD_RELEASE_MODEL_NAME)
print("Historical comparator:", analysis.LEGACY_MODEL_NAME)
print("Holdouts:", analysis.HOLDOUTS)
print("Excluded:", analysis.EXCLUDED_BATCHES)"""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Ejecución completa y partición experimental"))
    cells.append(
        nbformat.v4.new_code_cell(
            """result = analysis.run_analysis(n_starts=5, max_nfev=300, seed=20260812)

display(result["exclusions"])
display(result["inventory"][[
    "matrix", "experiment_code", "batch", "lot", "calibratable", "split",
    "acquisition_channel", "sensor_id", "sensor_zero_offset_sccm",
    "n_co2_hourly", "chemistry_first_h", "chemistry_last_h", "co2_first_h",
    "co2_last_h", "co2_peak_g_l_h", "median_setpoint_c",
    "n_artifacts_replaced", "n_left_censored_hourly", "note"
]])
display(result["natural_nutrient_pulses"][[
    "batch", "calendar_t_h", "density_crossing_t_h", "density_bracket_width_h",
    "model_pulse_time_h", "timing_source", "amount_N_kg_m3",
    "calendar_product", "calendar_dose", "excluded_from_co2_analysis"
]])
fig = analysis.plot_process_timeline_alignment(
    result["inventory"], result["nutrient_pulses"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """La separación se hace por fermentación completa. Los pulsos LAB aportan 0.14 kg N m⁻³; el pulso inicial se ignora porque el YAN inicial ya está en t=0. En todos los LAB con cruce disponible, el tiempo del pulso se obtiene por interpolación del cruce de densidad 1040 g/L. Un bracket químico >24 h se conserva, pero se etiqueta como incertidumbre de tiempo. En DOE-F0X se usa el tiempo de proceso existente, sin recalcularlo por densidad ni calendario."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Filtrado y sensibilidad del sensor"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["qc_summary"][[
    "matrix", "experiment_code", "batch", "acquisition_channel", "sensor_id",
    "sensor_zero_offset_sccm", "sensor_zero_offset_g_l_h", "median_setpoint_c",
    "n_artifacts_replaced", "n_sampling_window_transients",
    "n_other_short_transients", "negative_signed_fraction",
    "n_left_censored_hourly", "n_cold_low_sensitivity_hourly",
    "smoothing_roughness_ratio", "smoothing_integral_ratio", "smoothing_peak_ratio"
]])
display(result["sensor_zero_offsets"])
fig = analysis.plot_sensor_zero_correction(
    result["sensor_qc"], result["sensor_zero_offsets"], save=False
)
plt.show()
fig = analysis.plot_data_overview(
    result["observations"], result["nutrient_pulses"], save=False
)
plt.show()
fig = analysis.plot_sensor_filter_examples(
    result["sensor_qc"], result["sampling_schedule"],
    result["nutrient_pulses"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """El primer gráfico comprueba la corrección de cero antes de cualquier reconstrucción. El offset se estima por corrida porque los canales muestran deriva entre campañas. Luego, las excursiones asociadas a muestreo se reconstruyen sólo cuando el evento completo está respaldado por niveles pre/post consistentes; se aplica mediana robusta de 3 h y Savitzky–Golay cuadrático de 5 h. Los puntos bajo el umbral local aportan una penalización unilateral si el modelo excede el límite."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Temperatura y evidencia química de activación"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["driver_diagnostics"][[
    "matrix", "batch", "temperature_used_min_c", "temperature_used_mean_c",
    "temperature_used_max_c", "chemical_activity_lower_h",
    "chemical_activity_upper_h", "chemical_activity_center_h",
    "initial_o2_saturation_base_mg_l", "base_qprod_peak_g_l_h",
    "csat_base_min_g_l", "csat_base_max_g_l"
]])
display(result["temperature_alignment"])
fig = analysis.plot_temperature_profiles(
    result["observations"], result["temperature_inputs"],
    result["nutrient_pulses"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """La temperatura medida/reconstruida entra tanto a los drivers cinéticos como a la solubilidad de CO₂; el setpoint queda como referencia. La tabla permite verificar, lote por lote, el intervalo químico que condiciona el encendido de la fuente."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Calibración por matriz"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["fit_parameters"])
best_starts = (
    result["fit_starts"].sort_values("wsse_equal_batch")
    .groupby("matrix", as_index=False).first()
)
display(best_starts[[
    "matrix", "success", "nfev", "wsse_equal_batch", "kCO2_release_h",
    "CO2sat_scale", "O2_qmax_mg_gdw_h", "O2_initial_scale",
    "pulse_t_rise_h", "pulse_activity_gain",
    "chem_activation_start_fraction", "chem_activation_duration_fraction",
    "matrix_gain"
]])
fig = analysis.plot_parameter_comparison(result["fit_parameters"], save=False)
plt.show()
for matrix in ("synthetic", "natural"):
    fig = analysis.plot_calibration_overlays(result["predictions"], matrix, save=False)
    plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """`active_bound=True` indica que el dato sólo acota el parámetro en el borde permitido. Debe leerse como tensión estructural o falta de identificabilidad, no como una estimación interior resuelta. `chem_activation_start_fraction` ubica el inicio dentro del bracket químico y `chem_activation_duration_fraction` escala la duración por el ancho de ese bracket. `pulse_activity_gain=1` significa ausencia de modulación sobre la actividad preexistente; >1 es boost y <1 es atenuación."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Estimabilidad de inicio y duración"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["jacobian_identifiability"])
display(
    result["local_parameter_correlations"].query(
        "parameter_1 in ['chem_activation_start_fraction', 'chem_activation_duration_fraction']"
    )[["matrix", "parameter_1", "parameter_2", "local_log_parameter_correlation"]]
    .sort_values(["matrix", "parameter_1", "local_log_parameter_correlation"])
)
display(result["loo_activation_summary"])
display(result["loo_activation_estimates"][[
    "matrix", "omitted_batch", "success", "nfev",
    "chem_activation_start_fraction", "chem_activation_duration_fraction"
]])
display(result["activation_profiles"])
fig = analysis.plot_activation_identifiability(
    result["activation_profiles"], result["loo_activation_summary"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """El Jacobiano evalúa sensibilidad local en log-parámetros; un número de condición alto indica direcciones compensables. Los perfiles fijan uno de los dos parámetros y reoptimizan todos los demás. La franja naranjo muestra cuánto cambia la estimación al retirar una fermentación completa. La línea horizontal de 5 % es sólo un umbral descriptivo de sensibilidad del objetivo, no un intervalo de confianza de verosimilitud."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## ¿Se corrigió la cola inicial?"))
    cells.append(
        nbformat.v4.new_code_cell(
            """native_transition = result["source_transition_diagnostics"].query("native_matrix_fit").copy()
display(native_transition[[
    "target_matrix", "experiment_code", "batch", "median_setpoint_c",
    "cold_operation", "chemical_activity_lower_h", "chemical_activity_upper_h",
    "source_10pct_peak_onset_h", "o2_below_anaerobic_halfpoint_h",
    "observed_onset_h", "predicted_onset_h", "onset_delay_h",
    "O2_qmax_mg_gdw_h", "O2_initial_scale", "n_pulse_time_h",
    "pulse_t_rise_h", "pulse_activity_gain",
    "chem_activation_start_fraction", "chem_activation_duration_fraction",
    "pulse_utilization_complete_h",
    "observed_postpulse_peak_h", "predicted_postpulse_peak_h",
    "postpulse_peak_delay_h"
]])
fig = analysis.plot_initial_release_comparison(
    result["predictions"], result["previous_predictions"],
    result["batch_metrics"], result["previous_batch_metrics"], save=False
)
plt.show()
fig = analysis.plot_onset_model_comparison(
    result["batch_metrics"], result["previous_batch_metrics"], save=False
)
plt.show()
fig = analysis.plot_pulse_response_comparison(
    result["predictions"], result["previous_predictions"],
    result["batch_metrics"], result["previous_batch_metrics"],
    result["source_transition_diagnostics"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """El primer gráfico es el control visual principal: amplía el inicio de cada lote y compara liberación por umbral (naranjo) con liberación continua (azul). La línea horizontal marca 0.005 g L⁻¹ h⁻¹; el texto cuantifica la primera emisión y la duración 2–10 %. El segundo resume el error de inicio. El tercero amplía cada lote pulsado: línea magenta = adición; línea morada = término de la disponibilidad gradual estimada."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Validación retenida y transferencia cruzada"))
    cells.append(
        nbformat.v4.new_code_cell(
            """validation_columns = [
    "scenario", "calibration_matrix", "target_matrix", "experiment_code", "batch",
    "n", "n_total", "n_left_censored", "rmse_g_l_h", "nrmse_peak", "bias_g_l_h",
    "correlation", "r2", "integral_ratio_pred_over_observed_lower_bound",
    "observed_onset_h", "predicted_onset_h", "onset_delay_h",
    "predicted_first_emission_0p005_h", "observed_visual_rise_duration_h",
    "predicted_visual_rise_duration_h", "visual_rise_duration_error_h",
    "pulse_time_h", "observed_postpulse_peak_h", "predicted_postpulse_peak_h",
    "postpulse_peak_delay_h"
]
display(result["validation"][validation_columns])
display(result["model_comparison"])
display(result["filter_impact"])
fig = analysis.plot_validation(result["predictions"], result["validation"], save=False)
plt.show()
fig = analysis.plot_validation_model_comparison(
    result["predictions"], result["previous_predictions"],
    result["validation"], result["previous_validation"], save=False
)
plt.show()"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """`model_comparison` usa exactamente los mismos holdouts para la liberación por umbral anterior y la liberación continua nueva, conservando la misma respuesta finita al pulso. Valores negativos en las columnas `*_change_continuous_minus_threshold_release` significan mejora. `filter_impact` compara los ajustes crudo y filtrado sobre el mismo soporte cuantificable."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Diagnóstico por fermentación"))
    cells.append(
        nbformat.v4.new_code_cell(
            """display(result["batch_metrics"][[
    "calibration_matrix", "target_matrix", "experiment_code", "batch", "role",
    "n", "n_left_censored", "rmse_g_l_h", "nrmse_peak", "bias_g_l_h",
    "correlation", "integral_ratio_pred_over_observed_lower_bound",
    "observed_onset_h", "predicted_onset_h", "onset_delay_h", "pulse_time_h",
    "predicted_first_emission_0p005_h", "observed_visual_rise_duration_h",
    "predicted_visual_rise_duration_h", "visual_rise_duration_error_h",
    "observed_postpulse_peak_h", "predicted_postpulse_peak_h",
    "postpulse_peak_delay_h"
]].sort_values(["calibration_matrix", "target_matrix", "batch"]))"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """## Experimento pendiente para distinguir sensor de biología

El umbral dependiente de temperatura sigue siendo operacional, no un LOD/LOQ certificado. Para identificarlo se requiere un ensayo a 15, 18 y 21 °C con el mismo medio y volumen, referencia independiente de CO₂ y escalones de 0.01–0.30 g CO₂ L⁻¹ h⁻¹. En cada temperatura deben estimarse LOD, LOQ, sesgo, repetibilidad y tiempo de respuesta.

La corrección por percentil inicial también es operacional. Debe reemplazarse por blancos de gas sin fermentación medidos al comienzo y al final de cada corrida en cada canal. Eso separaría offset, deriva y flujo biológico temprano sin depender de una ventana temporal elegida.

La validación biológica complementaria debe iniciar con medición química frecuente durante las primeras 48–72 h, especialmente en frío, para estrechar el intervalo de activación. Sin ese muestreo, la química sólo identifica un intervalo y no un instante exacto."""
        )
    )
    cells.append(nbformat.v4.new_markdown_cell("## Resumen reproducible"))
    cells.append(
        nbformat.v4.new_code_cell(
            """for row in result["validation"].itertuples(index=False):
    direction = f"{row.calibration_matrix} → {row.target_matrix}"
    print(
        f"{direction:22s} | {row.experiment_code:7s} | "
        f"RMSE={row.rmse_g_l_h:.3f} | r={row.correlation:.2f} | "
        f"onset error={row.onset_delay_h:.1f} h | censored={row.n_left_censored}"
    )

active = result["fit_parameters"].query("active_bound")
if not active.empty:
    print("\\nAdvertencia: parámetros en límite")
    display(active)

comparison = result["model_comparison"]
print("\\nDiagnóstico local de estimabilidad:")
display(result["jacobian_identifiability"])
print("\\nEstabilidad leave-one-batch-out de inicio/duración:")
display(result["loo_activation_summary"])
print("\\nCambio medio de |error de inicio|, nuevo - anterior [h]:",
      comparison["abs_onset_error_change_continuous_minus_threshold_release"].mean())
print("Cambio medio de RMSE, nuevo - anterior [g/L/h]:",
      comparison["rmse_change_continuous_minus_threshold_release"].mean())
print("\\nCambio medio de |error de duración visual 2-10%|, nuevo - anterior [h]:",
      comparison["abs_visual_rise_duration_error_change_continuous_minus_threshold_release"].mean())
print("\\nCambio medio de |error del máximo post-pulso|, nuevo - anterior [h]:",
      comparison["abs_postpulse_peak_error_change_continuous_minus_threshold_release"].mean())
print("\\nArtefactos guardados en:", analysis.RESULTS_DIR.relative_to(ROOT))"""
        )
    )
    cells.append(
        nbformat.v4.new_markdown_cell(
            """## Artefactos

La ejecución guarda offsets por sensor/corrida, señal antes y después de corregir cero, observaciones filtradas, soporte químico, pulsos, temperatura, parámetros y multistart, perfiles de objetivo, Jacobiano, leave-one-batch-out, predicciones, métricas, validaciones, figuras y manifiesto JSON. El notebook ejecutado conserva todas las tablas y gráficos visibles."""
        )
    )
    nb["cells"] = cells
    nbformat.write(nb, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 1800) -> None:
    runtime_dir = RESULTS_DIR / ".jupyter_runtime"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime_dir)
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
    parser.add_argument("--create-only", action="store_true", help="Create the notebook without executing it.")
    parser.add_argument("--timeout", type=int, default=1800)
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
