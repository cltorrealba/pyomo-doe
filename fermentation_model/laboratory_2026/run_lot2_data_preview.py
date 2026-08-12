from __future__ import annotations

import argparse
import json
import os
import re
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

from shared.paths import DATA_DIR, LABORATORY_2026_RESULTS_DIR


LOT2_DIR = DATA_DIR / "Laboratorio 2026" / "MBDoE_2026" / "DOE_Lote_2"
Y15_DIR = LOT2_DIR / "Datos Y15 lote 2"
ETHANOL_WORKBOOK = LOT2_DIR / "DOE_Lote_2_ethanol_manual_entry.xlsx"
SCHEDULE_PATH = (
    SCRIPT_DIR
    / "archive"
    / "design_execution_copy_2026-06-09"
    / "rec"
    / "lot_2_schedule.csv"
)
RESULTS_DIR = LABORATORY_2026_RESULTS_DIR / "lot2_data_preview"
PROCESSED_DIR = RESULTS_DIR / "processed"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "fermentation_lot2_data_preview.ipynb"
EXECUTED_NOTEBOOK_PATH = NOTEBOOK_DIR / "fermentation_lot2_data_preview.executed.ipynb"

PROCESS_CONFIG = {
    "F1": {
        "lab": "LAB004",
        "candidate": "synthetic_glucose_rich_fructose_pulse",
        "folder": "DOE26-F04_F1_Lot2_synthetic_glucose_rich_fructose_pulse",
    },
    "F2": {
        "lab": "LAB005",
        "candidate": "synthetic_lit_SM410_24C_highN_strip",
        "folder": "DOE26-F05_F2_Lot2_synthetic_lit_SM410_24C_highN_strip",
    },
    "F3": {
        "lab": "LAB006",
        "candidate": "synthetic_viable_biomass_step",
        "folder": "DOE26-F06_F3_Lot2_synthetic_viable_biomass_step",
    },
}
LAB_TO_PROCESS = {config["lab"]: process for process, config in PROCESS_CONFIG.items()}

ANALYTE_MAP = {
    "GLUCOSE-320": "glucose_g_l",
    "GLU-FRU-320": "sugar_total_g_l",
    "AMMONIA": "ammonia_mg_l",
    "PAN": "pan_mg_l",
    "GLYCEROL": "glycerol_g_l",
    "PYRUVIC ACID": "pyruvic_acid_mg_l",
    "ACETALDEHID-1000": "acetaldehyde_mg_l",
    "ACETIC ACID": "acetic_acid_g_l",
}

BIOMASS_G_PER_CELL = 30e-12
MILLION_CELLS_ML_TO_G_L = 1e6 * 1000.0 * BIOMASS_G_PER_CELL
PERCENT_VV_TO_G_L = 7.89
INITIAL_VOLUME_ML = 2000.0
SMALL_SAMPLE_ML = 5.0
ETHANOL_SAMPLE_ML = 50.0
STANDARD_MOLAR_VOLUME_L_MOL = 22.414
CO2_MOLAR_MASS_G_MOL = 44.01
MAX_CO2_INTEGRATION_GAP_H = 0.5


def safe_numeric(values: pd.Series) -> pd.Series:
    if values.dtype == object:
        values = values.astype(str).str.strip().str.replace(",", ".", regex=False)
    return pd.to_numeric(values, errors="coerce")


def parse_description_datetime(values: pd.Series) -> pd.Series:
    cleaned = values.astype(str).str.strip().str.replace(r"[^0-9: /-]+$", "", regex=True)
    cleaned = cleaned.str.rstrip("-").str.strip()
    return pd.to_datetime(cleaned, format="mixed", dayfirst=True, errors="coerce")


def load_oculyze() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    repair_rows: list[dict[str, object]] = []
    for process, config in PROCESS_CONFIG.items():
        path = LOT2_DIR / config["folder"] / f"Oculyze_{config['lab']}.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        frame["source_file"] = str(path)
        frame["process"] = process
        frame["sample_id_raw"] = frame["Name"].astype(str).str.strip()
        frame["sample_sequence_raw"] = safe_numeric(
            frame["sample_id_raw"].str.extract(r"-(\d+)\s*$", expand=False)
        )
        frame["sample_collection_datetime_raw"] = parse_description_datetime(frame["Description"])
        frame["analysis_datetime"] = pd.to_datetime(
            frame["Date"].astype(str).str.strip() + " " + frame["Time"].astype(str).str.strip(),
            format="mixed",
            errors="coerce",
        )
        frame = frame.sort_values(["sample_collection_datetime_raw", "analysis_datetime"]).reset_index(drop=True)

        if process in {"F1", "F2"}:
            expected = [1] + list(range(3, 19))
        else:
            expected = [1] + list(range(3, 17))
        if len(frame) != len(expected):
            raise RuntimeError(
                f"Unexpected Oculyze row count for {process}: {len(frame)} rows; expected {len(expected)}"
            )
        frame["sample_sequence"] = expected
        frame["sample_id"] = [f"DOE-{config['lab']}-{sequence}" for sequence in expected]
        frame["sample_id_repaired"] = frame["sample_id_raw"].ne(frame["sample_id"])

        duplicate_time = frame["sample_collection_datetime_raw"].duplicated(keep=False)
        frame["sample_collection_datetime"] = frame["sample_collection_datetime_raw"]
        frame["time_source"] = "oculyze_description"
        for idx in frame.index[duplicate_time]:
            previous = frame.loc[: idx - 1, "sample_collection_datetime"].dropna()
            if previous.empty:
                continue
            same_time_before = (
                frame.loc[: idx - 1, "sample_collection_datetime_raw"]
                .eq(frame.loc[idx, "sample_collection_datetime_raw"])
                .any()
            )
            if same_time_before:
                inferred = frame.loc[idx, "analysis_datetime"].floor("30min") - pd.Timedelta(hours=1)
                if inferred <= previous.iloc[-1]:
                    inferred = previous.iloc[-1] + pd.Timedelta(hours=4)
                frame.loc[idx, "sample_collection_datetime"] = inferred
                frame.loc[idx, "time_source"] = "duplicate_description_repaired_from_analysis_order"

        for idx in frame.index[frame["sample_id_repaired"] | duplicate_time]:
            reasons = []
            if bool(frame.loc[idx, "sample_id_repaired"]):
                reasons.append("sequence_id_repair")
            if bool(duplicate_time.loc[idx]):
                reasons.append(
                    "duplicate_description_time_repaired"
                    if frame.loc[idx, "sample_collection_datetime"]
                    != frame.loc[idx, "sample_collection_datetime_raw"]
                    else "duplicate_time_context"
                )
            repair_rows.append(
                {
                    "process": process,
                    "sample_id_raw": frame.loc[idx, "sample_id_raw"],
                    "sample_id": frame.loc[idx, "sample_id"],
                    "raw_datetime": frame.loc[idx, "sample_collection_datetime_raw"],
                    "used_datetime": frame.loc[idx, "sample_collection_datetime"],
                    "repair_reason": ";".join(reasons),
                    "record_changed": bool(frame.loc[idx, "sample_id_repaired"])
                    or frame.loc[idx, "sample_collection_datetime"]
                    != frame.loc[idx, "sample_collection_datetime_raw"],
                }
            )

        for column in ("Concentration", "Viability", "Budding Index", "Density", "Temperature"):
            frame[column] = safe_numeric(frame[column])
        frame["x_total_g_l"] = frame["Concentration"] * MILLION_CELLS_ML_TO_G_L
        frame["x_viable_g_l"] = frame["x_total_g_l"] * frame["Viability"] / 100.0
        frame["x_dead_g_l"] = (frame["x_total_g_l"] - frame["x_viable_g_l"]).clip(lower=0.0)
        frame["oculyze_concentration_million_cells_ml"] = frame["Concentration"]
        frame["viability_percent"] = frame["Viability"]
        frame["budding_index_percent"] = frame["Budding Index"]
        frame["density_oculyze"] = frame["Density"]
        frame["temperature_oculyze_c"] = frame["Temperature"]
        frames.append(frame)

    if not frames:
        raise FileNotFoundError("No Oculyze reports were found for Lot 2")
    oculyze = pd.concat(frames, ignore_index=True)
    t0_by_process = (
        oculyze.sort_values(["process", "sample_sequence"])
        .groupby("process")["sample_collection_datetime"]
        .first()
    )
    oculyze["t0"] = oculyze["process"].map(t0_by_process)
    oculyze["t_h"] = (
        oculyze["sample_collection_datetime"] - oculyze["t0"]
    ).dt.total_seconds() / 3600.0
    keep = [
        "process",
        "sample_id",
        "sample_id_raw",
        "sample_sequence",
        "sample_id_repaired",
        "sample_collection_datetime_raw",
        "sample_collection_datetime",
        "analysis_datetime",
        "time_source",
        "t0",
        "t_h",
        "oculyze_concentration_million_cells_ml",
        "viability_percent",
        "budding_index_percent",
        "density_oculyze",
        "temperature_oculyze_c",
        "x_total_g_l",
        "x_viable_g_l",
        "x_dead_g_l",
        "source_file",
    ]
    return oculyze[keep].sort_values(["process", "sample_sequence"]), pd.DataFrame(repair_rows)


def load_y15() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for path in sorted(Y15_DIR.glob("*.txt")):
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["sample_id", "analyte", "replicate", "value_raw", "unit", "analysis_datetime"],
            dtype=str,
        )
        frame["source_file"] = path.name
        frames.append(frame)
    if not frames:
        raise FileNotFoundError(f"No Y15 text files found under {Y15_DIR}")
    raw = pd.concat(frames, ignore_index=True)
    for column in ("sample_id", "analyte", "replicate", "unit"):
        raw[column] = raw[column].astype(str).str.strip()
    raw["value"] = safe_numeric(raw["value_raw"])
    raw["analysis_datetime"] = pd.to_datetime(raw["analysis_datetime"], dayfirst=True, errors="coerce")
    raw["lab"] = raw["sample_id"].str.extract(r"DOE-(LAB\d+)-", expand=False)
    raw["process"] = raw["lab"].map(LAB_TO_PROCESS)
    raw["sample_sequence"] = safe_numeric(raw["sample_id"].str.extract(r"-(\d+)\s*$", expand=False))
    raw["mapped_column"] = raw["analyte"].map(ANALYTE_MAP)

    duplicate_key = ["sample_id", "analyte", "replicate", "value", "unit", "analysis_datetime"]
    exact_duplicate = raw.duplicated(duplicate_key, keep=False)
    canonical = raw.drop_duplicates(duplicate_key).copy()
    conflict = (
        canonical.groupby(["sample_id", "analyte"])["value"]
        .nunique(dropna=True)
        .rename("n_distinct_values")
        .reset_index()
    )
    conflict = conflict[conflict["n_distinct_values"] > 1].copy()
    canonical = canonical.sort_values("analysis_datetime").drop_duplicates(
        ["sample_id", "analyte"], keep="last"
    )

    mapped = canonical[canonical["mapped_column"].notna()].copy()
    wide = mapped.pivot(index=["process", "sample_id", "sample_sequence"], columns="mapped_column", values="value").reset_index()
    timing = (
        canonical.groupby(["process", "sample_id", "sample_sequence"], dropna=False)
        .agg(
            y15_analysis_datetime_start=("analysis_datetime", "min"),
            y15_analysis_datetime_end=("analysis_datetime", "max"),
            y15_source_files=("source_file", lambda values: "; ".join(sorted(set(values)))),
        )
        .reset_index()
    )
    wide = wide.merge(timing, on=["process", "sample_id", "sample_sequence"], how="left")
    for column in ANALYTE_MAP.values():
        if column not in wide:
            wide[column] = np.nan
    wide["fructose_g_l_raw"] = wide["sugar_total_g_l"] - wide["glucose_g_l"]
    wide["glucose_for_model_g_l"] = wide["glucose_g_l"].clip(lower=0.0)
    wide["fructose_for_model_g_l"] = wide["fructose_g_l_raw"].clip(lower=0.0)
    wide["ammonia_for_model_mg_l"] = wide["ammonia_mg_l"].clip(lower=0.0)
    wide["pan_for_model_mg_l"] = wide["pan_mg_l"].clip(lower=0.0)
    wide["yan_calc_mg_l_raw"] = wide[["ammonia_mg_l", "pan_mg_l"]].sum(axis=1, min_count=2)
    wide["yan_for_model_physical_mg_l"] = wide[
        ["ammonia_for_model_mg_l", "pan_for_model_mg_l"]
    ].sum(axis=1, min_count=2)
    wide["negative_censored_for_model"] = (
        wide[["glucose_g_l", "fructose_g_l_raw", "ammonia_mg_l", "pan_mg_l", "pyruvic_acid_mg_l"]]
        .lt(0.0)
        .any(axis=1)
    )
    duplicate_summary = pd.DataFrame(
        {
            "metric": ["raw_rows", "exact_duplicate_rows", "canonical_sample_analyte_rows", "conflicting_sample_analytes"],
            "value": [len(raw), int(exact_duplicate.sum()), len(canonical), len(conflict)],
        }
    )
    return wide.sort_values(["process", "sample_sequence"]), duplicate_summary, conflict


def load_ethanol() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for process in PROCESS_CONFIG:
        frame = pd.read_excel(ETHANOL_WORKBOOK, sheet_name=process, dtype=object)
        frame["source_sheet"] = process
        frames.append(frame)
    ethanol = pd.concat(frames, ignore_index=True)
    ethanol["process"] = ethanol["process"].astype(str).str.strip()
    ethanol["sample_id"] = ethanol["sample_id"].astype(str).str.strip()
    ethanol["sample_sequence"] = safe_numeric(ethanol["sample_sequence"])
    for column in ("ethanol_percent_vv_rep1", "ethanol_percent_vv_rep2", "ethanol_g_l_manual"):
        if column not in ethanol:
            ethanol[column] = np.nan
        ethanol[column + "_num"] = safe_numeric(ethanol[column])
    reps = ethanol[["ethanol_percent_vv_rep1_num", "ethanol_percent_vv_rep2_num"]]
    ethanol["ethanol_percent_vv_final_calc"] = reps.mean(axis=1, skipna=True)
    ethanol.loc[reps.notna().sum(axis=1).eq(0), "ethanol_percent_vv_final_calc"] = np.nan
    ethanol["ethanol_g_l_from_percent_vv"] = ethanol["ethanol_percent_vv_final_calc"] * PERCENT_VV_TO_G_L
    ethanol["ethanol_g_l_for_model"] = ethanol["ethanol_g_l_manual_num"].combine_first(
        ethanol["ethanol_g_l_from_percent_vv"]
    )
    ethanol["ethanol_measured"] = ethanol["ethanol_g_l_for_model"].notna()
    ethanol["sample_volume_ml"] = np.where(
        ethanol["ethanol_measured"], ETHANOL_SAMPLE_ML, SMALL_SAMPLE_ML
    )
    ethanol["ethanol_measurement_datetime"] = pd.to_datetime(
        ethanol.get("ethanol_measurement_datetime"), errors="coerce", dayfirst=True
    )
    keep = [
        "process",
        "sample_id",
        "sample_sequence",
        "ethanol_percent_vv_rep1_num",
        "ethanol_percent_vv_rep2_num",
        "ethanol_percent_vv_final_calc",
        "ethanol_g_l_manual_num",
        "ethanol_g_l_from_percent_vv",
        "ethanol_g_l_for_model",
        "ethanol_measured",
        "sample_volume_ml",
        "ethanol_measurement_datetime",
        "ethanol_method",
        "operator",
        "qc_flag",
        "notes",
    ]
    return ethanol[[column for column in keep if column in ethanol]].sort_values(
        ["process", "sample_sequence"]
    )


def build_sample_table(
    y15: pd.DataFrame, ethanol: pd.DataFrame, oculyze: pd.DataFrame
) -> pd.DataFrame:
    sample = ethanol.merge(
        y15,
        on=["process", "sample_id", "sample_sequence"],
        how="left",
        validate="one_to_one",
    )
    oculyze_merge = oculyze.drop(columns=["source_file"], errors="ignore")
    sample = sample.merge(
        oculyze_merge,
        on=["process", "sample_id", "sample_sequence"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_oculyze"),
    )
    t0_by_process = (
        oculyze.sort_values(["process", "sample_sequence"])
        .groupby("process")["sample_collection_datetime"]
        .first()
    )
    sample["t0"] = sample["process"].map(t0_by_process)

    sample["sample_datetime_used"] = sample["sample_collection_datetime"]
    sample["sample_time_source"] = sample["time_source"].fillna("")
    for process, group in sample.groupby("process"):
        known = group.dropna(subset=["sample_collection_datetime", "y15_analysis_datetime_start"])
        lags = (
            known["y15_analysis_datetime_start"] - known["sample_collection_datetime"]
        ).dt.total_seconds() / 60.0
        plausible = lags[(lags >= 0.0) & (lags <= 240.0)]
        median_lag_min = float(plausible.median()) if not plausible.empty else 30.0
        missing_idx = group.index[
            group["sample_datetime_used"].isna() & group["y15_analysis_datetime_start"].notna()
        ]
        sample.loc[missing_idx, "sample_datetime_used"] = (
            sample.loc[missing_idx, "y15_analysis_datetime_start"]
            - pd.to_timedelta(median_lag_min, unit="min")
        )
        sample.loc[missing_idx, "sample_time_source"] = (
            "y15_analysis_start_minus_process_median_lag"
        )
    sample["t_h"] = (
        sample["sample_datetime_used"] - sample["t0"]
    ).dt.total_seconds() / 3600.0
    sample["sample_observed"] = (
        sample["ethanol_measured"].fillna(False).astype(bool)
        | sample["y15_analysis_datetime_start"].notna()
        | sample["x_total_g_l"].notna()
    )
    sample["time_inferred"] = sample["sample_observed"] & ~sample[
        "sample_time_source"
    ].eq("oculyze_description")
    sample["sample_volume_ml"] = np.select(
        [sample["ethanol_measured"].fillna(False).astype(bool), sample["sample_observed"]],
        [ETHANOL_SAMPLE_ML, SMALL_SAMPLE_ML],
        default=0.0,
    )

    sample = sample.sort_values(["process", "t_h", "sample_sequence"]).reset_index(drop=True)
    sample["cumulative_sample_volume_ml"] = sample.groupby("process")["sample_volume_ml"].cumsum()
    sample["reactor_volume_after_sample_ml"] = INITIAL_VOLUME_ML - sample["cumulative_sample_volume_ml"]

    sample["ethanol_increment_g_l"] = sample.groupby("process")["ethanol_g_l_for_model"].diff()
    sample["known_sugar_input_before_sample_g_l"] = np.where(
        sample["process"].eq("F1") & sample["t_h"].ge(70.0), 35.0, 0.0
    )
    sample["sugar_consumed_since_first_g_l"] = (
        sample.groupby("process")["sugar_total_g_l"].transform("first")
        + sample["known_sugar_input_before_sample_g_l"]
        - sample["sugar_total_g_l"]
    )
    sample["ethanol_from_sugar_upper_g_l"] = 0.55 * sample["sugar_consumed_since_first_g_l"].clip(lower=0.0)
    sample["ethanol_mass_balance_screen"] = np.where(
        sample["ethanol_measured"]
        & sample["ethanol_g_l_for_model"].gt(sample["ethanol_from_sugar_upper_g_l"] + 15.0),
        "review_time_mapping_or_measurement",
        "pass_or_not_tested",
    )
    sample["ethanol_use_for_model"] = sample["ethanol_measured"] & sample[
        "ethanol_mass_balance_screen"
    ].eq("pass_or_not_tested")
    return sample


def load_temperature_and_co2(
    t0_by_process: dict[str, pd.Timestamp], sample: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    temperature_frames: list[pd.DataFrame] = []
    co2_frames: list[pd.DataFrame] = []
    inventory_rows: list[dict[str, object]] = []
    for process, config in PROCESS_CONFIG.items():
        folder = LOT2_DIR / config["folder"]
        t0 = pd.Timestamp(t0_by_process[process])
        temp_parts = []
        for path in sorted(folder.glob(f"Temp_{process}*.csv")):
            try:
                frame = pd.read_csv(path)
            except (FileNotFoundError, OSError) as exc:
                inventory_rows.append({"process": process, "kind": "temperature", "path": str(path), "status": f"unreadable:{exc}"})
                continue
            frame["source_file"] = str(path)
            temp_parts.append(frame)
            inventory_rows.append({"process": process, "kind": "temperature", "path": str(path), "status": "loaded", "n_rows": len(frame)})
        if temp_parts:
            temp = pd.concat(temp_parts, ignore_index=True)
            temp["timestamp"] = pd.to_datetime(temp["timestamp"], errors="coerce")
            for column in ("T", "SP", "nutricion_activa"):
                if column in temp:
                    temp[column] = safe_numeric(temp[column])
            temp = temp.dropna(subset=["timestamp", "T"]).sort_values("timestamp").drop_duplicates("timestamp")
            temp["process"] = process
            temp["t_h"] = (temp["timestamp"] - t0).dt.total_seconds() / 3600.0
            aggregation = {column: "median" for column in ("T", "SP") if column in temp}
            if "nutricion_activa" in temp:
                aggregation["nutricion_activa"] = "max"
            sampled = (
                temp.set_index("timestamp")
                .resample("10min")
                .agg(aggregation)
                .dropna(subset=["T"])
                .reset_index()
            )
            sampled["process"] = process
            sampled["t_h"] = (sampled["timestamp"] - t0).dt.total_seconds() / 3600.0
            temperature_frames.append(sampled)

        co2_parts = []
        filtered_paths = sorted(folder.glob(f"CO2_FILT_{process}*.csv"))
        raw_paths = sorted(folder.glob(f"CO2_{process}*.csv"))
        loaded_periods: list[tuple[pd.Timestamp, pd.Timestamp]] = []
        for path in filtered_paths:
            try:
                frame = pd.read_csv(path)
            except (FileNotFoundError, OSError) as exc:
                inventory_rows.append({"process": process, "kind": "co2_filtered", "path": str(path), "status": f"unreadable:{exc}"})
                continue
            frame["timestamp"] = pd.to_datetime(frame["timestamp"], errors="coerce")
            frame["flow_filt_sccm"] = safe_numeric(frame["flow_filt_sccm"])
            frame = frame.dropna(subset=["timestamp", "flow_filt_sccm"])
            frame["co2_source"] = "instrument_filtered"
            frame["source_file"] = str(path)
            co2_parts.append(frame)
            if not frame.empty:
                loaded_periods.append((frame["timestamp"].min(), frame["timestamp"].max()))
            inventory_rows.append({"process": process, "kind": "co2_filtered", "path": str(path), "status": "loaded", "n_rows": len(frame)})
        for path in raw_paths:
            try:
                raw = pd.read_csv(path)
            except (FileNotFoundError, OSError) as exc:
                inventory_rows.append({"process": process, "kind": "co2_raw", "path": str(path), "status": f"unreadable:{exc}"})
                continue
            raw["timestamp"] = pd.to_datetime(raw["timestamp"], errors="coerce")
            raw_column = next(
                (
                    column
                    for column in ("flow_raw_sccm", "flow_sccm", "flow_filt_sccm")
                    if column in raw
                ),
                None,
            )
            if raw_column not in raw:
                continue
            raw[raw_column] = safe_numeric(raw[raw_column])
            raw = raw.dropna(subset=["timestamp", raw_column]).sort_values("timestamp")
            if loaded_periods:
                covered = np.zeros(len(raw), dtype=bool)
                for start, end in loaded_periods:
                    covered |= raw["timestamp"].between(start, end).to_numpy()
                raw = raw.loc[~covered].copy()
            if raw.empty:
                continue
            raw["flow_filt_sccm"] = raw[raw_column].rolling(61, center=True, min_periods=5).median()
            raw["flow_filt_sccm"] = raw["flow_filt_sccm"].interpolate(limit_direction="both")
            raw["co2_source"] = "raw_rolling_median_fallback"
            raw["source_file"] = str(path)
            co2_parts.append(raw)
            inventory_rows.append({"process": process, "kind": "co2_raw_fallback", "path": str(path), "status": "loaded_uncovered_rows", "n_rows": len(raw)})
        if co2_parts:
            co2 = pd.concat(co2_parts, ignore_index=True)
            co2 = co2.sort_values(["timestamp", "co2_source"]).drop_duplicates("timestamp", keep="first")
            co2["process"] = process
            co2["t_h"] = (co2["timestamp"] - t0).dt.total_seconds() / 3600.0
            sampled = (
                co2.set_index("timestamp")
                .resample("10min")
                .agg(
                    flow_filt_sccm=("flow_filt_sccm", "median"),
                    co2_source=("co2_source", lambda values: ";".join(sorted(set(values)))),
                )
                .dropna(subset=["flow_filt_sccm"])
                .reset_index()
            )
            sampled["process"] = process
            sampled["t_h"] = (sampled["timestamp"] - t0).dt.total_seconds() / 3600.0
            process_samples = sample[sample["process"].eq(process)].sort_values("t_h")
            times = process_samples["t_h"].to_numpy(dtype=float)
            volumes = process_samples["reactor_volume_after_sample_ml"].to_numpy(dtype=float)
            positions = np.searchsorted(times, sampled["t_h"].to_numpy(dtype=float), side="right") - 1
            sampled["reactor_volume_ml_sampling_only"] = np.where(
                positions >= 0, volumes[np.maximum(positions, 0)], INITIAL_VOLUME_ML
            )
            sampled["flow_filt_sccm_physical"] = sampled["flow_filt_sccm"].clip(lower=0.0)
            sampled["flow_filt_sccm_v0_corrected"] = (
                sampled["flow_filt_sccm_physical"]
                * INITIAL_VOLUME_ML
                / sampled["reactor_volume_ml_sampling_only"]
            )
            sampled["co2_rate_equivalent_g_l_h"] = (
                sampled["flow_filt_sccm_v0_corrected"]
                * 60.0
                / 1000.0
                * CO2_MOLAR_MASS_G_MOL
                / STANDARD_MOLAR_VOLUME_L_MOL
                / (INITIAL_VOLUME_ML / 1000.0)
            )
            dt = sampled["t_h"].diff()
            increment = 0.5 * (
                sampled["co2_rate_equivalent_g_l_h"]
                + sampled["co2_rate_equivalent_g_l_h"].shift()
            ) * dt
            sampled["co2_increment_g_l"] = increment.where(
                dt.gt(0.0) & dt.le(MAX_CO2_INTEGRATION_GAP_H), 0.0
            ).fillna(0.0)
            sampled["co2_cumulative_observed_lower_bound_g_l"] = sampled["co2_increment_g_l"].cumsum()
            co2_frames.append(sampled)
    temperature = pd.concat(temperature_frames, ignore_index=True) if temperature_frames else pd.DataFrame()
    co2 = pd.concat(co2_frames, ignore_index=True) if co2_frames else pd.DataFrame()
    return temperature, co2, pd.DataFrame(inventory_rows)


def build_operational_inputs(
    sample: pd.DataFrame, temperature: pd.DataFrame
) -> pd.DataFrame:
    schedule = pd.read_csv(SCHEDULE_PATH)
    schedule = schedule[schedule["candidate"].isin([config["candidate"] for config in PROCESS_CONFIG.values()])]
    schedule = schedule[schedule["event"].eq("manual_pulse")].copy()
    candidate_to_process = {config["candidate"]: process for process, config in PROCESS_CONFIG.items()}
    schedule["process"] = schedule["candidate"].map(candidate_to_process)
    schedule["nominal_relative_time_h"] = safe_numeric(schedule["relative_time_h"])
    schedule["amount"] = safe_numeric(schedule["amount"])
    schedule["volume_ml"] = safe_numeric(schedule["volume_ml"])
    schedule["actual_or_model_time_h"] = schedule["nominal_relative_time_h"]
    schedule["evidence"] = "protocol_nominal_unverified"

    f1_active = temperature[
        temperature["process"].eq("F1")
        & pd.to_numeric(temperature.get("nutricion_activa"), errors="coerce").gt(0.0)
    ]
    f1_f = schedule["process"].eq("F1") & schedule["channel"].eq("F")
    if not f1_active.empty and f1_f.any():
        detected_time = float(f1_active["t_h"].iloc[0])
        schedule.loc[f1_f, "actual_or_model_time_h"] = detected_time
        schedule.loc[f1_f, "evidence"] = "controller_flag_plus_observed_fructose_jump"

    f3_x = schedule["process"].eq("F3") & schedule["channel"].eq("X")
    if f3_x.any():
        schedule.loc[f3_x, "evidence"] = "protocol_time_plus_observed_biomass_step"
    schedule["datetime_reanchored"] = schedule.apply(
        lambda row: sample.loc[sample["process"].eq(row["process"]), "t0"].iloc[0]
        + pd.to_timedelta(row["actual_or_model_time_h"], unit="h"),
        axis=1,
    )
    return schedule[
        [
            "process",
            "candidate",
            "channel",
            "amount",
            "volume_ml",
            "nominal_relative_time_h",
            "actual_or_model_time_h",
            "datetime_reanchored",
            "evidence",
            "details",
        ]
    ].sort_values(["process", "actual_or_model_time_h"])


def write_qc_tables(
    sample: pd.DataFrame,
    oculyze_repairs: pd.DataFrame,
    y15_duplicate_summary: pd.DataFrame,
    y15_conflicts: pd.DataFrame,
    temperature: pd.DataFrame,
    co2: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for process in PROCESS_CONFIG:
        s = sample[sample["process"].eq(process)]
        temp = temperature[temperature["process"].eq(process)] if not temperature.empty else pd.DataFrame()
        gas = co2[co2["process"].eq(process)] if not co2.empty else pd.DataFrame()
        rows.append(
            {
                "process": process,
                "n_sample_rows": len(s),
                "n_samples_with_observation_evidence": int(s["sample_observed"].sum()),
                "n_y15_samples": int(s["sugar_total_g_l"].notna().sum()),
                "n_ethanol": int(s["ethanol_measured"].sum()),
                "n_ethanol_used_after_qc": int(s["ethanol_use_for_model"].sum()),
                "n_biomass": int(s["x_viable_g_l"].notna().sum()),
                "n_time_inferred": int(s["time_inferred"].sum()),
                "n_negative_censored": int(s["negative_censored_for_model"].eq(True).sum()),
                "final_sampling_only_volume_ml": float(s["reactor_volume_after_sample_ml"].min()),
                "temperature_first_h": float(temp["t_h"].min()) if not temp.empty else np.nan,
                "temperature_last_h": float(temp["t_h"].max()) if not temp.empty else np.nan,
                "co2_first_h": float(gas["t_h"].min()) if not gas.empty else np.nan,
                "co2_last_h": float(gas["t_h"].max()) if not gas.empty else np.nan,
                "chemistry_last_h": float(s.loc[s["sugar_total_g_l"].notna(), "t_h"].max()),
            }
        )
    qc = pd.DataFrame(rows)
    qc.to_csv(PROCESSED_DIR / "lot2_data_qc_summary.csv", index=False)
    oculyze_repairs.to_csv(PROCESSED_DIR / "lot2_oculyze_repairs.csv", index=False)
    y15_duplicate_summary.to_csv(PROCESSED_DIR / "lot2_y15_duplicate_summary.csv", index=False)
    y15_conflicts.to_csv(PROCESSED_DIR / "lot2_y15_conflicts.csv", index=False)
    return qc


def add_pulse_lines(ax: plt.Axes, inputs: pd.DataFrame, process: str, channel: str) -> None:
    selected = inputs[inputs["process"].eq(process) & inputs["channel"].eq(channel)]
    for row in selected.itertuples(index=False):
        ax.axvline(float(row.actual_or_model_time_h), color="black", linestyle="--", linewidth=1)
        ax.annotate(
            f"{channel} +{float(row.amount):g}",
            xy=(float(row.actual_or_model_time_h), 0.98),
            xycoords=("data", "axes fraction"),
            xytext=(3, -3),
            textcoords="offset points",
            rotation=90,
            va="top",
            fontsize=8,
        )


def plot_process_panels(
    sample: pd.DataFrame,
    temperature: pd.DataFrame,
    co2: pd.DataFrame,
    inputs: pd.DataFrame,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    colors = {
        "G": "#0072B2",
        "F": "#D55E00",
        "total": "#4D4D4D",
        "N": "#009E73",
        "E": "#CC79A7",
        "Gly": "#E69F00",
        "X": "#56B4E9",
        "Xd": "#8C564B",
    }
    for process, config in PROCESS_CONFIG.items():
        s = sample[sample["process"].eq(process)].sort_values("t_h")
        temp = temperature[temperature["process"].eq(process)].sort_values("t_h")
        gas = co2[co2["process"].eq(process)].sort_values("t_h")
        fig, axes = plt.subplots(4, 2, figsize=(14, 16), sharex=False)
        ax = axes[0, 0]
        ax.plot(temp["t_h"], temp["T"], label="measured T", color="#D55E00")
        if "SP" in temp:
            ax.step(temp["t_h"], temp["SP"], where="post", label="setpoint", color="#0072B2")
        ax.set_ylabel("Temperature [degC]")
        ax.legend()

        ax = axes[0, 1]
        ax.plot(gas["t_h"], gas["flow_filt_sccm"], alpha=0.45, label="reported CO2_FILT")
        ax.plot(gas["t_h"], gas["flow_filt_sccm_v0_corrected"], label="volume-corrected", color="#009E73")
        ax.set_ylabel("CO2 flow [sccm]")
        ax.legend()

        ax = axes[1, 0]
        ax.plot(s["t_h"], s["glucose_g_l"], "o-", label="glucose", color=colors["G"])
        ax.plot(s["t_h"], s["fructose_g_l_raw"], "s-", label="fructose", color=colors["F"])
        ax.plot(s["t_h"], s["sugar_total_g_l"], ".-", label="total sugar", color=colors["total"])
        density = ax.twinx()
        density.plot(s["t_h"], s["density_oculyze"], "x:", label="density", color="#9467BD")
        ax.set_ylabel("Sugar [g/L]")
        density.set_ylabel("Oculyze density")
        add_pulse_lines(ax, inputs, process, "G")
        add_pulse_lines(ax, inputs, process, "F")
        handles, labels = ax.get_legend_handles_labels()
        handles2, labels2 = density.get_legend_handles_labels()
        ax.legend(handles + handles2, labels + labels2, fontsize=8)

        ax = axes[1, 1]
        ax.plot(s["t_h"], s["ammonia_mg_l"], "o-", label="ammonia", color="#0072B2")
        ax.plot(s["t_h"], s["pan_mg_l"], "s-", label="PAN", color="#E69F00")
        ax.plot(s["t_h"], s["yan_for_model_physical_mg_l"], ".-", label="YAN model-safe", color=colors["N"])
        add_pulse_lines(ax, inputs, process, "N")
        ax.axhline(0.0, color="black", linewidth=0.7)
        ax.set_ylabel("Nitrogen [mg/L]")
        ax.legend(fontsize=8)

        ax = axes[2, 0]
        used = s["ethanol_use_for_model"].fillna(False)
        ax.scatter(s.loc[used, "t_h"], s.loc[used, "ethanol_g_l_for_model"], label="ethanol used", color=colors["E"])
        ax.scatter(s.loc[s["ethanol_measured"] & ~used, "t_h"], s.loc[s["ethanol_measured"] & ~used, "ethanol_g_l_for_model"], marker="x", s=70, label="ethanol QC review", color="#D62728")
        ax.plot(s["t_h"], s["glycerol_g_l"], "o-", label="glycerol", color=colors["Gly"])
        add_pulse_lines(ax, inputs, process, "E")
        ax.set_ylabel("Concentration [g/L]")
        ax.legend(fontsize=8)

        ax = axes[2, 1]
        ax.plot(s["t_h"], s["x_viable_g_l"], "o-", label="viable biomass", color=colors["X"])
        ax.plot(s["t_h"], s["x_dead_g_l"], "s-", label="dead biomass", color=colors["Xd"])
        budding = ax.twinx()
        budding.plot(s["t_h"], s["budding_index_percent"], ".:", label="budding index", color="#009E73")
        add_pulse_lines(ax, inputs, process, "X")
        ax.set_ylabel("Biomass [g/L]")
        budding.set_ylabel("Budding index [%]")
        handles, labels = ax.get_legend_handles_labels()
        handles2, labels2 = budding.get_legend_handles_labels()
        ax.legend(handles + handles2, labels + labels2, fontsize=8)

        ax = axes[3, 0]
        ax.plot(s["t_h"], s["pyruvic_acid_mg_l"], "o-", label="pyruvate [mg/L]")
        ax.plot(s["t_h"], s["acetaldehyde_mg_l"], "s-", label="acetaldehyde [mg/L]")
        acetate = ax.twinx()
        acetate.plot(s["t_h"], s["acetic_acid_g_l"], ".-", label="acetic acid [g/L]", color="#D55E00")
        ax.set_ylabel("Pyruvate / acetaldehyde [mg/L]")
        acetate.set_ylabel("Acetic acid [g/L]")
        handles, labels = ax.get_legend_handles_labels()
        handles2, labels2 = acetate.get_legend_handles_labels()
        ax.legend(handles + handles2, labels + labels2, fontsize=8)

        ax = axes[3, 1]
        ax.step(s["t_h"], s["reactor_volume_after_sample_ml"], where="post", color="#4D4D4D")
        ax.set_ylabel("Volume after sampling [mL]")
        ax.set_ylim(0.0, INITIAL_VOLUME_ML * 1.05)
        ax.grid(True, alpha=0.2)

        for axis in axes.ravel():
            axis.set_xlabel("Elapsed time from first Oculyze sample [h]")
            axis.grid(True, alpha=0.2)
        fig.suptitle(f"Lot 2 {process}: {config['candidate']}", fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.98))
        fig.savefig(FIGURE_DIR / f"lot2_{process}_process_preview.png", dpi=170)
        plt.close(fig)


def write_report(qc: pd.DataFrame, sample: pd.DataFrame, repairs: pd.DataFrame) -> None:
    review = sample[sample["ethanol_mass_balance_screen"].ne("pass_or_not_tested")][
        ["process", "sample_id", "t_h", "sugar_total_g_l", "ethanol_g_l_for_model", "ethanol_mass_balance_screen"]
    ]
    report = f"""# Lot 2 data ingestion and QC

## Scope

This artifact integrates Y15 chemistry, manual ethanol entries, Oculyze biomass, online temperature, and online CO2 for the three synthetic-must MBDoE fermentations in Lot 2.

## Main conventions

- Time zero is the first Oculyze sample in each reactor.
- Total Oculyze concentration is interpreted as million cells/mL. Biomass uses 30 pg/cell.
- Negative low-range analytical readings are preserved in raw columns and censored at zero only in explicitly named model-safe columns.
- Ethanol is recomputed from the entered replicates as `% v/v * 7.89 = g/L`; cached Excel formula cells are ignored.
- Sampling removes 50 mL when ethanol was measured and 5 mL for another sample with analytical evidence. Planned rows without any measurement evidence remove 0 mL. The CO2 flow is normalized back to the initial 2 L volume.
- Planned input events are not silently treated as executed. The F1 fructose pulse and F3 biomass step have trajectory evidence; nitrogen pulse times remain protocol-based and unverified.

## QC summary

{qc.to_markdown(index=False)}

## Repaired Oculyze records

{repairs.to_markdown(index=False) if not repairs.empty else 'No repairs were required.'}

## Ethanol points held for review

The following points fail a deliberately conservative sugar-to-ethanol screening rule. They remain visible in plots and processed data but are excluded from the default calibration table until sample timing or measurement identity is confirmed.

{review.to_markdown(index=False) if not review.empty else 'No ethanol points were held for review.'}

## Interpretation boundary

The online process files end near 99 h, whereas chemistry continues to approximately 190 h. Cumulative CO2 after the last online record is therefore a lower bound, not a full-batch mass balance.
"""
    (RESULTS_DIR / "lot2_data_quality_report.md").write_text(report, encoding="utf-8")


def run_pipeline() -> dict[str, object]:
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    oculyze, repairs = load_oculyze()
    y15, duplicate_summary, conflicts = load_y15()
    ethanol = load_ethanol()
    sample = build_sample_table(y15, ethanol, oculyze)
    t0_by_process = {
        process: pd.Timestamp(group["t0"].iloc[0])
        for process, group in sample.groupby("process")
    }
    temperature, co2, process_inventory = load_temperature_and_co2(t0_by_process, sample)
    inputs = build_operational_inputs(sample, temperature)
    qc = write_qc_tables(sample, repairs, duplicate_summary, conflicts, temperature, co2)

    oculyze.to_csv(PROCESSED_DIR / "lot2_oculyze_processed.csv", index=False)
    y15.to_csv(PROCESSED_DIR / "lot2_y15_wide_processed.csv", index=False)
    ethanol.to_csv(PROCESSED_DIR / "lot2_ethanol_manual_processed.csv", index=False)
    sample.to_csv(PROCESSED_DIR / "lot2_samples_integrated.csv", index=False)
    temperature.to_csv(PROCESSED_DIR / "lot2_temperature_10min.csv", index=False)
    co2.to_csv(PROCESSED_DIR / "lot2_co2_volume_corrected_10min.csv", index=False)
    inputs.to_csv(PROCESSED_DIR / "lot2_operational_inputs.csv", index=False)
    process_inventory.to_csv(PROCESSED_DIR / "lot2_process_file_inventory.csv", index=False)
    plot_process_panels(sample, temperature, co2, inputs)
    write_report(qc, sample, repairs)

    summary = {
        "n_processes": int(sample["process"].nunique()),
        "n_sample_rows": int(len(sample)),
        "n_samples_with_observation_evidence": int(sample["sample_observed"].sum()),
        "n_y15_samples": int(sample["sugar_total_g_l"].notna().sum()),
        "n_ethanol_measured": int(sample["ethanol_measured"].sum()),
        "n_ethanol_used_after_qc": int(sample["ethanol_use_for_model"].sum()),
        "n_biomass_observations": int(sample["x_viable_g_l"].notna().sum()),
        "n_oculyze_id_or_time_repairs": int(repairs["record_changed"].sum()),
        "n_oculyze_duplicate_context_rows": int(len(repairs)),
        "n_y15_conflicts": int(len(conflicts)),
        "time_zero_definition": "first Oculyze sample per process",
        "biomass_conversion_g_per_cell": BIOMASS_G_PER_CELL,
        "ethanol_conversion_g_l_per_percent_vv": PERCENT_VV_TO_G_L,
    }
    (RESULTS_DIR / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def create_notebook() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    cells = [
        nbformat.v4.new_markdown_cell(
            """# MBDoE 2026 Lot 2: integrated data preview and QC

## tl;dr

This notebook reconstructs the three synthetic-must Lot 2 fermentations from raw laboratory and process files. It preserves raw analytical values, creates separately named model-safe values, and makes every timing or identifier repair auditable. The executed outputs below determine which observations can enter the next calibration and estimability comparison."""
        ),
        nbformat.v4.new_markdown_cell(
            r"""## Context & Methods

### Key assumptions

1. The first Oculyze sample defines $t=0$ for each reactor.
2. Oculyze concentration is total cells in million cells/mL. Using 30 pg/cell,

$$X_{\mathrm{total}}[\mathrm{g/L}] = C[10^6\,\mathrm{cells/mL}]\times 0.03.$$

Viable and dead biomass are

$$X = X_{\mathrm{total}}\frac{V}{100},\qquad X_d=X_{\mathrm{total}}-X.$$

3. Ethanol is recomputed from manual entries using $E[\mathrm{g/L}]=7.89\,E[\%\,v/v]$.
4. A measured-ethanol sample removes 50 mL; another sample with analytical evidence removes 5 mL. A planned row without any measurement evidence removes 0 mL.
5. CO2 flow is normalized to the initial 2 L volume with

$$q_{\mathrm{CO_2},0}=q_{\mathrm{CO_2}}\frac{V_0}{V(t)}.$$

Protocol events and inferred times are labeled; they are not presented as direct observations."""
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

from laboratory_2026 import run_lot2_data_preview as analysis

summary = analysis.run_pipeline()
summary"""
        ),
        nbformat.v4.new_markdown_cell("## Data"),
        nbformat.v4.new_code_cell(
            """processed = analysis.PROCESSED_DIR
qc = pd.read_csv(processed / "lot2_data_qc_summary.csv")
samples = pd.read_csv(processed / "lot2_samples_integrated.csv")
inputs = pd.read_csv(processed / "lot2_operational_inputs.csv")
repairs = pd.read_csv(processed / "lot2_oculyze_repairs.csv")
display(qc)
display(inputs)
display(repairs)"""
        ),
        nbformat.v4.new_markdown_cell(
            """### Data-quality rules

- Raw negative analytical values remain available. Only columns ending in `for_model` are censored at zero.
- Oculyze identifiers are repaired only when the row sequence makes the typo unambiguous; the raw identifier remains beside the repaired one.
- Missing collection times are estimated from Y15 analysis time minus the process-specific median analytical lag.
- The ethanol mass-balance flag is a screening test, not an automatic data correction."""
        ),
        nbformat.v4.new_code_cell(
            """columns = [
    "process", "sample_id", "sample_sequence", "sample_datetime_used", "t_h",
    "sample_time_source", "glucose_g_l", "fructose_g_l_raw", "yan_for_model_physical_mg_l",
    "ethanol_g_l_for_model", "ethanol_use_for_model", "glycerol_g_l", "x_viable_g_l",
    "x_dead_g_l", "reactor_volume_after_sample_ml"
]
display(samples[columns].head(24))"""
        ),
        nbformat.v4.new_markdown_cell("## Results"),
        nbformat.v4.new_code_cell(
            """for process in ("F1", "F2", "F3"):
    display(Image(filename=str(analysis.FIGURE_DIR / f"lot2_{process}_process_preview.png")))"""
        ),
        nbformat.v4.new_markdown_cell("## Takeaways"),
        nbformat.v4.new_code_cell(
            """display(Markdown((analysis.RESULTS_DIR / "lot2_data_quality_report.md").read_text(encoding="utf-8")))"""
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


def execute_notebook(timeout: int = 3600) -> None:
    notebook = nbformat.read(NOTEBOOK_PATH, as_version=4)
    with tempfile.TemporaryDirectory(prefix="lot2_jupyter_runtime_") as runtime:
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
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    create_notebook()
    if args.no_execute:
        run_pipeline()
    else:
        execute_notebook(timeout=args.timeout)
    print(EXECUTED_NOTEBOOK_PATH if not args.no_execute else NOTEBOOK_PATH)


if __name__ == "__main__":
    main()
