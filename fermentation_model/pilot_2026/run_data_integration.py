from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
import zipfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE = Path(__file__).resolve().parent
REPO = WORKSPACE.parents[1]
CONFIG_PATH = WORKSPACE / "data_integration_config.json"
RUN_ID_RE = re.compile(r"(26\d{3})")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_number(value: object) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, str):
        value = value.strip().replace(",", ".")
        if not value or value.upper() == "NQ":
            return np.nan
    return float(value)


def normalize_datetime(value: object, run_id: str) -> pd.Timestamp:
    if pd.isna(value):
        return pd.NaT
    if isinstance(value, (pd.Timestamp, np.datetime64)) or hasattr(value, "year"):
        timestamp = pd.Timestamp(value)
    else:
        text_value = str(value).strip()
        # The master workbook contains strings such as "19-04-2026 15:00 PM".
        # The 24-hour clock is authoritative; the redundant AM/PM suffix is invalid.
        text_value = re.sub(
            r"(\s(?:1[3-9]|2[0-3]):\d{2}(?::\d{2})?)\s*(?:AM|PM)$",
            r"\1",
            text_value,
            flags=re.IGNORECASE,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            timestamp = pd.to_datetime(text_value, dayfirst=True, errors="coerce")
    if pd.isna(timestamp):
        return pd.NaT
    if run_id in {"26157", "26158", "26159"}:
        if timestamp.year == 2026 and timestamp.day == 4 and 6 <= timestamp.month <= 12:
            timestamp = timestamp.replace(month=4, day=timestamp.month)
    return timestamp


def canonical_primary_sample_id(value: object) -> str:
    sample_id = str(value).strip()
    short = re.fullmatch(r"(26\d{3})-(\d{2})", sample_id)
    if short:
        return f"{short.group(1)}-P-{short.group(2)}"
    return sample_id


def atomic_savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{path.stem}.", suffix=path.suffix, dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        fig.savefig(temporary, dpi=170, bbox_inches="tight")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_master_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pd.read_excel(path, sheet_name=sheet_name, header=3)


def build_primary_and_windows(
    master: Path,
    runs: list[str],
    run_metadata: dict[str, dict[str, object]],
    initial_volume_l: float,
    yeast: str,
    output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = read_master_sheet(master, "06_Resultados_primarios")
    primary["experiment_id"] = (
        primary["Código muestra"].astype(str).str.extract(r"^(26\d{3})", expand=False)
    )
    primary = primary[primary["experiment_id"].isin(runs)].copy()
    primary["timestamp_original"] = primary["Fecha-hora"].astype(str)
    primary["timestamp"] = [
        normalize_datetime(value, run)
        for value, run in zip(primary["Fecha-hora"], primary["experiment_id"])
    ]
    primary["date_normalization"] = np.where(
        primary["timestamp_original"] != primary["timestamp"].astype(str),
        "parsed_or_day_month_corrected",
        "unchanged",
    )
    ethanol_column = "Cf Alcolyzer Real (% v/v)"
    raw_ethanol_column = "Cf Alcolyzer (% v/v)"
    primary["ethanol_real_original_percent_vv"] = primary[ethanol_column]
    artificial_negative = primary[ethanol_column].eq(-0.3) & primary[raw_ethanol_column].isna()
    primary.loc[artificial_negative, ethanol_column] = 0.0
    primary["ethanol_correction"] = np.where(
        artificial_negative,
        "owner_confirmed_negative_formula_artifact_set_to_zero",
        "none",
    )
    for field in ["lot", "reactor", "condition", "replicate", "protocol"]:
        primary[field] = primary["experiment_id"].map(
            {run: metadata[field] for run, metadata in run_metadata.items()}
        )
    primary["initial_volume_l"] = initial_volume_l
    primary["yeast"] = yeast
    primary = primary.sort_values(["experiment_id", "timestamp", "Código muestra"])

    rows = []
    for run, group in primary.groupby("experiment_id", sort=True):
        meta = run_metadata[run]
        rows.append(
            {
                "experiment_id": run,
                "lot": meta["lot"],
                "reactor": meta["reactor"],
                "condition": meta["condition"],
                "replicate": meta["replicate"],
                "protocol": meta["protocol"],
                "temperature_program": meta["temperature_program"],
                "sampling_start": group["timestamp"].min(),
                "sampling_end": group["timestamp"].max(),
                "primary_samples": int(group["timestamp"].notna().sum()),
                "mapping_status": "confirmed",
            }
        )
    windows = pd.DataFrame(rows).sort_values("experiment_id")
    if len(windows) != 9 or windows["experiment_id"].nunique() != 9:
        raise ValueError("Expected exactly nine Pilot 2026 sampling windows")
    primary["sampling_start"] = primary["experiment_id"].map(
        windows.set_index("experiment_id")["sampling_start"]
    )
    primary["time_h"] = (
        primary["timestamp"] - primary["sampling_start"]
    ).dt.total_seconds() / 3600.0
    primary.to_csv(output / "primary_results_qc.csv", index=False)
    return primary, windows


def finalize_primary_process_phase(
    primary: pd.DataFrame, windows: pd.DataFrame, output: Path
) -> pd.DataFrame:
    primary = primary.copy()
    lookup = windows.set_index("experiment_id")
    primary["active_end"] = primary["experiment_id"].map(lookup["active_end"])
    active_end = pd.to_datetime(primary["active_end"])
    sampling_end = pd.to_datetime(primary["experiment_id"].map(lookup["sampling_end"]))
    endpoint_is_active = active_end.ge(sampling_end)
    primary["calibration_include"] = np.where(
        endpoint_is_active,
        primary["timestamp"].le(active_end),
        primary["timestamp"].lt(active_end),
    )
    primary["process_phase"] = np.where(
        primary["calibration_include"], "active_process", "postprocess_cooling"
    )
    primary.to_csv(output / "primary_results_qc.csv", index=False)
    return primary


def export_temperature_source(raw_zip: Path, output: Path) -> pd.DataFrame:
    exporter = WORKSPACE / "export_temperature.mjs"
    dependency = WORKSPACE / "node_modules" / "mdb-reader"
    if not dependency.exists():
        raise RuntimeError(
            "Missing mdb-reader. Run `npm install` in fermentation_model/pilot_2026."
        )
    with tempfile.TemporaryDirectory(prefix="pilot2026_mdb_") as temp_dir:
        temp = Path(temp_dir)
        with zipfile.ZipFile(raw_zip) as archive:
            if archive.testzip() is not None:
                raise ValueError("Temperature-controller ZIP failed CRC validation")
            for member in archive.infolist():
                if member.filename.lower().endswith(".mdb"):
                    member.filename = Path(member.filename).name
                    archive.extract(member, temp)
        combined = temp / "temperature_all_tables.csv"
        subprocess.run(
            ["node", str(exporter), str(temp), str(combined)], check=True, capture_output=True
        )
        frame = pd.read_csv(combined, low_memory=False)
    return frame


def process_temperature(
    raw_zip: Path, windows: pd.DataFrame, runs: list[str], output: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    all_tables = export_temperature_source(raw_zip, output)
    all_tables["experiment_id"] = all_tables["source_file"].str.extract(
        RUN_ID_RE.pattern, expand=False
    )
    all_tables = all_tables[all_tables["experiment_id"].isin(runs)].copy()

    readings = all_tables[all_tables["source_table"].eq("LETTURE_CTA")].copy()
    lossless = readings[
        [
            "experiment_id",
            "source_file",
            "source_table",
            "source_row",
            "Data_ora",
            "Time",
            "Modo",
            "Temperatura",
            "Sonda1",
            "Sonda2",
            "IdEntita",
        ]
    ].rename(
        columns={
            "Data_ora": "timestamp_raw",
            "Modo": "mode",
            "Temperatura": "controller_temperature_c",
            "Sonda1": "sensor_1_c",
            "Sonda2": "sensor_2_c",
        }
    )
    lossless.to_csv(output / "temperature_controller_lossless_export.csv", index=False)

    readings["timestamp"] = (
        pd.to_datetime(readings["Data_ora"], utc=True, errors="coerce")
        .dt.tz_localize(None)
        .dt.floor("s")
    )
    readings = readings.sort_values(["experiment_id", "timestamp", "source_row"])
    duplicate_second_rows = int(readings.duplicated(["experiment_id", "timestamp"]).sum())
    readings = readings.drop_duplicates(["experiment_id", "timestamp"], keep="first")

    events = all_tables[all_tables["source_table"].eq("IMPOSTA_CTA")].copy()
    events["timestamp"] = pd.to_datetime(
        events["Data_ora"], utc=True, errors="coerce"
    ).dt.tz_localize(None)
    events["setpoint_c"] = (events["TempMin"] + events["TempMax"]) / 2.0
    events = events.sort_values(["experiment_id", "timestamp", "source_row"])

    windows_by_run = windows.set_index("experiment_id")
    active_end: dict[str, pd.Timestamp] = {}
    for run in runs:
        start = pd.Timestamp(windows_by_run.loc[run, "sampling_start"])
        end = pd.Timestamp(windows_by_run.loc[run, "sampling_end"])
        cooling = events[
            events["experiment_id"].eq(run)
            & events["timestamp"].ge(start + pd.Timedelta(hours=24))
            & events["Modo"].eq(2)
            & events["TempMin"].eq(9)
            & events["TempMax"].eq(11)
        ]
        active_end[run] = min(cooling["timestamp"].min(), end) if len(cooling) else end

    event_columns = [
        "experiment_id",
        "source_file",
        "source_row",
        "timestamp",
        "Modo",
        "TempMin",
        "TempMax",
        "setpoint_c",
        "CodiceOperatore",
    ]
    transitions = events[event_columns].rename(columns={"Modo": "mode"})
    transitions["sampling_start"] = transitions["experiment_id"].map(
        windows_by_run["sampling_start"]
    )
    transitions["sampling_end"] = transitions["experiment_id"].map(
        windows_by_run["sampling_end"]
    )
    transitions["active_end"] = transitions["experiment_id"].map(active_end)
    transitions["inside_sampling_window"] = transitions["timestamp"].between(
        transitions["sampling_start"], transitions["sampling_end"], inclusive="both"
    )
    transitions["calibration_include"] = transitions["timestamp"].lt(
        transitions["active_end"]
    ) & transitions["inside_sampling_window"]
    transitions.to_csv(output / "temperature_transitions_qc.csv", index=False)

    qc_frames = []
    for run in runs:
        run_readings = readings[readings["experiment_id"].eq(run)][
            [
                "timestamp",
                "source_file",
                "source_row",
                "Modo",
                "Temperatura",
                "Sonda1",
                "Sonda2",
            ]
        ].copy()
        run_events = events[events["experiment_id"].eq(run)][
            ["timestamp", "Modo", "TempMin", "TempMax", "setpoint_c"]
        ].copy()
        run_readings = pd.merge_asof(
            run_readings.sort_values("timestamp"),
            run_events.sort_values("timestamp"),
            on="timestamp",
            direction="backward",
            suffixes=("_reading", "_event"),
        )
        run_readings["experiment_id"] = run
        start = pd.Timestamp(windows_by_run.loc[run, "sampling_start"])
        end = pd.Timestamp(windows_by_run.loc[run, "sampling_end"])
        run_readings["inside_sampling_window"] = run_readings["timestamp"].between(
            start, end, inclusive="both"
        )
        run_readings["active_process"] = run_readings["timestamp"].lt(active_end[run]) & run_readings[
            "inside_sampling_window"
        ]
        qc_frames.append(run_readings)
    temperature_qc = pd.concat(qc_frames, ignore_index=True)
    keep = [
        "experiment_id",
        "timestamp",
        "source_file",
        "source_row",
        "Modo_reading",
        "controller_temperature_c",
        "Sonda1",
        "Sonda2",
        "setpoint_c",
        "TempMin",
        "TempMax",
        "inside_sampling_window",
        "active_process",
    ]
    temperature_qc = temperature_qc.rename(
        columns={"Temperatura": "controller_temperature_c", "Sonda1": "sensor_1_c", "Sonda2": "sensor_2_c"}
    )
    keep = [column.replace("Sonda1", "sensor_1_c").replace("Sonda2", "sensor_2_c") for column in keep]
    temperature_qc[keep].to_csv(output / "temperature_controller_qc.csv", index=False)

    summaries = []
    for run in runs:
        group = temperature_qc[
            temperature_qc["experiment_id"].eq(run) & temperature_qc["inside_sampling_window"]
        ]
        error = (group["sensor_1_c"] - group["setpoint_c"]).abs()
        summaries.append(
            {
                "experiment_id": run,
                "first_temperature_timestamp": group["timestamp"].min(),
                "last_temperature_timestamp": group["timestamp"].max(),
                "temperature_rows_in_sampling_window": len(group),
                "mean_absolute_setpoint_error_c": error.mean(),
                "p95_absolute_setpoint_error_c": error.quantile(0.95),
                "active_end": active_end[run],
                "active_end_source": "IMPOSTA_CTA_mode2_9_to_11C" if active_end[run] < pd.Timestamp(windows_by_run.loc[run, "sampling_end"]) else "last_primary_sample",
            }
        )
    temperature_summary = pd.DataFrame(summaries)
    temperature_summary.to_csv(output / "temperature_summary.csv", index=False)

    stats = {
        "raw_relevant_temperature_rows": int(len(lossless)),
        "duplicate_temperature_seconds_removed": duplicate_second_rows,
        "normalized_temperature_rows": int(len(readings)),
        "active_end": {run: active_end[run].isoformat() for run in runs},
    }
    return temperature_qc, transitions, temperature_summary, stats


def process_co2(
    raw_zip: Path,
    windows: pd.DataFrame,
    runs: list[str],
    active_end: dict[str, str],
    invalid_value: float,
    model_excluded_runs: set[str],
    output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    windows_by_run = windows.set_index("experiment_id")
    type_counts: Counter[str] = Counter()
    total_rows = 0
    rows_inside = 0
    sentinel_global = 0
    before_groups = 0
    after_rows = 0
    after_groups = 0
    out_of_order = 0
    maximum_gap = 0.0
    minute_frames = []
    summaries = []

    with zipfile.ZipFile(raw_zip) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise ValueError(f"CO2 ZIP failed CRC validation at {bad_member}")
        members = sorted(name for name in archive.namelist() if name.lower().endswith(".csv"))
        for member in members:
            match = RUN_ID_RE.search(member)
            if not match or match.group(1) not in runs:
                continue
            run = match.group(1)
            start = pd.Timestamp(windows_by_run.loc[run, "sampling_start"])
            end = pd.Timestamp(windows_by_run.loc[run, "sampling_end"])
            biological_end = pd.Timestamp(active_end[run])
            selected = []
            previous_timestamp = None
            run_type_counts: Counter[str] = Counter()
            run_sentinel_global = 0
            for chunk in pd.read_csv(archive.open(member), chunksize=250_000):
                total_rows += len(chunk)
                counts = chunk["row_type"].value_counts().to_dict()
                type_counts.update(counts)
                run_type_counts.update(counts)
                flow = chunk[chunk["row_type"].eq("flow")].copy()
                flow["timestamp"] = pd.to_datetime(flow["timestamp"], errors="raise")
                if previous_timestamp is not None and len(flow):
                    out_of_order += int(flow["timestamp"].iloc[0] < previous_timestamp)
                differences = flow["timestamp"].diff().dt.total_seconds().dropna()
                out_of_order += int(differences.lt(0).sum())
                if len(differences):
                    maximum_gap = max(maximum_gap, float(differences.max()))
                if len(flow):
                    previous_timestamp = flow["timestamp"].iloc[-1]
                raw_sentinel = int(flow["fmeasure"].eq(invalid_value).sum())
                sentinel_global += raw_sentinel
                run_sentinel_global += raw_sentinel
                inside = flow[flow["timestamp"].between(start, end, inclusive="both")][
                    ["timestamp", "fmeasure"]
                ]
                rows_inside += len(inside)
                selected.append(inside)

            raw_window = pd.concat(selected, ignore_index=True)
            raw_window["second"] = raw_window["timestamp"].dt.floor("s")
            raw_group_sizes = raw_window.groupby("second", sort=False).size()
            before_groups += len(raw_group_sizes)
            valid_raw = raw_window[raw_window["fmeasure"].ne(invalid_value)].copy()
            consolidated = valid_raw.groupby("second", sort=True)["fmeasure"].agg(
                flow_ln_min="mean", raw_records="size"
            )
            after_rows += len(valid_raw)
            after_groups += len(consolidated)

            second = consolidated.index
            # A controller cooling event starts postprocess at its event second.
            # If no cooling event occurs before the sampling endpoint, the final
            # primary-sample second remains part of the active sampling window.
            active_mask = (
                second <= biological_end
                if biological_end >= end
                else second < biological_end
            )
            postprocess_mask = ~active_mask
            positive = consolidated.loc[active_mask & consolidated["flow_ln_min"].gt(0)]
            if len(positive):
                first_positive = positive.index.min()
                last_positive = positive.index.max()
            else:
                first_positive = pd.NaT
                last_positive = pd.NaT
            state = pd.Series("valid_nonzero", index=second, dtype="object")
            zero = consolidated["flow_ln_min"].eq(0)
            if pd.notna(first_positive):
                edge_zero = zero & active_mask & ((second < first_positive) | (second > last_positive))
                intermediate_zero = zero & active_mask & ~edge_zero
            else:
                edge_zero = zero & active_mask
                intermediate_zero = pd.Series(False, index=second)
            state.loc[edge_zero] = "valid_observed_edge_zero"
            state.loc[intermediate_zero] = "invalid_intermediate_zero"
            state.loc[postprocess_mask] = "postprocess_cooling"

            seconds = consolidated.reset_index().rename(columns={"second": "timestamp"})
            seconds["experiment_id"] = run
            seconds["qc_state"] = state.to_numpy()
            seconds["artificial"] = False
            seconds["calibration_eligible"] = seconds["qc_state"].isin(
                ["valid_nonzero", "valid_observed_edge_zero"]
            )
            if run in model_excluded_runs:
                seconds["calibration_eligible"] = False

            artificial_initial = max(
                0, int((consolidated.index.min() - start).total_seconds())
            )
            artificial_final = max(
                0, int((end - consolidated.index.max()).total_seconds())
            )
            artificial_parts = []
            if artificial_initial:
                artificial_parts.append(
                    pd.DataFrame(
                        {
                            "timestamp": pd.date_range(
                                start, consolidated.index.min() - pd.Timedelta(seconds=1), freq="s"
                            ),
                            "flow_ln_min": 0.0,
                            "raw_records": 0,
                            "experiment_id": run,
                            "qc_state": "artificial_initial_zero",
                            "artificial": True,
                            "calibration_eligible": run not in model_excluded_runs,
                        }
                    )
                )
            if artificial_final:
                artificial_parts.append(
                    pd.DataFrame(
                        {
                            "timestamp": pd.date_range(
                                consolidated.index.max() + pd.Timedelta(seconds=1), end, freq="s"
                            ),
                            "flow_ln_min": 0.0,
                            "raw_records": 0,
                            "experiment_id": run,
                            "qc_state": "artificial_final_zero",
                            "artificial": True,
                            "calibration_eligible": run not in model_excluded_runs,
                        }
                    )
                )
            if artificial_parts:
                seconds = pd.concat([seconds, *artificial_parts], ignore_index=True)
            seconds = seconds.sort_values("timestamp")
            seconds["minute"] = seconds["timestamp"].dt.floor("min")
            seconds["flow_for_calibration"] = seconds["flow_ln_min"].where(
                seconds["calibration_eligible"]
            )
            seconds["observed_value"] = seconds["flow_ln_min"].where(~seconds["artificial"])

            base = seconds.groupby(["experiment_id", "minute"], as_index=False).agg(
                mean_flow_ln_min_calibration=("flow_for_calibration", "mean"),
                mean_flow_ln_min_observed=("observed_value", "mean"),
                consolidated_or_filled_seconds=("timestamp", "size"),
                raw_records_averaged=("raw_records", "sum"),
                calibration_seconds=("calibration_eligible", "sum"),
            )
            counts = (
                seconds.groupby(["experiment_id", "minute", "qc_state"])
                .size()
                .unstack(fill_value=0)
                .reset_index()
            )
            minute = base.merge(counts, on=["experiment_id", "minute"], how="left")
            minute["co2_model_include"] = run not in model_excluded_runs
            minute["co2_model_exclusion_reason"] = (
                "owner_confirmed_poor_signal_lot1"
                if run in model_excluded_runs
                else ""
            )
            minute["expected_seconds"] = [
                int(
                    (
                        min(end, timestamp + pd.Timedelta(seconds=59))
                        - max(start, timestamp)
                    ).total_seconds()
                )
                + 1
                for timestamp in minute["minute"]
            ]
            minute["unimputed_missing_seconds"] = (
                minute["expected_seconds"] - minute["consolidated_or_filled_seconds"]
            )
            minute_frames.append(minute)

            state_counts = state.value_counts().to_dict()
            summaries.append(
                {
                    "experiment_id": run,
                    "source_member": member,
                    "total_source_rows": sum(run_type_counts.values()),
                    "flow_rows_in_sampling_window": len(raw_window),
                    "same_second_groups_before_sentinel_filter": len(raw_group_sizes),
                    "same_second_repetitions_absorbed_before_sentinel_filter": len(raw_window) - len(raw_group_sizes),
                    "raw_2_621_global": run_sentinel_global,
                    "valid_raw_rows_after_sentinel_filter_in_window": len(valid_raw),
                    "consolidated_seconds_after_sentinel_filter": len(consolidated),
                    "repetitions_absorbed_after_sentinel_filter": len(valid_raw) - len(consolidated),
                    "sentinel_only_seconds_removed": len(raw_group_sizes) - len(consolidated),
                    "invalid_intermediate_zero_seconds": int(state_counts.get("invalid_intermediate_zero", 0)),
                    "valid_observed_edge_zero_seconds": int(state_counts.get("valid_observed_edge_zero", 0)),
                    "postprocess_observed_seconds": int(state_counts.get("postprocess_cooling", 0)),
                    "artificial_initial_zero_seconds": artificial_initial,
                    "artificial_final_zero_seconds": artificial_final,
                    "first_flow_timestamp": raw_window["timestamp"].min(),
                    "last_flow_timestamp": raw_window["timestamp"].max(),
                    "modbus_events": int(run_type_counts.get("error", 0)),
                    "co2_model_include": run not in model_excluded_runs,
                    "co2_model_exclusion_reason": (
                        "owner_confirmed_poor_signal_lot1"
                        if run in model_excluded_runs
                        else ""
                    ),
                }
            )

    minute_qc = pd.concat(minute_frames, ignore_index=True).sort_values(
        ["experiment_id", "minute"]
    )
    for state_column in [
        "valid_nonzero",
        "valid_observed_edge_zero",
        "invalid_intermediate_zero",
        "postprocess_cooling",
        "artificial_initial_zero",
        "artificial_final_zero",
    ]:
        if state_column not in minute_qc:
            minute_qc[state_column] = 0
    minute_qc.to_csv(output / "co2_minute_qc.csv", index=False)
    run_summary = pd.DataFrame(summaries).sort_values("experiment_id")
    run_summary.to_csv(output / "co2_run_qc_summary.csv", index=False)

    stats = {
        "massview_total_rows": total_rows,
        "massview_row_types": dict(type_counts),
        "flow_rows_in_sampling_windows": rows_inside,
        "same_second_groups_before_sentinel_filter": before_groups,
        "same_second_repetitions_absorbed_before_sentinel_filter": rows_inside - before_groups,
        "raw_2_621_observations_global": sentinel_global,
        "valid_raw_rows_after_sentinel_filter_in_windows": after_rows,
        "same_second_groups_after_sentinel_filter": after_groups,
        "same_second_repetitions_absorbed_after_sentinel_filter": after_rows - after_groups,
        "sentinel_only_seconds_removed": before_groups - after_groups,
        "out_of_order_flow_timestamps": out_of_order,
        "maximum_flow_gap_seconds": maximum_gap,
        "invalid_intermediate_zero_seconds": int(run_summary["invalid_intermediate_zero_seconds"].sum()),
        "valid_observed_edge_zero_seconds": int(run_summary["valid_observed_edge_zero_seconds"].sum()),
        "postprocess_observed_seconds": int(run_summary["postprocess_observed_seconds"].sum()),
        "artificial_initial_zero_seconds": int(run_summary["artificial_initial_zero_seconds"].sum()),
        "artificial_final_zero_seconds": int(run_summary["artificial_final_zero_seconds"].sum()),
        "modbus_events_26136": int(run_summary.loc[run_summary["experiment_id"].eq("26136"), "modbus_events"].sum()),
        "model_excluded_runs": sorted(model_excluded_runs),
    }
    return minute_qc, run_summary, stats


def process_operational_events(
    master: Path,
    windows: pd.DataFrame,
    runs: list[str],
    lot3_reference_map: dict[str, str],
    initial_volume_l: float,
    organic_yan_mg_per_mg: float,
    dap_yan_mg_per_mg: float,
    output: Path,
) -> pd.DataFrame:
    events = read_master_sheet(master, "03_Eventos_operacion")
    events["experiment_id"] = events["Código ensayo"].astype(str).str.extract(
        r"^(26\d{3})", expand=False
    )
    events = events[events["experiment_id"].isin(runs)].copy()
    events["timestamp_original"] = events["Fecha-hora"].astype(str)
    events["timestamp"] = [
        normalize_datetime(value, run)
        for value, run in zip(events["Fecha-hora"], events["experiment_id"])
    ]
    lookup = windows.set_index("experiment_id")
    events["event_origin"] = "master_workbook"
    events["reference_experiment_id"] = events["experiment_id"]
    events["event_timestamp_status"] = np.select(
        [
            pd.to_datetime(events["Fecha de ejecución"], errors="coerce").notna(),
            events["experiment_id"].isin({"26134", "26135", "26136"}),
        ],
        ["recorded_execution", "owner_confirmed_master_schedule"],
        default="master_schedule_execution_time_unconfirmed",
    )

    reconstructed = []
    for target, reference in lot3_reference_map.items():
        reference_events = events[events["experiment_id"].eq(reference)].copy()
        if reference_events.empty:
            raise ValueError(f"Missing Lot 1 reference events for {reference}")
        reference_start = pd.Timestamp(lookup.loc[reference, "sampling_start"])
        target_start = pd.Timestamp(lookup.loc[target, "sampling_start"])
        relative_time = reference_events["timestamp"] - reference_start
        reference_events["experiment_id"] = target
        reference_events["Código ensayo"] = target
        reference_events["Código evento"] = reference_events["Código evento"].astype(str).str.replace(
            reference, target, regex=False
        )
        reference_events["timestamp"] = target_start + relative_time
        reference_events["timestamp_original"] = reference_events["timestamp"].astype(str)
        reference_events["Fecha-hora"] = reference_events["timestamp"]
        reference_events["Fecha"] = reference_events["timestamp"].dt.date
        reference_events["Hora"] = reference_events["timestamp"].dt.time
        reference_events["Fecha de ejecución"] = pd.NaT
        reference_events["Estado"] = "Reconstruido confirmado"
        reference_events["event_origin"] = "reconstructed_from_lot1_relative_schedule"
        reference_events["reference_experiment_id"] = reference
        reference_events["event_timestamp_status"] = "owner_confirmed_relative_reconstruction"
        reconstructed.append(reference_events)
    if reconstructed:
        events = pd.concat([events, *reconstructed], ignore_index=True)

    events["sampling_start"] = events["experiment_id"].map(lookup["sampling_start"])
    events["sampling_end"] = events["experiment_id"].map(lookup["sampling_end"])
    events["active_end"] = events["experiment_id"].map(lookup["active_end"])
    events["relative_time_h"] = (
        events["timestamp"] - pd.to_datetime(events["sampling_start"])
    ).dt.total_seconds() / 3600.0
    events["inside_sampling_window"] = events["timestamp"].between(
        events["sampling_start"], events["sampling_end"], inclusive="both"
    )
    active_end = pd.to_datetime(events["active_end"])
    sampling_end = pd.to_datetime(events["sampling_end"])
    events["inside_active_process"] = np.where(
        active_end.ge(sampling_end),
        events["timestamp"].le(active_end),
        events["timestamp"].lt(active_end),
    )
    events["calibration_include"] = (
        events["inside_sampling_window"] & events["inside_active_process"]
    )

    nutrient = events["Tipo evento"].astype(str).str.contains("Pulso nutricional", na=False)
    doses = events["Dosis"].astype(str).str.extract(
        r"^\s*([0-9]+(?:[.,][0-9]+)?)\s*\+\s*([0-9]+(?:[.,][0-9]+)?)\s*$"
    )
    events["organic_product_g"] = pd.to_numeric(
        doses[0].str.replace(",", ".", regex=False), errors="coerce"
    ).where(nutrient)
    events["dap_product_g"] = pd.to_numeric(
        doses[1].str.replace(",", ".", regex=False), errors="coerce"
    ).where(nutrient)
    events["yan_added_mg_l"] = (
        events["organic_product_g"] * 1000.0 * organic_yan_mg_per_mg
        + events["dap_product_g"] * 1000.0 * dap_yan_mg_per_mg
    ) / initial_volume_l
    events["dose_parse_status"] = np.where(
        nutrient & events["organic_product_g"].notna() & events["dap_product_g"].notna(),
        "parsed_springferm_organic_plus_fda",
        np.where(nutrient, "unparsed", "not_applicable"),
    )
    events["model_input_note"] = np.where(
        events["Tipo evento"].astype(str).str.contains("Cambio temperatura", na=False),
        "use_measured_controller_setpoint_and_sensor_trace",
        np.where(
            events["Tipo evento"].astype(str).eq("Pulso nutricional 2"),
            "density_triggered_time_proxy",
            "event_table",
        ),
    )
    events = events.sort_values(["experiment_id", "timestamp", "Código evento"])
    events.to_csv(output / "operational_events_qc.csv", index=False)
    return events


def process_gc(
    master: Path,
    report: Path,
    runs: list[str],
    windows: pd.DataFrame,
    output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    pilot_gc_runs = set(runs) - {"26134", "26135", "26136"}
    window_lookup = windows.set_index("experiment_id")
    aromas = read_master_sheet(master, "09_Aromas_Cond")
    aromas = aromas[aromas["Código condensado / mix"].notna()].copy()
    aromas["mix_id"] = aromas["Código condensado / mix"].astype(str).str.strip()
    aromas["experiment_id"] = aromas["mix_id"].str.extract(r"^(26\d{3})", expand=False)
    aromas = aromas[aromas["experiment_id"].isin(runs)].copy()
    aromas["paired_wine_sample_id"] = aromas["Código muestra MEF"].map(
        canonical_primary_sample_id
    )
    aromas["timestamp_original"] = aromas["Fecha-hora"].astype(str)
    aromas["timestamp"] = [
        normalize_datetime(value, run)
        for value, run in zip(aromas["Fecha-hora"], aromas["experiment_id"])
    ]
    aromas["date_normalization"] = np.where(
        aromas["timestamp_original"] != aromas["timestamp"].astype(str),
        "parsed_or_day_month_corrected",
        "unchanged",
    )
    aromas["volume_a_ml"] = pd.to_numeric(aromas["Vol CA (mL)"], errors="coerce")
    aromas["volume_b_ml"] = pd.to_numeric(aromas["Vol CB (mL)"], errors="coerce")
    aromas["total_condensate_ml"] = aromas["volume_a_ml"].fillna(0) + aromas["volume_b_ml"].fillna(0)
    aromas["fraction_a"] = aromas["volume_a_ml"].fillna(0).div(aromas["total_condensate_ml"].replace(0, np.nan))
    aromas["fraction_b"] = aromas["volume_b_ml"].fillna(0).div(aromas["total_condensate_ml"].replace(0, np.nan))
    aromas = aromas.sort_values(["experiment_id", "timestamp", "mix_id"])
    aromas["capture_interval_start"] = aromas.groupby("experiment_id")["timestamp"].shift(1)
    aromas["capture_interval_start"] = aromas["capture_interval_start"].fillna(
        aromas["experiment_id"].map(window_lookup["sampling_start"])
    )
    aromas["capture_interval_end"] = aromas["timestamp"]
    aromas["capture_interval_h"] = (
        aromas["capture_interval_end"] - aromas["capture_interval_start"]
    ).dt.total_seconds() / 3600.0
    aromas["time_h"] = (
        aromas["timestamp"]
        - pd.to_datetime(aromas["experiment_id"].map(window_lookup["sampling_start"]))
    ).dt.total_seconds() / 3600.0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        raw = pd.read_excel(report, sheet_name="Hoja", header=None)
    sample_rows = raw.iloc[28:112]
    sample_map = pd.DataFrame(
        {
            "sample_number": sample_rows.iloc[:, 1],
            "original_sample_id": sample_rows.iloc[:, 2].astype(str).str.strip(),
            "lab_sample_id": sample_rows.iloc[:, 3].astype(str).str.strip(),
        }
    )
    sample_map = sample_map[sample_map["sample_number"].notna()]
    sample_map["experiment_id"] = sample_map["original_sample_id"].str.extract(
        r"^(26\d{3})", expand=False
    )
    sample_map["sample_id"] = sample_map["original_sample_id"].map(
        canonical_primary_sample_id
    )
    sample_map["sample_type"] = np.select(
        [
            sample_map["sample_id"].str.contains("-MIX-", na=False),
            sample_map["sample_id"].str.contains("-P-", na=False),
        ],
        ["condensate_mix", "wine_mef"],
        default="other",
    )

    analytes = [
        "Etil Acetato",
        "Isoamil Acetato",
        "Hexil Acetato",
        "Octanoato de Etilo",
        "Decanoato de Etilo",
        "Feniletil Acetato",
    ]
    result_rows = raw.iloc[117:201]
    result_table = pd.DataFrame({"lab_sample_id": result_rows.iloc[:, 2].astype(str).str.strip()})
    for offset, analyte in enumerate(analytes, start=3):
        result_table[analyte] = result_rows.iloc[:, offset].to_numpy()
    loq = {analyte: parse_number(raw.iloc[201, offset]) for offset, analyte in enumerate(analytes, start=3)}

    combined = sample_map.merge(result_table, on="lab_sample_id", how="left", validate="one_to_one")
    pilot_samples = combined[
        combined["experiment_id"].isin(pilot_gc_runs)
        & combined["sample_type"].isin({"condensate_mix", "wine_mef"})
    ].copy()
    mix_samples = pilot_samples[pilot_samples["sample_type"].eq("condensate_mix")].copy()
    wine_samples = pilot_samples[pilot_samples["sample_type"].eq("wine_mef")].copy()
    long = pilot_samples.melt(
        id_vars=[
            "sample_number",
            "original_sample_id",
            "sample_id",
            "sample_type",
            "lab_sample_id",
            "experiment_id",
        ],
        value_vars=analytes,
        var_name="analyte",
        value_name="reported_value",
    )
    long["vial_loq_ug_l"] = long["analyte"].map(loq)
    long["vial_concentration_ug_l"] = long["reported_value"].map(parse_number)
    long["result_status"] = np.select(
        [
            long["reported_value"].astype(str).str.strip().str.upper().eq("NQ"),
            long["vial_concentration_ug_l"].notna()
            & long["vial_concentration_ug_l"].lt(long["vial_loq_ug_l"]),
            long["vial_concentration_ug_l"].notna(),
        ],
        ["NQ", "below_loq", "quantified"],
        default="missing",
    )
    long["dilution_factor"] = np.where(long["sample_type"].eq("condensate_mix"), 1000, 1)
    long["sample_concentration_ug_l"] = (
        long["vial_concentration_ug_l"] * long["dilution_factor"]
    )
    long["sample_basis_loq_ug_l"] = long["vial_loq_ug_l"] * long["dilution_factor"]
    long["reported_basis"] = np.where(
        long["sample_type"].eq("condensate_mix"),
        "diluted_GC_vial",
        "undiluted_wine_sample",
    )
    long["preparation_stock_ml"] = np.where(
        long["sample_type"].eq("condensate_mix"), 50, np.nan
    )
    long["volume_transferred_to_analyst_ml"] = np.where(
        long["sample_type"].eq("condensate_mix"), 10, np.nan
    )
    long["model_observation_type"] = np.select(
        [
            long["result_status"].isin({"NQ", "below_loq"}),
            long["result_status"].eq("quantified"),
        ],
        ["left_censored", "observed"],
        default="missing",
    )
    long["censoring_lower_bound_ug_l"] = np.where(
        long["model_observation_type"].eq("left_censored"), 0.0, np.nan
    )
    long["censoring_upper_bound_ug_l"] = np.where(
        long["model_observation_type"].eq("left_censored"),
        long["sample_basis_loq_ug_l"],
        np.nan,
    )

    primary_results = read_master_sheet(master, "06_Resultados_primarios")
    primary_results["sample_id"] = primary_results["Código muestra"].astype(str).str.strip()
    primary_results["experiment_id"] = primary_results["sample_id"].str.extract(
        r"^(26\d{3})", expand=False
    )
    primary_results["sample_timestamp"] = [
        normalize_datetime(value, run)
        for value, run in zip(primary_results["Fecha-hora"], primary_results["experiment_id"])
    ]
    wine_timestamp = primary_results.set_index("sample_id")["sample_timestamp"]
    mix_timestamp = aromas.set_index("mix_id")["timestamp"]
    long["sample_timestamp"] = np.where(
        long["sample_type"].eq("condensate_mix"),
        long["sample_id"].map(mix_timestamp),
        long["sample_id"].map(wine_timestamp),
    )
    long["sample_timestamp"] = pd.to_datetime(long["sample_timestamp"])
    long["time_h"] = (
        long["sample_timestamp"]
        - pd.to_datetime(long["experiment_id"].map(window_lookup["sampling_start"]))
    ).dt.total_seconds() / 3600.0
    mix_volume = aromas.set_index("mix_id")["total_condensate_ml"]
    long["total_condensate_ml"] = long["sample_id"].map(mix_volume).where(
        long["sample_type"].eq("condensate_mix")
    )
    long["captured_mass_ug"] = (
        long["sample_concentration_ug_l"] * long["total_condensate_ml"] / 1000.0
    )
    long["captured_mass_upper_bound_ug"] = (
        long["censoring_upper_bound_ug_l"] * long["total_condensate_ml"] / 1000.0
    )
    long.to_csv(output / "gc_results_long_qc.csv", index=False)
    long[long["sample_type"].eq("wine_mef")].to_csv(
        output / "gc_wine_results_long_qc.csv", index=False
    )

    analyzed = aromas[aromas["experiment_id"].isin(pilot_gc_runs)].copy()
    analyzed = analyzed.merge(
        mix_samples[["sample_id", "lab_sample_id"]],
        left_on="mix_id",
        right_on="sample_id",
        how="left",
        validate="one_to_one",
    )
    analyzed["gc_prepared"] = analyzed["lab_sample_id"].notna()
    analyzed["gc_scope"] = "confirmed_lots_2_and_3_only"
    analyzed["paired_wine_gc_lab_sample_id"] = analyzed["paired_wine_sample_id"].map(
        wine_samples.set_index("sample_id")["lab_sample_id"]
    )
    analyzed["wine_pair_available"] = analyzed["paired_wine_gc_lab_sample_id"].notna()
    analyzed.to_csv(output / "gc_condensate_mix_qc.csv", index=False)

    pair_columns = [
        "experiment_id",
        "mix_id",
        "paired_wine_sample_id",
        "timestamp",
        "capture_interval_start",
        "capture_interval_end",
        "capture_interval_h",
        "time_h",
        "volume_a_ml",
        "volume_b_ml",
        "total_condensate_ml",
        "fraction_a",
        "fraction_b",
        "lab_sample_id",
        "paired_wine_gc_lab_sample_id",
        "wine_pair_available",
    ]
    pairs = analyzed[pair_columns].copy()
    pairs.to_csv(output / "gc_mix_wine_pairs_qc.csv", index=False)

    mix_long = long[long["sample_type"].eq("condensate_mix")][
        [
            "sample_id",
            "analyte",
            "result_status",
            "model_observation_type",
            "vial_concentration_ug_l",
            "sample_concentration_ug_l",
            "sample_basis_loq_ug_l",
            "captured_mass_ug",
            "captured_mass_upper_bound_ug",
        ]
    ].rename(
        columns={
            "sample_id": "mix_id_result",
            "result_status": "mix_result_status",
            "model_observation_type": "mix_model_observation_type",
            "vial_concentration_ug_l": "mix_vial_concentration_ug_l",
            "sample_concentration_ug_l": "mix_concentration_ug_l",
            "sample_basis_loq_ug_l": "mix_loq_ug_l",
        }
    )
    wine_long = long[long["sample_type"].eq("wine_mef")][
        [
            "sample_id",
            "analyte",
            "result_status",
            "model_observation_type",
            "sample_concentration_ug_l",
            "sample_basis_loq_ug_l",
        ]
    ].rename(
        columns={
            "sample_id": "wine_sample_id_result",
            "result_status": "wine_result_status",
            "model_observation_type": "wine_model_observation_type",
            "sample_concentration_ug_l": "wine_concentration_ug_l",
            "sample_basis_loq_ug_l": "wine_loq_ug_l",
        }
    )
    pair_long = pairs.merge(
        mix_long,
        left_on="mix_id",
        right_on="mix_id_result",
        how="left",
        validate="one_to_many",
    ).merge(
        wine_long,
        left_on=["paired_wine_sample_id", "analyte"],
        right_on=["wine_sample_id_result", "analyte"],
        how="left",
        validate="many_to_one",
    )
    pair_long.to_csv(output / "gc_mix_wine_pairs_long_qc.csv", index=False)

    stats = {
        "gc_report_samples_total": int(len(sample_map)),
        "gc_mix_samples": int(len(mix_samples)),
        "gc_wine_samples": int(len(wine_samples)),
        "gc_initial_wine_samples_without_mix": int(
            len(set(wine_samples["sample_id"]) - set(analyzed["paired_wine_sample_id"]))
        ),
        "gc_mix_wine_pairs": int(len(pairs)),
        "gc_mix_wine_pair_coverage": float(pairs["wine_pair_available"].mean()),
        "gc_analyte_results": int(len(long)),
        "gc_mix_analyte_results": int(long["sample_type"].eq("condensate_mix").sum()),
        "gc_wine_analyte_results": int(long["sample_type"].eq("wine_mef").sum()),
        "gc_paired_analyte_results": int(len(pair_long)),
        "gc_numeric_results": int(long["vial_concentration_ug_l"].notna().sum()),
        "gc_numeric_mix_results": int(
            (long["sample_type"].eq("condensate_mix") & long["vial_concentration_ug_l"].notna()).sum()
        ),
        "gc_status_counts": {
            f"{sample_type}:{analyte}:{status}": int(count)
            for (sample_type, analyte, status), count in long.groupby(
                ["sample_type", "analyte", "result_status"]
            ).size().items()
        },
        "lot1_gc_mix_samples_in_model_scope": 0,
    }
    return analyzed, long, stats


def make_figures(
    minute: pd.DataFrame,
    temperature: pd.DataFrame,
    gc_long: pd.DataFrame,
    windows: pd.DataFrame,
    runs: list[str],
    output: Path,
) -> None:
    lookup = windows.set_index("experiment_id")
    fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharey=False)
    for axis, run in zip(axes.flat, runs):
        group = minute[minute["experiment_id"].eq(run)].copy()
        hours = (group["minute"] - pd.Timestamp(lookup.loc[run, "sampling_start"])).dt.total_seconds() / 3600
        axis.plot(hours, group["mean_flow_ln_min_observed"], color="#4c78a8", lw=0.8, label="observed")
        invalid = group["invalid_intermediate_zero"].gt(0)
        axis.scatter(hours[invalid], np.zeros(invalid.sum()), color="#e45756", s=4, label="invalid zero")
        axis.set_title(f"{run} · {lookup.loc[run, 'reactor']} · {lookup.loc[run, 'protocol']}")
        axis.set_xlabel("Hours from first primary sample")
        axis.set_ylabel("CO₂ (Ln/min)")
        axis.grid(alpha=0.2)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2)
    fig.suptitle("Pilot 2026 CO₂ profiles and intermediate-zero QC", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.935))
    atomic_savefig(fig, output / "co2_profiles_qc.png")
    plt.close(fig)

    heat_source = minute.copy()
    heat_source["hour_from_start"] = [
        int((timestamp - pd.Timestamp(lookup.loc[run, "sampling_start"])).total_seconds() // 3600)
        for run, timestamp in zip(heat_source["experiment_id"], heat_source["minute"])
    ]
    heat = heat_source.pivot_table(
        index="experiment_id", columns="hour_from_start",
        values="invalid_intermediate_zero", aggfunc="sum", fill_value=0
    )
    heat = heat.reindex(runs).sort_index(axis=1)
    fig, axis = plt.subplots(figsize=(16, 4.8))
    image = axis.imshow(np.log1p(heat.to_numpy()), aspect="auto", cmap="magma")
    axis.set_yticks(range(len(heat.index)), heat.index)
    axis.set_xlabel("Hours from first primary sample")
    axis.set_ylabel("Experiment")
    axis.set_title("Intermediate zero seconds (log1p scale)")
    fig.colorbar(image, ax=axis, label="log(1 + invalid zero seconds)")
    fig.tight_layout()
    atomic_savefig(fig, output / "co2_intermediate_zero_heatmap.png")
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharey=True)
    for axis, run in zip(axes.flat, runs):
        group = temperature[
            temperature["experiment_id"].eq(run) & temperature["inside_sampling_window"]
        ].copy()
        hours = (group["timestamp"] - pd.Timestamp(lookup.loc[run, "sampling_start"])).dt.total_seconds() / 3600
        axis.plot(hours, group["sensor_1_c"], color="#4c78a8", lw=1, label="sensor")
        axis.step(hours, group["setpoint_c"], color="#f58518", lw=1.2, where="post", label="setpoint")
        axis.set_title(f"{run} · protocol {lookup.loc[run, 'protocol']}")
        axis.set_xlabel("Hours from first primary sample")
        axis.set_ylabel("Temperature (°C)")
        axis.grid(alpha=0.2)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2)
    fig.suptitle("Temperature sensor versus controller setpoint", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.935))
    atomic_savefig(fig, output / "temperature_sensor_vs_setpoint.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5), sharey=True)
    colors = ["#54a24b", "#eeca3b", "#e45756", "#bab0ac"]
    for axis, (sample_type, title) in zip(
        axes, [("wine_mef", "Wine MEF"), ("condensate_mix", "Condensate MIX")]
    ):
        status = (
            gc_long[gc_long["sample_type"].eq(sample_type)]
            .groupby(["analyte", "result_status"])
            .size()
            .unstack(fill_value=0)
            .reindex(columns=["quantified", "below_loq", "NQ", "missing"], fill_value=0)
        )
        status.plot.bar(stacked=True, ax=axis, color=colors, legend=False)
        axis.set_ylabel("Results")
        axis.set_xlabel("")
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=25)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, title="Status", loc="upper center", ncol=4)
    fig.suptitle("Pilot GC results by matrix and analytical status", y=1.02)
    fig.tight_layout()
    atomic_savefig(fig, output / "gc_condensate_result_status.png")
    plt.close(fig)


def build_manifest(raw_files: list[Path], output: Path) -> dict[str, object]:
    sources = [
        {"path": path.relative_to(REPO).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in raw_files
    ]
    outputs = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            outputs.append(
                {"path": path.relative_to(REPO).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
            )
    pipeline_path = Path(__file__).resolve()
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO, check=True, capture_output=True, text=True
    ).stdout.strip()
    git_status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    node_version = subprocess.run(
        ["node", "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    manifest = {
        "pipeline": {
            "path": pipeline_path.relative_to(REPO).as_posix(),
            "sha256": sha256(pipeline_path),
        },
        "config": {
            "path": CONFIG_PATH.relative_to(REPO).as_posix(),
            "sha256": sha256(CONFIG_PATH),
        },
        "execution_context": {
            "git_head_at_run": git_head,
            "git_worktree_dirty_at_run": bool(git_status),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "node": node_version,
            "packages": {
                name: importlib.metadata.version(name)
                for name in ["numpy", "pandas", "matplotlib", "openpyxl"]
            },
        },
        "sources": sources,
        "outputs": outputs,
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    raw_directory = REPO / config["raw_directory"]
    output = REPO / config["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    sources = {name: raw_directory / filename for name, filename in config["sources"].items()}
    missing = [str(path) for path in sources.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing Pilot 2026 sources: {missing}")
    runs = [item["experiment_id"] for item in config["runs"]]
    run_metadata = {item["experiment_id"]: item for item in config["runs"]}

    primary, windows = build_primary_and_windows(
        sources["master"],
        runs,
        run_metadata,
        config["rules"]["initial_volume_l"],
        config["rules"]["yeast"],
        output,
    )
    temperature, transitions, temperature_summary, temperature_stats = process_temperature(
        sources["temperature"], windows, runs, output
    )
    active_end = temperature_stats["active_end"]
    windows["active_end"] = windows["experiment_id"].map(active_end)
    windows["active_end_rule"] = np.where(
        pd.to_datetime(windows["active_end"]).lt(windows["sampling_end"]),
        "controller_mode2_setpoint_9_to_11C",
        "last_primary_sample",
    )
    windows["postprocess_excluded_from_calibration"] = True
    windows.to_csv(output / "process_windows.csv", index=False)
    primary = finalize_primary_process_phase(primary, windows, output)

    minute, co2_summary, co2_stats = process_co2(
        sources["co2"],
        windows,
        runs,
        active_end,
        config["rules"]["invalid_raw_value"],
        set(config["rules"]["co2_model_excluded_runs"]),
        output,
    )
    events = process_operational_events(
        sources["master"],
        windows,
        runs,
        config["rules"]["lot3_operational_reference_map"],
        config["rules"]["initial_volume_l"],
        config["rules"]["organic_product_yan_mg_per_mg"],
        config["rules"]["dap_product_yan_mg_per_mg"],
        output,
    )
    condensates, gc_long, gc_stats = process_gc(
        sources["master"], sources["gc"], runs, windows, output
    )

    c_transitions = transitions[
        transitions["experiment_id"].eq("26159")
        & transitions["calibration_include"]
        & transitions["mode"].eq(3)
    ]
    protocol_c_setpoints = sorted(c_transitions["setpoint_c"].dropna().unique().tolist())
    profile_c_reached_21_in_calibration_window = 21.0 in protocol_c_setpoints

    checks = {}
    computed = {
        "fermentations": len(windows),
        "massview_rows": co2_stats["massview_total_rows"],
        "flow_rows_in_sampling_windows": co2_stats["flow_rows_in_sampling_windows"],
        "same_second_groups_before_sentinel_filter": co2_stats["same_second_groups_before_sentinel_filter"],
        "same_second_repetitions_absorbed_before_sentinel_filter": co2_stats["same_second_repetitions_absorbed_before_sentinel_filter"],
        "raw_2_621_observations": co2_stats["raw_2_621_observations_global"],
        "modbus_events_26136": co2_stats["modbus_events_26136"],
        "out_of_order_flow_timestamps": co2_stats["out_of_order_flow_timestamps"],
        "maximum_flow_gap_seconds": co2_stats["maximum_flow_gap_seconds"],
    }
    for name, expected in config["expected_totals"].items():
        checks[name] = {"expected": expected, "observed": computed[name], "pass": computed[name] == expected}

    qc_summary = {
        "campaign_id": config["campaign_id"],
        "verdict": "processed_qc_ready_conditional_for_calibration",
        "calibration_gate": [
            "independent Ultra verification of code, tables, and figures",
            "confirm MassView factory normal reference temperature and pressure",
            "reconcile nominal 80 mg/L YAN per pulse with 90.435 mg/L derived from 116 g organic + 46 g FDA at 230 L",
        ],
        "expected_value_checks": checks,
        "co2": co2_stats,
        "temperature": temperature_stats,
        "gc": gc_stats,
        "primary_results": {
            "rows": int(len(primary)),
            "postprocess_rows_excluded": int((~primary["calibration_include"]).sum()),
            "ethanol_negative_formula_artifacts_set_to_zero": int(
                primary["ethanol_correction"].ne("none").sum()
            ),
        },
        "operational_events": {
            "rows": int(len(events)),
            "runs": int(events["experiment_id"].nunique()),
            "reconstructed_lot3_rows": int(
                events["event_origin"].eq("reconstructed_from_lot1_relative_schedule").sum()
            ),
            "outside_active_process_excluded": int((~events["calibration_include"]).sum()),
            "parsed_nutrition_rows": int(
                events["dose_parse_status"].eq("parsed_springferm_organic_plus_fda").sum()
            ),
            "derived_yan_per_pulse_mg_l": sorted(
                events.loc[
                    events["dose_parse_status"].eq("parsed_springferm_organic_plus_fda"),
                    "yan_added_mg_l",
                ]
                .dropna()
                .round(6)
                .unique()
                .tolist()
            ),
            "nominal_protocol_yan_per_pulse_mg_l": config["rules"][
                "nominal_protocol_yan_per_pulse_mg_l"
            ],
            "yan_reconciliation_status": "pending_owner_confirmation_of_80_vs_90_435_mg_l",
        },
        "profile_c": {
            "confirmed_protocol": "16_to_18_to_21C",
            "setpoints_inside_calibration_window": protocol_c_setpoints,
            "reached_21C_inside_calibration_window": profile_c_reached_21_in_calibration_window,
            "interpretation": "The planned 21 C segment was not observed inside the primary-sample/calibration window and must not be treated as executed calibration input.",
        },
        "events_outside_sampling_window_excluded": int(
            (~events["inside_sampling_window"]).sum()
        ),
        "mapping_status": "confirmed",
        "massview_unit": "Ln/min",
        "massview_second_normalization_applied": False,
        "owner_decisions_applied": [
            "Lot 2 day/month swaps normalized to April",
            "39 pilot wine GC samples retained: 33 paired to MIX and 6 initial baselines",
            "MIX observations modeled as interval accumulations",
            "NQ and below-LOQ results retained as left-censored observations",
            "Lot 3 operational events reconstructed from Lot 1 relative schedules",
            "eight -0.3 percent v/v ethanol formula artifacts set to zero",
            "Lot 1 CO2 retained for QC but excluded from model calibration",
            "Sonda1 is the authoritative executed-temperature measurement",
        ],
    }
    (output / "qc_summary.json").write_text(
        json.dumps(qc_summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )

    make_figures(minute, temperature, gc_long, windows, runs, output)
    readme = """# Pilot 2026 data integration results\n\nThis directory is regenerated by `python fermentation_model/pilot_2026/run_data_integration.py`.\n\nThe sampling window is defined by the first and last primary-result sample. Controller mode 2 with a 9-11 °C band marks cooling/postprocess, and observations at or after that transition are excluded when it precedes the final primary sample. Sonda1 is the authoritative executed-temperature measurement.\n\nMassView values remain in Ln/min without a second normalization. Raw 2.621 values are removed before same-second averaging; same-second records are averaged. Leading/trailing zeros are valid, intermediate zeros are invalid and edge fills are explicitly artificial. Lot 1 CO2 is retained for visual QC but `co2_model_include=false` because the owner confirmed that its signals are unsuitable for calibration.\n\nThe GC layer contains 39 pilot wine MEF samples and 33 condensate MIX samples. The six extra wine samples are the initial baselines without condensate; every MIX is paired 1:1 with its contemporaneous wine sample. MIX dates for Lot 2 are normalized to April. Reported MIX values describe the diluted vial and are multiplied by 1000; wine results are not multiplied. NQ and below-LOQ results remain left-censored. Captured aroma mass is represented over each collection interval as MIX concentration times recovered A+B volume.\n\nLot 3 operational events are owner-confirmed reconstructions obtained by shifting the matching Lot 1 schedule to the Lot 3 process start. Temperature calibration must use the measured controller trace; second nutrient-pulse times remain identified as density-triggered schedule proxies. Eight negative ethanol formula artifacts were set to zero according to the owner's instruction, while their original values and correction flags remain preserved.\n\nThe dataset remains conditional for calibration until the owner reviews the executed loading-QC notebook and the MassView factory normal reference is documented.\n"""
    (output / "README.md").write_text(readme, encoding="utf-8")
    ultra_prompt = """# Independent Ultra verification prompt\n\nAct as an independent scientific-data and code auditor. Re-run `fermentation_model/pilot_2026/run_data_integration.py`, inspect `data_integration_config.json`, verify every expected count in `qc_summary.json`, and execute `pilot_2026/notebooks/pilot_2026_data_loading_qc.ipynb`. Confirm hashes and ZIP CRCs. Verify: nine sampling windows; Lot 2 GC dates in April; 39 wine samples, 33 MIX, six initial wine baselines and 33/33 pairing; ×1000 only for MIX; interval aroma mass using A+B volume; NQ/below-LOQ preserved as left-censored; eight -0.3 ethanol artifacts changed to zero with the originals retained; 15 Lot 3 events reconstructed from the matching Lot 1 relative schedules; Lot 1 CO2 retained for QC but excluded from calibration; Sonda1 as measured temperature; and exclusion of cooling or unexecuted setpoint stages. Do not modify files. Report PASS, PASS CONDITIONAL, or FAIL with reproducible evidence.\n"""
    (output / "ULTRA_VERIFICATION_PROMPT.md").write_text(ultra_prompt, encoding="utf-8")

    manifest = build_manifest(list(sources.values()), output)
    failures = [name for name, result in checks.items() if not result["pass"]]
    if failures:
        raise AssertionError(f"Expected-value checks failed: {failures}")
    print(json.dumps({"verdict": qc_summary["verdict"], "checks": checks, "outputs": len(manifest["outputs"])}, indent=2))


if __name__ == "__main__":
    main()
