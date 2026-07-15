from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
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
        if run_id in {"26157", "26158", "26159"}:
            if timestamp.year == 2026 and timestamp.day == 4 and 6 <= timestamp.month <= 12:
                timestamp = timestamp.replace(month=4, day=timestamp.month)
        return timestamp
    return pd.to_datetime(value, dayfirst=True, errors="coerce")


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
    master: Path, runs: list[str], run_metadata: dict[str, dict[str, object]], output: Path
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
    primary = primary.sort_values(["experiment_id", "timestamp", "Código muestra"])
    primary.to_csv(output / "primary_results_qc.csv", index=False)

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
    return primary, windows


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
                            "calibration_eligible": True,
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
                            "calibration_eligible": True,
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
    }
    return minute_qc, run_summary, stats


def process_operational_events(
    master: Path, windows: pd.DataFrame, runs: list[str], output: Path
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
    events["sampling_start"] = events["experiment_id"].map(lookup["sampling_start"])
    events["sampling_end"] = events["experiment_id"].map(lookup["sampling_end"])
    events["inside_sampling_window"] = events["timestamp"].between(
        events["sampling_start"], events["sampling_end"], inclusive="both"
    )
    events["calibration_include"] = events["inside_sampling_window"]
    events.to_csv(output / "operational_events_qc.csv", index=False)
    return events


def process_gc(master: Path, report: Path, runs: list[str], output: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    aromas = read_master_sheet(master, "09_Aromas_Cond")
    aromas = aromas[aromas["Código condensado / mix"].notna()].copy()
    aromas["mix_id"] = aromas["Código condensado / mix"].astype(str).str.strip()
    aromas["experiment_id"] = aromas["mix_id"].str.extract(r"^(26\d{3})", expand=False)
    aromas = aromas[aromas["experiment_id"].isin(runs)].copy()
    aromas["volume_a_ml"] = pd.to_numeric(aromas["Vol CA (mL)"], errors="coerce")
    aromas["volume_b_ml"] = pd.to_numeric(aromas["Vol CB (mL)"], errors="coerce")
    aromas["total_condensate_ml"] = aromas["volume_a_ml"].fillna(0) + aromas["volume_b_ml"].fillna(0)
    aromas["fraction_a"] = aromas["volume_a_ml"].fillna(0).div(aromas["total_condensate_ml"].replace(0, np.nan))
    aromas["fraction_b"] = aromas["volume_b_ml"].fillna(0).div(aromas["total_condensate_ml"].replace(0, np.nan))

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
    mix_samples = combined[
        combined["original_sample_id"].str.contains("-MIX-", na=False)
        & combined["experiment_id"].isin({"26157", "26158", "26159", "26210", "26211", "26212"})
    ].copy()
    long = mix_samples.melt(
        id_vars=["sample_number", "original_sample_id", "lab_sample_id", "experiment_id"],
        value_vars=analytes,
        var_name="analyte",
        value_name="reported_value",
    )
    long["loq_ug_l"] = long["analyte"].map(loq)
    long["vial_concentration_ug_l"] = long["reported_value"].map(parse_number)
    long["result_status"] = np.select(
        [
            long["reported_value"].astype(str).str.strip().str.upper().eq("NQ"),
            long["vial_concentration_ug_l"].notna()
            & long["vial_concentration_ug_l"].lt(long["loq_ug_l"]),
            long["vial_concentration_ug_l"].notna(),
        ],
        ["NQ", "below_loq", "quantified"],
        default="missing",
    )
    long["dilution_factor"] = 1000
    long["mix_concentration_ug_l"] = long["vial_concentration_ug_l"] * long["dilution_factor"]
    long["reported_basis"] = "diluted_GC_vial"
    long["preparation_stock_ml"] = 50
    long["volume_transferred_to_analyst_ml"] = 10
    long.to_csv(output / "gc_results_long_qc.csv", index=False)

    analyzed = aromas[aromas["experiment_id"].isin({"26157", "26158", "26159", "26210", "26211", "26212"})].copy()
    analyzed = analyzed.merge(
        mix_samples[["original_sample_id", "lab_sample_id"]],
        left_on="mix_id",
        right_on="original_sample_id",
        how="left",
        validate="one_to_one",
    )
    analyzed["gc_prepared"] = analyzed["lab_sample_id"].notna()
    analyzed["gc_scope"] = "confirmed_lots_2_and_3_only"
    analyzed.to_csv(output / "gc_condensate_mix_qc.csv", index=False)

    stats = {
        "gc_mix_samples": int(len(mix_samples)),
        "gc_analyte_results": int(len(long)),
        "gc_numeric_results": int(long["vial_concentration_ug_l"].notna().sum()),
        "gc_status_counts": {
            f"{analyte}:{status}": int(count)
            for (analyte, status), count in long.groupby(["analyte", "result_status"]).size().items()
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

    status = gc_long.groupby(["analyte", "result_status"]).size().unstack(fill_value=0)
    status = status.reindex(columns=["quantified", "below_loq", "NQ", "missing"], fill_value=0)
    fig, axis = plt.subplots(figsize=(11, 5.5))
    status.plot.bar(stacked=True, ax=axis, color=["#54a24b", "#eeca3b", "#e45756", "#bab0ac"])
    axis.set_ylabel("Results")
    axis.set_xlabel("")
    axis.set_title("GC condensate results by analytical status")
    axis.legend(title="Status", loc="upper left", bbox_to_anchor=(1.01, 1.0))
    axis.tick_params(axis="x", rotation=25)
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
    manifest = {
        "pipeline": "fermentation_model/pilot_2026/run_data_integration.py",
        "config": "fermentation_model/pilot_2026/data_integration_config.json",
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

    primary, windows = build_primary_and_windows(sources["master"], runs, run_metadata, output)
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

    minute, co2_summary, co2_stats = process_co2(
        sources["co2"], windows, runs, active_end, config["rules"]["invalid_raw_value"], output
    )
    events = process_operational_events(sources["master"], windows, runs, output)
    condensates, gc_long, gc_stats = process_gc(sources["master"], sources["gc"], runs, output)

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
        ],
        "expected_value_checks": checks,
        "co2": co2_stats,
        "temperature": temperature_stats,
        "gc": gc_stats,
        "profile_c": {
            "confirmed_protocol": "16_to_18_to_21C",
            "setpoints_inside_calibration_window": protocol_c_setpoints,
            "reached_21C_inside_calibration_window": profile_c_reached_21_in_calibration_window,
            "interpretation": "The planned 21 C segment was not observed inside the primary-sample/calibration window and must not be treated as executed calibration input.",
        },
        "events_outside_sampling_window_excluded": int((~events["calibration_include"]).sum()),
        "mapping_status": "confirmed",
        "massview_unit": "Ln/min",
        "massview_second_normalization_applied": False,
    }
    (output / "qc_summary.json").write_text(
        json.dumps(qc_summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )

    make_figures(minute, temperature, gc_long, windows, runs, output)
    readme = """# Pilot 2026 data integration results\n\nThis directory is regenerated by `python fermentation_model/pilot_2026/run_data_integration.py`.\n\nThe sampling window is defined by the first and last primary-result sample. MassView values are retained in Ln/min and are not normalized a second time. Raw 2.621 values are removed before same-second averaging. Leading/trailing zeros are retained, intermediate zeros are invalid, and edge fills are explicit artificial zeros. Controller mode 2 with a 9-11 °C band marks postprocess cooling; its event second is excluded from dynamic calibration. GC report values are diluted-vial concentrations and numeric values are multiplied by 1000 to recover MIX concentration. Lot 1 condensates are outside GC model scope.\n\nThe dataset remains conditional for calibration until the independent Ultra audit and confirmation of the MassView factory normal reference.\n"""
    (output / "README.md").write_text(readme, encoding="utf-8")
    ultra_prompt = """# Independent Ultra verification prompt\n\nAct as an independent scientific-data and code auditor. Re-run `fermentation_model/pilot_2026/run_data_integration.py`, inspect `data_integration_config.json`, verify every expected count in `qc_summary.json`, and visually review all four PNG figures. Confirm the Lot 2 day/month correction, nine primary-sample windows, raw 2.621 removal before same-second means, zero classification, explicit artificial edge fills, controller-event active-end logic, protocol C 16→18→21 handling, exclusion of postprocess/out-of-window events, and GC ×1000 correction from diluted-vial results. Confirm hashes and ZIP CRCs. Do not modify files. Report PASS, PASS CONDITIONAL, or FAIL with reproducible evidence.\n"""
    (output / "ULTRA_VERIFICATION_PROMPT.md").write_text(ultra_prompt, encoding="utf-8")

    manifest = build_manifest(list(sources.values()), output)
    failures = [name for name, result in checks.items() if not result["pass"]]
    if failures:
        raise AssertionError(f"Expected-value checks failed: {failures}")
    print(json.dumps({"verdict": qc_summary["verdict"], "checks": checks, "outputs": len(manifest["outputs"])}, indent=2))


if __name__ == "__main__":
    main()
