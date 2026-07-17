from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_MODEL_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


INTEGRATION_RESULTS = PILOT_DIR / "results" / "data_integration_2026"
ADAPTIVE_RESULTS = PILOT_DIR / "results" / "adaptive_design_2026"
DEFAULT_CONFIG_PATH = ADAPTIVE_DIR / "model_dataset_config.json"
INTEGRATION_CONFIG_PATH = PILOT_DIR / "data_integration_config.json"

SOURCE_FILES = {
    "integration_config": INTEGRATION_CONFIG_PATH,
    "integration_qc": INTEGRATION_RESULTS / "qc_summary.json",
    "process_windows": INTEGRATION_RESULTS / "process_windows.csv",
    "primary_results": INTEGRATION_RESULTS / "primary_results_qc.csv",
    "temperature_controller": INTEGRATION_RESULTS / "temperature_controller_qc.csv",
    "co2_minute": INTEGRATION_RESULTS / "co2_minute_qc.csv",
    "operational_events": INTEGRATION_RESULTS / "operational_events_qc.csv",
    "gc_results": INTEGRATION_RESULTS / "gc_results_long_qc.csv",
    "gc_pairs": INTEGRATION_RESULTS / "gc_mix_wine_pairs_long_qc.csv",
}

EXCLUDED_CO2_RUNS = frozenset({"26134", "26135", "26136"})


@dataclass(frozen=True)
class AdapterConfig:
    campaign_id: str = "pilot_2026"
    adapter_schema_version: int = 2
    random_seed: int = 20260716
    sugar_observation_policy: str = "components_preferred"
    co2_bin_minutes: int = 10
    lag1_clip_lower: float = 0.0
    lag1_clip_upper: float = 0.999
    normal_molar_volume_l_per_mol: float = 22.414
    massview_reference_status: str = (
        "factory_normalized_signal_used_as_reported_no_second_conversion"
    )
    required_integration_verdict: str = (
        "processed_qc_ready_conditional_for_calibration"
    )

    def validate(self) -> None:
        if self.sugar_observation_policy != "components_preferred":
            raise ValueError(
                "Only components_preferred is implemented because total and individual "
                "sugars must not be double counted."
            )
        if self.co2_bin_minutes < 1:
            raise ValueError("co2_bin_minutes must be at least one")
        if not 0.0 <= self.lag1_clip_lower <= self.lag1_clip_upper < 1.0:
            raise ValueError("Invalid lag-1 autocorrelation clipping interval")
        if self.normal_molar_volume_l_per_mol <= 0.0:
            raise ValueError("normal_molar_volume_l_per_mol must be positive")

    @classmethod
    def from_json(cls, path: Path = DEFAULT_CONFIG_PATH) -> "AdapterConfig":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        co2 = raw["co2_downsampling"]
        clip = co2["lag1_autocorrelation_clip"]
        massview = raw["massview"]
        result = cls(
            campaign_id=str(raw["campaign_id"]),
            adapter_schema_version=int(raw["adapter_schema_version"]),
            random_seed=int(raw["random_seed"]),
            sugar_observation_policy=str(raw["sugar_observation_policy"]),
            co2_bin_minutes=int(co2["bin_minutes"]),
            lag1_clip_lower=float(clip[0]),
            lag1_clip_upper=float(clip[1]),
            normal_molar_volume_l_per_mol=float(
                massview["normal_molar_volume_l_per_mol"]
            ),
            massview_reference_status=str(massview["reference_status"]),
            required_integration_verdict=str(raw["required_integration_verdict"]),
        )
        result.validate()
        return result


@dataclass
class ModelDataset:
    primary_observations: pd.DataFrame
    temperature_inputs: pd.DataFrame
    co2_observations: pd.DataFrame
    operational_events: pd.DataFrame
    wine_aroma_observations: pd.DataFrame
    condensate_interval_observations: pd.DataFrame
    run_metadata: pd.DataFrame
    carbon_balance_diagnostic: pd.DataFrame
    qc: dict[str, Any]


PRIMARY_SPECS = {
    "Densidad": ("density", "kg/m^3", "identity"),
    "Brix": ("brix", "degree_Brix", "identity"),
    "pH": ("pH", "dimensionless", "identity"),
    "DO": ("dissolved_oxygen", "mg/L", "identity"),
    "Oculyze Concentration": (
        "cell_concentration",
        "10^6_cells/mL",
        "identity",
    ),
    "Oculyze Viability": ("cell_viability", "percent", "identity"),
    "Oculyze Budding": ("cell_budding", "percent", "identity"),
    "Y15 PAN": ("primary_amino_nitrogen", "mg/L_N", "identity"),
    "Y15 Amonio": ("ammonium_nitrogen", "mg/L_N", "identity"),
    "Y15 Glicerol": ("glycerol", "g/L", "identity"),
    "Y15 Acetico": ("acetic_acid", "g/L", "identity"),
    "Y15 Piruvic": ("pyruvic_acid", "mg/L", "identity"),
    "Y15 Acetaldehido": ("acetaldehyde", "mg/L", "identity"),
    "Cf Alcolyzer Real (% v/v)": ("ethanol", "percent_v/v", "identity"),
}


def _read_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _run_id(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.0$", "", regex=True)


def _require_sources(sources: dict[str, Path]) -> None:
    missing = [str(path) for path in sources.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing authoritative integration outputs: {missing}")


def _load_tables(sources: dict[str, Path]) -> dict[str, Any]:
    _require_sources(sources)
    return {
        "integration_config": json.loads(
            sources["integration_config"].read_text(encoding="utf-8")
        ),
        "integration_qc": json.loads(
            sources["integration_qc"].read_text(encoding="utf-8")
        ),
        "windows": pd.read_csv(
            sources["process_windows"],
            parse_dates=["sampling_start", "sampling_end", "active_end"],
        ),
        "primary": pd.read_csv(
            sources["primary_results"],
            parse_dates=["timestamp", "sampling_start", "active_end"],
        ),
        "temperature": pd.read_csv(
            sources["temperature_controller"], parse_dates=["timestamp"]
        ),
        "co2": pd.read_csv(
            sources["co2_minute"], parse_dates=["minute"], low_memory=False
        ),
        "events": pd.read_csv(
            sources["operational_events"],
            parse_dates=[
                "timestamp",
                "sampling_start",
                "sampling_end",
                "active_end",
            ],
        ),
        "gc": pd.read_csv(sources["gc_results"], parse_dates=["sample_timestamp"]),
        "gc_pairs": pd.read_csv(
            sources["gc_pairs"],
            parse_dates=["timestamp", "capture_interval_start", "capture_interval_end"],
        ),
    }


def _primary_observations(primary: pd.DataFrame) -> pd.DataFrame:
    primary = primary.copy()
    primary["experiment_id"] = _run_id(primary["experiment_id"])
    primary["calibration_include"] = _read_bool(primary["calibration_include"])
    # The model-ready table is a kinetic input, not a QC mirror.  Cooling and
    # post-process samples remain available in the authoritative integration
    # outputs, but must never enter parameter estimation implicitly.
    primary = primary[primary["calibration_include"]].copy()
    records: list[dict[str, Any]] = []
    common = [
        "experiment_id",
        "timestamp",
        "time_h",
        "process_phase",
        "calibration_include",
        "lot",
        "reactor",
        "condition",
        "replicate",
        "protocol",
        "yeast",
        "initial_volume_l",
    ]
    for row in primary.itertuples(index=False, name=None):
        values = dict(zip(primary.columns, row))
        base = {name: values[name] for name in common}
        base["sample_id"] = values["Código muestra"]
        for source, (state, unit, operator) in PRIMARY_SPECS.items():
            value = pd.to_numeric(pd.Series([values[source]]), errors="coerce").iloc[0]
            if pd.isna(value):
                continue
            records.append(
                {
                    **base,
                    "state": state,
                    "observed_value": float(value),
                    "unit": unit,
                    "observation_operator": operator,
                    "source_column": source,
                    "censoring_type": "observed",
                    "lower_bound": np.nan,
                    "upper_bound": np.nan,
                }
            )

        component_sources = []
        for source, state in (
            ("Y15 Glucosa", "glucose"),
            ("Y15 Fructosa", "fructose"),
        ):
            value = pd.to_numeric(pd.Series([values[source]]), errors="coerce").iloc[0]
            if pd.notna(value):
                component_sources.append(source)
                records.append(
                    {
                        **base,
                        "state": state,
                        "observed_value": float(value),
                        "unit": "g/L",
                        "observation_operator": "identity",
                        "source_column": source,
                        "censoring_type": "observed",
                        "lower_bound": np.nan,
                        "upper_bound": np.nan,
                    }
                )
        total = pd.to_numeric(
            pd.Series([values["Y15 Glucosa-Fructosa"]]), errors="coerce"
        ).iloc[0]
        if not component_sources and pd.notna(total):
            records.append(
                {
                    **base,
                    "state": "total_sugar",
                    "observed_value": float(total),
                    "unit": "g/L",
                    "observation_operator": "glucose + fructose",
                    "source_column": "Y15 Glucosa-Fructosa",
                    "censoring_type": "observed",
                    "lower_bound": np.nan,
                    "upper_bound": np.nan,
                }
            )

    observations = pd.DataFrame.from_records(records)
    grouped = observations.groupby(["sample_id", "state"], dropna=False).size()
    if grouped.max() > 1:
        raise ValueError("Duplicate primary state observations were generated")
    sugar = observations[observations["state"].isin({"glucose", "fructose", "total_sugar"})]
    states_by_sample = sugar.groupby("sample_id")["state"].agg(set)
    double_counted = states_by_sample[
        states_by_sample.apply(
            lambda states: "total_sugar" in states
            and bool(states.intersection({"glucose", "fructose"}))
        )
    ]
    if not double_counted.empty:
        raise ValueError(f"Sugar observations double counted: {double_counted.index.tolist()}")
    return observations.sort_values(["experiment_id", "timestamp", "state"]).reset_index(
        drop=True
    )


def _temperature_inputs(
    temperature: pd.DataFrame, windows: pd.DataFrame
) -> pd.DataFrame:
    temperature = temperature.copy()
    windows = windows.copy()
    temperature["experiment_id"] = _run_id(temperature["experiment_id"])
    windows["experiment_id"] = _run_id(windows["experiment_id"])
    temperature["inside_sampling_window"] = _read_bool(
        temperature["inside_sampling_window"]
    )
    temperature["active_process"] = _read_bool(temperature["active_process"])
    result = temperature[
        temperature["inside_sampling_window"] & temperature["active_process"]
    ].copy()
    starts = windows.set_index("experiment_id")["sampling_start"]
    result["time_h"] = (
        result["timestamp"] - result["experiment_id"].map(starts)
    ).dt.total_seconds() / 3600.0
    result["process_phase"] = "active_process"
    result["kinetic_include"] = True
    result["executed_temperature_c"] = pd.to_numeric(
        result["sensor_1_c"], errors="coerce"
    )
    result["commanded_setpoint_c"] = pd.to_numeric(result["setpoint_c"], errors="coerce")
    result["executed_temperature_operator"] = "Sonda1 sensor_1_c"
    result["command_operator"] = "controller setpoint"
    keep = [
        "experiment_id",
        "timestamp",
        "time_h",
        "process_phase",
        "kinetic_include",
        "executed_temperature_c",
        "commanded_setpoint_c",
        "executed_temperature_operator",
        "command_operator",
        "controller_temperature_c",
        "sensor_2_c",
        "Modo_reading",
        "TempMin",
        "TempMax",
        "source_file",
        "source_row",
    ]
    return result[keep].sort_values(["experiment_id", "timestamp"]).reset_index(drop=True)


def _lag1_effective_sample_size(values: pd.Series, config: AdapterConfig) -> tuple[float, float]:
    array = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    n = len(array)
    if n < 3 or np.nanstd(array) <= 1e-12:
        return 0.0, float(max(n, 1))
    rho = float(np.corrcoef(array[:-1], array[1:])[0, 1])
    if not np.isfinite(rho):
        rho = 0.0
    rho = float(np.clip(rho, config.lag1_clip_lower, config.lag1_clip_upper))
    neff = float(np.clip(n * (1.0 - rho) / (1.0 + rho), 1.0, n))
    return rho, neff


def _co2_observations(co2: pd.DataFrame, config: AdapterConfig) -> pd.DataFrame:
    co2 = co2.copy()
    co2["experiment_id"] = _run_id(co2["experiment_id"])
    sampling_starts = co2.groupby("experiment_id")["minute"].min()
    co2["co2_model_include"] = _read_bool(co2["co2_model_include"])
    valid = co2[
        co2["co2_model_include"]
        & pd.to_numeric(co2["mean_flow_ln_min_calibration"], errors="coerce").notna()
    ].copy()
    reentered = set(valid["experiment_id"]).intersection(EXCLUDED_CO2_RUNS)
    if reentered:
        raise ValueError(f"Excluded Lot 1 CO2 runs re-entered calibration: {sorted(reentered)}")

    frames: list[pd.DataFrame] = []
    for run, group in valid.groupby("experiment_id", sort=True):
        group = group.sort_values("minute").copy()
        group["flow"] = pd.to_numeric(
            group["mean_flow_ln_min_calibration"], errors="coerce"
        )
        rho, neff_minutes = _lag1_effective_sample_size(group["flow"], config)
        # co2_minute_qc spans the authoritative primary-sample window, including
        # explicitly masked/artificial leading minutes.  Model time must start
        # at that window boundary, not at the first usable sensor value.
        start = pd.Timestamp(sampling_starts.loc[run])
        group["relative_minute"] = (
            group["minute"] - start
        ).dt.total_seconds() / 60.0
        group["bin_index"] = np.floor(
            group["relative_minute"] / config.co2_bin_minutes
        ).astype(int)
        aggregated = group.groupby("bin_index", sort=True).agg(
            timestamp=("minute", "first"),
            time_h=("relative_minute", lambda values: float(values.mean()) / 60.0),
            observed_flow_ln_min=("flow", "mean"),
            minutes_represented=("flow", "size"),
            calibration_seconds=("calibration_seconds", "sum"),
            valid_nonzero_seconds=("valid_nonzero", "sum"),
            valid_observed_edge_zero_seconds=("valid_observed_edge_zero", "sum"),
            unimputed_missing_seconds=("unimputed_missing_seconds", "sum"),
        )
        n_bins = len(aggregated)
        aggregated["experiment_id"] = run
        aggregated["unit"] = "Ln/min"
        aggregated["observation_operator"] = "MassView reported normalized flow"
        aggregated["downsampling_method"] = (
            "fixed_relative_time_bin_mean_with_lag1_effective_sample_size"
        )
        aggregated["downsampling_bin_minutes"] = config.co2_bin_minutes
        aggregated["source_valid_minutes"] = len(group)
        aggregated["lag1_autocorrelation"] = rho
        aggregated["effective_sample_size_minutes"] = neff_minutes
        aggregated["likelihood_weight"] = neff_minutes / max(n_bins, 1)
        aggregated["sigma_multiplier"] = math.sqrt(max(n_bins, 1) / neff_minutes)
        aggregated["weight_basis"] = (
            "provisional_signal_lag1_only; calibration_must_replace_with_residual_ESS"
        )
        frames.append(aggregated.reset_index(drop=True))
    if not frames:
        raise ValueError("No usable Pilot 2026 CO2 observations remain after masks")
    result = pd.concat(frames, ignore_index=True)
    return result.sort_values(["experiment_id", "timestamp"]).reset_index(drop=True)


def _density_trigger_intervals(events: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    events = events.copy()
    primary = primary.copy()
    events["experiment_id"] = _run_id(events["experiment_id"])
    primary["experiment_id"] = _run_id(primary["experiment_id"])
    events["timing_uncertain"] = events["model_input_note"].eq(
        "density_triggered_time_proxy"
    )
    events["timing_interval_start"] = pd.NaT
    events["timing_interval_end"] = pd.NaT
    events["timing_uncertainty_basis"] = "recorded_or_owner_confirmed_event_time"
    for index, event in events[events["timing_uncertain"]].iterrows():
        run = str(event["experiment_id"])
        samples = primary[primary["experiment_id"].eq(run)].copy()
        samples["density"] = pd.to_numeric(samples["Densidad"], errors="coerce")
        samples = samples.dropna(subset=["density", "timestamp"]).sort_values("timestamp")
        above = samples[samples["density"] > 1040.0]
        below = samples[samples["density"] <= 1040.0]
        if not above.empty and not below.empty:
            end_candidates = below[below["timestamp"] >= above["timestamp"].min()]
            if not end_candidates.empty:
                end = end_candidates["timestamp"].iloc[0]
                start_candidates = above[above["timestamp"] <= end]
                if not start_candidates.empty:
                    start = start_candidates["timestamp"].iloc[-1]
                    events.at[index, "timing_interval_start"] = start
                    events.at[index, "timing_interval_end"] = end
                    events.at[index, "timing_uncertainty_basis"] = (
                        "observed_primary_density_bracket_1040_kg_m3"
                    )
                    continue
        events.at[index, "timing_interval_start"] = event["sampling_start"]
        events.at[index, "timing_interval_end"] = event["active_end"]
        events.at[index, "timing_uncertainty_basis"] = (
            "unresolved_within_active_process_requires_latent_event_time_prior"
        )
    starts = pd.to_datetime(events["sampling_start"])
    events["timing_interval_start_h"] = (
        pd.to_datetime(events["timing_interval_start"]) - starts
    ).dt.total_seconds() / 3600.0
    events["timing_interval_end_h"] = (
        pd.to_datetime(events["timing_interval_end"]) - starts
    ).dt.total_seconds() / 3600.0
    events["calibration_include"] = _read_bool(events["calibration_include"])
    events = events[events["calibration_include"]].copy()
    return events.sort_values(["experiment_id", "timestamp"]).reset_index(drop=True)


def _wine_aroma_observations(gc: pd.DataFrame) -> pd.DataFrame:
    gc = gc.copy()
    gc["experiment_id"] = _run_id(gc["experiment_id"])
    wine = gc[gc["sample_type"].eq("wine_mef")].copy()
    wine["reported_numeric_value"] = pd.to_numeric(
        wine["sample_concentration_ug_l"], errors="coerce"
    )
    wine["observed_value"] = wine["reported_numeric_value"].where(
        wine["model_observation_type"].eq("observed"), np.nan
    )
    wine["lower_bound"] = pd.to_numeric(
        wine["censoring_lower_bound_ug_l"], errors="coerce"
    )
    wine["upper_bound"] = pd.to_numeric(
        wine["censoring_upper_bound_ug_l"], errors="coerce"
    )
    wine["unit"] = "ug/L_wine"
    wine["observation_operator"] = "instantaneous liquid concentration"
    keep = [
        "experiment_id",
        "sample_id",
        "lab_sample_id",
        "sample_timestamp",
        "time_h",
        "analyte",
        "observed_value",
        "reported_numeric_value",
        "unit",
        "observation_operator",
        "model_observation_type",
        "lower_bound",
        "upper_bound",
        "result_status",
        "sample_basis_loq_ug_l",
        "reported_value",
    ]
    return wine[keep].sort_values(
        ["experiment_id", "sample_timestamp", "analyte"]
    ).reset_index(drop=True)


def _condensate_interval_observations(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = pairs.copy()
    pairs["experiment_id"] = _run_id(pairs["experiment_id"])
    pairs["reported_numeric_value"] = pd.to_numeric(
        pairs["captured_mass_ug"], errors="coerce"
    )
    pairs["observed_value"] = pairs["reported_numeric_value"].where(
        pairs["mix_model_observation_type"].eq("observed"), np.nan
    )
    pairs["lower_bound"] = np.where(
        pairs["mix_model_observation_type"].eq("left_censored"), 0.0, np.nan
    )
    pairs["upper_bound"] = pd.to_numeric(
        pairs["captured_mass_upper_bound_ug"], errors="coerce"
    )
    pairs["unit"] = "ug_captured"
    pairs["observation_operator"] = (
        "integral(capture_rate, capture_interval_start, capture_interval_end)"
    )
    keep = [
        "experiment_id",
        "mix_id",
        "lab_sample_id",
        "analyte",
        "timestamp",
        "time_h",
        "capture_interval_start",
        "capture_interval_end",
        "capture_interval_h",
        "volume_a_ml",
        "volume_b_ml",
        "total_condensate_ml",
        "observed_value",
        "reported_numeric_value",
        "unit",
        "observation_operator",
        "mix_model_observation_type",
        "lower_bound",
        "upper_bound",
        "mix_result_status",
        "mix_concentration_ug_l",
        "mix_loq_ug_l",
        "paired_wine_sample_id",
    ]
    return pairs[keep].sort_values(
        ["experiment_id", "capture_interval_end", "analyte"]
    ).reset_index(drop=True)


def _run_metadata(
    windows: pd.DataFrame,
    primary: pd.DataFrame,
    integration_config: dict[str, Any],
    model_config_raw: dict[str, Any],
) -> pd.DataFrame:
    windows = windows.copy()
    windows["experiment_id"] = _run_id(windows["experiment_id"])
    primary = primary.copy()
    primary["experiment_id"] = _run_id(primary["experiment_id"])
    primary_metadata = (
        primary.groupby("experiment_id", sort=True)
        .agg(initial_volume_l=("initial_volume_l", "first"), yeast=("yeast", "first"))
        .reset_index()
    )
    result = windows.merge(primary_metadata, on="experiment_id", how="left", validate="1:1")
    defaults = model_config_raw["metadata_defaults"]
    result["campaign"] = integration_config["campaign_id"]
    result["storage_history"] = defaults["storage_history"]
    result["capture_system"] = defaults["capture_system"]
    result["co2_sensor"] = defaults["co2_sensor"]
    result["temperature_sensor"] = defaults["temperature_sensor"]
    result["temperature_command"] = defaults["temperature_command"]
    result["sensor_metadata_status"] = (
        "signal_identity_confirmed; factory-normalized Ln/min used as reported"
    )
    return result


def _carbon_balance(
    primary: pd.DataFrame,
    co2: pd.DataFrame,
    config: AdapterConfig,
) -> pd.DataFrame:
    """Compute a transparent, incomplete observed-carbon recovery diagnostic.

    This is not used as a likelihood term.  Biomass carbon, unmeasured products,
    sampling losses, and invalid/missing CO2 intervals are not closed, and the
    MassView normal reference is still unverified.
    """

    primary = primary.copy()
    co2 = co2.copy()
    primary["experiment_id"] = _run_id(primary["experiment_id"])
    co2["experiment_id"] = _run_id(co2["experiment_id"])
    primary["calibration_include"] = _read_bool(primary["calibration_include"])
    active = primary[primary["calibration_include"]].sort_values(
        ["experiment_id", "timestamp"]
    )
    rows: list[dict[str, Any]] = []

    def finite(row: pd.Series, column: str, scale: float = 1.0) -> float:
        value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
        return float(value) * scale if pd.notna(value) else 0.0

    def carbon(row: pd.Series) -> float:
        glucose = max(finite(row, "Y15 Glucosa"), 0.0) * 6.0 / 180.156
        fructose = max(finite(row, "Y15 Fructosa"), 0.0) * 6.0 / 180.156
        ethanol_g_l = max(finite(row, "Cf Alcolyzer Real (% v/v)"), 0.0) * 7.89
        ethanol = ethanol_g_l * 2.0 / 46.068
        glycerol = max(finite(row, "Y15 Glicerol"), 0.0) * 3.0 / 92.094
        acetate = max(finite(row, "Y15 Acetico"), 0.0) * 2.0 / 60.052
        pyruvate = max(finite(row, "Y15 Piruvic", 1e-3), 0.0) * 3.0 / 88.062
        acetaldehyde = max(finite(row, "Y15 Acetaldehido", 1e-3), 0.0) * 2.0 / 44.053
        return glucose + fructose + ethanol + glycerol + acetate + pyruvate + acetaldehyde

    for run, group in active.groupby("experiment_id", sort=True):
        first = group.iloc[0]
        last = group.iloc[-1]
        initial_carbon = carbon(first)
        final_liquid_carbon = carbon(last)
        flow = co2[
            co2["experiment_id"].eq(run)
            & _read_bool(co2["co2_model_include"])
        ].copy()
        flow_values = pd.to_numeric(flow["mean_flow_ln_min_calibration"], errors="coerce")
        co2_mol = float(flow_values.dropna().sum()) / config.normal_molar_volume_l_per_mol
        volume_l = float(first["initial_volume_l"])
        co2_carbon_mol_l = co2_mol / volume_l if volume_l > 0 else np.nan
        accounted = final_liquid_carbon + co2_carbon_mol_l
        recovery = accounted / initial_carbon if initial_carbon > 0 else np.nan
        rows.append(
            {
                "experiment_id": run,
                "initial_observed_carbon_mol_C_l": initial_carbon,
                "final_observed_liquid_carbon_mol_C_l": final_liquid_carbon,
                "co2_observed_carbon_mol_C_l": co2_carbon_mol_l,
                "partial_observed_carbon_recovery": recovery,
                "co2_model_scope": "calibration" if run not in EXCLUDED_CO2_RUNS else "QC_only",
                "massview_reference_status": config.massview_reference_status,
                "diagnostic_status": (
                    "not_evaluable" if run in EXCLUDED_CO2_RUNS else "diagnostic_only_incomplete"
                ),
                "limitations": (
                    "excludes biomass carbon, unmeasured products, sampling losses, and "
                    "masked/missing CO2 intervals; normal reference conditions unverified"
                ),
            }
        )
    return pd.DataFrame(rows)


def _quality_summary(dataset: ModelDataset, integration_qc: dict[str, Any]) -> dict[str, Any]:
    primary = dataset.primary_observations
    co2 = dataset.co2_observations
    events = dataset.operational_events
    wine = dataset.wine_aroma_observations
    condensate = dataset.condensate_interval_observations
    sugar_sets = (
        primary[primary["state"].isin({"glucose", "fructose", "total_sugar"})]
        .groupby("sample_id")["state"]
        .agg(set)
    )
    no_sugar_double_count = not sugar_sets.apply(
        lambda states: "total_sugar" in states
        and bool(states.intersection({"glucose", "fructose"}))
    ).any()
    excluded_reentered = sorted(set(co2["experiment_id"]).intersection(EXCLUDED_CO2_RUNS))
    censored_wine = wine["model_observation_type"].eq("left_censored")
    censored_condensate = condensate["mix_model_observation_type"].eq("left_censored")
    checks = {
        "nine_runs": dataset.run_metadata["experiment_id"].nunique() == 9,
        "cooling_formally_excluded_from_kinetics": bool(
            dataset.temperature_inputs["process_phase"].eq("active_process").all()
            and dataset.temperature_inputs["kinetic_include"].all()
            and primary["calibration_include"].all()
            and events["calibration_include"].all()
        ),
        "sugar_not_double_counted": no_sugar_double_count,
        "excluded_lot1_co2_absent": not excluded_reentered,
        "co2_is_downsampled": len(co2) < integration_qc["co2"]["flow_rows_in_sampling_windows"],
        "co2_weights_positive": bool((co2["likelihood_weight"] > 0.0).all()),
        "wine_censoring_bounds_complete": bool(
            wine.loc[censored_wine, "upper_bound"].notna().all()
        ),
        "condensate_censoring_bounds_complete": bool(
            condensate.loc[censored_condensate, "upper_bound"].notna().all()
        ),
        "condensate_intervals_positive": bool(
            (condensate["capture_interval_h"] > 0.0).all()
        ),
        "density_trigger_uncertainty_explicit": bool(
            events.loc[events["timing_uncertain"], "timing_interval_start"].notna().all()
            and events.loc[events["timing_uncertain"], "timing_interval_end"].notna().all()
        ),
        "carbon_diagnostic_present": len(dataset.carbon_balance_diagnostic) == 9,
    }
    return {
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "excluded_co2_runs_reentered": excluded_reentered,
        "row_counts": {
            "primary_observations": len(primary),
            "temperature_inputs": len(dataset.temperature_inputs),
            "co2_observations": len(co2),
            "operational_events": len(events),
            "wine_aroma_observations": len(wine),
            "condensate_interval_observations": len(condensate),
        },
        "scientific_limitations": [
            "MassView factory reference temperature and pressure do not affect the "
            "reported Ln/min signal; no second normalization is applied.",
            "The 90.435 mg/L mass-derived YAN value remains auditable in integration "
            "tables; calibration uses the owner-confirmed fixed 80 mg/L per historical pulse.",
            "Storage history is not present in the authoritative integration outputs.",
            "Carbon recovery is diagnostic-only and is not a closed elemental balance.",
        ],
    }


def load_model_dataset(
    config: AdapterConfig | None = None,
    *,
    sources: dict[str, Path] | None = None,
    config_path: Path = DEFAULT_CONFIG_PATH,
) -> ModelDataset:
    config = config or AdapterConfig.from_json(config_path)
    config.validate()
    source_paths = dict(SOURCE_FILES if sources is None else sources)
    tables = _load_tables(source_paths)
    verdict = tables["integration_qc"].get("verdict")
    if verdict != config.required_integration_verdict:
        raise ValueError(
            f"Unexpected integration verdict {verdict!r}; expected "
            f"{config.required_integration_verdict!r}"
        )
    model_config_raw = json.loads(Path(config_path).read_text(encoding="utf-8"))
    dataset = ModelDataset(
        primary_observations=_primary_observations(tables["primary"]),
        temperature_inputs=_temperature_inputs(tables["temperature"], tables["windows"]),
        co2_observations=_co2_observations(tables["co2"], config),
        operational_events=_density_trigger_intervals(tables["events"], tables["primary"]),
        wine_aroma_observations=_wine_aroma_observations(tables["gc"]),
        condensate_interval_observations=_condensate_interval_observations(
            tables["gc_pairs"]
        ),
        run_metadata=_run_metadata(
            tables["windows"],
            tables["primary"],
            tables["integration_config"],
            model_config_raw,
        ),
        carbon_balance_diagnostic=_carbon_balance(tables["primary"], tables["co2"], config),
        qc={},
    )
    dataset.qc = _quality_summary(dataset, tables["integration_qc"])
    if dataset.qc["verdict"] != "PASS":
        raise ValueError(f"Model-ready dataset QC failed: {dataset.qc['checks']}")
    return dataset


def write_model_dataset(
    dataset: ModelDataset,
    *,
    config: AdapterConfig,
    config_path: Path,
    result_root: Path = ADAPTIVE_RESULTS,
    sources: dict[str, Path] | None = None,
) -> Path:
    source_paths = dict(SOURCE_FILES if sources is None else sources)
    config_payload = asdict(config)
    run_dir = create_immutable_run_directory(result_root, "model_dataset", config_payload)
    outputs: list[Path] = []
    table_map = {
        "primary_observations.csv": dataset.primary_observations,
        "temperature_inputs.csv": dataset.temperature_inputs,
        "co2_observations.csv": dataset.co2_observations,
        "operational_events.csv": dataset.operational_events,
        "wine_aroma_observations.csv": dataset.wine_aroma_observations,
        "condensate_observations.csv": dataset.condensate_interval_observations,
        "run_metadata.csv": dataset.run_metadata,
        "carbon_balance_diagnostic.csv": dataset.carbon_balance_diagnostic,
    }
    code_paths = [Path(__file__), ADAPTIVE_DIR / "run_artifacts.py", Path(config_path)]
    try:
        for name, table in table_map.items():
            path = run_dir / name
            table.to_csv(path, index=False)
            outputs.append(path)
        copied_config = run_dir / "adapter_config.json"
        write_json(copied_config, config_payload)
        outputs.append(copied_config)
        qc_path = run_dir / "dataset_qc.json"
        write_json(qc_path, dataset.qc)
        outputs.append(qc_path)
        manifest = build_manifest(
            run_dir=run_dir,
            stage="model_dataset",
            config=config_payload,
            sources=source_paths,
            code_paths=code_paths,
            random_seeds=[config.random_seed],
            status="completed",
            convergence={"solver": "not_applicable", "adapter_qc": dataset.qc["verdict"]},
            gate={
                "name": "phase_A_model_ready_adapter",
                "verdict": "PASS",
                "profiles_for_physical_execution": False,
            },
            outputs=outputs,
        )
        write_json(run_dir / "run_manifest.json", manifest)
    except Exception as error:
        failure = build_manifest(
            run_dir=run_dir,
            stage="model_dataset",
            config=config_payload,
            sources=source_paths,
            code_paths=code_paths,
            random_seeds=[config.random_seed],
            status="failed",
            convergence={
                "solver": "not_applicable",
                "exception_type": type(error).__name__,
                "exception": str(error),
            },
            gate={
                "name": "phase_A_model_ready_adapter",
                "verdict": "FAIL",
                "profiles_for_physical_execution": False,
            },
            outputs=outputs,
        )
        write_json(run_dir / "run_manifest.json", failure)
        raise
    return run_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the Pilot 2026 model-ready dataset")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--result-root", type=Path, default=ADAPTIVE_RESULTS)
    args = parser.parse_args()
    config = AdapterConfig.from_json(args.config)
    dataset = load_model_dataset(config, config_path=args.config)
    run_dir = write_model_dataset(
        dataset,
        config=config,
        config_path=args.config,
        result_root=args.result_root,
    )
    print(
        json.dumps(
            {
                "verdict": dataset.qc["verdict"],
                "run_directory": str(run_dir),
                "row_counts": dataset.qc["row_counts"],
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
