from __future__ import annotations

"""Recalibrate Pilot-2026 aromas with physical gas-liquid transfer.

The corrected loss operator converts emitted CO2 mass rate into gas volumetric
turnover with the ideal-gas density and applies an NTU transfer efficiency before
the equilibrium UNIFAC partition.  It is compared against the previous empirical
mass-rate scale using the same released-rCO2 forcing, observations and errors.
"""

import argparse
import copy
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from nbclient import NotebookClient


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_DIR = SCRIPT_DIR.parent
ROOT_DIR = FERMENTATION_DIR.parent
for path in (FERMENTATION_DIR, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from laboratory_2026 import run_co2_matrix_cross_validation_2026 as lab_co2  # noqa: E402
from pilot_2026 import run_co2_solubility_cross_lot_validation_2026 as pilot_co2  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_calibration  # noqa: E402


MODEL_DATASET_DIR = pilot_co2.MODEL_DATASET_DIR
PRIMARY_CALIBRATION_DIR = pilot_co2.PILOT_CALIBRATION_DIR
CO2_RESULTS_DIR = pilot_co2.RESULTS_DIR
AROMA_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "aroma_calibration_config.json"
PRIMARY_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "calibration_config.json"
HISTORICAL_AROMA_DIR = (
    SCRIPT_DIR
    / "results"
    / "adaptive_design_2026"
    / "aroma_calibration"
    / "20260717T161527Z_e927a5"
)
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_dynamic_transfer_recalibration_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_dynamic_transfer_recalibration.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_dynamic_transfer_recalibration.executed.ipynb"
)

RUNS = tuple(pilot_co2.VALID_RUNS)
HOLDOUTS = frozenset(pilot_co2.HOLDOUTS.values())
CALIBRATION_RUNS = frozenset(run for run in RUNS if run not in HOLDOUTS)
SOURCE_LABELS = {
    "released_co2_model": "Escala empírica anterior",
    "dynamic_transfer_model": "Transferencia gas-líquido corregida",
}
SPECIES_LABELS = {
    "ethyl_acetate": "Acetato de etilo",
    "ethyl_octanoate": "Octanoato de etilo",
    "isoamyl_acetate": "Acetato de isoamilo",
}
COLORS = {
    "released_co2_model": "#326FA8",
    "dynamic_transfer_model": "#D97732",
    "observed": "#27313B",
    "holdout": "#B45A7A",
    "grid": "#D8DEE5",
}
RANDOM_SEED = 20260820


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)!r}")


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default) + "\n",
        encoding="utf-8",
    )


def _load_co2_parameters() -> dict[str, dict[str, float]]:
    table = pd.read_csv(CO2_RESULTS_DIR / "fit_parameters.csv")
    fits: dict[str, dict[str, float]] = {}
    for matrix, group in table.groupby("calibration_matrix"):
        fits[str(matrix)] = dict(
            zip(group["parameter"].astype(str), group["estimate"].astype(float))
        )
    missing = set(pilot_co2.MATRICES) - set(fits)
    if missing:
        raise RuntimeError(f"Missing validated CO2 parameter sets: {sorted(missing)}")
    return fits


def _build_released_co2_overrides() -> tuple[
    dict[str, tuple[np.ndarray, np.ndarray]], pd.DataFrame, pd.DataFrame
]:
    co2_parameters = _load_co2_parameters()
    batches, theta, tables = pilot_co2._pilot_upstream_context()
    pulses = pilot_co2.build_pulse_schedule(tables.primary, tables.events)
    batches = pilot_co2.override_pilot_pulses(batches, pulses)
    sampling = pilot_co2.load_sampling_schedule(tables.primary)
    normalization = pilot_co2.load_sensor_normalization()
    observations, _, _, _ = pilot_co2.load_co2_observations(
        tables.metadata,
        tables.temperature,
        sampling,
        pulses,
        normalization,
    )
    theta_by_matrix = {matrix: dict(theta) for matrix in pilot_co2.MATRICES}
    cache, _ = lab_co2.build_driver_cache(batches, theta_by_matrix, observations)

    overrides: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    rows: list[dict[str, Any]] = []
    for run in RUNS:
        matrix = pilot_co2.LOT_BY_RUN[run]
        params = co2_parameters[matrix]
        qgas, qprod, o2, phi_ana = lab_co2.raw_qgas_grid_prediction(
            cache[run],
            params["kCO2_release_h"],
            params["CO2sat_scale"],
            params["O2_qmax_mg_gdw_h"],
            params["O2_initial_scale"],
            chemistry_aligned=True,
            nitrogen_boost_transition=True,
            pulse_t_rise_h=params["pulse_t_rise_h"],
            pulse_activity_gain=params["pulse_activity_gain"],
            continuous_release=True,
            bounded_chemical_activation=True,
            chem_activation_start_fraction=params["chem_activation_start_fraction"],
            chem_activation_duration_fraction=params["chem_activation_duration_fraction"],
        )
        released = np.maximum(float(params["matrix_gain"]) * qgas, 0.0)
        overrides[run] = (cache[run].time_h.copy(), released)
        for time_h, rate, produced, oxygen, anaerobic in zip(
            cache[run].time_h, released, qprod, o2, phi_ana
        ):
            rows.append(
                {
                    "experiment_id": run,
                    "matrix": matrix,
                    "time_h": float(time_h),
                    "released_rco2_g_l_h": float(rate),
                    "internal_qprod_g_l_h": float(produced),
                    "dissolved_o2_mg_l": float(oxygen),
                    "anaerobic_fraction": float(anaerobic),
                    "pulse_time_h": float(cache[run].n_pulse_time_h),
                    "co2_fit_role": "holdout" if run in HOLDOUTS else "calibration",
                }
            )
    return overrides, pd.DataFrame(rows), pulses


def _build_aroma_forcings(
    tables: pilot_calibration.CalibrationTables,
    config: dict[str, Any],
    primary_config: dict[str, Any],
    partitions: dict[str, dict[str, float]],
    released_overrides: dict[str, tuple[np.ndarray, np.ndarray]],
) -> tuple[
    dict[str, dict[str, dict[str, aroma.AromaForcing]]],
    dict[str, pd.DataFrame],
    dict[str, pd.DataFrame],
    pd.DataFrame,
]:
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]] = {
        source: {} for source in SOURCE_LABELS
    }
    wine_tables: dict[str, pd.DataFrame] = {}
    condensate_tables: dict[str, pd.DataFrame] = {}
    forcing_rows: list[dict[str, Any]] = []
    for species, analyte_label in config["priority_analytes"].items():
        empirical, wine, condensate, _ = aroma.build_forcings(
            tables,
            primary_config,
            PRIMARY_CALIBRATION_DIR,
            species,
            analyte_label,
            partitions[species],
            co2_rate_override_by_run=released_overrides,
        )
        dynamic, _, _, _ = aroma.build_forcings(
            tables,
            primary_config,
            PRIMARY_CALIBRATION_DIR,
            species,
            analyte_label,
            partitions[species],
            co2_rate_override_by_run=released_overrides,
            loss_model=aroma.LOSS_MODEL_DYNAMIC_TRANSFER,
        )
        forcing_sets["released_co2_model"][species] = empirical
        forcing_sets["dynamic_transfer_model"][species] = dynamic
        wine_tables[species] = wine
        condensate_tables[species] = condensate
        for source, forcings in (
            ("released_co2_model", empirical),
            ("dynamic_transfer_model", dynamic),
        ):
            for run, forcing in forcings.items():
                for time_h, rco2, loss_basis, qgas, partition_k, temp_c, ethanol, sugar in zip(
                    forcing.time_h,
                    forcing.co2_rate_g_l_h,
                    forcing.loss_basis_h_inv,
                    forcing.gas_turnover_h_inv,
                    forcing.partition_basis_l_g,
                    forcing.temperature_c,
                    forcing.ethanol_g_l,
                    forcing.total_sugar_g_l,
                ):
                    forcing_rows.append(
                        {
                            "source": source,
                            "species": species,
                            "experiment_id": run,
                            "matrix": pilot_co2.LOT_BY_RUN[run],
                            "time_h": float(time_h),
                            "rco2_g_l_h": float(rco2),
                            "loss_basis": float(loss_basis),
                            "gas_turnover_h_inv": float(qgas),
                            "partition_k_gas_over_liquid": float(partition_k),
                            "temperature_c": float(temp_c),
                            "ethanol_g_l": float(ethanol),
                            "total_sugar_g_l": float(sugar),
                            "log_k_temperature_term": float(
                                partitions[species]["temp_slope"] * (temp_c - 20.0)
                            ),
                            "log_k_ethanol_term": float(
                                partitions[species]["ethanol_slope"] * (ethanol - 50.0)
                            ),
                            "log_k_sugar_term": float(
                                partitions[species]["sugar_slope"] * (sugar - 100.0)
                            ),
                            "role": "holdout" if run in HOLDOUTS else "calibration",
                        }
                    )
    return forcing_sets, wine_tables, condensate_tables, pd.DataFrame(forcing_rows)


def _subset_inputs(
    forcings: dict[str, aroma.AromaForcing],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    runs: frozenset[str],
) -> tuple[dict[str, aroma.AromaForcing], pd.DataFrame, pd.DataFrame]:
    subset_forcings = {run: value for run, value in forcings.items() if run in runs}
    subset_wine = wine[wine["experiment_id"].astype(str).isin(runs)].copy()
    subset_condensate = condensate[
        condensate["experiment_id"].astype(str).isin(runs)
    ].copy()
    return subset_forcings, subset_wine, subset_condensate


def _fit_all(
    config: dict[str, Any],
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]],
    wine_tables: dict[str, pd.DataFrame],
    condensate_tables: dict[str, pd.DataFrame],
) -> dict[tuple[str, str, str], aroma.AromaFit]:
    fits: dict[tuple[str, str, str], aroma.AromaFit] = {}
    for source in SOURCE_LABELS:
        for fit_scope, run_set in (
            ("calibration_only", CALIBRATION_RUNS),
            ("all_data", frozenset(RUNS)),
        ):
            fit_config = copy.deepcopy(config)
            if fit_scope == "calibration_only":
                fit_config["optimization"]["sobol_multistarts"] = 8
                fit_config["optimization"]["profile_grid_log_offsets"] = [0.0]
            elif source == "released_co2_model":
                fit_config["optimization"]["sobol_multistarts"] = 8
                fit_config["optimization"]["profile_grid_log_offsets"] = [0.0]
            fit_config["optimization"]["seed"] = RANDOM_SEED + (
                1000 if source == "dynamic_transfer_model" else 2000
            ) + (0 if fit_scope == "calibration_only" else 100)
            for species in config["priority_analytes"]:
                forcings, wine, condensate = _subset_inputs(
                    forcing_sets[source][species],
                    wine_tables[species],
                    condensate_tables[species],
                    run_set,
                )
                fits[(source, fit_scope, species)] = aroma.fit_species(
                    species, forcings, wine, condensate, fit_config
                )
    return fits


def _collect_parameters(
    fits: dict[tuple[str, str, str], aroma.AromaFit]
) -> pd.DataFrame:
    rows = []
    for (source, fit_scope, species), fit in fits.items():
        table = fit.parameter_table.copy()
        table.insert(0, "fit_scope", fit_scope)
        table.insert(0, "forcing_source", source)
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def _collect_fit_diagnostics(
    fits: dict[tuple[str, str, str], aroma.AromaFit]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    validation_rows = []
    multistarts = []
    profiles = []
    covariance = []
    for (source, fit_scope, species), fit in fits.items():
        validation_rows.append(
            {
                "forcing_source": source,
                "fit_scope": fit_scope,
                "species": species,
                "objective": fit.validation["objective"],
                "objective_per_observation": fit.validation[
                    "objective_per_observation"
                ],
                "active_bound_fraction": fit.validation["active_bound_fraction"],
                "loss_separately_identified": fit.validation[
                    "loss_separately_identified"
                ],
                "weak_parameters": ";".join(fit.validation["weak_parameters"]),
            }
        )
        multistarts.append(
            fit.multistart_summary.assign(
                forcing_source=source, fit_scope=fit_scope
            )
        )
        profiles.append(
            fit.profiles.assign(forcing_source=source, fit_scope=fit_scope)
        )
        covariance.append(
            fit.covariance.rename_axis("row_parameter")
            .reset_index()
            .assign(species=species, forcing_source=source, fit_scope=fit_scope)
        )
    return (
        pd.DataFrame(validation_rows),
        pd.concat(multistarts, ignore_index=True),
        pd.concat(profiles, ignore_index=True),
        pd.concat(covariance, ignore_index=True),
    )


def _prediction_inventory(
    fits: dict[tuple[str, str, str], aroma.AromaFit],
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]],
    wine_tables: dict[str, pd.DataFrame],
    condensate_tables: dict[str, pd.DataFrame],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    wine_frames = []
    condensate_frames = []
    for (source, fit_scope, species), fit in fits.items():
        wine_pred, cond_pred = aroma.prediction_tables(
            forcing_sets[source][species],
            wine_tables[species],
            condensate_tables[species],
            fit.log_values,
            config,
        )
        first_time = wine_pred.groupby("experiment_id")["time_h"].transform("min")
        wine_pred["initialization_point"] = np.isclose(wine_pred["time_h"], first_time)
        wine_pred["role"] = wine_pred["experiment_id"].map(
            lambda run: "holdout" if str(run) in HOLDOUTS else "calibration"
        )
        cond_pred["role"] = cond_pred["experiment_id"].map(
            lambda run: "holdout" if str(run) in HOLDOUTS else "calibration"
        )
        for frame in (wine_pred, cond_pred):
            frame.insert(0, "species", species)
            frame.insert(0, "fit_scope", fit_scope)
            frame.insert(0, "forcing_source", source)
        wine_frames.append(wine_pred)
        condensate_frames.append(cond_pred)
    return pd.concat(wine_frames, ignore_index=True), pd.concat(
        condensate_frames, ignore_index=True
    )


def _transfer_diagnostic_inventory(
    fits: dict[tuple[str, str, str], aroma.AromaFit],
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]],
) -> pd.DataFrame:
    rows = []
    for source in SOURCE_LABELS:
        for species in SPECIES_LABELS:
            fit = fits[(source, "all_data", species)]
            for run, forcing in forcing_sets[source][species].items():
                diagnostics = aroma.transfer_diagnostics(forcing, fit.log_values)
                diagnostics.insert(0, "role", "holdout" if run in HOLDOUTS else "calibration")
                diagnostics.insert(0, "matrix", pilot_co2.LOT_BY_RUN[run])
                diagnostics.insert(0, "experiment_id", run)
                diagnostics.insert(0, "species", species)
                diagnostics.insert(0, "forcing_source", source)
                rows.append(diagnostics)
    return pd.concat(rows, ignore_index=True)


def _metric_row(
    frame: pd.DataFrame,
    prediction_column: str,
    *,
    domain: str,
) -> dict[str, Any]:
    observed = frame[frame["status"].eq("observed")].copy()
    if domain == "wine":
        observed = observed[~observed["initialization_point"]]
    error = (
        observed[prediction_column].to_numpy(dtype=float)
        - observed["observed_or_upper_bound"].to_numpy(dtype=float)
    )
    observed_values = observed["observed_or_upper_bound"].to_numpy(dtype=float)
    observed_mean = float(np.mean(observed_values)) if len(observed_values) else np.nan
    rmse = float(np.sqrt(np.mean(error**2))) if len(error) else np.nan
    censored = frame[~frame["status"].eq("observed")].copy()
    return {
        "n_observed_scored": int(len(observed)),
        "observed_mean": observed_mean,
        "rmse": rmse,
        "nrmse_over_mean": rmse / observed_mean
        if np.isfinite(observed_mean) and observed_mean > 0.0
        else np.nan,
        "mae": float(np.mean(np.abs(error))) if len(error) else np.nan,
        "bias": float(np.mean(error)) if len(error) else np.nan,
        "n_censored": int(len(censored)),
        "censor_compliance_fraction": float(
            (censored[prediction_column] <= censored["observed_or_upper_bound"]).mean()
        )
        if len(censored)
        else np.nan,
    }


def _score_predictions(
    wine_predictions: pd.DataFrame, condensate_predictions: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    group_columns = ["forcing_source", "fit_scope", "species", "role"]
    for keys, frame in wine_predictions.groupby(group_columns, sort=True):
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "domain": "wine",
                **_metric_row(frame, "predicted_ug_l", domain="wine"),
            }
        )
    for keys, frame in condensate_predictions.groupby(group_columns, sort=True):
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "domain": "condensate",
                **_metric_row(
                    frame, "predicted_captured_ug", domain="condensate"
                ),
            }
        )
    return pd.DataFrame(rows)


def _forcing_exposure(forcing_grid: pd.DataFrame) -> pd.DataFrame:
    representative = forcing_grid[
        forcing_grid["species"].eq("ethyl_acetate")
        & forcing_grid["source"].eq("dynamic_transfer_model")
    ]
    rows = []
    for run, group in representative.groupby("experiment_id", sort=True):
        group = group.sort_values("time_h")
        rows.append(
            {
                "experiment_id": run,
                "matrix": group["matrix"].iloc[0],
                "role": group["role"].iloc[0],
                "integrated_rco2_g_l": float(
                    np.trapz(group["rco2_g_l_h"], group["time_h"])
                ),
                "peak_rco2_g_l_h": float(group["rco2_g_l_h"].max()),
                "integrated_gas_turnover": float(
                    np.trapz(group["gas_turnover_h_inv"], group["time_h"])
                ),
                "peak_gas_turnover_h_inv": float(group["gas_turnover_h_inv"].max()),
            }
        )
    return pd.DataFrame(rows)


def _style_axis(ax: plt.Axes) -> None:
    ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)


def _save_figure(fig: plt.Figure, name: str, *, top: float = 0.95) -> Path:
    path = FIGURE_DIR / name
    fig.tight_layout(rect=(0.0, 0.0, 1.0, top))
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_rco2_forcing(forcing_grid: pd.DataFrame, pulses: pd.DataFrame) -> Path:
    representative = forcing_grid[
        forcing_grid["species"].eq("ethyl_acetate")
        & forcing_grid["source"].eq("dynamic_transfer_model")
    ]
    pulse_lookup = pulses.set_index("batch")["pulse_time_h"].to_dict()
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    right_axes = []
    for panel_index, (ax, run) in enumerate(zip(axes.flat, RUNS)):
        group = representative[representative["experiment_id"].eq(run)].sort_values("time_h")
        ax.plot(
            group["time_h"],
            group["rco2_g_l_h"],
            color=COLORS["released_co2_model"],
            lw=2.0,
            label="rCO₂ másico",
        )
        ax_right = ax.twinx()
        right_axes.append(ax_right)
        ax_right.plot(
            group["time_h"],
            group["gas_turnover_h_inv"],
            color=COLORS["dynamic_transfer_model"],
            lw=1.5,
            label="Qgas/VL",
        )
        ax_right.spines["top"].set_visible(False)
        if panel_index % 3 == 2:
            ax_right.set_ylabel("Qgas/VL (h⁻¹)")
        pulse = pulse_lookup.get(run, np.nan)
        if np.isfinite(pulse):
            ax.axvline(pulse, color="#B45A7A", ls="--", lw=1.2)
        ax.set_title(f"{run} · {pilot_co2.LOT_BY_RUN[run]}")
        ax.set_xlabel("Tiempo (h)")
        _style_axis(ax)
    axes[0, 0].set_ylabel("rCO2 (g L⁻¹ h⁻¹)")
    axes[1, 0].set_ylabel("rCO2 (g L⁻¹ h⁻¹)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    right_handles, right_labels = right_axes[0].get_legend_handles_labels()
    fig.legend(
        handles + right_handles,
        labels + right_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Conversión de rCO₂ másico a recambio volumétrico", y=0.995, fontsize=14)
    return _save_figure(fig, "gas_flow_conversion.png", top=0.89)


def plot_exposure(exposure: pd.DataFrame) -> Path:
    view = exposure.set_index("experiment_id").reindex(RUNS)
    x = np.arange(len(view))
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.bar(x, view["integrated_gas_turnover"], color=COLORS["dynamic_transfer_model"])
    ax.set_xticks(x, view.index)
    ax.set_ylabel("Recambios volumétricos acumulados (L gas/L líquido)")
    ax.set_xlabel("Fermentación")
    ax.set_title("Exposición acumulada al flujo gaseoso corregido")
    _style_axis(ax)
    return _save_figure(fig, "integrated_rco2_exposure.png")


def plot_holdout_wine(wine_predictions: pd.DataFrame) -> Path:
    view = wine_predictions[
        wine_predictions["fit_scope"].eq("calibration_only")
        & wine_predictions["role"].eq("holdout")
    ]
    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=False)
    for row_index, species in enumerate(SPECIES_LABELS):
        for col_index, run in enumerate(sorted(HOLDOUTS)):
            ax = axes[row_index, col_index]
            for source in SOURCE_LABELS:
                group = view[
                    view["species"].eq(species)
                    & view["experiment_id"].eq(run)
                    & view["forcing_source"].eq(source)
                ].sort_values("time_h")
                ax.plot(
                    group["time_h"],
                    group["predicted_ug_l"],
                    color=COLORS[source],
                    lw=2.0,
                    label=SOURCE_LABELS[source],
                )
            observed = view[
                view["species"].eq(species)
                & view["experiment_id"].eq(run)
                & view["forcing_source"].eq("released_co2_model")
            ]
            direct = observed[observed["status"].eq("observed")]
            censored = observed[~observed["status"].eq("observed")]
            ax.scatter(
                direct["time_h"],
                direct["observed_or_upper_bound"],
                color=COLORS["observed"],
                s=25,
                zorder=3,
                label="Observado",
            )
            if not censored.empty:
                ax.scatter(
                    censored["time_h"],
                    censored["observed_or_upper_bound"],
                    facecolors="none",
                    edgecolors=COLORS["observed"],
                    marker="v",
                    s=35,
                    label="Límite superior",
                )
            ax.set_title(f"{SPECIES_LABELS[species]} · {run}")
            ax.set_xlabel("Tiempo (h)")
            ax.set_ylabel("Aroma en vino (µg L⁻¹)")
            _style_axis(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=3,
        frameon=False,
    )
    fig.suptitle("Validación aroma holdout: ajuste sin 26158 ni 26211", y=0.995, fontsize=14)
    return _save_figure(fig, "holdout_wine_predictions.png", top=0.90)


def plot_holdout_metrics(metrics: pd.DataFrame) -> Path:
    view = metrics[
        metrics["fit_scope"].eq("calibration_only")
        & metrics["role"].eq("holdout")
        & metrics["domain"].eq("wine")
    ]
    pivot = view.pivot(index="species", columns="forcing_source", values="rmse").reindex(
        SPECIES_LABELS
    )
    x = np.arange(len(pivot))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 4.8))
    for offset, source in zip((-width / 2, width / 2), SOURCE_LABELS):
        ax.bar(
            x + offset,
            pivot[source],
            width,
            color=COLORS[source],
            label=SOURCE_LABELS[source],
        )
    ax.set_xticks(x, [SPECIES_LABELS[item] for item in pivot.index])
    ax.set_ylabel("RMSE holdout vino (µg L⁻¹)")
    ax.set_title("Desempeño cruzado del forcing de volatilización")
    ax.legend(frameon=False)
    _style_axis(ax)
    return _save_figure(fig, "holdout_wine_rmse_comparison.png")


def plot_condensate_predictions(condensate_predictions: pd.DataFrame) -> Path:
    view = condensate_predictions[
        condensate_predictions["forcing_source"].eq("dynamic_transfer_model")
        & condensate_predictions["fit_scope"].eq("calibration_only")
        & condensate_predictions["role"].eq("holdout")
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    for ax, species in zip(axes, SPECIES_LABELS):
        group = view[view["species"].eq(species)]
        direct = group[group["status"].eq("observed")]
        censored = group[~group["status"].eq("observed")]
        ax.scatter(
            direct["observed_or_upper_bound"],
            direct["predicted_captured_ug"],
            color=COLORS["dynamic_transfer_model"],
            label="Observado",
        )
        ax.scatter(
            censored["observed_or_upper_bound"],
            censored["predicted_captured_ug"],
            facecolors="none",
            edgecolors=COLORS["holdout"],
            marker="v",
            label="Límite superior",
        )
        maximum = max(
            float(group["observed_or_upper_bound"].max()) if len(group) else 1.0,
            float(group["predicted_captured_ug"].max()) if len(group) else 1.0,
            1.0,
        )
        ax.plot([0, maximum], [0, maximum], color="#7B8794", ls="--", lw=1)
        ax.set_xlim(0, maximum * 1.05)
        ax.set_ylim(0, maximum * 1.05)
        ax.set_title(SPECIES_LABELS[species])
        ax.set_xlabel("Observado o límite (µg)")
        ax.set_ylabel("Predicho capturado (µg)")
        _style_axis(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Condensado holdout con transferencia gas-líquido", y=0.995, fontsize=14)
    return _save_figure(fig, "holdout_condensate_predictions.png", top=0.87)


def plot_parameters(parameters: pd.DataFrame) -> Path:
    view = parameters[parameters["fit_scope"].eq("all_data")]
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.6))
    for ax, parameter in zip(axes[:2], aroma.FORMATION_PARAMETER_NAMES):
        group = view[view["parameter"].eq(parameter)]
        pivot = group.pivot(index="species", columns="forcing_source", values="estimate").reindex(
            SPECIES_LABELS
        )
        x = np.arange(len(pivot))
        width = 0.36
        for offset, source in zip((-width / 2, width / 2), SOURCE_LABELS):
            ax.bar(
                x + offset,
                pivot[source],
                width,
                color=COLORS[source],
                label=SOURCE_LABELS[source],
            )
        ax.set_xticks(x, [SPECIES_LABELS[item].replace(" ", "\n", 1) for item in pivot.index])
        ax.set_yscale("log")
        ax.set_title(parameter.replace("_", " "))
        _style_axis(ax)
    ax = axes[2]
    transfer = view[
        view["forcing_source"].eq("dynamic_transfer_model")
        & view["parameter"].eq("mass_transfer_kla_ref_h_inv")
    ].set_index("species").reindex(SPECIES_LABELS)
    x = np.arange(len(transfer))
    ax.bar(x, transfer["estimate"], color=COLORS["dynamic_transfer_model"])
    ax.set_xticks(x, [SPECIES_LABELS[item].replace(" ", "\n", 1) for item in transfer.index])
    ax.set_yscale("log")
    ax.set_ylabel("h⁻¹")
    ax.set_title("kLa de referencia")
    _style_axis(ax)
    ax = axes[3]
    multiplier = view[
        view["forcing_source"].eq("dynamic_transfer_model")
        & view["parameter"].eq("ethanol_kla_multiplier_per_10_g_l")
    ].set_index("species").reindex(SPECIES_LABELS)
    x = np.arange(len(multiplier))
    ax.bar(x, multiplier["estimate"], color=COLORS["dynamic_transfer_model"])
    ax.axhline(1.0, color=COLORS["grid"], lw=1.0)
    ax.set_xticks(x, [SPECIES_LABELS[item].replace(" ", "\n", 1) for item in multiplier.index])
    ax.set_yscale("log")
    ax.set_ylabel("factor por 10 g/L")
    ax.set_title("Efecto etanol sobre kLa")
    _style_axis(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.94),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("Parámetros recalibrados con todos los ensayos", y=0.995, fontsize=14)
    return _save_figure(fig, "aroma_parameter_comparison.png", top=0.86)


def plot_profiles(profiles: pd.DataFrame) -> Path:
    view = profiles[
        profiles["forcing_source"].eq("dynamic_transfer_model")
        & profiles["fit_scope"].eq("all_data")
    ]
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharex=False)
    for row_index, species in enumerate(SPECIES_LABELS):
        for col_index, parameter in enumerate(aroma.DYNAMIC_TRANSFER_PARAMETER_NAMES):
            ax = axes[row_index, col_index]
            group = view[
                view["species"].eq(species) & view["parameter"].eq(parameter)
            ].sort_values("fixed_value")
            ax.plot(group["fixed_value"], group["delta_objective"], marker="o")
            ax.axhline(3.84, color=COLORS["holdout"], ls="--", lw=1)
            ax.set_xscale("log")
            if row_index == 0:
                ax.set_title(parameter.replace("_", " "))
            if col_index == 0:
                ax.set_ylabel(f"{SPECIES_LABELS[species]}\nΔ objetivo")
            ax.set_xlabel("Valor fijo")
            _style_axis(ax)
    fig.suptitle("Perfiles de identificabilidad · transferencia gas-líquido", y=0.995, fontsize=14)
    return _save_figure(fig, "release_forcing_profile_likelihood.png", top=0.94)


def plot_loss_dynamics(transfer: pd.DataFrame, pulses: pd.DataFrame) -> Path:
    view = transfer[transfer["species"].eq("ethyl_octanoate")]
    pulse_lookup = pulses.set_index("batch")["pulse_time_h"].to_dict()
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for ax, run in zip(axes.flat, RUNS):
        for source in SOURCE_LABELS:
            group = view[
                view["forcing_source"].eq(source)
                & view["experiment_id"].eq(run)
            ].sort_values("time_h")
            ax.plot(
                group["time_h"],
                group["loss_coefficient_h_inv"],
                color=COLORS[source],
                lw=2.0,
                label=SOURCE_LABELS[source],
            )
        pulse = pulse_lookup.get(run, np.nan)
        if np.isfinite(pulse):
            ax.axvline(pulse, color=COLORS["holdout"], ls="--", lw=1.1)
        ax.set_title(f"{run} · {pilot_co2.LOT_BY_RUN[run]}")
        ax.set_xlabel("Tiempo (h)")
        _style_axis(ax)
    axes[0, 0].set_ylabel("λ de pérdida (h⁻¹)")
    axes[1, 0].set_ylabel("λ de pérdida (h⁻¹)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=2, frameon=False)
    fig.suptitle("Octanoato de etilo · dinámica efectiva de volatilización", y=0.995, fontsize=14)
    return _save_figure(fig, "ethyl_octanoate_loss_dynamics.png", top=0.89)


def plot_partition_decomposition(forcing_grid: pd.DataFrame, pulses: pd.DataFrame) -> Path:
    view = forcing_grid[
        forcing_grid["source"].eq("dynamic_transfer_model")
        & forcing_grid["species"].eq("ethyl_octanoate")
    ]
    pulse_lookup = pulses.set_index("batch")["pulse_time_h"].to_dict()
    styles = {
        "K total": ("partition_k_gas_over_liquid", "#27313B"),
        "Etanol": ("log_k_ethanol_term", "#326FA8"),
        "Temperatura": ("log_k_temperature_term", "#D97732"),
        "Azúcar": ("log_k_sugar_term", "#5A9367"),
    }
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for ax, run in zip(axes.flat, RUNS):
        group = view[view["experiment_id"].eq(run)].sort_values("time_h")
        for label, (column, color) in styles.items():
            values = group[column].to_numpy(dtype=float)
            relative = values / values[0] if label == "K total" else np.exp(values - values[0])
            ax.plot(group["time_h"], relative, color=color, lw=1.8, label=label)
        pulse = pulse_lookup.get(run, np.nan)
        if np.isfinite(pulse):
            ax.axvline(pulse, color=COLORS["holdout"], ls="--", lw=1.1)
        ax.axhline(1.0, color=COLORS["grid"], lw=0.8)
        ax.set_title(f"{run} · {pilot_co2.LOT_BY_RUN[run]}")
        ax.set_xlabel("Tiempo (h)")
        _style_axis(ax)
    axes[0, 0].set_ylabel("Factor relativo al inicio")
    axes[1, 0].set_ylabel("Factor relativo al inicio")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4, frameon=False)
    fig.suptitle("Octanoato de etilo · descomposición de K(T, etanol, azúcar)", y=0.995, fontsize=14)
    return _save_figure(fig, "ethyl_octanoate_partition_decomposition.png", top=0.89)


def _summary_payload(
    metrics: pd.DataFrame,
    exposure: pd.DataFrame,
    validation: pd.DataFrame,
    condensate_predictions: pd.DataFrame,
) -> dict[str, Any]:
    holdout = metrics[
        metrics["fit_scope"].eq("calibration_only")
        & metrics["role"].eq("holdout")
    ]
    pivot = holdout.pivot(
        index=["domain", "species"], columns="forcing_source", values="rmse"
    )
    comparisons = {}
    for (domain, species), row in pivot.iterrows():
        old = float(row["released_co2_model"])
        new = float(row["dynamic_transfer_model"])
        comparisons[f"{domain}__{species}"] = {
            "dynamic_transfer_rmse": new,
            "empirical_mass_basis_rmse": old,
            "relative_change": (new - old) / old if old > 0.0 else np.nan,
        }
    identified = validation[
        validation["forcing_source"].eq("dynamic_transfer_model")
        & validation["fit_scope"].eq("all_data")
    ].set_index("species")["loss_separately_identified"].to_dict()
    mix03 = condensate_predictions[
        condensate_predictions["fit_scope"].eq("all_data")
        & condensate_predictions["species"].eq("ethyl_octanoate")
        & condensate_predictions["experiment_id"].eq("26157")
        & condensate_predictions["mix_id"].eq("26157-MIX-03")
    ]
    mix03_payload = {
        str(row.forcing_source): float(row.predicted_captured_ug)
        for row in mix03.itertuples(index=False)
    }
    if not mix03.empty:
        mix03_payload["observed_captured_ug"] = float(
            mix03.iloc[0]["observed_or_upper_bound"]
        )
    return {
        "holdout_rmse_comparison": comparisons,
        "median_integrated_gas_turnover": float(exposure["integrated_gas_turnover"].median()),
        "loss_identified_with_dynamic_transfer": identified,
        "ethyl_octanoate_26157_mix03": mix03_payload,
        "holdout_initialization_policy": (
            "The first wine observation initializes each trajectory and is excluded from RMSE."
        ),
        "interpretation": (
            "Both models use the same fixed released-rCO2 curve. The corrected model converts "
            "mass rate to gas turnover and estimates kLa(E); its NTU efficiency remains bounded "
            "between zero and one while K(T, ethanol, sugar) supplies thermodynamics."
        ),
    }


def run_analysis() -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = aroma.load_json(AROMA_CONFIG_PATH)
    primary_config = aroma.load_json(PRIMARY_CONFIG_PATH)
    partitions, partition_provenance = aroma.load_partition_surrogates(config, ROOT_DIR)
    tables = pilot_calibration.load_tables(MODEL_DATASET_DIR)
    released_overrides, released_grid, pulses = _build_released_co2_overrides()
    forcing_sets, wine_tables, condensate_tables, forcing_grid = _build_aroma_forcings(
        tables,
        config,
        primary_config,
        partitions,
        released_overrides,
    )
    fits = _fit_all(config, forcing_sets, wine_tables, condensate_tables)
    parameters = _collect_parameters(fits)
    validation, multistarts, profiles, covariance = _collect_fit_diagnostics(fits)
    wine_predictions, condensate_predictions = _prediction_inventory(
        fits,
        forcing_sets,
        wine_tables,
        condensate_tables,
        config,
    )
    transfer_diagnostics = _transfer_diagnostic_inventory(fits, forcing_sets)
    metrics = _score_predictions(wine_predictions, condensate_predictions)
    exposure = _forcing_exposure(forcing_grid)
    release_gate = aroma.evaluate_gate(
        [fits[("dynamic_transfer_model", "all_data", species)] for species in config["priority_analytes"]]
    )
    summary = _summary_payload(
        metrics, exposure, validation, condensate_predictions
    )

    frames = {
        "rco2_release_driver_grid.csv": released_grid,
        "aroma_forcing_grid.csv": forcing_grid,
        "rco2_exposure_summary.csv": exposure,
        "aroma_parameter_estimates.csv": parameters,
        "aroma_fit_validation.csv": validation,
        "aroma_multistart_summary.csv": multistarts,
        "aroma_profile_likelihood.csv": profiles,
        "aroma_covariance_log_space.csv": covariance,
        "aroma_wine_predictions.csv": wine_predictions,
        "aroma_condensate_predictions.csv": condensate_predictions,
        "aroma_prediction_metrics.csv": metrics,
        "aroma_transfer_diagnostics.csv": transfer_diagnostics,
        "nutrient_pulse_schedule.csv": pulses,
    }
    for name, frame in frames.items():
        frame.to_csv(RESULTS_DIR / name, index=False)

    figure_paths = [
        plot_rco2_forcing(forcing_grid, pulses),
        plot_exposure(exposure),
        plot_holdout_wine(wine_predictions),
        plot_holdout_metrics(metrics),
        plot_condensate_predictions(condensate_predictions),
        plot_parameters(parameters),
        plot_profiles(profiles),
        plot_loss_dynamics(transfer_diagnostics, pulses),
        plot_partition_decomposition(forcing_grid, pulses),
    ]

    _write_json(RESULTS_DIR / "calibration_summary.json", summary)
    _write_json(RESULTS_DIR / "aroma_calibration_gate.json", release_gate)
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis": "pilot_2026_aroma_dynamic_transfer_recalibration",
            "model": lab_co2.MODEL_NAME,
            "rco2_definition": "predicted emitted CO2 mass rate after matrix_gain",
            "gas_flow_conversion": "Qgas/VL = rCO2 / ideal-gas CO2 density at reactor temperature and 1 atm",
            "transfer_model": "kLa(E)=kLa_ref*mE^((E-50)/10); lambda=K(T,E,S)*(Qgas/VL)*(1-exp(-kLa(E)/(Qgas/VL)))",
            "runs": list(RUNS),
            "calibration_runs": sorted(CALIBRATION_RUNS),
            "holdout_runs": sorted(HOLDOUTS),
            "priority_analytes": config["priority_analytes"],
            "sources": {
                "aroma_config": {
                    "path": str(AROMA_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(AROMA_CONFIG_PATH),
                },
                "co2_fit_parameters": {
                    "path": str((CO2_RESULTS_DIR / "fit_parameters.csv").relative_to(ROOT_DIR)),
                    "sha256": _sha256(CO2_RESULTS_DIR / "fit_parameters.csv"),
                },
                "model_dataset": str(MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
                "primary_calibration": str(PRIMARY_CALIBRATION_DIR.relative_to(ROOT_DIR)),
                "historical_aroma_calibration": str(HISTORICAL_AROMA_DIR.relative_to(ROOT_DIR)),
                "partition_surrogate": partition_provenance,
            },
            "outputs": [str(path.relative_to(RESULTS_DIR)) for path in figure_paths]
            + sorted(frames),
            "gate": release_gate,
            "summary": summary,
        },
    )
    return {
        "parameters": parameters,
        "validation": validation,
        "metrics": metrics,
        "exposure": exposure,
        "wine_predictions": wine_predictions,
        "condensate_predictions": condensate_predictions,
        "profiles": profiles,
        "transfer_diagnostics": transfer_diagnostics,
        "summary": summary,
        "gate": release_gate,
        "figures": figure_paths,
    }


def load_results() -> dict[str, Any]:
    """Load a completed analysis for fast notebook rendering/QA."""

    required = {
        "parameters": "aroma_parameter_estimates.csv",
        "validation": "aroma_fit_validation.csv",
        "metrics": "aroma_prediction_metrics.csv",
        "exposure": "rco2_exposure_summary.csv",
        "wine_predictions": "aroma_wine_predictions.csv",
        "condensate_predictions": "aroma_condensate_predictions.csv",
        "profiles": "aroma_profile_likelihood.csv",
        "transfer_diagnostics": "aroma_transfer_diagnostics.csv",
    }
    missing = [name for name in required.values() if not (RESULTS_DIR / name).exists()]
    if missing:
        raise RuntimeError(f"Cannot reuse incomplete aroma results: {missing}")
    result = {key: pd.read_csv(RESULTS_DIR / name) for key, name in required.items()}
    result["summary"] = json.loads(
        (RESULTS_DIR / "calibration_summary.json").read_text(encoding="utf-8")
    )
    result["gate"] = json.loads(
        (RESULTS_DIR / "aroma_calibration_gate.json").read_text(encoding="utf-8")
    )
    result["figures"] = sorted(FIGURE_DIR.glob("*.png"))
    return result


def create_notebook() -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    notebook = nbformat.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3"}
    notebook["cells"] = [
        nbformat.v4.new_markdown_cell(
            "# Piloto 2026 · aromas con transferencia gas–líquido"
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
import os
sys.path.insert(0, str(ROOT / "fermentation_model"))
from pilot_2026 import run_aroma_co2_release_recalibration_2026 as analysis
result = analysis.load_results() if os.environ.get("PILOT_AROMA_REUSE_RESULTS") == "1" else analysis.run_analysis()
print("Resultados:", analysis.RESULTS_DIR.relative_to(ROOT))"""
        ),
        nbformat.v4.new_markdown_cell("## tl;dr"),
        nbformat.v4.new_code_cell(
            """summary = result["summary"]
display(summary["holdout_rmse_comparison"])
print("Mediana de recambios gaseosos acumulados:", round(summary["median_integrated_gas_turnover"], 4))
print("Pérdida identificada:", summary["loss_identified_with_dynamic_transfer"])
print("26157 MIX-03:", summary["ethyl_octanoate_26157_mix03"])
print("Gate:", result["gate"]["verdict"])"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contexto y métodos

- Se conservan formación aromática, partición UNIFAC, captura, errores y censura.
- Ambos comparadores usan el mismo rCO₂ emitido predicho por `solubility_o2_nitrogen_boost_continuous_release`.
- La corrección convierte masa de CO₂ a `Qgas/VL` con gas ideal y usa `1-exp(-kLa(E)/(Qgas/VL))`.
- `kLa(E)=kLa_ref*mE**((E-50)/10)` permite contrastar directamente el efecto adicional del etanol.
- `K(T, etanol, azúcar)` sigue aportando la dependencia termodinámica con la composición.
- `26158` y `26211` son holdouts tanto para CO₂ como para aromas.
- La primera muestra de vino inicializa la trayectoria y no entra al RMSE de validación; las muestras siguientes sí.
- Se estiman `mass_transfer_kla_ref_h_inv` y el multiplicador de etanol, con eficiencia NTU acotada entre cero y uno."""
        ),
        nbformat.v4.new_markdown_cell("## Datos y forcing rCO₂"),
        nbformat.v4.new_code_cell(
            """display(result["exposure"].round(4))
display(Image(filename=analysis.FIGURE_DIR / "gas_flow_conversion.png"))
display(Image(filename=analysis.FIGURE_DIR / "integrated_rco2_exposure.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Resultados de calibración y validación"),
        nbformat.v4.new_code_cell(
            """display(result["parameters"].query("fit_scope == 'all_data'").round(5))
display(result["metrics"].query("fit_scope == 'calibration_only' and role == 'holdout'").round(4))
display(Image(filename=analysis.FIGURE_DIR / "holdout_wine_predictions.png"))
display(Image(filename=analysis.FIGURE_DIR / "holdout_wine_rmse_comparison.png"))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / "holdout_condensate_predictions.png"))
display(Image(filename=analysis.FIGURE_DIR / "aroma_parameter_comparison.png"))
display(Image(filename=analysis.FIGURE_DIR / "ethyl_octanoate_loss_dynamics.png"))
display(Image(filename=analysis.FIGURE_DIR / "ethyl_octanoate_partition_decomposition.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Identificabilidad"),
        nbformat.v4.new_code_cell(
            """display(result["validation"].query("forcing_source == 'dynamic_transfer_model' and fit_scope == 'all_data'").round(4))
display(Image(filename=analysis.FIGURE_DIR / "release_forcing_profile_likelihood.png"))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

La comparación relevante es el desempeño holdout con el mismo rCO₂, cambiando solamente el operador de volatilización. Una mejora no implica por sí sola identificabilidad completa: `kLa` puede saturarse cuando la corriente gaseosa ya sale cerca del equilibrio. Los límites de detección se incorporan como censura, no como ceros."""
        ),
        nbformat.v4.new_code_cell(
            """print("Notebook ejecutado sin errores.")
print("Figuras embebidas:", len(result["figures"]))
print("Artefactos:", analysis.RESULTS_DIR.relative_to(ROOT))"""
        ),
    ]
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    import tempfile

    runtime_dir = Path(tempfile.gettempdir()) / "pilot26_aroma_jupyter"
    runtime_dir.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime_dir)
    os.environ["IPYTHONDIR"] = str(runtime_dir / "ipython")
    os.environ["JUPYTER_ALLOW_INSECURE_WRITES"] = "1"
    os.environ["PILOT_AROMA_REUSE_RESULTS"] = "1"
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
