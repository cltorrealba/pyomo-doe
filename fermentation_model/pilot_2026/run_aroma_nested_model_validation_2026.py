from __future__ import annotations

"""Validate nested aroma production/capture models on pilot fermentations 2026.

The analysis performs complete leave-one-reactor-run-out validation for two
quantified target aromas.  Model selection is decided by pre-registered gates,
with a separate sensitivity analysis for the flagged 26211-P-12 wine sample.
"""

import argparse
import copy
import hashlib
import json
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

from pilot_2026 import run_aroma_co2_release_recalibration_2026 as released  # noqa: E402
from pilot_2026 import run_aroma_equilibrium_partition_test_2026 as equilibrium  # noqa: E402
from pilot_2026 import run_co2_solubility_cross_lot_validation_2026 as pilot_co2  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_nested_models as nested  # noqa: E402
from pilot_2026.adaptive_design import pilot_calibration  # noqa: E402


MODEL_DATASET_DIR = pilot_co2.MODEL_DATASET_DIR
PRIMARY_CALIBRATION_DIR = pilot_co2.PILOT_CALIBRATION_DIR
AROMA_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "aroma_calibration_config.json"
EQUILIBRIUM_CONFIG_PATH = (
    SCRIPT_DIR / "adaptive_design" / "aroma_equilibrium_partition_test_config.json"
)
NESTED_CONFIG_PATH = (
    SCRIPT_DIR / "adaptive_design" / "aroma_nested_validation_config.json"
)
PRIMARY_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "calibration_config.json"
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_nested_model_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_nested_model_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_nested_model_validation.executed.ipynb"
)

SPECIES_LABELS = {
    "ethyl_octanoate": "Octanoato de etilo",
    "isoamyl_acetate": "Acetato de isoamilo",
}
MODEL_LABELS = {
    nested.BASELINE: "Equilibrio basal",
    nested.DELAYED: "Respuesta pospulso",
    nested.RESERVOIR: "Reservorio de línea",
    nested.COMBINED: "Pospulso + reservorio",
    nested.CAPTURE: "Eficiencia de captura ajustada",
    nested.DELAYED_CAPTURE: "Pospulso + captura ajustada",
    nested.RESERVOIR_CAPTURE: "Reservorio + captura ajustada",
    nested.FULL_CAPTURE: "Pospulso + reservorio + captura",
    nested.ETHANOL_CAPTURE: "Captura dependiente de etanol",
    nested.ETHANOL_CAPTURE_RESERVOIR: "Captura por etanol + reservorio",
}
COLORS = {
    nested.BASELINE: "#4C78A8",
    nested.DELAYED: "#E6863B",
    nested.RESERVOIR: "#54A24B",
    nested.COMBINED: "#8E6C8A",
    nested.CAPTURE: "#2A9D8F",
    nested.DELAYED_CAPTURE: "#F4A261",
    nested.RESERVOIR_CAPTURE: "#6A994E",
    nested.FULL_CAPTURE: "#7B2CBF",
    nested.ETHANOL_CAPTURE: "#E76F51",
    nested.ETHANOL_CAPTURE_RESERVOIR: "#5A189A",
    "observed": "#20262E",
    "temperature": "#C44E52",
    "pulse": "#A51C30",
    "grid": "#D9DEE5",
}
RUNS = tuple(str(run) for run in pilot_co2.VALID_RUNS)


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _subset(
    forcings: dict[str, aroma.AromaForcing],
    wine: pd.DataFrame,
    condensate: pd.DataFrame,
    runs: set[str],
) -> tuple[dict[str, aroma.AromaForcing], pd.DataFrame, pd.DataFrame]:
    return (
        {run: forcing for run, forcing in forcings.items() if run in runs},
        wine[wine["experiment_id"].astype(str).isin(runs)].copy(),
        condensate[condensate["experiment_id"].astype(str).isin(runs)].copy(),
    )


def _build_inputs() -> tuple[
    dict[str, dict[str, aroma.AromaForcing]],
    dict[str, pd.DataFrame],
    dict[str, pd.DataFrame],
    dict[str, float],
    pd.DataFrame,
    pd.DataFrame,
    pilot_calibration.CalibrationTables,
]:
    aroma_config = aroma.load_json(AROMA_CONFIG_PATH)
    equilibrium_config = aroma.load_json(EQUILIBRIUM_CONFIG_PATH)
    primary_config = aroma.load_json(PRIMARY_CONFIG_PATH)
    unifac, _ = aroma.load_partition_surrogates(aroma_config, ROOT_DIR)
    morakul, _ = equilibrium._literature_partition_models(
        equilibrium_config, unifac
    )
    tables = pilot_calibration.load_tables(MODEL_DATASET_DIR)
    overrides, rco2_grid, pulses = released._build_released_co2_overrides()
    pulse_lookup = {
        str(row.batch): float(row.pulse_time_h)
        for row in pulses.itertuples(index=False)
        if str(row.batch) in RUNS
    }
    missing_pulses = set(RUNS) - set(pulse_lookup)
    if missing_pulses:
        raise RuntimeError(f"Missing nutrient pulse times: {sorted(missing_pulses)}")

    forcing_sets: dict[str, dict[str, aroma.AromaForcing]] = {}
    wine_tables: dict[str, pd.DataFrame] = {}
    condensate_tables: dict[str, pd.DataFrame] = {}
    for species, analyte_label in aroma_config["priority_analytes"].items():
        if species not in SPECIES_LABELS:
            continue
        forcings, wine, condensate, _ = aroma.build_forcings(
            tables,
            primary_config,
            PRIMARY_CALIBRATION_DIR,
            species,
            analyte_label,
            morakul[species],
            co2_rate_override_by_run=overrides,
            loss_model=aroma.LOSS_MODEL_EQUILIBRIUM,
        )
        forcing_sets[species] = {
            run: forcing for run, forcing in forcings.items() if run in RUNS
        }
        wine_tables[species] = wine[wine["experiment_id"].astype(str).isin(RUNS)].copy()
        condensate_tables[species] = condensate[
            condensate["experiment_id"].astype(str).isin(RUNS)
        ].copy()
    return (
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        pulses,
        rco2_grid,
        tables,
    )


def _data_quality_audit(tables: pilot_calibration.CalibrationTables) -> tuple[pd.DataFrame, pd.DataFrame]:
    wine = tables.wine_aroma.copy()
    wine = wine[
        wine["experiment_id"].astype(str).isin(RUNS)
        & wine["model_observation_type"].eq("observed")
    ].copy()
    rows: list[dict[str, Any]] = []
    for (run, analyte), group in wine.groupby(["experiment_id", "analyte"], sort=True):
        group = group.sort_values("time_h").reset_index(drop=True)
        for index in range(1, len(group) - 1):
            previous = group.iloc[index - 1]
            current = group.iloc[index]
            following = group.iloc[index + 1]
            previous_value = float(previous["observed_value"])
            current_value = float(current["observed_value"])
            following_value = float(following["observed_value"])
            if min(previous_value, current_value, following_value) <= 0.0:
                continue
            fraction = (
                float(current["time_h"]) - float(previous["time_h"])
            ) / max(float(following["time_h"]) - float(previous["time_h"]), 1e-12)
            expected = float(
                np.exp(
                    np.log(previous_value)
                    + fraction * (np.log(following_value) - np.log(previous_value))
                )
            )
            rows.append(
                {
                    "experiment_id": str(run),
                    "sample_id": str(current["sample_id"]),
                    "time_h": float(current["time_h"]),
                    "analyte": str(analyte),
                    "observed_value": current_value,
                    "log_linear_neighbor_expectation": expected,
                    "fold_vs_neighbor_expectation": current_value / expected,
                }
            )
    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby(["experiment_id", "sample_id", "time_h"], as_index=False)
        .agg(
            analytes_compared=("analyte", "size"),
            median_fold=("fold_vs_neighbor_expectation", "median"),
            maximum_fold=("fold_vs_neighbor_expectation", "max"),
        )
        .sort_values("median_fold", ascending=False)
        .reset_index(drop=True)
    )
    summary["flagged_26211_P12"] = summary["sample_id"].eq("26211-P-12")
    return detail, summary


def _fit_config(
    base: dict[str, Any], model_index: int, species_index: int, fold_index: int
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    config["optimization"]["seed"] = (
        int(base["optimization"]["seed"])
        + 1000 * model_index
        + 100 * species_index
        + fold_index
    )
    return config


def _decorate_predictions(
    frame: pd.DataFrame,
    *,
    species: str,
    model_variant: str,
    left_out_run: str,
    analysis_policy: str,
) -> pd.DataFrame:
    result = frame.copy()
    result.insert(0, "analysis_policy", analysis_policy)
    result.insert(1, "left_out_run", left_out_run)
    result.insert(2, "model_variant", model_variant)
    result.insert(3, "species", species)
    return result


def _run_loro(
    config: dict[str, Any],
    models: list[str],
    forcing_sets: dict[str, dict[str, aroma.AromaForcing]],
    wine_tables: dict[str, pd.DataFrame],
    condensate_tables: dict[str, pd.DataFrame],
    pulse_lookup: dict[str, float],
    *,
    analysis_policy: str,
    exclude_flagged_wine: bool,
) -> dict[str, pd.DataFrame]:
    wine_predictions: list[pd.DataFrame] = []
    condensate_predictions: list[pd.DataFrame] = []
    parameter_rows: list[pd.DataFrame] = []
    validation_rows: list[dict[str, Any]] = []
    simulation_rows: list[dict[str, Any]] = []

    for model_index, model_variant in enumerate(models):
        for species_index, species in enumerate(SPECIES_LABELS):
            wine_all = wine_tables[species].copy()
            if exclude_flagged_wine:
                flagged = set(config["data_quality_policy"]["flagged_wine_sample_ids"])
                wine_all = wine_all[~wine_all["sample_id"].astype(str).isin(flagged)].copy()
            for fold_index, left_out in enumerate(RUNS):
                training_runs = set(RUNS) - {left_out}
                train_forcing, train_wine, train_condensate = _subset(
                    forcing_sets[species], wine_all, condensate_tables[species], training_runs
                )
                fit = nested.fit_nested_species(
                    species,
                    model_variant,
                    train_forcing,
                    pulse_lookup,
                    train_wine,
                    train_condensate,
                    _fit_config(config, model_index, species_index, fold_index),
                )
                parameters = fit.parameter_table.copy()
                parameters.insert(0, "analysis_policy", analysis_policy)
                parameters.insert(1, "left_out_run", left_out)
                parameter_rows.append(parameters)
                validation_rows.append(
                    {
                        "analysis_policy": analysis_policy,
                        "left_out_run": left_out,
                        "model_variant": model_variant,
                        "species": species,
                        **fit.validation,
                    }
                )

                test_forcing, test_wine, test_condensate = _subset(
                    forcing_sets[species], wine_all, condensate_tables[species], {left_out}
                )
                wine_prediction, condensate_prediction, simulations = (
                    nested.prediction_tables_nested(
                        test_forcing,
                        pulse_lookup,
                        test_wine,
                        test_condensate,
                        fit.log_values,
                        model_variant,
                        config,
                    )
                )
                wine_predictions.append(
                    _decorate_predictions(
                        wine_prediction,
                        species=species,
                        model_variant=model_variant,
                        left_out_run=left_out,
                        analysis_policy=analysis_policy,
                    )
                )
                condensate_predictions.append(
                    _decorate_predictions(
                        condensate_prediction,
                        species=species,
                        model_variant=model_variant,
                        left_out_run=left_out,
                        analysis_policy=analysis_policy,
                    )
                )
                simulation = simulations[left_out]
                simulation_rows.append(
                    {
                        "analysis_policy": analysis_policy,
                        "left_out_run": left_out,
                        "model_variant": model_variant,
                        "species": species,
                        "maximum_relative_mass_balance_error": float(
                            np.max(simulation.relative_mass_balance_error)
                        ),
                        "terminal_line_inventory_ug": float(
                            simulation.line_inventory_ug[-1]
                        ),
                        "terminal_cumulative_production_ug": float(
                            simulation.cumulative_production_ug[-1]
                        ),
                        "terminal_captured_ug": float(simulation.captured_ug[-1]),
                    }
                )
    return {
        "wine": pd.concat(wine_predictions, ignore_index=True),
        "condensate": pd.concat(condensate_predictions, ignore_index=True),
        "parameters": pd.concat(parameter_rows, ignore_index=True),
        "validation": pd.DataFrame(validation_rows),
        "simulation": pd.DataFrame(simulation_rows),
    }


def _observed_errors(
    frame: pd.DataFrame, prediction_column: str, *, wine: bool
) -> tuple[np.ndarray, np.ndarray]:
    observed = frame[frame["status"].eq("observed")].copy()
    if wine:
        observed = observed[~observed["initialization_point"]]
    actual = observed["observed_or_upper_bound"].to_numpy(dtype=float)
    error = observed[prediction_column].to_numpy(dtype=float) - actual
    return error, actual


def _score_loro(wine: pd.DataFrame, condensate: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    group_columns = ["analysis_policy", "model_variant", "species"]
    for keys, group in wine.groupby(group_columns, sort=True):
        error, actual = _observed_errors(group, "predicted_ug_l", wine=True)
        mean = float(np.mean(actual)) if len(actual) else np.nan
        rmse = float(np.sqrt(np.mean(error**2))) if len(error) else np.nan
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "domain": "wine",
                "n_observed_scored": int(len(error)),
                "observed_mean": mean,
                "rmse": rmse,
                "nrmse_over_mean": rmse / mean if mean > 0.0 else np.nan,
                "mae": float(np.mean(np.abs(error))) if len(error) else np.nan,
                "bias": float(np.mean(error)) if len(error) else np.nan,
            }
        )
    for keys, group in condensate.groupby(group_columns, sort=True):
        error, actual = _observed_errors(
            group, "predicted_captured_ug", wine=False
        )
        mean = float(np.mean(actual)) if len(actual) else np.nan
        rmse = float(np.sqrt(np.mean(error**2))) if len(error) else np.nan
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "domain": "condensate",
                "n_observed_scored": int(len(error)),
                "observed_mean": mean,
                "rmse": rmse,
                "nrmse_over_mean": rmse / mean if mean > 0.0 else np.nan,
                "mae": float(np.mean(np.abs(error))) if len(error) else np.nan,
                "bias": float(np.mean(error)) if len(error) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _runwise_condensate_sse(condensate: pd.DataFrame) -> pd.DataFrame:
    observed = condensate[condensate["status"].eq("observed")].copy()
    observed["squared_error"] = (
        observed["predicted_captured_ug"]
        - observed["observed_or_upper_bound"]
    ) ** 2
    return (
        observed.groupby(
            ["analysis_policy", "model_variant", "species", "experiment_id"],
            as_index=False,
        )
        .agg(n_observed=("squared_error", "size"), squared_error=("squared_error", "sum"))
    )


def _gate(
    config: dict[str, Any],
    metrics: pd.DataFrame,
    run_sse: pd.DataFrame,
    parameters: pd.DataFrame,
    validation: pd.DataFrame,
    simulation: pd.DataFrame,
) -> tuple[dict[str, Any], pd.DataFrame]:
    policy = config["validation_policy"]
    primary_metrics = metrics[metrics["analysis_policy"].eq("primary")]
    primary_sse = run_sse[run_sse["analysis_policy"].eq("primary")]
    comparison_rows: list[dict[str, Any]] = []
    candidate_payload: dict[str, Any] = {}
    comparison_baseline = str(
        policy.get("comparison_baseline_model", nested.BASELINE)
    )
    baseline_metrics = primary_metrics[
        primary_metrics["model_variant"].eq(comparison_baseline)
    ].set_index(["species", "domain"])

    for candidate in config["candidate_models"]:
        if candidate == comparison_baseline:
            continue
        checks: list[bool] = []
        species_payload: dict[str, Any] = {}
        for species in SPECIES_LABELS:
            baseline_cond = float(
                baseline_metrics.loc[(species, "condensate"), "nrmse_over_mean"]
            )
            baseline_wine = float(
                baseline_metrics.loc[(species, "wine"), "nrmse_over_mean"]
            )
            candidate_metrics = primary_metrics[
                primary_metrics["model_variant"].eq(candidate)
                & primary_metrics["species"].eq(species)
            ].set_index("domain")
            candidate_cond = float(
                candidate_metrics.loc["condensate", "nrmse_over_mean"]
            )
            candidate_wine = float(candidate_metrics.loc["wine", "nrmse_over_mean"])
            cond_reduction = 1.0 - candidate_cond / baseline_cond
            wine_degradation = candidate_wine / baseline_wine - 1.0

            base_run = primary_sse[
                primary_sse["model_variant"].eq(comparison_baseline)
                & primary_sse["species"].eq(species)
            ][["experiment_id", "squared_error"]].rename(
                columns={"squared_error": "baseline_squared_error"}
            )
            candidate_run = primary_sse[
                primary_sse["model_variant"].eq(candidate)
                & primary_sse["species"].eq(species)
            ][["experiment_id", "squared_error"]].rename(
                columns={"squared_error": "candidate_squared_error"}
            )
            paired = base_run.merge(candidate_run, on="experiment_id", how="inner")
            fraction_runs_better = float(
                (paired["candidate_squared_error"] < paired["baseline_squared_error"]).mean()
            )
            cond_pass = bool(
                cond_reduction
                >= float(policy["minimum_condensate_nrmse_reduction_fraction"])
            )
            wine_pass = bool(
                wine_degradation
                <= float(policy["maximum_wine_nrmse_degradation_fraction"])
            )
            run_pass = bool(
                fraction_runs_better
                >= float(
                    policy["minimum_fraction_runs_with_lower_condensate_squared_error"]
                )
            )
            checks.extend([cond_pass, wine_pass, run_pass])
            species_payload[species] = {
                "baseline_condensate_nrmse": baseline_cond,
                "candidate_condensate_nrmse": candidate_cond,
                "condensate_nrmse_reduction_fraction": cond_reduction,
                "condensate_reduction_pass": cond_pass,
                "baseline_wine_nrmse": baseline_wine,
                "candidate_wine_nrmse": candidate_wine,
                "wine_nrmse_degradation_fraction": wine_degradation,
                "wine_guardrail_pass": wine_pass,
                "fraction_runs_lower_condensate_squared_error": fraction_runs_better,
                "run_consistency_pass": run_pass,
            }
            comparison_rows.append(
                {
                    "model_variant": candidate,
                    "species": species,
                    **species_payload[species],
                }
            )

        parameter_subset = parameters[
            parameters["analysis_policy"].eq("primary")
            & parameters["model_variant"].eq(candidate)
        ]
        active_fraction = float(parameter_subset["active_bound"].mean())
        mass_error = float(
            simulation[
                simulation["analysis_policy"].eq("primary")
                & simulation["model_variant"].eq(candidate)
            ]["maximum_relative_mass_balance_error"].max()
        )
        fit_subset = validation[
            validation["analysis_policy"].eq("primary")
            & validation["model_variant"].eq(candidate)
        ]
        convergence_pass = bool(fit_subset["success"].all())
        active_pass = bool(
            active_fraction <= float(policy["maximum_active_bound_fraction"])
        )
        mass_pass = bool(
            mass_error <= float(policy["maximum_relative_mass_balance_error"])
        )
        checks.extend([convergence_pass, active_pass, mass_pass])
        candidate_payload[candidate] = {
            "species": species_payload,
            "all_fits_converged": convergence_pass,
            "active_bound_fraction": active_fraction,
            "active_bound_pass": active_pass,
            "maximum_relative_mass_balance_error": mass_error,
            "mass_balance_pass": mass_pass,
            "passes_all_gates": bool(all(checks)),
        }

    passing = [
        candidate
        for candidate, payload in candidate_payload.items()
        if payload["passes_all_gates"]
    ]
    passing.sort(
        key=lambda candidate: (
            len(nested.PARAMETER_NAMES_BY_MODEL[candidate]),
            config["candidate_models"].index(candidate),
        )
    )
    selected = passing[0] if passing else None
    condensate_metric = primary_metrics[
        primary_metrics["domain"].eq("condensate")
    ]
    diagnostic_best = str(
        condensate_metric.groupby("model_variant")["nrmse_over_mean"].mean().idxmin()
    )
    verdict = "PASS" if selected is not None else "NO_VALID_MODEL"
    return (
        {
            "verdict": verdict,
            "selected_model": selected,
            "diagnostic_best_model": diagnostic_best,
            "selection_rule": policy["selection_rule"],
            "comparison_baseline_model": comparison_baseline,
            "candidate_checks": candidate_payload,
            "external_validation_status": "retrospective_internal_only",
        },
        pd.DataFrame(comparison_rows),
    )


def _fit_all_data(
    config: dict[str, Any],
    forcing_sets: dict[str, dict[str, aroma.AromaForcing]],
    wine_tables: dict[str, pd.DataFrame],
    condensate_tables: dict[str, pd.DataFrame],
    pulse_lookup: dict[str, float],
) -> dict[str, pd.DataFrame]:
    parameter_rows: list[pd.DataFrame] = []
    validation_rows: list[dict[str, Any]] = []
    wine_rows: list[pd.DataFrame] = []
    condensate_rows: list[pd.DataFrame] = []
    state_rows: list[dict[str, Any]] = []
    for model_index, model_variant in enumerate(config["candidate_models"]):
        for species_index, species in enumerate(SPECIES_LABELS):
            fit = nested.fit_nested_species(
                species,
                model_variant,
                forcing_sets[species],
                pulse_lookup,
                wine_tables[species],
                condensate_tables[species],
                _fit_config(config, model_index, species_index, 99),
            )
            parameters = fit.parameter_table.copy()
            parameters.insert(0, "fit_scope", "all_data")
            parameter_rows.append(parameters)
            validation_rows.append(
                {
                    "fit_scope": "all_data",
                    "model_variant": model_variant,
                    "species": species,
                    **fit.validation,
                }
            )
            wine = fit.wine_predictions.copy()
            wine.insert(0, "model_variant", model_variant)
            wine.insert(1, "species", species)
            wine_rows.append(wine)
            condensate = fit.condensate_predictions.copy()
            condensate.insert(0, "model_variant", model_variant)
            condensate.insert(1, "species", species)
            condensate_rows.append(condensate)
            for run, forcing in forcing_sets[species].items():
                simulation = nested.simulate_nested(
                    forcing, pulse_lookup[run], fit.log_values, model_variant
                )
                for index, time_h in enumerate(simulation.time_h):
                    state_rows.append(
                        {
                            "model_variant": model_variant,
                            "species": species,
                            "experiment_id": run,
                            "time_h": float(time_h),
                            "liquid_ug_l": float(simulation.liquid_ug_l[index]),
                            "liquid_inventory_ug": float(
                                simulation.liquid_ug_l[index] * forcing.volume_l
                            ),
                            "volume_l": float(forcing.volume_l),
                            "captured_cumulative_ug": float(simulation.captured_ug[index]),
                            "line_inventory_ug": float(simulation.line_inventory_ug[index]),
                            "cumulative_production_ug": float(
                                simulation.cumulative_production_ug[index]
                            ),
                            "cumulative_volatilized_ug": float(
                                simulation.cumulative_volatilized_ug[index]
                            ),
                            "cumulative_drained_ug": float(
                                simulation.cumulative_drained_ug[index]
                            ),
                            "cumulative_unrecovered_ug": float(
                                simulation.cumulative_drained_ug[index]
                                - simulation.captured_ug[index]
                            ),
                            "production_rate_ug_l_h": float(
                                simulation.production_rate_ug_l_h[index]
                            ),
                            "pulse_activation": float(
                                simulation.pulse_activation[index]
                            ),
                            "capture_efficiency_fraction": float(
                                simulation.capture_efficiency_fraction[index]
                            ),
                            "relative_mass_balance_error": float(
                                simulation.relative_mass_balance_error[index]
                            ),
                            "temperature_c": float(forcing.temperature_c[index]),
                            "rco2_g_l_h": float(forcing.co2_rate_g_l_h[index]),
                            "ethanol_g_l": float(forcing.ethanol_g_l[index]),
                        }
                    )
    return {
        "parameters": pd.concat(parameter_rows, ignore_index=True),
        "validation": pd.DataFrame(validation_rows),
        "wine": pd.concat(wine_rows, ignore_index=True),
        "condensate": pd.concat(condensate_rows, ignore_index=True),
        "states": pd.DataFrame(state_rows),
    }


def mass_balance_summary(states: pd.DataFrame) -> pd.DataFrame:
    """Summarize biological production, transfer and observed recovery."""

    rows: list[dict[str, Any]] = []
    group_columns = ["model_variant", "species", "experiment_id"]
    for keys, group in states.groupby(group_columns, sort=True):
        group = group.sort_values("time_h")
        first = group.iloc[0]
        last = group.iloc[-1]
        initial = float(first["liquid_inventory_ug"])
        produced = float(last["cumulative_production_ug"])
        liquid = float(last["liquid_inventory_ug"])
        line = float(last["line_inventory_ug"])
        drained = float(last["cumulative_drained_ug"])
        recovered = float(last["captured_cumulative_ug"])
        available = initial + produced
        rows.append(
            {
                **dict(zip(group_columns, keys)),
                "terminal_time_h": float(last["time_h"]),
                "initial_liquid_inventory_ug": initial,
                "biological_production_ug": produced,
                "terminal_liquid_inventory_ug": liquid,
                "terminal_line_inventory_ug": line,
                "emitted_or_drained_ug": drained,
                "recovered_condensate_ug": recovered,
                "emitted_but_unrecovered_ug": drained - recovered,
                "mass_closure_error_ug": available - liquid - line - drained,
                "volatile_fraction_of_available_mass": drained / available
                if available > 0.0
                else np.nan,
                "effective_recovery_fraction_of_drained_mass": recovered / drained
                if drained > 0.0
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _parameter_stability(parameters: pd.DataFrame) -> pd.DataFrame:
    primary = parameters[parameters["analysis_policy"].eq("primary")]
    return (
        primary.groupby(["model_variant", "species", "parameter"], as_index=False)
        .agg(
            folds=("estimate", "size"),
            median_estimate=("estimate", "median"),
            minimum_estimate=("estimate", "min"),
            maximum_estimate=("estimate", "max"),
            log_standard_deviation=(
                "estimate",
                lambda values: float(np.std(np.log(np.asarray(values, dtype=float)), ddof=1)),
            ),
            active_bound_fraction=("active_bound", "mean"),
        )
    )


def _style_axis(axis: plt.Axes) -> None:
    axis.grid(color=COLORS["grid"], linewidth=0.6, alpha=0.8)
    axis.spines[["top", "right"]].set_visible(False)


def _save(fig: plt.Figure, name: str, *, top: float = 0.94) -> Path:
    path = FIGURE_DIR / name
    fig.tight_layout(rect=(0.0, 0.0, 1.0, top))
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_anomaly(summary: pd.DataFrame) -> Path:
    top = summary.head(12).sort_values("median_fold")
    labels = [f"{run}\n{sample}" for run, sample in zip(top["experiment_id"], top["sample_id"])]
    colors = [
        COLORS["pulse"] if flagged else "#778899"
        for flagged in top["flagged_26211_P12"]
    ]
    fig, axis = plt.subplots(figsize=(10, 5))
    axis.barh(labels, top["median_fold"], color=colors)
    axis.axvline(1.0, color="#333333", linewidth=1.0)
    axis.set_xlabel("Mediana observado / interpolación log-lineal de vecinos")
    axis.set_title("Auditoría multianalito de picos químicos internos")
    _style_axis(axis)
    return _save(fig, "01_sample_anomaly_audit.png", top=0.9)


def plot_cv_metrics(metrics: pd.DataFrame) -> Path:
    primary = metrics[metrics["analysis_policy"].eq("primary")]
    model_variants = list(dict.fromkeys(primary["model_variant"].astype(str)))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    for axis, domain in zip(axes, ("wine", "condensate")):
        subset = primary[primary["domain"].eq(domain)]
        x = np.arange(len(SPECIES_LABELS))
        width = 0.8 / max(len(model_variants), 1)
        center = (len(model_variants) - 1.0) / 2.0
        for index, model_variant in enumerate(model_variants):
            model = subset[subset["model_variant"].eq(model_variant)].set_index("species")
            values = [float(model.loc[species, "nrmse_over_mean"]) for species in SPECIES_LABELS]
            axis.bar(
                x + (index - center) * width,
                values,
                width,
                label=MODEL_LABELS.get(model_variant, model_variant),
                color=COLORS.get(model_variant, "#4C78A8"),
            )
        axis.set_xticks(x, [SPECIES_LABELS[species] for species in SPECIES_LABELS])
        axis.set_ylabel("NRMSE / media observada")
        axis.set_title("Vino" if domain == "wine" else "Condensado")
        _style_axis(axis)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("Validación cruzada leave-one-run-out")
    return _save(fig, "02_loro_nrmse_comparison.png")


def plot_drivers(states: pd.DataFrame, pulse_lookup: dict[str, float]) -> Path:
    representative_model = (
        nested.BASELINE
        if states["model_variant"].eq(nested.BASELINE).any()
        else str(states["model_variant"].iloc[0])
    )
    representative = states[
        states["model_variant"].eq(representative_model)
        & states["species"].eq("ethyl_octanoate")
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for axis, run in zip(axes.flat, RUNS):
        group = representative[
            representative["experiment_id"].astype(str).eq(run)
        ].sort_values("time_h")
        axis.plot(
            group["time_h"],
            group["rco2_g_l_h"],
            color=COLORS.get(representative_model, "#4C78A8"),
        )
        axis.axvline(pulse_lookup[run], color=COLORS["pulse"], linestyle="--", linewidth=1.1)
        twin = axis.twinx()
        twin.plot(group["time_h"], group["temperature_c"], color=COLORS["temperature"], alpha=0.7)
        twin.set_ylim(13.0, 23.0)
        if run not in {"26159", "26212"}:
            twin.set_yticklabels([])
        axis.set_title(run)
        axis.set_xlabel("Tiempo (h)")
        _style_axis(axis)
    axes[0, 0].set_ylabel("rCO₂ emitido (g L⁻¹ h⁻¹)")
    axes[1, 0].set_ylabel("rCO₂ emitido (g L⁻¹ h⁻¹)")
    fig.text(0.99, 0.5, "Temperatura (°C)", rotation=90, va="center", ha="right")
    fig.suptitle("Forzantes validados: CO₂, temperatura y pulso nutricional")
    return _save(fig, "03_rco2_temperature_pulse_drivers.png")


def plot_liquid(
    states: pd.DataFrame,
    wine: pd.DataFrame,
    pulse_lookup: dict[str, float],
    model_variant: str,
    species: str,
) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for axis, run in zip(axes.flat, RUNS):
        curve = states[
            states["model_variant"].eq(model_variant)
            & states["species"].eq(species)
            & states["experiment_id"].astype(str).eq(run)
        ].sort_values("time_h")
        points = wine[
            wine["model_variant"].eq(model_variant)
            & wine["species"].eq(species)
            & wine["experiment_id"].astype(str).eq(run)
        ]
        axis.plot(curve["time_h"], curve["liquid_ug_l"], color=COLORS[model_variant], linewidth=1.8)
        observed = points[points["status"].eq("observed")]
        axis.scatter(
            observed["time_h"],
            observed["observed_or_upper_bound"],
            color=COLORS["observed"],
            s=20,
            zorder=3,
        )
        axis.axvline(pulse_lookup[run], color=COLORS["pulse"], linestyle="--", linewidth=1.0)
        twin = axis.twinx()
        twin.plot(curve["time_h"], curve["temperature_c"], color=COLORS["temperature"], alpha=0.35)
        twin.set_ylim(13.0, 23.0)
        twin.set_yticks([])
        axis.set_title(run)
        axis.set_xlabel("Tiempo (h)")
        _style_axis(axis)
    axes[0, 0].set_ylabel("Concentración en vino (µg/L)")
    axes[1, 0].set_ylabel("Concentración en vino (µg/L)")
    fig.suptitle(f"{SPECIES_LABELS[species]} — {MODEL_LABELS[model_variant]}")
    return _save(fig, f"04_liquid_{species}.png")


def plot_condensate(
    condensate: pd.DataFrame,
    pulse_lookup: dict[str, float],
    model_variant: str,
    species: str,
) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey=True)
    for axis, run in zip(axes.flat, RUNS):
        points = condensate[
            condensate["model_variant"].eq(model_variant)
            & condensate["species"].eq(species)
            & condensate["experiment_id"].astype(str).eq(run)
        ].sort_values("interval_end_h")
        observed = points[points["status"].eq("observed")]
        if observed.empty:
            axis.text(
                0.5,
                0.5,
                "Sin muestra de\ncondensado",
                ha="center",
                va="center",
                transform=axis.transAxes,
                color="#5F6B76",
                fontsize=10,
            )
            axis.axvline(
                pulse_lookup[run], color=COLORS["pulse"], linestyle="--", linewidth=1.0
            )
            axis.set_title(run)
            axis.set_xlabel("Tiempo medio de intervalo (h)")
            _style_axis(axis)
            continue
        center = 0.5 * (observed["interval_start_h"] + observed["interval_end_h"])
        observed_mass = observed["observed_or_upper_bound"].to_numpy(dtype=float)
        predicted_mass = observed["predicted_captured_ug"].to_numpy(dtype=float)
        axis.scatter(
            center,
            observed_mass,
            marker="s",
            s=42,
            color="#4B5563",
            label="Observado",
            zorder=4,
        )
        axis.scatter(
            center,
            predicted_mass,
            marker="o",
            s=42,
            color=COLORS[model_variant],
            label="Predicho",
            zorder=4,
        )
        for x_value, observed_value, predicted_value in zip(
            center, observed_mass, predicted_mass
        ):
            if observed_value > 0.0 and predicted_value > 0.0:
                axis.plot(
                    [x_value, x_value],
                    [observed_value, predicted_value],
                    color="#C7CDD4",
                    linewidth=1.2,
                    zorder=2,
                )
        positive = np.concatenate(
            [observed_mass[observed_mass > 0.0], predicted_mass[predicted_mass > 0.0]]
        )
        if positive.size:
            axis.set_yscale("log")
            axis.set_ylim(max(float(positive.min()) / 1.8, 1e-3), float(positive.max()) * 1.8)
        ratio = np.divide(
            predicted_mass,
            observed_mass,
            out=np.full_like(predicted_mass, np.nan),
            where=observed_mass > 0.0,
        )
        finite_ratio = ratio[np.isfinite(ratio)]
        if finite_ratio.size:
            axis.text(
                0.97,
                0.95,
                f"mediana P/O = {np.median(finite_ratio):.1f}×",
                ha="right",
                va="top",
                transform=axis.transAxes,
                fontsize=8,
                color="#374151",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75},
            )
        axis.axvline(pulse_lookup[run], color=COLORS["pulse"], linestyle="--", linewidth=1.0)
        axis.set_title(run)
        axis.set_xlabel("Tiempo medio de intervalo (h)")
        _style_axis(axis)
    axes[0, 0].set_ylabel("Masa por intervalo (µg, escala log)")
    axes[1, 0].set_ylabel("Masa por intervalo (µg, escala log)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        axes[0, 0].legend(handles, labels, frameon=False, fontsize=8)
    fig.suptitle(f"{SPECIES_LABELS[species]} — condensado, {MODEL_LABELS[model_variant]}")
    return _save(fig, f"05_condensate_{species}.png")


def plot_parameter_stability(stability: pd.DataFrame, model_variant: str) -> Path:
    subset = stability[stability["model_variant"].eq(model_variant)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for axis, species in zip(axes, SPECIES_LABELS):
        group = subset[subset["species"].eq(species)].reset_index(drop=True)
        x = np.arange(len(group))
        median = group["median_estimate"].to_numpy(dtype=float)
        lower = median - group["minimum_estimate"].to_numpy(dtype=float)
        upper = group["maximum_estimate"].to_numpy(dtype=float) - median
        axis.errorbar(x, median, yerr=np.vstack([lower, upper]), fmt="o", color=COLORS[model_variant], capsize=4)
        axis.set_yscale("log")
        axis.set_xticks(x, group["parameter"], rotation=35, ha="right")
        axis.set_title(SPECIES_LABELS[species])
        axis.set_ylabel("Estimación por fold (escala log)")
        _style_axis(axis)
    fig.suptitle(f"Estabilidad leave-one-run-out — {MODEL_LABELS[model_variant]}")
    return _save(fig, "06_parameter_stability.png")


def run_analysis() -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = aroma.load_json(NESTED_CONFIG_PATH)
    (
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        pulses,
        rco2_grid,
        tables,
    ) = _build_inputs()
    anomaly_detail, anomaly_summary = _data_quality_audit(tables)

    primary = _run_loro(
        config,
        list(config["candidate_models"]),
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        analysis_policy="primary",
        exclude_flagged_wine=False,
    )
    primary_metrics = _score_loro(primary["wine"], primary["condensate"])
    primary_sse = _runwise_condensate_sse(primary["condensate"])
    provisional_gate, _ = _gate(
        config,
        primary_metrics,
        primary_sse,
        primary["parameters"],
        primary["validation"],
        primary["simulation"],
    )
    sensitivity_models = list(
        dict.fromkeys([nested.BASELINE, provisional_gate["diagnostic_best_model"]])
    )
    sensitivity = _run_loro(
        config,
        sensitivity_models,
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        analysis_policy="exclude_26211_P12_wine",
        exclude_flagged_wine=True,
    )

    cv_wine = pd.concat([primary["wine"], sensitivity["wine"]], ignore_index=True)
    cv_condensate = pd.concat(
        [primary["condensate"], sensitivity["condensate"]], ignore_index=True
    )
    cv_parameters = pd.concat(
        [primary["parameters"], sensitivity["parameters"]], ignore_index=True
    )
    cv_validation = pd.concat(
        [primary["validation"], sensitivity["validation"]], ignore_index=True
    )
    cv_simulation = pd.concat(
        [primary["simulation"], sensitivity["simulation"]], ignore_index=True
    )
    metrics = _score_loro(cv_wine, cv_condensate)
    run_sse = _runwise_condensate_sse(cv_condensate)
    gate, comparison = _gate(
        config,
        metrics,
        run_sse,
        cv_parameters,
        cv_validation,
        cv_simulation,
    )
    all_data = _fit_all_data(
        config,
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
    )
    stability = _parameter_stability(cv_parameters)
    display_model = gate["selected_model"] or gate["diagnostic_best_model"]

    frames = {
        "data_quality_anomaly_detail.csv": anomaly_detail,
        "data_quality_sample_summary.csv": anomaly_summary,
        "nutrient_pulse_schedule.csv": pulses,
        "validated_rco2_grid.csv": rco2_grid,
        "loro_wine_predictions.csv": cv_wine,
        "loro_condensate_predictions.csv": cv_condensate,
        "loro_parameter_estimates.csv": cv_parameters,
        "loro_fit_validation.csv": cv_validation,
        "loro_simulation_diagnostics.csv": cv_simulation,
        "loro_metrics.csv": metrics,
        "loro_runwise_condensate_sse.csv": run_sse,
        "validation_gate_comparison.csv": comparison,
        "parameter_stability.csv": stability,
        "all_data_parameter_estimates.csv": all_data["parameters"],
        "all_data_fit_validation.csv": all_data["validation"],
        "all_data_wine_predictions.csv": all_data["wine"],
        "all_data_condensate_predictions.csv": all_data["condensate"],
        "all_data_state_trajectories.csv": all_data["states"],
    }
    for filename, frame in frames.items():
        frame.to_csv(RESULTS_DIR / filename, index=False)

    figures = [
        plot_anomaly(anomaly_summary),
        plot_cv_metrics(metrics),
        plot_drivers(all_data["states"], pulse_lookup),
    ]
    for species in SPECIES_LABELS:
        figures.append(
            plot_liquid(
                all_data["states"],
                all_data["wine"],
                pulse_lookup,
                display_model,
                species,
            )
        )
        figures.append(
            plot_condensate(
                all_data["condensate"], pulse_lookup, display_model, species
            )
        )
    figures.append(plot_parameter_stability(stability, display_model))

    flagged = anomaly_summary[anomaly_summary["flagged_26211_P12"]]
    flagged_payload = (
        flagged.iloc[0].to_dict() if not flagged.empty else {"not_found": True}
    )
    summary = {
        "verdict": gate["verdict"],
        "selected_model": gate["selected_model"],
        "diagnostic_best_model": gate["diagnostic_best_model"],
        "display_model": display_model,
        "data_quality": {
            "primary_policy": config["data_quality_policy"]["primary"],
            "sensitivity_policy": config["data_quality_policy"]["sensitivity"],
            "flagged_sample": flagged_payload,
        },
        "interpretation": {
            "activation_peak": config["interpretation_limits"]["activation_peak"],
            "line_reservoir": config["interpretation_limits"]["line_reservoir"],
            "external_validation": config["interpretation_limits"]["external_validation"],
        },
        "literature_basis": {
            "transfer": "Morakul et al. 2011, DOI 10.1016/j.procbio.2011.01.034; Mouret et al. 2014, DOI 10.1016/j.foodres.2014.02.044",
            "nitrogen_response": "Seguinot et al. 2018, DOI 10.1016/j.fm.2018.04.005",
            "compound_specific_export": "Saerens et al. 2008, DOI 10.1128/AEM.01616-07",
        },
    }
    _write_json(RESULTS_DIR / "validation_gate.json", gate)
    _write_json(RESULTS_DIR / "analysis_summary.json", summary)
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis_id": config["analysis_id"],
            "runs": list(RUNS),
            "species": list(SPECIES_LABELS),
            "candidate_models": list(config["candidate_models"]),
            "primary_data_policy": config["data_quality_policy"]["primary"],
            "sensitivity_data_policy": config["data_quality_policy"]["sensitivity"],
            "sources": {
                "nested_config": {
                    "path": str(NESTED_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(NESTED_CONFIG_PATH),
                },
                "aroma_config": {
                    "path": str(AROMA_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(AROMA_CONFIG_PATH),
                },
                "equilibrium_config": {
                    "path": str(EQUILIBRIUM_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(EQUILIBRIUM_CONFIG_PATH),
                },
                "primary_config": {
                    "path": str(PRIMARY_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(PRIMARY_CONFIG_PATH),
                },
                "model_dataset": str(MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
            },
            "outputs": sorted(frames)
            + ["validation_gate.json", "analysis_summary.json"]
            + [str(path.relative_to(RESULTS_DIR)) for path in figures],
            "gate": gate,
        },
    )
    return {
        **frames,
        "gate": gate,
        "summary": summary,
        "figures": figures,
    }


def load_results() -> dict[str, Any]:
    files = {
        "anomaly_summary": "data_quality_sample_summary.csv",
        "metrics": "loro_metrics.csv",
        "comparison": "validation_gate_comparison.csv",
        "stability": "parameter_stability.csv",
        "parameters": "all_data_parameter_estimates.csv",
        "fit_validation": "all_data_fit_validation.csv",
        "wine": "all_data_wine_predictions.csv",
        "condensate": "all_data_condensate_predictions.csv",
        "states": "all_data_state_trajectories.csv",
        "pulses": "nutrient_pulse_schedule.csv",
    }
    result = {
        name: pd.read_csv(RESULTS_DIR / filename) for name, filename in files.items()
    }
    result["gate"] = aroma.load_json(RESULTS_DIR / "validation_gate.json")
    result["summary"] = aroma.load_json(RESULTS_DIR / "analysis_summary.json")
    result["figures"] = sorted(FIGURE_DIR.glob("*.png"))
    return result


def create_notebook(summary: dict[str, Any] | None = None) -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    if summary is None and (RESULTS_DIR / "analysis_summary.json").exists():
        summary = aroma.load_json(RESULTS_DIR / "analysis_summary.json")
    verdict = summary["verdict"] if summary else "PENDING"
    selected = summary["selected_model"] if summary else None
    notebook = nbformat.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3"}
    notebook["cells"] = [
        nbformat.v4.new_markdown_cell(
            "# Piloto 2026 — validación anidada del modelo de aromas"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto pre-registrado:** `{verdict}`.  
**Modelo seleccionado:** `{selected}`.

Se contrastan cuatro explicaciones anidadas con leave-one-reactor-run-out: equilibrio basal, respuesta biológica gradual después del pulso, reservorio de línea y ambos mecanismos. Un modelo sólo se acepta si mejora ≥30 % el NRMSE de condensado para ambos compuestos, no empeora >10 % el vino, mejora al menos la mitad de los reactores, conserva masa y no depende excesivamente de límites paramétricos."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contexto y supuestos

- El balance químico en vino y el condensado se modelan como observaciones distintas del mismo proceso.
- `rCO₂`, temperatura y pulsos se heredan del modelo piloto validado.
- La partición gas/líquido usa Morakul/Mouret y cambia con etanol y temperatura.
- La respuesta pospulso es gradual y específica por compuesto; no se interpreta automáticamente como cinética pura de nitrógeno porque el pulso y la transición térmica están confundidos en el protocolo A.
- El reservorio conserva masa, pero su constante sólo representa holdup físico si coincide con mediciones del tren de captura.
- Todos los seis reactores informaron el diseño: esta es validación interna retrospectiva, no confirmación prospectiva.
- 26211-P-12 se conserva en el primario y se excluye sólo en sensibilidad."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import os
import sys
from IPython.display import display, Image

ROOT = Path.cwd()
while ROOT != ROOT.parent and not (ROOT / 'fermentation_model').exists():
    ROOT = ROOT.parent
if not (ROOT / 'fermentation_model').exists():
    raise RuntimeError('Execute from the repository or a descendant directory')
sys.path.insert(0, str(ROOT / 'fermentation_model'))
from pilot_2026 import run_aroma_nested_model_validation_2026 as analysis
result = analysis.load_results() if os.environ.get('PILOT_AROMA_REUSE_RESULTS') == '1' else analysis.run_analysis()
print('Resultados:', analysis.RESULTS_DIR.relative_to(ROOT))
print('Veredicto:', result['gate']['verdict'])
print('Modelo seleccionado:', result['gate']['selected_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Calidad de datos"),
        nbformat.v4.new_code_cell(
            """display(result['anomaly_summary'].head(12).round(3))
display(Image(filename=analysis.FIGURE_DIR / '01_sample_anomaly_audit.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Forzantes del proceso"),
        nbformat.v4.new_code_cell(
            """display(result['pulses'][['batch', 'pulse_time_h', 'timing_source']].round(3))
display(Image(filename=analysis.FIGURE_DIR / '03_rco2_temperature_pulse_drivers.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Resultados de validación cruzada"),
        nbformat.v4.new_code_cell(
            """display(result['metrics'].round(4))
display(result['comparison'].round(4))
display(Image(filename=analysis.FIGURE_DIR / '02_loro_nrmse_comparison.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Curvas ajustadas y pérdida por intervalo"),
        nbformat.v4.new_code_cell(
            """for species in analysis.SPECIES_LABELS:
    display(Image(filename=analysis.FIGURE_DIR / f'04_liquid_{species}.png'))
    display(Image(filename=analysis.FIGURE_DIR / f'05_condensate_{species}.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Estabilidad e identificabilidad práctica"),
        nbformat.v4.new_code_cell(
            """display(result['stability'].round(4))
display(result['parameters'].round(5))
display(result['fit_validation'].round(5))
display(Image(filename=analysis.FIGURE_DIR / '06_parameter_stability.png'))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- `PASS` significa que una extensión cumplió todos los guardrails internos para ambos compuestos; no equivale a validación prospectiva.
- `NO_VALID_MODEL` significa que la evidencia actual no separa de forma transferible producción, transferencia y captura bajo los criterios definidos.
- Una respuesta pospulso útil sólo para acetato de isoamilo es biológicamente plausible, pero no satisface un modelo común si octanoato no cruza el mismo gate.
- Un reservorio con τ de decenas de horas debe tratarse como término fenomenológico hasta medir holdup, eficiencia del condensador y aroma en gas.
- El siguiente experimento discriminante es medir simultáneamente vino, gas de salida y condensado alrededor del pulso."""
        ),
        nbformat.v4.new_code_cell(
            """assert result['gate']['verdict'] in {'PASS', 'NO_VALID_MODEL'}
assert len(result['figures']) == 8
assert result['fit_validation']['maximum_relative_mass_balance_error'].max() <= 1e-8
print('Notebook ejecutado sin errores.')
print('Figuras embebidas:', len(result['figures']))"""
        ),
    ]
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    import tempfile

    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_nested_jupyter"
    runtime.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime)
    os.environ["IPYTHONDIR"] = str(runtime / "ipython")
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
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    if args.create_only:
        create_notebook()
    else:
        result = run_analysis()
        create_notebook(result["summary"])
        if not args.analysis_only:
            execute_notebook(timeout=args.timeout)
    print(NOTEBOOK_PATH)
    if not args.create_only and not args.analysis_only:
        print(EXECUTED_NOTEBOOK_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
