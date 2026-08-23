from __future__ import annotations

"""Test equilibrium aroma stripping with UNIFAC and Morakul partition models.

The experiment keeps the validated emitted-rCO2 forcing, aroma observations,
capture efficiencies, censoring and calibration/holdout split fixed. It compares:

1. equilibrium with the existing UNIFAC surrogate and no transfer parameter;
2. equilibrium with published Morakul/Mouret k_i(T, ethanol) coefficients;
3. the UNIFAC model with a liquid-side kLa and the corrected NTU basis.

The previous gas-side-NTU result is loaded only as a historical benchmark.
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

from pilot_2026 import run_aroma_co2_release_recalibration_2026 as previous  # noqa: E402
from pilot_2026 import run_co2_solubility_cross_lot_validation_2026 as pilot_co2  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_calibration  # noqa: E402


MODEL_DATASET_DIR = pilot_co2.MODEL_DATASET_DIR
PRIMARY_CALIBRATION_DIR = pilot_co2.PILOT_CALIBRATION_DIR
AROMA_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "aroma_calibration_config.json"
TEST_CONFIG_PATH = (
    SCRIPT_DIR / "adaptive_design" / "aroma_equilibrium_partition_test_config.json"
)
PRIMARY_CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "calibration_config.json"
PREVIOUS_RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_dynamic_transfer_recalibration_2026"
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_equilibrium_partition_test_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_equilibrium_partition_test.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_equilibrium_partition_test.executed.ipynb"
)

RUNS = tuple(pilot_co2.VALID_RUNS)
HOLDOUTS = frozenset(pilot_co2.HOLDOUTS.values())
CALIBRATION_RUNS = frozenset(run for run in RUNS if run not in HOLDOUTS)
MODEL_LABELS = {
    "equilibrium_unifac": "Equilibrio - UNIFAC",
    "equilibrium_morakul": "Equilibrio - Morakul",
    "liquid_kla_unifac": "kLa líquido - UNIFAC",
    "historical_gas_ntu": "NTU previo (referencia)",
}
SPECIES_LABELS = {
    "ethyl_acetate": "Acetato de etilo",
    "ethyl_octanoate": "Octanoato de etilo",
    "isoamyl_acetate": "Acetato de isoamilo",
}
COLORS = {
    "equilibrium_unifac": "#326FA8",
    "equilibrium_morakul": "#D97732",
    "liquid_kla_unifac": "#5E8C61",
    "historical_gas_ntu": "#8A7AA8",
    "observed": "#222A33",
    "censored": "#B45A7A",
    "pulse": "#CC3D3D",
    "grid": "#D8DEE5",
}
RANDOM_SEED = 20260821


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


def _literature_partition_models(
    test_config: dict[str, Any],
    unifac: dict[str, dict[str, float]],
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    models: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    for species, payload in test_config["published_partition_models"].items():
        model = dict(payload)
        model["trap_efficiency"] = float(unifac[species]["trap_efficiency"])
        models[species] = model
        rows.append(
            {
                "species": species,
                "F1": float(payload["F1"]),
                "F2_l_g": float(payload["F2_l_g"]),
                "F3_kj_mol": float(payload["F3_kj_mol"]),
                "F4_kj_l_mol_g": float(payload["F4_kj_l_mol_g"]),
                "reference_temperature_k": float(
                    payload["reference_temperature_k"]
                ),
                "source_matrix": payload["source_matrix"],
                "parameter_source": payload["parameter_source"],
                "doi": payload["doi"],
                "trap_efficiency": model["trap_efficiency"],
            }
        )
    return models, pd.DataFrame(rows)


def _build_forcing_sets(
    tables: pilot_calibration.CalibrationTables,
    aroma_config: dict[str, Any],
    primary_config: dict[str, Any],
    unifac: dict[str, dict[str, float]],
    morakul: dict[str, dict[str, Any]],
    released_overrides: dict[str, tuple[np.ndarray, np.ndarray]],
) -> tuple[
    dict[str, dict[str, dict[str, aroma.AromaForcing]]],
    dict[str, pd.DataFrame],
    dict[str, pd.DataFrame],
    pd.DataFrame,
]:
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]] = {
        model: {} for model in MODEL_LABELS if model != "historical_gas_ntu"
    }
    wine_tables: dict[str, pd.DataFrame] = {}
    condensate_tables: dict[str, pd.DataFrame] = {}
    rows: list[dict[str, Any]] = []

    for species, analyte_label in aroma_config["priority_analytes"].items():
        variants = (
            (
                "equilibrium_unifac",
                unifac[species],
                aroma.LOSS_MODEL_EQUILIBRIUM,
                "unifac_surrogate",
            ),
            (
                "equilibrium_morakul",
                morakul[species],
                aroma.LOSS_MODEL_EQUILIBRIUM,
                "morakul_mouret_empirical",
            ),
            (
                "liquid_kla_unifac",
                unifac[species],
                aroma.LOSS_MODEL_DYNAMIC_TRANSFER,
                "unifac_surrogate",
            ),
        )
        for variant, partition, loss_model, partition_family in variants:
            forcings, wine, condensate, _ = aroma.build_forcings(
                tables,
                primary_config,
                PRIMARY_CALIBRATION_DIR,
                species,
                analyte_label,
                partition,
                co2_rate_override_by_run=released_overrides,
                loss_model=loss_model,
            )
            forcing_sets[variant][species] = forcings
            wine_tables[species] = wine
            condensate_tables[species] = condensate
            for run, forcing in forcings.items():
                for values in zip(
                    forcing.time_h,
                    forcing.co2_rate_g_l_h,
                    forcing.gas_turnover_h_inv,
                    forcing.partition_basis_l_g,
                    forcing.temperature_c,
                    forcing.ethanol_g_l,
                    forcing.total_sugar_g_l,
                ):
                    time_h, rco2, qgas, partition_k, temp_c, ethanol, sugar = values
                    rows.append(
                        {
                            "model_variant": variant,
                            "partition_family": partition_family,
                            "species": species,
                            "experiment_id": run,
                            "matrix": pilot_co2.LOT_BY_RUN[run],
                            "role": "holdout" if run in HOLDOUTS else "calibration",
                            "time_h": float(time_h),
                            "rco2_g_l_h": float(rco2),
                            "gas_turnover_h_inv": float(qgas),
                            "partition_k_gas_over_liquid": float(partition_k),
                            "temperature_c": float(temp_c),
                            "ethanol_g_l": float(ethanol),
                            "total_sugar_g_l": float(sugar),
                        }
                    )
    return forcing_sets, wine_tables, condensate_tables, pd.DataFrame(rows)


def _fit_all(
    aroma_config: dict[str, Any],
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]],
    wine_tables: dict[str, pd.DataFrame],
    condensate_tables: dict[str, pd.DataFrame],
) -> dict[tuple[str, str, str], aroma.AromaFit]:
    fits: dict[tuple[str, str, str], aroma.AromaFit] = {}
    for model_index, variant in enumerate(forcing_sets):
        for fit_scope, run_set in (
            ("calibration_only", CALIBRATION_RUNS),
            ("all_data", frozenset(RUNS)),
        ):
            fit_config = copy.deepcopy(aroma_config)
            fit_config["optimization"]["sobol_multistarts"] = (
                8 if variant == "liquid_kla_unifac" else 12
            )
            if fit_scope == "calibration_only":
                fit_config["optimization"]["profile_grid_log_offsets"] = [0.0]
            fit_config["optimization"]["seed"] = (
                RANDOM_SEED + model_index * 1000 + (100 if fit_scope == "all_data" else 0)
            )
            for species in aroma_config["priority_analytes"]:
                forcings, wine, condensate = previous._subset_inputs(
                    forcing_sets[variant][species],
                    wine_tables[species],
                    condensate_tables[species],
                    run_set,
                )
                fits[(variant, fit_scope, species)] = aroma.fit_species(
                    species, forcings, wine, condensate, fit_config
                )
    return fits


def _collect_parameters(
    fits: dict[tuple[str, str, str], aroma.AromaFit]
) -> pd.DataFrame:
    frames = []
    for (variant, fit_scope, species), fit in fits.items():
        frame = fit.parameter_table.copy()
        frame.insert(0, "fit_scope", fit_scope)
        frame.insert(0, "model_variant", variant)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _collect_fit_validation(
    fits: dict[tuple[str, str, str], aroma.AromaFit]
) -> pd.DataFrame:
    rows = []
    for (variant, fit_scope, species), fit in fits.items():
        rows.append(
            {
                "model_variant": variant,
                "fit_scope": fit_scope,
                "species": species,
                "objective": fit.validation["objective"],
                "objective_per_observation": fit.validation[
                    "objective_per_observation"
                ],
                "converged_multistarts": fit.validation["converged_multistarts"],
                "active_bound_fraction": fit.validation["active_bound_fraction"],
                "weak_parameters": ";".join(fit.validation["weak_parameters"]),
                "loss_separately_identified": fit.validation[
                    "loss_separately_identified"
                ],
                "loss_parameter_fixed_by_physics": fit.validation[
                    "loss_parameter_fixed_by_physics"
                ],
                "loss_model": fit.validation["loss_model"],
            }
        )
    return pd.DataFrame(rows)


def _transfer_inventory(
    fits: dict[tuple[str, str, str], aroma.AromaFit],
    forcing_sets: dict[str, dict[str, dict[str, aroma.AromaForcing]]],
) -> pd.DataFrame:
    frames = []
    for variant in forcing_sets:
        for species in SPECIES_LABELS:
            fit = fits[(variant, "all_data", species)]
            for run, forcing in forcing_sets[variant][species].items():
                frame = aroma.transfer_diagnostics(forcing, fit.log_values)
                frame.insert(0, "role", "holdout" if run in HOLDOUTS else "calibration")
                frame.insert(0, "matrix", pilot_co2.LOT_BY_RUN[run])
                frame.insert(0, "experiment_id", run)
                frame.insert(0, "species", species)
                frame.insert(0, "model_variant", variant)
                frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _historical_metrics() -> pd.DataFrame:
    path = PREVIOUS_RESULTS_DIR / "aroma_prediction_metrics.csv"
    if not path.exists():
        return pd.DataFrame()
    metrics = pd.read_csv(path)
    metrics = metrics[metrics["forcing_source"].eq("dynamic_transfer_model")].copy()
    metrics = metrics.rename(columns={"forcing_source": "model_variant"})
    metrics["model_variant"] = "historical_gas_ntu"
    return metrics


def _with_historical_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    historical = _historical_metrics()
    if historical.empty:
        return metrics.copy()
    common = [column for column in metrics.columns if column in historical.columns]
    return pd.concat([metrics, historical[common]], ignore_index=True)


def _metric_lookup(
    metrics: pd.DataFrame,
    variant: str,
    species: str,
    domain: str,
    metric: str,
) -> float:
    row = metrics[
        metrics["model_variant"].eq(variant)
        & metrics["fit_scope"].eq("calibration_only")
        & metrics["role"].eq("holdout")
        & metrics["species"].eq(species)
        & metrics["domain"].eq(domain)
    ]
    if row.empty:
        return float("nan")
    return float(row.iloc[0][metric])


def _validation_summary(
    metrics_with_history: pd.DataFrame,
    validation: pd.DataFrame,
    test_config: dict[str, Any],
) -> tuple[dict[str, Any], pd.DataFrame]:
    policy = test_config["validation_policy"]
    benchmark = "historical_gas_ntu"
    rows: list[dict[str, Any]] = []
    for variant in ("equilibrium_unifac", "equilibrium_morakul", "liquid_kla_unifac"):
        for species in SPECIES_LABELS:
            for domain in ("wine", "condensate"):
                current_rmse = _metric_lookup(
                    metrics_with_history, variant, species, domain, "rmse"
                )
                benchmark_rmse = _metric_lookup(
                    metrics_with_history, benchmark, species, domain, "rmse"
                )
                rows.append(
                    {
                        "model_variant": variant,
                        "species": species,
                        "domain": domain,
                        "rmse": current_rmse,
                        "benchmark_rmse": benchmark_rmse,
                        "relative_change_vs_historical": (
                            (current_rmse - benchmark_rmse) / benchmark_rmse
                            if np.isfinite(current_rmse)
                            and np.isfinite(benchmark_rmse)
                            and benchmark_rmse > 0.0
                            else np.nan
                        ),
                    }
                )
    comparison = pd.DataFrame(rows)
    candidate_rows = []
    for variant in ("equilibrium_unifac", "equilibrium_morakul", "liquid_kla_unifac"):
        cond_changes = comparison[
            comparison["model_variant"].eq(variant)
            & comparison["domain"].eq("condensate")
            & comparison["species"].isin(policy["primary_condensate_species"])
        ]["relative_change_vs_historical"].dropna()
        wine_changes = comparison[
            comparison["model_variant"].eq(variant)
            & comparison["domain"].eq("wine")
        ]["relative_change_vs_historical"].dropna()
        cond_pass = bool(
            len(cond_changes) == len(policy["primary_condensate_species"])
            and np.all(
                cond_changes
                <= -float(policy["minimum_condensate_nrmse_reduction_fraction"])
            )
        )
        wine_pass = bool(
            len(wine_changes) == len(SPECIES_LABELS)
            and np.all(
                wine_changes
                <= float(policy["maximum_wine_rmse_degradation_fraction"])
            )
        )
        fit_rows = validation[
            validation["model_variant"].eq(variant)
            & validation["fit_scope"].eq("calibration_only")
        ]
        computational_pass = bool(
            len(fit_rows) == len(SPECIES_LABELS)
            and (fit_rows["converged_multistarts"] >= 6).all()
        )
        candidate_rows.append(
            {
                "model_variant": variant,
                "condensate_target_pass": cond_pass,
                "wine_guardrail_pass": wine_pass,
                "computational_pass": computational_pass,
                "mean_condensate_relative_change": float(cond_changes.mean())
                if len(cond_changes)
                else np.nan,
                "maximum_wine_relative_change": float(wine_changes.max())
                if len(wine_changes)
                else np.nan,
            }
        )
    candidates = pd.DataFrame(candidate_rows)
    passing = candidates[
        candidates["condensate_target_pass"]
        & candidates["wine_guardrail_pass"]
        & candidates["computational_pass"]
    ]
    if len(passing):
        verdict = "PASS"
        selected = str(
            passing.sort_values("mean_condensate_relative_change").iloc[0][
                "model_variant"
            ]
        )
    elif candidates["computational_pass"].all():
        verdict = "INFORMATIVE_NEGATIVE"
        selected = None
    else:
        verdict = "COMPUTATIONAL_FAILURE"
        selected = None
    return (
        {
            "verdict": verdict,
            "selected_model": selected,
            "criteria": policy,
            "candidate_checks": candidates.to_dict(orient="records"),
            "interpretation": (
                "PASS requires at least 30% lower holdout condensate error for both "
                "observed ester species, no more than 10% wine degradation, and stable "
                "optimization. INFORMATIVE_NEGATIVE means the physical hypothesis was "
                "tested successfully but did not meet predictive gates."
            ),
        },
        comparison,
    )


def _style_axis(ax: plt.Axes) -> None:
    ax.grid(color=COLORS["grid"], lw=0.6, alpha=0.7)
    ax.spines[["top", "right"]].set_visible(False)


def _save_figure(fig: plt.Figure, name: str, *, top: float = 0.94) -> Path:
    path = FIGURE_DIR / name
    fig.tight_layout(rect=(0.0, 0.0, 1.0, top))
    fig.savefig(path, dpi=175, bbox_inches="tight")
    plt.close(fig)
    return path


def _pulse_by_run(pulses: pd.DataFrame) -> dict[str, float]:
    return {
        str(row.batch): float(row.pulse_time_h)
        for row in pulses.itertuples(index=False)
    }


def plot_partition_trajectories(
    forcing_grid: pd.DataFrame, pulses: pd.DataFrame
) -> Path:
    pulse_map = _pulse_by_run(pulses)
    runs = ("26157", "26211")
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=False)
    for row_index, species in enumerate(SPECIES_LABELS):
        for col_index, run in enumerate(runs):
            ax = axes[row_index, col_index]
            for variant in ("equilibrium_unifac", "equilibrium_morakul"):
                group = forcing_grid[
                    forcing_grid["model_variant"].eq(variant)
                    & forcing_grid["species"].eq(species)
                    & forcing_grid["experiment_id"].eq(run)
                ].sort_values("time_h")
                ax.plot(
                    group["time_h"],
                    group["partition_k_gas_over_liquid"],
                    color=COLORS[variant],
                    lw=2,
                    label=MODEL_LABELS[variant],
                )
            ax.axvline(
                pulse_map[run], color=COLORS["pulse"], ls=":", lw=1.3, label="Pulso N"
            )
            ax.set_yscale("log")
            ax.set_title(f"{SPECIES_LABELS[species]} - {run}")
            ax.set_xlabel("Tiempo de proceso (h)")
            ax.set_ylabel("K = Cgas / Clíquido")
            _style_axis(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Comparación de partición: UNIFAC vs Morakul/Mouret", y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.972),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, "partition_trajectories.png", top=0.91)


def plot_holdout_wine(
    predictions: pd.DataFrame, pulses: pd.DataFrame
) -> Path:
    pulse_map = _pulse_by_run(pulses)
    runs = tuple(sorted(HOLDOUTS))
    view = predictions[
        predictions["fit_scope"].eq("calibration_only")
        & predictions["role"].eq("holdout")
    ]
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=False)
    for row_index, species in enumerate(SPECIES_LABELS):
        for col_index, run in enumerate(runs):
            ax = axes[row_index, col_index]
            group = view[
                view["species"].eq(species)
                & view["experiment_id"].eq(run)
            ]
            observed = group[
                group["model_variant"].eq("equilibrium_unifac")
                & group["status"].eq("observed")
            ]
            ax.scatter(
                observed["time_h"],
                observed["observed_or_upper_bound"],
                color=COLORS["observed"],
                s=32,
                zorder=4,
                label="Química",
            )
            for variant in ("equilibrium_unifac", "equilibrium_morakul", "liquid_kla_unifac"):
                model = group[group["model_variant"].eq(variant)].sort_values("time_h")
                ax.plot(
                    model["time_h"],
                    model["predicted_ug_l"],
                    color=COLORS[variant],
                    lw=1.8,
                    marker="o",
                    ms=3,
                    label=MODEL_LABELS[variant],
                )
            ax.axvline(
                pulse_map[run], color=COLORS["pulse"], ls=":", lw=1.3, label="Pulso N"
            )
            ax.set_title(f"{SPECIES_LABELS[species]} - {run}")
            ax.set_xlabel("Tiempo de proceso (h)")
            ax.set_ylabel("Concentración en vino (ug/L)")
            _style_axis(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Validación holdout en vino", y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.972),
        ncol=5,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, "holdout_wine_predictions.png", top=0.90)


def plot_holdout_condensate(
    predictions: pd.DataFrame, pulses: pd.DataFrame
) -> Path:
    pulse_map = _pulse_by_run(pulses)
    runs = tuple(sorted(HOLDOUTS))
    view = predictions[
        predictions["fit_scope"].eq("calibration_only")
        & predictions["role"].eq("holdout")
    ].copy()
    view["midpoint_h"] = 0.5 * (
        view["interval_start_h"] + view["interval_end_h"]
    )
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=False)
    for row_index, species in enumerate(SPECIES_LABELS):
        for col_index, run in enumerate(runs):
            ax = axes[row_index, col_index]
            group = view[
                view["species"].eq(species)
                & view["experiment_id"].eq(run)
            ]
            observed = group[group["model_variant"].eq("equilibrium_unifac")]
            direct = observed[observed["status"].eq("observed")]
            censored = observed[~observed["status"].eq("observed")]
            ax.scatter(
                direct["midpoint_h"],
                direct["observed_or_upper_bound"],
                color=COLORS["observed"],
                s=34,
                zorder=4,
                label="Observado",
            )
            ax.scatter(
                censored["midpoint_h"],
                censored["observed_or_upper_bound"],
                facecolors="none",
                edgecolors=COLORS["censored"],
                marker="v",
                s=38,
                label="Límite superior",
            )
            for variant in ("equilibrium_unifac", "equilibrium_morakul", "liquid_kla_unifac"):
                model = group[group["model_variant"].eq(variant)].sort_values("midpoint_h")
                ax.plot(
                    model["midpoint_h"],
                    model["predicted_captured_ug"],
                    color=COLORS[variant],
                    lw=1.8,
                    marker="o",
                    ms=3,
                    label=MODEL_LABELS[variant],
                )
            ax.axvline(
                pulse_map[run], color=COLORS["pulse"], ls=":", lw=1.3, label="Pulso N"
            )
            ax.set_yscale("symlog", linthresh=1.0)
            ax.set_title(f"{SPECIES_LABELS[species]} - {run}")
            ax.set_xlabel("Tiempo medio del MIX (h)")
            ax.set_ylabel("Masa capturada (ug)")
            _style_axis(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Validación holdout en condensado y pulso nutricional", y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.972),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, "holdout_condensate_predictions.png", top=0.86)


def plot_holdout_metrics(metrics: pd.DataFrame) -> Path:
    view = metrics[
        metrics["fit_scope"].eq("calibration_only")
        & metrics["role"].eq("holdout")
    ]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    for ax, domain in zip(axes, ("wine", "condensate")):
        domain_view = view[view["domain"].eq(domain)]
        pivot = domain_view.pivot(
            index="species", columns="model_variant", values="nrmse_over_mean"
        ).reindex(SPECIES_LABELS)
        variants = [variant for variant in MODEL_LABELS if variant in pivot.columns]
        x = np.arange(len(pivot))
        width = 0.8 / max(len(variants), 1)
        for index, variant in enumerate(variants):
            ax.bar(
                x - 0.4 + width / 2 + index * width,
                pivot[variant],
                width,
                color=COLORS[variant],
                label=MODEL_LABELS[variant],
            )
        ax.set_xticks(x, [SPECIES_LABELS[item].replace(" ", "\n", 1) for item in pivot.index])
        ax.set_ylabel("NRMSE / media observada")
        ax.set_title("Vino" if domain == "wine" else "Condensado")
        _style_axis(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle("Desempeño holdout frente al NTU histórico", y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=4,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, "holdout_metric_comparison.png", top=0.82)


def plot_effective_loss(
    transfer: pd.DataFrame, pulses: pd.DataFrame
) -> Path:
    pulse_map = _pulse_by_run(pulses)
    runs = ("26157", "26211")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, run in zip(axes, runs):
        view = transfer[
            transfer["species"].eq("ethyl_octanoate")
            & transfer["experiment_id"].eq(run)
        ]
        for variant in ("equilibrium_unifac", "equilibrium_morakul", "liquid_kla_unifac"):
            group = view[view["model_variant"].eq(variant)].sort_values("time_h")
            ax.plot(
                group["time_h"],
                group["loss_coefficient_h_inv"],
                color=COLORS[variant],
                lw=2,
                label=MODEL_LABELS[variant],
            )
        ax.axvline(
            pulse_map[run], color=COLORS["pulse"], ls=":", lw=1.3, label="Pulso N"
        )
        ax.set_title(f"Octanoato de etilo - {run}")
        ax.set_xlabel("Tiempo de proceso (h)")
        ax.set_ylabel("Coeficiente de pérdida (h-1)")
        _style_axis(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle("Pérdida efectiva después de corregir la base NTU", y=0.995)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=4,
        frameon=False,
        fontsize=9,
    )
    return _save_figure(fig, "ethyl_octanoate_effective_loss.png", top=0.80)


def _summary_payload(
    metrics: pd.DataFrame,
    validation_gate: dict[str, Any],
    condensate_predictions: pd.DataFrame,
) -> dict[str, Any]:
    holdout = metrics[
        metrics["fit_scope"].eq("calibration_only")
        & metrics["role"].eq("holdout")
    ]
    metric_records = holdout[
        ["model_variant", "species", "domain", "rmse", "nrmse_over_mean", "bias"]
    ].to_dict(orient="records")
    mix03 = condensate_predictions[
        condensate_predictions["fit_scope"].eq("all_data")
        & condensate_predictions["species"].eq("ethyl_octanoate")
        & condensate_predictions["mix_id"].eq("26157-MIX-03")
    ]
    mix03_payload = {
        str(row.model_variant): float(row.predicted_captured_ug)
        for row in mix03.itertuples(index=False)
    }
    if len(mix03):
        mix03_payload["observed_captured_ug"] = float(
            mix03.iloc[0]["observed_or_upper_bound"]
        )
    return {
        "verdict": validation_gate["verdict"],
        "selected_model": validation_gate["selected_model"],
        "holdout_metrics": metric_records,
        "ethyl_octanoate_26157_mix03": mix03_payload,
        "holdout_initialization_policy": (
            "The first wine observation initializes each trajectory and is excluded from RMSE."
        ),
        "caveat": (
            "Morakul parameters were identified in natural must for ethyl and isoamyl "
            "acetate, and in synthetic pilot fermentations for ethyl octanoate. Matrix "
            "transfer is tested here, not assumed."
        ),
    }


def run_analysis() -> dict[str, Any]:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    aroma_config = aroma.load_json(AROMA_CONFIG_PATH)
    test_config = aroma.load_json(TEST_CONFIG_PATH)
    primary_config = aroma.load_json(PRIMARY_CONFIG_PATH)
    unifac, unifac_provenance = aroma.load_partition_surrogates(
        aroma_config, ROOT_DIR
    )
    morakul, coefficient_table = _literature_partition_models(test_config, unifac)
    tables = pilot_calibration.load_tables(MODEL_DATASET_DIR)
    released_overrides, released_grid, pulses = previous._build_released_co2_overrides()
    forcing_sets, wine_tables, condensate_tables, forcing_grid = _build_forcing_sets(
        tables,
        aroma_config,
        primary_config,
        unifac,
        morakul,
        released_overrides,
    )
    fits = _fit_all(aroma_config, forcing_sets, wine_tables, condensate_tables)
    parameters = _collect_parameters(fits)
    validation = _collect_fit_validation(fits)
    wine_predictions, condensate_predictions = previous._prediction_inventory(
        fits, forcing_sets, wine_tables, condensate_tables, aroma_config
    )
    transfer = _transfer_inventory(fits, forcing_sets)
    metrics = previous._score_predictions(wine_predictions, condensate_predictions)
    wine_predictions = wine_predictions.rename(columns={"forcing_source": "model_variant"})
    condensate_predictions = condensate_predictions.rename(columns={"forcing_source": "model_variant"})
    metrics = metrics.rename(columns={"forcing_source": "model_variant"})
    metrics_with_history = _with_historical_metrics(metrics)
    gate, comparison = _validation_summary(
        metrics_with_history, validation, test_config
    )
    summary = _summary_payload(metrics_with_history, gate, condensate_predictions)

    frames = {
        "partition_model_coefficients.csv": coefficient_table,
        "rco2_release_driver_grid.csv": released_grid,
        "aroma_forcing_grid.csv": forcing_grid,
        "aroma_parameter_estimates.csv": parameters,
        "aroma_fit_validation.csv": validation,
        "aroma_wine_predictions.csv": wine_predictions,
        "aroma_condensate_predictions.csv": condensate_predictions,
        "aroma_prediction_metrics.csv": metrics_with_history,
        "aroma_transfer_diagnostics.csv": transfer,
        "model_vs_historical_comparison.csv": comparison,
        "nutrient_pulse_schedule.csv": pulses,
    }
    for name, frame in frames.items():
        frame.to_csv(RESULTS_DIR / name, index=False)

    figure_paths = [
        plot_partition_trajectories(forcing_grid, pulses),
        plot_holdout_wine(wine_predictions, pulses),
        plot_holdout_condensate(condensate_predictions, pulses),
        plot_holdout_metrics(metrics_with_history),
        plot_effective_loss(transfer, pulses),
    ]

    _write_json(RESULTS_DIR / "validation_gate.json", gate)
    _write_json(RESULTS_DIR / "calibration_summary.json", summary)
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis": "pilot_2026_aroma_equilibrium_partition_test",
            "runs": list(RUNS),
            "calibration_runs": sorted(CALIBRATION_RUNS),
            "holdout_runs": sorted(HOLDOUTS),
            "models": test_config["model_variants"],
            "rco2_definition": "validated predicted emitted CO2 mass rate",
            "sources": {
                "aroma_config": {
                    "path": str(AROMA_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(AROMA_CONFIG_PATH),
                },
                "test_config": {
                    "path": str(TEST_CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(TEST_CONFIG_PATH),
                },
                "unifac_surrogate": unifac_provenance,
                "morakul_2011_doi": "10.1016/j.procbio.2011.01.034",
                "mouret_2014_doi": "10.1016/j.foodres.2014.02.044",
                "model_dataset": str(MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
                "previous_benchmark": str(PREVIOUS_RESULTS_DIR.relative_to(ROOT_DIR)),
            },
            "outputs": sorted(frames)
            + [str(path.relative_to(RESULTS_DIR)) for path in figure_paths],
            "gate": gate,
            "summary": summary,
        },
    )
    return {
        "parameters": parameters,
        "validation": validation,
        "wine_predictions": wine_predictions,
        "condensate_predictions": condensate_predictions,
        "metrics": metrics_with_history,
        "comparison": comparison,
        "coefficients": coefficient_table,
        "forcing_grid": forcing_grid,
        "transfer": transfer,
        "pulses": pulses,
        "gate": gate,
        "summary": summary,
        "figures": figure_paths,
    }


def load_results() -> dict[str, Any]:
    required = {
        "parameters": "aroma_parameter_estimates.csv",
        "validation": "aroma_fit_validation.csv",
        "wine_predictions": "aroma_wine_predictions.csv",
        "condensate_predictions": "aroma_condensate_predictions.csv",
        "metrics": "aroma_prediction_metrics.csv",
        "comparison": "model_vs_historical_comparison.csv",
        "coefficients": "partition_model_coefficients.csv",
        "forcing_grid": "aroma_forcing_grid.csv",
        "transfer": "aroma_transfer_diagnostics.csv",
        "pulses": "nutrient_pulse_schedule.csv",
    }
    result = {key: pd.read_csv(RESULTS_DIR / name) for key, name in required.items()}
    result["gate"] = aroma.load_json(RESULTS_DIR / "validation_gate.json")
    result["summary"] = aroma.load_json(RESULTS_DIR / "calibration_summary.json")
    result["figures"] = sorted(FIGURE_DIR.glob("*.png"))
    return result


def create_notebook(summary: dict[str, Any] | None = None) -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    if summary is None and (RESULTS_DIR / "calibration_summary.json").exists():
        summary = aroma.load_json(RESULTS_DIR / "calibration_summary.json")
    verdict = summary["verdict"] if summary else "PENDIENTE"
    mix03 = summary.get("ethyl_octanoate_26157_mix03", {}) if summary else {}
    notebook = nbformat.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["metadata"]["language_info"] = {"name": "python", "version": "3"}
    notebook["cells"] = [
        nbformat.v4.new_markdown_cell(
            "# Piloto 2026 - prueba M0 de pérdida aromática en equilibrio"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto ejecutado:** `{verdict}`.

Esta prueba elimina el parámetro libre de transferencia para contrastar dos modelos de partición y, en paralelo, corrige la base NTU del modelo con `kLa` líquido. El caso diagnóstico 26157-MIX-03 queda registrado en la salida ejecutada: `{json.dumps(mix03, ensure_ascii=False)}`.

El criterio de éxito exige mejorar al menos 30 % el error holdout de condensado para octanoato e isoamilo sin degradar más de 10 % el ajuste de vino."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Context & Methods

### Key Assumptions

- La química líquida inicializa cada trayectoria; la primera muestra no entra al RMSE.
- Los holdouts completos son 26158 y 26211.
- El `rCO2` emitido proviene del modelo de CO2 validado previamente.
- `K=Cgas/Clíquido`; el modelo de equilibrio usa `lambda=K*(Qgas/VL)`.
- Morakul/Mouret representa composición mediante etanol y temperatura. Los coeficientes de acetatos provienen de mosto natural y los de octanoato de mosto sintético piloto.
- El tren de captura conserva las eficiencias históricas; por ello esta prueba aísla partición/transferencia, no valida todavía el condensador.
- Los puntos censurados se tratan como límites, no como ceros.

Fuentes primarias: Morakul et al. 2011, DOI 10.1016/j.procbio.2011.01.034; Mouret et al. 2014, DOI 10.1016/j.foodres.2014.02.044."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import os
import sys
from IPython.display import display, Image

ROOT = Path.cwd()
while ROOT != ROOT.parent and not (ROOT / "fermentation_model").exists():
    ROOT = ROOT.parent
if not (ROOT / "fermentation_model").exists():
    raise RuntimeError("Execute from the repository or a descendant directory")
sys.path.insert(0, str(ROOT / "fermentation_model"))
from pilot_2026 import run_aroma_equilibrium_partition_test_2026 as analysis
result = analysis.load_results() if os.environ.get("PILOT_AROMA_REUSE_RESULTS") == "1" else analysis.run_analysis()
print("Resultados:", analysis.RESULTS_DIR.relative_to(ROOT))
print("Veredicto:", result["gate"]["verdict"])"""
        ),
        nbformat.v4.new_markdown_cell("## Data"),
        nbformat.v4.new_code_cell(
            """display(result["coefficients"])
display(result["pulses"][["batch", "pulse_time_h", "timing_source"]].round(3))
display(Image(filename=analysis.FIGURE_DIR / "partition_trajectories.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Results"),
        nbformat.v4.new_code_cell(
            """holdout = result["metrics"].query("fit_scope == 'calibration_only' and role == 'holdout'")
display(holdout[["model_variant", "species", "domain", "n_observed_scored", "rmse", "nrmse_over_mean", "bias"]].round(4))
display(result["comparison"].round(4))
display(Image(filename=analysis.FIGURE_DIR / "holdout_metric_comparison.png"))"""
        ),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / "holdout_wine_predictions.png"))
display(Image(filename=analysis.FIGURE_DIR / "holdout_condensate_predictions.png"))"""
        ),
        nbformat.v4.new_markdown_cell("## Transfer and identifiability checks"),
        nbformat.v4.new_code_cell(
            """display(result["validation"].round(5))
display(result["parameters"].query("fit_scope == 'all_data'").round(6))
display(Image(filename=analysis.FIGURE_DIR / "ethyl_octanoate_effective_loss.png"))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- Si equilibrio-Morakul mejora condensado y respeta el guardrail de vino, el siguiente paso es el estado de headspace/trampa.
- Si el modelo sigue fallando en MIX-04 o sobrepasa el techo calculado con química interpolada, el bloqueo está en producción temporal, captura o resolución de muestreo, no en un `kLa` adicional.
- Un `kLa` corregido que caiga al límite inferior indica que el optimizador intenta reconstruir artificialmente la penalización histórica.
- Ningún resultado de acetato de etilo permite declarar identificabilidad de pérdida porque todos sus condensados están censurados."""
        ),
        nbformat.v4.new_code_cell(
            """assert result["gate"]["verdict"] in {"PASS", "INFORMATIVE_NEGATIVE", "COMPUTATIONAL_FAILURE"}
assert len(result["figures"]) == 5
print("Notebook ejecutado sin errores.")
print("Figuras embebidas:", len(result["figures"]))
print("Gate:", result["gate"]["verdict"])"""
        ),
    ]
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    import tempfile

    runtime_dir = Path(tempfile.gettempdir()) / "pilot26_aroma_equilibrium_jupyter"
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
    if args.create_only:
        create_notebook()
    else:
        result = run_analysis()
        create_notebook(result["summary"])
        execute_notebook(timeout=args.timeout)
    print(NOTEBOOK_PATH)
    if not args.create_only:
        print(EXECUTED_NOTEBOOK_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
