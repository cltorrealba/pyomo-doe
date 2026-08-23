from __future__ import annotations

"""Calibrate pilot aroma production only to wine and audit condensate closure.

The combined A+B capture fraction is fixed at one by engineering assumption.
Condensate observations are excluded from the optimizer and are scored only as
an independent closure diagnostic in leave-one-reactor-run-out validation.
"""

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")

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

from pilot_2026 import run_aroma_assumed_complete_capture_validation_2026 as complete  # noqa: E402
from pilot_2026 import run_aroma_nested_model_validation_2026 as base  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_nested_models as nested  # noqa: E402


CONFIG_PATH = (
    SCRIPT_DIR
    / "adaptive_design"
    / "aroma_wine_calibrated_complete_capture_config.json"
)
RESULTS_DIR = (
    SCRIPT_DIR / "results" / "aroma_wine_calibrated_complete_capture_2026"
)
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_PATH = (
    SCRIPT_DIR
    / "notebooks"
    / "pilot_2026_aroma_wine_calibrated_complete_capture.ipynb"
)
EXECUTED_NOTEBOOK_PATH = NOTEBOOK_PATH.with_name(
    "pilot_2026_aroma_wine_calibrated_complete_capture.executed.ipynb"
)

MODEL_LABELS = {
    nested.ASSUMED_COMPLETE: "Producción basal — ajuste sólo a vino",
    nested.ASSUMED_COMPLETE_DELAYED: "Pospulso — ajuste sólo a vino",
}
COLORS = {
    nested.ASSUMED_COMPLETE: "#4C78A8",
    nested.ASSUMED_COMPLETE_DELAYED: "#E6863B",
}


def _configure_base_output() -> None:
    base.RESULTS_DIR = RESULTS_DIR
    base.FIGURE_DIR = FIGURE_DIR
    base.MODEL_LABELS.update(MODEL_LABELS)
    base.COLORS.update(COLORS)


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
        json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default)
        + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _closure_summary(condensate_predictions: pd.DataFrame) -> pd.DataFrame:
    observed = condensate_predictions[
        condensate_predictions["analysis_policy"].eq("primary")
        & condensate_predictions["status"].eq("observed")
    ].copy()
    rows: list[dict[str, Any]] = []
    for keys, group in observed.groupby(["model_variant", "species"], sort=True):
        observed_total = float(group["observed_or_upper_bound"].sum())
        predicted_total = float(group["predicted_captured_ug"].sum())
        rows.append(
            {
                "model_variant": keys[0],
                "species": keys[1],
                "n_intervals": int(len(group)),
                "observed_total_ug": observed_total,
                "predicted_total_ug": predicted_total,
                "predicted_to_observed_ratio": (
                    predicted_total / observed_total
                    if observed_total > 0.0
                    else np.nan
                ),
                "closure_difference_ug": predicted_total - observed_total,
                "intervals_overpredicted_fraction": float(
                    (
                        group["predicted_captured_ug"]
                        > group["observed_or_upper_bound"]
                    ).mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def _validation_gate(
    config: dict[str, Any],
    metrics: pd.DataFrame,
    parameters: pd.DataFrame,
    validation: pd.DataFrame,
    simulation: pd.DataFrame,
) -> tuple[dict[str, Any], pd.DataFrame]:
    policy = config["validation_policy"]
    primary_metrics = metrics[metrics["analysis_policy"].eq("primary")]
    primary_parameters = parameters[parameters["analysis_policy"].eq("primary")]
    primary_validation = validation[validation["analysis_policy"].eq("primary")]
    primary_simulation = simulation[simulation["analysis_policy"].eq("primary")]
    baseline = nested.ASSUMED_COMPLETE
    baseline_wine = primary_metrics[
        primary_metrics["model_variant"].eq(baseline)
        & primary_metrics["domain"].eq("wine")
    ].set_index("species")["nrmse_over_mean"]
    candidates: dict[str, Any] = {}
    rows: list[dict[str, Any]] = []
    passing: list[str] = []
    for model_variant in config["candidate_models"]:
        model_metrics = primary_metrics[
            primary_metrics["model_variant"].eq(model_variant)
        ].set_index(["species", "domain"])
        species_payload: dict[str, Any] = {}
        checks: list[bool] = []
        for species in base.SPECIES_LABELS:
            wine_nrmse = float(
                model_metrics.loc[(species, "wine"), "nrmse_over_mean"]
            )
            condensate_nrmse = float(
                model_metrics.loc[(species, "condensate"), "nrmse_over_mean"]
            )
            wine_threshold = float(
                policy["maximum_wine_nrmse_over_mean"][species]
            )
            condensate_threshold = float(
                policy["maximum_condensate_nrmse_over_mean"][species]
            )
            wine_pass = bool(wine_nrmse <= wine_threshold)
            condensate_pass = bool(condensate_nrmse <= condensate_threshold)
            delayed_reduction = (
                1.0 - wine_nrmse / float(baseline_wine.loc[species])
                if model_variant != baseline
                else 0.0
            )
            delayed_support_pass = bool(
                model_variant == baseline
                or delayed_reduction
                >= float(
                    policy["minimum_delayed_wine_nrmse_reduction_fraction"]
                )
            )
            checks.extend([wine_pass, condensate_pass, delayed_support_pass])
            species_payload[species] = {
                "wine_nrmse": wine_nrmse,
                "wine_threshold": wine_threshold,
                "wine_pass": wine_pass,
                "independent_condensate_nrmse": condensate_nrmse,
                "condensate_threshold": condensate_threshold,
                "condensate_closure_pass": condensate_pass,
                "wine_nrmse_reduction_vs_baseline": delayed_reduction,
                "delayed_support_pass": delayed_support_pass,
            }
            rows.append(
                {
                    "model_variant": model_variant,
                    "species": species,
                    **species_payload[species],
                }
            )
        model_parameters = primary_parameters[
            primary_parameters["model_variant"].eq(model_variant)
        ]
        active_fraction = float(model_parameters["active_bound"].mean())
        active_pass = bool(
            active_fraction <= float(policy["maximum_active_bound_fraction"])
        )
        model_validation = primary_validation[
            primary_validation["model_variant"].eq(model_variant)
        ]
        convergence_pass = bool(model_validation["success"].all())
        domain_pass = bool(model_validation["calibration_domains"].eq("wine").all())
        max_mass_error = float(
            primary_simulation[
                primary_simulation["model_variant"].eq(model_variant)
            ]["maximum_relative_mass_balance_error"].max()
        )
        mass_pass = bool(
            max_mass_error <= float(policy["maximum_relative_mass_balance_error"])
        )
        checks.extend([active_pass, convergence_pass, domain_pass, mass_pass])
        passes = bool(all(checks))
        if passes:
            passing.append(model_variant)
        candidates[model_variant] = {
            "species": species_payload,
            "calibration_domains": "wine",
            "condensate_used_in_calibration": False,
            "active_bound_fraction": active_fraction,
            "active_bound_pass": active_pass,
            "all_fits_converged": convergence_pass,
            "wine_only_domain_pass": domain_pass,
            "maximum_relative_mass_balance_error": max_mass_error,
            "mass_balance_pass": mass_pass,
            "passes_all_gates": passes,
        }
    passing.sort(key=lambda model: len(nested.PARAMETER_NAMES_BY_MODEL[model]))
    selected = passing[0] if passing else None
    score = primary_metrics.copy()
    score["threshold"] = score.apply(
        lambda row: float(
            policy[
                "maximum_wine_nrmse_over_mean"
                if row["domain"] == "wine"
                else "maximum_condensate_nrmse_over_mean"
            ][row["species"]]
        ),
        axis=1,
    )
    score["normalized_to_gate"] = score["nrmse_over_mean"] / score["threshold"]
    diagnostic_best = str(
        score.groupby("model_variant")["normalized_to_gate"].mean().idxmin()
    )
    return (
        {
            "verdict": "PASS" if selected else "NO_VALID_MODEL",
            "selected_model": selected,
            "diagnostic_best_model": diagnostic_best,
            "selection_rule": policy["selection_rule"],
            "candidate_checks": candidates,
            "external_validation_status": "retrospective_internal_only",
        },
        pd.DataFrame(rows),
    )


def _frames_from_run(
    primary: dict[str, pd.DataFrame],
    sensitivity: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return tuple(
        pd.concat([primary[key], sensitivity[key]], ignore_index=True)
        for key in ("wine", "condensate", "parameters", "validation", "simulation")
    )


def run_analysis() -> dict[str, Any]:
    _configure_base_output()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    config = aroma.load_json(CONFIG_PATH)
    (
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        pulses,
        rco2_grid,
        tables,
    ) = base._build_inputs()
    anomaly_detail, anomaly_summary = base._data_quality_audit(tables)
    models = list(config["candidate_models"])
    primary = base._run_loro(
        config,
        models,
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        analysis_policy="primary",
        exclude_flagged_wine=False,
    )
    sensitivity = base._run_loro(
        config,
        models,
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        analysis_policy="exclude_26211_P12_wine",
        exclude_flagged_wine=True,
    )
    (
        cv_wine,
        cv_condensate,
        cv_parameters,
        cv_validation,
        cv_simulation,
    ) = _frames_from_run(primary, sensitivity)
    metrics = base._score_loro(cv_wine, cv_condensate)
    gate, comparison = _validation_gate(
        config, metrics, cv_parameters, cv_validation, cv_simulation
    )
    all_data = base._fit_all_data(
        config,
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
    )
    stability = base._parameter_stability(cv_parameters)
    mass_summary = base.mass_balance_summary(all_data["states"])
    capture_audit = complete._capture_contract_audit(
        all_data["states"], mass_summary
    )
    closure = _closure_summary(cv_condensate)
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
        "validation_gate_comparison.csv": comparison,
        "parameter_stability.csv": stability,
        "condensate_closure_summary.csv": closure,
        "all_data_parameter_estimates.csv": all_data["parameters"],
        "all_data_fit_validation.csv": all_data["validation"],
        "all_data_wine_predictions.csv": all_data["wine"],
        "all_data_condensate_predictions.csv": all_data["condensate"],
        "all_data_state_trajectories.csv": all_data["states"],
        "all_data_mass_balance_summary.csv": mass_summary,
        "capture_contract_audit.csv": capture_audit,
    }
    for filename, frame in frames.items():
        frame.to_csv(RESULTS_DIR / filename, index=False)

    figures = _make_figures(
        anomaly_summary,
        metrics,
        all_data,
        stability,
        pulse_lookup,
        display_model,
    )
    summary = _summary(config, gate, capture_audit, closure, display_model)
    _write_outputs(config, gate, summary, frames, figures)
    return {**frames, "gate": gate, "summary": summary, "figures": figures}


def _make_figures(
    anomaly_summary: pd.DataFrame,
    metrics: pd.DataFrame,
    all_data: dict[str, pd.DataFrame],
    stability: pd.DataFrame,
    pulse_lookup: dict[str, float],
    display_model: str,
) -> list[Path]:
    figures = [
        base.plot_anomaly(anomaly_summary),
        base.plot_cv_metrics(metrics),
        base.plot_drivers(all_data["states"], pulse_lookup),
    ]
    for species in base.SPECIES_LABELS:
        figures.append(
            base.plot_liquid(
                all_data["states"],
                all_data["wine"],
                pulse_lookup,
                display_model,
                species,
            )
        )
        figures.append(
            base.plot_condensate(
                all_data["condensate"], pulse_lookup, display_model, species
            )
        )
    figures.append(base.plot_parameter_stability(stability, display_model))
    return figures


def _summary(
    config: dict[str, Any],
    gate: dict[str, Any],
    capture_audit: pd.DataFrame,
    closure: pd.DataFrame,
    display_model: str,
) -> dict[str, Any]:
    display_closure = closure[closure["model_variant"].eq(display_model)]
    return {
        "verdict": gate["verdict"],
        "selected_model": gate["selected_model"],
        "diagnostic_best_model": gate["diagnostic_best_model"],
        "display_model": display_model,
        "calibration_contract": config["calibration_contract"],
        "capture_contract": config["capture_contract"],
        "capture_contract_audit": capture_audit.iloc[0].to_dict(),
        "diagnostic_closure_ratio": {
            str(row.species): float(row.predicted_to_observed_ratio)
            for row in display_closure.itertuples(index=False)
        },
        "interpretation_limits": config["interpretation_limits"],
        "literature_basis": config["literature_basis"],
    }


def _write_outputs(
    config: dict[str, Any],
    gate: dict[str, Any],
    summary: dict[str, Any],
    frames: dict[str, pd.DataFrame],
    figures: list[Path],
) -> None:
    _write_json(RESULTS_DIR / "validation_gate.json", gate)
    _write_json(RESULTS_DIR / "analysis_summary.json", summary)
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis_id": config["analysis_id"],
            "runs": list(base.RUNS),
            "species": list(base.SPECIES_LABELS),
            "candidate_models": list(config["candidate_models"]),
            "calibration_contract": config["calibration_contract"],
            "capture_contract": config["capture_contract"],
            "sources": {
                "config": {
                    "path": str(CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": _sha256(CONFIG_PATH),
                },
                "nested_model": {
                    "path": str(
                        Path(nested.__file__).resolve().relative_to(ROOT_DIR)
                    ),
                    "sha256": _sha256(Path(nested.__file__).resolve()),
                },
                "model_dataset": str(base.MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
            },
            "outputs": sorted(frames)
            + ["validation_gate.json", "analysis_summary.json"]
            + [str(path.relative_to(RESULTS_DIR)) for path in figures],
            "gate": gate,
        },
    )


def load_results() -> dict[str, Any]:
    files = {
        "metrics": "loro_metrics.csv",
        "comparison": "validation_gate_comparison.csv",
        "stability": "parameter_stability.csv",
        "closure": "condensate_closure_summary.csv",
        "parameters": "all_data_parameter_estimates.csv",
        "fit_validation": "all_data_fit_validation.csv",
        "wine": "all_data_wine_predictions.csv",
        "condensate": "all_data_condensate_predictions.csv",
        "states": "all_data_state_trajectories.csv",
        "mass": "all_data_mass_balance_summary.csv",
        "capture_audit": "capture_contract_audit.csv",
    }
    result = {
        name: pd.read_csv(RESULTS_DIR / filename) for name, filename in files.items()
    }
    result["gate"] = aroma.load_json(RESULTS_DIR / "validation_gate.json")
    result["summary"] = aroma.load_json(RESULTS_DIR / "analysis_summary.json")
    return result


def finalize_existing_results() -> dict[str, Any]:
    _configure_base_output()
    config = aroma.load_json(CONFIG_PATH)
    result = load_results()
    parameters = pd.read_csv(RESULTS_DIR / "loro_parameter_estimates.csv")
    validation = pd.read_csv(RESULTS_DIR / "loro_fit_validation.csv")
    simulation = pd.read_csv(RESULTS_DIR / "loro_simulation_diagnostics.csv")
    gate, comparison = _validation_gate(
        config, result["metrics"], parameters, validation, simulation
    )
    comparison.to_csv(RESULTS_DIR / "validation_gate_comparison.csv", index=False)
    display_model = gate["selected_model"] or gate["diagnostic_best_model"]
    anomaly_summary = pd.read_csv(
        RESULTS_DIR / "data_quality_sample_summary.csv"
    )
    pulses = pd.read_csv(RESULTS_DIR / "nutrient_pulse_schedule.csv")
    pulse_lookup = {
        str(row.batch): float(row.pulse_time_h)
        for row in pulses.itertuples(index=False)
        if str(row.batch) in base.RUNS
    }
    all_data = {
        "parameters": result["parameters"],
        "wine": result["wine"],
        "condensate": result["condensate"],
        "states": result["states"],
    }
    figures = _make_figures(
        anomaly_summary,
        result["metrics"],
        all_data,
        result["stability"],
        pulse_lookup,
        display_model,
    )
    summary = _summary(
        config, gate, result["capture_audit"], result["closure"], display_model
    )
    frames = {
        path.name: pd.DataFrame()
        for path in RESULTS_DIR.iterdir()
        if path.is_file() and path.suffix == ".csv"
    }
    _write_outputs(config, gate, summary, frames, figures)
    return {**result, "gate": gate, "summary": summary, "figures": figures}


def create_notebook() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    notebook = nbformat.v4.new_notebook()
    notebook["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook["cells"] = [
        nbformat.v4.new_markdown_cell(
            """# Aromas piloto 2026 — producción calibrada sólo con vino

## tl;dr

Los parámetros biológicos se estiman exclusivamente desde la química líquida. El condensado no participa en el objetivo y se utiliza como prueba independiente del contrato `eta_AB=1` más partición de Morakul/Mouret."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contexto y métodos

### Supuestos clave

- Captura conjunta A+B fijada a 1 por diseño, no estimada.
- Producción calibrada únicamente contra vino.
- Condensado completamente excluido del optimizador.
- Transferencia fijada como `K(E,T)·QCO2`.
- Validación leave-one-reactor-run-out en los seis procesos."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import sys
from IPython.display import Image, display

ROOT = Path.cwd()
if not (ROOT / 'fermentation_model').exists():
    ROOT = ROOT.parents[1]
sys.path.insert(0, str(ROOT / 'fermentation_model'))
from pilot_2026 import run_aroma_wine_calibrated_complete_capture_2026 as analysis

result = analysis.load_results()
print('Veredicto:', result['gate']['verdict'])
print('Modelo diagnóstico:', result['gate']['diagnostic_best_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Resultados"),
        nbformat.v4.new_code_cell(
            """primary = result['metrics'].query("analysis_policy == 'primary'")
display(primary.round(4))
display(result['comparison'].round(4))
display(Image(filename=analysis.FIGURE_DIR / '02_loro_nrmse_comparison.png'))"""
        ),
        nbformat.v4.new_code_cell(
            """display(result['closure'].round(4))
display(result['capture_audit'])
display(result['parameters'].round(5))"""
        ),
        nbformat.v4.new_code_cell(
            """for name in [
    '03_rco2_temperature_pulse_drivers.png',
    '04_liquid_ethyl_octanoate.png',
    '04_liquid_isoamyl_acetate.png',
    '05_condensate_ethyl_octanoate.png',
    '05_condensate_isoamyl_acetate.png',
    '06_parameter_stability.png',
]:
    display(Image(filename=analysis.FIGURE_DIR / name))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Conclusiones

El ajuste del vino evalúa la estructura biológica sin contaminación del bloque de captura. El condensado evalúa por separado si las pérdidas predichas son compatibles con captura A+B completa. Un fallo de cierre no identifica eficiencia física: rechaza el contrato conjunto bajo la estructura y datos actuales."""
        ),
        nbformat.v4.new_code_cell(
            """assert result['fit_validation']['calibration_domains'].eq('wine').all()
assert result['capture_audit']['fixed_total_capture_fraction'].eq(1.0).all()
assert result['capture_audit']['maximum_absolute_unrecovered_ug'].max() < 1e-8
assert result['fit_validation']['maximum_relative_mass_balance_error'].max() <= 1e-8
print('Notebook ejecutado sin errores; condensado excluido del ajuste.')"""
        ),
    ]
    nbformat.validate(notebook)
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_wine_only_jupyter"
    runtime.mkdir(parents=True, exist_ok=True)
    os.environ["JUPYTER_RUNTIME_DIR"] = str(runtime)
    os.environ["IPYTHONDIR"] = str(runtime / "ipython")
    os.environ["JUPYTER_ALLOW_INSECURE_WRITES"] = "1"
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
    nbformat.validate(notebook)
    nbformat.write(notebook, EXECUTED_NOTEBOOK_PATH)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-notebook", action="store_true")
    parser.add_argument("--reuse-results", action="store_true")
    args = parser.parse_args()
    result = finalize_existing_results() if args.reuse_results else run_analysis()
    create_notebook()
    if not args.skip_notebook:
        execute_notebook()
    print(json.dumps(result["summary"], indent=2, ensure_ascii=False))
    print(NOTEBOOK_PATH)
    if not args.skip_notebook:
        print(EXECUTED_NOTEBOOK_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
