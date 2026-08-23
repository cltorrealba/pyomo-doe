from __future__ import annotations

"""Validate the pilot aroma model under the imposed A+B complete-capture contract.

The available observations do not identify physical condenser efficiency.  This
analysis therefore fixes total A+B capture to one, retains the published
Morakul/Mouret gas-liquid loss driver, and asks whether the resulting joint
production/loss model predicts held-out wine and combined-condensate data.
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

from pilot_2026 import run_aroma_nested_model_validation_2026 as base  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_nested_models as nested  # noqa: E402


CONFIG_PATH = (
    SCRIPT_DIR
    / "adaptive_design"
    / "aroma_assumed_complete_capture_validation_config.json"
)
RESULTS_DIR = (
    SCRIPT_DIR / "results" / "aroma_assumed_complete_capture_validation_2026"
)
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_PATH = (
    SCRIPT_DIR
    / "notebooks"
    / "pilot_2026_aroma_assumed_complete_capture_validation.ipynb"
)
EXECUTED_NOTEBOOK_PATH = NOTEBOOK_PATH.with_name(
    "pilot_2026_aroma_assumed_complete_capture_validation.executed.ipynb"
)

MODEL_LABELS = {
    nested.ASSUMED_COMPLETE: "Captura completa + equilibrio",
    nested.ASSUMED_COMPLETE_RESERVOIR: "Captura completa + reservorio",
    nested.ASSUMED_COMPLETE_DELAYED: "Captura completa + respuesta pospulso",
    nested.ASSUMED_COMPLETE_COMBINED: (
        "Captura completa + pospulso + reservorio"
    ),
}
COLORS = {
    nested.ASSUMED_COMPLETE: "#4C78A8",
    nested.ASSUMED_COMPLETE_RESERVOIR: "#54A24B",
    nested.ASSUMED_COMPLETE_DELAYED: "#E6863B",
    nested.ASSUMED_COMPLETE_COMBINED: "#8E6C8A",
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


def _capture_contract_audit(
    states: pd.DataFrame, mass_summary: pd.DataFrame
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "contract": "combined_A_plus_B_complete_capture",
                "fixed_total_capture_fraction": 1.0,
                "minimum_simulated_capture_fraction": float(
                    states["capture_efficiency_fraction"].min()
                ),
                "maximum_simulated_capture_fraction": float(
                    states["capture_efficiency_fraction"].max()
                ),
                "maximum_absolute_unrecovered_ug": float(
                    mass_summary["emitted_but_unrecovered_ug"].abs().max()
                ),
                "maximum_absolute_mass_closure_error_ug": float(
                    mass_summary["mass_closure_error_ug"].abs().max()
                ),
                "physical_efficiency_identified": False,
                "interpretation": (
                    "eta_AB=1 is imposed by design assumption; it is not estimated "
                    "from the combined condensate observations."
                ),
            }
        ]
    )


def _comparison_to_previous_model(
    current_metrics: pd.DataFrame, current_model: str
) -> pd.DataFrame:
    previous_path = (
        SCRIPT_DIR
        / "results"
        / "aroma_final_nested_validation_2026"
        / "loro_metrics.csv"
    )
    if not previous_path.exists():
        return pd.DataFrame()
    previous = pd.read_csv(previous_path)
    previous = previous[
        previous["analysis_policy"].eq("primary")
        & previous["model_variant"].eq("ethanol_capture_plus_line_reservoir")
    ][["species", "domain", "nrmse_over_mean"]].rename(
        columns={"nrmse_over_mean": "previous_effective_model_nrmse"}
    )
    current = current_metrics[
        current_metrics["analysis_policy"].eq("primary")
        & current_metrics["model_variant"].eq(current_model)
    ][["species", "domain", "nrmse_over_mean"]].rename(
        columns={"nrmse_over_mean": "complete_capture_model_nrmse"}
    )
    comparison = previous.merge(current, on=["species", "domain"], how="inner")
    comparison["relative_nrmse_change_fraction"] = (
        comparison["complete_capture_model_nrmse"]
        / comparison["previous_effective_model_nrmse"]
        - 1.0
    )
    comparison.insert(0, "complete_capture_model", current_model)
    return comparison


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

    primary = base._run_loro(
        config,
        list(config["candidate_models"]),
        forcing_sets,
        wine_tables,
        condensate_tables,
        pulse_lookup,
        analysis_policy="primary",
        exclude_flagged_wine=False,
    )
    primary_metrics = base._score_loro(primary["wine"], primary["condensate"])
    primary_sse = base._runwise_condensate_sse(primary["condensate"])
    provisional_gate, _ = base._gate(
        config,
        primary_metrics,
        primary_sse,
        primary["parameters"],
        primary["validation"],
        primary["simulation"],
    )
    baseline_model = str(
        config["validation_policy"]["comparison_baseline_model"]
    )
    sensitivity_models = list(
        dict.fromkeys(
            [baseline_model, provisional_gate["diagnostic_best_model"]]
        )
    )
    sensitivity = base._run_loro(
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
    metrics = base._score_loro(cv_wine, cv_condensate)
    run_sse = base._runwise_condensate_sse(cv_condensate)
    gate, comparison = base._gate(
        config,
        metrics,
        run_sse,
        cv_parameters,
        cv_validation,
        cv_simulation,
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
    capture_audit = _capture_contract_audit(all_data["states"], mass_summary)
    display_model = gate["selected_model"] or gate["diagnostic_best_model"]
    previous_comparison = _comparison_to_previous_model(metrics, display_model)

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
        "all_data_mass_balance_summary.csv": mass_summary,
        "capture_contract_audit.csv": capture_audit,
        "comparison_to_previous_effective_model.csv": previous_comparison,
    }
    for filename, frame in frames.items():
        frame.to_csv(RESULTS_DIR / filename, index=False)

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

    summary = {
        "verdict": gate["verdict"],
        "selected_model": gate["selected_model"],
        "diagnostic_best_model": gate["diagnostic_best_model"],
        "display_model": display_model,
        "capture_contract": config["capture_contract"],
        "capture_contract_audit": capture_audit.iloc[0].to_dict(),
        "data_quality": {
            "primary_policy": config["data_quality_policy"]["primary"],
            "sensitivity_policy": config["data_quality_policy"]["sensitivity"],
        },
        "interpretation_limits": config["interpretation_limits"],
        "literature_basis": config["literature_basis"],
    }
    _write_json(RESULTS_DIR / "validation_gate.json", gate)
    _write_json(RESULTS_DIR / "analysis_summary.json", summary)
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis_id": config["analysis_id"],
            "runs": list(base.RUNS),
            "species": list(base.SPECIES_LABELS),
            "candidate_models": list(config["candidate_models"]),
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
    return {**frames, "gate": gate, "summary": summary, "figures": figures}


def load_results() -> dict[str, Any]:
    files = {
        "metrics": "loro_metrics.csv",
        "comparison": "validation_gate_comparison.csv",
        "stability": "parameter_stability.csv",
        "parameters": "all_data_parameter_estimates.csv",
        "fit_validation": "all_data_fit_validation.csv",
        "wine": "all_data_wine_predictions.csv",
        "condensate": "all_data_condensate_predictions.csv",
        "states": "all_data_state_trajectories.csv",
        "mass": "all_data_mass_balance_summary.csv",
        "capture_audit": "capture_contract_audit.csv",
        "previous_comparison": "comparison_to_previous_effective_model.csv",
    }
    result = {
        name: pd.read_csv(RESULTS_DIR / filename) for name, filename in files.items()
    }
    result["gate"] = aroma.load_json(RESULTS_DIR / "validation_gate.json")
    result["summary"] = aroma.load_json(RESULTS_DIR / "analysis_summary.json")
    return result


def finalize_existing_results() -> dict[str, Any]:
    """Build gates, figures and metadata from completed calibration CSVs."""

    _configure_base_output()
    config = aroma.load_json(CONFIG_PATH)
    metrics = pd.read_csv(RESULTS_DIR / "loro_metrics.csv")
    run_sse = pd.read_csv(RESULTS_DIR / "loro_runwise_condensate_sse.csv")
    parameters = pd.read_csv(RESULTS_DIR / "loro_parameter_estimates.csv")
    validation = pd.read_csv(RESULTS_DIR / "loro_fit_validation.csv")
    simulation = pd.read_csv(RESULTS_DIR / "loro_simulation_diagnostics.csv")
    gate, comparison = base._gate(
        config, metrics, run_sse, parameters, validation, simulation
    )
    comparison.to_csv(RESULTS_DIR / "validation_gate_comparison.csv", index=False)
    display_model = gate["selected_model"] or gate["diagnostic_best_model"]
    all_data = {
        "parameters": pd.read_csv(RESULTS_DIR / "all_data_parameter_estimates.csv"),
        "wine": pd.read_csv(RESULTS_DIR / "all_data_wine_predictions.csv"),
        "condensate": pd.read_csv(
            RESULTS_DIR / "all_data_condensate_predictions.csv"
        ),
        "states": pd.read_csv(RESULTS_DIR / "all_data_state_trajectories.csv"),
    }
    anomaly_summary = pd.read_csv(
        RESULTS_DIR / "data_quality_sample_summary.csv"
    )
    stability = pd.read_csv(RESULTS_DIR / "parameter_stability.csv")
    mass_summary = pd.read_csv(
        RESULTS_DIR / "all_data_mass_balance_summary.csv"
    )
    capture_audit = pd.read_csv(RESULTS_DIR / "capture_contract_audit.csv")
    previous_comparison = _comparison_to_previous_model(metrics, display_model)
    previous_comparison.to_csv(
        RESULTS_DIR / "comparison_to_previous_effective_model.csv", index=False
    )
    pulses = pd.read_csv(RESULTS_DIR / "nutrient_pulse_schedule.csv")
    pulse_lookup = {
        str(row.batch): float(row.pulse_time_h)
        for row in pulses.itertuples(index=False)
        if str(row.batch) in base.RUNS
    }
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

    summary = {
        "verdict": gate["verdict"],
        "selected_model": gate["selected_model"],
        "diagnostic_best_model": gate["diagnostic_best_model"],
        "display_model": display_model,
        "capture_contract": config["capture_contract"],
        "capture_contract_audit": capture_audit.iloc[0].to_dict(),
        "data_quality": {
            "primary_policy": config["data_quality_policy"]["primary"],
            "sensitivity_policy": config["data_quality_policy"]["sensitivity"],
        },
        "interpretation_limits": config["interpretation_limits"],
        "literature_basis": config["literature_basis"],
    }
    _write_json(RESULTS_DIR / "validation_gate.json", gate)
    _write_json(RESULTS_DIR / "analysis_summary.json", summary)
    output_files = sorted(
        path.name for path in RESULTS_DIR.iterdir() if path.is_file()
    )
    _write_json(
        RESULTS_DIR / "analysis_manifest.json",
        {
            "analysis_id": config["analysis_id"],
            "runs": list(base.RUNS),
            "species": list(base.SPECIES_LABELS),
            "candidate_models": list(config["candidate_models"]),
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
            "outputs": output_files
            + [str(path.relative_to(RESULTS_DIR)) for path in figures],
            "gate": gate,
        },
    )
    return {
        "gate": gate,
        "summary": summary,
        "metrics": metrics,
        "parameters": all_data["parameters"],
        "mass": mass_summary,
        "capture_audit": capture_audit,
        "previous_comparison": previous_comparison,
        "figures": figures,
    }


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
            """# Aromas piloto 2026 — captura A+B completa asumida

## tl;dr

Este notebook prueba explícitamente el contrato de diseño `eta_AB = 1`. La eficiencia física no se estima: toda masa que abandona el estado de línea se asigna al condensado combinado A+B. La validación leave-one-run-out determina si producción, partición de Mouret/Morakul y captura completa son conjuntamente compatibles con los datos."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contexto y métodos

### Supuestos clave

- A y B se observan químicamente como un único `MIX`; sólo sus volúmenes están separados.
- La captura total A+B se fija en 1 por supuesto de diseño, no por medición.
- La volatilización usa el equilibrio `K(E,T)·QCO2` de Morakul/Mouret.
- El condensado participa en el ajuste como medición de toda la pérdida gaseosa bajo este contrato.
- La validación deja fuera una fermentación completa por fold."""
        ),
        nbformat.v4.new_code_cell(
            """from pathlib import Path
import sys
import pandas as pd
from IPython.display import Image, display

ROOT = Path.cwd()
if not (ROOT / 'fermentation_model').exists():
    ROOT = ROOT.parents[1]
sys.path.insert(0, str(ROOT / 'fermentation_model'))
from pilot_2026 import run_aroma_assumed_complete_capture_validation_2026 as analysis

result = analysis.load_results()
print('Resultados:', analysis.RESULTS_DIR.relative_to(ROOT))
print('Veredicto:', result['gate']['verdict'])
print('Modelo diagnóstico:', result['gate']['diagnostic_best_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Resultados"),
        nbformat.v4.new_code_cell(
            """primary_metrics = result['metrics'].query("analysis_policy == 'primary'")
display(primary_metrics.round(4))
display(Image(filename=analysis.FIGURE_DIR / '02_loro_nrmse_comparison.png'))"""
        ),
        nbformat.v4.new_code_cell(
            """display(result['capture_audit'])
display(result['previous_comparison'].round(4))
display(result['parameters'].round(5))"""
        ),
        nbformat.v4.new_code_cell(
            """display(result['mass'].round(4))
assert result['capture_audit']['fixed_total_capture_fraction'].eq(1.0).all()
assert result['capture_audit']['minimum_simulated_capture_fraction'].eq(1.0).all()
assert result['capture_audit']['maximum_simulated_capture_fraction'].eq(1.0).all()
assert result['capture_audit']['maximum_absolute_unrecovered_ug'].max() < 1e-8
assert result['capture_audit']['maximum_absolute_mass_closure_error_ug'].max() < 1e-6"""
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

El veredicto se interpreta sobre el contrato conjunto, no sobre la eficiencia del equipo. Si el condensado no cierra, los datos actuales no permiten asignar la discrepancia exclusivamente a producción, partición o captura completa. El resultado correcto es aceptar o rechazar la compatibilidad del supuesto, manteniéndolo visible."""
        ),
    ]
    nbformat.validate(notebook)
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 3600) -> None:
    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_complete_capture_jupyter"
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
