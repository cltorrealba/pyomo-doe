from __future__ import annotations

"""Validate the condensate observation fraction for pilot-2026 aromas."""

import argparse
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

from pilot_2026 import run_aroma_nested_model_validation_2026 as base  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_calibration as aroma  # noqa: E402
from pilot_2026.adaptive_design import pilot_aroma_nested_models as nested  # noqa: E402


CONFIG_PATH = (
    SCRIPT_DIR
    / "adaptive_design"
    / "aroma_capture_efficiency_validation_config.json"
)
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_capture_efficiency_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_capture_efficiency_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_capture_efficiency_validation.executed.ipynb"
)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=base._json_default)
        + "\n",
        encoding="utf-8",
    )


def _plot_cv(metrics: pd.DataFrame, models: list[str]) -> Path:
    primary = metrics[metrics["analysis_policy"].eq("primary")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for axis, domain in zip(axes, ("wine", "condensate")):
        subset = primary[primary["domain"].eq(domain)]
        x = np.arange(len(base.SPECIES_LABELS))
        for index, model_variant in enumerate(models):
            table = subset[subset["model_variant"].eq(model_variant)].set_index("species")
            values = [
                float(table.loc[species, "nrmse_over_mean"])
                for species in base.SPECIES_LABELS
            ]
            axis.bar(
                x + (index - 0.5) * 0.32,
                values,
                width=0.32,
                color=base.COLORS[model_variant],
                label=base.MODEL_LABELS[model_variant],
            )
        axis.set_xticks(
            x, [base.SPECIES_LABELS[species] for species in base.SPECIES_LABELS]
        )
        axis.set_ylabel("NRMSE / media observada")
        axis.set_title("Vino" if domain == "wine" else "Condensado")
        base._style_axis(axis)
    axes[1].legend(frameon=False)
    fig.suptitle("Impacto fuera de muestra de estimar la fracción de captura")
    return base._save(fig, "02_capture_efficiency_loro_nrmse.png")


def _plot_efficiency(parameters: pd.DataFrame, model_variant: str) -> Path:
    subset = parameters[
        parameters["analysis_policy"].eq("primary")
        & parameters["model_variant"].eq(model_variant)
        & parameters["parameter"].eq("capture_efficiency_fraction")
    ].copy()
    fig, axis = plt.subplots(figsize=(10, 5))
    positions = np.arange(len(base.RUNS))
    width = 0.34
    for species_index, species in enumerate(base.SPECIES_LABELS):
        group = subset[subset["species"].eq(species)].set_index("left_out_run")
        values = [float(group.loc[run, "estimate"]) for run in base.RUNS]
        axis.bar(
            positions + (species_index - 0.5) * width,
            values,
            width=width,
            color=("#2A9D8F" if species_index == 0 else "#E9C46A"),
            label=base.SPECIES_LABELS[species],
        )
    axis.axhline(0.85, color="#6C757D", linestyle="--", linewidth=1.0, label="Nominal EO 0,85")
    axis.axhline(0.88, color="#9A8C98", linestyle=":", linewidth=1.0, label="Nominal IAA 0,88")
    axis.set_yscale("log")
    axis.set_xticks(positions, [f"sin {run}" for run in base.RUNS])
    axis.set_ylabel("Fracción efectiva recuperada")
    axis.set_title("Estabilidad leave-one-run-out de la observación del condensado")
    axis.legend(frameon=False, ncol=2, fontsize=8)
    base._style_axis(axis)
    return base._save(fig, "07_capture_efficiency_by_fold.png")


def run_analysis() -> dict[str, Any]:
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
    sensitivity = base._run_loro(
        config,
        list(config["candidate_models"]),
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
    balance_summary = base.mass_balance_summary(all_data["states"])
    display_model = gate["selected_model"] or gate["diagnostic_best_model"]
    target_model = str(config["candidate_models"][-1])
    volume_context = tables.metadata[
        tables.metadata["experiment_id"].astype(str).isin(base.RUNS)
    ][["experiment_id", "initial_volume_l"]].copy()
    volume_context["nominal_capture_efficiency_ethyl_octanoate"] = 0.85
    volume_context["nominal_capture_efficiency_isoamyl_acetate"] = 0.88

    frames = {
        "data_quality_anomaly_detail.csv": anomaly_detail,
        "data_quality_sample_summary.csv": anomaly_summary,
        "reactor_volume_and_nominal_capture_context.csv": volume_context,
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
        "all_data_mass_balance_summary.csv": balance_summary,
    }
    for filename, frame in frames.items():
        frame.to_csv(RESULTS_DIR / filename, index=False)

    original_figure_dir = base.FIGURE_DIR
    base.FIGURE_DIR = FIGURE_DIR
    try:
        figures = [
            base.plot_anomaly(anomaly_summary),
            _plot_cv(metrics, list(config["candidate_models"])),
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
        figures.append(_plot_efficiency(cv_parameters, target_model))
    finally:
        base.FIGURE_DIR = original_figure_dir

    fitted_efficiency = all_data["parameters"][
        all_data["parameters"]["model_variant"].eq(target_model)
        & all_data["parameters"]["parameter"].eq("capture_efficiency_fraction")
    ][["species", "estimate", "std_log_local", "active_bound"]]
    summary = {
        "verdict": gate["verdict"],
        "selected_model": gate["selected_model"],
        "diagnostic_best_model": gate["diagnostic_best_model"],
        "display_model": display_model,
        "target_model": target_model,
        "reactor_volume_l": sorted(volume_context["initial_volume_l"].astype(float).unique()),
        "nominal_capture_efficiency": {
            "ethyl_octanoate": 0.85,
            "isoamyl_acetate": 0.88,
        },
        "fitted_capture_efficiency_all_data": fitted_efficiency.set_index("species")[
            "estimate"
        ].astype(float).to_dict(),
        "interpretation_limits": config["interpretation_limits"],
        "literature_basis": {
            "production_loss_separation": "Seguinot et al. 2018, DOI 10.1016/j.fm.2018.04.005",
            "gas_liquid_partition": "Morakul et al. 2011, DOI 10.1016/j.procbio.2011.01.034; Mouret et al. 2014, DOI 10.1016/j.foodres.2014.02.044",
            "condensate_is_time_varying_matrix": "Humaj et al. 2024, DOI 10.3390/fermentation10040206",
        },
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
            "post_hoc_status": config["interpretation_limits"]["post_hoc_status"],
            "sources": {
                "config": {
                    "path": str(CONFIG_PATH.relative_to(ROOT_DIR)),
                    "sha256": base._sha256(CONFIG_PATH),
                },
                "model_dataset": str(base.MODEL_DATASET_DIR.relative_to(ROOT_DIR)),
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
        "volume_context": "reactor_volume_and_nominal_capture_context.csv",
        "metrics": "loro_metrics.csv",
        "comparison": "validation_gate_comparison.csv",
        "stability": "parameter_stability.csv",
        "parameters": "all_data_parameter_estimates.csv",
        "fit_validation": "all_data_fit_validation.csv",
        "wine": "all_data_wine_predictions.csv",
        "condensate": "all_data_condensate_predictions.csv",
        "mass_balance": "all_data_mass_balance_summary.csv",
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
            "# Piloto 2026 — validación de la observación del condensado"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto:** `{verdict}`.  
**Modelo seleccionado:** `{selected}`.

La prueba corrige una incompatibilidad de escala: los fermentadores contienen 230 L, mientras que 85–88 % de captura era una hipótesis nominal, no una recuperación medida. Se estima una fracción efectiva por compuesto, manteniendo fijos la producción biológica, el equilibrio gas/líquido y el balance de masa."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Alcance y caveats

- La fracción ajustada representa todo el operador entre aroma emitido y masa recuperada: condensador, línea, trampa y pérdidas no observadas.
- No debe llamarse eficiencia termodinámica del condensador sin calibración con estándar gaseoso.
- El test fue motivado después de observar el fallo de los modelos de retardo y la escala de 230 L; es exploratorio.
- Leave-one-run-out mide transferencia interna retrospectiva. Aún se necesita una corrida prospectiva sellada.
- 26211-P-12 se conserva en el primario y sólo se excluye en sensibilidad."""
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
from pilot_2026 import run_aroma_capture_efficiency_validation_2026 as analysis
result = analysis.load_results() if os.environ.get('PILOT_AROMA_REUSE_RESULTS') == '1' else analysis.run_analysis()
print('Veredicto:', result['gate']['verdict'])
print('Modelo:', result['gate']['selected_model'])
print('Eficiencias ajustadas:', result['summary']['fitted_capture_efficiency_all_data'])"""
        ),
        nbformat.v4.new_markdown_cell("## Escala y calidad de datos"),
        nbformat.v4.new_code_cell(
            """display(result['volume_context'])
display(result['anomaly_summary'].head(10).round(3))
display(Image(filename=analysis.FIGURE_DIR / '01_sample_anomaly_audit.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Validación leave-one-run-out"),
        nbformat.v4.new_code_cell(
            """display(result['metrics'].round(4))
display(result['comparison'].round(4))
display(Image(filename=analysis.FIGURE_DIR / '02_capture_efficiency_loro_nrmse.png'))
display(Image(filename=analysis.FIGURE_DIR / '07_capture_efficiency_by_fold.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Forzantes y curvas"),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / '03_rco2_temperature_pulse_drivers.png'))
for species in analysis.base.SPECIES_LABELS:
    display(Image(filename=analysis.FIGURE_DIR / f'04_liquid_{species}.png'))
    display(Image(filename=analysis.FIGURE_DIR / f'05_condensate_{species}.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Parámetros e identificabilidad"),
        nbformat.v4.new_code_cell(
            """display(result['stability'].round(5))
display(result['parameters'].round(6))
display(result['fit_validation'].round(5))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- Si el modelo pasa, existe una corrección predictiva interna defendible, pero la fracción de captura sigue siendo efectiva y específica de este montaje.
- Si la fracción cambia fuertemente entre folds o cae en límites, la observación de condensado no es transferible con un único valor constante.
- El experimento confirmatorio debe medir simultáneamente concentración líquida, gas de salida y recuperación de trampa con un estándar, registrando volumen y %EtOH del condensado."""
        ),
        nbformat.v4.new_code_cell(
            """assert result['gate']['verdict'] in {'PASS', 'NO_VALID_MODEL'}
assert len(result['figures']) == 8
assert result['fit_validation']['maximum_relative_mass_balance_error'].max() <= 1e-8
print('Notebook ejecutado sin errores; figuras:', len(result['figures']))"""
        ),
    ]
    nbformat.write(notebook, NOTEBOOK_PATH)


def execute_notebook(timeout: int = 1200) -> None:
    import tempfile

    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_capture_jupyter"
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
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--create-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=1200)
    args = parser.parse_args()
    if args.create_only:
        create_notebook()
    else:
        result = run_analysis()
        create_notebook(result["summary"])
        if not args.analysis_only:
            execute_notebook(timeout=args.timeout)
    print(NOTEBOOK_PATH)
    if not args.analysis_only and not args.create_only:
        print(EXECUTED_NOTEBOOK_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
