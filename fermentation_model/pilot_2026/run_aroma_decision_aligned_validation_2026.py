from __future__ import annotations

"""Decision-aligned validation of pilot-2026 aroma recovery models."""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

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

from pilot_2026 import run_aroma_ethanol_capture_validation_2026 as ethanol  # noqa: E402


CONFIG_PATH = (
    SCRIPT_DIR / "adaptive_design" / "aroma_decision_aligned_validation_config.json"
)
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_decision_aligned_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_decision_aligned_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_decision_aligned_validation.executed.ipynb"
)


def _configure() -> None:
    ethanol.CONFIG_PATH = CONFIG_PATH
    ethanol.RESULTS_DIR = RESULTS_DIR
    ethanol.FIGURE_DIR = FIGURE_DIR
    ethanol.NOTEBOOK_PATH = NOTEBOOK_PATH
    ethanol.EXECUTED_NOTEBOOK_PATH = EXECUTED_NOTEBOOK_PATH


def _repeatability_benchmark(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_model = "ethanol_dependent_capture_efficiency"
    observed = predictions[
        predictions["analysis_policy"].eq("primary")
        & predictions["model_variant"].eq(target_model)
        & predictions["status"].eq("observed")
    ].copy()
    observed["mix_number"] = observed["mix_id"].str.extract(r"MIX-(\d+)").astype(int)
    pairs = (("26157", "26158"), ("26210", "26211"))
    rows: list[dict[str, Any]] = []
    for species, group in observed.groupby("species"):
        for first, second in pairs:
            first_table = group[group["experiment_id"].astype(str).eq(first)][
                ["mix_number", "observed_or_upper_bound"]
            ].rename(columns={"observed_or_upper_bound": "first_value"})
            second_table = group[group["experiment_id"].astype(str).eq(second)][
                ["mix_number", "observed_or_upper_bound"]
            ].rename(columns={"observed_or_upper_bound": "second_value"})
            paired = first_table.merge(second_table, on="mix_number", how="inner")
            for row in paired.itertuples(index=False):
                rows.extend(
                    [
                        {
                            "species": species,
                            "target_run": first,
                            "replicate_predictor_run": second,
                            "mix_number": int(row.mix_number),
                            "observed_ug": float(row.first_value),
                            "replicate_predicted_ug": float(row.second_value),
                        },
                        {
                            "species": species,
                            "target_run": second,
                            "replicate_predictor_run": first,
                            "mix_number": int(row.mix_number),
                            "observed_ug": float(row.second_value),
                            "replicate_predicted_ug": float(row.first_value),
                        },
                    ]
                )
    detail = pd.DataFrame(rows)
    summary_rows = []
    for species, group in detail.groupby("species"):
        error = group["replicate_predicted_ug"] - group["observed_ug"]
        rmse = float(np.sqrt(np.mean(error**2)))
        observed_mean = float(group["observed_ug"].mean())
        fold = np.maximum(
            group["replicate_predicted_ug"] / group["observed_ug"],
            group["observed_ug"] / group["replicate_predicted_ug"],
        )
        summary_rows.append(
            {
                "species": species,
                "n_directional_pairs": int(len(group)),
                "rmse_ug": rmse,
                "observed_mean_ug": observed_mean,
                "replicate_as_predictor_nrmse": rmse / observed_mean,
                "median_absolute_fold_difference": float(np.median(fold)),
            }
        )
    return detail, pd.DataFrame(summary_rows)


def run_analysis() -> dict[str, Any]:
    _configure()
    result = ethanol.run_analysis()
    detail, repeatability = _repeatability_benchmark(
        result["loro_condensate_predictions.csv"]
    )
    detail.to_csv(RESULTS_DIR / "matched_replicate_benchmark_detail.csv", index=False)
    repeatability.to_csv(
        RESULTS_DIR / "matched_replicate_benchmark_summary.csv", index=False
    )
    summary_path = RESULTS_DIR / "analysis_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["matched_replicate_benchmark"] = repeatability.set_index("species")[
        "replicate_as_predictor_nrmse"
    ].astype(float).to_dict()
    ethanol.capture._write_json(summary_path, summary)
    manifest_path = RESULTS_DIR / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for filename in (
        "matched_replicate_benchmark_detail.csv",
        "matched_replicate_benchmark_summary.csv",
    ):
        if filename not in manifest["outputs"]:
            manifest["outputs"].append(filename)
    manifest["outputs"] = sorted(manifest["outputs"])
    manifest["matched_replicate_benchmark"] = summary[
        "matched_replicate_benchmark"
    ]
    ethanol.capture._write_json(manifest_path, manifest)
    result["summary"] = summary
    result["repeatability_detail"] = detail
    result["repeatability"] = repeatability
    return result


def load_results() -> dict[str, Any]:
    _configure()
    result = ethanol.load_results()
    result["repeatability_detail"] = pd.read_csv(
        RESULTS_DIR / "matched_replicate_benchmark_detail.csv"
    )
    result["repeatability"] = pd.read_csv(
        RESULTS_DIR / "matched_replicate_benchmark_summary.csv"
    )
    result["summary"] = json.loads(
        (RESULTS_DIR / "analysis_summary.json").read_text(encoding="utf-8")
    )
    return result


def create_notebook(summary: dict[str, Any] | None = None) -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    if summary is None and (RESULTS_DIR / "analysis_summary.json").exists():
        summary = json.loads(
            (RESULTS_DIR / "analysis_summary.json").read_text(encoding="utf-8")
        )
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
            "# Piloto 2026 — validación de aromas alineada con el criterio de decisión"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto:** `{verdict}`.  
**Modelo seleccionado:** `{selected}`.

Esta iteración no agrega parámetros: compara captura constante y dependiente de etanol usando una ponderación homoscedástica por dominio. El objetivo de calibración queda así alineado con el NRMSE usado en el gate."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Caveats predefinidos

- La ponderación por media del dominio es una función de decisión, no una estimación de error analítico.
- La eficiencia sigue siendo un operador efectivo de recuperación; no una propiedad termodinámica medida.
- El benchmark de réplicas cuantifica el piso experimental de predicción.
- El análisis es retrospectivo y post hoc; cualquier `PASS` requiere una campaña confirmatoria."""
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
from pilot_2026 import run_aroma_decision_aligned_validation_2026 as analysis
result = analysis.load_results() if os.environ.get('PILOT_AROMA_REUSE_RESULTS') == '1' else analysis.run_analysis()
print('Veredicto:', result['gate']['verdict'])
print('Seleccionado:', result['gate']['selected_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Reproducibilidad experimental"),
        nbformat.v4.new_code_cell(
            """display(result['repeatability'].round(4))
display(result['repeatability_detail'].round(3))
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
for species in analysis.ethanol.capture.base.SPECIES_LABELS:
    display(Image(filename=analysis.FIGURE_DIR / f'04_liquid_{species}.png'))
    display(Image(filename=analysis.FIGURE_DIR / f'05_condensate_{species}.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Estabilidad paramétrica"),
        nbformat.v4.new_code_cell(
            """display(result['stability'].round(5))
display(result['parameters'].round(6))
display(result['fit_validation'].round(5))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- El gate compara la extensión dependiente de etanol contra la captura constante bajo la misma función de decisión.
- Si el error se aproxima al benchmark de réplicas, más complejidad no puede justificarse sin mejorar el protocolo de medición.
- El ensayo confirmatorio debe incluir estándar gaseoso, concentración en gas de salida, volumen/%EtOH por MIX y duplicados de condensado."""
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

    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_decision_jupyter"
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
