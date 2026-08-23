from __future__ import annotations

"""Final nested validation for pilot-2026 aroma production and recovery."""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import nbformat
from nbclient import NotebookClient


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_DIR = SCRIPT_DIR.parent
ROOT_DIR = FERMENTATION_DIR.parent
for path in (FERMENTATION_DIR, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from pilot_2026 import run_aroma_decision_aligned_validation_2026 as decision  # noqa: E402


CONFIG_PATH = SCRIPT_DIR / "adaptive_design" / "aroma_final_nested_validation_config.json"
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_final_nested_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_final_nested_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_final_nested_validation.executed.ipynb"
)


def _configure() -> None:
    decision.CONFIG_PATH = CONFIG_PATH
    decision.RESULTS_DIR = RESULTS_DIR
    decision.FIGURE_DIR = FIGURE_DIR
    decision.NOTEBOOK_PATH = NOTEBOOK_PATH
    decision.EXECUTED_NOTEBOOK_PATH = EXECUTED_NOTEBOOK_PATH


def run_analysis() -> dict[str, Any]:
    _configure()
    return decision.run_analysis()


def load_results() -> dict[str, Any]:
    _configure()
    return decision.load_results()


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
            "# Piloto 2026 — validación final anidada del modelo de aromas"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto:** `{verdict}`.  
**Modelo seleccionado:** `{selected}`.

La comparación final contiene: recuperación constante, recuperación dependiente de etanol y esta última con un reservorio de línea que conserva masa. La calibración se pondera de acuerdo con el NRMSE del gate y se valida dejando fuera un reactor completo."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Contrato de interpretación

- Producción: dos rendimientos ligados al consumo de azúcar; no se agrega aquí otro estado de pulso porque no fue transferible en la prueba anterior.
- Pérdida: partición Morakul/Mouret dependiente de etanol y temperatura, forzada por el rCO₂ validado.
- Observación: fracción efectiva de recuperación dependiente de etanol.
- Memoria: inventario de línea de primer orden; su τ sólo es físico si coincide con una medición de residencia.
- Validación: interna retrospectiva; una campaña prospectiva sellada sigue siendo obligatoria."""
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
from pilot_2026 import run_aroma_final_nested_validation_2026 as analysis
result = analysis.load_results() if os.environ.get('PILOT_AROMA_REUSE_RESULTS') == '1' else analysis.run_analysis()
print('Veredicto:', result['gate']['verdict'])
print('Modelo seleccionado:', result['gate']['selected_model'])
print('Mejor diagnóstico:', result['gate']['diagnostic_best_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Calidad y límite de reproducibilidad"),
        nbformat.v4.new_code_cell(
            """display(result['repeatability'].round(4))
display(result['anomaly_summary'].head(10).round(3))
display(Image(filename=analysis.FIGURE_DIR / '01_sample_anomaly_audit.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Validación cruzada y estabilidad"),
        nbformat.v4.new_code_cell(
            """display(result['metrics'].round(4))
display(result['comparison'].round(4))
display(result['stability'].round(5))
display(Image(filename=analysis.FIGURE_DIR / '02_capture_efficiency_loro_nrmse.png'))
display(Image(filename=analysis.FIGURE_DIR / '07_capture_efficiency_by_fold.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Forzantes y curvas finales"),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / '03_rco2_temperature_pulse_drivers.png'))
for species in analysis.decision.ethanol.capture.base.SPECIES_LABELS:
    display(Image(filename=analysis.FIGURE_DIR / f'04_liquid_{species}.png'))
    display(Image(filename=analysis.FIGURE_DIR / f'05_condensate_{species}.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Parámetros finales"),
        nbformat.v4.new_code_cell(
            """display(result['parameters'].round(6))
display(result['fit_validation'].round(6))"""
        ),
        nbformat.v4.new_markdown_cell("## Separación producción–transferencia–captura"),
        nbformat.v4.new_code_cell(
            """selected = result['gate']['selected_model'] or result['gate']['diagnostic_best_model']
display(result['mass_balance'].query('model_variant == @selected').round(4))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- `PASS` indica transferibilidad interna bajo todos los guardrails; no convierte η o τ en parámetros físicos medidos.
- `NO_VALID_MODEL` indica que la variabilidad de condensado y/o la inestabilidad paramétrica impiden una conclusión predictiva con estos datos.
- La confirmación mínima mide simultáneamente vino, gas de salida, volumen y %EtOH de cada MIX, recuperación de estándar gaseoso y duplicados."""
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

    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_final_jupyter"
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
