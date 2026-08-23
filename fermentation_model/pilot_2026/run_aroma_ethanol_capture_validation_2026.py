from __future__ import annotations

"""Validate an ethanol-dependent condensate recovery operator."""

import argparse
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

from pilot_2026 import run_aroma_capture_efficiency_validation_2026 as capture  # noqa: E402


CONFIG_PATH = (
    SCRIPT_DIR
    / "adaptive_design"
    / "aroma_ethanol_capture_validation_config.json"
)
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_ethanol_capture_validation_2026"
FIGURE_DIR = RESULTS_DIR / "figures"
NOTEBOOK_DIR = SCRIPT_DIR / "notebooks"
NOTEBOOK_PATH = NOTEBOOK_DIR / "pilot_2026_aroma_ethanol_capture_validation.ipynb"
EXECUTED_NOTEBOOK_PATH = (
    NOTEBOOK_DIR / "pilot_2026_aroma_ethanol_capture_validation.executed.ipynb"
)


def _configure_capture_module() -> None:
    capture.CONFIG_PATH = CONFIG_PATH
    capture.RESULTS_DIR = RESULTS_DIR
    capture.FIGURE_DIR = FIGURE_DIR
    capture.NOTEBOOK_PATH = NOTEBOOK_PATH
    capture.EXECUTED_NOTEBOOK_PATH = EXECUTED_NOTEBOOK_PATH


def run_analysis() -> dict[str, Any]:
    _configure_capture_module()
    return capture.run_analysis()


def load_results() -> dict[str, Any]:
    _configure_capture_module()
    return capture.load_results()


def create_notebook(summary: dict[str, Any] | None = None) -> None:
    NOTEBOOK_DIR.mkdir(parents=True, exist_ok=True)
    if summary is None and (RESULTS_DIR / "analysis_summary.json").exists():
        summary = capture.aroma.load_json(RESULTS_DIR / "analysis_summary.json")
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
            "# Piloto 2026 — recuperación aromática dependiente de etanol"
        ),
        nbformat.v4.new_markdown_cell(
            f"""## tl;dr

**Veredicto:** `{verdict}`.  
**Modelo seleccionado:** `{selected}`.

Se compara una fracción efectiva constante con `eta(E)=eta50·m^((E−50)/10)`. El cambio de etanol modifica sólo el operador de recuperación del condensado; la producción, la partición gas/líquido y el balance de masa permanecen iguales."""
        ),
        nbformat.v4.new_markdown_cell(
            """## Alcance y supuestos

- El etanol del vino es un proxy de la composición del vapor/condensado; no se midió etanol en cada MIX.
- `eta(E)` agrega un parámetro respecto del modelo de captura constante y se evalúa contra éste, no contra el supuesto nominal 0,85–0,88.
- La muestra 26211-P-12 se conserva en el análisis primario.
- El análisis es post hoc e internamente validado por reactor; requiere confirmación prospectiva y calibración de trampa."""
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
from pilot_2026 import run_aroma_ethanol_capture_validation_2026 as analysis
result = analysis.load_results() if os.environ.get('PILOT_AROMA_REUSE_RESULTS') == '1' else analysis.run_analysis()
print('Veredicto:', result['gate']['verdict'])
print('Comparador:', result['gate']['comparison_baseline_model'])
print('Seleccionado:', result['gate']['selected_model'])"""
        ),
        nbformat.v4.new_markdown_cell("## Validación cruzada"),
        nbformat.v4.new_code_cell(
            """display(result['metrics'].round(4))
display(result['comparison'].round(4))
display(Image(filename=analysis.FIGURE_DIR / '02_capture_efficiency_loro_nrmse.png'))
display(Image(filename=analysis.FIGURE_DIR / '07_capture_efficiency_by_fold.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## CO₂, temperatura, pulso y curvas"),
        nbformat.v4.new_code_cell(
            """display(Image(filename=analysis.FIGURE_DIR / '03_rco2_temperature_pulse_drivers.png'))
for species in analysis.capture.base.SPECIES_LABELS:
    display(Image(filename=analysis.FIGURE_DIR / f'04_liquid_{species}.png'))
    display(Image(filename=analysis.FIGURE_DIR / f'05_condensate_{species}.png'))"""
        ),
        nbformat.v4.new_markdown_cell("## Estabilidad"),
        nbformat.v4.new_code_cell(
            """display(result['stability'].round(5))
display(result['parameters'].round(6))
display(result['fit_validation'].round(5))"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- Una mejora transferible apoyaría que el cambio de composición del condensado es parte del error temporal.
- Un multiplicador en el límite o inestable entre folds implica que etanol sólo está actuando como proxy de tiempo.
- Incluso con `PASS`, la validación confirmatoria debe medir etanol/volumen de cada MIX, concentración en gas de salida y recuperación de un estándar gaseoso."""
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

    runtime = Path(tempfile.gettempdir()) / "pilot26_aroma_ethanol_capture_jupyter"
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
