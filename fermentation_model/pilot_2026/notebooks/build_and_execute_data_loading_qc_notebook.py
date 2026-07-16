from __future__ import annotations

import base64
import contextlib
import io
import json
import os
import traceback
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE / "pilot_2026_data_loading_qc.ipynb"


def markdown(source: str) -> dict[str, object]:
    return {"cell_type": "markdown", "metadata": {}, "source": source.strip() + "\n"}


def code(source: str) -> dict[str, object]:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.strip() + "\n",
    }


def build_notebook() -> dict[str, object]:
    cells = [
        markdown(
            """
# Piloto 2026 — verificación de carga y calidad de datos

## tl;dr

Este notebook verifica la integración reproducible de nueve fermentaciones piloto. La carga corregida conserva 39 muestras de vino para GC —seis baselines iniciales y 33 muestras emparejadas con 33 MIX—, trata NQ y bajo-LOQ como censura izquierda, representa cada MIX como una acumulación por intervalo y aplica ×1.000 sólo a los condensados. Los tres perfiles de CO₂ del Lote 1 se muestran para control visual, pero quedan fuera de la calibración. Ocho artefactos de etanol `−0,3 % v/v` fueron corregidos a cero conservando el valor original.

**Estado esperado:** carga reproducible y apta para revisión del responsable experimental; calibración todavía condicional a esa revisión y a documentar la referencia normal de MassView.
"""
        ),
        markdown(
            """
## Contexto y métodos

### Supuestos confirmados

- Ventana de proceso: primera a última muestra primaria.
- Temperatura ejecutada: `Sonda1`; el setpoint es la orden del controlador.
- Protocolo C: 16→18→21 °C; si el vino se secó antes, el tramo no ejecutado se excluye.
- MassView: Ln/min, sin segunda normalización.
- Ceros de CO₂: válidos al inicio/final e inválidos entre la primera y última señal positiva.
- Lote 3: misma operación relativa que Lote 1, desplazada a su fecha de inicio.
- MIX: concentración del vial diluido ×1.000 y masa acumulada `C_MIX × (V_A + V_B)`.
"""
        ),
        code(
            """
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-whitegrid")
pd.set_option("display.max_columns", 30)
pd.set_option("display.width", 180)

def find_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "fermentation_model" / "pilot_2026").exists():
            return candidate
    raise FileNotFoundError("No se encontró la raíz del repositorio pyomo-doe")

REPO = find_repo_root(Path.cwd().resolve())
RESULTS = REPO / "fermentation_model" / "pilot_2026" / "results" / "data_integration_2026"
print(f"Repositorio: {REPO}")
print(f"Resultados:  {RESULTS}")
"""
        ),
        markdown("## Datos"),
        code(
            """
qc = json.loads((RESULTS / "qc_summary.json").read_text(encoding="utf-8"))
windows = pd.read_csv(RESULTS / "process_windows.csv", parse_dates=["sampling_start", "sampling_end", "active_end"])
primary = pd.read_csv(RESULTS / "primary_results_qc.csv", parse_dates=["timestamp", "sampling_start", "active_end"])
temperature = pd.read_csv(RESULTS / "temperature_controller_qc.csv", parse_dates=["timestamp"])
temperature_summary = pd.read_csv(RESULTS / "temperature_summary.csv")
co2 = pd.read_csv(RESULTS / "co2_minute_qc.csv", parse_dates=["minute"], low_memory=False)
co2_summary = pd.read_csv(RESULTS / "co2_run_qc_summary.csv")
events = pd.read_csv(RESULTS / "operational_events_qc.csv", parse_dates=["timestamp", "sampling_start", "sampling_end", "active_end"])
gc = pd.read_csv(RESULTS / "gc_results_long_qc.csv", parse_dates=["sample_timestamp"])
gc_pairs = pd.read_csv(RESULTS / "gc_mix_wine_pairs_qc.csv", parse_dates=["timestamp", "capture_interval_start", "capture_interval_end"])
gc_pairs_long = pd.read_csv(RESULTS / "gc_mix_wine_pairs_long_qc.csv", parse_dates=["timestamp", "capture_interval_start", "capture_interval_end"])

tables = {
    "process_windows": windows,
    "primary_results": primary,
    "temperature": temperature,
    "co2_minute": co2,
    "operational_events": events,
    "gc_results_long": gc,
    "gc_pairs": gc_pairs,
    "gc_pairs_long": gc_pairs_long,
}
print(pd.DataFrame({"rows": {name: len(table) for name, table in tables.items()},
                    "columns": {name: len(table.columns) for name, table in tables.items()}}).to_string())
"""
        ),
        code(
            """
expected_checks = pd.DataFrame(qc["expected_value_checks"]).T
print("Controles de ingestión reproducibles")
print(expected_checks.to_string())
assert expected_checks["pass"].all()
assert len(windows) == 9 and windows["experiment_id"].nunique() == 9
assert qc["gc"]["gc_mix_wine_pair_coverage"] == 1.0
assert qc["primary_results"]["ethanol_negative_formula_artifacts_set_to_zero"] == 8
print("\\nPASS: conteos históricos y reglas nuevas reconciliados.")
"""
        ),
        code(
            """
window_view = windows[["experiment_id", "lot", "reactor", "protocol", "sampling_start", "active_end", "sampling_end", "primary_samples"]].copy()
window_view["duration_h"] = (window_view["sampling_end"] - window_view["sampling_start"]).dt.total_seconds() / 3600
window_view["active_h"] = (window_view["active_end"] - window_view["sampling_start"]).dt.total_seconds() / 3600
print(window_view.round(2).to_string(index=False))
"""
        ),
        markdown("## Resultados"),
        markdown("### 1. Resultados primarios y corrección de etanol"),
        code(
            """
core_columns = [
    "Densidad", "Brix", "Oculyze Concentration", "Oculyze Viability",
    "Y15 Glucosa", "Y15 Fructosa", "Y15 PAN", "Y15 Amonio",
    "Y15 Glicerol", "Y15 Acetico", "Y15 Piruvic", "Y15 Acetaldehido",
    "Cf Alcolyzer Real (% v/v)", "pH", "DO",
]
coverage = primary.groupby("experiment_id")[core_columns].apply(lambda frame: frame.notna().mean()).T
fig, ax = plt.subplots(figsize=(13, 6))
image = ax.imshow(coverage.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="viridis")
ax.set_yticks(range(len(coverage.index)), coverage.index)
ax.set_xticks(range(len(coverage.columns)), coverage.columns, rotation=45, ha="right")
ax.set_title("Cobertura de resultados primarios por fermentación")
fig.colorbar(image, ax=ax, label="Fracción no nula")
fig.tight_layout()
plt.show()

corrected = primary[primary["ethanol_correction"].ne("none")][
    ["Código muestra", "timestamp", "ethanol_real_original_percent_vv", "Cf Alcolyzer Real (% v/v)", "process_phase"]
]
print("Correcciones de etanol (original → procesado):")
print(corrected.to_string(index=False))
"""
        ),
        markdown("### 2. Temperatura medida frente a setpoint"),
        code(
            """
window_lookup = windows.set_index("experiment_id")
runs = windows["experiment_id"].astype(str).tolist()
fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharey=True)
for ax, run in zip(axes.flat, runs):
    frame = temperature[(temperature["experiment_id"].astype(str).eq(run)) & temperature["inside_sampling_window"]].copy()
    t0 = window_lookup.loc[int(run) if windows["experiment_id"].dtype.kind in "iu" else run, "sampling_start"]
    hours = (frame["timestamp"] - t0).dt.total_seconds() / 3600
    ax.plot(hours, frame["sensor_1_c"], lw=1.2, label="Sonda1")
    ax.step(hours, frame["setpoint_c"], where="post", lw=1.1, label="Setpoint")
    active_h = (pd.Timestamp(window_lookup.loc[int(run) if windows["experiment_id"].dtype.kind in "iu" else run, "active_end"]) - t0).total_seconds() / 3600
    ax.axvline(active_h, color="black", ls="--", lw=0.8, alpha=0.7)
    ax.set_title(f"{run} · protocolo {window_lookup.loc[int(run) if windows['experiment_id'].dtype.kind in 'iu' else run, 'protocol']}")
    ax.set_xlabel("h desde primera muestra")
    ax.set_ylabel("°C")
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2)
fig.suptitle("Temperatura ejecutada; línea discontinua = término activo", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.925))
plt.show()
print(temperature_summary[["experiment_id", "temperature_rows_in_sampling_window", "mean_absolute_setpoint_error_c", "p95_absolute_setpoint_error_c", "active_end_source"]].round(3).to_string(index=False))
"""
        ),
        markdown("### 3. CO₂: calidad visual y alcance de calibración"),
        code(
            """
print(co2_summary[["experiment_id", "co2_model_include", "invalid_intermediate_zero_seconds", "valid_observed_edge_zero_seconds", "postprocess_observed_seconds", "co2_model_exclusion_reason"]].to_string(index=False))
fig, axes = plt.subplots(3, 3, figsize=(16, 11))
for ax, run in zip(axes.flat, runs):
    frame = co2[co2["experiment_id"].astype(str).eq(run)].copy()
    t0 = pd.Timestamp(window_lookup.loc[int(run) if windows["experiment_id"].dtype.kind in "iu" else run, "sampling_start"])
    hours = (frame["minute"] - t0).dt.total_seconds() / 3600
    ax.plot(hours, frame["mean_flow_ln_min_observed"], lw=0.8, color="#4c78a8")
    invalid = frame["invalid_intermediate_zero"].fillna(0).gt(0)
    ax.scatter(hours[invalid], np.zeros(invalid.sum()), s=3, color="#e45756")
    included = bool(frame["co2_model_include"].iloc[0])
    ax.set_title(f"{run} · {'MODELO' if included else 'QC SOLAMENTE'}")
    ax.set_xlabel("h desde primera muestra")
    ax.set_ylabel("CO₂ (Ln/min)")
fig.suptitle("Perfiles de CO₂; rojo = minutos con ceros intermedios inválidos", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97))
plt.show()
"""
        ),
        markdown("### 4. Eventos operacionales"),
        code(
            """
event_counts = events.groupby(["experiment_id", "event_origin"]).size().unstack(fill_value=0)
print(event_counts.to_string())
print("\\nEventos excluidos del proceso activo:")
print(events[~events["calibration_include"]][["experiment_id", "timestamp", "Tipo evento", "Producto/Nutriente", "event_origin"]].to_string(index=False))
print("\\nNutrición derivada de 116 g orgánico + 46 g FDA por pulso:")
print(events[events["dose_parse_status"].eq("parsed_springferm_organic_plus_fda")][["experiment_id", "relative_time_h", "organic_product_g", "dap_product_g", "yan_added_mg_l", "event_timestamp_status"]].round(3).to_string(index=False))

fig, ax = plt.subplots(figsize=(13, 5))
event_types = {name: idx for idx, name in enumerate(sorted(events["Tipo evento"].dropna().unique()))}
for run_index, run in enumerate(runs):
    frame = events[events["experiment_id"].astype(str).eq(run)]
    for _, row in frame.iterrows():
        color = plt.cm.tab10(event_types[row["Tipo evento"]] % 10)
        ax.scatter(row["relative_time_h"], run_index, color=color, s=55, marker="o" if row["calibration_include"] else "x")
ax.set_yticks(range(len(runs)), runs)
ax.set_xlabel("h desde primera muestra")
ax.set_ylabel("Fermentación")
ax.set_title("Eventos operacionales; × = fuera del proceso activo")
fig.tight_layout()
plt.show()
"""
        ),
        markdown("### 5. GC: baseline, emparejamiento, censura y masa capturada"),
        code(
            """
gc_summary = pd.DataFrame({
    "samples": gc.groupby("sample_type")["sample_id"].nunique(),
    "analyte_results": gc.groupby("sample_type").size(),
    "numeric": gc.groupby("sample_type")["vial_concentration_ug_l"].apply(lambda x: x.notna().sum()),
    "left_censored": gc.groupby("sample_type")["model_observation_type"].apply(lambda x: x.eq("left_censored").sum()),
})
print(gc_summary.to_string())
print(f"\\nPares MIX–vino: {len(gc_pairs)}; cobertura: {gc_pairs['wine_pair_available'].mean():.1%}")
print(f"Baselines de vino sin MIX: {qc['gc']['gc_initial_wine_samples_without_mix']}")
print(gc_pairs[["experiment_id", "mix_id", "paired_wine_sample_id", "timestamp", "capture_interval_h", "total_condensate_ml"]].head(12).round(2).to_string(index=False))

status = gc.groupby(["sample_type", "analyte", "result_status"]).size().rename("n").reset_index()
fig, axes = plt.subplots(1, 2, figsize=(16, 5.5), sharey=True)
for ax, sample_type in zip(axes, ["wine_mef", "condensate_mix"]):
    pivot = status[status["sample_type"].eq(sample_type)].pivot(index="analyte", columns="result_status", values="n").fillna(0)
    pivot = pivot.reindex(columns=["quantified", "below_loq", "NQ", "missing"], fill_value=0)
    pivot.plot.bar(stacked=True, ax=ax, legend=False)
    ax.set_title(sample_type)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=25)
axes[0].set_ylabel("Resultados")
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4)
fig.suptitle("Estado analítico por matriz", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.89))
plt.show()
"""
        ),
        code(
            """
analytes = sorted(gc_pairs_long["analyte"].unique())
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ax, analyte in zip(axes.flat, analytes):
    frame = gc_pairs_long[
        gc_pairs_long["analyte"].eq(analyte)
        & gc_pairs_long["wine_concentration_ug_l"].notna()
        & gc_pairs_long["captured_mass_ug"].notna()
        & gc_pairs_long["wine_concentration_ug_l"].gt(0)
        & gc_pairs_long["captured_mass_ug"].gt(0)
    ]
    if len(frame):
        ax.scatter(frame["wine_concentration_ug_l"], frame["captured_mass_ug"], alpha=0.75)
        ax.set_xscale("log")
        ax.set_yscale("log")
    else:
        ax.text(0.5, 0.5, "Sin pares cuantificados", ha="center", va="center", transform=ax.transAxes)
    ax.set_title(analyte)
    ax.set_xlabel("Concentración en vino (µg/L)")
    ax.set_ylabel("Masa capturada intervalo (µg)")
fig.suptitle("Correspondencia exploratoria vino–condensado (sólo pares cuantificados)", y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97))
plt.show()
"""
        ),
        markdown("### 6. Controles finales automáticos"),
        code(
            """
checks = {
    "9 fermentaciones": len(windows) == 9,
    "todas las fechas MIX Lote 2 están en abril": gc_pairs.loc[gc_pairs["experiment_id"].isin([26157, 26158, 26159]), "timestamp"].dt.month.eq(4).all(),
    "39 vinos GC": gc.loc[gc["sample_type"].eq("wine_mef"), "sample_id"].nunique() == 39,
    "33 MIX GC": gc.loc[gc["sample_type"].eq("condensate_mix"), "sample_id"].nunique() == 33,
    "33 pares MIX–vino": len(gc_pairs) == 33 and gc_pairs["wine_pair_available"].all(),
    "6 baselines sin MIX": qc["gc"]["gc_initial_wine_samples_without_mix"] == 6,
    "8 etanoles corregidos a cero": len(corrected) == 8 and corrected["Cf Alcolyzer Real (% v/v)"].eq(0).all(),
    "15 eventos reconstruidos Lote 3": events["event_origin"].eq("reconstructed_from_lot1_relative_schedule").sum() == 15,
    "CO2 Lote 1 fuera del modelo": not co2_summary.loc[co2_summary["experiment_id"].isin([26134, 26135, 26136]), "co2_model_include"].any(),
    "NQ/bajo LOQ permanecen censurados": gc.loc[gc["result_status"].isin(["NQ", "below_loq"]), "model_observation_type"].eq("left_censored").all(),
}
check_table = pd.Series(checks, name="PASS").to_frame()
print(check_table.to_string())
assert check_table["PASS"].all(), check_table[~check_table["PASS"]]
print("\\nVEREDICTO DEL NOTEBOOK: PASS para carga y reglas del responsable experimental.")
"""
        ),
        markdown(
            """
## Takeaways

- La carga es internamente coherente y reproduce todos los conteos auditados.
- Las fechas GC del Lote 2 quedan en abril y los 33 MIX tienen correspondencia 1:1 con vino; las seis muestras iniciales de vino permanecen como baseline.
- Las observaciones censuradas siguen disponibles para una likelihood censurada y no fueron convertidas a cero.
- El CO₂ de Lote 1 puede inspeccionarse, pero no se entregará como observación de calibración.
- La temperatura de calibración será `Sonda1`; los perfiles se diseñarán como setpoints y deberán pasar por un modelo del controlador.
- Los eventos del Lote 3 están reconstruidos de forma trazable desde Lote 1. El segundo pulso conserva la advertencia de que su hora es un proxy ligado a densidad.
- La dosis `116 + 46 g` equivale a 90,435 mg/L YAN por pulso a 230 L con los factores confirmados; falta reconciliarla con el nominal de 80 mg/L antes de calibrar nutrición.
- Siguiente gate: revisión visual del responsable y construcción del adaptador de estas tablas al modelo jerárquico previo al MBDoE adaptativo.
"""
        ),
    ]
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def execute_notebook(notebook: dict[str, object]) -> None:
    namespace: dict[str, object] = {"__name__": "__main__"}
    execution_count = 0
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        execution_count += 1
        outputs: list[dict[str, object]] = []
        stdout = io.StringIO()
        stderr = io.StringIO()
        images: list[str] = []
        original_show = plt.show

        def capture_show(*_args: object, **_kwargs: object) -> None:
            for number in plt.get_fignums():
                figure = plt.figure(number)
                buffer = io.BytesIO()
                figure.savefig(buffer, format="png", dpi=125, bbox_inches="tight")
                images.append(base64.b64encode(buffer.getvalue()).decode("ascii"))
            plt.close("all")

        plt.show = capture_show
        try:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exec(compile(cell["source"], str(NOTEBOOK), "exec"), namespace)
            if plt.get_fignums():
                capture_show()
        except Exception as exc:
            trace = traceback.format_exc().splitlines()
            cell["outputs"] = [
                {
                    "output_type": "error",
                    "ename": type(exc).__name__,
                    "evalue": str(exc),
                    "traceback": trace,
                }
            ]
            cell["execution_count"] = execution_count
            raise
        finally:
            plt.show = original_show
        if stdout.getvalue():
            outputs.append(
                {"output_type": "stream", "name": "stdout", "text": stdout.getvalue()}
            )
        if stderr.getvalue():
            outputs.append(
                {"output_type": "stream", "name": "stderr", "text": stderr.getvalue()}
            )
        outputs.extend(
            {
                "output_type": "display_data",
                "metadata": {},
                "data": {"image/png": image},
            }
            for image in images
        )
        cell["outputs"] = outputs
        cell["execution_count"] = execution_count


def validate_notebook(notebook: dict[str, object]) -> None:
    assert notebook["nbformat"] == 4
    assert isinstance(notebook["cells"], list) and notebook["cells"]
    for cell in notebook["cells"]:
        assert cell["cell_type"] in {"markdown", "code"}
        assert isinstance(cell["source"], str)
        if cell["cell_type"] == "code":
            assert cell["execution_count"] is not None
            assert isinstance(cell["outputs"], list)


def main() -> None:
    notebook = build_notebook()
    execute_notebook(notebook)
    validate_notebook(notebook)
    NOTEBOOK.write_text(json.dumps(notebook, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"Built and executed {NOTEBOOK}")


if __name__ == "__main__":
    main()
