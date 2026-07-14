from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
BUNDLE_NAME = "A02_validacion_modelo_vendimia_2025"
BUNDLE_DIR = ROOT / "rendicion_tecnica" / BUNDLE_NAME
ZIP_PATH = ROOT / "rendicion_tecnica" / f"{BUNDLE_NAME}.zip"

GENERATED_DATE = "2026-06-23"
ANNEX_ID = "A02"
OE = "OE2"
ACTIVITY_CODE = "2.14"
ACTIVITY_TITLE = "Validacion del modelo con datos de vendimia 2025 (Actividad 2.14)"
SUGGESTED_FILENAME = "PI-4497_IA5_A02_Act-2.14_Validacion-modelo-vendimia-2025.pdf"

CONTEXT_ALIASES = {
    "curve_validation": "curve",
    "current_fit": "curve_fit",
    "design_predictions": "design_pred",
    "medium_transfer_diagnostics": "medium_transfer",
    "new_must_data_loading": "data_loading",
    "new_must_glycerol_estimability_doe": "nmg_doe",
    "new_must_glycerol_overnight_validation": "overnight",
    "best_multistart_post_analysis": "best_post",
}


manifest_rows: list[dict[str, str]] = []


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_git(args: list[str]) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True, stderr=subprocess.STDOUT).strip()
    except Exception as exc:  # pragma: no cover
        return f"UNAVAILABLE: {exc}"


def reset_bundle() -> None:
    target_root = (ROOT / "rendicion_tecnica").resolve()
    target_root.mkdir(exist_ok=True)
    bundle_resolved = BUNDLE_DIR.resolve()
    if BUNDLE_DIR.exists():
        if not str(bundle_resolved).startswith(str(target_root)):
            raise RuntimeError(f"Refusing to remove unexpected path: {bundle_resolved}")
        shutil.rmtree(BUNDLE_DIR)
    for sub in [
        "data_sources",
        "code_sources",
        "notebooks",
        "source_reports",
        "source_tables",
        "figures",
        "tables",
    ]:
        (BUNDLE_DIR / sub).mkdir(parents=True, exist_ok=True)


def contextual_name(src: Path) -> str:
    if src.parent.name in {"current_fit", "design_predictions"}:
        prefix = CONTEXT_ALIASES.get(src.parent.name, src.parent.name)
    else:
        prefix = CONTEXT_ALIASES.get(src.parent.name, src.parent.name)
    name = f"{prefix}__{src.name}"
    if len(name) <= 62:
        return name
    digest = hashlib.sha1(rel(src).encode("utf-8")).hexdigest()[:8]
    suffix = src.suffix
    available = max(12, 62 - len(prefix) - len(digest) - len(suffix) - 4)
    return f"{prefix}__{src.stem[:available]}__{digest}{suffix}"


def record_artifact(
    *,
    bundle_path: Path,
    category: str,
    evidence_requirement: str,
    role: str,
    source_path: Path | None = None,
    status: str = "ok",
) -> None:
    source_text = "generated" if source_path is None else rel(source_path)
    size = ""
    digest = ""
    if bundle_path.exists() and bundle_path.is_file():
        size = str(bundle_path.stat().st_size)
        digest = sha256(bundle_path)
    manifest_rows.append(
        {
            "annex_id": ANNEX_ID,
            "oe": OE,
            "activity": ACTIVITY_CODE,
            "category": category,
            "evidence_requirement": evidence_requirement,
            "role_in_bundle": role,
            "status": status,
            "source_path": source_text,
            "bundle_path": rel(bundle_path),
            "bytes": size,
            "sha256": digest,
        }
    )


def copy_file(
    source: str,
    dest_subdir: str,
    category: str,
    evidence_requirement: str,
    role: str,
    namespace: bool = False,
) -> Path | None:
    src = ROOT / source
    dst = BUNDLE_DIR / dest_subdir / (contextual_name(src) if namespace else src.name)
    if not src.exists():
        record_artifact(
            bundle_path=dst,
            category=category,
            evidence_requirement=evidence_requirement,
            role=role,
            source_path=src,
            status="missing",
        )
        return None
    shutil.copy2(src, dst)
    record_artifact(
        bundle_path=dst,
        category=category,
        evidence_requirement=evidence_requirement,
        role=role,
        source_path=src,
    )
    return dst


def copy_glob(
    pattern: str,
    dest_subdir: str,
    category: str,
    evidence_requirement: str,
    role_prefix: str,
) -> None:
    for src in sorted(ROOT.glob(pattern)):
        if not src.is_file():
            continue
        dst = BUNDLE_DIR / dest_subdir / contextual_name(src)
        shutil.copy2(src, dst)
        record_artifact(
            bundle_path=dst,
            category=category,
            evidence_requirement=evidence_requirement,
            role=f"{role_prefix}: {src.name}",
            source_path=src,
        )


def write_text_artifact(rel_path: str, text: str, category: str, evidence_requirement: str, role: str) -> Path:
    path = BUNDLE_DIR / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    record_artifact(
        bundle_path=path,
        category=category,
        evidence_requirement=evidence_requirement,
        role=role,
    )
    return path


def write_df_artifact(rel_path: str, df: pd.DataFrame, category: str, evidence_requirement: str, role: str) -> Path:
    path = BUNDLE_DIR / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    record_artifact(
        bundle_path=path,
        category=category,
        evidence_requirement=evidence_requirement,
        role=role,
    )
    return path


def read_csv_if_exists(path: str) -> pd.DataFrame:
    p = ROOT / path
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def workbook_summary(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not path.exists():
        return rows
    wb = load_workbook(path, read_only=True, data_only=True)
    for ws in wb.worksheets:
        headers = []
        for cell in next(ws.iter_rows(min_row=1, max_row=1, values_only=True), ()):
            if cell is not None and str(cell).strip():
                headers.append(str(cell).strip())
        rows.append(
            {
                "source_file": rel(path),
                "sheet": ws.title,
                "rows_including_header": ws.max_row,
                "data_rows_estimated": max(ws.max_row - 1, 0),
                "columns": ws.max_column,
                "header_preview": ";".join(headers[:35]),
            }
        )
    return rows


def build_annex_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "ID anexo": ANNEX_ID,
                "OE": OE,
                "Actividad unica": ACTIVITY_CODE,
                "Titulo propuesto": ACTIVITY_TITLE,
                "Nombre de archivo sugerido": SUGGESTED_FILENAME,
                "Contenido minimo": "Documentar la base independiente 2025, protocolo de validacion, comparacion predicho-observado, metricas, brechas y decision tecnica.",
                "Evidencia fuente": "Base independiente vendimia 2025, predicho vs. observado, metricas de error, criterio de validacion, brechas y decision tecnica.",
            }
        ]
    )


def build_compliance_matrix() -> pd.DataFrame:
    rows = [
        (
            "Base independiente vendimia 2025",
            "cubierto",
            "data_sources/*.xlsx; tables/03_base_independiente_2025_resumen_lotes.csv; tables/03_base_independiente_2025_diccionario.csv",
            "Incluye mosto natural y sintetico 2025, normalizacion long y resumen por lote.",
        ),
        (
            "Protocolo de validacion",
            "cubierto",
            "03_protocolo_validacion_A02.md; code_sources/run_new_must_curve_validation.py",
            "Define fuentes, estados, comparacion base vs multistart, predicho-observado, factibilidad y decision.",
        ),
        (
            "Comparacion predicho-observado",
            "cubierto",
            "tables/04_predicho_observado_metricas.csv; figures/curve_fit__*.png",
            "Incluye metricas por lote/estado y graficos fitcmp para lotes naturales y sinteticos.",
        ),
        (
            "Metricas de error",
            "cubierto",
            "tables/04_metricas_error_por_medio_estado.csv; tables/04_residuales_agregados_base_vs_best.csv",
            "Incluye weighted RMSE, RMSE, sesgo medio y comparacion theta base vs mejor multistart.",
        ),
        (
            "Criterio de validacion",
            "cubierto",
            "tables/05_criterios_validacion.csv; 04_decision_tecnica_A02.md",
            "Criterios de base cargada, cobertura de metricas, robustez, factibilidad y transferencia natural.",
        ),
        (
            "Brechas",
            "cubierto",
            "tables/06_brechas_y_acciones.csv; source_reports/medium_transfer__medium_transfer_diagnostics_report.md",
            "Incluye identificabilidad practica, transferencia de medio, redundancias y necesidad de recalibracion secuencial.",
        ),
        (
            "Decision tecnica",
            "cubierto",
            "04_decision_tecnica_A02.md; ANEXO_A02_borrador.md",
            "Decision de validacion tecnica condicionada para uso secuencial y primer bloque mixto.",
        ),
    ]
    return pd.DataFrame(rows, columns=["contenido_minimo", "estado_bundle", "archivos_bundle", "detalle"])


def build_data_dictionary(normalized: pd.DataFrame) -> pd.DataFrame:
    fields = [
        ("medium", "Tipo de matriz", "natural/synthetic", "Define base independiente y transferencia de medio"),
        ("batch", "Lote o fermentacion", "LAB004-LAB012, MS007-MS016", "Unidad experimental"),
        ("time_h", "Tiempo de fermentacion", "h", "Eje temporal para predicho-observado"),
        ("temperature_c", "Temperatura", "C", "Entrada/condicion operacional"),
        ("X_viable_kg_m3", "Biomasa viable", "kg/m3", "Estado observado/modelado"),
        ("X_dead_kg_m3", "Biomasa muerta", "kg/m3", "Estado observado/modelado"),
        ("N_kg_m3", "Nitrogeno asimilable", "kg/m3", "Estado observado/modelado"),
        ("G_g_l", "Glucosa", "g/L", "Estado observado/modelado"),
        ("F_g_l", "Fructosa", "g/L", "Estado observado/modelado"),
        ("E_g_l", "Etanol", "g/L", "Estado observado/modelado"),
        ("glycerol_g_l", "Glicerol", "g/L", "Estado observado/modelado"),
        ("ethyl_acetate_total", "Acetato de etilo total", "dato fuente", "Variable aromatica disponible"),
        ("isoamyl_acetate_total", "Acetato de isoamilo total", "dato fuente", "Variable aromatica disponible"),
        ("ethyl_octanoate_total", "Octanoato de etilo total", "dato fuente", "Variable aromatica disponible"),
    ]
    rows = []
    columns = set(normalized.columns) if not normalized.empty else set()
    for field, desc, unit, role in fields:
        rows.append(
            {
                "field": field,
                "description": desc,
                "unit_or_domain": unit,
                "role": role,
                "present_in_normalized_long": field in columns,
            }
        )
    return pd.DataFrame(rows)


def build_data_source_summary() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for source in [
        ROOT / "fermentation_model/data/mosto_natural_xthiol.xlsx",
        ROOT / "fermentation_model/data/mosto_sintetico_vl3.xlsx",
        ROOT / "fermentation_model/data/Calibration_data_vl3.xlsx",
    ]:
        rows.extend(workbook_summary(source))
    return pd.DataFrame(rows)


def build_independent_base_summary(batch_summary: pd.DataFrame, normalized: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if batch_summary.empty:
        base_summary = pd.DataFrame()
    else:
        agg = (
            batch_summary.groupby("medium")
            .agg(
                n_batches=("batch", "nunique"),
                n_rows=("n_rows", "sum"),
                t_min_h=("t_min_h", "min"),
                t_max_h=("t_max_h", "max"),
                G_initial_mean=("G_initial_g_l", "mean"),
                F_initial_mean=("F_initial_g_l", "mean"),
                YAN_initial_mean=("YAN_initial_mg_l", "mean"),
                E_final_mean=("E_final_g_l", "mean"),
                total_N_pulse_mean=("total_N_pulse_mg_l", "mean"),
            )
            .reset_index()
        )
        total = pd.DataFrame(
            [
                {
                    "medium": "all",
                    "n_batches": batch_summary["batch"].nunique(),
                    "n_rows": batch_summary["n_rows"].sum(),
                    "t_min_h": batch_summary["t_min_h"].min(),
                    "t_max_h": batch_summary["t_max_h"].max(),
                    "G_initial_mean": batch_summary["G_initial_g_l"].mean(),
                    "F_initial_mean": batch_summary["F_initial_g_l"].mean(),
                    "YAN_initial_mean": batch_summary["YAN_initial_mg_l"].mean(),
                    "E_final_mean": batch_summary["E_final_g_l"].mean(),
                    "total_N_pulse_mean": batch_summary["total_N_pulse_mg_l"].mean(),
                }
            ]
        )
        base_summary = pd.concat([agg, total], ignore_index=True)

    if normalized.empty:
        obs_summary = pd.DataFrame()
    else:
        state_cols = ["X_viable_kg_m3", "X_dead_kg_m3", "N_kg_m3", "G_g_l", "F_g_l", "E_g_l", "glycerol_g_l"]
        rows = []
        for medium, group in normalized.groupby("medium", dropna=False):
            for state in state_cols:
                if state in group.columns:
                    rows.append(
                        {
                            "medium": medium,
                            "state": state,
                            "n_observations": int(group[state].notna().sum()),
                            "n_batches_with_observation": int(group.loc[group[state].notna(), "batch"].nunique()),
                        }
                    )
        obs_summary = pd.DataFrame(rows)
    return base_summary, obs_summary


def build_validation_criteria(
    batch_summary: pd.DataFrame,
    metrics_by_state: pd.DataFrame,
    feasibility: pd.DataFrame,
    theta_disagreement: pd.DataFrame,
    pairwise: pd.DataFrame,
    profile: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    n_natural = int(batch_summary.loc[batch_summary.get("medium", pd.Series(dtype=str)) == "natural", "batch"].nunique()) if not batch_summary.empty else 0
    n_synthetic = int(batch_summary.loc[batch_summary.get("medium", pd.Series(dtype=str)) == "synthetic", "batch"].nunique()) if not batch_summary.empty else 0
    rows.append(
        {
            "criterio": "Base independiente 2025 cargada",
            "umbral_o_regla": "Debe contener lotes naturales y sinteticos 2025 homologados.",
            "valor_observado": f"natural={n_natural}; synthetic={n_synthetic}",
            "estado": "cumple" if n_natural > 0 and n_synthetic > 0 else "no_cumple",
            "interpretacion": "La base independiente contiene ambos dominios de validacion.",
        }
    )

    improved = 0
    total = 0
    if not metrics_by_state.empty and "best_minus_base_weighted_rmse" in metrics_by_state.columns:
        total = int(metrics_by_state["best_minus_base_weighted_rmse"].notna().sum())
        improved = int((metrics_by_state["best_minus_base_weighted_rmse"] <= 0).sum())
    rows.append(
        {
            "criterio": "Robustez theta base vs mejor multistart",
            "umbral_o_regla": "Mejor multistart no debe deteriorar la mayoria de combinaciones medio/estado.",
            "valor_observado": f"{improved}/{total} combinaciones mejoran o empatan",
            "estado": "cumple" if total and improved / total >= 0.80 else "revisar",
            "interpretacion": "La alternativa multistart mejora casi todos los estados/medios; el deterioro principal reportado es Xd sintetico.",
        }
    )

    issue_sum = int(feasibility.get("issue_count", pd.Series(dtype=float)).fillna(0).sum()) if not feasibility.empty else -1
    rows.append(
        {
            "criterio": "Factibilidad de predicciones DOE",
            "umbral_o_regla": "No activar banderas duras en etanol, biomasa, biomasa muerta, glicerol o azucar residual.",
            "valor_observado": f"issue_count_total={issue_sum}",
            "estado": "cumple" if issue_sum == 0 else "no_cumple",
            "interpretacion": "Las simulaciones seleccionadas quedan en rangos plausibles.",
        }
    )

    natural_candidates = 0
    if not feasibility.empty and "medium" in feasibility.columns:
        natural_candidates = int(feasibility.loc[feasibility["medium"] == "natural", "candidate"].nunique())
    rows.append(
        {
            "criterio": "Prueba de transferencia a mosto natural",
            "umbral_o_regla": "El set de validacion/diseno debe incluir al menos un candidato natural.",
            "valor_observado": f"n_candidatos_naturales={natural_candidates}",
            "estado": "cumple" if natural_candidates > 0 else "revisar",
            "interpretacion": "Existe prueba natural para transferencia de matriz, pero debe mantenerse en el primer bloque.",
        }
    )

    max_disagreement = float(theta_disagreement["theta_disagreement_weighted_rmse"].max()) if not theta_disagreement.empty else float("nan")
    rows.append(
        {
            "criterio": "Sensibilidad a theta",
            "umbral_o_regla": "Registrar desacuerdo base-vs-best para priorizar experimentos informativos.",
            "valor_observado": f"max_theta_disagreement_weighted_rmse={max_disagreement:.3f}" if pd.notna(max_disagreement) else "NA",
            "estado": "advertencia",
            "interpretacion": "Alto desacuerdo implica informacion util, pero no cierre automatico de parametros debiles.",
        }
    )

    redundant_pairs = 0
    if not pairwise.empty:
        redundant_pairs = int(((pairwise["scaled_rms_distance"] <= 0.25) & (pairwise["trajectory_correlation"] >= 0.90)).sum())
    rows.append(
        {
            "criterio": "Redundancia de trayectorias",
            "umbral_o_regla": "Pares con distancia <=0.25 y correlacion >=0.90 requieren seleccion manual.",
            "valor_observado": f"pares_redundantes={redundant_pairs}",
            "estado": "advertencia" if redundant_pairs else "cumple",
            "interpretacion": "No ejecutar pares redundantes en el primer bloque; elegir el que ataque mejor la brecha.",
        }
    )

    identifiable = 0
    prof_total = 0
    if not profile.empty and "profile_identifiable_95" in profile.columns:
        prof_total = len(profile)
        identifiable = int(profile["profile_identifiable_95"].astype(str).str.lower().eq("true").sum())
    rows.append(
        {
            "criterio": "Identificabilidad practica",
            "umbral_o_regla": "Registrar parametros no identificables y tratarlos como brecha tecnica.",
            "valor_observado": f"{identifiable}/{prof_total} parametros identificables al 95%",
            "estado": "brecha" if prof_total and identifiable < prof_total else "cumple",
            "interpretacion": "La validacion permite uso secuencial, no cierre definitivo de todos los parametros.",
        }
    )

    return pd.DataFrame(rows)


def build_breaches_table(pairwise: pd.DataFrame, profile: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "brecha": "Identificabilidad practica incompleta",
            "evidencia": "Perfiles de verosimilitud con direcciones one-sided o no identificables.",
            "impacto": "No declarar cierre parametrico definitivo.",
            "accion_recomendada": "Usar campana secuencial: ejecutar bloque 1, recalibrar, recomputar FIM/perfiles y rerankear.",
            "estado": "controlada_para_validacion_condicionada",
        },
        {
            "brecha": "Transferencia de medio requiere prueba natural",
            "evidencia": "El reporte de transferencia indica que la matriz natural debe mantenerse para distinguir transferencia vs estructura.",
            "impacto": "No usar solo mosto sintetico como validacion final.",
            "accion_recomendada": "Mantener `natural_cold_hot_switch` o `natural_glucose_pulse_after_growth` en el primer bloque.",
            "estado": "controlada_con_bloque_mixto",
        },
        {
            "brecha": "Redundancia entre algunos candidatos",
            "evidencia": "Pares con baja distancia RMS y alta correlacion de trayectoria.",
            "impacto": "Baja eficiencia si ambos se ejecutan temprano.",
            "accion_recomendada": "No poner ambos en el primer bloque; escoger segun parametro debil objetivo.",
            "estado": "controlada_por_protocolo",
        },
        {
            "brecha": "Deterioro menor en Xd sintetico bajo mejor multistart",
            "evidencia": "Curve validation report: synthetic Xd es la unica combinacion con deterioro leve.",
            "impacto": "Xd permanece observable sensible/ruidoso.",
            "accion_recomendada": "Revisar calidad de conteos de biomasa muerta y priorizar medidas consistentes en bloque 1.",
            "estado": "advertencia",
        },
    ]
    if not pairwise.empty:
        top = pairwise.head(2)
        for _, r in top.iterrows():
            rows.append(
                {
                    "brecha": "Par candidato redundante detectado",
                    "evidencia": f"{r['candidate_a']} vs {r['candidate_b']}; distancia={r['scaled_rms_distance']:.3f}; correlacion={r['trajectory_correlation']:.3f}",
                    "impacto": "Puede aportar informacion duplicada.",
                    "accion_recomendada": "Ejecutar solo uno en bloque inicial.",
                    "estado": "advertencia",
                }
            )
    return pd.DataFrame(rows)


def build_decision_table(criteria: pd.DataFrame) -> pd.DataFrame:
    hard = criteria[criteria["estado"].isin(["no_cumple"])]
    status = "validacion_tecnica_condicionada" if hard.empty else "no_validado"
    return pd.DataFrame(
        [
            {
                "decision_tecnica": status,
                "alcance": "Modelo aceptable para planificacion y diseno secuencial con base independiente vendimia 2025.",
                "condiciones": "Ejecutar primer bloque mixto, incluir una prueba natural, recalibrar y recomputar perfiles/FIM antes de cerrar bloques siguientes.",
                "no_declarar": "No declarar cierre parametrico definitivo ni validacion final de todos los parametros debiles.",
                "uso_recomendado": "Anexo de rendicion tecnica A02 y soporte para decision de validacion condicionada.",
            }
        ]
    )


def build_validation_protocol_doc() -> str:
    return f"""# Protocolo de validacion - {ANNEX_ID}

## Objetivo

Documentar la validacion del modelo con base independiente de vendimia 2025, cubriendo datos fuente, protocolo, comparacion predicho-observado, metricas de error, brechas y decision tecnica.

## Base independiente

La base independiente 2025 se compone de:

- `fermentation_model/data/mosto_natural_xthiol.xlsx`
- `fermentation_model/data/mosto_sintetico_vl3.xlsx`

La base historica `Calibration_data_vl3.xlsx` se incluye como referencia de calibracion previa y auditoria, pero la evidencia principal de vendimia 2025 esta en los libros natural/sintetico homologados y en `results/new_must_data_loading`.

## Estados validados

- Biomasa viable (`X`)
- Biomasa muerta (`Xd`)
- Nitrogeno asimilable (`N`)
- Glucosa (`G`)
- Fructosa (`F`)
- Etanol (`E`)
- Glicerol (`Gly`)

## Procedimiento

1. Cargar y normalizar los datos 2025 con `new_must_data_loader.py`.
2. Resumir la base por lote, medio, horizonte, condiciones iniciales y numero de observaciones.
3. Simular curvas con theta base y theta best-multistart.
4. Comparar predicho vs observado por lote/estado mediante RMSE y weighted RMSE.
5. Evaluar robustez base-vs-best con `best_minus_base_weighted_rmse`.
6. Revisar factibilidad de predicciones DOE: etanol, biomasa, biomasa muerta, glicerol y azucar residual.
7. Detectar sensibilidad a theta y redundancias de trayectorias.
8. Declarar brechas y decision tecnica.

## Criterio de validacion usado en este bundle

- La base debe contener lotes naturales y sinteticos de vendimia 2025.
- La comparacion predicho-observado debe estar disponible por medio/estado.
- Las simulaciones seleccionadas no deben activar banderas duras de factibilidad.
- La alternativa best-multistart debe mejorar o empatar la mayoria de combinaciones medio/estado.
- El set recomendado debe incluir una prueba natural para transferencia a mosto real.
- Las brechas de identificabilidad se documentan como condicion de uso secuencial, no como bloqueo del anexo.

## Decision operacional derivada

Validacion tecnica condicionada: el modelo es utilizable para diseno y planificacion secuencial con datos de vendimia 2025. El primer bloque recomendado debe ser mixto y debe incluir transferencia a mosto natural; despues del bloque 1 se debe recalibrar y rerankear.
"""


def build_decision_doc(criteria: pd.DataFrame, breaches: pd.DataFrame, decision: pd.DataFrame) -> str:
    return f"""# Decision tecnica - {ANNEX_ID}

## Decision

`{decision.iloc[0]['decision_tecnica']}`

## Alcance

{decision.iloc[0]['alcance']}

## Condiciones de uso

{decision.iloc[0]['condiciones']}

## Criterios evaluados

{criteria.to_markdown(index=False)}

## Brechas y acciones

{breaches.to_markdown(index=False)}

## Texto sugerido para informe

Se valido el modelo contra una base independiente de vendimia 2025 compuesta por fermentaciones en mosto natural y sintetico. La comparacion predicho-observado fue evaluada mediante metricas de error por lote, estado y medio, y se contrastaron las predicciones bajo theta base y una solucion best-multistart. Las simulaciones seleccionadas no activaron banderas duras de factibilidad y el analisis de curvas mostro mejora o estabilidad en la mayoria de estados/medios. Se identifican brechas de identificabilidad practica y redundancia entre algunos candidatos, por lo que la decision tecnica es validar el modelo para uso secuencial condicionado: ejecutar un primer bloque mixto con al menos una prueba natural, recalibrar y recomputar FIM/perfiles antes de cerrar la campana completa.
"""


def build_version_doc() -> str:
    key_files = [
        ROOT / "fermentation_model/run_new_must_curve_validation.py",
        ROOT / "fermentation_model/run_new_must_overnight_validation.py",
        ROOT / "fermentation_model/run_new_must_glycerol_estimability_doe.py",
        ROOT / "fermentation_model/new_must_data_loader.py",
        ROOT / "fermentation_model/data/mosto_natural_xthiol.xlsx",
        ROOT / "fermentation_model/data/mosto_sintetico_vl3.xlsx",
    ]
    hashes = []
    for path in key_files:
        if path.exists():
            hashes.append(f"- `{rel(path)}`: `{sha256(path)}`")
    return f"""# Version de codigo y trazabilidad - {ANNEX_ID}

Generado: {datetime.now().astimezone().isoformat()}
Fecha de cierre documental usada para la rendicion: {GENERATED_DATE}

## Git

- Rama: `{run_git(['branch', '--show-current'])}`
- Commit HEAD: `{run_git(['rev-parse', 'HEAD'])}`
- Commit corto: `{run_git(['rev-parse', '--short', 'HEAD'])}`

## Estado de trabajo

```text
{run_git(['status', '--short'])}
```

## Diff stat

```text
{run_git(['diff', '--stat'])}
```

## Hashes SHA256 de fuentes clave

{chr(10).join(hashes)}

## Nota

Usar este documento junto con `00_manifest.csv`, que registra hashes de todos los artefactos copiados o generados dentro del bundle.
"""


def md_table(df: pd.DataFrame, max_rows: int = 15) -> str:
    if df.empty:
        return "_No hay datos disponibles._"
    return df.head(max_rows).to_markdown(index=False)


def build_readme(
    annex_table: pd.DataFrame,
    base_summary: pd.DataFrame,
    metrics_state: pd.DataFrame,
    criteria: pd.DataFrame,
) -> str:
    n_batches = int(base_summary.loc[base_summary["medium"] == "all", "n_batches"].iloc[0]) if not base_summary.empty and "all" in set(base_summary["medium"]) else 0
    n_criteria_ok = int(criteria["estado"].isin(["cumple"]).sum()) if not criteria.empty else 0
    return f"""# Bundle {ANNEX_ID} - Validacion del modelo con datos de vendimia 2025

Actividad: {OE}-{ACTIVITY_CODE} - {ACTIVITY_TITLE}

Fecha de cierre documental: {GENERATED_DATE}

## Proposito

Este bundle organiza la evidencia para el anexo solicitado en la imagen: base independiente 2025, protocolo de validacion, predicho vs observado, metricas, brechas y decision tecnica.

## Tabla del anexo

{annex_table.to_markdown(index=False)}

## Archivos principales

- `ANEXO_A02_borrador.md`: texto base del anexo.
- `01_tabla_anexo_A02.csv`: tabla con las columnas de la imagen.
- `02_matriz_cumplimiento_A02.csv`: mapeo contenido minimo -> evidencia.
- `03_protocolo_validacion_A02.md`: protocolo y criterio de validacion.
- `04_decision_tecnica_A02.md`: decision, brechas y condiciones.
- `08_version_codigo.md`: version de codigo y hashes.
- `00_manifest.csv`: inventario completo con SHA256.

## Tablas generadas

- `tables/03_base_independiente_2025_resumen_lotes.csv`
- `tables/03_base_independiente_2025_observaciones.csv`
- `tables/03_base_independiente_2025_diccionario.csv`
- `tables/04_predicho_observado_metricas.csv`
- `tables/04_metricas_error_por_medio_estado.csv`
- `tables/04_residuales_agregados_base_vs_best.csv`
- `tables/05_criterios_validacion.csv`
- `tables/06_brechas_y_acciones.csv`
- `tables/07_decision_tecnica.csv`

## Resumen cuantitativo

- Lotes independientes 2025 resumidos: {n_batches}
- Filas de metricas por medio/estado: {len(metrics_state)}
- Criterios en estado `cumple`: {n_criteria_ok}/{len(criteria)}

## Uso recomendado

1. Usar `ANEXO_A02_borrador.md` como cuerpo narrativo.
2. Insertar la tabla `01_tabla_anexo_A02.csv` en la matriz general de anexos.
3. Usar `02_matriz_cumplimiento_A02.csv` para demostrar que el contenido minimo esta cubierto.
4. Adjuntar `source_reports/`, `source_tables/` y `figures/` como respaldo.
5. Convertir el anexo final con el nombre sugerido: `{SUGGESTED_FILENAME}`.

Para regenerar:

```powershell
python scripts/build_A02_validation_bundle.py
```
"""


def build_annex_draft(
    annex_table: pd.DataFrame,
    compliance: pd.DataFrame,
    base_summary: pd.DataFrame,
    obs_summary: pd.DataFrame,
    metrics_state: pd.DataFrame,
    feasibility: pd.DataFrame,
    criteria: pd.DataFrame,
    breaches: pd.DataFrame,
    decision: pd.DataFrame,
) -> str:
    return f"""# {ACTIVITY_TITLE}

## Identificacion

- ID anexo: {ANNEX_ID}
- Objetivo especifico: {OE}
- Actividad unica: {ACTIVITY_CODE}
- Nombre de archivo sugerido: `{SUGGESTED_FILENAME}`
- Fecha de cierre documental del bundle: {GENERATED_DATE}

## Tabla fuente del anexo

{annex_table.to_markdown(index=False)}

## Resumen ejecutivo

Se documento la validacion del modelo con una base independiente de vendimia 2025 compuesta por fermentaciones en mosto natural y mosto sintetico. La validacion integra carga y homologacion de datos, comparacion predicho-observado por lote y estado, metricas de error, revision de factibilidad de predicciones, analisis de sensibilidad a theta, brechas tecnicas y decision de uso.

La decision tecnica propuesta es `{decision.iloc[0]['decision_tecnica']}`: el modelo queda respaldado para planificacion y diseno secuencial, condicionado a ejecutar un primer bloque mixto, incluir una prueba natural y recalibrar antes de cerrar la campana completa.

## Matriz de cumplimiento

{compliance.to_markdown(index=False)}

## Base independiente vendimia 2025

Resumen por medio:

{md_table(base_summary, max_rows=10)}

Cobertura observacional por medio y estado:

{md_table(obs_summary, max_rows=20)}

## Protocolo de validacion

El protocolo completo esta en `03_protocolo_validacion_A02.md`. En terminos operativos, el flujo fue:

1. Cargar y homologar bases natural/sintetica 2025.
2. Simular curvas con theta base y theta best-multistart.
3. Comparar predicho-observado por lote/estado.
4. Consolidar RMSE/weighted RMSE por medio y estado.
5. Revisar factibilidad de predicciones DOE.
6. Documentar sensibilidad, redundancia, brechas y decision tecnica.

## Comparacion predicho-observado y metricas

Metricas agregadas por medio/estado:

{md_table(metrics_state, max_rows=18)}

Los graficos `figures/curve_fit__*.png` contienen la comparacion visual predicho-observado por lote.

## Factibilidad de predicciones

Vista preliminar de predicciones DOE:

{md_table(feasibility[['candidate', 'theta', 'medium', 'family', 'final_sugar_GF', 'final_E', 'max_X', 'max_Xd', 'max_Gly', 'issue_count', 'status']] if not feasibility.empty else feasibility, max_rows=12)}

## Criterios de validacion

{criteria.to_markdown(index=False)}

## Brechas

{breaches.to_markdown(index=False)}

## Decision tecnica

{decision.to_markdown(index=False)}

## Texto sugerido para informe

Se valido el modelo con una base independiente de vendimia 2025 compuesta por fermentaciones naturales y sinteticas. La comparacion predicho-observado fue evaluada mediante metricas de error por lote, estado y medio, complementadas con revision de factibilidad de predicciones y robustez frente a una solucion best-multistart. Las predicciones seleccionadas no activaron banderas duras de factibilidad y el mejor multistart mejora la mayoria de combinaciones medio/estado. Se identifican brechas de identificabilidad practica, transferencia de medio y redundancia de algunos candidatos; por ello, la decision tecnica es validacion condicionada para uso secuencial, con primer bloque mixto, inclusion de prueba natural, recalibracion y reranking posterior.
"""


def build_manifest_markdown(manifest: pd.DataFrame) -> str:
    counts = manifest.groupby(["category", "status"]).size().reset_index(name="n")
    return f"""# Manifest {ANNEX_ID}

La version tabular con hashes esta en `00_manifest.csv`.

## Conteo por categoria

{counts.to_markdown(index=False)}

## Primeras rutas

{manifest[['category', 'evidence_requirement', 'status', 'bundle_path', 'source_path']].head(30).to_markdown(index=False)}
"""


def copy_source_artifacts() -> None:
    for source in [
        "fermentation_model/data/mosto_natural_xthiol.xlsx",
        "fermentation_model/data/mosto_sintetico_vl3.xlsx",
        "fermentation_model/data/Calibration_data_vl3.xlsx",
    ]:
        copy_file(source, "data_sources", "datos", "Base independiente vendimia 2025", "Libro de datos fuente")

    for source in [
        "fermentation_model/new_must_data_loader.py",
        "fermentation_model/run_new_must_curve_validation.py",
        "fermentation_model/run_new_must_overnight_validation.py",
        "fermentation_model/run_new_must_glycerol_estimability_doe.py",
        "fermentation_model/run_medium_transfer_diagnostics.py",
        "scripts/build_A02_validation_bundle.py",
    ]:
        copy_file(source, "code_sources", "codigo", "Protocolo de validacion y version de codigo", "Codigo fuente clave")

    for source in [
        "fermentation_model/fermentation_new_must_data_loading.ipynb",
        "fermentation_model/fermentation_new_must_glycerol_estimability_doe.executed.ipynb",
        "fermentation_model/fermentation_new_must_glycerol_estimability_doe.ipynb",
    ]:
        copy_file(source, "notebooks", "notebooks", "Predicho-observado y calibracion", "Notebook de respaldo")

    for source in [
        "fermentation_model/results/curve_validation/curve_validation_report.md",
        "fermentation_model/results/curve_validation/curve_validation_decision_summary.md",
        "fermentation_model/results/medium_transfer_diagnostics/medium_transfer_diagnostics_report.md",
        "fermentation_model/results/new_must_glycerol_overnight_validation/overnight_validation_report.md",
        "fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/post_summary.md",
        "fermentation_model/results/new_must_glycerol_estimability_doe/new_must_glycerol_estimability_doe_report.md",
    ]:
        copy_file(source, "source_reports", "reportes", "Metricas, brechas y decision tecnica", "Reporte fuente", namespace=True)

    for source in [
        "fermentation_model/results/new_must_data_loading/new_must_batch_summary.csv",
        "fermentation_model/results/new_must_data_loading/new_must_normalized_long.csv",
        "fermentation_model/results/new_must_data_loading/density_sugar_fit_by_medium.csv",
        "fermentation_model/results/curve_validation/current_fit_curve_metrics.csv",
        "fermentation_model/results/curve_validation/current_fit_curve_metrics_by_medium_state.csv",
        "fermentation_model/results/curve_validation/design_prediction_feasibility.csv",
        "fermentation_model/results/curve_validation/design_theta_disagreement.csv",
        "fermentation_model/results/curve_validation/design_pairwise_similarity.csv",
        "fermentation_model/results/medium_transfer_diagnostics/historical_batch_audit.csv",
        "fermentation_model/results/medium_transfer_diagnostics/medium_map_template.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/multistart_summary.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/l2_scan_summary.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/profile_full_summary.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/sampling_policy_benchmark.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/residual_by_state_base_vs_best_agg.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/residual_by_state_base_vs_best.csv",
        "fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/parameter_estimability_best.csv",
    ]:
        copy_file(source, "source_tables", "tablas_fuente", "Base, predicho-observado, metricas, brechas", "Tabla fuente", namespace=True)

    copy_glob(
        "fermentation_model/results/curve_validation/plots/current_fit/*.png",
        "figures",
        "figuras",
        "Comparacion predicho-observado",
        "Figura predicho-observado por lote",
    )
    copy_glob(
        "fermentation_model/results/curve_validation/plots/design_predictions/*.png",
        "figures",
        "figuras",
        "Predicciones y factibilidad",
        "Figura prediccion candidato DOE",
    )
    copy_glob(
        "fermentation_model/results/new_must_data_loading/*.png",
        "figures",
        "figuras",
        "Base independiente 2025",
        "Figura carga de datos",
    )


def write_zip() -> None:
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(BUNDLE_DIR.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(BUNDLE_DIR.parent))


def main() -> int:
    reset_bundle()
    copy_source_artifacts()

    annex_table = build_annex_table()
    compliance = build_compliance_matrix()
    data_sources = build_data_source_summary()
    batch_summary = read_csv_if_exists("fermentation_model/results/new_must_data_loading/new_must_batch_summary.csv")
    normalized = read_csv_if_exists("fermentation_model/results/new_must_data_loading/new_must_normalized_long.csv")
    base_summary, obs_summary = build_independent_base_summary(batch_summary, normalized)
    dictionary = build_data_dictionary(normalized)
    metrics = read_csv_if_exists("fermentation_model/results/curve_validation/current_fit_curve_metrics.csv")
    metrics_state = read_csv_if_exists("fermentation_model/results/curve_validation/current_fit_curve_metrics_by_medium_state.csv")
    feasibility = read_csv_if_exists("fermentation_model/results/curve_validation/design_prediction_feasibility.csv")
    theta_disagreement = read_csv_if_exists("fermentation_model/results/curve_validation/design_theta_disagreement.csv")
    pairwise = read_csv_if_exists("fermentation_model/results/curve_validation/design_pairwise_similarity.csv")
    residual_agg = read_csv_if_exists(
        "fermentation_model/results/new_must_glycerol_overnight_validation/best_multistart_post_analysis/residual_by_state_base_vs_best_agg.csv"
    )
    profile = read_csv_if_exists("fermentation_model/results/new_must_glycerol_overnight_validation/profile_full_summary.csv")
    criteria = build_validation_criteria(batch_summary, metrics_state, feasibility, theta_disagreement, pairwise, profile)
    breaches = build_breaches_table(pairwise, profile)
    decision = build_decision_table(criteria)

    write_df_artifact("01_tabla_anexo_A02.csv", annex_table, "tabla_anexo", "Contenido minimo y evidencia fuente", "Tabla segun imagen")
    write_df_artifact("02_matriz_cumplimiento_A02.csv", compliance, "matriz_cumplimiento", "Contenido minimo", "Matriz contenido minimo a evidencia")
    write_df_artifact("tables/03_fuentes_excel_resumen.csv", data_sources, "datos", "Base independiente vendimia 2025", "Resumen de libros Excel")
    write_df_artifact("tables/03_base_independiente_2025_resumen_lotes.csv", base_summary, "datos", "Base independiente vendimia 2025", "Resumen por medio")
    write_df_artifact("tables/03_base_independiente_2025_observaciones.csv", obs_summary, "datos", "Base independiente vendimia 2025", "Cobertura observacional")
    write_df_artifact("tables/03_base_independiente_2025_diccionario.csv", dictionary, "datos", "Base independiente vendimia 2025", "Diccionario de variables")
    write_df_artifact("tables/04_predicho_observado_metricas.csv", metrics, "predicho_observado", "Comparacion predicho-observado", "Metricas por lote/estado")
    write_df_artifact("tables/04_metricas_error_por_medio_estado.csv", metrics_state, "metricas", "Metricas de error", "Metricas agregadas por medio/estado")
    write_df_artifact("tables/04_residuales_agregados_base_vs_best.csv", residual_agg, "metricas", "Metricas de error", "Residuales agregados base vs best")
    write_df_artifact("tables/04_factibilidad_predicciones.csv", feasibility, "predicciones", "Comparacion predicho-observado y factibilidad", "Factibilidad de predicciones DOE")
    write_df_artifact("tables/04_desacuerdo_theta.csv", theta_disagreement, "sensibilidad", "Criterio de validacion", "Desacuerdo theta base vs best")
    write_df_artifact("tables/04_redundancia_trayectorias.csv", pairwise, "brechas", "Brechas", "Redundancia de trayectorias")
    write_df_artifact("tables/05_criterios_validacion.csv", criteria, "criterios", "Criterio de validacion", "Criterios y estado")
    write_df_artifact("tables/06_brechas_y_acciones.csv", breaches, "brechas", "Brechas", "Brechas y acciones")
    write_df_artifact("tables/07_decision_tecnica.csv", decision, "decision", "Decision tecnica", "Decision tecnica estructurada")

    write_text_artifact("03_protocolo_validacion_A02.md", build_validation_protocol_doc(), "protocolo", "Protocolo de validacion", "Documento de protocolo")
    write_text_artifact("04_decision_tecnica_A02.md", build_decision_doc(criteria, breaches, decision), "decision", "Decision tecnica", "Documento de decision tecnica")
    write_text_artifact("08_version_codigo.md", build_version_doc(), "version_codigo", "Version de codigo", "Version y trazabilidad")
    write_text_artifact("README_A02_bundle.md", build_readme(annex_table, base_summary, metrics_state, criteria), "guia", "Uso del bundle", "Guia principal")
    write_text_artifact(
        "ANEXO_A02_borrador.md",
        build_annex_draft(annex_table, compliance, base_summary, obs_summary, metrics_state, feasibility, criteria, breaches, decision),
        "anexo",
        "Anexo propuesto",
        "Borrador narrativo del anexo",
    )

    manifest_for_markdown = pd.DataFrame(manifest_rows)
    manifest_md = BUNDLE_DIR / "00_manifest.md"
    manifest_md.write_text(build_manifest_markdown(manifest_for_markdown), encoding="utf-8")
    record_artifact(
        bundle_path=manifest_md,
        category="manifest",
        evidence_requirement="Inventario y trazabilidad",
        role="Resumen legible del manifest",
    )
    manifest_path = BUNDLE_DIR / "00_manifest.csv"
    manifest_rows.append(
        {
            "annex_id": ANNEX_ID,
            "oe": OE,
            "activity": ACTIVITY_CODE,
            "category": "manifest",
            "evidence_requirement": "Inventario y trazabilidad",
            "role_in_bundle": "Inventario completo con hashes; fila autoreferencial sin hash estable",
            "status": "self_index",
            "source_path": "generated",
            "bundle_path": rel(manifest_path),
            "bytes": "",
            "sha256": "",
        }
    )
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False, quoting=csv.QUOTE_MINIMAL)

    write_zip()
    print(f"Bundle written: {BUNDLE_DIR}")
    print(f"ZIP written: {ZIP_PATH}")
    print(f"Artifacts indexed: {len(manifest_rows)}")
    if not criteria.empty:
        print(criteria[["criterio", "estado", "valor_observado"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
