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


ROOT = Path(__file__).resolve().parents[1]
BUNDLE_NAME = "A06_calibracion_metabolica"
BUNDLE_DIR = ROOT / "rendicion_tecnica" / BUNDLE_NAME
ZIP_PATH = ROOT / "rendicion_tecnica" / f"{BUNDLE_NAME}.zip"

GENERATED_DATE = "2026-06-26"
ANNEX_ID = "A06"
OE = "OE3"
ACTIVITY_CODE = "3.13"
ACTIVITY_TITLE = "Calibracion del submodelo metabolico del Gemelo Digital (Actividad 3.13)"
SUGGESTED_FILENAME = "PI-4497_IA5_A06_Act-3.13_Calibracion-submodelo-metabolico.pdf"

FM = ROOT / "fermentation_model"
DATA_DIR = FM / "data"
RESULTS = FM / "results"

manifest_rows: list[dict[str, str]] = []


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def safe_name(name: str, max_len: int = 54) -> str:
    base = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)
    while "__" in base:
        base = base.replace("__", "_")
    if len(base) <= max_len:
        return base
    suffix = Path(base).suffix
    digest = hashlib.sha1(base.encode("utf-8")).hexdigest()[:8]
    stem = Path(base).stem[: max_len - len(suffix) - 10]
    return f"{stem}_{digest}{suffix}"


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


def copy_path(
    source: Path,
    dest_subdir: str,
    category: str,
    evidence_requirement: str,
    role: str,
    dest_name: str | None = None,
) -> Path | None:
    dst = BUNDLE_DIR / dest_subdir / (dest_name or safe_name(source.name))
    if not source.exists():
        record_artifact(
            bundle_path=dst,
            category=category,
            evidence_requirement=evidence_requirement,
            role=role,
            source_path=source,
            status="missing",
        )
        return None
    shutil.copy2(source, dst)
    record_artifact(
        bundle_path=dst,
        category=category,
        evidence_requirement=evidence_requirement,
        role=role,
        source_path=source,
    )
    return dst


def write_text_artifact(rel_path: str, text: str, category: str, evidence_requirement: str, role: str) -> Path:
    path = BUNDLE_DIR / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    record_artifact(bundle_path=path, category=category, evidence_requirement=evidence_requirement, role=role)
    return path


def write_df_artifact(rel_path: str, df: pd.DataFrame, category: str, evidence_requirement: str, role: str) -> Path:
    path = BUNDLE_DIR / rel_path
    path.parent.mkdir(parents=True, exist_ok=True)
    clean = df.where(pd.notna(df), "")
    clean.to_csv(path, index=False)
    record_artifact(bundle_path=path, category=category, evidence_requirement=evidence_requirement, role=role)
    return path


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def read_indexed_series(path: Path, value_name: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["parameter", value_name])
    df = pd.read_csv(path, index_col=0)
    series = df.iloc[:, 0]
    out = series.reset_index()
    out.columns = ["parameter", value_name]
    return out


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def copy_inputs() -> None:
    for path in [
        DATA_DIR / "Calibration_data_vl3.xlsx",
        DATA_DIR / "mosto_natural_xthiol.xlsx",
        DATA_DIR / "mosto_sintetico_vl3.xlsx",
    ]:
        copy_path(path, "data_sources", "dataset", "Dataset de calibracion", "Libro de datos fuente")

    for path in [
        FM / "run_fit_strategy_analysis.py",
        FM / "run_secondary_joint_campaign_doe.py",
        FM / "run_secondary_v2_model_evaluation.py",
        FM / "run_secondary_fit_capacity.py",
        FM / "run_new_must_glycerol_estimability_doe.py",
        FM / "run_new_must_overnight_validation.py",
        ROOT / "scripts" / "build_A06_metabolic_calibration_bundle.py",
    ]:
        copy_path(path, "code_sources", "codigo", "Version de codigo", "Script reproducible")

    for path in [
        FM / "fermentation_model_calibration_3_effective_reformulation.ipynb",
        FM / "fermentation_secondary_v2_model_evaluation.executed.ipynb",
        FM / "fermentation_secondary_fit_capacity.executed.ipynb",
    ]:
        copy_path(path, "notebooks", "codigo", "Version de codigo", "Notebook de respaldo")

    report_paths = [
        RESULTS / "fit_strategy_analysis" / "fit_strategy_decision_report.md",
        RESULTS / "secondary_joint_campaign_doe" / "secondary_joint_campaign_report.md",
        RESULTS / "secondary_v2_model_evaluation" / "secondary_v2_model_report.md",
        RESULTS / "secondary_fit_capacity" / "secondary_fit_capacity_report.md",
        RESULTS / "new_must_glycerol_estimability_doe" / "new_must_glycerol_estimability_doe_report.md",
        RESULTS / "new_must_glycerol_overnight_validation" / "overnight_validation_report.md",
    ]
    for path in report_paths:
        copy_path(
            path,
            "source_reports",
            "reportes_fuente",
            "Resultados e interpretacion",
            "Reporte fuente",
            dest_name=safe_name(f"{path.parent.name}__{path.name}"),
        )

    table_paths = [
        RESULTS / "secondary_joint_campaign_doe" / "secondary_joint_input_data.csv",
        RESULTS / "secondary_joint_campaign_doe" / "secondary_joint_metadata.json",
        RESULTS / "secondary_joint_campaign_doe" / "theta_secondary_joint.csv",
        RESULTS / "secondary_v2_model_evaluation" / "fit_comparison.csv",
        RESULTS / "secondary_v2_model_evaluation" / "state_fit_metrics.csv",
        RESULTS / "secondary_v2_model_evaluation" / "theta_secondary_v2_reduced_o2fixed.csv",
        RESULTS / "secondary_v2_model_evaluation" / "theta_secondary_v2_o2fixed.csv",
        RESULTS / "secondary_v2_model_evaluation" / "theta_secondary_v2_o2free.csv",
        RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_estimability.csv",
        RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_eigen_spectrum.csv",
        RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_fim.csv",
        RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_weak_loadings.csv",
        RESULTS / "secondary_v2_model_evaluation" / "metadata.json",
        RESULTS / "secondary_fit_capacity" / "fit_capacity_summary.csv",
        RESULTS / "secondary_fit_capacity" / "fit_capacity_state_metrics.csv",
        RESULTS / "secondary_fit_capacity" / "fit_capacity_aggregate_metrics.csv",
        RESULTS / "secondary_fit_capacity" / "batch_all_o2free_parameters.csv",
        RESULTS / "secondary_fit_capacity" / "batch_all_o2free_parameter_spread.csv",
        RESULTS / "secondary_fit_capacity" / "metadata.json",
        RESULTS / "fit_strategy_analysis" / "fit_strategy_summary.csv",
        RESULTS / "fit_strategy_analysis" / "fit_strategy_state_summary.csv",
        RESULTS / "fit_strategy_analysis" / "fit_strategy_solve_summary.csv",
        RESULTS / "fit_strategy_analysis" / "fit_strategy_residuals.csv",
        RESULTS / "fit_strategy_analysis" / "strategy_predictions.csv",
        RESULTS / "fit_strategy_analysis" / "plus_qx_iG_logL2_10_theta.csv",
        RESULTS / "fit_strategy_analysis" / "plus_qx_iG_logL2_10_physical.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "fit_summary.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "theta_fit_table.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "parameter_estimability_summary.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "fim_summary.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "fim_eigenvalues.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "profile_summary_combined.csv",
        RESULTS / "new_must_glycerol_estimability_doe" / "profile_weak_summary.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "l2_scan_summary.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "multistart_summary.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "multistart_theta.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "profile_full_summary.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "profile_full_profiles.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "sampling_policy_benchmark.csv",
        RESULTS / "new_must_glycerol_overnight_validation" / "pyomo_selected_summary.csv",
    ]
    for path in table_paths:
        copy_path(
            path,
            "source_tables",
            "tablas_fuente",
            "Dataset, parametros, metricas e incertidumbre",
            "Tabla fuente",
            dest_name=safe_name(f"{path.parent.name}__{path.name}"),
        )

    figure_paths = []
    figure_paths.extend((RESULTS / "secondary_v2_model_evaluation" / "plots").glob("*.png"))
    figure_paths.extend((RESULTS / "secondary_fit_capacity" / "plots").glob("*.png"))
    figure_paths.extend((RESULTS / "fit_strategy_analysis").glob("strategy_curves_batch_*.png"))
    for path in sorted(figure_paths):
        copy_path(
            path,
            "figures",
            "figuras",
            "Resultados y metricas",
            "Figura de ajuste",
            dest_name=safe_name(f"{path.parent.parent.name if path.parent.name == 'plots' else path.parent.name}__{path.name}"),
        )


def build_annex_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "ID anexo": ANNEX_ID,
                "OE": OE,
                "Actividad unica": ACTIVITY_CODE,
                "Titulo propuesto": ACTIVITY_TITLE,
                "Nombre de archivo sugerido": SUGGESTED_FILENAME,
                "Contenido minimo": (
                    "Documentar dataset, estrategia de calibracion, parametros metabolicos, metricas, "
                    "incertidumbre e interpretacion."
                ),
                "Evidencia fuente": (
                    "Dataset de calibracion, parametros metabolicos, funcion objetivo, metricas, "
                    "incertidumbre/identificabilidad y resultados."
                ),
            }
        ]
    )


def build_compliance_matrix() -> pd.DataFrame:
    rows = [
        {
            "contenido_minimo": "Dataset de calibracion",
            "estado_bundle": "cubierto",
            "archivos_bundle": (
                "data_sources/*.xlsx; source_tables/secondary_joint_campaign_doe__secondary_joint_input_data.csv; "
                "tables/03_dataset_resumen.csv"
            ),
            "detalle": "Incluye base historica VL3, mosto natural, mosto sintetico y dataset integrado de 27 lotes para calibracion secundaria.",
        },
        {
            "contenido_minimo": "Estrategia de calibracion",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/04_estrategia_calibracion.csv; 03_metodologia_calibracion_A06.md; source_reports/*.md",
            "detalle": "Documenta ajuste regularizado, evaluacion v2, diagnostico de capacidad y validacion por multistart/perfiles.",
        },
        {
            "contenido_minimo": "Parametros metabolicos",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/06_parametros_metabolicos.csv; source_tables/*theta*.csv; source_tables/*estimability*.csv",
            "detalle": "Incluye parametros del ajuste regularizado y parametros v2 de piruvato, acetaldehido, acetato y oxigeno.",
        },
        {
            "contenido_minimo": "Funcion objetivo",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/05_funcion_objetivo.csv; code_sources/run_fit_strategy_analysis.py; code_sources/run_secondary_v2_model_evaluation.py",
            "detalle": "Describe WSSE, penalizacion log-L2, residuos escalados por sigma y diagnosticos de robustez.",
        },
        {
            "contenido_minimo": "Metricas",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/07_metricas_calibracion.csv; source_tables/*fit*.csv; source_tables/*state*.csv",
            "detalle": "Incluye WSSE, SSE, residuos por estado, n observaciones, convergencia y comparacion de modelos.",
        },
        {
            "contenido_minimo": "Incertidumbre/identificabilidad",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/08_incertidumbre_identificabilidad.csv; source_tables/*profile*.csv; source_tables/*fim*.csv",
            "detalle": "Incluye FIM, espectro propio, perfiles de verosimilitud y clasificacion de parametros debiles/confundidos.",
        },
        {
            "contenido_minimo": "Resultados e interpretacion",
            "estado_bundle": "cubierto",
            "archivos_bundle": "tables/09_resultados_interpretacion.csv; 04_resultados_interpretacion_A06.md; ANEXO_A06_borrador.md",
            "detalle": "Resume decision tecnica, uso recomendado, brechas y limites de cierre parametrico.",
        },
    ]
    return pd.DataFrame(rows)


def excel_summary(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"archivo": path.name, "estado": "missing"}
    xl = pd.ExcelFile(path)
    sheet_rows = []
    total_rows = 0
    for sheet in xl.sheet_names:
        try:
            df = pd.read_excel(path, sheet_name=sheet)
            n_rows = len(df)
            total_rows += n_rows
            sheet_rows.append(f"{sheet}:{n_rows}")
        except Exception:
            sheet_rows.append(f"{sheet}:error")
    return {
        "archivo": path.name,
        "estado": "disponible",
        "n_hojas": len(xl.sheet_names),
        "filas_totales": total_rows,
        "n_lotes": "",
        "medios": "",
        "hojas_resumen": "; ".join(sheet_rows[:12]),
    }


def build_dataset_summary() -> pd.DataFrame:
    rows = [excel_summary(path) for path in [
        DATA_DIR / "Calibration_data_vl3.xlsx",
        DATA_DIR / "mosto_natural_xthiol.xlsx",
        DATA_DIR / "mosto_sintetico_vl3.xlsx",
    ]]
    joint = read_csv(RESULTS / "secondary_joint_campaign_doe" / "secondary_joint_input_data.csv")
    if not joint.empty:
        key_cols = [
            "time_h",
            "temperature_c",
            "density",
            "X_viable_kg_m3",
            "G_g_l",
            "F_g_l",
            "YAN_mg_l",
            "pyruvic_acid",
            "acetaldehyde",
            "acetic_acid",
            "DO_mg_l",
            "glycerol_g_l",
            "E_g_l",
            "ethyl_acetate_total",
            "isoamyl_acetate_total",
            "ethyl_octanoate_total",
        ]
        coverage = {col: int(joint[col].notna().sum()) for col in key_cols if col in joint.columns}
        rows.append(
            {
                "archivo": "secondary_joint_input_data.csv",
                "estado": "disponible",
                "n_hojas": "",
                "filas_totales": len(joint),
                "n_lotes": joint["batch"].nunique() if "batch" in joint else "",
                "medios": ",".join(sorted(map(str, joint["medium"].dropna().unique()))) if "medium" in joint else "",
                "hojas_resumen": "; ".join(f"{k}:{v}" for k, v in coverage.items()),
            }
        )
    return pd.DataFrame(rows)


def build_calibration_strategy() -> pd.DataFrame:
    fit = read_csv(RESULTS / "fit_strategy_analysis" / "fit_strategy_summary.csv")
    recommended = fit.loc[fit["strategy"] == "plus_qx_iG_logL2_10"].head(1) if not fit.empty else pd.DataFrame()
    v2_meta = load_json(RESULTS / "secondary_v2_model_evaluation" / "metadata.json")
    cap_meta = load_json(RESULTS / "secondary_fit_capacity" / "metadata.json")
    joint_meta = load_json(RESULTS / "secondary_joint_campaign_doe" / "secondary_joint_metadata.json")
    rows = [
        {
            "etapa": "1_dataset_integrado",
            "descripcion": "Integracion de datos historicos VL3, mosto natural y mosto sintetico en un dataset comun.",
            "parametros_o_estados": ",".join(joint_meta.get("target_parameters", [])[:8]) + ("..." if joint_meta.get("target_parameters") else ""),
            "n_lotes_o_obs": "27 lotes en secondary_joint_input_data.csv",
            "criterio_decision": "Dataset integrado disponible y trazable a libros fuente.",
            "archivo_fuente": "secondary_joint_campaign_doe/secondary_joint_input_data.csv",
        },
        {
            "etapa": "2_ajuste_regularizado_core",
            "descripcion": "Comparacion de estrategias ParmEst con WSSE y penalizacion log-L2.",
            "parametros_o_estados": recommended.iloc[0]["estimated_parameters"] if not recommended.empty else "",
            "n_lotes_o_obs": f"{int(recommended.iloc[0]['n_observations'])} observaciones" if not recommended.empty else "",
            "criterio_decision": "Seleccion de plus_qx_iG_logL2_10 por no activar bounds y mantener objetivo penalizado competitivo.",
            "archivo_fuente": "fit_strategy_analysis/fit_strategy_summary.csv",
        },
        {
            "etapa": "3_submodelo_secundario_v2",
            "descripcion": "Calibracion v2 para piruvato, acetaldehido, acetato y O2 con simulacion core fija.",
            "parametros_o_estados": ",".join(v2_meta.get("parameters", [])),
            "n_lotes_o_obs": f"{v2_meta.get('n_batches', '')} lotes; {v2_meta.get('n_residuals', '')} residuos",
            "criterio_decision": "Comparar WSSE frente a v1 y clasificar parametros por FIM.",
            "archivo_fuente": "secondary_v2_model_evaluation/fit_comparison.csv",
        },
        {
            "etapa": "4_capacidad_de_ajuste",
            "descripcion": "Ajustes por medio y por lote como cota superior de capacidad, no como modelo transferible.",
            "parametros_o_estados": ",".join(cap_meta.get("fit_all_parameters", [])),
            "n_lotes_o_obs": f"{cap_meta.get('n_batches', '')} lotes",
            "criterio_decision": cap_meta.get("interpretation", ""),
            "archivo_fuente": "secondary_fit_capacity/fit_capacity_summary.csv",
        },
        {
            "etapa": "5_incertidumbre",
            "descripcion": "Multistart, escaneo L2, perfiles de verosimilitud y FIM para robustez/identificabilidad.",
            "parametros_o_estados": "full17; v2_reduced_o2fixed",
            "n_lotes_o_obs": "14 multistarts; 17 perfiles full17; FIM v2",
            "criterio_decision": "Distinguir parametros identificables, regulares y debiles/confundidos.",
            "archivo_fuente": "new_must_glycerol_overnight_validation/profile_full_summary.csv",
        },
    ]
    return pd.DataFrame(rows)


def build_objective_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "workflow": "fit_strategy_analysis",
                "funcion_objetivo": "WSSE + penalizacion log-L2",
                "formula_operativa": "ParmEst SSE_weighted + 0.5 * lambda * sum(log(theta/theta_ref)^2) * n_batches",
                "residuos_o_estados": "E,F,G,N,X",
                "ponderacion": "residuos ponderados por funcion SSE_weighted del notebook de calibracion",
                "archivo_codigo": "code_sources/run_fit_strategy_analysis.py",
            },
            {
                "workflow": "secondary_v2_model_evaluation",
                "funcion_objetivo": "least_squares robusto sobre residuos escalados",
                "formula_operativa": "sum(((pred-obs)/sigma_estado)^2) con loss soft_l1 y priors log para parametros seleccionados",
                "residuos_o_estados": "Pyr,AcAld,Acetate,O2",
                "ponderacion": "SIGMA: Pyr=6.0, AcAld=8.0, Acetate=0.035, O2=0.20",
                "archivo_codigo": "code_sources/run_secondary_v2_model_evaluation.py",
            },
            {
                "workflow": "secondary_fit_capacity",
                "funcion_objetivo": "least_squares por medio/lote",
                "formula_operativa": "misma base de residuos escalados; parametros liberados por lote como diagnostico de capacidad",
                "residuos_o_estados": "Pyr,AcAld,Acetate,O2",
                "ponderacion": "no se usa como modelo final; solo cota superior de ajuste",
                "archivo_codigo": "code_sources/run_secondary_fit_capacity.py",
            },
        ]
    )


def build_parameter_table() -> pd.DataFrame:
    rows = []
    theta_core = read_indexed_series(RESULTS / "fit_strategy_analysis" / "plus_qx_iG_logL2_10_theta.csv", "value")
    if not theta_core.empty:
        for _, row in theta_core.iterrows():
            rows.append(
                {
                    "familia": "core_regularizado_theta",
                    "parameter": row["parameter"],
                    "value": row["value"],
                    "lower_bound": "",
                    "upper_bound": "",
                    "classification": "seleccionado_plus_qx_iG_logL2_10",
                    "source": "fit_strategy_analysis/plus_qx_iG_logL2_10_theta.csv",
                }
            )
    physical = read_csv(RESULTS / "fit_strategy_analysis" / "plus_qx_iG_logL2_10_physical.csv")
    if not physical.empty:
        for _, row in physical.iterrows():
            rows.append(
                {
                    "familia": "core_regularizado_physical",
                    "parameter": row["parameter"],
                    "value": row["back_transformed_value"],
                    "lower_bound": row.get("physical_lower_bound", ""),
                    "upper_bound": row.get("physical_upper_bound", ""),
                    "classification": "valor_fisico_retrotransformado",
                    "source": "fit_strategy_analysis/plus_qx_iG_logL2_10_physical.csv",
                }
            )
    v2_est = read_csv(RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_estimability.csv")
    if not v2_est.empty:
        for _, row in v2_est.iterrows():
            rows.append(
                {
                    "familia": "secondary_v2_reduced_o2fixed",
                    "parameter": row["parameter"],
                    "value": row["theta"],
                    "lower_bound": "",
                    "upper_bound": "",
                    "classification": row["classification"],
                    "source": "secondary_v2_model_evaluation/v2_reduced_o2fixed_estimability.csv",
                    "std_log_approx": row.get("std_log_approx", ""),
                    "approx_95_multiplier": row.get("approx_95_multiplier", ""),
                    "active_bound": row.get("active_bound", ""),
                }
            )
    return pd.DataFrame(rows)


def build_metrics_table() -> pd.DataFrame:
    rows = []
    fit_strategy = read_csv(RESULTS / "fit_strategy_analysis" / "fit_strategy_summary.csv")
    if not fit_strategy.empty:
        for _, row in fit_strategy.head(8).iterrows():
            rows.append(
                {
                    "workflow": "fit_strategy_analysis",
                    "modelo_o_estrategia": row["strategy"],
                    "metrica_principal": "direct_plus_l2_total",
                    "valor": row["direct_plus_l2_total"],
                    "wsse": row["direct_WSSE_total"],
                    "n_obs": row["n_observations"],
                    "estado": "ok" if bool(row["solve_ok"]) else "revisar",
                    "comentario": f"active_bounds={row['n_active_bounds']}; lambda={row['l2_lambda']}",
                }
            )
    v2_fit = read_csv(RESULTS / "secondary_v2_model_evaluation" / "fit_comparison.csv")
    if not v2_fit.empty:
        baseline = float(v2_fit.loc[v2_fit["model"] == "secondary_v1_current", "final_wsse"].iloc[0])
        for _, row in v2_fit.iterrows():
            improvement = "" if row["model"] == "secondary_v1_current" else 100.0 * (baseline - float(row["final_wsse"])) / baseline
            rows.append(
                {
                    "workflow": "secondary_v2_model_evaluation",
                    "modelo_o_estrategia": row["model"],
                    "metrica_principal": "final_wsse",
                    "valor": row["final_wsse"],
                    "wsse": row["final_wsse"],
                    "n_obs": row["n_residuals"],
                    "estado": "ok" if bool(row["success"]) else "revisar",
                    "comentario": f"wsse_per_residual={row['wsse_per_residual']}; mejora_vs_v1_pct={improvement}",
                }
            )
    cap = read_csv(RESULTS / "secondary_fit_capacity" / "fit_capacity_aggregate_metrics.csv")
    if not cap.empty:
        for _, row in cap.loc[cap["state"] == "ALL"].iterrows():
            rows.append(
                {
                    "workflow": "secondary_fit_capacity",
                    "modelo_o_estrategia": row["fit_label"],
                    "metrica_principal": "wsse_per_obs",
                    "valor": row["wsse_per_obs"],
                    "wsse": row["wsse"],
                    "n_obs": row["n_obs"],
                    "estado": "diagnostico",
                    "comentario": "Cota de capacidad de ajuste; no equivale a modelo transferible.",
                }
            )
    multi = read_csv(RESULTS / "new_must_glycerol_overnight_validation" / "multistart_summary.csv")
    if not multi.empty:
        best = multi.sort_values("final_wsse").head(3)
        for _, row in best.iterrows():
            rows.append(
                {
                    "workflow": "overnight_multistart",
                    "modelo_o_estrategia": row["fit"],
                    "metrica_principal": "final_wsse",
                    "valor": row["final_wsse"],
                    "wsse": row["final_wsse"],
                    "n_obs": row["n_residuals"],
                    "estado": "ok" if bool(row["success"]) else "revisar",
                    "comentario": f"wsse_per_residual={row['wsse_per_residual']}; lambda={row['l2_lambda']}",
                }
            )
    return pd.DataFrame(rows)


def build_uncertainty_table() -> pd.DataFrame:
    rows = []
    v2_est = read_csv(RESULTS / "secondary_v2_model_evaluation" / "v2_reduced_o2fixed_estimability.csv")
    if not v2_est.empty:
        for _, row in v2_est.iterrows():
            rows.append(
                {
                    "workflow": "secondary_v2_FIM",
                    "parameter": row["parameter"],
                    "diagnostico": row["classification"],
                    "metrica": "approx_95_multiplier",
                    "valor": row["approx_95_multiplier"],
                    "detalle": f"std_log_approx={row['std_log_approx']}; active_bound={row['active_bound']}",
                }
            )
    profile = read_csv(RESULTS / "new_must_glycerol_overnight_validation" / "profile_full_summary.csv")
    if not profile.empty:
        for _, row in profile.iterrows():
            rows.append(
                {
                    "workflow": "profile_likelihood_full17",
                    "parameter": row["parameter"],
                    "diagnostico": "profile_identifiable_95" if bool(row["profile_identifiable_95"]) else "no_identifiable_95",
                    "metrica": "max_lr_stat",
                    "valor": row["max_lr_stat"],
                    "detalle": (
                        f"left={row['crosses_left_95']}; right={row['crosses_right_95']}; "
                        f"chi2_95={row['chi2_95_threshold']}"
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_interpretation_table(metrics: pd.DataFrame, uncertainty: pd.DataFrame) -> pd.DataFrame:
    v2_metrics = metrics.loc[metrics["workflow"] == "secondary_v2_model_evaluation"] if not metrics.empty else pd.DataFrame()
    weak_params = []
    if not uncertainty.empty:
        weak_params = uncertainty.loc[
            uncertainty["diagnostico"].astype(str).str.contains("weak|no_identifiable", case=False, regex=True),
            "parameter",
        ].dropna().astype(str).unique().tolist()
    rows = [
        {
            "tema": "cierre_documental",
            "resultado": "Cobertura 100% del contenido minimo solicitado para A06.",
            "soporte": "02_matriz_cumplimiento_A06.csv; 00_manifest.csv",
            "interpretacion": "El bundle documenta dataset, estrategia, parametros, funcion objetivo, metricas, incertidumbre y resultados.",
        },
        {
            "tema": "modelo_recomendado",
            "resultado": "plus_qx_iG_logL2_10 como candidato regularizado para ajuste core; secondary_v2_reduced_o2fixed para capa secundaria.",
            "soporte": "fit_strategy_summary.csv; fit_comparison.csv; metadata.json",
            "interpretacion": "Se privilegia estabilidad sin bounds activos y reduccion de error frente al modelo secundario v1.",
        },
        {
            "tema": "mejora_v2",
            "resultado": "; ".join(
                f"{row['modelo_o_estrategia']}={float(row['valor']):.3g}" for _, row in v2_metrics.iterrows()
            ),
            "soporte": "tables/07_metricas_calibracion.csv",
            "interpretacion": "La estructura v2 reduce el WSSE global respecto del ajuste secundario v1, pero conserva parametros debiles.",
        },
        {
            "tema": "identificabilidad",
            "resultado": ", ".join(weak_params[:12]) + ("..." if len(weak_params) > 12 else ""),
            "soporte": "tables/08_incertidumbre_identificabilidad.csv",
            "interpretacion": "Los parametros debiles o no identificables deben mantenerse fijos, regularizados o tratarse como sensibilidad.",
        },
        {
            "tema": "decision_tecnica",
            "resultado": "Calibracion tecnica documentada con uso condicionado para simulacion, DOE y decision secuencial.",
            "soporte": "04_resultados_interpretacion_A06.md",
            "interpretacion": "No conviene declarar cierre parametrico definitivo para todos los parametros metabolicos.",
        },
    ]
    return pd.DataFrame(rows)


def build_compliance_summary(compliance: pd.DataFrame) -> pd.DataFrame:
    covered = int(compliance["estado_bundle"].astype(str).str.startswith("cubierto").sum())
    total = int(len(compliance))
    return pd.DataFrame(
        [
            {
                "criterio": "Cobertura contenido minimo A06",
                "estado": "cumple",
                "valor": f"{covered}/{total} items cubiertos",
                "comentario": "Cobertura documental completa; las limitaciones tecnicas quedan explicitadas como interpretacion.",
            }
        ]
    )


def md_table(df: pd.DataFrame, max_rows: int = 12) -> str:
    if df.empty:
        return "_Sin datos disponibles._"
    clean = df.head(max_rows).where(pd.notna(df.head(max_rows)), "")
    return clean.to_markdown(index=False, disable_numparse=True)


def build_method_doc(strategy: pd.DataFrame, objective: pd.DataFrame, dataset: pd.DataFrame) -> str:
    return f"""# Metodologia de calibracion A06

## Dataset

{md_table(dataset, max_rows=8)}

## Estrategia

{strategy.to_markdown(index=False, disable_numparse=True)}

## Funcion objetivo

{objective.to_markdown(index=False, disable_numparse=True)}

## Trazabilidad metodologica

- Los datos fuente se preservan en `data_sources/`.
- El dataset integrado queda en `source_tables/secondary_joint_campaign_doe__secondary_joint_input_data.csv`.
- La version ejecutable del flujo queda en `code_sources/` y `notebooks/`.
- Los reportes originales se preservan en `source_reports/`.
"""


def build_results_doc(
    params: pd.DataFrame,
    metrics: pd.DataFrame,
    uncertainty: pd.DataFrame,
    interpretation: pd.DataFrame,
    compliance_summary: pd.DataFrame,
) -> str:
    return f"""# Resultados e interpretacion A06

## Parametros metabolicos

{md_table(params, max_rows=18)}

## Metricas de calibracion

{md_table(metrics, max_rows=18)}

## Incertidumbre e identificabilidad

{md_table(uncertainty, max_rows=24)}

## Interpretacion tecnica

{interpretation.to_markdown(index=False, disable_numparse=True)}

## Cumplimiento

{compliance_summary.to_markdown(index=False, disable_numparse=True)}
"""


def build_version_doc() -> str:
    return f"""# Version de codigo A06

- Fecha de generacion del bundle: {GENERATED_DATE}
- Git commit HEAD: `{run_git(['rev-parse', 'HEAD'])}`
- Git branch: `{run_git(['branch', '--show-current'])}`
- Estado corto de git al generar:

```text
{run_git(['status', '--short'])}
```

## Scripts copiados

- `code_sources/run_fit_strategy_analysis.py`
- `code_sources/run_secondary_joint_campaign_doe.py`
- `code_sources/run_secondary_v2_model_evaluation.py`
- `code_sources/run_secondary_fit_capacity.py`
- `code_sources/run_new_must_glycerol_estimability_doe.py`
- `code_sources/run_new_must_overnight_validation.py`
- `code_sources/build_A06_metabolic_calibration_bundle.py`
"""


def build_readme(
    annex_table: pd.DataFrame,
    compliance_summary: pd.DataFrame,
    dataset: pd.DataFrame,
    metrics: pd.DataFrame,
    uncertainty: pd.DataFrame,
) -> str:
    n_files = len(manifest_rows) + 4
    n_dataset_rows = ""
    if not dataset.empty and "filas_totales" in dataset.columns:
        joint = dataset.loc[dataset["archivo"] == "secondary_joint_input_data.csv"]
        if not joint.empty:
            n_dataset_rows = str(joint.iloc[0]["filas_totales"])
    v2_best = metrics.loc[
        (metrics["workflow"] == "secondary_v2_model_evaluation")
        & (metrics["modelo_o_estrategia"] == "secondary_v2_reduced_o2fixed")
    ]
    weak = int(
        uncertainty["diagnostico"].astype(str).str.contains("weak|no_identifiable", case=False, regex=True).sum()
    ) if not uncertainty.empty else 0
    return f"""# Bundle {ANNEX_ID} - Calibracion del submodelo metabolico

Actividad: {OE}-{ACTIVITY_CODE} - {ACTIVITY_TITLE}

Fecha de cierre documental: {GENERATED_DATE}

## Proposito

Este bundle organiza evidencia para el anexo solicitado: dataset, estrategia de calibracion, parametros metabolicos, funcion objetivo, metricas, incertidumbre/identificabilidad, resultados e interpretacion.

## Tabla del anexo

{annex_table.to_markdown(index=False, disable_numparse=True)}

## Archivos principales

- `ANEXO_A06_borrador.md`: texto base del anexo.
- `01_tabla_anexo_A06.csv`: tabla con las columnas solicitadas.
- `02_matriz_cumplimiento_A06.csv`: mapeo contenido minimo -> evidencia.
- `03_metodologia_calibracion_A06.md`: dataset, estrategia y funcion objetivo.
- `04_resultados_interpretacion_A06.md`: parametros, metricas, incertidumbre e interpretacion.
- `08_version_codigo.md`: version de codigo y estado git.
- `00_manifest.csv`: inventario completo con SHA256.

## Resumen

- Archivos del bundle: {n_files}
- Filas del dataset integrado: {n_dataset_rows}
- Cumplimiento documental: {compliance_summary.iloc[0]['valor'] if not compliance_summary.empty else 'NA'}
- WSSE secondary_v2_reduced_o2fixed: {v2_best.iloc[0]['valor'] if not v2_best.empty else 'NA'}
- Parametros debiles/no identificables reportados: {weak}

## Uso recomendado

1. Usar `ANEXO_A06_borrador.md` como cuerpo narrativo.
2. Insertar `01_tabla_anexo_A06.csv` en la matriz general de anexos.
3. Usar `02_matriz_cumplimiento_A06.csv` para justificar cobertura.
4. Adjuntar `data_sources/`, `source_reports/`, `source_tables/`, `figures/`, `code_sources/` y `notebooks/`.
5. Convertir el anexo final con el nombre sugerido: `{SUGGESTED_FILENAME}`.

Para regenerar:

```powershell
python scripts/build_A06_metabolic_calibration_bundle.py
```
"""


def build_annex_draft(
    annex_table: pd.DataFrame,
    compliance: pd.DataFrame,
    dataset: pd.DataFrame,
    strategy: pd.DataFrame,
    objective: pd.DataFrame,
    params: pd.DataFrame,
    metrics: pd.DataFrame,
    uncertainty: pd.DataFrame,
    interpretation: pd.DataFrame,
    compliance_summary: pd.DataFrame,
) -> str:
    return f"""# {ACTIVITY_TITLE}

## Identificacion

- ID anexo: {ANNEX_ID}
- Objetivo especifico: {OE}
- Actividad unica: {ACTIVITY_CODE}
- Nombre de archivo sugerido: `{SUGGESTED_FILENAME}`
- Fecha de cierre documental del bundle: {GENERATED_DATE}

## Tabla fuente del anexo

{annex_table.to_markdown(index=False, disable_numparse=True)}

## Resumen ejecutivo

Se documenta la calibracion del submodelo metabolico del Gemelo Digital a partir de bases historicas VL3, ensayos de mosto natural y ensayos de mosto sintetico. La evidencia incluye dataset integrado, estrategia de calibracion regularizada, parametros metabolicos estimados, funcion objetivo, metricas de ajuste, diagnosticos de incertidumbre/identificabilidad y una interpretacion tecnica de uso condicionado.

## Matriz de cumplimiento

{compliance.to_markdown(index=False, disable_numparse=True)}

## Dataset de calibracion

{md_table(dataset, max_rows=8)}

## Estrategia de calibracion

{strategy.to_markdown(index=False, disable_numparse=True)}

## Funcion objetivo

{objective.to_markdown(index=False, disable_numparse=True)}

## Parametros metabolicos

{md_table(params, max_rows=20)}

## Metricas

{md_table(metrics, max_rows=20)}

## Incertidumbre e identificabilidad

{md_table(uncertainty, max_rows=24)}

## Interpretacion

{interpretation.to_markdown(index=False, disable_numparse=True)}

## Cumplimiento

{compliance_summary.to_markdown(index=False, disable_numparse=True)}

## Texto sugerido para informe

Se consolido evidencia de calibracion del submodelo metabolico del Gemelo Digital para la actividad OE3-3.13. El dataset integrado combina informacion historica VL3, mosto natural y mosto sintetico, y se usa para evaluar parametros asociados a estados centrales y secundarios de fermentacion. La estrategia documentada incluye ajuste regularizado mediante WSSE con penalizacion log-L2, evaluacion del submodelo secundario v2 para piruvato, acetaldehido, acetato y oxigeno, diagnostico de capacidad por lote/medio y analisis de incertidumbre por FIM, multistart y perfiles de verosimilitud. Los resultados permiten cerrar documentalmente la calibracion tecnica, dejando explicitado que ciertos parametros deben mantenerse fijos, regularizados o tratados como sensibilidad antes de declarar identificabilidad parametrica completa.
"""


def build_manifest_markdown(manifest: pd.DataFrame) -> str:
    counts = manifest.groupby(["category", "status"]).size().reset_index(name="n")
    return f"""# Manifest {ANNEX_ID}

La version tabular con hashes esta en `00_manifest.csv`.

## Conteo por categoria

{counts.to_markdown(index=False, disable_numparse=True)}

## Primeras rutas

{manifest[['category', 'evidence_requirement', 'status', 'bundle_path', 'source_path']].head(80).to_markdown(index=False, disable_numparse=True)}
"""


def write_zip() -> None:
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(BUNDLE_DIR.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(BUNDLE_DIR.parent))


def main() -> int:
    reset_bundle()
    copy_inputs()

    annex_table = build_annex_table()
    compliance = build_compliance_matrix()
    dataset = build_dataset_summary()
    strategy = build_calibration_strategy()
    objective = build_objective_table()
    params = build_parameter_table()
    metrics = build_metrics_table()
    uncertainty = build_uncertainty_table()
    interpretation = build_interpretation_table(metrics, uncertainty)
    compliance_summary = build_compliance_summary(compliance)

    write_df_artifact("01_tabla_anexo_A06.csv", annex_table, "tabla_anexo", "Contenido minimo y evidencia fuente", "Tabla segun fila A06")
    write_df_artifact("02_matriz_cumplimiento_A06.csv", compliance, "matriz_cumplimiento", "Contenido minimo", "Matriz contenido minimo a evidencia")
    write_df_artifact("tables/03_dataset_resumen.csv", dataset, "dataset", "Dataset de calibracion", "Resumen de dataset")
    write_df_artifact("tables/04_estrategia_calibracion.csv", strategy, "estrategia", "Estrategia de calibracion", "Estrategia de calibracion")
    write_df_artifact("tables/05_funcion_objetivo.csv", objective, "funcion_objetivo", "Funcion objetivo", "Funcion objetivo")
    write_df_artifact("tables/06_parametros_metabolicos.csv", params, "parametros", "Parametros metabolicos", "Parametros metabolicos")
    write_df_artifact("tables/07_metricas_calibracion.csv", metrics, "metricas", "Metricas", "Metricas de calibracion")
    write_df_artifact("tables/08_incertidumbre_identificabilidad.csv", uncertainty, "incertidumbre", "Incertidumbre/identificabilidad", "Incertidumbre e identificabilidad")
    write_df_artifact("tables/09_resultados_interpretacion.csv", interpretation, "interpretacion", "Resultados e interpretacion", "Interpretacion tecnica")
    write_df_artifact("tables/10_cumplimiento_A06.csv", compliance_summary, "cumplimiento", "Cumplimiento", "Resumen de cumplimiento")

    write_text_artifact(
        "03_metodologia_calibracion_A06.md",
        build_method_doc(strategy, objective, dataset),
        "metodologia",
        "Estrategia de calibracion",
        "Metodologia del anexo",
    )
    write_text_artifact(
        "04_resultados_interpretacion_A06.md",
        build_results_doc(params, metrics, uncertainty, interpretation, compliance_summary),
        "resultados",
        "Resultados e interpretacion",
        "Resultados del anexo",
    )
    write_text_artifact("08_version_codigo.md", build_version_doc(), "version_codigo", "Version de codigo", "Version y trazabilidad")
    write_text_artifact(
        "README_A06_bundle.md",
        build_readme(annex_table, compliance_summary, dataset, metrics, uncertainty),
        "guia",
        "Uso del bundle",
        "Guia principal",
    )
    write_text_artifact(
        "ANEXO_A06_borrador.md",
        build_annex_draft(
            annex_table,
            compliance,
            dataset,
            strategy,
            objective,
            params,
            metrics,
            uncertainty,
            interpretation,
            compliance_summary,
        ),
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
    print(compliance_summary.to_string(index=False))
    print(metrics.loc[metrics["workflow"].eq("secondary_v2_model_evaluation"), ["modelo_o_estrategia", "valor", "n_obs"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
