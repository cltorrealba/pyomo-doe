from __future__ import annotations

import hashlib
import math
import re
import shutil
import sys
import textwrap
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = SCRIPT_DIR / "results"
BUNDLE_NAME = "A04_vino_condensados_aroma_bundle_20260701"
BUNDLE_DIR = RESULTS_DIR / BUNDLE_NAME
FIG_DIR = BUNDLE_DIR / "figures"
TABLE_DIR = BUNDLE_DIR / "tables"
TRACE_DIR = BUNDLE_DIR / "source_trace"
TEMPLATE_DIR = BUNDLE_DIR / "templates"
SNAPSHOT_DIR = BUNDLE_DIR / "source_snapshot"

LAB_SYNTHETIC = SCRIPT_DIR / "data" / "Laboratorio 2025" / "mosto_sintetico_vl3.xlsx"
PILOT_WORKBOOK = SCRIPT_DIR / "data" / "Piloto 2025" / "Calibration_data_vl3.xlsx"
A07_BUNDLE = SCRIPT_DIR / "pilot_2025" / "bundles" / "A07_aroma_nlp_calibration"
A15_BUNDLE = SCRIPT_DIR / "pilot_2025" / "bundles" / "A15_pilot_scalability_2025"
AROMA_DOE = RESULTS_DIR / "aroma_joint_campaign_doe"
OPERATIONAL_DOE = RESULTS_DIR / "final_operational_doe_volume_constrained"

TARGET_AROMAS = [
    "ethyl_acetate",
    "isoamyl_acetate",
    "ethyl_octanoate",
]

AROMA_TOTAL_COLUMNS = [
    "bencil_alcohol_total",
    "benzaldehido_total",
    "decanoato_de_etilo_total",
    "hexil_acetate_total",
    "isoamil_acetate_total",
    "octanoate_de_etilo_total",
    "phenylethylacetate_total",
    "Ethyl_Acetate_total",
]

AROMA_CONDENSATE_COLUMNS = [
    "bencil_alcohol_condensado",
    "benzaldehido_condensado",
    "decanoato_de_etilo_condensado",
    "hexil_acetate_condensado",
    "isoamil_acetate_condensado",
    "octanoate_de_etilo_condensado",
    "phenylethylacetate_condensado",
    "Ethyl_Acetate_condensado",
]

PROCESS_COLUMNS = [
    "densidad",
    "Viability",
    "Oculyze Concentration",
    "GLUCOSE",
    "FRUCTOSE",
    "PAN",
    "AMMONIA",
    "YAN",
    "GLYCEROL",
    "PYRUVIC ACID",
    "ACETALDEHIDO",
    "ETANOL",
    "pulso_nut",
]

SOURCE_KEYWORDS = [
    "muestra",
    "muestras",
    "aroma",
    "aromas",
    "condens",
    "gc",
    "cromat",
    "almacen",
    "volatil",
    "headspace",
    "sPME".lower(),
]

HIGH_VALUE_SOURCE_KEYWORDS = [
    "muestra",
    "muestras",
    "muestreo",
    "aroma",
    "aromas",
    "condens",
    "gc",
    "cromat",
    "almacen",
    "volatil",
    "headspace",
    "spme",
    "protocolo",
    "analisis",
    "análisis",
]


def resolve_one_drive_root() -> Path | None:
    for parent in SCRIPT_DIR.resolve().parents:
        if parent.name.startswith("OneDrive") and "Concha" in parent.name and "Toro" in parent.name:
            return parent
    candidates = list(Path.home().glob("OneDrive - Vi*a Concha y Toro S.A"))
    if candidates:
        exact = [path for path in candidates if "ñ" in path.name]
        return exact[0] if exact else candidates[0]
    for candidate in Path.home().glob("OneDrive*"):
        if "Concha" in candidate.name and "Toro" in candidate.name:
            return candidate
    return None


def discover_project_paths() -> dict[str, Path | None]:
    one_drive = resolve_one_drive_root()
    paths: dict[str, Path | None] = {
        "one_drive": one_drive,
        "project_a04_dir": None,
        "reference_pdf": None,
        "gc_protocol_docx": None,
        "pilot_2025_docx": None,
    }
    if one_drive is None:
        return paths

    project_base = one_drive / "Documentos" / "Proyectos I+D" / "PI-4497" / "Ley I+D" / "Rendición Técnica 6"
    if not project_base.exists():
        ley_dir = one_drive / "Documentos" / "Proyectos I+D" / "PI-4497" / "Ley I+D"
        matches = list(ley_dir.glob("Rendici*cnica 6")) if ley_dir.exists() else []
        project_base = matches[0] if matches else project_base

    anexos = project_base / "Anexos"
    paths["project_a04_dir"] = anexos / "A04 - Actividad 3.11"
    if project_base.exists():
        pdf_matches = list(project_base.rglob("CII_VCT_Trabajo_de_Titulo_final.pdf"))
        docx_2026_matches = list(project_base.rglob("Ensayos I+D Piloto 2026*.docx"))
        docx_2025_matches = list(project_base.rglob("Ensayos I+D Piloto 2025*.docx"))
        paths["reference_pdf"] = pdf_matches[0] if pdf_matches else None
        paths["gc_protocol_docx"] = docx_2026_matches[0] if docx_2026_matches else None
        paths["pilot_2025_docx"] = docx_2025_matches[0] if docx_2025_matches else None
    return paths


def prepare_dirs() -> None:
    if BUNDLE_DIR.exists():
        resolved_bundle = BUNDLE_DIR.resolve()
        resolved_results = RESULTS_DIR.resolve()
        if resolved_results not in resolved_bundle.parents:
            raise RuntimeError(f"Refusing to refresh unexpected bundle path: {resolved_bundle}")
        shutil.rmtree(BUNDLE_DIR)
    for directory in (BUNDLE_DIR, FIG_DIR, TABLE_DIR, TRACE_DIR, TEMPLATE_DIR, SNAPSHOT_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def fmt(value: object, ndigits: int = 3) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, str):
        return value
    value = float(value)
    if abs(value) >= 10000 or (abs(value) > 0 and abs(value) < 0.001):
        return f"{value:.{ndigits}e}"
    return f"{value:.{ndigits}f}"


def md_table(frame: pd.DataFrame, max_rows: int | None = None, columns: list[str] | None = None) -> str:
    if frame.empty:
        return "_No data available._"
    out = frame.copy()
    if columns is not None:
        out = out[[col for col in columns if col in out.columns]]
    if max_rows is not None:
        out = out.head(max_rows)
    for id_col in ["batch", "source_id", "file_name", "source", "location", "option", "dataset_or_campaign"]:
        if id_col in out.columns:
            out[id_col] = out[id_col].astype(str)
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(fmt)
    return out.to_markdown(index=False, disable_numparse=True)


def to_numeric(series: pd.Series) -> pd.Series:
    if series.empty:
        return pd.Series(dtype=float)
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    return pd.to_numeric(
        series.astype(str)
        .str.replace(",", ".", regex=False)
        .str.replace("%", "", regex=False)
        .str.strip(),
        errors="coerce",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_manifest(paths: dict[str, Path | None]) -> pd.DataFrame:
    rows = []
    for key, role in [
        ("reference_pdf", "Camilo thesis / synthetic campaign source"),
        ("gc_protocol_docx", "GC and 2026 pilot analytical protocol"),
        ("pilot_2025_docx", "2025 pilot protocol draft, complementary"),
    ]:
        path = paths.get(key)
        if path is None or not path.exists():
            rows.append(
                {
                    "source_id": key,
                    "role": role,
                    "file_name": "not_found",
                    "path": "not_found",
                    "size_mb": np.nan,
                    "modified": "not_found",
                    "sha256": "not_found",
                }
            )
            continue
        rows.append(
            {
                "source_id": key,
                "role": role,
                "file_name": path.name,
                "path": str(path),
                "size_mb": path.stat().st_size / (1024 * 1024),
                "modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                "sha256": sha256_file(path),
            }
        )
    return pd.DataFrame(rows)


def normalize_text(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def keyword_snippet(text: str, width: int = 360, keywords: list[str] | None = None) -> str | None:
    keywords = keywords or SOURCE_KEYWORDS
    clean = normalize_text(text)
    lower = clean.lower()
    positions = [lower.find(keyword) for keyword in keywords if lower.find(keyword) >= 0]
    if not positions:
        return None
    pos = min(positions)
    start = max(0, pos - width // 2)
    end = min(len(clean), pos + width // 2)
    snippet = clean[start:end]
    if start > 0:
        snippet = "..." + snippet
    if end < len(clean):
        snippet += "..."
    return snippet


def snippet_score(text: str) -> int:
    lower = text.lower()
    score = 0
    weights = {
        "muestra": 5,
        "muestreo": 5,
        "gc": 5,
        "cromat": 5,
        "condens": 5,
        "almacen": 4,
        "spme": 4,
        "headspace": 4,
        "protocolo": 3,
        "aroma": 2,
        "volatil": 2,
        "analisis": 2,
        "análisis": 2,
    }
    for key, weight in weights.items():
        if key in lower:
            score += weight
    return score


def extract_pdf_snippets(pdf_path: Path | None, max_rows: int = 14) -> pd.DataFrame:
    if pdf_path is None or not pdf_path.exists():
        return pd.DataFrame(
            [{"source": "reference_pdf", "location": "not_found", "keyword_context": "PDF source not found"}]
        )
    try:
        from pypdf import PdfReader
    except Exception as exc:  # pragma: no cover - environment dependent
        return pd.DataFrame(
            [{"source": pdf_path.name, "location": "import_error", "keyword_context": f"pypdf not available: {exc}"}]
        )
    rows = []
    try:
        reader = PdfReader(str(pdf_path))
        for idx, page in enumerate(reader.pages, start=1):
            text = page.extract_text() or ""
            snippet = keyword_snippet(text, keywords=HIGH_VALUE_SOURCE_KEYWORDS)
            if snippet:
                rows.append(
                    {
                        "source": pdf_path.name,
                        "location": f"page {idx}",
                        "relevance_score": snippet_score(snippet),
                        "keyword_context": snippet,
                    }
                )
    except Exception as exc:
        rows.append({"source": pdf_path.name, "location": "extract_error", "relevance_score": 0, "keyword_context": str(exc)})
    out = pd.DataFrame(rows)
    if not out.empty and "relevance_score" in out:
        out = out.sort_values(["relevance_score", "location"], ascending=[False, True]).head(max_rows)
    return out


def docx_text_blocks(docx_path: Path) -> list[tuple[str, str]]:
    import xml.etree.ElementTree as ET

    blocks: list[tuple[str, str]] = []
    with zipfile.ZipFile(docx_path) as zf:
        xml_data = zf.read("word/document.xml")
    root = ET.fromstring(xml_data)
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    idx = 0
    for paragraph in root.iterfind(".//w:p", ns):
        texts = [node.text for node in paragraph.iterfind(".//w:t", ns) if node.text]
        text = normalize_text("".join(texts))
        if text:
            idx += 1
            blocks.append((f"paragraph {idx}", text))
    return blocks


def extract_docx_snippets(docx_path: Path | None, max_rows: int = 18) -> pd.DataFrame:
    if docx_path is None or not docx_path.exists():
        return pd.DataFrame(
            [{"source": "gc_protocol_docx", "location": "not_found", "keyword_context": "DOCX source not found"}]
        )
    rows = []
    try:
        for location, text in docx_text_blocks(docx_path):
            snippet = keyword_snippet(text, keywords=HIGH_VALUE_SOURCE_KEYWORDS)
            if snippet:
                rows.append(
                    {
                        "source": docx_path.name,
                        "location": location,
                        "relevance_score": snippet_score(snippet),
                        "keyword_context": snippet,
                    }
                )
    except Exception as exc:
        rows.append({"source": docx_path.name, "location": "extract_error", "relevance_score": 0, "keyword_context": str(exc)})
    out = pd.DataFrame(rows)
    if not out.empty and "relevance_score" in out:
        out = out.sort_values(["relevance_score", "location"], ascending=[False, True]).head(max_rows)
    return out


def relevant_columns(columns: list[object], candidates: list[str]) -> list[str]:
    column_map = {str(col).strip().lower(): str(col) for col in columns}
    selected = []
    for candidate in candidates:
        key = candidate.lower()
        if key in column_map:
            selected.append(column_map[key])
    return selected


def workbook_campaign_summary(workbook: Path, label: str, sheet_pattern: str) -> pd.DataFrame:
    if not workbook.exists():
        return pd.DataFrame()
    xls = pd.ExcelFile(workbook)
    rows = []
    pattern = re.compile(sheet_pattern)
    for sheet in xls.sheet_names:
        if not pattern.match(str(sheet)):
            continue
        df = pd.read_excel(workbook, sheet_name=sheet)
        df.columns = [str(col).strip() for col in df.columns]
        total_cols = relevant_columns(list(df.columns), AROMA_TOTAL_COLUMNS)
        cond_cols = relevant_columns(list(df.columns), AROMA_CONDENSATE_COLUMNS)
        process_cols = relevant_columns(list(df.columns), PROCESS_COLUMNS)
        time_col = "time" if "time" in df.columns else "t" if "t" in df.columns else None
        temp_col = "temperature" if "temperature" in df.columns else "temperatura" if "temperatura" in df.columns else None
        yan_col = "YAN" if "YAN" in df.columns else None
        glucose_col = "GLUCOSE" if "GLUCOSE" in df.columns else None
        fructose_col = "FRUCTOSE" if "FRUCTOSE" in df.columns else None
        ethanol_col = "ETANOL" if "ETANOL" in df.columns else None
        pulse_col = "pulso_nut" if "pulso_nut" in df.columns else None

        time = to_numeric(df[time_col]) if time_col else pd.Series(dtype=float)
        temperature = to_numeric(df[temp_col]) if temp_col else pd.Series(dtype=float)
        glucose = to_numeric(df[glucose_col]) if glucose_col else pd.Series(dtype=float)
        fructose = to_numeric(df[fructose_col]) if fructose_col else pd.Series(dtype=float)
        ethanol = to_numeric(df[ethanol_col]) if ethanol_col else pd.Series(dtype=float)
        yan = to_numeric(df[yan_col]) if yan_col else pd.Series(dtype=float)
        pulse = to_numeric(df[pulse_col]) if pulse_col else pd.Series(dtype=float)

        n_aroma_total = int(df[total_cols].notna().sum().sum()) if total_cols else 0
        n_aroma_cond = int(df[cond_cols].notna().sum().sum()) if cond_cols else 0
        rows.append(
            {
                "campaign": label,
                "batch": str(sheet),
                "n_rows": int(len(df)),
                "t_min_h": float(time.min()) if len(time.dropna()) else np.nan,
                "t_max_h": float(time.max()) if len(time.dropna()) else np.nan,
                "temperature_min_c": float(temperature.min()) if len(temperature.dropna()) else np.nan,
                "temperature_max_c": float(temperature.max()) if len(temperature.dropna()) else np.nan,
                "S_initial_g_l": float((glucose + fructose).dropna().iloc[0]) if len((glucose + fructose).dropna()) else np.nan,
                "S_final_g_l": float((glucose + fructose).dropna().iloc[-1]) if len((glucose + fructose).dropna()) else np.nan,
                "YAN_initial_mg_l": float(yan.dropna().iloc[0]) if len(yan.dropna()) else np.nan,
                "E_final_raw": float(ethanol.dropna().iloc[-1]) if len(ethanol.dropna()) else np.nan,
                "n_process_obs": int(df[process_cols].notna().sum().sum()) if process_cols else 0,
                "n_aroma_total_obs": n_aroma_total,
                "n_aroma_condensate_obs": n_aroma_cond,
                "n_aroma_timepoints": int(df[total_cols + cond_cols].notna().any(axis=1).sum()) if (total_cols or cond_cols) else 0,
                "n_positive_n_pulse_rows": int((pulse > 0).sum()) if len(pulse) else 0,
            }
        )
    return pd.DataFrame(rows)


def aggregate_coverage(lab_summary: pd.DataFrame, pilot_summary: pd.DataFrame, operational_protocol: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if not lab_summary.empty:
        n_batches = lab_summary["batch"].nunique()
        rows.append(
            {
                "dataset_or_campaign": "synthetic_vl3_camilos_archived_wines",
                "status": "existing process and wine chemistry; aroma GC pending",
                "n_fermentations": n_batches,
                "process_observations": int(lab_summary["n_process_obs"].sum()),
                "liquid_or_total_aroma_observations_current": int(lab_summary["n_aroma_total_obs"].sum()),
                "condensate_aroma_observations_current": int(lab_summary["n_aroma_condensate_obs"].sum()),
                "potential_terminal_target_aroma_observations": int(n_batches * len(TARGET_AROMAS)),
                "co2_sensor_points": 0,
                "model_use": "secondary validation / archive; weak for dynamic aroma calibration",
            }
        )
    if not pilot_summary.empty:
        rows.append(
            {
                "dataset_or_campaign": "pilot_2025_natural_must",
                "status": "existing pilot data with aroma and condensate",
                "n_fermentations": int(pilot_summary["batch"].nunique()),
                "process_observations": int(pilot_summary["n_process_obs"].sum()) if "n_process_obs" in pilot_summary else np.nan,
                "liquid_or_total_aroma_observations_current": int(pilot_summary["n_aroma_total_obs"].sum()),
                "condensate_aroma_observations_current": int(pilot_summary["n_aroma_condensate_obs"].sum()),
                "potential_terminal_target_aroma_observations": np.nan,
                "co2_sensor_points": int(pilot_summary.get("n_co2_sensor_points", pd.Series([0])).sum()),
                "model_use": "calibration evidence for aroma synthesis-loss and gas/liquid closure",
            }
        )
    if not operational_protocol.empty:
        if "full_samples" in operational_protocol:
            full_series = pd.to_numeric(operational_protocol["full_samples"], errors="coerce").fillna(0)
        else:
            full_series = pd.Series([0] * len(operational_protocol), dtype=float)
        full_samples = int(full_series.sum())
        n_exp = int(operational_protocol["candidate"].nunique()) if "candidate" in operational_protocol else len(operational_protocol)
        rows.append(
            {
                "dataset_or_campaign": "mbdoe_2026_operational_campaign",
                "status": "planned / partly under execution; designed for weak model directions",
                "n_fermentations": n_exp,
                "process_observations": np.nan,
                "liquid_or_total_aroma_observations_current": 0,
                "condensate_aroma_observations_current": 0,
                "potential_terminal_target_aroma_observations": int(n_exp * len(TARGET_AROMAS)),
                "planned_dynamic_target_aroma_observations": int(full_samples * len(TARGET_AROMAS)),
                "co2_sensor_points": "online",
                "model_use": "primary GC budget target for calibrating aroma layer and validating MPCC constraints",
            }
        )
    return pd.DataFrame(rows)


def build_priority_matrix() -> pd.DataFrame:
    rows = [
        {
            "option": "Analyze stored synthetic VL3/Camilo wines",
            "dynamic_aroma_trajectory": 0,
            "condensate_mass_closure": 0,
            "online_co2_stripping_driver": 0,
            "mbdoe_weak_direction_targeting": 0,
            "matrix_relevance_to_wine": 1,
            "existing_process_context": 2,
            "budget_efficiency_for_model": 1,
            "interpretation": "Useful as archive/independent endpoint check, but low leverage for synthesis-loss identifiability.",
        },
        {
            "option": "Use pilot 2025 natural aroma/condensate data",
            "dynamic_aroma_trajectory": 2,
            "condensate_mass_closure": 2,
            "online_co2_stripping_driver": 1,
            "mbdoe_weak_direction_targeting": 1,
            "matrix_relevance_to_wine": 2,
            "existing_process_context": 2,
            "budget_efficiency_for_model": 2,
            "interpretation": "Best existing calibration base; not fully optimized for current weak directions.",
        },
        {
            "option": "Prioritize MBDoE 2026 optimal campaign samples",
            "dynamic_aroma_trajectory": 2,
            "condensate_mass_closure": 2,
            "online_co2_stripping_driver": 2,
            "mbdoe_weak_direction_targeting": 3,
            "matrix_relevance_to_wine": 1,
            "existing_process_context": 2,
            "budget_efficiency_for_model": 3,
            "interpretation": "Highest marginal information per GC run for model calibration and future MPCC use.",
        },
    ]
    out = pd.DataFrame(rows)
    score_cols = [col for col in out.columns if col not in {"option", "interpretation"}]
    out["total_score"] = out[score_cols].sum(axis=1)
    return out


def build_templates() -> dict[str, pd.DataFrame]:
    return {
        "sample_traceability_template": pd.DataFrame(
            columns=[
                "sample_id",
                "campaign",
                "lot",
                "reactor_id",
                "treatment_id",
                "medium",
                "relative_time_h",
                "sample_datetime",
                "sample_type",
                "matrix",
                "volume_ml",
                "storage_temperature_c",
                "preservative_or_headspace",
                "operator",
                "chain_of_custody_step",
                "freezer_box_position",
                "gc_batch_id",
                "analysis_datetime",
                "use_in_calibration",
                "comments",
            ]
        ),
        "gc_qaqc_template": pd.DataFrame(
            [
                ("calibration_curve", "each GC batch", "R2 >= 0.99 or SOP limit", "repeat curve or flag batch"),
                ("internal_standard_response", "each sample", "within control chart", "dilute/reinject/flag"),
                ("method_blank", "beginning and end", "below LOQ", "subtract or repeat"),
                ("duplicate", ">= 1 per batch", "RSD within method threshold", "flag uncertainty"),
                ("matrix_spike", "representative wine and condensate", "70-130% recovery", "matrix correction or repeat"),
                ("carryover_check", "after high standard/sample", "below LOQ", "clean/inject blank/repeat"),
                ("sample_temperature_log", "from sampling to GC", "documented cold chain", "flag traceability deviation"),
            ],
            columns=["qa_qc_item", "frequency", "acceptance_criterion", "action_if_failed"],
        ),
        "aroma_analysis_request_template": pd.DataFrame(
            [
                ("ethyl_acetate", "wine", "dynamic sample", "ug/L or mg/L", "target MPCC state"),
                ("isoamyl_acetate", "wine", "dynamic sample", "ug/L", "target MPCC state"),
                ("ethyl_octanoate", "wine", "dynamic sample", "ug/L", "target MPCC state"),
                ("ethyl_acetate", "condensate", "terminal integral", "ug or ug/L equivalent", "loss closure"),
                ("isoamyl_acetate", "condensate", "terminal integral", "ug or ug/L equivalent", "loss closure"),
                ("ethyl_octanoate", "condensate", "terminal integral", "ug or ug/L equivalent", "loss closure"),
            ],
            columns=["compound", "matrix", "sample_role", "preferred_unit", "model_role"],
        ),
    }


def plot_coverage(coverage: pd.DataFrame) -> None:
    if coverage.empty:
        return
    labels = coverage["dataset_or_campaign"].str.replace("_", "\n")
    liquid = pd.to_numeric(coverage["liquid_or_total_aroma_observations_current"], errors="coerce").fillna(0)
    cond = pd.to_numeric(coverage["condensate_aroma_observations_current"], errors="coerce").fillna(0)
    planned = pd.to_numeric(coverage.get("planned_dynamic_target_aroma_observations", pd.Series([0] * len(coverage))), errors="coerce").fillna(0)
    terminal = pd.to_numeric(coverage["potential_terminal_target_aroma_observations"], errors="coerce").fillna(0)
    x = np.arange(len(coverage))
    width = 0.2
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - 1.5 * width, liquid, width, label="current liquid/total aroma obs", color="#3B82F6")
    ax.bar(x - 0.5 * width, cond, width, label="current condensate aroma obs", color="#10B981")
    ax.bar(x + 0.5 * width, planned, width, label="planned dynamic target obs", color="#F59E0B")
    ax.bar(x + 1.5 * width, terminal, width, label="potential terminal target obs", color="#6B7280")
    ax.set_ylabel("observation count")
    ax.set_title("Aroma analytical information currently available vs planned")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_01_aroma_coverage_current_vs_planned.png", dpi=180)
    plt.close(fig)


def plot_priority_matrix(priority: pd.DataFrame) -> None:
    score_cols = [
        "dynamic_aroma_trajectory",
        "condensate_mass_closure",
        "online_co2_stripping_driver",
        "mbdoe_weak_direction_targeting",
        "matrix_relevance_to_wine",
        "existing_process_context",
        "budget_efficiency_for_model",
    ]
    data = priority.set_index("option")[score_cols]
    fig, ax = plt.subplots(figsize=(12, 4.2))
    im = ax.imshow(data.to_numpy(dtype=float), cmap="YlGnBu", vmin=0, vmax=3, aspect="auto")
    ax.set_xticks(np.arange(len(score_cols)))
    ax.set_xticklabels([col.replace("_", "\n") for col in score_cols], fontsize=8)
    ax.set_yticks(np.arange(len(data.index)))
    ax.set_yticklabels(data.index, fontsize=9)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data.iloc[i, j]:.0f}", ha="center", va="center", color="#111827", fontsize=9)
    ax.set_title("Decision support score for GC budget allocation")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="score")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_02_gc_budget_priority_matrix.png", dpi=180)
    plt.close(fig)


def plot_model_selection(aroma_model_selection: pd.DataFrame) -> None:
    if aroma_model_selection.empty or "selection_score" not in aroma_model_selection:
        return
    df = aroma_model_selection.sort_values("selection_score", ascending=True).copy()
    colors = np.where(df.get("selected", False).astype(bool), "#10B981", "#6B7280")
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh(df["model"], df["selection_score"], color=colors)
    ax.set_xlabel("selection score (lower is better)")
    ax.set_title("Pilot 2025 aroma model benchmark")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_03_aroma_model_selection_score.png", dpi=180)
    plt.close(fig)


def plot_estimability_current_after(current: pd.DataFrame, after: pd.DataFrame) -> None:
    if current.empty or after.empty:
        return
    target = [
        "k_EA_growth",
        "k_EA_stationary",
        "k_IAA_growth",
        "k_IAA_stationary",
        "k_EO_growth",
        "k_EO_stationary",
        "alpha_EA_loss",
        "alpha_IAA_loss",
        "alpha_EO_loss",
        "k_EA_XE_Nlim",
    ]
    cur = current[current["parameter"].isin(target)][["parameter", "approx_95_multiplier"]].rename(
        columns={"approx_95_multiplier": "current"}
    )
    aft = after[after["parameter"].isin(target)][["parameter", "approx_95_multiplier"]].rename(
        columns={"approx_95_multiplier": "after_campaign"}
    )
    df = cur.merge(aft, on="parameter", how="outer").set_index("parameter").reindex(target)
    df = df.apply(pd.to_numeric, errors="coerce")
    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(df.index))
    width = 0.38
    ax.bar(x - width / 2, df["current"].clip(upper=1e4), width, label="current", color="#EF4444")
    ax.bar(x + width / 2, df["after_campaign"].clip(upper=1e4), width, label="after MBDoE campaign", color="#10B981")
    ax.set_yscale("log")
    ax.set_ylabel("approx. 95% multiplier (log scale, clipped at 1e4)")
    ax.set_title("Aroma-layer parameter uncertainty before and after campaign")
    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=45, ha="right", fontsize=8)
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_04_aroma_estimability_current_after.png", dpi=180)
    plt.close(fig)


def plot_var_reduction(reduction: pd.DataFrame) -> None:
    if reduction.empty:
        return
    df = reduction.copy()
    if "campaign_var_reduction" not in df:
        return
    df = df[df["group"].isin(["synthesis", "fermentation"])].copy()
    colors = df["group"].map({"synthesis": "#8B5CF6", "fermentation": "#3B82F6"}).fillna("#6B7280")
    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(df))
    ax.bar(x, df["campaign_var_reduction"], color=colors)
    ax.axhline(0.9, color="#111827", lw=1, ls="--", label="90% variance reduction")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("variance reduction")
    ax.set_title("Expected variance reduction from selected aroma MBDoE campaign")
    ax.set_xticks(x)
    ax.set_xticklabels(df["parameter"], rotation=60, ha="right", fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_05_mbdoe_parameter_variance_reduction.png", dpi=180)
    plt.close(fig)


def plot_dopt_benchmark(dopt: pd.DataFrame) -> None:
    if dopt.empty:
        return
    metrics = ["logdet", "condition_number", "trace_inv", "mean_var_reduction", "worst_var_reduction"]
    available = [metric for metric in metrics if metric in dopt.columns]
    if not available:
        return
    fig, axes = plt.subplots(1, len(available), figsize=(3.2 * len(available), 4))
    if len(available) == 1:
        axes = [axes]
    for ax, metric in zip(axes, available):
        values = pd.to_numeric(dopt[metric], errors="coerce")
        ax.bar(dopt["design"], values, color=["#10B981", "#6B7280", "#6B7280"][: len(dopt)])
        ax.set_title(metric.replace("_", " "))
        if metric == "condition_number":
            ax.set_yscale("log")
        ax.tick_params(axis="x", labelrotation=65, labelsize=8)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Hybrid selected design vs pure D-opt benchmark", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_06_hybrid_vs_dopt_benchmark.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def copy_reference_figures() -> None:
    figure_map = {
        A07_BUNDLE / "figures" / "fig_01_aroma_model_selection.png": "source_A07_fig_01_aroma_model_selection.png",
        A07_BUNDLE / "figures" / "fig_04_final_fit_relative_rmse.png": "source_A07_fig_04_final_fit_relative_rmse.png",
        A07_BUNDLE / "figures" / "fig_05_global_estimability_current.png": "source_A07_fig_05_global_estimability_current.png",
        A07_BUNDLE / "figures" / "fig_06_global_estimability_after_campaign.png": "source_A07_fig_06_global_estimability_after_campaign.png",
        A07_BUNDLE / "figures" / "fig_09_global_selected_campaign_inputs.png": "source_A07_fig_09_global_selected_campaign_inputs.png",
        A07_BUNDLE / "figures" / "fig_15_final_fit_representative_25170.png": "source_A07_fig_15_final_fit_representative_25170.png",
        A15_BUNDLE / "figures" / "fig_12_pilot_aroma_model_selection.png": "source_A15_fig_12_pilot_aroma_model_selection.png",
        AROMA_DOE / "aroma_eigen_spectrum.png": "source_aroma_doe_eigen_spectrum.png",
        AROMA_DOE / "aroma_weak_eigendirection_loadings.png": "source_aroma_doe_weak_direction_loadings.png",
    }
    for src, name in figure_map.items():
        copy_if_exists(src, FIG_DIR / name)


def write_tables(tables: dict[str, pd.DataFrame]) -> None:
    for name, table in tables.items():
        if isinstance(table, pd.DataFrame):
            table.to_csv(TABLE_DIR / f"{name}.csv", index=False)


def write_templates(templates: dict[str, pd.DataFrame]) -> None:
    for name, table in templates.items():
        table.to_csv(TEMPLATE_DIR / f"{name}.csv", index=False)


def write_source_trace(paths: dict[str, Path | None], snippets: pd.DataFrame, source_files: pd.DataFrame) -> None:
    source_files.to_csv(TRACE_DIR / "source_file_manifest.csv", index=False)
    snippets.to_csv(TRACE_DIR / "source_document_snippets.csv", index=False)
    data_sources = pd.DataFrame(
        [
            {
                "artifact": "laboratory synthetic VL3 workbook",
                "path": str(LAB_SYNTHETIC),
                "role": "existing synthetic/wine chemistry campaign; aroma GC not executed in workbook",
            },
            {
                "artifact": "pilot 2025 calibration workbook",
                "path": str(PILOT_WORKBOOK),
                "role": "existing natural-must pilot data with aroma total, condensate and some CO2",
            },
            {
                "artifact": "A07 aroma NLP calibration bundle",
                "path": str(A07_BUNDLE),
                "role": "aroma model selection, calibration, estimability and MBDoE evidence",
            },
            {
                "artifact": "A15 pilot scalability bundle",
                "path": str(A15_BUNDLE),
                "role": "pilot transferability and scale-up evidence",
            },
            {
                "artifact": "aroma joint campaign DOE results",
                "path": str(AROMA_DOE),
                "role": "selected MBDoE campaign, variance reduction, eigenanalysis and D-opt benchmark",
            },
            {
                "artifact": "operational volume-constrained DOE protocol",
                "path": str(OPERATIONAL_DOE),
                "role": "operational sampling schedule, full-sample volume and reactor-volume feasibility",
            },
        ]
    )
    data_sources.to_csv(TRACE_DIR / "repository_evidence_index.csv", index=False)


def build_claims_limitations(priority: pd.DataFrame, dopt: pd.DataFrame) -> str:
    dopt_text = "not available"
    if not dopt.empty:
        hybrid = dopt[dopt["design"].eq("hybrid_current")]
        greedy = dopt[dopt["design"].eq("dopt_greedy")]
        if not hybrid.empty and not greedy.empty:
            h = hybrid.iloc[0]
            g = greedy.iloc[0]
            dopt_text = (
                f"hybrid logdet={h['logdet']:.3f}, trace_inv={h['trace_inv']:.3f}, "
                f"worst_var_reduction={h['worst_var_reduction']:.3f}; pure D-opt logdet={g['logdet']:.3f}, "
                f"trace_inv={g['trace_inv']:.3f}, worst_var_reduction={g['worst_var_reduction']:.3f}"
            )
    top_option = priority.sort_values("total_score", ascending=False).iloc[0]["option"] if not priority.empty else "MBDoE campaign"
    return f"""# Claims y limitaciones A04

## Claims soportados

1. El set sintetico VL3/Camilo contiene informacion de fermentacion y vino suficiente para trazabilidad y comparacion quimica general, pero no contiene mediciones aromaticas GC ejecutadas en el workbook de trabajo.
2. Los datos piloto 2025 de mosto natural si contienen mediciones de aromas totales y condensado, por lo que son la base existente mas fuerte para ajustar una capa dinamica de sintesis-perdida aromatica.
3. La estructura seleccionada para aroma, `ea_ethanol_nlimited`, fue elegida por benchmark determinista de modelos y supera a la estructura base por selection score/BIC/WSSE en piloto.
4. La campana MBDoE 2026 apunta directamente a direcciones debiles del modelo. En el benchmark disponible: {dopt_text}.
5. Con presupuesto GC limitado, la opcion prioritaria es: `{top_option}`. Esta decision combina evidencia FIM/MBDoE, disponibilidad de CO2, cierre por condensado y valor marginal para el MPCC.

## Limitaciones que deben declararse

1. A04 todavia no debe afirmar resultados aromaticos finales de la campana sintetica VL3/Camilo si las muestras permanecen almacenadas y no analizadas por GC.
2. Las puntuaciones de priorizacion de presupuesto GC son una matriz de decision documentada, no un criterio de optimalidad Pyomo.DoE por si solo.
3. La FIM usada en los resultados MBDoE depende del modelo nominal y de escalas de error asumidas; debe actualizarse con datos reales de la campana ejecutada.
4. Ethyl acetate sigue siendo el marcador mas desafiante, con sesgo/errores mayores y parametros parcialmente confundidos; se recomienda reportarlo explicitamente.
5. El condensado terminal entrega una restriccion integral de perdida volatilizada, pero no reemplaza mediciones GC online de fase gas.
6. La transferencia desde mosto sintetico a mosto natural requiere validacion; por eso se conserva evidencia piloto natural y al menos un brazo natural en protocolos posteriores.
"""


def build_report(
    tables: dict[str, pd.DataFrame],
    snippets: pd.DataFrame,
    source_files: pd.DataFrame,
    created_at: str,
) -> str:
    coverage = tables["coverage_summary"]
    priority = tables["gc_budget_priority_matrix"]
    lab_summary = tables["lab_synthetic_summary"]
    pilot_summary = tables["pilot_summary"]
    aroma_model = tables["a07_aroma_model_selection"]
    final_fit = tables["a07_final_fit_metrics"]
    reduction = tables["aroma_campaign_parameter_reduction"]
    dopt = tables["aroma_dopt_benchmark_summary"]
    operational = tables["operational_campaign_protocol"]

    selected_model = "not available"
    if not aroma_model.empty and "selected" in aroma_model:
        selected_rows = aroma_model[aroma_model["selected"].astype(bool)]
        if not selected_rows.empty:
            row = selected_rows.iloc[0]
            selected_model = (
                f"{row['model']} (data WSSE={fmt(row['data_wsse'])}, BIC={fmt(row['bic'])}, "
                f"selection score={fmt(row['selection_score'])})"
            )

    synthesis_reduction = reduction[reduction["group"].eq("synthesis")] if not reduction.empty and "group" in reduction else pd.DataFrame()
    min_synth_reduction = (
        float(synthesis_reduction["campaign_var_reduction"].min()) if not synthesis_reduction.empty else np.nan
    )
    mean_synth_reduction = (
        float(synthesis_reduction["campaign_var_reduction"].mean()) if not synthesis_reduction.empty else np.nan
    )

    report = f"""# A04 - Analisis de vino y condensados de aroma (Actividad 3.11)

Generated: {created_at}

## 1. Objetivo del bundle

Este bundle entrega evidencia tecnica para sostener A04 bajo una situacion operacional concreta: existe una campana sintetica VL3/Camilo con buena quimica de vino y trazabilidad de fermentacion, pero sin analisis aromatico GC ejecutado en la base de trabajo. Las muestras de vino fueron consideradas/almacenadas para analisis posterior, pero el presupuesto GC debe priorizarse.

La decision tecnica documentada aqui es priorizar el analisis de muestras provenientes de la campana experimental optima basada en modelo, y usar las muestras historicas de Camilo como material de validacion secundaria si queda presupuesto analitico. La razon no es que la campana historica carezca de valor, sino que su diseno no fue construido para excitar las direcciones debiles del modelo de aromas ni para cerrar el balance liquido-condensado con CO2.

## 2. Evidencia fuente y estado de trazabilidad

Fuentes primarias indexadas:

{md_table(source_files, columns=["source_id", "role", "file_name", "size_mb", "modified"], max_rows=10)}

Los extractos documentales se guardan en `source_trace/source_document_snippets.csv`. Estos extractos no sustituyen el documento original; solo fijan ubicaciones y contexto para auditoria.

Extractos relevantes:

{md_table(snippets, max_rows=8)}

## 3. Disponibilidad de datos

La base `mosto_sintetico_vl3.xlsx` contiene multiples fermentaciones sinteticas con temperatura, densidad, biomasa, glucosa, fructosa, nitrogeno, glicerol, piruvato, acetaldehido, etanol y pulsos. En el workbook actual, las columnas de aroma no aportan observaciones GC utiles para calibrar la capa aromatica.

La base piloto 2025 (`Calibration_data_vl3.xlsx`) contiene mosto natural con mediciones de aromas totales y condensado. Esta es la evidencia existente mas fuerte para modelar particion liquido-gas y perdida volatilizada.

Resumen comparativo:

{md_table(coverage, max_rows=10)}

Resumen sintetico VL3/Camilo:

{md_table(lab_summary, max_rows=12, columns=["batch", "n_rows", "t_max_h", "temperature_min_c", "temperature_max_c", "S_initial_g_l", "S_final_g_l", "YAN_initial_mg_l", "n_process_obs", "n_aroma_total_obs", "n_aroma_condensate_obs"])}

Resumen piloto 2025:

{md_table(pilot_summary, max_rows=12, columns=["batch", "n_rows", "t_max_h", "temperature_min_c", "temperature_max_c", "S_initial_g_l", "S_final_g_l", "YAN_initial_mg_l", "n_aroma_total_obs", "n_aroma_condensate_obs", "n_co2_sensor_points"])}

## 4. Modelo dinamico usado para justificar la decision

La capa de fermentacion primaria mantiene estados de biomasa viable, biomasa muerta, nitrogeno asimilable, glucosa, fructosa, etanol y glicerol:

$$
\\frac{{d x}}{{dt}} = f(x, u, \\theta), \\quad
x = [X, X_d, N, G, F, E, Gly]^T
$$

La capa secundaria incluye piruvato, acetaldehido y acido acetico como proxies quimicos ligados al metabolismo central. Para A04, el componente critico es la capa aromatica. Para cada aroma objetivo \\(i\\):

$$
\\frac{{d A_i^{{liq}}}}{{dt}} = r_i^{{prod}} - r_i^{{loss}}
$$

$$
\\frac{{d A_i^{{cond}}}}{{dt}} = r_i^{{loss}}
$$

$$
A_i^{{total}} = A_i^{{liq}} + A_i^{{cond}}
$$

La produccion se representa como funcion de la actividad fermentativa y del estado de limitacion de nitrogeno:

$$
r_i^{{prod}} =
\\left(k_{{i,g}}\\,\\phi_N + k_{{i,s}}(1-\\phi_N)\\right) q_S
$$

$$
\\phi_N = \\frac{{N}}{{N + K_N}}
$$

Para ethyl acetate se mantiene una extension dependiente de etanol, biomasa y limitacion de nitrogeno:

$$
r_{{EA}}^{{prod}} =
r_{{EA}}^{{phase}} +
k_{{EA,XE}} \\frac{{X E}}{{K_E + E}} +
k_{{EA,XE,Nlim}} \\frac{{X E}}{{K_E + E}}\\frac{{K_N}}{{K_N + N}}
$$

La perdida hacia condensado se modela como stripping impulsado por CO2 y particion gas-liquido:

$$
r_i^{{loss}} =
\\alpha_i^{{loss}}\\,K_i^{{LG}}(T,E,G,F)\\,q_{{gas}}\\,A_i^{{liq}}
$$

Donde \\(K_i^{{LG}}\\) se fija principalmente desde literatura/UNIFAC y \\(\\alpha_i^{{loss}}\\) absorbe desviaciones de transferencia no capturadas por el equilibrio idealizado.

## 5. Calibracion, benchmarking y calidad de ajuste

El benchmark piloto selecciono la estructura aromatica:

`{selected_model}`

Tabla de seleccion de modelos:

{md_table(aroma_model, max_rows=8, columns=["model", "n_parameters", "data_wsse", "wsse_per_data_residual", "bic", "active_bound_count", "selected", "selection_score"])}

Metricas integradas del ajuste piloto:

{md_table(final_fit, max_rows=24, columns=["group", "state", "species", "state_pool", "n", "relative_rmse", "relative_bias", "corr"])}

Interpretacion:

- Isoamyl acetate total/retained queda en rango util para modelamiento dinamico.
- Ethyl octanoate queda razonablemente informativo, aunque condensado mantiene incertidumbre.
- Ethyl acetate es el marcador mas desafiante y justifica una estructura enriquecida con etanol, biomasa y limitacion de nitrogeno.
- Las muestras de condensado son importantes porque convierten la volatilizacion desde un termino libre a una restriccion integral observable.

## 6. Estimabilidad y MBDoE

La FIM se construye a partir de sensibilidades locales:

$$
F(\\theta, d) =
J(\\theta, d)^T W J(\\theta, d)
$$

con \\(J\\) la matriz de sensibilidades de salidas medidas frente a parametros y \\(W\\) la matriz de pesos/varianzas de medicion. La campana seleccionada maximiza informacion manteniendo balance entre:

- D-optimalidad: maximizar \\(\\log \\det(F)\\).
- Robustez de direccion debil: evitar que el menor autovalor quede excesivamente bajo.
- Menor incertidumbre posterior aproximada: reducir \\(\\mathrm{{trace}}(F^{{-1}})\\) y varianzas diagonales.

El diseno MBDoE de aromas predice reducciones de varianza altas en parametros de sintesis. Para los parametros de sintesis aromatica, la reduccion minima es {fmt(min_synth_reduction)} y la reduccion media es {fmt(mean_synth_reduction)}.

Reduccion esperada por parametro:

{md_table(reduction, max_rows=30)}

Benchmark hibrido vs D-opt puro:

{md_table(dopt, max_rows=5)}

La lectura tecnica es que el D-opt puro obtiene un logdet levemente mayor, pero el diseno hibrido seleccionado conserva mejor peor reduccion de varianza y menor `trace_inv`, por lo que es mas defendible bajo incertidumbre operacional y presupuesto GC limitado.

## 7. Campana operacionalizada para analisis GC

La campana operacional con restriccion de volumen traduce el diseno a fermentaciones, muestras completas y muestras pequenas. Las muestras completas contienen etanol + aromas + panel completo; las muestras pequenas reducen carga de volumen cuando no se requiere GC completo.

{md_table(operational, max_rows=12, columns=["campaign_order", "lot", "candidate", "medium", "family", "horizon_h", "temperature_segments", "full_samples", "small_samples", "planned_sample_volume_ml", "rationale"])}

Esta operacionalizacion es clave para A04 porque evita gastar GC en muestras que no resuelven direcciones debiles del modelo. En particular, la campana MBDoE integra:

- trayectoria liquida de aromas en vino/mosto;
- condensado terminal como restriccion integral de volatilizacion;
- CO2 online como driver de stripping;
- covariables de metabolismo central que conectan el modelo ODE con MPCC.

## 8. Decision de priorizacion analitica

La matriz siguiente no reemplaza la FIM; documenta el criterio experimental/analitico usado para traducir FIM, trazabilidad y presupuesto a una decision operable.

{md_table(priority, max_rows=10)}

Conclusion de decision:

1. Analizar primero muestras de la campana MBDoE operacionalizada.
2. Usar piloto 2025 como base de calibracion/benchmark ya disponible.
3. Mantener muestras historicas sinteticas VL3/Camilo como validacion secundaria o respaldo de endpoint si hay presupuesto GC remanente.

## 9. QA/QC minimo recomendado para A04

{md_table(tables["gc_qaqc_plan"], max_rows=12)}

El bundle incluye plantillas editables:

- `templates/sample_traceability_template.csv`
- `templates/gc_qaqc_template.csv`
- `templates/aroma_analysis_request_template.csv`

## 10. Resultados y figuras incluidas

Figuras generadas directamente en este bundle:

- `figures/fig_01_aroma_coverage_current_vs_planned.png`
- `figures/fig_02_gc_budget_priority_matrix.png`
- `figures/fig_03_aroma_model_selection_score.png`
- `figures/fig_04_aroma_estimability_current_after.png`
- `figures/fig_05_mbdoe_parameter_variance_reduction.png`
- `figures/fig_06_hybrid_vs_dopt_benchmark.png`

Figuras fuente copiadas desde A07/A15/DOE:

- seleccion de modelo aroma piloto;
- RMSE relativo del ajuste integrado;
- estimabilidad antes/despues de campana;
- inputs de campana seleccionada;
- ajuste representativo piloto 25170;
- espectro de autovalores y direccion debil del diseno de aromas.

## 11. Conclusion

Para A04, la evidencia mas honesta es separar tres niveles:

1. **Trazabilidad y vino sintetico existente:** existe material fermentativo/quimico y muestras almacenadas, pero sin GC aromatico ejecutado.
2. **Evidencia aroma ya calibrable:** piloto 2025 aporta aromas totales, condensado y CO2 parcial; por eso sostiene la estructura dinamica de aroma.
3. **Decision de gasto GC:** la campana MBDoE 2026 entrega mayor informacion marginal para parametros aromaticos y para el uso futuro del modelo dentro del MPCC.

Por lo tanto, la postergacion o repriorizacion de analisis GC sobre muestras historicas no es una brecha metodologica; es una decision tecnica basada en estimabilidad, cierre de balance y uso eficiente de presupuesto analitico.
"""
    return report


def write_readme(created_at: str) -> str:
    return f"""# README - Bundle A04 vino y condensados de aroma

Generated: {created_at}

## Contenido principal

- `A04_report.md`: reporte tecnico autoexplicativo para sostener A04.
- `CLAIMS_Y_LIMITACIONES_A04.md`: claims permitidos, limitaciones y cautelas de interpretacion.
- `tables/coverage_summary.csv`: comparacion compacta entre campana sintetica VL3/Camilo, piloto 2025 y campana MBDoE operacional.
- `tables/gc_budget_priority_matrix.csv`: matriz de decision para priorizar presupuesto GC.
- `tables/a07_aroma_model_selection.csv`: benchmark de modelos aromaticos en piloto.
- `tables/aroma_campaign_parameter_reduction.csv`: reduccion esperada de varianza por parametro.
- `tables/aroma_dopt_benchmark_summary.csv`: benchmark hibrido vs D-opt puro.
- `figures/`: figuras generadas y figuras fuente copiadas desde bundles A07/A15/DOE.
- `source_trace/`: manifiesto de fuentes, rutas originales, hashes y extractos documentales.
- `templates/`: plantillas CSV para trazabilidad de muestras, solicitud GC y QA/QC.

## Lectura recomendada

1. Revisar `A04_report.md` completo.
2. Usar `CLAIMS_Y_LIMITACIONES_A04.md` para redactar el anexo sin sobreafirmar resultados GC no ejecutados.
3. Abrir `figures/fig_01_aroma_coverage_current_vs_planned.png` y `figures/fig_02_gc_budget_priority_matrix.png` para explicar la decision analitica.
4. Respaldar el argumento cuantitativo con `figures/fig_05_mbdoe_parameter_variance_reduction.png` y `figures/fig_06_hybrid_vs_dopt_benchmark.png`.

## Mensaje central

La campana historica sintetica tiene valor como trazabilidad y validacion secundaria, pero la campana MBDoE operacionalizada es la prioridad racional para gasto GC porque combina trayectoria de aromas, condensado terminal, CO2 y excitacion de parametros debiles.
"""


def archive_bundle(project_a04_dir: Path | None) -> Path:
    zip_path = RESULTS_DIR / f"{BUNDLE_NAME}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(BUNDLE_DIR.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(BUNDLE_DIR.parent))
    if project_a04_dir is not None:
        project_a04_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(zip_path, project_a04_dir / zip_path.name)
        for name in ["A04_report.md", "README_BUNDLE.md", "CLAIMS_Y_LIMITACIONES_A04.md"]:
            copy_if_exists(BUNDLE_DIR / name, project_a04_dir / name)
    return zip_path


def main() -> None:
    created_at = datetime.now().isoformat(timespec="seconds")
    prepare_dirs()
    project_paths = discover_project_paths()

    lab_summary = workbook_campaign_summary(LAB_SYNTHETIC, "synthetic_vl3_camilos_archived_wines", r"^MS\d+")
    pilot_summary = workbook_campaign_summary(PILOT_WORKBOOK, "pilot_2025_natural_must", r"^\d+")
    a15_pilot_summary = read_csv(A15_BUNDLE / "tables" / "03_pilot_batch_summary.csv")
    if not a15_pilot_summary.empty:
        pilot_summary = a15_pilot_summary.copy()
    operational = read_csv(OPERATIONAL_DOE / "campaign_protocol.csv")
    coverage = aggregate_coverage(lab_summary, pilot_summary, operational)

    source_files = source_manifest(project_paths)
    snippets = pd.concat(
        [
            extract_pdf_snippets(project_paths.get("reference_pdf")),
            extract_docx_snippets(project_paths.get("gc_protocol_docx")),
            extract_docx_snippets(project_paths.get("pilot_2025_docx"), max_rows=8),
        ],
        ignore_index=True,
    )
    if not snippets.empty and "relevance_score" in snippets:
        snippets = snippets.sort_values(["relevance_score", "source", "location"], ascending=[False, True, True])

    priority = build_priority_matrix()
    templates = build_templates()

    tables = {
        "lab_synthetic_summary": lab_summary,
        "pilot_summary": pilot_summary,
        "coverage_summary": coverage,
        "source_document_snippets": snippets,
        "source_file_manifest": source_files,
        "gc_budget_priority_matrix": priority,
        "operational_campaign_protocol": operational,
        "operational_schedule": read_csv(OPERATIONAL_DOE / "operational_schedule.csv"),
        "a07_aroma_model_selection": read_csv(A07_BUNDLE / "tables" / "01_aroma_model_selection.csv"),
        "a07_final_fit_metrics": read_csv(A07_BUNDLE / "tables" / "06_final_fit_metrics.csv"),
        "a07_estimability_current": read_csv(A07_BUNDLE / "tables" / "07_global_parameter_estimability_current.csv"),
        "a07_estimability_after_campaign": read_csv(A07_BUNDLE / "tables" / "08_global_parameter_estimability_after_campaign.csv"),
        "a07_target_parameters": pd.DataFrame({"source": [str(A07_BUNDLE / "tables" / "15_target_parameters.json")]}),
        "a15_transferability_decision_matrix": read_csv(A15_BUNDLE / "tables" / "22_transferability_decision_matrix.csv"),
        "aroma_campaign_selected": read_csv(AROMA_DOE / "aroma_campaign_selected.csv"),
        "aroma_campaign_parameter_reduction": read_csv(AROMA_DOE / "aroma_campaign_parameter_reduction.csv"),
        "aroma_dopt_benchmark_summary": read_csv(AROMA_DOE / "aroma_dopt_benchmark_summary.csv"),
        "aroma_eigen_summary": read_csv(AROMA_DOE / "aroma_eigen_summary.csv"),
        "gc_qaqc_plan": templates["gc_qaqc_template"],
    }

    write_tables(tables)
    write_templates(templates)
    write_source_trace(project_paths, snippets, source_files)

    plot_coverage(coverage)
    plot_priority_matrix(priority)
    plot_model_selection(tables["a07_aroma_model_selection"])
    plot_estimability_current_after(tables["a07_estimability_current"], tables["a07_estimability_after_campaign"])
    plot_var_reduction(tables["aroma_campaign_parameter_reduction"])
    plot_dopt_benchmark(tables["aroma_dopt_benchmark_summary"])
    copy_reference_figures()

    report = build_report(tables, snippets, source_files, created_at)
    (BUNDLE_DIR / "A04_report.md").write_text(report, encoding="utf-8")
    (BUNDLE_DIR / "README_BUNDLE.md").write_text(write_readme(created_at), encoding="utf-8")
    (BUNDLE_DIR / "CLAIMS_Y_LIMITACIONES_A04.md").write_text(
        build_claims_limitations(priority, tables["aroma_dopt_benchmark_summary"]), encoding="utf-8"
    )
    copy_if_exists(Path(__file__), SNAPSHOT_DIR / Path(__file__).name)

    manifest_rows = []
    for path in sorted(BUNDLE_DIR.rglob("*")):
        if path.is_file():
            manifest_rows.append(
                {
                    "relative_path": str(path.relative_to(BUNDLE_DIR)).replace("\\", "/"),
                    "size_bytes": path.stat().st_size,
                }
            )
    pd.DataFrame(manifest_rows).to_csv(BUNDLE_DIR / "bundle_manifest.csv", index=False)

    zip_path = archive_bundle(project_paths.get("project_a04_dir"))
    print(f"[done] bundle_dir={BUNDLE_DIR}")
    print(f"[done] zip={zip_path}")
    if project_paths.get("project_a04_dir") is not None:
        print(f"[done] copied_to={project_paths['project_a04_dir']}")


if __name__ == "__main__":
    main()
