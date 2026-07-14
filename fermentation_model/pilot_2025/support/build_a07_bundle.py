from __future__ import annotations

import json
import math
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PILOT_DIR = SCRIPT_DIR.parent if SCRIPT_DIR.name == "support" else SCRIPT_DIR
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
RESULTS = PILOT_DIR / "results" / "co2_solubility_integrated_doe"
AROMA_RESULTS = PILOT_DIR / "results" / "aroma_model_selection_doe"
GLOBAL_RESULTS = PILOT_DIR / "results" / "global_state_model_selection_doe"
BUNDLE_ROOT = PILOT_DIR / "bundles" / "A07_aroma_nlp_calibration"
ZIP_PATH = PILOT_DIR / "bundles" / "A07_aroma_nlp_calibration.zip"


def read_text(path: Path, default: str = "") -> str:
    if path.exists():
        return path.read_text(encoding="utf-8").strip()
    return default


def safe_read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, **kwargs)


def md_table(df: pd.DataFrame, columns: list[str] | None = None, n: int | None = None) -> str:
    if df.empty:
        return "_No disponible._"
    out = df.copy()
    if columns:
        out = out[[col for col in columns if col in out.columns]]
    if n is not None:
        out = out.head(n)
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda x: f"{x:.4g}" if pd.notna(x) else "")
    widths = [max(len(str(col)), *(len(str(v)) for v in out[col].astype(str))) for col in out.columns]
    header = "| " + " | ".join(str(col).ljust(widths[i]) for i, col in enumerate(out.columns)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(widths))) + " |"
    rows = [
        "| " + " | ".join(str(row[col]).ljust(widths[i]) for i, col in enumerate(out.columns)) + " |"
        for _, row in out.iterrows()
    ]
    return "\n".join([header, sep, *rows])


def ensure_clean_bundle() -> tuple[Path, Path, Path, Path]:
    bundles_dir = PILOT_DIR / "bundles"
    bundles_dir.mkdir(parents=True, exist_ok=True)
    if BUNDLE_ROOT.exists():
        resolved = BUNDLE_ROOT.resolve()
        allowed = bundles_dir.resolve()
        if allowed not in resolved.parents and resolved != allowed:
            raise RuntimeError(f"Refusing to remove unexpected path: {BUNDLE_ROOT}")
        shutil.rmtree(BUNDLE_ROOT)
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    figures = BUNDLE_ROOT / "figures"
    tables = BUNDLE_ROOT / "tables"
    source = BUNDLE_ROOT / "source_snapshot"
    for folder in (figures, tables, source):
        folder.mkdir(parents=True, exist_ok=True)
    return figures, tables, source, BUNDLE_ROOT


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def load_core_tables() -> dict[str, pd.DataFrame | dict | str]:
    co2_model = read_text(RESULTS / "selected_co2_model.txt", "not_available")
    aroma_model = read_text(RESULTS / "selected_aroma_model_inherited.txt", "not_available")
    secondary_model = read_text(RESULTS / "selected_secondary_model_inherited.txt", "not_available")
    tables: dict[str, pd.DataFrame | dict | str] = {
        "co2_model": co2_model,
        "aroma_model": aroma_model,
        "secondary_model": secondary_model,
        "co2_selection": safe_read_csv(RESULTS / "co2_model_selection_summary.csv"),
        "co2_fit_selected": safe_read_csv(RESULTS / f"fit_co2_{co2_model}.csv"),
        "co2_metrics": safe_read_csv(RESULTS / "co2_metrics_selected.csv"),
        "co2_params": safe_read_csv(RESULTS / "co2_params_selected.csv", header=None, names=["parameter", "value"]),
        "fit_integrated": safe_read_csv(RESULTS / "fit_integrated_selected.csv"),
        "final_metrics": safe_read_csv(RESULTS / "final_fit_metrics.csv"),
        "initial_metrics": safe_read_csv(RESULTS / "initial_fit_metrics.csv"),
        "aroma_selection": safe_read_csv(AROMA_RESULTS / "model_selection_summary.csv"),
        "secondary_selection": safe_read_csv(GLOBAL_RESULTS / "secondary_model_selection_summary.csv"),
        "global_estimability": safe_read_csv(GLOBAL_RESULTS / "parameter_estimability_global_current.csv"),
        "global_estimability_after": safe_read_csv(GLOBAL_RESULTS / "parameter_estimability_global_current_plus_campaign.csv"),
        "global_campaign": safe_read_csv(GLOBAL_RESULTS / "selected_campaign_hybrid_global_selected.csv"),
        "co2_o2_estimability": safe_read_csv(RESULTS / "co2_o2_parameter_estimability_current.csv"),
        "co2_o2_estimability_after": safe_read_csv(RESULTS / "co2_o2_parameter_estimability_current_plus_campaign.csv"),
        "co2_o2_campaign": safe_read_csv(RESULTS / "co2_o2_selected_campaign_hybrid.csv"),
        "co2_o2_ranking": safe_read_csv(RESULTS / "co2_o2_candidate_ranking.csv"),
    }
    fim_metrics = {}
    for key, path in {
        "global_current": GLOBAL_RESULTS / "fim_metrics_global_current.json",
        "co2_o2_current": RESULTS / "co2_o2_fim_metrics_current.json",
    }.items():
        if path.exists():
            fim_metrics[key] = json.loads(path.read_text(encoding="utf-8"))
    tables["fim_metrics"] = fim_metrics
    return tables


def write_curated_tables(tables: dict[str, pd.DataFrame | dict | str], out: Path) -> None:
    mapping = {
        "01_aroma_model_selection.csv": tables["aroma_selection"],
        "02_secondary_model_selection.csv": tables["secondary_selection"],
        "03_co2_model_selection.csv": tables["co2_selection"],
        "04_co2_selected_fit_multistart.csv": tables["co2_fit_selected"],
        "05_integrated_aroma_fit.csv": tables["fit_integrated"],
        "06_final_fit_metrics.csv": tables["final_metrics"],
        "07_global_parameter_estimability_current.csv": tables["global_estimability"],
        "08_global_parameter_estimability_after_campaign.csv": tables["global_estimability_after"],
        "09_co2_parameter_estimability_current.csv": tables["co2_o2_estimability"],
        "10_co2_parameter_estimability_after_campaign.csv": tables["co2_o2_estimability_after"],
        "11_selected_global_campaign.csv": tables["global_campaign"],
        "12_selected_co2_campaign.csv": tables["co2_o2_campaign"],
        "13_theta_selected_integrated.csv": safe_read_csv(RESULTS / "theta_selected_integrated.csv", header=None, names=["parameter", "value"]),
        "14_co2_selected_parameters.csv": tables["co2_params"],
    }
    for name, table in mapping.items():
        if isinstance(table, pd.DataFrame) and not table.empty:
            table.to_csv(out / name, index=False)
    target_path = RESULTS / "target_parameters.json"
    if target_path.exists():
        copy_if_exists(target_path, out / "15_target_parameters.json")
    fim_metrics = tables.get("fim_metrics", {})
    if isinstance(fim_metrics, dict):
        (out / "16_fim_metrics_summary.json").write_text(json.dumps(fim_metrics, indent=2), encoding="utf-8")


def display_label(row: pd.Series) -> str:
    state = str(row.get("state", "") or "").strip()
    species = str(row.get("species", "") or "").strip()
    pool = str(row.get("pool", "") or "").strip()
    if state and state != "nan":
        return state
    if species and species != "nan":
        return f"{species}:{pool}" if pool and pool != "nan" else species
    return str(row.name)


def plot_model_selection(df: pd.DataFrame, label_col: str, value_col: str, selected_col: str, title: str, path: Path) -> None:
    if df.empty or value_col not in df.columns:
        return
    data = df.copy()
    data[value_col] = pd.to_numeric(data[value_col], errors="coerce")
    data = data.dropna(subset=[value_col]).sort_values(value_col)
    if data.empty:
        return
    labels = data[label_col].astype(str).tolist()
    colors = ["tab:green" if str(v).lower() == "true" else "tab:blue" for v in data.get(selected_col, pd.Series(False, index=data.index))]
    fig, ax = plt.subplots(figsize=(11, max(4, 0.35 * len(data))))
    ax.barh(labels, data[value_col], color=colors)
    ax.set_xlabel(value_col)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_secondary_selection(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        return
    data = df.copy()
    for col in ("bic", "selection_score"):
        data[col] = pd.to_numeric(data[col], errors="coerce")
    data = data.dropna(subset=["selection_score"]).sort_values("selection_score")
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.barh(data["model"].astype(str), data["selection_score"], color="tab:purple")
    ax.set_xlabel("selection score")
    ax.set_title("Secondary-state structural model selection")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fit_metrics(df: pd.DataFrame, path: Path) -> None:
    if df.empty or "relative_rmse" not in df.columns:
        return
    data = df.copy()
    data["relative_rmse"] = pd.to_numeric(data["relative_rmse"], errors="coerce")
    data = data.dropna(subset=["relative_rmse"])
    data["label"] = data.apply(display_label, axis=1)
    data = data.sort_values(["group", "relative_rmse"])
    colors = data["group"].map({"core": "tab:blue", "secondary": "tab:orange", "aroma": "tab:green", "co2": "tab:red"}).fillna("tab:gray")
    fig, ax = plt.subplots(figsize=(12, max(5, 0.32 * len(data))))
    ax.barh(data["group"].astype(str) + " | " + data["label"], data["relative_rmse"], color=colors)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("relative RMSE")
    ax.set_title("Final model fit quality by observed state/pool")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_estimability(df: pd.DataFrame, path: Path, title: str) -> None:
    if df.empty:
        return
    data = df.copy()
    data["std_log_approx"] = pd.to_numeric(data.get("std_log_approx"), errors="coerce")
    data = data.dropna(subset=["std_log_approx"])
    data["std_plot"] = data["std_log_approx"].clip(lower=1e-4, upper=1e2)
    colors = data["classification"].map(
        {
            "well_estimated": "tab:green",
            "moderate": "tab:blue",
            "weak_but_actionable": "tab:orange",
            "weak_or_confounded": "tab:red",
        }
    ).fillna("tab:gray")
    fig, ax = plt.subplots(figsize=(12, max(5, 0.28 * len(data))))
    ax.barh(data["parameter"].astype(str), data["std_plot"], color=colors)
    ax.set_xscale("log")
    ax.set_xlabel("approx. std(log theta), clipped at 1e2")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_segments(value: str) -> list[float]:
    values: list[float] = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values


def parse_pulses(value: str) -> list[tuple[float, float]]:
    pulses: list[tuple[float, float]] = []
    for item in str(value).split(";"):
        item = item.strip()
        if not item or ":" not in item:
            continue
        t, amount = item.split(":", 1)
        try:
            pulses.append((float(t.lower().replace("h", "").strip()), float(amount.strip())))
        except ValueError:
            continue
    return pulses


def plot_campaign(df: pd.DataFrame, path: Path, title: str) -> None:
    if df.empty:
        return
    data = df.copy().sort_values("campaign_order") if "campaign_order" in df.columns else df.copy()
    n = len(data)
    fig, axes = plt.subplots(n, 1, figsize=(12, max(3, 1.8 * n)), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, (_, row) in zip(axes, data.iterrows()):
        horizon = float(row.get("horizon_h", 300.0) or 300.0)
        seg = parse_segments(row.get("temperature_segments", ""))
        if not seg:
            seg = [20.0]
        edges = np.linspace(0.0, horizon, len(seg) + 1)
        x_vals = []
        y_vals = []
        for i, temp in enumerate(seg):
            x_vals.extend([edges[i], edges[i + 1]])
            y_vals.extend([temp, temp])
        ax.plot(x_vals, y_vals, color="tab:red", linewidth=2)
        for t, amount in parse_pulses(row.get("N_pulses_kg_m3", "")):
            ax.axvline(t, color="tab:blue", linestyle="--", linewidth=1.2)
            ax.text(t, max(seg) + 0.4, f"N {amount:g}", rotation=90, va="bottom", ha="center", fontsize=8, color="tab:blue")
        ax.set_ylabel("T (C)")
        ax.set_ylim(min(seg) - 2, max(seg) + 3)
        ax.grid(alpha=0.25)
        ax.set_title(f"{int(float(row.get('campaign_order', 0) or 0))}. {row.get('candidate', '')}", loc="left", fontsize=9)
    axes[-1].set_xlabel("time (h)")
    fig.suptitle(title, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def create_figures(tables: dict[str, pd.DataFrame | dict | str], figures: Path) -> None:
    copy_if_exists(AROMA_RESULTS / "plots" / "model_selection_comparison.png", figures / "fig_01_aroma_model_selection.png")
    plot_secondary_selection(tables["secondary_selection"], figures / "fig_02_secondary_model_selection.png")  # type: ignore[arg-type]
    plot_model_selection(
        tables["co2_selection"],  # type: ignore[arg-type]
        "mode",
        "selection_score",
        "selected",
        "CO2 model benchmark and selection",
        figures / "fig_03_co2_model_selection.png",
    )
    plot_fit_metrics(tables["final_metrics"], figures / "fig_04_final_fit_relative_rmse.png")  # type: ignore[arg-type]
    plot_estimability(
        tables["global_estimability"],  # type: ignore[arg-type]
        figures / "fig_05_global_estimability_current.png",
        "Current-data practical estimability for global model parameters",
    )
    plot_estimability(
        tables["global_estimability_after"],  # type: ignore[arg-type]
        figures / "fig_06_global_estimability_after_campaign.png",
        "Practical estimability after selected global MBDoE campaign",
    )
    plot_estimability(
        tables["co2_o2_estimability"],  # type: ignore[arg-type]
        figures / "fig_07_co2_estimability_current.png",
        "Current-data estimability for CO2/O2 solubility block",
    )
    plot_estimability(
        tables["co2_o2_estimability_after"],  # type: ignore[arg-type]
        figures / "fig_08_co2_estimability_after_campaign.png",
        "CO2/O2 estimability after selected MBDoE campaign",
    )
    plot_campaign(tables["global_campaign"], figures / "fig_09_global_selected_campaign_inputs.png", "Selected global MBDoE campaign")  # type: ignore[arg-type]
    plot_campaign(tables["co2_o2_campaign"], figures / "fig_10_co2_o2_selected_campaign_inputs.png", "Selected CO2/O2-focused MBDoE campaign")  # type: ignore[arg-type]

    selected_existing = {
        RESULTS / "plots" / "data" / "co2_curated_25170.png": figures / "fig_11_co2_data_curation_25170.png",
        RESULTS / "plots" / "data" / "co2_curated_25171.png": figures / "fig_12_co2_data_curation_25171.png",
        RESULTS / "plots" / "co2_benchmark" / "co2_benchmark_25170.png": figures / "fig_13_co2_benchmark_25170.png",
        RESULTS / "plots" / "co2_benchmark" / "co2_benchmark_25171.png": figures / "fig_14_co2_benchmark_25171.png",
        RESULTS / "plots" / "final_fit" / "final_fit_25170.png": figures / "fig_15_final_fit_representative_25170.png",
        RESULTS / "plots" / "final_fit" / "final_fit_25171.png": figures / "fig_16_final_fit_representative_25171.png",
        RESULTS / "plots" / "fim" / "co2_o2_current_relative_eigenvalues.png": figures / "fig_17_co2_o2_current_relative_eigenvalues.png",
        RESULTS / "plots" / "fim" / "co2_o2_current_plus_campaign_relative_eigenvalues.png": figures / "fig_18_co2_o2_after_campaign_relative_eigenvalues.png",
        GLOBAL_RESULTS / "plots" / "fim" / "global_current_relative_eigenvalues.png": figures / "fig_19_global_current_relative_eigenvalues.png",
        GLOBAL_RESULTS / "plots" / "fim" / "global_current_plus_campaign_relative_eigenvalues.png": figures / "fig_20_global_after_campaign_relative_eigenvalues.png",
    }
    for src, dst in selected_existing.items():
        copy_if_exists(src, dst)


def copy_source_snapshot(source: Path) -> None:
    sources = [
        PILOT_DIR / "run_pilot_2025_co2_solubility_integrated_doe.py",
        PILOT_DIR / "pilot_2025_co2_solubility_integrated_doe.ipynb",
        SCRIPT_DIR / "run_pilot_2025_aroma_model_selection_doe.py",
        SCRIPT_DIR / "run_pilot_2025_global_model_selection_doe.py",
        SCRIPT_DIR / "run_pilot_2025_co2_stripping_benchmark.py",
        SCRIPT_DIR / "pilot_2025_data_loader.py",
        FERMENTATION_MODEL_DIR / "run_new_must_glycerol_estimability_doe.py",
        FERMENTATION_MODEL_DIR / "run_secondary_joint_campaign_doe.py",
    ]
    manifest = []
    for src in sources:
        if src.exists():
            dst = source / src.name
            copy_if_exists(src, dst)
            manifest.append(
                {
                    "bundle_file": f"source_snapshot/{dst.name}",
                    "source_path": str(src),
                    "size_bytes": src.stat().st_size,
                }
            )
    pd.DataFrame(manifest).to_csv(source / "source_manifest.csv", index=False)


def counts_by_class(df: pd.DataFrame) -> str:
    if df.empty or "classification" not in df.columns:
        return "_No disponible._"
    counts = df["classification"].value_counts().rename_axis("classification").reset_index(name="count")
    return md_table(counts)


def best_row(df: pd.DataFrame, selected_col: str = "selected", sort_col: str = "selection_score") -> pd.Series:
    if df.empty:
        return pd.Series(dtype=object)
    if selected_col in df.columns:
        selected = df[df[selected_col].astype(str).str.lower().eq("true")]
        if not selected.empty:
            return selected.iloc[0]
    if sort_col in df.columns:
        data = df.copy()
        data[sort_col] = pd.to_numeric(data[sort_col], errors="coerce")
        data = data.dropna(subset=[sort_col]).sort_values(sort_col)
        if not data.empty:
            return data.iloc[0]
    return df.iloc[0]


def create_report(tables: dict[str, pd.DataFrame | dict | str], root: Path) -> None:
    co2_model = str(tables["co2_model"])
    aroma_model = str(tables["aroma_model"])
    secondary_model = str(tables["secondary_model"])
    co2_row = best_row(tables["co2_selection"])  # type: ignore[arg-type]
    aroma_selection = tables["aroma_selection"] if isinstance(tables["aroma_selection"], pd.DataFrame) else pd.DataFrame()
    secondary_selection = tables["secondary_selection"] if isinstance(tables["secondary_selection"], pd.DataFrame) else pd.DataFrame()
    final_metrics = tables["final_metrics"] if isinstance(tables["final_metrics"], pd.DataFrame) else pd.DataFrame()
    fit_integrated = tables["fit_integrated"] if isinstance(tables["fit_integrated"], pd.DataFrame) else pd.DataFrame()
    co2_fit = tables["co2_fit_selected"] if isinstance(tables["co2_fit_selected"], pd.DataFrame) else pd.DataFrame()
    co2_params = tables["co2_params"] if isinstance(tables["co2_params"], pd.DataFrame) else pd.DataFrame()
    global_estimability = tables["global_estimability"] if isinstance(tables["global_estimability"], pd.DataFrame) else pd.DataFrame()
    global_after = tables["global_estimability_after"] if isinstance(tables["global_estimability_after"], pd.DataFrame) else pd.DataFrame()
    co2_o2_after = tables["co2_o2_estimability_after"] if isinstance(tables["co2_o2_estimability_after"], pd.DataFrame) else pd.DataFrame()
    global_campaign = tables["global_campaign"] if isinstance(tables["global_campaign"], pd.DataFrame) else pd.DataFrame()
    co2_campaign = tables["co2_o2_campaign"] if isinstance(tables["co2_o2_campaign"], pd.DataFrame) else pd.DataFrame()

    integrated_success = ""
    if not fit_integrated.empty and "success" in fit_integrated.columns:
        integrated_success = str(fit_integrated.iloc[0].get("success", ""))

    report = f"""# A07 - Optimización NLP para calibración de síntesis aromática

Actividad 3.14. Bundle generado el {datetime.now().strftime("%Y-%m-%d %H:%M")}.

## 1. Resumen ejecutivo

Se construyó y documentó un pipeline determinista para calibrar y validar un modelo dinámico de fermentación alcohólica con síntesis y pérdida de aromas. El objetivo técnico fue transformar datos piloto de mosto natural, sensores de CO2 y mediciones aromáticas en un modelo ODE calibrable que pueda transferirse posteriormente como restricción dinámica dentro del MPCC.

La estructura vigente seleccionada es:

- Fermentación primaria con biomasa viable/muerta, nitrógeno asimilable, glucosa, fructosa, etanol y glicerol.
- Capa secundaria con piruvato, acetaldehído, acetato y oxígeno latente.
- Aromas: ethyl acetate, isoamyl acetate y ethyl octanoate, separados en retenido líquido y fracción volatilizada/condensada.
- CO2: modelo seleccionado `{co2_model}`, con transición macroscópica de O2 y acumulación/liberación de CO2 disuelto.

Modelos seleccionados:

- Aromas: `{aroma_model}`.
- Secundarios: `{secondary_model}`.
- CO2: `{co2_model}`.

El benchmark de CO2 seleccionó `{co2_model}` con score `{co2_row.get("selection_score", "NA")}` y WSSE `{co2_row.get("data_wsse", "NA")}`. La selección aromática favoreció una estructura donde ethyl acetate depende de etanol, biomasa y limitación de nitrógeno, no solo de una tasa de consumo de azúcar.

## 2. Evidencia fuente incluida

El bundle contiene:

- `A07_technical_report.md`: este informe técnico.
- `tables/`: tablas curadas con selección de modelos, parámetros, métricas de ajuste, estimabilidad y campaña MBDoE.
- `figures/`: figuras seleccionadas para anexar directamente.
- `source_snapshot/`: copia compacta de scripts/notebook fuente que implementan el simulador, calibración, benchmarking, FIM y MBDoE.

No se incluyen todos los CSV o figuras intermedias para evitar un paquete excesivo. Los resultados completos siguen en `fermentation_model/pilot_2025/results/`.

## 3. Datos y validación de lectura

Se procesaron planillas piloto de mosto natural y sensores CO2 online. La curación relevante fue:

- CO2 de 25150 y 25151 se excluye.
- 25170 se usa como set de CO2 más confiable.
- 25171 se reancla a un nuevo tiempo cero en el punto post-reinóculo.
- Aromas `xxx_total` se interpretan como retenido líquido más condensado equivalente.
- Aromas `xxx_condensado` se interpretan como fracción volatilizada acumulada, expresada en concentración equivalente de mosto.

La lectura se verificó gráficamente en `figures/fig_11_*` y `figures/fig_12_*`.

## 4. Modelo dinámico de simulación

### 4.1 Estados primarios

El núcleo de fermentación usa:

\\[
x = [X, X_d, N, G, F, E, Gly]^T
\\]

donde `X` es biomasa viable, `X_d` biomasa muerta, `N` nitrógeno asimilable, `G` glucosa, `F` fructosa, `E` etanol y `Gly` glicerol.

\\[
\\frac{{dX}}{{dt}}=(\\mu-k_d)X+u_X
\\]

\\[
\\frac{{dX_d}}{{dt}}=k_dX
\\]

\\[
\\frac{{dN}}{{dt}}=-q_N f_N(T,N)X+u_N
\\]

\\[
\\frac{{dG}}{{dt}}=-\\left(q_{{XG}}f_N+q_{{EG}}f_G+m\\frac{{G}}{{G+F}}\\right)X+u_G
\\]

\\[
\\frac{{dF}}{{dt}}=-\\left(q_{{XF}}f_N+q_{{EF}}f_F+m\\frac{{F}}{{G+F}}\\right)X+u_F
\\]

\\[
\\frac{{dE}}{{dt}}=(\\beta_G f_G+\\beta_F f_F)X+u_E
\\]

\\[
\\frac{{dGly}}{{dt}}=(\\gamma_G f_G+\\gamma_F f_F)X
\\]

Los factores `f_N`, `f_G` y `f_F` combinan limitación tipo Monod, corrección térmica Arrhenius, inhibición por etanol e inhibición de fructosa por glucosa. Los pulsos operacionales entran como funciones gaussianas suavizadas:

\\[
u_j(t)=\\sum_k \\frac{{a_{{j,k}}}}{{\\sqrt{{\\pi}}w}}\\exp\\left[-\\left(\\frac{{t-t_k}}{{w}}\\right)^2\\right]
\\]

### 4.2 Estados secundarios

La capa secundaria seleccionada (`{secondary_model}`) representa especies que explican desviaciones de aromas y metabolismo lateral:

\\[
z=[Pyr, AcAld, Acetate, O_2]^T
\\]

La estructura usa términos de producción por consumo de azúcar, conversión piruvato-acetaldehído, reducción/oxidación de acetaldehído, generación/asimilación de acetato y efecto latente de O2. El benchmark de estructuras secundarias se resume en `figures/fig_02_secondary_model_selection.png` y `tables/02_secondary_model_selection.csv`.

### 4.3 CO2 disuelto y flujo gaseoso

El CO2 producido se vincula estequiométricamente a producción de etanol:

\\[
r_{{CO2,base}}=\\frac{{44.01}}{{2\\cdot 46.07}}r_E
\\]

En el modelo seleccionado, la producción efectiva se modula por una transición macroscópica de O2:

\\[
r_{{CO2,eff}}=
r_{{CO2,base}}\\left[f_C+(1-f_C)\\phi_{{ana}}(O_2)\\right]+r_{{CO2,resp}}
\\]

\\[
\\phi_{{ana}}(O_2)=\\frac{{K_{{ana}}^n}}{{K_{{ana}}^n+O_2^n}}
\\]

El CO2 disuelto se representa con:

\\[
\\frac{{dC_{{CO2,L}}}}{{dt}}=r_{{CO2,eff}}-r_{{CO2,gas}}
\\]

\\[
r_{{CO2,gas}}=k_{{rel}}\\max(C_{{CO2,L}}-C^*_{{CO2}},0)
\\]

La saturación efectiva:

\\[
C^*_{{CO2}}=
s_{{CO2}}1.69\\exp[-0.032(T-20)]\\exp(0.0016E)\\exp[-0.0012(G+F)]
\\]

Parámetros CO2 seleccionados:

{md_table(co2_params)}

### 4.4 Síntesis y volatilización aromática

Para cada aroma `i`:

\\[
\\frac{{dA_i^{{liq}}}}{{dt}}=r_i^{{prod}}-r_i^{{loss}}
\\]

\\[
\\frac{{dA_i^{{cond}}}}{{dt}}=r_i^{{loss}}
\\]

\\[
A_i^{{total}}=A_i^{{liq}}+A_i^{{cond}}
\\]

La producción base separa fase de crecimiento y fase estacionaria:

\\[
r_i^{{phase}}=
\\left(k_{{i,g}}\\phi_N+k_{{i,s}}(1-\\phi_N)\\right)q_S
\\]

\\[
\\phi_N=\\frac{{N}}{{N+K_N}}
\\]

Para ethyl acetate, el modelo seleccionado agrega una fuente dependiente de etanol/biomasa y limitación de nitrógeno:

\\[
r_{{EA}}=
r_{{EA}}^{{phase}}
+k_{{EA,XE}}X\\frac{{E}}{{K_E+E}}
+k_{{EA,XE,Nlim}}X\\frac{{E}}{{K_E+E}}\\frac{{K_N}}{{K_N+N}}
\\]

La pérdida aromática no se trata como un `k_vol` libre genérico. Se usa:

\\[
r_i^{{loss}}=
\\alpha_{{i,loss}}K_i^{{LG}}(T,E,G,F)q_{{gas}}A_i^{{liq}}
\\]

donde `K_i^{{LG}}` proviene de presión de vapor Antoine y coeficientes de actividad UNIFAC con mezcla agua-etanol-azúcar equivalente.

## 5. Formulación NLP de calibración

La calibración se formula como un problema de mínimos cuadrados no lineal con restricciones de caja:

\\[
\\min_{{\\theta}} \\sum_{{j=1}}^M \\rho\\left(
\\frac{{\\hat y_j(\\theta)-y_j}}{{\\sigma_j}}
\\right)^2
\\]

sujeto a:

\\[
\\dot x=f(x,z,A,C_{{CO2}},u,T,\\theta)
\\]

\\[
\\theta_L \\leq \\theta \\leq \\theta_U
\\]

\\[
x(t)\\geq 0,\\quad z(t)\\geq 0,\\quad A(t)\\geq 0
\\]

La optimización se ejecutó en espacio logarítmico para parámetros positivos:

\\[
\\eta=\\log(\\theta)
\\]

El solver numérico fue `scipy.optimize.least_squares` con:

- método `trf` (trust-region reflective),
- restricciones de caja,
- pérdida robusta `soft_l1`,
- `f_scale=2.0`,
- multistart para modelos estructurales y CO2,
- penalización adicional por parámetros en borde durante selección de modelo.

La selección estructural usó:

\\[
AICc = n\\log(WSSE/n)+2p+\\frac{{2p(p+1)}}{{n-p-1}}
\\]

\\[
BIC = n\\log(WSSE/n)+p\\log(n)
\\]

\\[
Score = BIC + P_{{bounds}} + P_{{adequacy}}
\\]

## 6. Benchmarking de modelos

### 6.1 Aroma

La comparación de modelos de ethyl acetate mostró que la estructura `{aroma_model}` entrega el mejor compromiso entre ajuste y complejidad. La evidencia principal está en `figures/fig_01_aroma_model_selection.png`.

{md_table(aroma_selection, ["model", "n_parameters", "data_wsse", "bic", "active_bound_count", "ea_retained_relative_rmse"], 8)}

### 6.2 Secundarios

{md_table(secondary_selection, ["model", "data_wsse", "n_parameters", "bic", "active_bound_count", "state_penalty", "selection_score"], 8)}

### 6.3 CO2

{md_table(tables["co2_selection"], ["mode", "mechanistic_class", "data_wsse", "bic", "active_bound_count", "selection_penalty", "selection_score", "selected"], 10)}

La estructura `old_lag_threshold` ajusta razonablemente, pero queda penalizada por ser empírica y por tocar borde de parámetro. El modelo `{co2_model}` fue seleccionado porque mejora la interpretación fenomenológica del retraso: transición O2 + acumulación de CO2 disuelto + liberación gaseosa.

## 7. Convergencia y robustez numérica

### CO2

{md_table(co2_fit, ["mode", "start", "success", "status", "nfev", "data_wsse", "active_bound_count"], 10)}

### Ajuste integrado de aromas

{md_table(fit_integrated, ["mode", "start", "success", "status", "nfev", "initial_objective_wsse", "final_objective_wsse", "active_bound_count"], 10)}

Interpretación: el bloque CO2 seleccionado converge satisfactoriamente. El ajuste integrado aromático actual quedó como reconciliación/validación con parámetros heredados de la selección aromática y CO2 fijo; el campo `success={integrated_success}` debe leerse como evidencia de que esta pasada no fue una reoptimización global larga. Para una versión final cerrada del anexo, conviene ejecutar un multistart integrado más largo si se quiere declarar convergencia completa de todos los parámetros aromáticos bajo el CO2 seleccionado.

## 8. Resultados de ajuste

Las métricas finales se reportan como RMSE relativo, sesgo relativo y correlación por estado o pool observado.

{md_table(final_metrics, ["group", "state", "species", "pool", "n", "relative_rmse", "relative_bias", "corr"], 30)}

Figuras principales:

- `figures/fig_04_final_fit_relative_rmse.png`
- `figures/fig_15_final_fit_representative_25170.png`
- `figures/fig_16_final_fit_representative_25171.png`

Resultados generales:

- Estados primarios como etanol y glicerol quedan bien representados.
- Glucosa/fructosa y nitrógeno presentan errores relativos mayores, esperables por heterogeneidad de mosto natural y medición puntual.
- Isoamyl acetate total y retained quedan en rango útil para modelamiento MPCC; condensado mantiene incertidumbre.
- Ethyl acetate sigue siendo el aroma más desafiante, por eso se seleccionó una estructura enriquecida con etanol, biomasa y limitación de N.
- CO2 mejora al incorporar transición O2/solubilidad, especialmente para el inicio efectivo de fermentación.

## 9. Estimabilidad y FIM

La FIM se calculó por diferencias finitas en espacio log-paramétrico:

\\[
J_{{ij}}\\approx
\\frac{{r_i(\\theta_j e^h)-r_i(\\theta_j e^{{-h}})}}{{2h}}
\\]

\\[
F=J^TJ
\\]

La covarianza aproximada:

\\[
\\Sigma_{{\\log\\theta}}\\approx F^{{-1}}
\\]

Clasificación con datos actuales:

{counts_by_class(global_estimability)}

Después de la campaña MBDoE global:

{counts_by_class(global_after)}

Para el bloque CO2/O2, la campaña propuesta transforma ambos parámetros CO2 principales en bien estimados:

{md_table(co2_o2_after, ["parameter", "theta", "std_log_approx", "approx_95_multiplier", "classification"])}

Figuras:

- `figures/fig_05_global_estimability_current.png`
- `figures/fig_06_global_estimability_after_campaign.png`
- `figures/fig_07_co2_estimability_current.png`
- `figures/fig_08_co2_estimability_after_campaign.png`
- `figures/fig_17_*` a `fig_20_*` para espectros de eigenvalores.

## 10. MBDoE

El diseño experimental usa FIM aditiva:

\\[
F_{{total}}=F_{{data}}+\\sum_k F_{{candidate,k}}
\\]

El criterio usado es híbrido:

\\[
\\Phi_{{hybrid}}=
\\log\\det(F)
-2\\left|\\log\\left(\\frac{{\\lambda_{{min}}}}{{\\lambda_{{max}}}}\\right)\\right|
-0.05\\log(\\operatorname{{tr}}(F^{{-1}}))
\\]

Este criterio conserva D-optimalidad como base, pero penaliza direcciones débiles.

Campaña global seleccionada:

{md_table(global_campaign, ["campaign_order", "candidate", "family", "temperature_segments", "N_pulses_kg_m3", "score"], 10)}

Campaña enfocada CO2/O2 seleccionada:

{md_table(co2_campaign, ["campaign_order", "candidate", "family", "temperature_segments", "N_pulses_kg_m3", "score"], 10)}

Las figuras `figures/fig_09_*` y `figures/fig_10_*` muestran los perfiles de temperatura y pulsos de N.

## 11. Conclusiones

1. Se dispone de un simulador ODE documentado y calibrado parcialmente que conecta fermentación primaria, secundarios, CO2 y aromas.
2. La síntesis aromática no queda adecuadamente representada por una cinética de fase simple para ethyl acetate; la estructura seleccionada incorpora dependencia de etanol/biomasa y limitación de N.
3. La volatilización se modela por partición gas-líquido y stripping impulsado por CO2, no como una pérdida empírica aislada.
4. El retraso de CO2 se explica mejor con transición O2 + solubilidad de CO2 que con un lag arbitrario.
5. La FIM muestra parámetros bien estimables y parámetros aún confundidos; el MBDoE propuesto apunta a mejorar direcciones débiles.
6. Para uso MPCC inmediato, conviene transferir como robustos los parámetros bien estimados y usar bounds/priors fuertes en parámetros confounded.

## 12. Pasos siguientes

- Ejecutar un multistart integrado largo para cerrar formalmente convergencia NLP de aromas bajo el CO2 seleccionado.
- Recolectar CO2 online y condensado final en nuevos ensayos diseñados.
- Mantener ethyl acetate, isoamyl acetate y ethyl octanoate como salidas MPCC, pero tratar ethyl acetate con mayor incertidumbre.
- Revisar si el MPCC debe usar todos los parámetros libres o una versión reducida con parámetros poco estimables fijados.
- Usar la campaña MBDoE como base para priorizar nuevos experimentos con mosto natural.
"""
    (root / "A07_technical_report.md").write_text(report, encoding="utf-8")


def create_readme(root: Path) -> None:
    readme = """# Bundle A07 - Optimización NLP para calibración de síntesis aromática

Este paquete resume la evidencia técnica para el anexo A07 de la Actividad 3.14.

## Elementos principales

- `A07_technical_report.md`: informe autoexplicativo con metodología, formulación matemática, resultados y conclusiones.
- `figures/fig_01_*` a `fig_10_*`: figuras de más alto valor para el anexo.
- `tables/01_*` a `tables/14_*`: tablas curadas para trazabilidad de modelo, solver, parámetros, métricas, estimabilidad y MBDoE.

## Elementos auxiliares

- `figures/fig_11_*` a `fig_20_*`: evidencia complementaria de curación CO2, benchmark, ajustes representativos y eigenvalores.
- `source_snapshot/`: scripts/notebook fuente que implementan el pipeline. No están pensados como lectura principal, sino como respaldo reproducible.

## Lectura recomendada

1. Leer `A07_technical_report.md`.
2. Revisar `figures/fig_01_aroma_model_selection.png`, `fig_03_co2_model_selection.png`, `fig_04_final_fit_relative_rmse.png` y `fig_09_global_selected_campaign_inputs.png`.
3. Usar las tablas `03`, `06`, `07`, `11` y `14` para extraer valores concretos.
"""
    (root / "README_bundle.md").write_text(readme, encoding="utf-8")


def create_bundle_manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "relative_path": str(path.relative_to(root)).replace("\\", "/"),
                    "size_bytes": path.stat().st_size,
                }
            )
    pd.DataFrame(rows).to_csv(root / "bundle_manifest.csv", index=False)


def zip_bundle(root: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(root.parent))


def main() -> None:
    figures, tables_dir, source, root = ensure_clean_bundle()
    tables = load_core_tables()
    write_curated_tables(tables, tables_dir)
    create_figures(tables, figures)
    copy_source_snapshot(source)
    create_report(tables, root)
    create_readme(root)
    create_bundle_manifest(root)
    zip_bundle(root, ZIP_PATH)
    print(root)
    print(ZIP_PATH)


if __name__ == "__main__":
    main()
