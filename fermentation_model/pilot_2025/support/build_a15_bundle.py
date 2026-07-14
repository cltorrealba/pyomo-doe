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
REPO_ROOT = FERMENTATION_MODEL_DIR.parent

LAB_DATA = FERMENTATION_MODEL_DIR / "data" / "Laboratorio 2025" / "mosto_sintetico_vl3.xlsx"
PILOT_DATA = FERMENTATION_MODEL_DIR / "data" / "Piloto 2025" / "Calibration_data_vl3.xlsx"

LAB_LOADING = FERMENTATION_MODEL_DIR / "results" / "new_must_data_loading"
LAB_DOE = FERMENTATION_MODEL_DIR / "results" / "new_must_glycerol_estimability_doe"
LAB_OVERNIGHT = FERMENTATION_MODEL_DIR / "results" / "new_must_glycerol_overnight_validation"
MEDIUM_TRANSFER = FERMENTATION_MODEL_DIR / "results" / "medium_transfer_diagnostics"

PILOT_CO2 = PILOT_DIR / "results" / "co2_solubility_integrated_doe"
PILOT_GLOBAL = PILOT_DIR / "results" / "global_state_model_selection_doe"
PILOT_AROMA = PILOT_DIR / "results" / "aroma_model_selection_doe"

BUNDLE_ROOT = PILOT_DIR / "bundles" / "A15_pilot_scalability_2025"
ZIP_PATH = PILOT_DIR / "bundles" / "A15_pilot_scalability_2025.zip"


def read_text(path: Path, default: str = "") -> str:
    if path.exists():
        return path.read_text(encoding="utf-8").strip()
    return default


def safe_csv(path: Path, **kwargs) -> pd.DataFrame:
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
        out[col] = out[col].map(format_cell)
    widths = [max(len(str(col)), *(len(str(v)) for v in out[col].astype(str))) for col in out.columns]
    header = "| " + " | ".join(str(col).ljust(widths[i]) for i, col in enumerate(out.columns)) + " |"
    sep = "| " + " | ".join("-" * widths[i] for i in range(len(widths))) + " |"
    rows = [
        "| " + " | ".join(str(row[col]).ljust(widths[i]) for i, col in enumerate(out.columns)) + " |"
        for _, row in out.iterrows()
    ]
    return "\n".join([header, sep, *rows])


def format_cell(value: object) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.4g}"
    text = str(value)
    return text if len(text) <= 80 else text[:77] + "..."


def ensure_clean_bundle() -> tuple[Path, Path, Path, Path, Path]:
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
    data = BUNDLE_ROOT / "source_data"
    for folder in (figures, tables, source, data):
        folder.mkdir(parents=True, exist_ok=True)
    return figures, tables, source, data, BUNDLE_ROOT


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def first_valid(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.iloc[0]) if not values.empty else float("nan")


def last_valid(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.iloc[-1]) if not values.empty else float("nan")


def count_valid(series: pd.Series) -> int:
    return int(pd.to_numeric(series, errors="coerce").notna().sum())


def workbook_inventory(path: Path, dataset: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    xl = pd.ExcelFile(path)
    rows = []
    for sheet in xl.sheet_names:
        try:
            df = pd.read_excel(path, sheet_name=sheet, nrows=3)
            rows.append({"dataset": dataset, "workbook": path.name, "sheet": sheet, "n_preview_rows": len(df), "n_columns": len(df.columns)})
        except Exception as err:
            rows.append({"dataset": dataset, "workbook": path.name, "sheet": sheet, "n_preview_rows": 0, "n_columns": 0, "error": str(err)})
    return pd.DataFrame(rows)


def pilot_batch_summary(pilot_data: pd.DataFrame, co2_down: pd.DataFrame) -> pd.DataFrame:
    if pilot_data.empty:
        return pd.DataFrame()
    aroma_cols = [col for col in pilot_data.columns if col.endswith("_total")]
    cond_cols = [col for col in pilot_data.columns if col.endswith("_condensate")]
    rows = []
    for batch, group in pilot_data.groupby("batch", sort=True):
        row = {
            "scale": "pilot",
            "medium": "natural",
            "batch": batch,
            "n_rows": len(group),
            "t_min_h": float(pd.to_numeric(group["time_h"], errors="coerce").min()),
            "t_max_h": float(pd.to_numeric(group["time_h"], errors="coerce").max()),
            "temperature_min_c": float(pd.to_numeric(group["temperature_c"], errors="coerce").min()),
            "temperature_max_c": float(pd.to_numeric(group["temperature_c"], errors="coerce").max()),
            "temperature_mean_c": float(pd.to_numeric(group["temperature_c"], errors="coerce").mean()),
            "S_initial_g_l": first_valid(group.get("S_GF_g_l", pd.Series(dtype=float))),
            "S_final_g_l": last_valid(group.get("S_GF_g_l", pd.Series(dtype=float))),
            "YAN_initial_mg_l": first_valid(group.get("YAN_mg_l", pd.Series(dtype=float))),
            "E_initial_g_l": first_valid(group.get("E_g_l", pd.Series(dtype=float))),
            "E_final_g_l": last_valid(group.get("E_g_l", pd.Series(dtype=float))),
            "n_positive_n_pulse_rows": int((pd.to_numeric(group.get("N_pulse_kg_m3", pd.Series(dtype=float)), errors="coerce").fillna(0.0) > 0).sum()),
            "n_aroma_total_obs": int(sum(count_valid(group[col]) for col in aroma_cols)),
            "n_aroma_condensate_obs": int(sum(count_valid(group[col]) for col in cond_cols)),
        }
        if not co2_down.empty and "batch" in co2_down.columns:
            row["n_co2_sensor_points"] = int(co2_down[co2_down["batch"].astype(str).eq(str(batch))].shape[0])
        else:
            row["n_co2_sensor_points"] = 0
        rows.append(row)
    return pd.DataFrame(rows)


def lab_synthetic_summary(batch_summary: pd.DataFrame) -> pd.DataFrame:
    if batch_summary.empty:
        return pd.DataFrame()
    out = batch_summary[batch_summary["medium"].astype(str).eq("synthetic")].copy()
    out["scale"] = "laboratory"
    out["S_initial_g_l"] = pd.to_numeric(out.get("G_initial_g_l"), errors="coerce") + pd.to_numeric(out.get("F_initial_g_l"), errors="coerce")
    return out


def aggregate_lab_residuals(residuals: pd.DataFrame) -> pd.DataFrame:
    if residuals.empty:
        return pd.DataFrame()
    data = residuals.copy()
    for col in ("rmse", "weighted_rmse", "wsse", "n"):
        data[col] = pd.to_numeric(data[col], errors="coerce")
    return (
        data.groupby(["fit", "medium", "state"], as_index=False)
        .agg(
            n_batches=("batch", "nunique"),
            n_obs=("n", "sum"),
            mean_rmse=("rmse", "mean"),
            median_rmse=("rmse", "median"),
            mean_weighted_rmse=("weighted_rmse", "mean"),
            total_wsse=("wsse", "sum"),
        )
        .sort_values(["fit", "medium", "state"])
    )


def load_tables() -> dict[str, object]:
    pilot_data = safe_csv(PILOT_CO2 / "pilot_model_data_used.csv")
    co2_down = safe_csv(PILOT_CO2 / "co2_downsampled.csv")
    lab_batch = safe_csv(LAB_LOADING / "new_must_batch_summary.csv")
    lab_resid = safe_csv(LAB_DOE / "residual_by_state.csv")
    tables: dict[str, object] = {
        "workbook_inventory": pd.concat(
            [workbook_inventory(LAB_DATA, "laboratory_synthetic_vl3"), workbook_inventory(PILOT_DATA, "pilot_2025")],
            ignore_index=True,
        ),
        "lab_batch_summary": lab_synthetic_summary(lab_batch),
        "all_lab_batch_summary": lab_batch,
        "pilot_batch_summary": pilot_batch_summary(pilot_data, co2_down),
        "pilot_model_data": pilot_data,
        "density_fit": safe_csv(LAB_LOADING / "density_sugar_fit_by_medium.csv"),
        "lab_fit_summary": safe_csv(LAB_DOE / "fit_summary.csv"),
        "lab_residual_by_state": lab_resid,
        "lab_residual_agg": aggregate_lab_residuals(lab_resid),
        "lab_fim_summary": safe_csv(LAB_DOE / "fim_summary.csv"),
        "lab_estimability": safe_csv(LAB_DOE / "parameter_estimability_summary.csv"),
        "lab_profile_summary": safe_csv(LAB_DOE / "profile_summary_combined.csv"),
        "lab_campaign": safe_csv(LAB_DOE / "selected_campaign_hybrid.csv"),
        "lab_candidate_ranking": safe_csv(LAB_DOE / "candidate_ranking.csv"),
        "pilot_final_metrics": safe_csv(PILOT_CO2 / "final_fit_metrics.csv"),
        "pilot_co2_selection": safe_csv(PILOT_CO2 / "co2_model_selection_summary.csv"),
        "pilot_co2_params": safe_csv(PILOT_CO2 / "co2_params_selected.csv", header=None, names=["parameter", "value"]),
        "pilot_co2_curation": safe_csv(PILOT_CO2 / "co2_curation_decisions.csv"),
        "pilot_aroma_selection": safe_csv(PILOT_AROMA / "model_selection_summary.csv"),
        "pilot_secondary_selection": safe_csv(PILOT_GLOBAL / "secondary_model_selection_summary.csv"),
        "pilot_estimability": safe_csv(PILOT_GLOBAL / "parameter_estimability_global_current.csv"),
        "pilot_estimability_after": safe_csv(PILOT_GLOBAL / "parameter_estimability_global_current_plus_campaign.csv"),
        "pilot_campaign": safe_csv(PILOT_GLOBAL / "selected_campaign_hybrid_global_selected.csv"),
        "pilot_co2_campaign": safe_csv(PILOT_CO2 / "co2_o2_selected_campaign_hybrid.csv"),
        "medium_transfer_audit": safe_csv(MEDIUM_TRANSFER / "historical_batch_audit.csv"),
        "lab_flags": read_lab_flags(),
        "pilot_incidents": pilot_incidents(),
        "selected_co2_model": read_text(PILOT_CO2 / "selected_co2_model.txt", "not_available"),
        "selected_aroma_model": read_text(PILOT_CO2 / "selected_aroma_model_inherited.txt", "not_available"),
        "selected_secondary_model": read_text(PILOT_CO2 / "selected_secondary_model_inherited.txt", "not_available"),
    }
    return tables


def read_lab_flags() -> pd.DataFrame:
    if not LAB_DATA.exists():
        return pd.DataFrame()
    try:
        return pd.read_excel(LAB_DATA, sheet_name="Flags_calidad")
    except Exception:
        return pd.DataFrame()


def pilot_incidents() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "item": "CO2 25150/25151",
                "type": "sensor_data_exclusion",
                "description": "Los archivos CO2 de 25150 y 25151 fueron marcados como no usables para calibración.",
                "model_action": "Excluir del ajuste CO2; mantener datos offline de composición.",
            },
            {
                "item": "25171",
                "type": "process_incident",
                "description": "Inóculo inicial con baja/no viabilidad; el proceso útil se reancla al punto post reinóculo.",
                "model_action": "Usar el punto ING25-SB012-Pre reinóculo como t=0 efectivo.",
            },
            {
                "item": "Configuración física piloto",
                "type": "metadata_gap",
                "description": "El workbook no codifica volumen de estanque, geometría, agitación ni detalle de equipamiento.",
                "model_action": "Reportar como brecha documental; el análisis usa variables operacionales medidas.",
            },
        ]
    )


def write_tables(tables: dict[str, object], out: Path) -> None:
    mapping = {
        "01_workbook_inventory.csv": tables["workbook_inventory"],
        "02_laboratory_synthetic_batch_summary.csv": tables["lab_batch_summary"],
        "03_pilot_batch_summary.csv": tables["pilot_batch_summary"],
        "04_density_sugar_fit.csv": tables["density_fit"],
        "05_lab_calibration_fit_summary.csv": tables["lab_fit_summary"],
        "06_lab_residuals_by_state_aggregated.csv": tables["lab_residual_agg"],
        "07_lab_fim_summary.csv": tables["lab_fim_summary"],
        "08_lab_parameter_estimability.csv": tables["lab_estimability"],
        "09_lab_profile_summary.csv": tables["lab_profile_summary"],
        "10_lab_selected_campaign_hybrid.csv": tables["lab_campaign"],
        "11_pilot_final_fit_metrics.csv": tables["pilot_final_metrics"],
        "12_pilot_co2_model_selection.csv": tables["pilot_co2_selection"],
        "13_pilot_secondary_model_selection.csv": tables["pilot_secondary_selection"],
        "14_pilot_aroma_model_selection.csv": tables["pilot_aroma_selection"],
        "15_pilot_parameter_estimability.csv": tables["pilot_estimability"],
        "16_pilot_selected_campaign_global.csv": tables["pilot_campaign"],
        "17_pilot_selected_campaign_co2.csv": tables["pilot_co2_campaign"],
        "18_lab_quality_flags.csv": tables["lab_flags"],
        "19_pilot_incidents_and_actions.csv": tables["pilot_incidents"],
        "20_pilot_co2_curation_decisions.csv": tables["pilot_co2_curation"],
        "21_pilot_co2_selected_parameters.csv": tables["pilot_co2_params"],
    }
    for name, table in mapping.items():
        if isinstance(table, pd.DataFrame) and not table.empty:
            table.to_csv(out / name, index=False)
    transfer_matrix(tables).to_csv(out / "22_transferability_decision_matrix.csv", index=False)


def transfer_matrix(tables: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "layer": "primary_kinetics",
                "lab_synthetic_evidence": "CCD synthetic VL3 excites sugar, temperature and nitrogen directions.",
                "pilot_evidence": "Pilot natural must reproduces primary trends but with matrix and scale deviations.",
                "transferability": "partial_to_good",
                "action": "Use shared structure; keep recalibrated pilot priors for MPCC.",
            },
            {
                "layer": "glycerol",
                "lab_synthetic_evidence": "Glycerol state added and fitted in mixed lab dataset.",
                "pilot_evidence": "Pilot final relative RMSE for Gly is low compared with other primary states.",
                "transferability": "good",
                "action": "Retain glycerol state in scale-up model.",
            },
            {
                "layer": "nitrogen",
                "lab_synthetic_evidence": "Synthetic low/high YAN informs qN and sN but sN remains weak.",
                "pilot_evidence": "YAN residuals remain high and natural must has remanent YAN behavior.",
                "transferability": "limited",
                "action": "Use stronger priors; consider ammonium/PAN split in future.",
            },
            {
                "layer": "CO2",
                "lab_synthetic_evidence": "Not present in synthetic VL3 model stage.",
                "pilot_evidence": "Online pilot CO2 required dissolved CO2 and O2-transition structure.",
                "transferability": "requires_pilot_extension",
                "action": "Use pilot-specific CO2 model for MPCC constraints.",
            },
            {
                "layer": "aromas",
                "lab_synthetic_evidence": "Not directly calibrated in synthetic VL3 stage.",
                "pilot_evidence": "Pilot aroma/condensate data enabled EA/IAA/EO synthesis-loss layer.",
                "transferability": "requires_pilot_extension",
                "action": "Use pilot aroma structure; keep ethyl acetate uncertainty explicit.",
            },
            {
                "layer": "DOE",
                "lab_synthetic_evidence": "Synthetic designs ranked high for separating sugar and weak kinetic directions.",
                "pilot_evidence": "Pilot designs favor natural-must temperature and N perturbations.",
                "transferability": "complementary",
                "action": "Use synthetic for mechanism separation; pilot for scale/matrix validation.",
            },
        ]
    )


def create_figures(tables: dict[str, object], figures: Path) -> None:
    plot_dataset_overview(tables, figures / "fig_01_dataset_coverage_and_operational_space.png")
    plot_lab_vs_pilot_kpis(tables, figures / "fig_02_lab_vs_pilot_initial_conditions_kpis.png")
    copy_if_exists(LAB_LOADING / "density_sugar_regression.png", figures / "fig_03_density_sugar_regression_lab.png")
    copy_if_exists(LAB_LOADING / "operational_inputs_by_medium.png", figures / "fig_04_lab_operational_inputs_by_medium.png")
    plot_lab_fit_metrics(tables, figures / "fig_05_lab_transfer_fit_metrics.png")
    plot_pilot_fit_metrics(tables, figures / "fig_06_pilot_fit_metrics.png")
    plot_fim_summary(tables, figures / "fig_07_lab_fim_transfer_summary.png")
    plot_estimability_counts(tables, figures / "fig_08_estimability_counts_lab_and_pilot.png")
    plot_campaign(tables["lab_campaign"], figures / "fig_09_lab_mbdoe_selected_campaign.png", "Laboratory/synthetic MBDoE selected campaign")
    plot_campaign(tables["pilot_campaign"], figures / "fig_10_pilot_mbdoe_selected_campaign.png", "Pilot natural-must MBDoE selected campaign")
    plot_model_selection(tables["pilot_co2_selection"], "mode", "selection_score", "selected", "Pilot CO2 model selection", figures / "fig_11_pilot_co2_model_selection.png")
    plot_model_selection(tables["pilot_aroma_selection"], "model", "bic", None, "Pilot aroma model selection by BIC", figures / "fig_12_pilot_aroma_model_selection.png")
    copy_if_exists(LAB_DOE / "plots" / "fit_mixed_full17_l2_synthetic_MS007.png", figures / "fig_13_lab_synthetic_fit_MS007.png")
    copy_if_exists(LAB_DOE / "plots" / "fit_mixed_full17_l2_synthetic_MS009.png", figures / "fig_14_lab_synthetic_fit_MS009.png")
    copy_if_exists(PILOT_CO2 / "plots" / "final_fit" / "final_fit_25170.png", figures / "fig_15_pilot_final_fit_25170.png")
    copy_if_exists(PILOT_CO2 / "plots" / "final_fit" / "final_fit_25171.png", figures / "fig_16_pilot_final_fit_25171.png")
    copy_if_exists(PILOT_CO2 / "plots" / "co2_benchmark" / "co2_benchmark_25170.png", figures / "fig_17_pilot_co2_benchmark_25170.png")
    copy_if_exists(PILOT_CO2 / "plots" / "co2_benchmark" / "co2_benchmark_25171.png", figures / "fig_18_pilot_co2_benchmark_25171.png")
    copy_if_exists(LAB_DOE / "plots" / "fim_eigen_spectra.png", figures / "fig_19_lab_fim_eigen_spectra.png")
    copy_if_exists(PILOT_GLOBAL / "plots" / "fim" / "global_current_relative_eigenvalues.png", figures / "fig_20_pilot_global_current_relative_eigenvalues.png")


def plot_dataset_overview(tables: dict[str, object], path: Path) -> None:
    lab = tables["lab_batch_summary"]
    pilot = tables["pilot_batch_summary"]
    if not isinstance(lab, pd.DataFrame) or not isinstance(pilot, pd.DataFrame) or lab.empty or pilot.empty:
        return
    rows = [
        {"dataset": "Lab synthetic VL3", "batches": lab["batch"].nunique(), "rows": len(lab), "max_h": pd.to_numeric(lab["t_max_h"], errors="coerce").max()},
        {"dataset": "Pilot 2025", "batches": pilot["batch"].nunique(), "rows": int(pd.to_numeric(pilot["n_rows"], errors="coerce").sum()), "max_h": pd.to_numeric(pilot["t_max_h"], errors="coerce").max()},
    ]
    data = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, col, title in zip(axes, ["batches", "rows", "max_h"], ["Batches", "Rows/samples", "Max horizon (h)"]):
        ax.bar(data["dataset"], data[col], color=["tab:blue", "tab:green"])
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=20)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Dataset coverage used for transferability/scalability analysis")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_lab_vs_pilot_kpis(tables: dict[str, object], path: Path) -> None:
    lab = tables["lab_batch_summary"]
    pilot = tables["pilot_batch_summary"]
    if not isinstance(lab, pd.DataFrame) or not isinstance(pilot, pd.DataFrame) or lab.empty or pilot.empty:
        return
    frames = []
    frames.append(
        lab.assign(dataset="Lab synthetic VL3").rename(
            columns={"temperature_initial_c": "temperature_initial_c", "S_initial_g_l": "S_initial_g_l"}
        )[["dataset", "batch", "S_initial_g_l", "temperature_initial_c", "t_max_h"]]
    )
    frames.append(
        pilot.assign(dataset="Pilot 2025").rename(columns={"temperature_mean_c": "temperature_initial_c"})[
            ["dataset", "batch", "S_initial_g_l", "temperature_initial_c", "t_max_h"]
        ]
    )
    data = pd.concat(frames, ignore_index=True)
    metrics = [("S_initial_g_l", "Initial G+F (g/L)"), ("temperature_initial_c", "Temperature KPI (C)"), ("t_max_h", "Process horizon (h)")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (col, title) in zip(axes, metrics):
        groups = [pd.to_numeric(data[data["dataset"].eq(ds)][col], errors="coerce").dropna() for ds in ["Lab synthetic VL3", "Pilot 2025"]]
        ax.boxplot(groups, tick_labels=["Lab synthetic", "Pilot"], patch_artist=True)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Operational-space comparison")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_lab_fit_metrics(tables: dict[str, object], path: Path) -> None:
    agg = tables["lab_residual_agg"]
    if not isinstance(agg, pd.DataFrame) or agg.empty:
        return
    data = agg[agg["fit"].astype(str).eq("mixed_full17_l2")].copy()
    if data.empty:
        return
    data["mean_weighted_rmse"] = pd.to_numeric(data["mean_weighted_rmse"], errors="coerce")
    pivot = data.pivot_table(index="state", columns="medium", values="mean_weighted_rmse", aggfunc="mean")
    pivot = pivot.reindex(["X", "Xd", "N", "G", "F", "E", "Gly"])
    fig, ax = plt.subplots(figsize=(10, 4.8))
    pivot.plot(kind="bar", ax=ax)
    ax.set_ylabel("mean weighted RMSE")
    ax.set_title("Laboratory mixed model fit by state and medium")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def label_pilot_metric(row: pd.Series) -> str:
    state = str(row.get("state", "") or "").strip()
    species = str(row.get("species", "") or "").strip()
    pool = str(row.get("pool", "") or "").strip()
    if state and state != "nan":
        return state
    if species and species != "nan":
        return f"{species}:{pool}" if pool and pool != "nan" else species
    return "aroma"


def plot_pilot_fit_metrics(tables: dict[str, object], path: Path) -> None:
    metrics = tables["pilot_final_metrics"]
    if not isinstance(metrics, pd.DataFrame) or metrics.empty:
        return
    data = metrics.copy()
    data["relative_rmse"] = pd.to_numeric(data["relative_rmse"], errors="coerce")
    data = data.dropna(subset=["relative_rmse"])
    data["label"] = data.apply(label_pilot_metric, axis=1)
    data = data.sort_values(["group", "relative_rmse"])
    colors = data["group"].map({"core": "tab:blue", "secondary": "tab:orange", "aroma": "tab:green", "co2": "tab:red"}).fillna("tab:gray")
    fig, ax = plt.subplots(figsize=(11, max(4, 0.3 * len(data))))
    ax.barh(data["group"].astype(str) + " | " + data["label"], data["relative_rmse"], color=colors)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("relative RMSE")
    ax.set_title("Pilot 2025 calibrated model fit quality")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fim_summary(tables: dict[str, object], path: Path) -> None:
    fim = tables["lab_fim_summary"]
    if not isinstance(fim, pd.DataFrame) or fim.empty:
        return
    data = fim.copy()
    for col in ("logdet", "condition_number", "trace_inv"):
        data[col] = pd.to_numeric(data[col], errors="coerce")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].bar(data["analysis"], data["logdet"], color="tab:blue")
    axes[0].set_title("FIM logdet")
    axes[1].bar(data["analysis"], data["condition_number"], color="tab:orange")
    axes[1].set_yscale("log")
    axes[1].set_title("Condition number")
    axes[2].bar(data["analysis"], data["trace_inv"], color="tab:green")
    axes[2].set_yscale("log")
    axes[2].set_title("trace(inv(F))")
    for ax in axes:
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Laboratory natural/synthetic/mixed information diagnostics")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_estimability_counts(tables: dict[str, object], path: Path) -> None:
    lab = tables["lab_estimability"]
    pilot = tables["pilot_estimability"]
    frames = []
    if isinstance(lab, pd.DataFrame) and not lab.empty:
        subset = lab[lab["analysis"].astype(str).eq("mixed_current_full_l2")].copy()
        if not subset.empty:
            frames.append(subset.assign(dataset="Lab mixed full L2")[["dataset", "classification"]])
    if isinstance(pilot, pd.DataFrame) and not pilot.empty:
        frames.append(pilot.assign(dataset="Pilot global")[["dataset", "classification"]])
    if not frames:
        return
    data = pd.concat(frames, ignore_index=True)
    counts = data.value_counts(["dataset", "classification"]).reset_index(name="count")
    pivot = counts.pivot(index="dataset", columns="classification", values="count").fillna(0)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("parameter count")
    ax.set_title("Practical estimability classification")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_model_selection(obj: object, label_col: str, value_col: str, selected_col: str | None, title: str, path: Path) -> None:
    if not isinstance(obj, pd.DataFrame) or obj.empty or value_col not in obj.columns:
        return
    data = obj.copy()
    data[value_col] = pd.to_numeric(data[value_col], errors="coerce")
    data = data.dropna(subset=[value_col]).sort_values(value_col)
    if data.empty or label_col not in data.columns:
        return
    colors = "tab:blue"
    if selected_col and selected_col in data.columns:
        colors = ["tab:green" if str(v).lower() == "true" else "tab:blue" for v in data[selected_col]]
    fig, ax = plt.subplots(figsize=(11, max(4, 0.35 * len(data))))
    ax.barh(data[label_col].astype(str), data[value_col], color=colors)
    ax.set_xlabel(value_col)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_segments(value: object) -> list[float]:
    values = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values


def parse_pulses(value: object) -> list[tuple[float, float]]:
    pulses = []
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


def plot_campaign(obj: object, path: Path, title: str) -> None:
    if not isinstance(obj, pd.DataFrame) or obj.empty:
        return
    data = obj.copy()
    if "campaign_order" in data.columns:
        data["campaign_order_num"] = pd.to_numeric(data["campaign_order"], errors="coerce")
        data = data.sort_values("campaign_order_num")
    data = data.head(9)
    n = len(data)
    fig, axes = plt.subplots(n, 1, figsize=(12, max(3, 1.65 * n)), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, (_, row) in zip(axes, data.iterrows()):
        horizon = float(row.get("horizon_h", 216.0) or 216.0)
        seg = parse_segments(row.get("temperature_segments", ""))
        if not seg:
            seg = [20.0]
        edges = np.linspace(0.0, horizon, len(seg) + 1)
        xs: list[float] = []
        ys: list[float] = []
        for idx, temp in enumerate(seg):
            xs.extend([edges[idx], edges[idx + 1]])
            ys.extend([temp, temp])
        ax.plot(xs, ys, color="tab:red", linewidth=2)
        for t, amount in parse_pulses(row.get("N_pulses_kg_m3", "")):
            ax.axvline(t, color="tab:blue", linestyle="--", linewidth=1.2)
            ax.text(t, max(seg) + 0.4, f"N {amount:g}", rotation=90, va="bottom", ha="center", fontsize=8, color="tab:blue")
        ax.set_ylabel("T (C)")
        ax.set_ylim(min(seg) - 2, max(seg) + 3)
        ax.grid(alpha=0.25)
        order = row.get("campaign_order", "")
        ax.set_title(f"{order}. {row.get('candidate', '')}", loc="left", fontsize=9)
    axes[-1].set_xlabel("time (h)")
    fig.suptitle(title, y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_source_snapshot(source: Path, data_dir: Path) -> None:
    source_files = [
        SCRIPT_DIR / "build_a15_bundle.py",
        PILOT_DIR / "run_pilot_2025_co2_solubility_integrated_doe.py",
        SCRIPT_DIR / "pilot_2025_data_loader.py",
        SCRIPT_DIR / "run_pilot_2025_calibration_estimability.py",
        SCRIPT_DIR / "run_pilot_2025_global_model_selection_doe.py",
        SCRIPT_DIR / "run_pilot_2025_aroma_model_selection_doe.py",
        FERMENTATION_MODEL_DIR / "run_new_must_data_loading.py",
        FERMENTATION_MODEL_DIR / "run_new_must_glycerol_estimability_doe.py",
        FERMENTATION_MODEL_DIR / "run_secondary_joint_campaign_doe.py",
    ]
    manifest = []
    for src in source_files:
        if src.exists():
            dst = source / src.name
            copy_if_exists(src, dst)
            manifest.append({"bundle_file": f"source_snapshot/{dst.name}", "source_path": str(src), "size_bytes": src.stat().st_size})
    for src in (LAB_DATA, PILOT_DATA):
        if src.exists():
            dst = data_dir / src.name
            copy_if_exists(src, dst)
            manifest.append({"bundle_file": f"source_data/{dst.name}", "source_path": str(src), "size_bytes": src.stat().st_size})
    pd.DataFrame(manifest).to_csv(source / "source_manifest.csv", index=False)


def class_counts(df: pd.DataFrame, analysis: str | None = None) -> str:
    if df.empty or "classification" not in df.columns:
        return "_No disponible._"
    data = df.copy()
    if analysis and "analysis" in data.columns:
        data = data[data["analysis"].astype(str).eq(analysis)]
    if data.empty:
        return "_No disponible._"
    counts = data["classification"].value_counts().rename_axis("classification").reset_index(name="count")
    return md_table(counts)


def create_report(tables: dict[str, object], root: Path) -> None:
    lab_summary = tables["lab_batch_summary"] if isinstance(tables["lab_batch_summary"], pd.DataFrame) else pd.DataFrame()
    pilot_summary = tables["pilot_batch_summary"] if isinstance(tables["pilot_batch_summary"], pd.DataFrame) else pd.DataFrame()
    density = tables["density_fit"] if isinstance(tables["density_fit"], pd.DataFrame) else pd.DataFrame()
    lab_fit = tables["lab_fit_summary"] if isinstance(tables["lab_fit_summary"], pd.DataFrame) else pd.DataFrame()
    lab_fim = tables["lab_fim_summary"] if isinstance(tables["lab_fim_summary"], pd.DataFrame) else pd.DataFrame()
    lab_est = tables["lab_estimability"] if isinstance(tables["lab_estimability"], pd.DataFrame) else pd.DataFrame()
    pilot_metrics = tables["pilot_final_metrics"] if isinstance(tables["pilot_final_metrics"], pd.DataFrame) else pd.DataFrame()
    pilot_co2 = tables["pilot_co2_selection"] if isinstance(tables["pilot_co2_selection"], pd.DataFrame) else pd.DataFrame()
    pilot_sec = tables["pilot_secondary_selection"] if isinstance(tables["pilot_secondary_selection"], pd.DataFrame) else pd.DataFrame()
    pilot_aroma = tables["pilot_aroma_selection"] if isinstance(tables["pilot_aroma_selection"], pd.DataFrame) else pd.DataFrame()
    pilot_est = tables["pilot_estimability"] if isinstance(tables["pilot_estimability"], pd.DataFrame) else pd.DataFrame()
    lab_campaign = tables["lab_campaign"] if isinstance(tables["lab_campaign"], pd.DataFrame) else pd.DataFrame()
    pilot_campaign = tables["pilot_campaign"] if isinstance(tables["pilot_campaign"], pd.DataFrame) else pd.DataFrame()
    incidents = tables["pilot_incidents"] if isinstance(tables["pilot_incidents"], pd.DataFrame) else pd.DataFrame()
    flags = tables["lab_flags"] if isinstance(tables["lab_flags"], pd.DataFrame) else pd.DataFrame()
    transfer = transfer_matrix(tables)

    selected_co2 = str(tables.get("selected_co2_model", "not_available"))
    selected_aroma = str(tables.get("selected_aroma_model", "not_available"))
    selected_secondary = str(tables.get("selected_secondary_model", "not_available"))

    report = f"""# A15 - Validación de escalabilidad en bodega piloto vendimia 2025

Actividad 7.18. Bundle generado el {datetime.now().strftime("%Y-%m-%d %H:%M")}.

## 1. Resumen ejecutivo

Este anexo documenta la validación de escalabilidad desde ensayos de laboratorio con mosto sintético VL3 hacia vinificaciones piloto de vendimia 2025. El foco no es solo mostrar curvas de ajuste, sino evaluar si la estructura cinética/metabólica desarrollada en laboratorio puede transferirse a mosto natural y escala piloto, y qué extensiones fueron necesarias.

Resultado principal:

- La estructura primaria de fermentación es transferible de forma parcial a buena: biomasa, azúcar, nitrógeno, etanol y glicerol comparten una arquitectura ODE común.
- La transferencia directa desde mosto sintético no es suficiente para describir todo el piloto: se requirieron extensiones para CO2, transición de oxígeno, secundarios y aromas.
- Los datos piloto 2025 permitieron validar escalabilidad operacional y revelar brechas reales: CO2 online, condensado aromático, reinóculo 25171, y exclusión de sensores CO2 25150/25151.
- El enfoque recomendado es mixto: laboratorio sintético para separar mecanismos e identificabilidad; piloto natural para validar matriz, escala y outputs aromáticos.

Modelos piloto seleccionados:

- CO2: `{selected_co2}`.
- Secundarios: `{selected_secondary}`.
- Aromas: `{selected_aroma}`.

## 2. Evidencia fuente incluida

El bundle contiene:

- `A15_technical_report.md`: informe autoexplicativo.
- `README_bundle.md`: guía del contenido del paquete.
- `tables/`: tablas curadas de batches, KPIs, selección de modelos, FIM, estimabilidad, incidencias y campañas.
- `figures/`: figuras seleccionadas para anexar.
- `source_data/`: workbooks fuente compactos usados como evidencia (`mosto_sintetico_vl3.xlsx` y `Calibration_data_vl3.xlsx`).
- `source_snapshot/`: scripts fuente mínimos para trazabilidad.

## 3. Configuración experimental y datasets

### 3.1 Laboratorio sintético VL3

Fuente: `fermentation_model/data/Laboratorio 2025/mosto_sintetico_vl3.xlsx`.

La planilla incluye 10 cubadas sintéticas MS007-MS016, hojas de resumen, datos homologados, diseño CCD, flags de calidad y diccionario de mapeo. El diseño VL3 cubre principalmente variación de temperatura y composición inicial de azúcares/nutrientes, con mosto sintético como matriz controlada.

Resumen de batches sintéticos:

{md_table(lab_summary, ["batch", "n_rows", "t_max_h", "temperature_initial_c", "temperature_min_c", "temperature_max_c", "G_initial_g_l", "F_initial_g_l", "S_initial_g_l"], 12)}

### 3.2 Piloto vendimia 2025

Fuente: `fermentation_model/data/Piloto 2025/Calibration_data_vl3.xlsx`.

La planilla piloto contiene 8 vinificaciones naturales: 25026, 25027, 25085, 25086, 25150, 25151, 25170 y 25171. Incluye temperatura, densidad, biomasa viable, peso seco, glucosa, fructosa, PAN, amonio, YAN, glicerol, piruvato, acetaldehído, etanol, pulsos de nutriente, aromas totales y aromas en condensado.

Resumen de batches piloto usados por el pipeline:

{md_table(pilot_summary, ["batch", "n_rows", "t_max_h", "temperature_min_c", "temperature_max_c", "S_initial_g_l", "S_final_g_l", "YAN_initial_mg_l", "E_final_g_l", "n_aroma_total_obs", "n_co2_sensor_points"], 12)}

### 3.3 Brechas de configuración física

El workbook piloto no codifica completamente la configuración física de bodega piloto: volumen útil, geometría, agitación, transferencia térmica ni hardware de sensor. Por tanto, el análisis de escalabilidad se basa en variables operacionales medidas y en evidencia de proceso, no en un balance físico completo de estanque.

Esto debe completarse manualmente en la bitácora oficial del proyecto si se requiere trazabilidad de hardware.

## 4. Validación de datos y homologación

Se homologaron variables entre laboratorio y piloto:

- Tiempo de fermentación.
- Temperatura operacional.
- Glucosa y fructosa.
- Azúcar total a partir de densidad.
- Biomasa viable.
- Nitrógeno asimilable.
- Etanol.
- Glicerol.
- Pulsos de nutriente.
- Estados secundarios y aromas en piloto.

La regresión densidad-azúcar mostró alta consistencia dentro de cada medio:

{md_table(density, ["group", "source", "intercept", "slope", "r2", "rmse_g_l", "n_points"])}

Figura clave: `figures/fig_03_density_sugar_regression_lab.png`.

Flags/incidencias de laboratorio:

{md_table(flags)}

Incidencias piloto:

{md_table(incidents)}

## 5. Modelo base transferido desde laboratorio

El modelo primario transferido usa:

\\[
x=[X,X_d,N,G,F,E,Gly]^T
\\]

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

Los factores cinéticos incluyen dependencia de temperatura, saturación por sustrato, inhibición por etanol e interacción glucosa/fructosa. Esta estructura fue primero evaluada en laboratorio y después reutilizada como base para piloto.

## 6. Calibración laboratorio y evaluación medio sintético/natural

La calibración de laboratorio comparó fits por medio y un fit mixto:

{md_table(lab_fit, ["fit", "n_batches", "mediums", "n_parameters", "success", "nfev", "final_wsse", "wsse_per_residual", "l2_lambda"], 10)}

Interpretación:

- El ajuste sintético reducido tuvo menor WSSE/residual que el natural reducido.
- El ajuste mixto completo con L2 permitió una estructura común con más parámetros sin perder estabilidad.
- La matriz natural introduce desviaciones que no aparecen en sintético puro.
- Por tanto, el mosto sintético es útil para identificación mecanística, pero no debe ser el único prior para piloto.

Figuras:

- `figures/fig_05_lab_transfer_fit_metrics.png`
- `figures/fig_13_lab_synthetic_fit_MS007.png`
- `figures/fig_14_lab_synthetic_fit_MS009.png`

## 7. FIM e identificabilidad en laboratorio

La FIM se calculó en espacio log-paramétrico:

\\[
J_{{ij}}\\approx\\frac{{r_i(\\theta_j e^h)-r_i(\\theta_j e^{{-h}})}}{{2h}}
\\]

\\[
F=J^TJ
\\]

Resumen FIM:

{md_table(lab_fim, ["analysis", "n_residuals", "logdet", "min_eigenvalue", "condition_number", "trace_inv", "rank_1e-8"])}

Estimabilidad del modelo mixto laboratorio:

{class_counts(lab_est, "mixed_current_full_l2")}

Conclusión de laboratorio:

- La mezcla sintético/natural mejora el contenido de información frente a cada medio por separado.
- Persisten direcciones débiles: principalmente saturación de nitrógeno, mantenimiento y algunos términos acoplados a fructosa.
- La campaña MBDoE de laboratorio prioriza diseños sintéticos porque permiten perturbar composición inicial y separar mecanismos que en mosto natural están correlacionados.

## 8. Validación piloto y extensiones requeridas

Al transferir el modelo a piloto, se mantuvo el núcleo primario pero se extendió el simulador:

1. Capa secundaria: piruvato, acetaldehído, acetato y O2 latente.
2. Aromas: ethyl acetate, isoamyl acetate y ethyl octanoate.
3. Volatilización: partición gas-líquido Antoine+UNIFAC y stripping por flujo CO2.
4. CO2: estado de CO2 disuelto y transición macroscópica de oxígeno.

Benchmark de modelos secundarios piloto:

{md_table(pilot_sec, ["model", "data_wsse", "n_parameters", "bic", "active_bound_count", "state_penalty", "selection_score"], 8)}

Benchmark aromático piloto:

{md_table(pilot_aroma, ["model", "n_parameters", "data_wsse", "bic", "active_bound_count", "ea_retained_relative_rmse"], 8)}

Benchmark CO2 piloto:

{md_table(pilot_co2, ["mode", "mechanistic_class", "data_wsse", "bic", "active_bound_count", "selection_penalty", "selection_score", "selected"], 10)}

## 9. KPIs y resultados piloto

Métricas finales del modelo calibrado piloto:

{md_table(pilot_metrics, ["group", "state", "species", "pool", "n", "relative_rmse", "relative_bias", "corr"], 30)}

Figuras:

- `figures/fig_06_pilot_fit_metrics.png`
- `figures/fig_15_pilot_final_fit_25170.png`
- `figures/fig_16_pilot_final_fit_25171.png`
- `figures/fig_17_pilot_co2_benchmark_25170.png`
- `figures/fig_18_pilot_co2_benchmark_25171.png`

Interpretación de escalabilidad:

- Glicerol y etanol muestran transferencia razonable de estructura primaria.
- Azúcares y nitrógeno son más sensibles a matriz natural, medición y eventos operacionales.
- CO2 online no es explicable por la cinética primaria sola: requiere solubilidad, transición de O2 y liberación gas-líquido.
- Aromas requieren estructura propia y datos piloto; no son extrapolables desde mosto sintético sin medición aromática.

## 10. Diseño experimental y escalabilidad

Campaña MBDoE laboratorio:

{md_table(lab_campaign, ["campaign_order", "candidate", "family", "medium", "temperature_segments", "score", "rationale"], 12)}

Campaña MBDoE piloto natural:

{md_table(pilot_campaign, ["campaign_order", "candidate", "family", "medium", "temperature_segments", "N_pulses_kg_m3", "score"], 10)}

Lectura técnica:

- Laboratorio sintético maximiza separación de mecanismos por manipulación de composición inicial, azúcar, biomasa y etanol.
- Piloto natural prioriza temperatura, N y condiciones compatibles con mosto real.
- La escalabilidad se valida combinando ambas escalas: sintético para identificabilidad y piloto para validez de matriz/proceso.

## 11. Matriz de decisión de transferibilidad

{md_table(transfer)}

## 12. Conclusiones

1. La estructura primaria desarrollada con mosto sintético VL3 es útil y parcialmente transferible a escala piloto.
2. La transferencia directa no es suficiente para los outputs de bodega piloto: CO2, O2, aromas y condensado requieren extensiones fenomenológicas.
3. La validación piloto 2025 entregó evidencia real de escalabilidad operacional: múltiples vinificaciones, tratamientos térmicos/nutricionales, datos aromáticos, CO2 online y eventos operativos.
4. El modelo escalado debe usarse con parámetros piloto recalibrados y no con parámetros sintéticos puros.
5. La campaña MBDoE sugiere continuar con un approach mixto: sintético para excitar direcciones poco identificables y piloto natural para validar transferencia.

## 13. Pasos a seguir

- Completar bitácora física de bodega piloto: volumen, geometría, sensores, configuración térmica y protocolo de muestreo.
- Ejecutar nuevos ensayos piloto con CO2 online confiable y condensado final.
- Recalibrar el modelo integrado tras los ensayos MBDoE.
- Separar YAN en amonio/PAN si se busca explicar remanentes de nitrógeno.
- Definir qué parámetros quedan libres en MPCC y cuáles se fijan por robustez.
"""
    (root / "A15_technical_report.md").write_text(report, encoding="utf-8")


def create_readme(root: Path) -> None:
    readme = """# Bundle A15 - Validación de escalabilidad en bodega piloto vendimia 2025

Este paquete resume la evidencia técnica para el anexo A15, Actividad 7.18.

## Elementos principales

- `A15_technical_report.md`: informe técnico autoexplicativo.
- `figures/fig_01_*` a `fig_12_*`: figuras de mayor valor para el anexo.
- `tables/02_*`, `03_*`, `05_*`, `07_*`, `11_*`, `22_*`: tablas clave para batches, KPIs, calibración, FIM, piloto y decisión de transferibilidad.

## Elementos complementarios

- `figures/fig_13_*` a `fig_20_*`: ajustes representativos, benchmark CO2 y espectros FIM.
- `source_data/`: workbooks fuente compactos incluidos como respaldo.
- `source_snapshot/`: scripts usados para generar resultados y trazabilidad.

## Lectura sugerida

1. Leer `A15_technical_report.md`.
2. Revisar `fig_01`, `fig_02`, `fig_05`, `fig_06`, `fig_09` y `fig_10`.
3. Usar `tables/22_transferability_decision_matrix.csv` para conclusiones ejecutivas.
"""
    (root / "README_bundle.md").write_text(readme, encoding="utf-8")


def create_manifest(root: Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file():
            rows.append({"relative_path": str(path.relative_to(root)).replace("\\", "/"), "size_bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(root / "bundle_manifest.csv", index=False)


def zip_bundle(root: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(root.parent))


def main() -> None:
    figures, tables_dir, source, data_dir, root = ensure_clean_bundle()
    tables = load_tables()
    write_tables(tables, tables_dir)
    create_figures(tables, figures)
    write_source_snapshot(source, data_dir)
    create_report(tables, root)
    create_readme(root)
    create_manifest(root)
    zip_bundle(root, ZIP_PATH)
    print(root)
    print(ZIP_PATH)


if __name__ == "__main__":
    main()
