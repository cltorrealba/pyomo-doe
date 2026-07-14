from __future__ import annotations

import math
import re
import shutil
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
RESULTS_DIR = SCRIPT_DIR / "results"
BUNDLE_NAME = "A12_vinificaciones_protocolos_aromas_kpi_bundle_20260630"
BUNDLE_DIR = RESULTS_DIR / BUNDLE_NAME
FIG_DIR = BUNDLE_DIR / "figures"
TABLE_DIR = BUNDLE_DIR / "tables"
TEMPLATE_DIR = BUNDLE_DIR / "templates"
TRACE_DIR = BUNDLE_DIR / "source_trace"

A09_ROOT = Path(
    r"c:\Users\ctorrealba\OneDrive - Viña Concha y Toro S.A\Documentos\Proyectos I+D\PI-4497\Ley I+D\Rendición Técnica 6\Anexos\A09 - Actividad 4.5-4.7"
)
A09_BUNDLE = A09_ROOT / "A09_optimizacion_dfba_kkt_mcdm_bundle_20260629"
A09_PDF = A09_ROOT / "Anexo A09.pdf"
PROJECT_A12_DIR = A09_ROOT.parent / "A12 - Actividad 5.6-5.8"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_new_must_glycerol_estimability_doe as core


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def as_float(value: object) -> float:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return np.nan
    if isinstance(value, (int, float, np.number)):
        return float(value)
    text = str(value).strip().replace(",", ".")
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return float(match.group(0)) if match else np.nan


def isoamyl_to_mg_l(value: object) -> float:
    value_float = as_float(value)
    if not np.isfinite(value_float):
        return np.nan
    text = str(value).lower()
    if "mg/l" in text:
        return value_float
    if value_float < 0.1:
        return value_float * 1000.0
    return value_float


def fmt(value: object, ndigits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, str):
        return value
    value = float(value)
    if abs(value) >= 10000 or (abs(value) > 0 and abs(value) < 0.001):
        return f"{value:.{ndigits}e}"
    return f"{value:.{ndigits}f}"


def md_table(frame: pd.DataFrame, max_rows: int | None = None) -> str:
    if frame.empty:
        return "_No data available._"
    out = frame.copy()
    if max_rows is not None:
        out = out.head(max_rows)
    for col in out.columns:
        if pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(fmt)
    return out.to_markdown(index=False)


def prepare_dirs() -> None:
    if BUNDLE_DIR.exists():
        # Refresh files in place while preserving directory permissions.
        for path in sorted(BUNDLE_DIR.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                try:
                    path.rmdir()
                except OSError:
                    pass
    for directory in (BUNDLE_DIR, FIG_DIR, TABLE_DIR, TEMPLATE_DIR, TRACE_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def build_compact_tables() -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}

    tables["validation_metrics"] = read_csv(
        RESULTS_DIR / "A08_independent_validation_bundle_2026-06-29" / "tables" / "synthetic_validation_metrics_by_state.csv"
    )
    tables["curve_metrics"] = read_csv(RESULTS_DIR / "curve_validation" / "current_fit_curve_metrics_by_medium_state.csv")
    tables["fit_summary"] = read_csv(RESULTS_DIR / "new_must_glycerol_estimability_doe" / "fit_summary.csv")
    tables["profile_likelihood"] = read_csv(
        RESULTS_DIR / "new_must_glycerol_overnight_validation" / "profile_full_summary.csv"
    )
    tables["volume_policy"] = read_csv(RESULTS_DIR / "final_operational_doe_volume_constrained" / "policy_summary.csv")
    tables["selected_campaign"] = read_csv(RESULTS_DIR / "final_operational_doe_volume_constrained" / "selected_campaign.csv")
    tables["campaign_protocol"] = read_csv(RESULTS_DIR / "final_operational_doe_volume_constrained" / "campaign_protocol.csv")
    tables["operational_schedule"] = read_csv(RESULTS_DIR / "final_operational_doe_volume_constrained" / "operational_schedule.csv")
    tables["volume_audit"] = read_csv(RESULTS_DIR / "final_operational_doe_volume_constrained" / "volume_audit.csv")
    tables["stock_assumptions"] = read_csv(
        RESULTS_DIR / "final_operational_doe_volume_constrained" / "stock_concentration_assumptions.csv"
    )
    tables["lot1_sampling"] = read_csv(RESULTS_DIR / "lot1_express_optimal_sampling" / "express_optimal_sampling_summary.csv")
    tables["aroma_selected"] = read_csv(RESULTS_DIR / "aroma_joint_campaign_doe" / "aroma_campaign_selected.csv")
    tables["aroma_parameter_reduction"] = read_csv(
        RESULTS_DIR / "aroma_joint_campaign_doe" / "aroma_campaign_parameter_reduction.csv"
    )
    tables["aroma_dopt_benchmark"] = read_csv(RESULTS_DIR / "aroma_joint_campaign_doe" / "aroma_dopt_benchmark_summary.csv")
    tables["secondary_fit"] = read_csv(RESULTS_DIR / "secondary_v2_model_evaluation" / "fit_comparison.csv")
    tables["a09_key_results"] = read_csv(A09_BUNDLE / "tablas" / "RESUMEN_RESULTADOS_CLAVE_A09.csv")
    tables["a09_mpcc_summary"] = read_csv(A09_BUNDLE / "tablas" / "optimized_control_summary.csv")
    tables["a09_policy_selection"] = read_csv(A09_BUNDLE / "tablas" / "dynamic_control_policy_selection.csv")
    tables["a09_pareto_summary"] = read_csv(A09_BUNDLE / "tablas" / "simplified_mpcc_pareto_summary.csv")

    # A12 status matrix.
    tables["activity_status"] = pd.DataFrame(
        [
            {
                "activity": "5.6 Vinifications: optimal and commercial protocols",
                "current_status": "In progress / technically gated",
                "evidence_in_bundle": "DOE operational protocol, volume audit, MPCC preliminary candidates",
                "claim_level": "Protocol design and readiness; not final in-vivo optimal-vs-commercial comparison",
            },
            {
                "activity": "5.7 Aroma markers in wine and condensates",
                "current_status": "Analytical design and MBDoE completed; final campaign data pending",
                "evidence_in_bundle": "Aroma DOE, target markers, terminal condensate balance, QA/QC template",
                "claim_level": "Analytical plan and model-based prioritization; not final marker statistics",
            },
            {
                "activity": "5.8 In-vivo KPI comparison",
                "current_status": "KPI framework prepared; experimental comparison pending calibrated protocols",
                "evidence_in_bundle": "KPI schema, MPCC preliminary objective metrics, DOE metrics",
                "claim_level": "Framework and preliminary model-based indicators; not final in-vivo KPI delta",
            },
        ]
    )

    tables["kpi_schema"] = pd.DataFrame(
        [
            ("fermentation_time_h", "time to residual sugar gate or endpoint", "h", "shorter is better under quality constraints"),
            ("final_residual_sugar_g_l", "G+F at endpoint", "g/L", "primary dryness/productivity metric"),
            ("ethanol_g_l", "ethanol at endpoint", "g/L", "yield and product specification"),
            ("glycerol_g_l", "glycerol at endpoint", "g/L", "quality/metabolic side product"),
            ("co2_cumulative_g", "integrated CO2 sensor signal", "g or mol", "fermentation rate and mass balance"),
            ("yan_consumed_mg_l", "initial YAN minus terminal YAN", "mg/L", "nitrogen use efficiency"),
            ("nutrient_added_mg_l", "sum of nitrogen additions", "mg/L", "input cost/process intervention"),
            ("temperature_energy_proxy", "mean squared deviation from ambient/reference", "dimensionless", "thermal operation cost"),
            ("aroma_wine_ug_l", "target aroma concentration in wine", "ug/L", "chemical quality marker"),
            ("aroma_condensate_ug", "terminal trapped aroma mass", "ug", "volatilized-loss closure"),
            ("aroma_retention_ratio", "wine mass / total produced mass", "dimensionless", "aroma retention efficiency"),
            ("volume_removed_ml", "liquid removed by samples", "mL", "sampling burden/operational feasibility"),
        ],
        columns=["kpi", "definition", "unit", "interpretation"],
    )

    tables["qaqc_plan"] = pd.DataFrame(
        [
            ("chemical_calibration", "external calibration curve plus internal standard where applicable", "R2 >= 0.99 or method SOP threshold"),
            ("blank_control", "solvent/process blank for aroma and condensate analysis", "below LOQ or documented correction"),
            ("duplicate_sample", "at least one duplicate per analytical batch or critical timepoint", "RSD <= method threshold"),
            ("spike_recovery", "matrix spike for representative wine/condensate", "70-130% unless method-specific criterion applies"),
            ("mass_balance", "wine + condensate + residual/gas closure for target aromatics where possible", "qualitative/quantitative closure reported"),
            ("sensor_qc", "CO2 zero/span check and timestamp synchronization", "documented before/after run"),
            ("sample_traceability", "reactor, lot, treatment, time, volume, preservative, storage", "complete chain-of-custody fields"),
        ],
        columns=["qa_qc_item", "procedure", "acceptance_or_action"],
    )

    tables["execution_log_template"] = pd.DataFrame(
        columns=[
            "lot",
            "reactor_id",
            "treatment_id",
            "planned_event",
            "planned_datetime",
            "actual_datetime",
            "relative_time_h",
            "event_type",
            "sample_type",
            "volume_removed_ml",
            "control_channel",
            "stock_concentration",
            "stock_volume_ml",
            "target_dose",
            "actual_dose",
            "temperature_setpoint_c",
            "measured_temperature_c",
            "co2_sensor_status",
            "incident_flag",
            "incident_description",
            "operator",
            "use_in_fit",
        ]
    )

    tables["raw_data_inventory"] = pd.DataFrame(
        [
            ("mosto_sintetico_vl3.xlsx", "historical/current synthetic must data", "used for calibration/validation", "not copied; summarized by normalized tables"),
            ("mosto_natural_xthiol.xlsx", "natural must transfer data", "used for model comparison and transfer diagnosis", "not copied; summarized by normalized tables"),
            ("A09 bundle/PDF", "MPCC, KKT, MCDM and preliminary protocol optimization", "cross-reference for optimal protocol work", "selected docs/tables/figures copied"),
            ("new in-vivo lot logs", "actual vinification execution records", "pending ingestion", "templates included"),
            ("aroma wine/condensate reports", "GC or equivalent analytical reports", "pending ingestion", "QA/QC and schema included"),
        ],
        columns=["source", "description", "role", "bundle_handling"],
    )

    for name, frame in tables.items():
        if not frame.empty or name.endswith("_template"):
            frame.to_csv(TABLE_DIR / f"{name}.csv", index=False)
    return tables


def read_theta(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, index_col=0)
    return {str(k): float(v) for k, v in frame.iloc[:, 0].to_dict().items()}


def make_model_predicted_observed_figure() -> None:
    theta_path = RESULTS_DIR / "new_must_glycerol_overnight_validation" / "theta_multistart_07.csv"
    if not theta_path.exists():
        theta_path = RESULTS_DIR / "new_must_glycerol_overnight_validation" / "theta_reference.csv"
    if not theta_path.exists():
        fallback = RESULTS_DIR / "curve_validation" / "plots" / "current_fit" / "fitcmp_synthetic_MS007.png"
        copy_if_exists(fallback, FIG_DIR / "fig02_model_validation_predicted_observed.png")
        return
    theta = read_theta(theta_path)
    data = core.load_normalized_data()
    batches = [batch for batch in core.make_batches(data) if batch.medium == "synthetic"]
    rows = []
    for batch in batches:
        sim = core.simulate(batch, theta, batch.time)
        if sim is None:
            continue
        for state in core.STATE_NAMES:
            obs = np.asarray(batch.observations[state], dtype=float)
            pred = sim.loc[batch.time, state].to_numpy(dtype=float)
            mask = np.isfinite(obs)
            for observed, predicted in zip(obs[mask], pred[mask]):
                rows.append({"state": state, "observed": float(observed), "predicted": float(predicted)})
    predobs = pd.DataFrame(rows)
    predobs.to_csv(TABLE_DIR / "model_predicted_observed_synthetic_long.csv", index=False)
    if predobs.empty:
        return
    states = list(core.STATE_NAMES)
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    axes = axes.ravel()
    for ax, state in zip(axes, states):
        group = predobs[predobs["state"].eq(state)]
        if group.empty:
            ax.axis("off")
            continue
        x = group["observed"].to_numpy(dtype=float)
        y = group["predicted"].to_numpy(dtype=float)
        lo = min(np.nanmin(x), np.nanmin(y))
        hi = max(np.nanmax(x), np.nanmax(y))
        pad = 0.05 * (hi - lo if hi > lo else 1.0)
        ax.scatter(x, y, s=16, alpha=0.75, color="#1f77b4", edgecolor="none")
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="black", lw=1)
        ax.set_title(state)
        ax.set_xlabel("Observed")
        ax.set_ylabel("Predicted")
    for ax in axes[len(states) :]:
        ax.axis("off")
    fig.suptitle("Predicted vs observed, synthetic must validation", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig02_model_validation_predicted_observed.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def make_selected_campaign_montage() -> None:
    images = sorted((RESULTS_DIR / "vc_doe_nb_plots").glob("doe_*.png"))
    if not images:
        images = sorted((RESULTS_DIR / "final_operational_doe_volume_constrained" / "plots").glob("volume_profile_*.png"))
    if not images:
        return
    selected = images[:9]
    ncols = 3
    nrows = int(math.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.6 * nrows))
    axes = np.asarray(axes).ravel()
    for ax, path in zip(axes, selected):
        img = plt.imread(path)
        ax.imshow(img)
        ax.set_title(path.stem.replace("_", " ")[:45])
        ax.axis("off")
    for ax in axes[len(selected) :]:
        ax.axis("off")
    fig.suptitle("Selected campaign input and sampling designs", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig08_selected_campaign_designs_montage.png", dpi=180)
    plt.close(fig)


def make_figures(tables: dict[str, pd.DataFrame]) -> None:
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    status = tables["activity_status"]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    y = np.arange(len(status))
    colors = ["#fdae6b", "#9ecae1", "#c7e9c0"]
    ax.barh(y, [0.72, 0.65, 0.55], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(["5.6 protocols", "5.7 aroma", "5.8 KPI"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Evidence maturity, qualitative")
    ax.set_title("A12 current evidence status: design-ready, final in-vivo comparison pending")
    for idx, text in enumerate(status["current_status"]):
        ax.text(0.02, idx, text, va="center", ha="left", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig01_a12_activity_status.png", dpi=180)
    plt.close(fig)

    make_model_predicted_observed_figure()

    policy = tables["volume_policy"]
    if not policy.empty:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        names = policy["policy"].astype(str).to_list()
        axes[0].bar(names, pd.to_numeric(policy["campaign_logdet"]), color="#3182bd")
        axes[0].set_title("Campaign logdet")
        axes[1].bar(names, pd.to_numeric(policy["campaign_trace_inv"]), color="#756bb1")
        axes[1].set_title("Trace(inv FIM)")
        axes[2].bar(names, pd.to_numeric(policy["new_param_mean_var_reduction"]), color="#31a354")
        axes[2].set_title("Mean variance reduction")
        for ax in axes:
            ax.tick_params(axis="x", rotation=25)
        fig.suptitle("Operational DOE policy benchmark")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig03_doe_policy_benchmark.png", dpi=180)
        plt.close(fig)

    volume = tables["volume_audit"]
    if not volume.empty:
        fig, ax = plt.subplots(figsize=(10, 4.8))
        order = volume["campaign_order"].astype(str)
        final_volume = pd.to_numeric(volume["final_volume_ml"], errors="coerce")
        removed = pd.to_numeric(volume["net_volume_removed_ml"], errors="coerce")
        ax.bar(order, final_volume, label="Final volume", color="#74c476")
        ax.bar(order, removed, bottom=final_volume, label="Removed/input net volume", color="#fd8d3c", alpha=0.8)
        ax.axhline(1100, color="black", ls="--", lw=1, label="Minimum final volume target")
        ax.set_xlabel("Campaign order")
        ax.set_ylabel("mL")
        ax.set_title("2 L reactor volume audit for proposed campaign")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig04_volume_audit.png", dpi=180)
        plt.close(fig)

    a09 = tables["a09_key_results"].copy()
    if not a09.empty:
        a09["sugar_g_l_num"] = a09["azucar_final_g_L"].map(as_float)
        a09["isoamyl_mg_l_num"] = a09["isoamyl_acetate"].map(isoamyl_to_mg_l)
        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        for category, group in a09.groupby("categoria"):
            ax.scatter(group["sugar_g_l_num"], group["isoamyl_mg_l_num"], s=70, label=category, alpha=0.85)
            for _, row in group.iterrows():
                label = str(row["candidato_o_etapa"])[:18]
                ax.annotate(label, (row["sugar_g_l_num"], row["isoamyl_mg_l_num"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel("Terminal sugar reported by A09 (g/L)")
        ax.set_ylabel("Isoamyl acetate, normalized when possible (mg/L)")
        ax.set_title("A09 preliminary MPCC/sequential tradeoff: not final in-vivo protocol")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig05_a09_preliminary_mpcc_tradeoff.png", dpi=180)
        plt.close(fig)

    aroma_red = tables["aroma_parameter_reduction"]
    if not aroma_red.empty:
        red = aroma_red.copy()
        red["campaign_var_reduction"] = pd.to_numeric(red["campaign_var_reduction"], errors="coerce")
        fig, ax = plt.subplots(figsize=(10, 6.5))
        colors = red["group"].map({"fermentation": "#3182bd", "synthesis": "#e6550d"}).fillna("#969696")
        ax.barh(red["parameter"], red["campaign_var_reduction"], color=colors)
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("Expected variance reduction")
        ax.set_title("Aroma/fermentation DOE: expected parameter variance reduction")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig06_aroma_parameter_reduction.png", dpi=180)
        plt.close(fig)

    aroma_bench = tables["aroma_dopt_benchmark"]
    if not aroma_bench.empty:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        x = aroma_bench["design"].astype(str)
        axes[0].bar(x, pd.to_numeric(aroma_bench["logdet"], errors="coerce"), color="#3182bd")
        axes[0].set_title("logdet")
        axes[1].bar(x, pd.to_numeric(aroma_bench["condition_number"], errors="coerce"), color="#756bb1")
        axes[1].set_title("condition number")
        axes[2].bar(x, pd.to_numeric(aroma_bench["trace_inv"], errors="coerce"), color="#31a354")
        axes[2].set_title("trace(inv)")
        for ax in axes:
            ax.tick_params(axis="x", rotation=25)
        fig.suptitle("Aroma DOE benchmark: hybrid vs D-opt variants")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig07_aroma_doe_benchmark.png", dpi=180)
        plt.close(fig)

    make_selected_campaign_montage()

    # Copy selected source figures from A09 and model pipeline.
    source_figures = {
        "source_a09_optimized_control_profiles.png": A09_BUNDLE / "figuras" / "optimized_control_profiles.png",
        "source_a09_pareto_sugar_aroma.svg": A09_BUNDLE / "figuras" / "pareto_sugar_aroma.svg",
        "source_a09_pareto_nutrient_aroma.svg": A09_BUNDLE / "figuras" / "pareto_nutrient_aroma.svg",
        "source_a09_control_movement_overview.svg": A09_BUNDLE / "figuras" / "control_movement_overview.svg",
        "source_aroma_eigen_spectrum.png": RESULTS_DIR / "aroma_joint_campaign_doe" / "aroma_eigen_spectrum.png",
        "source_secondary_fit_synthetic_MS016.png": RESULTS_DIR / "secondary_v2_model_evaluation" / "plots" / "secondary_v2_fit_synthetic_MS016.png",
    }
    for name, src_path in source_figures.items():
        copy_if_exists(src_path, FIG_DIR / name)


def copy_source_trace() -> None:
    source_files = {
        "A09_Anexo_pdf_reference.pdf": A09_PDF,
        "A09_README_BUNDLE.md": A09_BUNDLE / "README_BUNDLE.md",
        "A09_CLAIMS_Y_LIMITACIONES.md": A09_BUNDLE / "CLAIMS_Y_LIMITACIONES.md",
        "A09_METODOLOGIA_FORMULACION.md": A09_BUNDLE / "documentos" / "METODOLOGIA_FORMULACION.md",
        "A09_RESULTADOS_Y_TRAZABILIDAD.md": A09_BUNDLE / "documentos" / "RESULTADOS_Y_TRAZABILIDAD.md",
        "A09_METODO_MULTICRITERIO.md": A09_BUNDLE / "documentos" / "METODO_MULTICRITERIO.md",
        "A08_validation_report.md": RESULTS_DIR / "A08_independent_validation_bundle_2026-06-29" / "A08_report.md",
        "curve_validation_decision_summary.md": RESULTS_DIR / "curve_validation" / "curve_validation_decision_summary.md",
        "volume_constrained_doe_report.md": RESULTS_DIR / "final_operational_doe_volume_constrained" / "volume_constrained_doe_report.md",
        "aroma_campaign_report.md": RESULTS_DIR / "aroma_joint_campaign_doe" / "aroma_campaign_report.md",
        "secondary_v2_model_report.md": RESULTS_DIR / "secondary_v2_model_evaluation" / "secondary_v2_model_report.md",
        "lot1_express_optimal_sampling_report.md": RESULTS_DIR / "lot1_express_optimal_sampling" / "express_optimal_sampling_report.md",
    }
    for name, src in source_files.items():
        copy_if_exists(src, TRACE_DIR / name)


def write_templates() -> None:
    # CSV templates are also stored in tables for direct spreadsheet opening.
    for name in ("execution_log_template", "qaqc_plan", "kpi_schema"):
        src = TABLE_DIR / f"{name}.csv"
        if src.exists():
            copy_if_exists(src, TEMPLATE_DIR / f"{name}.csv")


def build_report(tables: dict[str, pd.DataFrame]) -> None:
    fit_summary = tables["fit_summary"]
    fit_short = (
        fit_summary[["fit", "n_batches", "mediums", "n_parameters", "final_wsse", "n_residuals", "wsse_per_residual", "l2_lambda"]]
        if not fit_summary.empty
        else fit_summary
    )
    profile = tables["profile_likelihood"]
    profile_short = (
        profile[["parameter", "max_lr_stat", "crosses_left_95", "crosses_right_95", "profile_identifiable_95"]]
        if not profile.empty
        else profile
    )
    validation = tables["validation_metrics"]
    validation_short = (
        validation[validation.get("theta", pd.Series(dtype=str)).eq("best_multistart_07")][
            ["state", "n", "rmse", "mae", "bias", "weighted_rmse", "r2"]
        ]
        if not validation.empty and "theta" in validation.columns
        else validation
    )
    campaign_protocol = tables["campaign_protocol"]
    protocol_short = (
        campaign_protocol[
            ["campaign_order", "lot", "start_datetime", "candidate", "medium", "family", "horizon_h"]
        ]
        if not campaign_protocol.empty
        else campaign_protocol
    )
    selected = tables["selected_campaign"]
    selected_short_cols = [
        c
        for c in [
            "campaign_order",
            "candidate",
            "family",
            "medium",
            "horizon_h",
            "campaign_logdet",
            "new_param_mean_var_reduction",
            "new_param_worst_var_reduction",
            "rationale",
        ]
        if c in selected.columns
    ]
    selected_short = selected[selected.get("objective", pd.Series(dtype=str)).eq("hybrid")][selected_short_cols] if selected_short_cols and not selected.empty else selected
    a09_key = tables["a09_key_results"]
    a09_key_short = a09_key.copy()
    if not a09_key_short.empty:
        a09_key_short["isoamyl_acetate_mg_l_interpreted"] = a09_key_short["isoamyl_acetate"].map(isoamyl_to_mg_l)
    mpcc = tables["a09_mpcc_summary"]
    policy = tables["volume_policy"]
    aroma_selected = tables["aroma_selected"]
    aroma_red = tables["aroma_parameter_reduction"]
    aroma_bench = tables["aroma_dopt_benchmark"]
    secondary = tables["secondary_fit"]

    report = r"""# A12 - Vinificaciones, aromas y KPIs in vivo

**Actividades:** 5.6, 5.7 y 5.8  
**Titulo integrado:** Vinificaciones con protocolos optimos y comerciales + analisis de marcadores aromaticos en vino y condensados + evaluacion de indicadores clave de eficiencia in vivo  
**Fecha de generacion:** @@DATE@@  
**Bundle:** `@@BUNDLE@@`

## 1. Proposito y criterio de reporte

Este bundle consolida la evidencia tecnica disponible para el Anexo A12. El objetivo es documentar el avance metodologico y experimental hacia la comparacion in vivo de protocolos optimos versus comerciales, incluyendo analisis de aromas y KPIs de eficiencia.

La diferencia clave respecto de un anexo experimental cerrado es que los protocolos optimos realistas aun no deben declararse finales. La razon tecnica es que el MPCC/dFBA/KKT disponible en A09 fue construido sobre modelo nominal y objetivos provisionales. Antes de ejecutar una comparacion in vivo definitiva, se requirio robustecer la calibracion del modelo de fermentacion, evaluar identificabilidad practica, incorporar glicerinol, biomasa muerta, CO2, estados secundarios y capa de aromas, y disenar una campana MBDoE que mejore la estimabilidad de parametros.

Por tanto, este anexo debe defender tres puntos:

1. Existe metodologia y software para generar politicas dinamicas optimas.
2. Existe una campana experimental operacionalizable para recalibrar y validar el gemelo antes de cerrar protocolos optimos.
3. La postergacion de la comparacion in vivo final no es una omision, sino una decision tecnica para evitar artefactos de optimizacion.

## 2. Estado de cumplimiento por actividad

@@TABLE_ACTIVITY_STATUS@@

Figura asociada: `figures/fig01_a12_activity_status.png`.

## 3. Evidencia fuente consolidada

La evidencia usada proviene de cuatro capas:

- Datos historicos y nuevos de mosto sintetico/natural procesados para calibracion ODE.
- Validacion independiente A08 del Gemelo Digital.
- Diseno experimental MBDoE operacionalizado para nuevos ensayos.
- Bundle A09 de optimizacion dFBA/KKT/MPCC/MCDM.

Archivos trazables se incluyen en `source_trace/`. Las tablas compactas estan en `tables/` y las figuras principales en `figures/`.

## 4. Modelo de simulacion fermentativa

El modelo operativo principal es un sistema ODE:

$$
x(t)=\left[X, X_d, N, G, F, E, Gly\right]^T
$$

con estados de biomasa viable, biomasa muerta, nitrogeno asimilable, glucosa, fructosa, etanol y glicerol. La dinamica es:

$$
\frac{dX}{dt}=(\mu-k_d)X+u_X(t)
$$

$$
\frac{dX_d}{dt}=k_dX
$$

$$
\frac{dN}{dt}=-q_N a_{\mu}(T)\frac{N}{N+K_N(T)}X+u_N(t)
$$

$$
\frac{dG}{dt}=-\left[q_{XG}a_{\mu}(T)\frac{N}{N+K_N(T)}
+q_{EG}a_{\beta}(T)\frac{G}{G+K_G(T)}I_E(E)
+m(T)\frac{G}{G+F}\right]X+u_G(t)
$$

$$
\frac{dF}{dt}=-\left[q_{XF}a_{\mu}(T)\frac{N}{N+K_N(T)}
+q_{EF}a_{\beta}(T)\frac{F}{F+K_F(T)}I_G(G)I_E(E)
+m(T)\frac{F}{G+F}\right]X+u_F(t)
$$

$$
\frac{dE}{dt}=(\beta_G+\beta_F)X+u_E(t)
$$

$$
\frac{dGly}{dt}=\left[\gamma_{G0}a_{\beta}(T)\frac{G}{G+K_G(T)}I_E(E)
+\gamma_{F0}a_{\beta}(T)\frac{F}{F+K_F(T)}I_G(G)I_E(E)\right]X
$$

con inhibiciones:

$$
I_E(E)=\frac{1}{1+i_E(T)E},\qquad I_G(G)=\frac{1}{1+i_G(T)G}
$$

Los pulsos operacionales se modelan como entradas suaves:

$$
u_j(t)=\sum_k \Delta C_{j,k}\frac{\exp[-((t-t_{j,k})/w)^2]}{\sqrt{\pi}w}
$$

## 5. Calibracion, validacion y estimabilidad

La calibracion se formulo como minimos cuadrados ponderados:

$$
\min_{\theta} J(\theta)=\sum_{b,i,k}\left(\frac{\hat{y}_{b,i}(t_k;\theta)-y_{b,i,k}}{\sigma_i(y_{b,i,k})}\right)^2
+\lambda\sum_{p\in P_R}\left[\log\left(\frac{\theta_p}{\theta_{p,ref}}\right)\right]^2
$$

donde:

$$
\sigma_i(y)=\max(\sigma_{i,floor},r_i|y|)
$$

La estimabilidad se evaluo con FIM, perfiles de likelihood y aproximacion bayesiana/Laplace:

$$
F(\theta)=S(\theta)^T W S(\theta),\qquad
LR_p(\theta_p)=J(\theta_p,\hat{\theta}_{-p})-J(\hat{\theta})
$$

El umbral usado fue:

$$
LR_p>\chi^2_{1,0.95}=3.841
$$

**Resumen de ajustes**

@@TABLE_FIT@@

**Metricas de validacion sintetica**

@@TABLE_VALIDATION@@

Figura central de validacion: `figures/fig02_model_validation_predicted_observed.png`.

**Perfil likelihood**

@@TABLE_PROFILE@@

## 6. Modelo dFBA/KKT/MPCC y decision multicriterio

El trabajo A09 implemento una formulacion simultanea dFBA/KKT/MPCC. El simulador secuencial resuelve:

$$
\frac{dx}{dt}=f_{kin}(x,u,p)+Bv
$$

con un problema FBA local:

$$
\max_v c^Tv
$$

$$
S v=0,\qquad l(x,u,p)\leq v\leq q(x,u,p)
$$

La reformulacion KKT usa duales \alpha,\beta y estacionariedad:

$$
c-S^T\lambda-\alpha+\beta=0
$$

con complementariedad:

$$
\alpha_i(v_i-l_i)=0,\qquad \beta_i(q_i-v_i)=0
$$

En el MPCC se usa relajacion tipo Scholtes:

$$
0\leq \alpha_i(v_i-l_i)\leq \tau,\qquad
0\leq \beta_i(q_i-v_i)\leq \tau
$$

El objetivo dinamico A09 incluyo productividad/azucar residual, aroma, nutrientes, energia termica y suavidad:

$$
\min J=w_sJ_{sugar}+w_aJ_{aroma}+w_nJ_{nutrient}+w_eJ_{energy}+w_{sm}J_{smooth}
$$

El bundle A09 demostro que el pipeline puede generar candidatos dinamicos y frentes tipo Pareto, pero sus propios claims indican que no debe presentarse como decision industrial final ni optimo global. La politica MPCC promovible principal fue local, numericamente sana y generada bajo modelo nominal.

**Resultados A09 clave**

@@TABLE_A09_KEY@@

Figura asociada: `figures/fig05_a09_preliminary_mpcc_tradeoff.png`.

**Candidato MPCC MVP**

@@TABLE_MPCC@@

Figura fuente: `figures/source_a09_optimized_control_profiles.png`.

## 7. Protocolo experimental propuesto para A12

La campana experimental actual se estructura como 9 fermentaciones en 3 lotes de 3 reactores. No corresponde aun a la comparacion final "protocolo optimo vs comercial", sino a la campana necesaria para recalibrar y seleccionar con mayor confianza esos protocolos.

**Lotes y tratamientos propuestos**

@@TABLE_PROTOCOL@@

**Campana seleccionada por MBDoE**

@@TABLE_SELECTED@@

El criterio principal fue D-optimalidad:

$$
\Phi_D(d)=\log\det(F(\theta,d)+\epsilon I)
$$

con criterios complementarios:

$$
\Phi_A(d)=tr(F^{-1})
$$

La politica con restriccion de volumen seleccionada fue `full14_small6`, que combina muestras completas y pequenas para mantener volumen final viable en reactores de 2 L.

**Benchmark de politicas de muestreo/volumen**

@@TABLE_POLICY@@

Figuras:

- `figures/fig03_doe_policy_benchmark.png`
- `figures/fig04_volume_audit.png`
- `figures/fig08_selected_campaign_designs_montage.png`

**Stocks y restricciones operacionales**

@@TABLE_STOCKS@@

## 8. Analisis de aromas, condensados y QA/QC

La capa de aromas considera, como minimo, etil acetato, isoamyl acetate y ethyl octanoate en el modelo, y deja trazabilidad para extender a bencil alcohol, benzaldehido, decanoato de etilo, hexil acetato, phenylethyl acetate y otros marcadores disponibles.

La formulacion usada separa:

- produccion liquida aparente dependiente de fase metabolica;
- perdida gas-liquido fijada por literatura/UNIFAC o surrogate validado;
- CO2 online como senal de tasa de fermentacion;
- condensado terminal como restriccion de cierre de masa, no como dinamica gas en linea.

El balance conceptual por aroma `a` es:

$$
\frac{dC_{a,L}}{dt}=r_{a,prod}(x,\theta_a,T)-r_{a,loss}(C_{a,L},T,F_{CO2})
$$

$$
\frac{dM_{a,cond}}{dt}=\eta_{trap}\,r_{a,loss}(C_{a,L},T,F_{CO2})V_L
$$

Los parametros de particion se fijan por literatura para evitar confundir equilibrio gas-liquido con sintesis biologica. Los parametros libres principales quedan asociados a sintesis por fase.

**Aroma DOE seleccionado**

@@TABLE_AROMA_SELECTED@@

**Benchmark aroma D-opt**

@@TABLE_AROMA_BENCH@@

**Reduccion esperada de varianza**

@@TABLE_AROMA_REDUCTION@@

Figuras:

- `figures/fig06_aroma_parameter_reduction.png`
- `figures/fig07_aroma_doe_benchmark.png`
- `figures/source_aroma_eigen_spectrum.png`

**Plan QA/QC**

@@TABLE_QAQC@@

## 9. KPIs de eficiencia para comparacion in vivo

Los KPIs preparados para A12 separan eficiencia fermentativa, eficiencia operacional, calidad quimica y retencion aromatica.

@@TABLE_KPI@@

La comparacion final debe reportarse por tratamiento, lote y replica:

$$
\Delta KPI = KPI_{optimal}-KPI_{commercial}
$$

con normalizacion multicriterio:

$$
z_j(d)=\frac{KPI_j(d)-KPI_j^{min}}{KPI_j^{max}-KPI_j^{min}}
$$

y score agregado:

$$
Score(d)=\sum_j w_j z_j(d)
$$

Los pesos deben quedar definidos antes de mirar resultados finales para evitar sesgo posterior.

## 10. Bitacoras, incidencias y datos crudos

Este bundle incluye plantillas en `templates/`:

- `execution_log_template.csv`
- `qaqc_plan.csv`
- `kpi_schema.csv`

El inventario de fuentes queda en:

@@TABLE_INVENTORY@@

Las fotos, cromatogramas/reportes analiticos originales, bitacoras firmadas y datos crudos de la nueva campana no se incluyen porque aun no estan cerrados en el repositorio actual. Deben anexarse en la version final posterior a ejecucion.

## 11. Interpretacion y decision tecnica

La decision recomendada es **no cerrar todavia protocolos optimos realistas para ejecucion comparativa final**. La evidencia justifica postergar la actividad in vivo final hasta completar:

1. Ingestion de datos reales del primer lote.
2. Recalibracion del modelo con tiempos, pulsos y condiciones iniciales ejecutadas.
3. Recalculo de FIM/perfiles para identificar parametros aun debiles.
4. Recalculo de politicas MPCC sobre el modelo calibrado.
5. Seleccion multicriterio congelada antes de ejecutar comparacion optimo vs comercial.

Esta decision evita que la optimizacion MPCC produzca artefactos derivados de parametros no identificables, estados secundarios mal ajustados o una capa aromatica aun no validada in vivo.

## 12. Conclusiones

- El Gemelo Digital ya cuenta con una base ODE calibrada/validada, diagnostico de estimabilidad y diseno MBDoE operacionalizable.
- El pipeline A09 demuestra capacidad de optimizacion dinamica dFBA/KKT/MPCC/MCDM, pero aun bajo modelo nominal.
- La campana propuesta de 9 fermentaciones es defendible como etapa de robustecimiento previo a A12 final.
- La actividad 5.7 queda metodologicamente cubierta por el plan de aromas, CO2 y condensado, aunque los resultados analiticos finales siguen pendientes.
- La actividad 5.8 queda cubierta como marco KPI y preparacion de comparacion; el resultado cuantitativo in vivo requiere ejecutar y cerrar la campana.

## 13. Contenido del bundle

- `A12_report.md`: este reporte principal.
- `README_BUNDLE.md`: guia de lectura.
- `CLAIMS_Y_LIMITACIONES_A12.md`: claims defendibles y no defendibles.
- `figures/`: figuras compactas de validacion, DOE, aromas y MPCC.
- `tables/`: tablas agregadas y schemas.
- `templates/`: bitacora, QA/QC y KPIs.
- `source_trace/`: documentos fuente A08/A09 y reportes tecnicos seleccionados.
"""
    replacements = {
        "@@DATE@@": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "@@BUNDLE@@": BUNDLE_NAME,
        "@@TABLE_ACTIVITY_STATUS@@": md_table(tables["activity_status"]),
        "@@TABLE_FIT@@": md_table(fit_short),
        "@@TABLE_VALIDATION@@": md_table(validation_short),
        "@@TABLE_PROFILE@@": md_table(profile_short),
        "@@TABLE_A09_KEY@@": md_table(a09_key_short),
        "@@TABLE_MPCC@@": md_table(mpcc),
        "@@TABLE_PROTOCOL@@": md_table(protocol_short),
        "@@TABLE_SELECTED@@": md_table(selected_short),
        "@@TABLE_POLICY@@": md_table(policy),
        "@@TABLE_STOCKS@@": md_table(tables["stock_assumptions"]),
        "@@TABLE_AROMA_SELECTED@@": md_table(aroma_selected),
        "@@TABLE_AROMA_BENCH@@": md_table(aroma_bench),
        "@@TABLE_AROMA_REDUCTION@@": md_table(aroma_red),
        "@@TABLE_QAQC@@": md_table(tables["qaqc_plan"]),
        "@@TABLE_KPI@@": md_table(tables["kpi_schema"]),
        "@@TABLE_INVENTORY@@": md_table(tables["raw_data_inventory"]),
    }
    for marker, value in replacements.items():
        report = report.replace(marker, value)
    (BUNDLE_DIR / "A12_report.md").write_text(report, encoding="utf-8")


def build_claims_and_readme() -> None:
    claims = """# Claims y limitaciones A12

## Claims respaldados

- Se dispone de un pipeline de calibracion, validacion, estimabilidad y MBDoE para el modelo de fermentacion.
- Se dispone de una campana experimental operacionalizable de 9 fermentaciones, 3 lotes, con restriccion de volumen y muestreo.
- Se dispone de una formulacion MPCC/dFBA/KKT/MCDM preliminar documentada en A09.
- Se dispone de un plan de analisis de aromas en vino y condensado, con CO2 online y condensado terminal como cierre de masa parcial.
- Se dispone de schemas de KPI, QA/QC y bitacora para ejecutar A12 de forma trazable.

## Claims parciales

- Los protocolos optimos existen como candidatos metodologicos o preliminares, no como protocolos realistas finales.
- La comparacion optimo vs comercial esta disenada conceptualmente, pero no cerrada con resultados in vivo.
- La capa aromatica permite disenar experimentos y definir mediciones, pero aun no entrega estadistica final por tratamiento.

## Claims que no deben hacerse

- No afirmar que las vinificaciones A12 finales ya compararon protocolo optimo versus comercial.
- No afirmar que el candidato MPCC A09 es optimo global o recomendacion industrial final.
- No afirmar que los marcadores aromaticos en vino/condensado ya fueron estadisticamente comparados entre protocolos finales.
- No afirmar eficiencia in vivo superior hasta tener datos ejecutados, QA/QC aprobado y analisis estadistico.

## Mensaje defendible

El avance A12 queda tecnicamente respaldado como transicion desde optimizacion nominal hacia validacion experimental robusta. La postergacion del cierre in vivo final se justifica porque ejecutar protocolos optimos no recalibrados podria validar artefactos del modelo, no mejoras reales de proceso.
"""
    (BUNDLE_DIR / "CLAIMS_Y_LIMITACIONES_A12.md").write_text(claims, encoding="utf-8")

    readme = f"""# Bundle A12 - Vinificaciones, aromas y KPIs

Fecha: {datetime.now().strftime("%Y-%m-%d %H:%M")}

Este bundle compacta evidencia para el Anexo A12:

- Actividad 5.6: Vinificaciones con protocolos optimos y comerciales.
- Actividad 5.7: Analisis de marcadores aromaticos en vino y condensados.
- Actividad 5.8: Indicadores clave de eficiencia para politicas optimas y comerciales in vivo.

## Lectura recomendada

1. `A12_report.md`: reporte principal autocontenido.
2. `CLAIMS_Y_LIMITACIONES_A12.md`: que se puede declarar y que no.
3. `figures/fig01_a12_activity_status.png`: estado de avance por actividad.
4. `figures/fig02_model_validation_predicted_observed.png`: evidencia de validacion del modelo.
5. `figures/fig05_a09_preliminary_mpcc_tradeoff.png`: conexion con MPCC/A09.
6. `figures/fig06_aroma_parameter_reduction.png`: evidencia de diseno de aromas.

## Evidencia principal

- `tables/campaign_protocol.csv`: lotes y tratamientos propuestos.
- `tables/selected_campaign.csv`: campana MBDoE seleccionada.
- `tables/volume_policy.csv` y `tables/volume_audit.csv`: factibilidad operacional.
- `tables/aroma_selected.csv`: campana aromatica MBDoE.
- `tables/kpi_schema.csv`: KPIs de evaluacion in vivo.
- `templates/`: bitacora, QA/QC y schema de KPIs.

## Complementario

- `source_trace/`: reportes fuente A08/A09 y documentos metodologicos seleccionados.
- `figures/source_*`: figuras fuente copiadas desde A09 o diagnosticos del modelo.
- `tables/a09_*`: tablas resumidas del trabajo MPCC/KKT/MCDM.

El archivo zip queda junto a la carpeta del bundle.
"""
    (BUNDLE_DIR / "README_BUNDLE.md").write_text(readme, encoding="utf-8")


def build_manifest_and_zip() -> Path:
    rows = []
    for path in sorted(BUNDLE_DIR.rglob("*")):
        if path.is_file():
            rows.append({"relative_path": str(path.relative_to(BUNDLE_DIR)).replace("\\", "/"), "bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(BUNDLE_DIR / "manifest.csv", index=False)

    zip_path = RESULTS_DIR / f"{BUNDLE_NAME}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=8) as zf:
        for path in sorted(BUNDLE_DIR.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=str(Path(BUNDLE_NAME) / path.relative_to(BUNDLE_DIR)))

    PROJECT_A12_DIR.mkdir(parents=True, exist_ok=True)
    copy_if_exists(zip_path, PROJECT_A12_DIR / zip_path.name)
    copy_if_exists(BUNDLE_DIR / "A12_report.md", PROJECT_A12_DIR / "A12_report.md")
    copy_if_exists(BUNDLE_DIR / "README_BUNDLE.md", PROJECT_A12_DIR / "README_BUNDLE.md")
    return zip_path


def main() -> None:
    prepare_dirs()
    tables = build_compact_tables()
    make_figures(tables)
    copy_source_trace()
    write_templates()
    build_report(tables)
    build_claims_and_readme()
    zip_path = build_manifest_and_zip()
    print(f"[done] bundle: {BUNDLE_DIR}")
    print(f"[done] zip: {zip_path}")
    print(f"[done] project copy: {PROJECT_A12_DIR}")


if __name__ == "__main__":
    main()
