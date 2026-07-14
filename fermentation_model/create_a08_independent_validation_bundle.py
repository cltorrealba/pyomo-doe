from __future__ import annotations

import math
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
BUNDLE_DIR = RESULTS_DIR / "A08_independent_validation_bundle_2026-06-29"
FIG_DIR = BUNDLE_DIR / "figures"
TABLE_DIR = BUNDLE_DIR / "tables"
SRC_DIR = BUNDLE_DIR / "source_trace"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_new_must_glycerol_estimability_doe as core


def _read_indexed_theta(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, index_col=0)
    if frame.shape[1] != 1:
        raise ValueError(f"Expected one theta column in {path}")
    return {str(k): float(v) for k, v in frame.iloc[:, 0].to_dict().items()}


def _safe_read(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _format_float(value: float, ndigits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    value = float(value)
    if abs(value) >= 10000 or (abs(value) > 0 and abs(value) < 0.001):
        return f"{value:.{ndigits}e}"
    return f"{value:.{ndigits}f}"


def _markdown_table(frame: pd.DataFrame, max_rows: int | None = None) -> str:
    if frame.empty:
        return "_No data available._"
    out = frame.copy()
    if max_rows is not None:
        out = out.head(max_rows)
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]) or pd.api.types.is_numeric_dtype(out[col]):
            out[col] = out[col].map(lambda x: _format_float(x) if pd.notna(x) else "NA")
    return out.to_markdown(index=False)


def prepare_dirs() -> None:
    for directory in (BUNDLE_DIR, FIG_DIR, TABLE_DIR, SRC_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def build_prediction_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    data = core.load_normalized_data()
    batches = [b for b in core.make_batches(data) if b.medium == "synthetic"]
    theta_reference = _read_indexed_theta(RESULTS_DIR / "new_must_glycerol_overnight_validation" / "theta_reference.csv")
    theta_best = _read_indexed_theta(RESULTS_DIR / "new_must_glycerol_overnight_validation" / "theta_multistart_07.csv")
    theta_map = {
        "reference_mixed_full17_l2": theta_reference,
        "best_multistart_07": theta_best,
    }
    rows: list[dict[str, object]] = []
    for batch in batches:
        for theta_label, theta in theta_map.items():
            sim = core.simulate(batch, theta, batch.time)
            if sim is None:
                continue
            for state in core.STATE_NAMES:
                obs = np.asarray(batch.observations[state], dtype=float)
                pred = sim.loc[batch.time, state].to_numpy(dtype=float)
                mask = np.isfinite(obs)
                if not mask.any():
                    continue
                sigma = core.sigma_for_state(state, obs[mask])
                for t, observed, predicted, sig in zip(batch.time[mask], obs[mask], pred[mask], sigma):
                    residual = float(predicted - observed)
                    rows.append(
                        {
                            "set": "synthetic_independent_source",
                            "medium": batch.medium,
                            "batch": batch.batch,
                            "theta": theta_label,
                            "time_h": float(t),
                            "state": state,
                            "observed": float(observed),
                            "predicted": float(predicted),
                            "residual": residual,
                            "sigma": float(sig),
                            "standardized_residual": residual / float(sig),
                        }
                    )
    pred = pd.DataFrame(rows)
    metrics_rows: list[dict[str, object]] = []
    for (theta_label, state), group in pred.groupby(["theta", "state"], sort=True):
        obs = group["observed"].to_numpy(dtype=float)
        p = group["predicted"].to_numpy(dtype=float)
        residual = p - obs
        ss_res = float(np.sum(residual**2))
        ss_tot = float(np.sum((obs - np.mean(obs)) ** 2))
        metrics_rows.append(
            {
                "theta": theta_label,
                "state": state,
                "n": int(len(group)),
                "rmse": math.sqrt(float(np.mean(residual**2))),
                "mae": float(np.mean(np.abs(residual))),
                "bias": float(np.mean(residual)),
                "weighted_rmse": math.sqrt(float(np.mean(group["standardized_residual"].to_numpy(dtype=float) ** 2))),
                "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
            }
        )
    metrics = pd.DataFrame(metrics_rows)
    return pred, metrics


def make_figures(pred: pd.DataFrame, metrics: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    batch_summary = _safe_read(RESULTS_DIR / "new_must_data_loading" / "new_must_batch_summary.csv")
    if not batch_summary.empty:
        fig, ax = plt.subplots(figsize=(10, 4.8))
        plot_data = batch_summary.copy()
        plot_data["label"] = plot_data["medium"] + "/" + plot_data["batch"]
        colors = np.where(plot_data["medium"].eq("synthetic"), "#1f77b4", "#7f7f7f")
        ax.bar(plot_data["label"], pd.to_numeric(plot_data["n_rows"], errors="coerce"), color=colors)
        ax.set_ylabel("Number of time rows")
        ax.set_title("Data coverage by fermentation")
        ax.tick_params(axis="x", rotation=70)
        ax.legend(handles=[
            plt.Rectangle((0, 0), 1, 1, color="#1f77b4", label="Synthetic"),
            plt.Rectangle((0, 0), 1, 1, color="#7f7f7f", label="Natural"),
        ])
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig01_data_coverage.png", dpi=180)
        plt.close(fig)

    best = pred[pred["theta"].eq("best_multistart_07")].copy()
    if not best.empty:
        states = list(core.STATE_NAMES)
        fig, axes = plt.subplots(2, 4, figsize=(14, 7))
        axes = axes.ravel()
        for ax, state in zip(axes, states):
            g = best[best["state"].eq(state)]
            if g.empty:
                ax.axis("off")
                continue
            x = g["observed"].to_numpy(dtype=float)
            y = g["predicted"].to_numpy(dtype=float)
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
        fig.suptitle("Predicted vs observed, synthetic independent source, best theta", y=1.02)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig02_predicted_observed_synthetic.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 4.8))
        ordered = states
        data = [best.loc[best["state"].eq(state), "standardized_residual"].to_numpy(dtype=float) for state in ordered]
        ax.boxplot(data, labels=ordered, showfliers=False)
        ax.axhline(0.0, color="black", lw=1)
        ax.axhline(2.0, color="#d62728", lw=1, ls="--")
        ax.axhline(-2.0, color="#d62728", lw=1, ls="--")
        ax.set_ylabel("Standardized residual")
        ax.set_title("Residual structure by state, synthetic independent source")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig04_standardized_residuals_synthetic.png", dpi=180)
        plt.close(fig)

    curve = _safe_read(RESULTS_DIR / "curve_validation" / "current_fit_curve_metrics_by_medium_state.csv")
    curve = curve[curve.get("medium", pd.Series(dtype=str)).eq("synthetic")] if not curve.empty else curve
    if not curve.empty:
        x = np.arange(len(curve))
        width = 0.35
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.bar(x - width / 2, pd.to_numeric(curve["base_weighted_rmse"]), width, label="Reference theta", color="#9ecae1")
        ax.bar(x + width / 2, pd.to_numeric(curve["best_weighted_rmse"]), width, label="Best multistart", color="#fdae6b")
        ax.set_xticks(x)
        ax.set_xticklabels(curve["state"], rotation=0)
        ax.set_ylabel("Weighted RMSE")
        ax.set_title("Synthetic-must curve fit: reference vs best multistart")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig03_synthetic_weighted_rmse_by_state.png", dpi=180)
        plt.close(fig)

    profile = _safe_read(RESULTS_DIR / "new_must_glycerol_overnight_validation" / "profile_full_summary.csv")
    if not profile.empty:
        profile["max_lr_stat"] = pd.to_numeric(profile["max_lr_stat"], errors="coerce")
        profile = profile.sort_values("max_lr_stat")
        colors = np.where(profile["profile_identifiable_95"].astype(str).str.lower().eq("true"), "#2ca02c", "#d62728")
        fig, ax = plt.subplots(figsize=(9, 6.5))
        ax.barh(profile["parameter"], profile["max_lr_stat"].clip(upper=500), color=colors)
        ax.axvline(3.841458820694124, color="black", lw=1, ls="--", label="95% chi-square threshold")
        ax.set_xlabel("Profile LR statistic, clipped at 500")
        ax.set_title("Profile-likelihood identifiability screen")
        ax.legend(loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig05_profile_identifiability.png", dpi=180)
        plt.close(fig)

    sampling = _safe_read(RESULTS_DIR / "new_must_glycerol_overnight_validation" / "sampling_policy_benchmark.csv")
    if not sampling.empty:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        policies = sampling["policy"].astype(str).to_list()
        axes[0].bar(policies, pd.to_numeric(sampling["final_logdet"]), color="#3182bd")
        axes[0].set_title("D-opt logdet")
        axes[1].bar(policies, pd.to_numeric(sampling["final_trace_inv"]), color="#756bb1")
        axes[1].set_title("A-opt trace(inv FIM)")
        axes[2].bar(policies, pd.to_numeric(sampling["weak_mean_var_reduction"]), color="#31a354")
        axes[2].set_title("Weak-dir variance reduction")
        for ax in axes:
            ax.tick_params(axis="x", rotation=25)
        fig.suptitle("Sampling policy benchmark")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig06_sampling_policy_benchmark.png", dpi=180)
        plt.close(fig)

    policy = _safe_read(RESULTS_DIR / "final_operational_doe_volume_constrained" / "policy_summary.csv")
    if not policy.empty:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        names = policy["policy"].astype(str).to_list()
        axes[0].bar(names, pd.to_numeric(policy["campaign_logdet"]), color="#3182bd")
        axes[0].set_title("Campaign logdet")
        axes[1].bar(names, pd.to_numeric(policy["campaign_trace_inv"]), color="#756bb1")
        axes[1].set_title("Campaign trace(inv)")
        axes[2].bar(names, pd.to_numeric(policy["new_param_mean_var_reduction"]), color="#31a354")
        axes[2].set_title("Mean variance reduction")
        for ax in axes:
            ax.tick_params(axis="x", rotation=25)
        fig.suptitle("Volume-constrained DOE policy benchmark")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig07_volume_constrained_doe_policy.png", dpi=180)
        plt.close(fig)

    doe_images = sorted((RESULTS_DIR / "vc_doe_nb_plots").glob("doe_*.png"))
    if doe_images:
        fig, axes = plt.subplots(3, 3, figsize=(12, 10))
        for ax, path in zip(axes.ravel(), doe_images[:9]):
            img = plt.imread(path)
            ax.imshow(img)
            ax.set_title(path.stem.replace("doe_", "DOE "))
            ax.axis("off")
        for ax in axes.ravel()[len(doe_images[:9]) :]:
            ax.axis("off")
        fig.suptitle("Selected campaign input and sampling designs", y=0.995)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "fig08_selected_campaign_designs_montage.png", dpi=180)
        plt.close(fig)


def copy_tables_and_sources(pred: pd.DataFrame, metrics: pd.DataFrame) -> None:
    pred.to_csv(TABLE_DIR / "synthetic_predicted_observed_long.csv", index=False)
    metrics.to_csv(TABLE_DIR / "synthetic_validation_metrics_by_state.csv", index=False)

    table_sources = {
        "data_batch_summary.csv": RESULTS_DIR / "new_must_data_loading" / "new_must_batch_summary.csv",
        "density_sugar_fit_by_medium.csv": RESULTS_DIR / "new_must_data_loading" / "density_sugar_fit_by_medium.csv",
        "fit_summary.csv": RESULTS_DIR / "new_must_glycerol_estimability_doe" / "fit_summary.csv",
        "theta_fit_table.csv": RESULTS_DIR / "new_must_glycerol_estimability_doe" / "theta_fit_table.csv",
        "curve_metrics_by_medium_state.csv": RESULTS_DIR / "curve_validation" / "current_fit_curve_metrics_by_medium_state.csv",
        "residual_by_state.csv": RESULTS_DIR / "new_must_glycerol_estimability_doe" / "residual_by_state.csv",
        "profile_likelihood_summary_full17.csv": RESULTS_DIR / "new_must_glycerol_overnight_validation" / "profile_full_summary.csv",
        "sampling_policy_benchmark.csv": RESULTS_DIR / "new_must_glycerol_overnight_validation" / "sampling_policy_benchmark.csv",
        "pyomo_doe_selected_summary.csv": RESULTS_DIR / "new_must_glycerol_overnight_validation" / "pyomo_selected_summary.csv",
        "volume_constrained_policy_summary.csv": RESULTS_DIR / "final_operational_doe_volume_constrained" / "policy_summary.csv",
        "volume_constrained_selected_campaign.csv": RESULTS_DIR / "final_operational_doe_volume_constrained" / "selected_campaign.csv",
        "volume_audit.csv": RESULTS_DIR / "final_operational_doe_volume_constrained" / "volume_audit.csv",
        "post_campaign_estimability.csv": RESULTS_DIR / "final_operational_doe_volume_constrained" / "post_campaign_estimability.csv",
        "secondary_v2_state_fit_metrics.csv": RESULTS_DIR / "secondary_v2_model_evaluation" / "state_fit_metrics.csv",
        "secondary_v2_fit_comparison.csv": RESULTS_DIR / "secondary_v2_model_evaluation" / "fit_comparison.csv",
        "lot1_pulse_timing_mbdoe_summary.csv": RESULTS_DIR / "lot1_pulse_timing_mbdoe" / "lot1_pulse_timing_mbdoe_summary.csv",
        "lot1_express_optimal_sampling_summary.csv": RESULTS_DIR / "lot1_express_optimal_sampling" / "express_optimal_sampling_summary.csv",
    }
    for name, src in table_sources.items():
        _copy_if_exists(src, TABLE_DIR / name)

    figure_sources = {
        "source_density_sugar_regression.png": RESULTS_DIR / "new_must_data_loading" / "density_sugar_regression.png",
        "source_batch_trajectory_synthetic_MS007.png": RESULTS_DIR / "new_must_data_loading" / "batch_trajectory_synthetic_MS007.png",
        "source_fitcmp_synthetic_MS007.png": RESULTS_DIR / "curve_validation" / "plots" / "current_fit" / "fitcmp_synthetic_MS007.png",
        "source_fitcmp_synthetic_MS009.png": RESULTS_DIR / "curve_validation" / "plots" / "current_fit" / "fitcmp_synthetic_MS009.png",
        "source_fitcmp_synthetic_MS016.png": RESULTS_DIR / "curve_validation" / "plots" / "current_fit" / "fitcmp_synthetic_MS016.png",
        "source_fim_eigen_spectra.png": RESULTS_DIR / "new_must_glycerol_estimability_doe" / "plots" / "fim_eigen_spectra.png",
    }
    for name, src in figure_sources.items():
        _copy_if_exists(src, FIG_DIR / name)

    source_reports = {
        "curve_validation_report.md": RESULTS_DIR / "curve_validation" / "curve_validation_report.md",
        "curve_validation_decision_summary.md": RESULTS_DIR / "curve_validation" / "curve_validation_decision_summary.md",
        "overnight_validation_report.md": RESULTS_DIR / "new_must_glycerol_overnight_validation" / "overnight_validation_report.md",
        "secondary_v2_model_report.md": RESULTS_DIR / "secondary_v2_model_evaluation" / "secondary_v2_model_report.md",
        "volume_constrained_doe_report.md": RESULTS_DIR / "final_operational_doe_volume_constrained" / "volume_constrained_doe_report.md",
        "lot1_pulse_timing_mbdoe_report.md": RESULTS_DIR / "lot1_pulse_timing_mbdoe" / "lot1_pulse_timing_mbdoe_report.md",
        "express_optimal_sampling_report.md": RESULTS_DIR / "lot1_express_optimal_sampling" / "express_optimal_sampling_report.md",
    }
    for name, src in source_reports.items():
        _copy_if_exists(src, SRC_DIR / name)


def build_report(pred: pd.DataFrame, metrics: pd.DataFrame) -> None:
    batch_summary = _safe_read(TABLE_DIR / "data_batch_summary.csv")
    fit_summary = _safe_read(TABLE_DIR / "fit_summary.csv")
    curve = _safe_read(TABLE_DIR / "curve_metrics_by_medium_state.csv")
    profile = _safe_read(TABLE_DIR / "profile_likelihood_summary_full17.csv")
    sampling = _safe_read(TABLE_DIR / "sampling_policy_benchmark.csv")
    doe_policy = _safe_read(TABLE_DIR / "volume_constrained_policy_summary.csv")
    selected_campaign = _safe_read(TABLE_DIR / "volume_constrained_selected_campaign.csv")
    secondary = _safe_read(TABLE_DIR / "secondary_v2_fit_comparison.csv")
    lot1 = _safe_read(TABLE_DIR / "lot1_pulse_timing_mbdoe_summary.csv")
    lot1_sampling = _safe_read(TABLE_DIR / "lot1_express_optimal_sampling_summary.csv")

    synthetic_batches = batch_summary[batch_summary.get("medium", pd.Series(dtype=str)).eq("synthetic")]
    synthetic_metrics = metrics[metrics["theta"].eq("best_multistart_07")].copy()
    if not synthetic_metrics.empty:
        synthetic_metrics = synthetic_metrics[["state", "n", "rmse", "mae", "bias", "weighted_rmse", "r2"]]
    synthetic_curve = curve[curve.get("medium", pd.Series(dtype=str)).eq("synthetic")] if not curve.empty else curve
    if not synthetic_curve.empty:
        synthetic_curve = synthetic_curve[["state", "n", "base_weighted_rmse", "best_weighted_rmse", "best_minus_base_weighted_rmse"]]
    fit_short = fit_summary[["fit", "n_batches", "mediums", "n_parameters", "final_wsse", "n_residuals", "wsse_per_residual", "l2_lambda"]] if not fit_summary.empty else fit_summary
    profile_short = profile[["parameter", "n_success", "max_lr_stat", "crosses_left_95", "crosses_right_95", "profile_identifiable_95"]] if not profile.empty else profile
    sampling_short = sampling[["policy", "final_logdet", "final_min_eigenvalue", "final_condition_number", "final_trace_inv", "weak_mean_var_reduction", "weak_worst_var_reduction"]] if not sampling.empty else sampling
    policy_short = doe_policy[["policy", "campaign_logdet", "campaign_min_relative_eigenvalue", "campaign_trace_inv", "new_param_mean_var_reduction", "new_param_worst_var_reduction"]] if not doe_policy.empty else doe_policy
    selected_short_cols = [c for c in ["rank", "candidate", "medium", "family", "campaign_step", "objective_score", "final_volume_ml"] if c in selected_campaign.columns]
    selected_short = selected_campaign[selected_short_cols] if selected_short_cols else selected_campaign.head(9)
    secondary_short = secondary[["model", "success", "final_wsse", "n_residuals", "wsse_per_residual"]] if not secondary.empty and "model" in secondary.columns else secondary

    synthetic_n_rows = int(pd.to_numeric(synthetic_batches.get("n_rows", pd.Series(dtype=float)), errors="coerce").sum()) if not synthetic_batches.empty else 0
    table_synthetic_coverage = _markdown_table(
        synthetic_batches[
            [
                "medium",
                "batch",
                "n_rows",
                "t_min_h",
                "t_max_h",
                "temperature_initial_c",
                "temperature_min_c",
                "temperature_max_c",
                "G_initial_g_l",
            ]
        ]
        if not synthetic_batches.empty
        else synthetic_batches
    )

    report = r"""# A08 - Validacion independiente del Gemelo Digital

**Actividad:** 3.15  
**Titulo:** Validacion independiente del Gemelo Digital  
**Fecha de generacion:** @@GEN_TIME@@  
**Bundle:** `@@BUNDLE_NAME@@`

## 1. Objetivo del anexo

Este anexo documenta la validacion independiente del Gemelo Digital de fermentacion, usando como evidencia principal un set de mosto sintetico (`MS007` a `MS016`) cargado desde las bases nuevas del proyecto. El objetivo fue verificar si el modelo ODE calibrado reproduce curvas observadas de biomasa viable, biomasa muerta, nitrogeno asimilable, glucosa, fructosa, etanol y glicerol, y si la estructura resultante es suficientemente robusta para avanzar a diseno experimental y posterior uso como restriccion dinamica en control predictivo.

La independencia se interpreta como independencia respecto de la base historica usada al inicio del desarrollo. El set sintetico fue luego usado para diagnostico, recalibracion y MBDoE; por tanto, no debe presentarse como holdout ciego definitivo. Para una validacion externa estricta, el primer lote experimental nuevo debe evaluarse como set prospectivo.

## 2. Evidencia fuente

El set sintetico contiene @@SYN_N_BATCHES@@ fermentaciones y @@SYN_N_ROWS@@ filas de tiempo-proceso. Las variables relevantes fueron normalizadas antes de modelar:

- `time_h`: tiempo de proceso.
- `temperature_c`: setpoint de temperatura.
- `X`: biomasa viable, convertida desde millones de celulas/mL usando 30 pg/cell.
- `Xd`: biomasa muerta, calculada como biomasa total Oculyze menos biomasa viable.
- `N`: YAN en kg/m3.
- `G`, `F`: glucosa y fructosa en g/L.
- `E`: etanol, convertido a g/L.
- `Gly`: glicerol en g/L.
- Pulsos nutricionales: tratados como entradas suaves en el balance de nitrogeno.

**Cobertura del set sintetico**

@@TABLE_SYNTHETIC_COVERAGE@@

Figura principal de cobertura: `figures/fig01_data_coverage.png`.  
Figura de limpieza de proxy densidad-azucar: `figures/source_density_sugar_regression.png`.

## 3. Modelo de simulacion

El modelo dinamico principal es un sistema ODE con estados:

$$
x(t)=\\left[X, X_d, N, G, F, E, Gly\\right]^T
$$

donde `X` es biomasa viable, `X_d` biomasa muerta, `N` nitrogeno asimilable, `G` glucosa, `F` fructosa, `E` etanol y `Gly` glicerol. La temperatura `T(t)` y los pulsos externos `u_j(t)` son entradas de proceso.

La dinamica usada para simulacion fue:

$$
\\frac{dX}{dt}=\\left(\\mu-k_d\\right)X+u_X(t)
$$

$$
\\frac{dX_d}{dt}=k_d X
$$

$$
\\frac{dN}{dt}=-q_N a_{\\mu}(T)\\frac{N}{N+K_N(T)}X+u_N(t)
$$

$$
\\frac{dG}{dt}=-\\left[q_{XG}a_{\\mu}(T)\\frac{N}{N+K_N(T)}
+q_{EG}a_{\\beta}(T)\\frac{G}{G+K_G(T)}I_E(E)
+m(T)\\frac{G}{G+F}\\right]X+u_G(t)
$$

$$
\\frac{dF}{dt}=-\\left[q_{XF}a_{\\mu}(T)\\frac{N}{N+K_N(T)}
+q_{EF}a_{\\beta}(T)\\frac{F}{F+K_F(T)}I_G(G)I_E(E)
+m(T)\\frac{F}{G+F}\\right]X+u_F(t)
$$

$$
\\frac{dE}{dt}=\\left[\\beta_G(T,G,E)+\\beta_F(T,F,G,E)\\right]X+u_E(t)
$$

$$
\\frac{dGly}{dt}=\\left[\\gamma_{G0}a_{\\beta}(T)\\frac{G}{G+K_G(T)}I_E(E)
+\\gamma_{F0}a_{\\beta}(T)\\frac{F}{F+K_F(T)}I_G(G)I_E(E)\\right]X
$$

con:

$$
\\mu=\\mu_0 a_{\\mu}(T)\\frac{N}{N+K_N(T)}
$$

$$
\\beta_G=\\beta_{G0}a_{\\beta}(T)\\frac{G}{G+K_G(T)}I_E(E)
$$

$$
\\beta_F=\\beta_{F0}a_{\\beta}(T)\\frac{F}{F+K_F(T)}I_G(G)I_E(E)
$$

$$
I_E(E)=\\frac{1}{1+i_E(T)E},\\qquad I_G(G)=\\frac{1}{1+i_G(T)G}
$$

Los factores de temperatura fueron parametrizados tipo Arrhenius usando constantes fijas del modelo historico. Los pulsos se representaron como gaussianas normalizadas:

$$
u_j(t)=\\sum_k \\Delta C_{j,k}\\frac{\\exp\\left[-\\left(\\frac{t-t_{j,k}}{w}\\right)^2\\right]}{\\sqrt{\\pi}w}
$$

El vector completo evaluado incluyo:

$$
\\theta=\\left[\\mu_0,s_N,q_N,q_{XG},q_{XF},\\beta_{G0},s_G,\\beta_{F0},s_F,q_{EG},q_{EF},i_G,i_E,K_{d0},m_0,\\gamma_{G0},\\gamma_{F0}\\right]
$$

## 4. Modelo de calibracion y validacion

La calibracion se formulo como minimos cuadrados ponderados en escala de observacion:

$$
\\min_{\\theta}\\; J(\\theta)=\\sum_{b\\in B}\\sum_{i\\in S}\\sum_{k\\in T_{b,i}}
\\left(\\frac{\\hat{y}_{b,i}(t_k;\\theta)-y_{b,i,k}}{\\sigma_i(y_{b,i,k})}\\right)^2
+\\lambda\\sum_{p\\in P_R}\\left[\\log\\left(\\frac{\\theta_p}{\\theta_{p,ref}}\\right)\\right]^2
$$

sujeto a:

$$
\\theta_p^{min}\\leq \\theta_p\\leq \\theta_p^{max}
$$

y a las ecuaciones ODE del modelo. La desviacion experimental se definio como:

$$
\\sigma_i(y)=\\max\\left(\\sigma_{i,floor}, r_i|y|\\right)
$$

Se evaluaron estrategias de ajuste natural, sintetico, mixto y mixto regularizado. La regularizacion L2 en log-parametros se uso para evitar soluciones numericamente mejores pero fisicamente poco transferibles en parametros debiles.

**Resumen de ajustes**

@@TABLE_FIT_SHORT@@

## 5. Predicho-observado y metricas de validacion sintetica

La evidencia principal de validacion se calculo simulando las fermentaciones sinteticas con el mejor conjunto `best_multistart_07` y comparando contra observaciones. La tabla resume metricas por estado.

@@TABLE_SYNTHETIC_METRICS@@

Figura predicho-observado: `figures/fig02_predicted_observed_synthetic.png`.  
Figura de residuos estandarizados: `figures/fig04_standardized_residuals_synthetic.png`.

La comparacion contra el theta de referencia muestra que el multistart mejora la mayoria de las curvas sinteticas, especialmente fructosa, etanol, glucosa y glicerol.

@@TABLE_SYNTHETIC_CURVE@@

Figura de benchmark de RMSE ponderado: `figures/fig03_synthetic_weighted_rmse_by_state.png`.

## 6. Analisis de residuos

Los residuos fueron evaluados como:

$$
r_{b,i,k}=\\frac{\\hat{y}_{b,i}(t_k;\\theta)-y_{b,i,k}}{\\sigma_i(y_{b,i,k})}
$$

El diagnostico muestra ajuste razonable para consumo de azucares, etanol y glicerol, con sesgos residuales todavia visibles en `N` y `X_d`. Esto es coherente con dos limitaciones observadas: el YAN tiende a presentar remanentes no consumidos bajo algunas condiciones sinteticas, y `X_d` depende de una resta entre biomasa total y viable, amplificando ruido experimental.

Los residuos punto a punto quedan disponibles en `tables/synthetic_predicted_observed_long.csv`, y las metricas agregadas en `tables/synthetic_validation_metrics_by_state.csv`.

## 7. Estimabilidad e identificabilidad practica

La estimabilidad se evaluo por aproximacion FIM, perfiles de likelihood y diagnostico bayesiano/Laplace. Para perfil likelihood se uso:

$$
LR_p(\\theta_p)=J(\\theta_p,\\hat{\\theta}_{-p}(\\theta_p))-J(\\hat{\\theta})
$$

con decision a 95% usando:

$$
LR_p>\\chi^2_{1,0.95}=3.841
$$

Un parametro se considero identificable por perfil si el perfil cruza el umbral a ambos lados del optimo.

**Resumen de perfiles**

@@TABLE_PROFILE_SHORT@@

Figura de perfiles agregados: `figures/fig05_profile_identifiability.png`.

Resultado principal: `qN`, `Kd0`, `gammaG0` y `gammaF0` aparecen identificables por perfil bajo el ajuste extendido. Parametros como `qXG`, `m0`, `betaG0` y varias direcciones unilaterales siguen siendo practicamente debiles o confusas con las condiciones experimentales disponibles.

## 8. MBDoE y campana experimental derivada

El diseno de experimentos se baso en la matriz de informacion de Fisher:

$$
F(\\theta,d)=\\sum_{k}\\left(\\frac{\\partial y(t_k;\\theta,d)}{\\partial \\theta}\\right)^T
W_k
\\left(\\frac{\\partial y(t_k;\\theta,d)}{\\partial \\theta}\\right)+F_{prior}
$$

El criterio principal fue D-optimalidad:

$$
\\Phi_D(d)=\\log\\det\\left(F(\\theta,d)+\\epsilon I\\right)
$$

y se reportaron criterios complementarios:

$$
\\Phi_A(d)=\\mathrm{tr}\\left(F^{-1}\\right)
$$

ademas de autovalores, condicionamiento y reduccion esperada de varianza en direcciones debiles.

**Benchmark de politica de muestreo**

@@TABLE_SAMPLING_SHORT@@

La politica `balanced` tuvo mayor `logdet`, menor `trace_inv` y mejor reduccion media de varianza en direcciones debiles. Posteriormente, el diseno fue auditado por restricciones operacionales y volumen de reactor.

**Benchmark DOE con restriccion de volumen**

@@TABLE_POLICY_SHORT@@

La politica seleccionada fue `full14_small6`, porque mantuvo el mejor compromiso entre informacion, restriccion de volumen y factibilidad de laboratorio.

**Campana seleccionada**

@@TABLE_SELECTED_SHORT@@

Figura de disenos seleccionados: `figures/fig08_selected_campaign_designs_montage.png`.

## 9. Extension secundaria y limitaciones

Se evaluo una capa secundaria para piruvato, acetaldehido, acetato y oxigeno disuelto. Esta capa no reemplaza el modelo principal del anexo, pero informa riesgos para etapas posteriores del Gemelo Digital y MPCC.

@@TABLE_SECONDARY_SHORT@@

La version `secondary_v2_reduced_o2fixed` mejoro respecto de la estructura secundaria inicial, pero algunos estados secundarios siguen mostrando error estructurado. Por ello, para control predictivo se recomienda usar esos estados con regularizacion/fijacion parcial y no exigir identificacion completa antes de disponer de nuevos datos con CO2, aromas y condensado.

## 10. Decision de avance

La decision tecnica es **avanzar con el Gemelo Digital calibrado como modelo operativo para diseno experimental y validacion prospectiva**, con las siguientes condiciones:

1. Usar el theta `best_multistart_07` como prior operativo alternativo y conservar el theta regularizado de referencia como comparador.
2. No declarar cerrada la identificabilidad practica de todos los parametros; actualizar perfiles y FIM tras el primer lote experimental.
3. Tratar `N`, `X_d`, `m0`, `qXG` y direcciones unilaterales como parametros/estados de riesgo.
4. Usar la campana DOE operacionalizada con politica `full14_small6`, registrando tiempos reales, volumen extraido, volumen de pulsos, CO2 online y condensado final.
5. Usar el primer lote nuevo como validacion prospectiva estricta.

## 11. Limitaciones

- El set sintetico es independiente respecto de la base historica inicial, pero fue incorporado al pipeline de recalibracion y diagnostico; no es holdout ciego final.
- El modelo no incluye aun una capa completa gas-liquido para CO2/aromas dentro de esta evidencia A08.
- La medicion de biomasa muerta deriva de diferencia entre biomasa total y viable, lo que aumenta incertidumbre.
- La representacion de nitrogeno total no separa completamente amonio/PAN ni remanentes no asimilables.
- La transferencia a mosto natural requiere validacion posterior, ya que micronutrientes y matriz pueden cambiar parametros aparentes.

## 12. Archivos principales del bundle

- `A08_report.md`: este reporte.
- `figures/fig02_predicted_observed_synthetic.png`: figura central predicho-observado.
- `figures/fig03_synthetic_weighted_rmse_by_state.png`: mejora referencia vs multistart.
- `figures/fig04_standardized_residuals_synthetic.png`: diagnostico de residuos.
- `figures/fig05_profile_identifiability.png`: parametros identificables/no identificables por perfil.
- `figures/fig06_sampling_policy_benchmark.png`: benchmark MBDoE de muestreo.
- `figures/fig08_selected_campaign_designs_montage.png`: campana DOE seleccionada.
- `tables/synthetic_validation_metrics_by_state.csv`: metricas principales.
- `tables/synthetic_predicted_observed_long.csv`: datos predicho-observado y residuos.
- `source_trace/`: reportes fuente usados para trazabilidad.
"""
    replacements = {
        "@@GEN_TIME@@": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "@@BUNDLE_NAME@@": BUNDLE_DIR.name,
        "@@SYN_N_BATCHES@@": str(len(synthetic_batches)),
        "@@SYN_N_ROWS@@": str(synthetic_n_rows),
        "@@TABLE_SYNTHETIC_COVERAGE@@": table_synthetic_coverage,
        "@@TABLE_FIT_SHORT@@": _markdown_table(fit_short),
        "@@TABLE_SYNTHETIC_METRICS@@": _markdown_table(synthetic_metrics),
        "@@TABLE_SYNTHETIC_CURVE@@": _markdown_table(synthetic_curve),
        "@@TABLE_PROFILE_SHORT@@": _markdown_table(profile_short),
        "@@TABLE_SAMPLING_SHORT@@": _markdown_table(sampling_short),
        "@@TABLE_POLICY_SHORT@@": _markdown_table(policy_short),
        "@@TABLE_SELECTED_SHORT@@": _markdown_table(selected_short),
        "@@TABLE_SECONDARY_SHORT@@": _markdown_table(secondary_short),
    }
    for marker, value in replacements.items():
        report = report.replace(marker, value)
    report = report.replace("\\\\", "\\")
    (BUNDLE_DIR / "A08_report.md").write_text(report, encoding="utf-8")


def build_readme() -> None:
    manifest_rows = []
    for path in sorted(BUNDLE_DIR.rglob("*")):
        if path.is_file():
            manifest_rows.append(
                {
                    "relative_path": str(path.relative_to(BUNDLE_DIR)).replace("\\", "/"),
                    "bytes": path.stat().st_size,
                }
            )
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(BUNDLE_DIR / "manifest.csv", index=False)

    readme = f"""# Bundle A08 - Validacion independiente del Gemelo Digital

Este bundle contiene la evidencia compacta para el anexo A08, Actividad 3.15.

## Archivo principal

- `A08_report.md`: reporte autocontenido con objetivo, set independiente, modelo matematico, calibracion, validacion predicho-observado, residuos, estimabilidad, MBDoE, limitaciones y decision de avance.

## Evidencia mas relevante

- `figures/fig02_predicted_observed_synthetic.png`: figura central de validacion predicho-observado para mosto sintetico.
- `figures/fig03_synthetic_weighted_rmse_by_state.png`: comparacion referencia vs mejor multistart.
- `figures/fig04_standardized_residuals_synthetic.png`: residuos estandarizados por estado.
- `figures/fig05_profile_identifiability.png`: lectura de identificabilidad practica por perfil likelihood.
- `figures/fig06_sampling_policy_benchmark.png`: criterio MBDoE para muestreo.
- `figures/fig07_volume_constrained_doe_policy.png`: benchmark operacional con restriccion de volumen.
- `figures/fig08_selected_campaign_designs_montage.png`: resumen visual de inputs/muestreos del diseno seleccionado.

## Tablas clave

- `tables/synthetic_validation_metrics_by_state.csv`: metricas por estado.
- `tables/synthetic_predicted_observed_long.csv`: predicho, observado y residuos.
- `tables/profile_likelihood_summary_full17.csv`: resultados de perfiles.
- `tables/volume_constrained_selected_campaign.csv`: campana experimental operacionalizada.

## Complementario

- `source_trace/`: copias de reportes fuente generados durante el pipeline.
- `manifest.csv`: listado de todos los archivos incluidos.

La version zip del bundle se genera junto a la carpeta: `{BUNDLE_DIR.name}.zip`.
"""
    (BUNDLE_DIR / "README.md").write_text(readme, encoding="utf-8")


def zip_bundle() -> Path:
    zip_path = RESULTS_DIR / f"{BUNDLE_DIR.name}.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=8) as zf:
        for path in sorted(BUNDLE_DIR.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=str(Path(BUNDLE_DIR.name) / path.relative_to(BUNDLE_DIR)))
    return zip_path


def main() -> None:
    prepare_dirs()
    pred, metrics = build_prediction_table()
    make_figures(pred, metrics)
    copy_tables_and_sources(pred, metrics)
    build_report(pred, metrics)
    build_readme()
    zip_path = zip_bundle()
    print(f"[done] bundle: {BUNDLE_DIR}")
    print(f"[done] zip: {zip_path}")


if __name__ == "__main__":
    main()
