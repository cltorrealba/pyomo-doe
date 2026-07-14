from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_new_must_glycerol_estimability_doe as model

BASE_DIR = SCRIPT_DIR / "results" / "new_must_glycerol_estimability_doe"
OVERNIGHT_DIR = SCRIPT_DIR / "results" / "new_must_glycerol_overnight_validation"
POST_DIR = OVERNIGHT_DIR / "best_multistart_post_analysis"
OUT_DIR = SCRIPT_DIR / "results" / "curve_validation"
PLOT_DIR = OUT_DIR / "plots"

STATE_UNITS = {
    "X": "kg/m3",
    "Xd": "kg/m3",
    "N": "kg/m3",
    "G": "g/L",
    "F": "g/L",
    "E": "g/L",
    "Gly": "g/L",
}
FEATURE_SCALES = {
    "X": 1.0,
    "Xd": 0.5,
    "N": 0.10,
    "G": 50.0,
    "F": 50.0,
    "E": 50.0,
    "Gly": 4.0,
}


def load_theta(path: Path) -> dict[str, float]:
    values = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
    return model.clip_theta({str(name): float(value) for name, value in values.items()})


def load_base_and_best_theta() -> tuple[dict[str, float], dict[str, float], str]:
    base_theta = load_theta(BASE_DIR / "theta_mixed_full17_l2.csv")
    multistart = pd.read_csv(OVERNIGHT_DIR / "multistart_summary.csv").sort_values("final_wsse")
    best_label = str(multistart.iloc[0]["fit"])
    best_theta = load_theta(OVERNIGHT_DIR / f"theta_{best_label}.csv")
    return base_theta, best_theta, best_label


def dense_time_for_batch(batch: model.BatchData, n: int = 320) -> np.ndarray:
    dense = np.linspace(float(batch.time[0]), float(batch.time[-1]), int(n))
    return np.array(sorted(set(np.round(np.concatenate([dense, batch.time]), 8))), dtype=float)


def dense_time_for_design(design: model.FutureDesign, sample_policy: str = "balanced", n: int = 360) -> np.ndarray:
    dense = np.linspace(0.0, float(design.horizon_h), int(n))
    samples = model.operational_sample_times(design.horizon_h, policy=sample_policy)
    pulse_times = []
    for schedule in design.pulses.values():
        pulse_times.extend([float(t) for t, amount in schedule if float(amount) > 0.0])
    return np.array(sorted(set(np.round(np.concatenate([dense, samples, np.array(pulse_times, dtype=float)]), 8))), dtype=float)


def write_fit_comparison_plots(
    batches: list[model.BatchData],
    base_theta: dict[str, float],
    best_theta: dict[str, float],
    best_label: str,
) -> pd.DataFrame:
    rows = []
    fit_dir = PLOT_DIR / "current_fit"
    fit_dir.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        time = dense_time_for_batch(batch)
        sim_base = model.simulate(batch, base_theta, time)
        sim_best = model.simulate(batch, best_theta, time)
        if sim_base is None or sim_best is None:
            rows.append({"medium": batch.medium, "batch": batch.batch, "status": "simulation_failed"})
            continue
        fig, axes = plt.subplots(4, 2, figsize=(12, 10), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes, model.STATE_NAMES):
            obs = np.asarray(batch.observations[state], dtype=float)
            mask = np.isfinite(obs)
            ax.plot(sim_base.index, sim_base[state], color="tab:blue", lw=2.0, label="base theta")
            ax.plot(sim_best.index, sim_best[state], color="tab:orange", lw=2.0, linestyle="--", label=best_label)
            ax.scatter(batch.time[mask], obs[mask], color="black", s=20, zorder=4, label="data")
            ax.set_title(f"{state} [{STATE_UNITS[state]}]")
            ax.grid(True, alpha=0.25)
            if mask.any():
                pred_base = np.interp(batch.time[mask], sim_base.index.to_numpy(dtype=float), sim_base[state].to_numpy(dtype=float))
                pred_best = np.interp(batch.time[mask], sim_best.index.to_numpy(dtype=float), sim_best[state].to_numpy(dtype=float))
                sigma = model.sigma_for_state(state, obs[mask])
                rows.append(
                    {
                        "medium": batch.medium,
                        "batch": batch.batch,
                        "state": state,
                        "n": int(mask.sum()),
                        "base_weighted_rmse": float(np.sqrt(np.mean(((pred_base - obs[mask]) / sigma) ** 2))),
                        "best_weighted_rmse": float(np.sqrt(np.mean(((pred_best - obs[mask]) / sigma) ** 2))),
                        "base_rmse": float(np.sqrt(np.mean((pred_base - obs[mask]) ** 2))),
                        "best_rmse": float(np.sqrt(np.mean((pred_best - obs[mask]) ** 2))),
                        "best_minus_base_weighted_rmse": float(
                            np.sqrt(np.mean(((pred_best - obs[mask]) / sigma) ** 2))
                            - np.sqrt(np.mean(((pred_base - obs[mask]) / sigma) ** 2))
                        ),
                        "status": "ok",
                    }
                )
        axes[-1].axis("off")
        axes[0].legend(loc="best", fontsize=8)
        fig.suptitle(f"{batch.medium}/{batch.batch}: current-data fit comparison", y=0.995)
        fig.tight_layout()
        fig.savefig(fit_dir / f"fitcmp_{batch.medium}_{batch.batch}.png", dpi=160)
        plt.close(fig)
    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "current_fit_curve_metrics.csv", index=False)
    if not table.empty and "state" in table.columns:
        summary = (
            table[table["status"].eq("ok")]
            .groupby(["medium", "state"], as_index=False)
            .agg(
                n=("n", "sum"),
                base_weighted_rmse=("base_weighted_rmse", "mean"),
                best_weighted_rmse=("best_weighted_rmse", "mean"),
                best_minus_base_weighted_rmse=("best_minus_base_weighted_rmse", "mean"),
            )
        )
        summary.to_csv(OUT_DIR / "current_fit_curve_metrics_by_medium_state.csv", index=False)
    return table


def selected_design_names() -> list[str]:
    names: list[str] = []
    for path in [
        BASE_DIR / "selected_campaign_hybrid.csv",
        OVERNIGHT_DIR / "selected_campaign_balanced.csv",
        POST_DIR / "sel_best.csv",
    ]:
        if path.exists():
            selected = pd.read_csv(path)
            for name in selected["candidate"].dropna().astype(str):
                if name not in names:
                    names.append(name)
    return names


def pulse_times(design: model.FutureDesign) -> list[tuple[str, float, float]]:
    rows = []
    for channel, schedule in design.pulses.items():
        for time_h, amount in schedule:
            if float(amount) > 0.0:
                rows.append((channel, float(time_h), float(amount)))
    return rows


def feasibility_issues(sim: pd.DataFrame, design: model.FutureDesign) -> list[str]:
    issues = []
    if float(sim["E"].max()) > 150.0:
        issues.append("ethanol_above_150_g_L")
    if float(sim["X"].max()) > 6.0:
        issues.append("viable_biomass_above_6_kg_m3")
    if float(sim["Xd"].max()) > 4.0:
        issues.append("dead_biomass_above_4_kg_m3")
    if float(sim["Gly"].max()) > 15.0:
        issues.append("glycerol_above_15_g_L")
    if float((sim["G"] + sim["F"]).iloc[-1]) > 35.0:
        issues.append("high_final_residual_sugar")
    if float(sim["G"].min()) <= 1e-6 and float(sim["F"].min()) <= 1e-6 and float(sim["E"].iloc[-1]) < 50.0:
        issues.append("low_ethanol_after_sugar_depletion")
    for state in model.STATE_NAMES:
        if not np.all(np.isfinite(sim[state].to_numpy(dtype=float))):
            issues.append(f"{state}_nonfinite")
    return issues


def write_design_prediction_plots(
    designs: dict[str, model.FutureDesign],
    names: list[str],
    base_theta: dict[str, float],
    best_theta: dict[str, float],
    best_label: str,
    sample_policy: str = "balanced",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pred_dir = PLOT_DIR / "design_predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    disagreement_rows = []
    for name in names:
        design = designs[name]
        time = dense_time_for_design(design, sample_policy=sample_policy)
        sim_base = model.simulate(design, base_theta, time)
        sim_best = model.simulate(design, best_theta, time)
        if sim_base is None or sim_best is None:
            rows.append({"candidate": name, "status": "simulation_failed"})
            continue
        samples = model.operational_sample_times(design.horizon_h, policy=sample_policy)
        pulses = pulse_times(design)
        for theta_label, sim in [("base", sim_base), (best_label, sim_best)]:
            issues = feasibility_issues(sim, design)
            rows.append(
                {
                    "candidate": name,
                    "theta": theta_label,
                    "medium": design.medium,
                    "family": design.family,
                    "horizon_h": design.horizon_h,
                    "final_G": float(sim["G"].iloc[-1]),
                    "final_F": float(sim["F"].iloc[-1]),
                    "final_sugar_GF": float((sim["G"] + sim["F"]).iloc[-1]),
                    "final_E": float(sim["E"].iloc[-1]),
                    "max_E": float(sim["E"].max()),
                    "max_X": float(sim["X"].max()),
                    "max_Xd": float(sim["Xd"].max()),
                    "max_Gly": float(sim["Gly"].max()),
                    "min_N": float(sim["N"].min()),
                    "issue_count": len(issues),
                    "issues": "; ".join(issues),
                    "status": "ok",
                }
            )
        pieces = []
        for state in model.STATE_NAMES:
            base_values = sim_base[state].to_numpy(dtype=float)
            best_values = sim_best[state].to_numpy(dtype=float)
            sigma = model.sigma_for_state(state, 0.5 * (base_values + best_values))
            pieces.append(((best_values - base_values) / sigma) ** 2)
        disagreement_rows.append(
            {
                "candidate": name,
                "medium": design.medium,
                "family": design.family,
                "theta_disagreement_weighted_rmse": float(np.sqrt(np.mean(np.concatenate(pieces)))),
            }
        )
        fig, axes = plt.subplots(4, 2, figsize=(12, 10), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes, model.STATE_NAMES):
            ax.plot(sim_base.index, sim_base[state], color="tab:blue", lw=2.0, label="base theta")
            ax.plot(sim_best.index, sim_best[state], color="tab:orange", lw=2.0, linestyle="--", label=best_label)
            for sample in samples:
                ax.axvline(sample, color="0.85", lw=0.6, zorder=0)
            for channel, time_h, _amount in pulses:
                ax.axvline(time_h, color="tab:red", lw=0.9, linestyle=":", alpha=0.7)
            ax.set_title(f"{state} [{STATE_UNITS[state]}]")
            ax.grid(True, alpha=0.25)
        axes[-1].axis("off")
        axes[0].legend(loc="best", fontsize=8)
        pulse_text = ", ".join(f"{ch}@{t:.0f}h" for ch, t, _amount in pulses) or "no pulses"
        fig.suptitle(f"{name}: predicted trajectories ({sample_policy}); {pulse_text}", y=0.995, fontsize=11)
        fig.tight_layout()
        fig.savefig(pred_dir / f"pred_{name}.png", dpi=160)
        plt.close(fig)
    feasibility = pd.DataFrame(rows)
    disagreement = pd.DataFrame(disagreement_rows).sort_values("theta_disagreement_weighted_rmse", ascending=False)
    feasibility.to_csv(OUT_DIR / "design_prediction_feasibility.csv", index=False)
    disagreement.to_csv(OUT_DIR / "design_theta_disagreement.csv", index=False)
    return feasibility, disagreement


def design_feature_vector(design: model.FutureDesign, theta: dict[str, float]) -> np.ndarray | None:
    fractions = np.linspace(0.0, 1.0, 21)
    time = fractions * float(design.horizon_h)
    sim = model.simulate(design, theta, time)
    if sim is None:
        return None
    pieces = []
    for state in model.STATE_NAMES:
        pieces.append(sim[state].to_numpy(dtype=float) / FEATURE_SCALES[state])
    return np.concatenate(pieces)


def write_redundancy_analysis(designs: dict[str, model.FutureDesign], names: list[str], theta: dict[str, float]) -> pd.DataFrame:
    vectors = {}
    for name in names:
        vector = design_feature_vector(designs[name], theta)
        if vector is not None and np.all(np.isfinite(vector)):
            vectors[name] = vector
    rows = []
    for idx, name_a in enumerate(vectors):
        for name_b in list(vectors)[idx + 1 :]:
            va = vectors[name_a]
            vb = vectors[name_b]
            distance = float(np.linalg.norm(va - vb) / math.sqrt(len(va)))
            corr = float(np.corrcoef(va, vb)[0, 1]) if np.std(va) > 0 and np.std(vb) > 0 else np.nan
            rows.append(
                {
                    "candidate_a": name_a,
                    "candidate_b": name_b,
                    "scaled_rms_distance": distance,
                    "trajectory_correlation": corr,
                }
            )
    table = pd.DataFrame(rows).sort_values("scaled_rms_distance")
    table.to_csv(OUT_DIR / "design_pairwise_similarity.csv", index=False)
    return table


def write_report(
    fit_metrics: pd.DataFrame,
    feasibility: pd.DataFrame,
    disagreement: pd.DataFrame,
    similarity: pd.DataFrame,
) -> None:
    fit_summary = pd.DataFrame()
    if not fit_metrics.empty and "state" in fit_metrics.columns:
        fit_summary = (
            fit_metrics[fit_metrics["status"].eq("ok")]
            .groupby(["medium", "state"], as_index=False)
            .agg(
                base_weighted_rmse=("base_weighted_rmse", "mean"),
                best_weighted_rmse=("best_weighted_rmse", "mean"),
                best_minus_base_weighted_rmse=("best_minus_base_weighted_rmse", "mean"),
            )
            .sort_values(["medium", "state"])
        )
    feasibility_issues = feasibility[feasibility["issue_count"].fillna(0).astype(float) > 0].copy()
    lines = [
        "# Curve-level validation for the new-must glycerol DOE",
        "",
        "## Purpose",
        "",
        "This validation checks whether the calibration and DOE conclusions survive curve-level inspection. It compares the base theta against the best overnight multistart theta, simulates the selected DOE candidates, checks feasibility flags, and quantifies trajectory redundancy.",
        "",
        "## Current-data fit comparison",
        "",
        fit_summary.to_markdown(index=False) if not fit_summary.empty else "_No fit summary available._",
        "",
        "Negative `best_minus_base_weighted_rmse` means the best multistart improves that state/medium relative to the base theta.",
        "",
        "## Predictive feasibility flags",
        "",
        feasibility_issues[
            [
                "candidate",
                "theta",
                "medium",
                "family",
                "final_sugar_GF",
                "final_E",
                "max_X",
                "max_Xd",
                "max_Gly",
                "issues",
            ]
        ].to_markdown(index=False)
        if not feasibility_issues.empty
        else "No hard feasibility flags were triggered by the selected candidate simulations.",
        "",
        "## Base-vs-best theta prediction disagreement",
        "",
        disagreement.head(20).to_markdown(index=False) if not disagreement.empty else "_No disagreement table available._",
        "",
        "High disagreement means the experiment is informative but also sensitive to the current local optimum.",
        "",
        "## Most similar selected designs",
        "",
        similarity.head(20).to_markdown(index=False) if not similarity.empty else "_No pairwise similarity table available._",
        "",
        "Low distance and high correlation indicate potential redundancy.",
        "",
        "## Interpretation",
        "",
        "- The design set should be kept only if the predictive curves are plausible and the first block contains at least one natural-must transfer test.",
        "- If an experiment triggers feasibility flags under either theta, treat it as a candidate for manual review rather than automatic execution.",
        "- If two candidates have very low pairwise distance, choose the one that better targets the weak profile-likelihood direction or is easier operationally.",
    ]
    (OUT_DIR / "curve_validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    data = model.load_normalized_data()
    batches = model.make_batches(data)
    base_theta, best_theta, best_label = load_base_and_best_theta()
    designs = {design.name: design for design in model.make_future_designs(data)}
    names = selected_design_names()
    print(f"[fit curves] {len(batches)} batches", flush=True)
    fit_metrics = write_fit_comparison_plots(batches, base_theta, best_theta, best_label)
    print(f"[design predictions] {len(names)} candidates", flush=True)
    feasibility, disagreement = write_design_prediction_plots(designs, names, base_theta, best_theta, best_label, sample_policy="balanced")
    print("[similarity]", flush=True)
    similarity = write_redundancy_analysis(designs, names, best_theta)
    write_report(fit_metrics, feasibility, disagreement, similarity)
    print(f"[done] {OUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
