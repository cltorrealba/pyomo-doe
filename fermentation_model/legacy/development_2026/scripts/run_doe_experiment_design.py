from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.contrib.doe import DesignOfExperiments

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_fit_strategy_analysis import load_notebook_context, series_to_theta

RESULTS_DIR = SCRIPT_DIR / "results" / "doe_experiment_design"
FINAL_THETA_PATH = SCRIPT_DIR / "results" / "identifiability_reduction" / "theta_final_identifiable.csv"

FINAL7 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF", "iG")
PRIMARY9 = FINAL7 + ("qXG", "qXF")
EXTENDED12 = PRIMARY9 + ("sN", "sG", "sF")
WEAK_PARAMETERS = ("qXG", "qXF", "iG")


@dataclass(frozen=True)
class DesignCandidate:
    name: str
    temperature_c: tuple[float, float, float, float]
    pulse_time_fractions: tuple[float, ...]
    pulse_amounts: tuple[float, ...]
    family: str


class FixedDesignExperiment:
    def __init__(
        self,
        ns: dict,
        batch_id: str,
        theta: dict[str, float],
        parameters: tuple[str, ...],
        temperature_c: tuple[float, ...],
        pulses: tuple[tuple[float, float], ...],
    ):
        self.ns = ns
        self.batch_id = str(batch_id)
        self.theta = dict(theta)
        self.parameters = list(parameters)
        self.temperature_c = tuple(float(v) for v in temperature_c)
        self.pulses = tuple((float(t), float(a)) for t, a in pulses)

    def get_labeled_model(self):
        ns = self.ns
        batch = ns["load_batch"](self.batch_id)
        experiment = ns["FermentationExperiment"](
            batch,
            theta_initial=self.theta,
            parameters_to_estimate=self.parameters,
            input_mode="design",
            temperature_segments=len(self.temperature_c),
            pulse_max_count=len(self.pulses),
            fix_design_inputs=False,
        )
        model = experiment.get_labeled_model()
        for segment, value in enumerate(self.temperature_c):
            model.T_set[segment].set_value(float(value))
        for slot, (time_h, amount) in enumerate(self.pulses):
            model.pulse_time[slot].set_value(float(time_h))
            model.pulse_amount[slot].set_value(float(amount))
        return model


class MeasuredExperiment:
    def __init__(self, ns: dict, batch_id: str, theta: dict[str, float], parameters: tuple[str, ...]):
        self.ns = ns
        self.batch_id = str(batch_id)
        self.theta = dict(theta)
        self.parameters = list(parameters)

    def get_labeled_model(self):
        ns = self.ns
        return ns["FermentationExperiment"](
            ns["load_batch"](self.batch_id),
            theta_initial=self.theta,
            parameters_to_estimate=self.parameters,
            input_mode="measured",
            fix_design_inputs=True,
        ).get_labeled_model()


def make_solver(ns: dict):
    solver = pyo.SolverFactory("ipopt", executable=ns["IPOPT_EXECUTABLE"])
    for key, value in ns["PARMEST_SOLVER_OPTIONS"].items():
        solver.options[key] = value
    solver.options["max_iter"] = max(int(solver.options.get("max_iter", 0) or 0), 3000)
    solver.options["tol"] = min(float(solver.options.get("tol", 1e-6)), 1e-6)
    return solver


def pulse_schedule_for_batch(ns: dict, batch_id: str, candidate: DesignCandidate, pulse_slots: int) -> tuple[tuple[float, float], ...]:
    batch = ns["load_batch"](batch_id)
    t0 = float(batch.time[0])
    tf = float(batch.time[-1])
    allowed = float(ns["DESIGN_PULSE_ALLOWED_FRACTION"]) * (tf - t0)
    fractions = list(candidate.pulse_time_fractions)[:pulse_slots]
    amounts = list(candidate.pulse_amounts)[:pulse_slots]
    while len(fractions) < pulse_slots:
        fractions.append(1.0)
    while len(amounts) < pulse_slots:
        amounts.append(0.0)
    schedule = []
    for frac, amount in zip(fractions, amounts):
        frac = min(max(float(frac), 0.0), 1.0)
        schedule.append((t0 + frac * allowed, min(max(float(amount), 0.0), float(ns["DESIGN_MAX_YAN_PER_PULSE"]))))
    schedule = sorted(schedule, key=lambda item: item[0])
    total = sum(amount for _time, amount in schedule)
    max_total = float(ns["DESIGN_MAX_TOTAL_YAN"])
    if total > max_total and total > 0:
        schedule = [(time, amount * max_total / total) for time, amount in schedule]
    return tuple(schedule)


def manual_candidates() -> list[DesignCandidate]:
    return [
        DesignCandidate("constant_20_split", (20, 20, 20, 20), (0.00, 0.25, 0.50, 0.75, 1.00), (0.04, 0.04, 0.04, 0.04, 0.04), "manual"),
        DesignCandidate("ramp_up_staggered", (15, 18, 22, 25), (0.05, 0.30, 0.60, 0.85), (0.05, 0.06, 0.06, 0.03), "manual"),
        DesignCandidate("ramp_down_frontload", (25, 22, 18, 15), (0.00, 0.15, 0.35, 0.65), (0.08, 0.06, 0.04, 0.02), "manual"),
        DesignCandidate("hot_start_cool_finish", (25, 25, 18, 15), (0.00, 0.20, 0.45), (0.10, 0.06, 0.04), "manual"),
        DesignCandidate("cool_start_hot_finish", (15, 15, 22, 25), (0.10, 0.35, 0.70), (0.06, 0.07, 0.06), "manual"),
        DesignCandidate("alternating_extremes", (25, 15, 25, 15), (0.00, 0.20, 0.45, 0.70), (0.05, 0.05, 0.05, 0.05), "manual"),
        DesignCandidate("late_identification", (18, 22, 25, 25), (0.25, 0.55, 0.90), (0.06, 0.07, 0.06), "manual"),
        DesignCandidate("early_biomass_boost", (23, 25, 20, 18), (0.00, 0.10, 0.25, 0.50), (0.08, 0.06, 0.04, 0.02), "manual"),
        DesignCandidate("low_temp_growth_separation", (15, 17, 20, 22), (0.00, 0.35, 0.70), (0.07, 0.07, 0.06), "manual"),
        DesignCandidate("high_temp_inhibition_probe", (25, 25, 25, 22), (0.00, 0.25, 0.50), (0.07, 0.07, 0.06), "manual"),
    ]


def random_candidates(n: int, seed: int) -> list[DesignCandidate]:
    rng = np.random.default_rng(seed)
    rows = []
    for idx in range(int(n)):
        if idx % 4 == 0:
            temps = tuple(float(v) for v in rng.choice([15, 18, 21, 25], size=4, replace=True))
        else:
            temps = tuple(float(v) for v in rng.uniform(15.0, 25.0, size=4))
        n_pulses = int(rng.integers(3, 6))
        fractions = tuple(sorted(float(v) for v in rng.uniform(0.0, 1.0, size=n_pulses)))
        total = float(rng.uniform(0.12, 0.20))
        weights = rng.dirichlet(np.ones(n_pulses))
        amounts = np.minimum(weights * total, 0.1)
        if amounts.sum() > 0:
            amounts *= min(total, 0.2) / amounts.sum()
        rows.append(
            DesignCandidate(
                f"lhs_{idx + 1:02d}",
                tuple(round(v, 3) for v in temps),
                tuple(round(v, 4) for v in fractions),
                tuple(round(float(v), 5) for v in amounts),
                "lhs",
            )
        )
    return rows


def fim_for_experiment(ns: dict, experiment, parameters: tuple[str, ...], step: float, solver) -> tuple[np.ndarray, np.ndarray]:
    doe = DesignOfExperiments(
        experiment=experiment,
        step=float(step),
        scale_nominal_param_value=True,
        solver=solver,
        tee=False,
    )
    fim = np.asarray(doe.compute_FIM(method="sequential"), dtype=float)
    jac = np.asarray(doe.seq_jac, dtype=float)
    if fim.shape != (len(parameters), len(parameters)):
        raise RuntimeError(f"Unexpected FIM shape {fim.shape}; expected {(len(parameters), len(parameters))}.")
    return 0.5 * (fim + fim.T), jac


def current_prior_fim(ns: dict, theta: dict[str, float], parameters: tuple[str, ...], batches: list[str], step: float, solver) -> np.ndarray:
    cache = RESULTS_DIR / f"prior_fim_{'_'.join(parameters)}_step{step:g}.csv"
    if cache.exists():
        return pd.read_csv(cache, index_col=0).loc[list(parameters), list(parameters)].to_numpy(dtype=float)
    fim_total = np.zeros((len(parameters), len(parameters)), dtype=float)
    rows = []
    for batch_id in batches:
        print(f"[DOE prior] batch={batch_id} parameters={len(parameters)}", flush=True)
        experiment = MeasuredExperiment(ns, batch_id, theta, parameters)
        fim, _jac = fim_for_experiment(ns, experiment, parameters, step, solver)
        fim_total += fim
        rows.append({"batch": batch_id, **fim_metrics(fim, parameters, prefix="batch_")})
    pd.DataFrame(fim_total, index=parameters, columns=parameters).to_csv(cache)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"prior_fim_{'_'.join(parameters)}_batch_metrics.csv", index=False)
    return fim_total


def stable_inverse(fim: np.ndarray, ridge_fraction: float = 1e-9) -> np.ndarray:
    fim = np.asarray(fim, dtype=float)
    scale = max(float(np.trace(fim)) / max(fim.shape[0], 1), 1.0)
    ridge = ridge_fraction * scale
    return np.linalg.pinv(fim + ridge * np.eye(fim.shape[0]))


def fim_metrics(fim: np.ndarray, parameters: tuple[str, ...], prefix: str = "") -> dict[str, float]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eig = np.linalg.eigvalsh(fim)
    max_eig = float(np.max(eig)) if eig.size else np.nan
    floor = max(max_eig * 1e-12, np.finfo(float).tiny) if np.isfinite(max_eig) and max_eig > 0 else np.finfo(float).tiny
    eig_pos = np.clip(eig, floor, None)
    cond = float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan
    out = {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0 else np.nan,
        f"{prefix}condition_number": cond,
        f"{prefix}trace": float(np.trace(fim)),
        f"{prefix}trace_inv": float(np.trace(stable_inverse(fim))),
    }
    for name in WEAK_PARAMETERS:
        if name in parameters:
            idx = parameters.index(name)
            out[f"{prefix}diag_{name}"] = float(fim[idx, idx])
    return out


def variance_reduction_metrics(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    combined_cov = stable_inverse(combined)
    out = {}
    ratios = []
    for name in WEAK_PARAMETERS:
        if name not in parameters:
            continue
        idx = parameters.index(name)
        before = float(prior_cov[idx, idx])
        after = float(combined_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        out[f"posterior_var_ratio_{name}"] = ratio
        out[f"posterior_var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        ratios.append(ratio)
    if ratios:
        out["weak_mean_var_reduction"] = float(1.0 - np.nanmean(ratios))
        out["weak_worst_var_reduction"] = float(1.0 - np.nanmax(ratios))
    return out


def evaluate_candidate(
    ns: dict,
    theta: dict[str, float],
    candidate: DesignCandidate,
    batch_id: str,
    parameters: tuple[str, ...],
    prior_fim: np.ndarray,
    step: float,
    solver,
    pulse_slots: int,
) -> tuple[dict, np.ndarray | None]:
    pulses = pulse_schedule_for_batch(ns, batch_id, candidate, pulse_slots=pulse_slots)
    row = {
        "candidate": candidate.name,
        "family": candidate.family,
        "batch_template": str(batch_id),
        "parameters": ", ".join(parameters),
        "temperature_c": ", ".join(f"{v:.3g}" for v in candidate.temperature_c),
        "pulse_times_h": ", ".join(f"{t:.3g}" for t, _a in pulses),
        "pulse_amounts_kg_m3": ", ".join(f"{a:.4g}" for _t, a in pulses),
        "total_pulse_kg_m3": float(sum(a for _t, a in pulses)),
        "status": "failed",
        "error": "",
    }
    try:
        experiment = FixedDesignExperiment(
            ns,
            batch_id=batch_id,
            theta=theta,
            parameters=parameters,
            temperature_c=candidate.temperature_c,
            pulses=pulses,
        )
        fim_new, _jac = fim_for_experiment(ns, experiment, parameters, step, solver)
        combined = prior_fim + fim_new
        row.update(fim_metrics(fim_new, parameters, prefix="new_"))
        row.update(fim_metrics(combined, parameters, prefix="combined_"))
        row.update(variance_reduction_metrics(prior_fim, combined, parameters))
        row["status"] = "ok"
        return row, fim_new
    except Exception as err:
        row["error"] = f"{type(err).__name__}: {err}"
        print(f"[DOE] failed {candidate.name} / {batch_id}: {row['error']}", flush=True)
        return row, None


def simulate_design(ns: dict, theta: dict[str, float], candidate_row: pd.Series, parameters: tuple[str, ...]):
    candidate = DesignCandidate(
        name=str(candidate_row["candidate"]),
        family=str(candidate_row.get("family", "selected")),
        temperature_c=tuple(float(v.strip()) for v in str(candidate_row["temperature_c"]).split(",")),
        pulse_time_fractions=(),
        pulse_amounts=(),
    )
    pulse_times = tuple(float(v.strip()) for v in str(candidate_row["pulse_times_h"]).split(","))
    pulse_amounts = tuple(float(v.strip()) for v in str(candidate_row["pulse_amounts_kg_m3"]).split(","))
    pulses = tuple(zip(pulse_times, pulse_amounts))
    batch_id = str(candidate_row["batch_template"])
    batch = ns["load_batch"](batch_id)
    experiment = ns["FermentationExperiment"](
        batch,
        theta_initial=theta,
        parameters_to_estimate=parameters,
        input_mode="design",
        temperature_segments=len(candidate.temperature_c),
        pulse_max_count=len(pulses),
        fix_design_inputs=False,
    )
    model = experiment.get_labeled_model()
    for segment, value in enumerate(candidate.temperature_c):
        model.T_set[segment].set_value(value)
        model.T_set[segment].fix()
    for slot, (time_h, amount) in enumerate(pulses):
        model.pulse_time[slot].set_value(time_h)
        model.pulse_amount[slot].set_value(amount)
        model.pulse_time[slot].fix()
        model.pulse_amount[slot].fix()
    result = ns["solve_dynamic_model"](model, tee=False)
    sim = ns["extract_simulation_results"](model)
    sim["candidate"] = candidate.name
    sim["batch_template"] = batch_id
    sim["status"] = str(result.solver.status)
    sim["termination"] = str(result.solver.termination_condition)
    return sim


def plot_design_simulations(ns: dict, simulations: pd.DataFrame, selected_rows: pd.DataFrame):
    if simulations.empty:
        return
    for stale in RESULTS_DIR.glob("doe_simulation_*.png"):
        stale.unlink()
    states = list(ns["STATE_LABELS"])
    for candidate_name, group in simulations.groupby("candidate"):
        batch_id = str(group["batch_template"].iloc[0])
        batch = ns["load_batch"](batch_id)
        fig, axes = plt.subplots(len(states), 1, figsize=(9, 11), sharex=True)
        for ax, state in zip(axes, states):
            measured = batch.measurements[state].dropna()
            ax.scatter(measured.index.astype(float), measured.values.astype(float), s=20, color="black", alpha=0.45, label="template data")
            ax.plot(group["t"], group[state], color="tab:green", label="DOE simulation")
            ax.set_ylabel(state)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("t [h]")
        row = selected_rows[selected_rows["candidate"].eq(candidate_name)].iloc[0]
        protocol_line = (
            f"T={row['temperature_c']} | pulses={row['pulse_times_h']} h / "
            f"{row['pulse_amounts_kg_m3']} kg/m3"
        )
        fig.suptitle(
            f"{candidate_name} | batch template {batch_id}\n"
            + "\n".join(textwrap.wrap(protocol_line, width=88)),
            y=0.997,
            fontsize=10,
        )
        axes[0].legend(loc="upper right", fontsize=8, frameon=True)
        fig.tight_layout(rect=(0, 0, 1, 0.91))
        fig.savefig(RESULTS_DIR / f"doe_simulation_{candidate_name}_{batch_id}.png", dpi=180)
        plt.close(fig)


def write_report(primary_summary: pd.DataFrame, final7_summary: pd.DataFrame, extended_summary: pd.DataFrame):
    lines = [
        "# Fermentation DOE design report",
        "",
        "Objective: propose temperature setpoint profiles and YAN pulse schedules that improve identification of the current free parameters and the weak practical-identifiability directions.",
        "",
        "Primary target set: `mu0`, `qN`, `betaG0`, `betaF0`, `qEG`, `qEF`, `iG`, `qXG`, `qXF`.",
        "The accepted calibration model estimates the first seven and keeps `qXG/qXF` fixed; this DOE intentionally targets `qXG/qXF` because they are the main correctable weak directions.",
        "",
        "Ranking criterion: combined current-data FIM plus proposed-experiment FIM, using weighted relative sensitivities from Pyomo DoE. Designs are ranked primarily by mean posterior variance reduction for `qXG`, `qXF`, and `iG`, then by combined FIM condition/minimum eigenvalue/logdet.",
        "",
        "## Recommended Primary Designs",
        "",
        primary_summary.head(8).to_markdown(index=False),
        "",
    ]
    if final7_summary is not None and not final7_summary.empty:
        lines.extend(["## Current Free-Parameter Support", "", final7_summary.head(6).to_markdown(index=False), ""])
    if extended_summary is not None and not extended_summary.empty:
        lines.extend(["## Extended Saturation Screening", "", extended_summary.head(6).to_markdown(index=False), ""])
    lines.extend(
        [
            "## Interpretation",
            "",
            "- High-information designs tend to combine temperature changes with split nutrient additions instead of a single static profile.",
            "- `qXG/qXF` need biomass/sugar-growth separation; therefore designs that alter biomass build-up before sugar depletion are preferred.",
            "- `iG` benefits from high-sugar/high-activity windows early in the run, then temperature changes that separate inhibition from fermentation-rate effects.",
            "- Designs are conditional on the initial conditions of the selected batch template. Re-run this script with a new template if the next fermentation has different starting sugar, ethanol, YAN, or biomass.",
        ]
    )
    (RESULTS_DIR / "doe_design_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Design fermentation experiments with Pyomo DoE FIM screening.")
    parser.add_argument("--templates", nargs="*", default=["25026", "25085", "25170"], help="Batch templates/initial conditions to screen.")
    parser.add_argument("--prior-batches", nargs="*", default=["25026", "25086", "25150", "25170"], help="Existing calibration batches used as prior FIM.")
    parser.add_argument("--random-candidates", type=int, default=12)
    parser.add_argument("--seed", type=int, default=86021)
    parser.add_argument("--step", type=float, default=1e-2)
    parser.add_argument("--pulse-slots", type=int, default=5)
    parser.add_argument("--max-candidates", type=int, default=0, help="Optional cap after manual + random generation; 0 means no cap.")
    parser.add_argument("--skip-final7", action="store_true")
    parser.add_argument("--skip-extended", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_notebook_context()
    theta = series_to_theta(pd.read_csv(FINAL_THETA_PATH, index_col=0).iloc[:, 0])
    solver = make_solver(ns)

    candidates = manual_candidates() + random_candidates(args.random_candidates, args.seed)
    if args.max_candidates and args.max_candidates > 0:
        candidates = candidates[: int(args.max_candidates)]
    candidate_manifest = pd.DataFrame(
        [
            {
                "candidate": c.name,
                "family": c.family,
                "temperature_c": ", ".join(map(str, c.temperature_c)),
                "pulse_time_fractions": ", ".join(map(str, c.pulse_time_fractions)),
                "pulse_amounts": ", ".join(map(str, c.pulse_amounts)),
            }
            for c in candidates
        ]
    )
    candidate_manifest.to_csv(RESULTS_DIR / "candidate_manifest.csv", index=False)

    print("[DOE] Computing/loading prior FIM for primary target set.", flush=True)
    prior_primary = current_prior_fim(ns, theta, PRIMARY9, list(args.prior_batches), args.step, solver)
    prior_primary_metrics = fim_metrics(prior_primary, PRIMARY9, prefix="prior_")
    (pd.DataFrame(prior_primary, index=PRIMARY9, columns=PRIMARY9)).to_csv(RESULTS_DIR / "prior_primary9_fim.csv")
    (pd.Series(prior_primary_metrics, name="prior_primary9_metrics")).to_csv(RESULTS_DIR / "prior_primary9_metrics.csv")

    primary_rows = []
    primary_fims = {}
    total = len(candidates) * len(args.templates)
    counter = 0
    for batch_id in args.templates:
        for candidate in candidates:
            counter += 1
            print(f"[DOE primary] {counter}/{total}: template={batch_id} candidate={candidate.name}", flush=True)
            row, fim = evaluate_candidate(
                ns,
                theta,
                candidate,
                str(batch_id),
                PRIMARY9,
                prior_primary,
                args.step,
                solver,
                args.pulse_slots,
            )
            primary_rows.append(row)
            if fim is not None:
                primary_fims[(row["candidate"], row["batch_template"])] = fim
    primary = pd.DataFrame(primary_rows)
    ok_primary = primary[primary["status"].eq("ok")].copy()
    if not ok_primary.empty:
        ok_primary = ok_primary.sort_values(
            [
                "weak_mean_var_reduction",
                "weak_worst_var_reduction",
                "combined_min_relative_eigenvalue",
                "combined_logdet",
            ],
            ascending=[False, False, False, False],
        )
    primary.to_csv(RESULTS_DIR / "doe_primary9_all_candidates.csv", index=False)
    ok_primary.to_csv(RESULTS_DIR / "doe_primary9_ranked.csv", index=False)

    selected_rows = ok_primary.head(5).copy()
    selected_rows.to_csv(RESULTS_DIR / "doe_recommended_designs.csv", index=False)
    for _, row in selected_rows.iterrows():
        key = (row["candidate"], row["batch_template"])
        if key in primary_fims:
            pd.DataFrame(primary_fims[key], index=PRIMARY9, columns=PRIMARY9).to_csv(
                RESULTS_DIR / f"fim_new_primary9_{row['candidate']}_{row['batch_template']}.csv"
            )

    final7_ranked = pd.DataFrame()
    if not args.skip_final7 and not selected_rows.empty:
        print("[DOE] Evaluating selected designs for final7 support.", flush=True)
        prior_final7 = current_prior_fim(ns, theta, FINAL7, list(args.prior_batches), args.step, solver)
        final7_rows = []
        selected_candidates = {row["candidate"]: row for _, row in selected_rows.iterrows()}
        candidate_by_name = {c.name: c for c in candidates}
        for _idx, selected in selected_rows.iterrows():
            candidate = candidate_by_name[selected["candidate"]]
            row, _fim = evaluate_candidate(
                ns,
                theta,
                candidate,
                str(selected["batch_template"]),
                FINAL7,
                prior_final7,
                args.step,
                solver,
                args.pulse_slots,
            )
            final7_rows.append(row)
        final7_ranked = pd.DataFrame(final7_rows)
        final7_ranked = final7_ranked[final7_ranked["status"].eq("ok")].sort_values(
            ["combined_min_relative_eigenvalue", "combined_logdet"], ascending=[False, False]
        )
        final7_ranked.to_csv(RESULTS_DIR / "doe_final7_selected_ranked.csv", index=False)

    extended_ranked = pd.DataFrame()
    if not args.skip_extended and not selected_rows.empty:
        print("[DOE] Evaluating selected designs for extended saturation screening.", flush=True)
        prior_extended = current_prior_fim(ns, theta, EXTENDED12, list(args.prior_batches), args.step, solver)
        extended_rows = []
        candidate_by_name = {c.name: c for c in candidates}
        for _idx, selected in selected_rows.head(4).iterrows():
            candidate = candidate_by_name[selected["candidate"]]
            row, _fim = evaluate_candidate(
                ns,
                theta,
                candidate,
                str(selected["batch_template"]),
                EXTENDED12,
                prior_extended,
                args.step,
                solver,
                args.pulse_slots,
            )
            extended_rows.append(row)
        extended_ranked = pd.DataFrame(extended_rows)
        extended_ranked = extended_ranked[extended_ranked["status"].eq("ok")].sort_values(
            ["combined_min_relative_eigenvalue", "combined_logdet"], ascending=[False, False]
        )
        extended_ranked.to_csv(RESULTS_DIR / "doe_extended12_selected_ranked.csv", index=False)

    simulations = []
    for _idx, row in selected_rows.head(3).iterrows():
        try:
            sim = simulate_design(ns, theta, row, PRIMARY9)
            simulations.append(sim)
        except Exception as err:
            print(f"[DOE] simulation failed for {row['candidate']}: {type(err).__name__}: {err}", flush=True)
    simulation_frame = pd.concat(simulations, ignore_index=True) if simulations else pd.DataFrame()
    simulation_frame.to_csv(RESULTS_DIR / "doe_selected_simulations.csv", index=False)
    if not selected_rows.empty and not simulation_frame.empty:
        plot_design_simulations(ns, simulation_frame, selected_rows)

    write_report(selected_rows, final7_ranked, extended_ranked)
    metadata = {
        "templates": list(args.templates),
        "prior_batches": list(args.prior_batches),
        "step": float(args.step),
        "pulse_slots": int(args.pulse_slots),
        "n_candidates": len(candidates),
        "target_primary": list(PRIMARY9),
        "target_final7": list(FINAL7),
        "target_extended12": list(EXTENDED12),
    }
    (RESULTS_DIR / "doe_run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("\nRecommended designs:")
    if selected_rows.empty:
        print("No successful designs.")
    else:
        print(
            selected_rows[
                [
                    "candidate",
                    "batch_template",
                    "temperature_c",
                    "pulse_times_h",
                    "pulse_amounts_kg_m3",
                    "weak_mean_var_reduction",
                    "weak_worst_var_reduction",
                    "combined_min_relative_eigenvalue",
                    "combined_condition_number",
                ]
            ].to_string(index=False)
        )
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
