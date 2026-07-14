from __future__ import annotations

import argparse
import math
import os
import sys
import types
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

os.environ.setdefault("MPLBACKEND", "Agg")

import nbformat
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
NOTEBOOK_PATH = SCRIPT_DIR / "fermentation_model_calibration_3_effective_reformulation.ipynb"
RESULTS_DIR = SCRIPT_DIR / "results" / "fit_strategy_analysis"


def _quiet_display(*_args, **_kwargs):
    return None


def load_notebook_context() -> dict:
    """Execute only the reusable setup cells from the calibration notebook."""
    os.chdir(REPO_ROOT)
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    module_name = "__fit_strategy_notebook_context__"
    module = types.ModuleType(module_name)
    module.__file__ = str(NOTEBOOK_PATH)
    sys.modules[module_name] = module
    ns: dict = module.__dict__

    def run_cell(index: int, source_transform: Callable[[str], str] | None = None) -> None:
        source = nb.cells[index].source
        if source_transform is not None:
            source = source_transform(source)
        exec(compile(source, f"{NOTEBOOK_PATH.name}:cell_{index}", "exec"), ns)

    run_cell(1)
    ns["display"] = _quiet_display
    for index in [3, 5, 8, 10, 12]:
        run_cell(index)

    def only_simulation_functions(source: str) -> str:
        marker = "\n\nmodel, solve_results, sim = simulate_batch(batch)"
        if marker in source:
            return source.split(marker, 1)[0]
        return source

    run_cell(14, only_simulation_functions)
    for index in [19, 21, 25, 27]:
        run_cell(index)
    return ns


@dataclass(frozen=True)
class FitStrategy:
    name: str
    estimated_parameters: tuple[str, ...]
    initial_source: str = "nominal"
    l2_lambda: float = 0.0
    l2_parameters: tuple[str, ...] = ()
    pso_seed: int | None = None

    @property
    def uses_l2(self) -> bool:
        return self.l2_lambda > 0.0 and bool(self.l2_parameters)

    @property
    def uses_pso(self) -> bool:
        return self.pso_seed is not None


def series_to_theta(theta_like) -> dict[str, float]:
    if isinstance(theta_like, pd.Series):
        data = theta_like.to_dict()
    elif isinstance(theta_like, pd.DataFrame):
        if theta_like.shape[1] == 1:
            data = theta_like.iloc[:, 0].to_dict()
        else:
            raise ValueError("Cannot convert multi-column DataFrame to theta.")
    else:
        data = dict(theta_like)
    return {str(k): float(v) for k, v in data.items()}


def reference_theta(ns: dict, source: str, batch_id: str = "25026") -> dict[str, float]:
    load_batch = ns["load_batch"]
    theta_for_batch = ns["_theta_for_batch"]
    clip_theta_to_bounds = ns["clip_theta_to_bounds"]
    default_theta = ns["DEFAULT_THETA"]

    theta = theta_for_batch(load_batch(batch_id), default_theta)
    if source == "reduced6":
        reduced_path = SCRIPT_DIR / "results" / "identifiability_reduction" / "theta_reduced6_nominal_fixed.csv"
        if reduced_path.exists():
            theta.update(series_to_theta(pd.read_csv(reduced_path, index_col=0).iloc[:, 0]))
    elif source != "nominal":
        raise ValueError(f"Unknown initial source: {source}")
    clipped, _ = clip_theta_to_bounds(theta)
    return {name: float(clipped[name]) for name in default_theta}


def complete_theta(ns: dict, theta_partial, initial_theta: dict[str, float]) -> dict[str, float]:
    theta = dict(initial_theta)
    theta.update(series_to_theta(theta_partial))
    clipped, _ = ns["clip_theta_to_bounds"](theta)
    return {name: float(clipped[name]) for name in ns["DEFAULT_THETA"]}


def make_l2_log_objective(ns: dict, reference: dict[str, float], l2_parameters: tuple[str, ...], l2_lambda: float):
    pyo = ns["pyo"]
    parmest = ns["parmest"]
    bounds = ns["PARAMETER_BOUNDS"]
    penalized = tuple(l2_parameters)

    def l2_log_sse_weighted(model):
        expr = parmest.SSE_weighted(model)
        penalty = 0.0
        for param in model.unknown_parameters:
            name = param.name
            if name not in penalized:
                continue
            ref = float(reference[name])
            lb, _ub = bounds[name]
            if ref <= 0.0 or lb <= 0.0:
                scale = max(abs(ref), 1e-8)
                penalty += ((param - ref) / scale) ** 2
            else:
                penalty += pyo.log(param / ref) ** 2
        return expr + 0.5 * float(l2_lambda) * penalty

    l2_log_sse_weighted.__name__ = f"SSE_weighted_logL2_lambda_{l2_lambda:g}"
    return l2_log_sse_weighted


def log_l2_penalty_total(
    theta: dict[str, float],
    reference: dict[str, float],
    l2_parameters: tuple[str, ...],
    l2_lambda: float,
    n_batches: int,
) -> float:
    if l2_lambda <= 0.0 or not l2_parameters:
        return 0.0
    total = 0.0
    for name in l2_parameters:
        value = float(theta[name])
        ref = float(reference[name])
        if value <= 0.0 or ref <= 0.0:
            scale = max(abs(ref), 1e-8)
            total += ((value - ref) / scale) ** 2
        else:
            total += math.log(value / ref) ** 2
    return 0.5 * float(l2_lambda) * total * int(n_batches)


def run_parmest_strategy(
    ns: dict,
    strategy: FitStrategy,
    batches: list[str],
    initial_theta: dict[str, float],
    tee: bool = False,
) -> tuple[float, dict[str, float], object, float]:
    if strategy.uses_l2:
        obj_function = make_l2_log_objective(
            ns,
            reference=initial_theta,
            l2_parameters=strategy.l2_parameters,
            l2_lambda=strategy.l2_lambda,
        )
    else:
        obj_function = "SSE_weighted"

    obj_value, theta_estimated, estimator = ns["run_parmest_estimation"](
        obj_function=obj_function,
        batch_ids=batches,
        tee=tee,
        theta_initial=initial_theta,
        parameters_to_estimate=list(strategy.estimated_parameters),
    )
    theta = complete_theta(ns, theta_estimated, initial_theta)
    penalty = log_l2_penalty_total(
        theta=theta,
        reference=initial_theta,
        l2_parameters=strategy.l2_parameters,
        l2_lambda=strategy.l2_lambda,
        n_batches=len(batches),
    )
    return float(obj_value), theta, estimator, penalty


def run_pso_then_polish(
    ns: dict,
    strategy: FitStrategy,
    batches: list[str],
    initial_theta: dict[str, float],
    pso_epoch: int,
    pso_pop_size: int,
) -> tuple[float, dict[str, float], object, float, dict]:
    pso_dir = RESULTS_DIR / "pso_work" / strategy.name
    pso_dir.mkdir(parents=True, exist_ok=True)
    ns["CUSTOM_PSO_BATCHES"] = list(batches)
    ns["CUSTOM_PSO_PARAMETERS"] = list(strategy.estimated_parameters)
    ns["CUSTOM_PSO_FIXED_PARAMETERS"] = set(ns["DEFAULT_THETA"]) - set(strategy.estimated_parameters)
    ns["CUSTOM_PSO_LOG_PARAMETERS"] = [
        name for name in strategy.estimated_parameters if ns["PARAMETER_BOUNDS"][name][0] > 0.0
    ]
    ns["CUSTOM_PSO_REFERENCE_THETA"] = dict(initial_theta)
    ns["CUSTOM_PSO_INCLUDE_WSSE_SEED"] = True
    ns["theta_WSSE"] = dict(initial_theta)
    ns["CUSTOM_PSO_CONFIG"] = {
        "epoch": int(pso_epoch),
        "pop_size": int(pso_pop_size),
        "w": 0.5,
        "c1": 1.5,
        "c2": 1.5,
        "seed": int(strategy.pso_seed),
        "verbose": False,
        "save_history": False,
        "relative_gap_threshold": 5e-4,
    }
    ns["CUSTOM_PSO_CHECK_INTERVAL"] = max(3, min(6, int(pso_epoch)))
    ns["CUSTOM_PSO_PATIENCE_CHECKS"] = 2
    ns["CUSTOM_PSO_PRINT_PROGRESS"] = True
    ns["CUSTOM_PSO_PROGRESS_EVERY_EVALUATIONS"] = max(1, int(pso_pop_size))
    ns["CUSTOM_PSO_RESULTS_DIR"] = pso_dir
    ns["CUSTOM_PSO_RESULT_PREFIX"] = strategy.name

    pso_result = ns["run_custom_pso"]()
    pso_theta = complete_theta(ns, pso_result["best_theta"], initial_theta)
    pd.Series(pso_theta, name=strategy.name).to_csv(RESULTS_DIR / f"{strategy.name}_pso_best_theta.csv")
    pso_result["progress"].to_csv(RESULTS_DIR / f"{strategy.name}_pso_progress.csv", index=False)

    polish_strategy = FitStrategy(
        name=f"{strategy.name}_polish",
        estimated_parameters=strategy.estimated_parameters,
        initial_source=strategy.initial_source,
        l2_lambda=strategy.l2_lambda,
        l2_parameters=strategy.l2_parameters,
        pso_seed=None,
    )
    obj_value, theta, estimator, penalty = run_parmest_strategy(
        ns,
        polish_strategy,
        batches=batches,
        initial_theta=pso_theta,
        tee=False,
    )
    return obj_value, theta, estimator, penalty, pso_result


def active_bound_summary(ns: dict, theta: dict[str, float], estimated_parameters: tuple[str, ...]) -> tuple[int, str]:
    active = []
    for name in estimated_parameters:
        value = float(theta[name])
        lb, ub = ns["PARAMETER_BOUNDS"][name]
        width = float(ub) - float(lb)
        if width <= 0.0:
            continue
        near_lower = (value - float(lb)) / width <= 1e-3
        near_upper = (float(ub) - value) / width <= 1e-3
        if near_lower:
            active.append(f"{name}:lower")
        elif near_upper:
            active.append(f"{name}:upper")
    return len(active), "; ".join(active)


def evaluate_strategy(
    ns: dict,
    strategy: FitStrategy,
    theta: dict[str, float],
    parmest_objective_total: float,
    l2_penalty_total: float,
    batches: list[str],
    runtime_s: float,
    extra: dict | None = None,
) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    objective_summary, residuals, solve_summary = ns["evaluate_fit_objectives"](theta, batch_ids=batches)
    active_count, active_names = active_bound_summary(ns, theta, strategy.estimated_parameters)
    state_summary = pd.DataFrame()
    if not residuals.empty:
        state_summary = (
            residuals.groupby("state")[["squared_error", "weighted_squared_error"]]
            .sum()
            .rename(columns={"squared_error": "SSE_raw_sum", "weighted_squared_error": "WSSE_raw_sum"})
        )
        state_summary["strategy"] = strategy.name
        state_summary = state_summary.reset_index().set_index(["strategy", "state"])
    row = {
        "strategy": strategy.name,
        "estimated_parameters": ", ".join(strategy.estimated_parameters),
        "n_estimated": len(strategy.estimated_parameters),
        "l2_lambda": float(strategy.l2_lambda),
        "l2_parameters": ", ".join(strategy.l2_parameters),
        "uses_pso": bool(strategy.uses_pso),
        "pso_seed": strategy.pso_seed if strategy.pso_seed is not None else np.nan,
        "ParmEst_objective_total": float(parmest_objective_total),
        "direct_WSSE_total": float(objective_summary["WSSE_raw_sum"]),
        "direct_SSE_total": float(objective_summary["SSE_raw_sum"]),
        "l2_penalty_total": float(l2_penalty_total),
        "direct_plus_l2_total": float(objective_summary["WSSE_raw_sum"]) + float(l2_penalty_total),
        "n_observations": int(objective_summary["n_observations"]),
        "n_active_bounds": int(active_count),
        "active_bounds": active_names,
        "runtime_s": float(runtime_s),
        "solve_ok": bool(
            not solve_summary.empty
            and solve_summary["termination"].astype(str).str.contains("optimal", case=False, na=False).all()
        ),
    }
    if extra:
        row.update(extra)
    return row, residuals, solve_summary, state_summary


def save_predictions(ns: dict, theta_by_strategy: dict[str, dict[str, float]], batches: list[str]) -> pd.DataFrame:
    records = []
    for strategy_name, theta in theta_by_strategy.items():
        for batch_id in batches:
            try:
                _model, solve_result, sim = ns["simulate_batch"](ns["load_batch"](batch_id), theta_initial=theta)
                sim = sim.copy()
                sim["strategy"] = strategy_name
                sim["batch"] = str(batch_id)
                sim["status"] = str(solve_result.solver.status)
                sim["termination"] = str(solve_result.solver.termination_condition)
                records.append(sim)
            except Exception as err:
                records.append(
                    pd.DataFrame(
                        [
                            {
                                "strategy": strategy_name,
                                "batch": str(batch_id),
                                "status": "failed",
                                "termination": f"{type(err).__name__}: {err}",
                            }
                        ]
                    )
                )
    predictions = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    predictions.to_csv(RESULTS_DIR / "strategy_predictions.csv", index=False)
    return predictions


def plot_strategy_curves(ns: dict, predictions: pd.DataFrame, selected_strategies: list[str], batches: list[str]) -> None:
    if predictions.empty:
        return
    import matplotlib.pyplot as plt

    state_labels = list(ns["STATE_LABELS"])
    palette = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
    colors = {name: palette[i % len(palette)] for i, name in enumerate(selected_strategies)}
    for batch_id in batches:
        batch = ns["load_batch"](batch_id)
        fig, axes = plt.subplots(len(state_labels), 1, figsize=(9, 11), sharex=True)
        for ax, state in zip(axes, state_labels):
            measured = batch.measurements[state].dropna()
            ax.scatter(measured.index.astype(float), measured.values.astype(float), s=22, color="black", label="data")
            for strategy_name in selected_strategies:
                subset = predictions[
                    predictions["strategy"].eq(strategy_name)
                    & predictions["batch"].astype(str).eq(str(batch_id))
                    & predictions[state].notna()
                ].sort_values("t")
                if subset.empty:
                    continue
                ax.plot(
                    subset["t"].astype(float),
                    subset[state].astype(float),
                    color=colors.get(strategy_name),
                    label=strategy_name,
                    alpha=0.9,
                )
            ax.set_ylabel(state)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("t [h]")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=min(3, len(labels)))
        fig.suptitle(f"Fit strategy comparison - batch {batch_id}", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(RESULTS_DIR / f"strategy_curves_batch_{batch_id}.png", dpi=180)
        plt.close(fig)


def default_strategies(include_pso: bool) -> list[FitStrategy]:
    reduced6 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF")
    plus_qx = reduced6 + ("qXG", "qXF")
    plus_ig = reduced6 + ("iG",)
    plus_qx_ig = reduced6 + ("qXG", "qXF", "iG")
    strategies = [
        FitStrategy("reduced6_local", reduced6, initial_source="reduced6"),
        FitStrategy("plus_qx_unregularized", plus_qx, initial_source="reduced6"),
        FitStrategy("plus_iG_unregularized", plus_ig, initial_source="reduced6"),
        FitStrategy("plus_qx_iG_unregularized", plus_qx_ig, initial_source="reduced6"),
        FitStrategy("plus_qx_iG_logL2_10", plus_qx_ig, initial_source="reduced6", l2_lambda=10.0, l2_parameters=("qXG", "qXF", "iG")),
        FitStrategy("plus_qx_iG_logL2_50", plus_qx_ig, initial_source="reduced6", l2_lambda=50.0, l2_parameters=("qXG", "qXF", "iG")),
        FitStrategy("plus_qx_iG_logL2_100", plus_qx_ig, initial_source="reduced6", l2_lambda=100.0, l2_parameters=("qXG", "qXF", "iG")),
        FitStrategy("plus_qx_logL2_50", plus_qx, initial_source="reduced6", l2_lambda=50.0, l2_parameters=("qXG", "qXF")),
        FitStrategy("plus_iG_logL2_50", plus_ig, initial_source="reduced6", l2_lambda=50.0, l2_parameters=("iG",)),
    ]
    if include_pso:
        strategies.append(
            FitStrategy(
                "plus_qx_iG_pso_pilot_logL2_50",
                plus_qx_ig,
                initial_source="reduced6",
                l2_lambda=50.0,
                l2_parameters=("qXG", "qXF", "iG"),
                pso_seed=321,
            )
        )
    return strategies


def choose_recommended_strategy(summary: pd.DataFrame) -> dict[str, str | float]:
    reduced_name = "reduced6_local"
    reduced_wsse = float(summary.loc[reduced_name, "direct_WSSE_total"]) if reduced_name in summary.index else np.nan
    valid = summary[summary["solve_ok"].fillna(False)].copy()
    no_active = valid[valid["n_active_bounds"].fillna(999).astype(float).eq(0.0)].copy()
    regularized = no_active[no_active["l2_lambda"].fillna(0.0).astype(float).gt(0.0)].copy()

    best_free = valid.sort_values(["direct_WSSE_total", "direct_plus_l2_total"], na_position="last").head(1)
    best_stable = no_active.sort_values(["direct_WSSE_total", "direct_plus_l2_total"], na_position="last").head(1)
    best_regularized = regularized.sort_values(["direct_plus_l2_total", "direct_WSSE_total"], na_position="last").head(1)

    recommended = reduced_name
    reason = "Only the reduced model satisfied the stability filters."
    if not best_regularized.empty:
        candidate_name = str(best_regularized.index[0])
        candidate = best_regularized.iloc[0]
        if np.isfinite(reduced_wsse) and float(candidate["direct_plus_l2_total"]) < reduced_wsse:
            recommended = candidate_name
            reason = (
                "Best log-L2 candidate with no active estimated-parameter bounds; "
                "its penalized objective is still below the reduced-model WSSE."
            )
        else:
            reason = (
                "Regularized candidates were stable, but their penalized objective did not improve on the reduced model."
            )

    return {
        "recommended": recommended,
        "reason": reason,
        "best_free": str(best_free.index[0]) if not best_free.empty else "",
        "best_stable": str(best_stable.index[0]) if not best_stable.empty else "",
        "best_regularized": str(best_regularized.index[0]) if not best_regularized.empty else "",
        "reduced_wsse": reduced_wsse,
    }


def write_decision_report(summary: pd.DataFrame, state_summary: pd.DataFrame, selected: list[str]) -> None:
    report_path = RESULTS_DIR / "fit_strategy_decision_report.md"
    decision = choose_recommended_strategy(summary)
    recommended = str(decision["recommended"])
    recommended_row = summary.loc[recommended]
    lines = [
        "# Fit strategy analysis",
        "",
        f"Recommended regularized curve-fit candidate: `{recommended}`.",
        "",
        str(decision["reason"]),
        "",
        f"- Best free fit: `{decision['best_free']}`.",
        f"- Best stable direct-WSSE fit: `{decision['best_stable']}`.",
        f"- Best stable log-L2 fit: `{decision['best_regularized']}`.",
        "",
        "Selection metric for the regularized candidate: direct weighted SSE on the data plus the explicit log-L2 penalty.",
        "Fits with active estimated-parameter bounds are reported as sensitivity fits, not as recommended estimable models.",
        "",
        "## Recommended follow-up profile set",
        "",
        f"Run profile-likelihood diagnostics for `{recommended}` and compare it against `reduced6_local`.",
        "If the released parameters do not cross the 95% profile threshold without relying on the prior penalty, keep them fixed for DOE and use this candidate only as a Bayesian regularized curve-fit scenario.",
        "",
        "Recommended candidate metrics:",
        "",
        recommended_row.to_frame("value").to_markdown(),
        "",
        "## Top strategies",
        "",
        summary.head(8).to_markdown(),
        "",
    ]
    if not state_summary.empty:
        lines.extend(["## State residuals for selected strategies", "", state_summary.loc[selected].to_markdown(), ""])
    report_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pso", action="store_true", help="Skip the PSO pilot candidate.")
    parser.add_argument("--pso-epoch", type=int, default=12, help="Epochs for the PSO pilot.")
    parser.add_argument("--pso-pop-size", type=int, default=6, help="Population size for the PSO pilot.")
    parser.add_argument(
        "--batches",
        nargs="*",
        default=["25026", "25086", "25150", "25170"],
        help="Calibration batches to use.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_notebook_context()
    batches = [str(batch) for batch in args.batches]
    strategies = default_strategies(include_pso=not args.skip_pso)

    summary_rows = []
    residual_tables = []
    solve_tables = []
    state_tables = []
    theta_by_strategy: dict[str, dict[str, float]] = {}

    import time

    for strategy in strategies:
        print(f"\n=== Strategy: {strategy.name} ===", flush=True)
        initial_theta = reference_theta(ns, strategy.initial_source, batch_id=batches[0])
        t0 = time.perf_counter()
        try:
            if strategy.uses_pso:
                obj_value, theta, _estimator, penalty, pso_result = run_pso_then_polish(
                    ns,
                    strategy=strategy,
                    batches=batches,
                    initial_theta=initial_theta,
                    pso_epoch=args.pso_epoch,
                    pso_pop_size=args.pso_pop_size,
                )
                extra = {
                    "pso_best_direct_WSSE": float(pso_result["best_objective"]),
                    "pso_evaluations": int(pso_result["n_evaluations"]),
                    "pso_successful_evaluations": int(pso_result["n_successful_evaluations"]),
                    "pso_failed_evaluations": int(pso_result["n_failed_evaluations"]),
                }
            else:
                obj_value, theta, _estimator, penalty = run_parmest_strategy(
                    ns,
                    strategy=strategy,
                    batches=batches,
                    initial_theta=initial_theta,
                    tee=False,
                )
                extra = {}
            runtime_s = time.perf_counter() - t0
            row, residuals, solve_summary, state_summary = evaluate_strategy(
                ns,
                strategy=strategy,
                theta=theta,
                parmest_objective_total=obj_value,
                l2_penalty_total=penalty,
                batches=batches,
                runtime_s=runtime_s,
                extra=extra,
            )
            summary_rows.append(row)
            theta_by_strategy[strategy.name] = theta
            residuals["strategy"] = strategy.name
            residual_tables.append(residuals)
            solve_summary = solve_summary.copy()
            solve_summary["strategy"] = strategy.name
            solve_tables.append(solve_summary.reset_index())
            if not state_summary.empty:
                state_tables.append(state_summary)
            pd.Series(theta, name=strategy.name).to_csv(RESULTS_DIR / f"{strategy.name}_theta.csv")
            ns["physical_parameter_table"](ns["load_batch"](batches[0]), theta).to_csv(
                RESULTS_DIR / f"{strategy.name}_physical.csv"
            )
            print(
                f"{strategy.name}: direct_WSSE={row['direct_WSSE_total']:.6g}, "
                f"direct+L2={row['direct_plus_l2_total']:.6g}, "
                f"active_bounds={row['n_active_bounds']}, runtime={runtime_s:.1f}s",
                flush=True,
            )
        except Exception as err:
            runtime_s = time.perf_counter() - t0
            summary_rows.append(
                {
                    "strategy": strategy.name,
                    "estimated_parameters": ", ".join(strategy.estimated_parameters),
                    "n_estimated": len(strategy.estimated_parameters),
                    "l2_lambda": float(strategy.l2_lambda),
                    "l2_parameters": ", ".join(strategy.l2_parameters),
                    "uses_pso": bool(strategy.uses_pso),
                    "pso_seed": strategy.pso_seed if strategy.pso_seed is not None else np.nan,
                    "ParmEst_objective_total": np.nan,
                    "direct_WSSE_total": np.nan,
                    "direct_SSE_total": np.nan,
                    "l2_penalty_total": np.nan,
                    "direct_plus_l2_total": np.nan,
                    "n_observations": np.nan,
                    "n_active_bounds": np.nan,
                    "active_bounds": "",
                    "runtime_s": runtime_s,
                    "solve_ok": False,
                    "error": f"{type(err).__name__}: {err}",
                }
            )
            print(f"{strategy.name}: failed after {runtime_s:.1f}s: {type(err).__name__}: {err}", flush=True)

    summary = pd.DataFrame(summary_rows).set_index("strategy")
    summary = summary.sort_values(["direct_plus_l2_total", "direct_WSSE_total"], na_position="last")
    summary.to_csv(RESULTS_DIR / "fit_strategy_summary.csv")

    residual_summary = pd.concat(residual_tables, ignore_index=True) if residual_tables else pd.DataFrame()
    residual_summary.to_csv(RESULTS_DIR / "fit_strategy_residuals.csv", index=False)
    solve_summary = pd.concat(solve_tables, ignore_index=True) if solve_tables else pd.DataFrame()
    solve_summary.to_csv(RESULTS_DIR / "fit_strategy_solve_summary.csv", index=False)
    state_summary = pd.concat(state_tables).sort_index() if state_tables else pd.DataFrame()
    state_summary.to_csv(RESULTS_DIR / "fit_strategy_state_summary.csv")

    selected = ["reduced6_local"]
    if not summary.empty:
        decision = choose_recommended_strategy(summary)
        for candidate in [str(decision["best_free"]), str(decision["recommended"])]:
            if candidate and candidate in theta_by_strategy and candidate not in selected:
                selected.append(candidate)
    predictions = save_predictions(ns, {name: theta_by_strategy[name] for name in selected if name in theta_by_strategy}, batches)
    plot_strategy_curves(ns, predictions, [name for name in selected if name in theta_by_strategy], batches)
    write_decision_report(summary, state_summary, [name for name in selected if name in theta_by_strategy])

    print("\nTop strategies:")
    print(summary.head(8)[["n_estimated", "direct_WSSE_total", "l2_penalty_total", "direct_plus_l2_total", "n_active_bounds", "solve_ok"]])
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
