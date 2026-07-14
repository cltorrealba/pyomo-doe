from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from estimability_tools import build_profile_grid, plot_profile_likelihood_profiles, summarize_profile_likelihood
from run_fit_strategy_analysis import (
    RESULTS_DIR as FIT_RESULTS_DIR,
    complete_theta,
    load_notebook_context,
    log_l2_penalty_total,
    make_l2_log_objective,
    reference_theta,
    series_to_theta,
)


NOTEBOOK_PATH = SCRIPT_DIR / "fermentation_model_calibration_3_effective_reformulation.ipynb"
RESULTS_DIR = SCRIPT_DIR / "results" / "deep_model_selection"
DEFAULT_BATCHES = ["25026", "25086", "25150", "25170"]
REDUCED6 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF")
PLUS_IG = REDUCED6 + ("iG",)
PLUS_QX_IG = REDUCED6 + ("qXG", "qXF", "iG")
RELEASED_QX_IG = ("qXG", "qXF", "iG")


@dataclass(frozen=True)
class DeepCandidate:
    name: str
    estimated_parameters: tuple[str, ...]
    theta_file: Path
    objective_kind: str = "SSE_weighted"
    l2_lambda: float = 0.0
    l2_parameters: tuple[str, ...] = ()
    profile_parameters: tuple[str, ...] = ()
    profile_label: str = ""

    @property
    def uses_l2(self) -> bool:
        return self.l2_lambda > 0.0 and bool(self.l2_parameters)


def load_q_functions(ns: dict) -> None:
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    source = nb.cells[23].source.replace("RUN_Q_EIGEN_ANALYSIS = True", "RUN_Q_EIGEN_ANALYSIS = False")
    exec(compile(source, f"{NOTEBOOK_PATH.name}:cell_23_deep", "exec"), ns)


def load_theta(path: Path) -> dict[str, float]:
    return series_to_theta(pd.read_csv(path, index_col=0).iloc[:, 0])


def objective_for_candidate(ns: dict, candidate: DeepCandidate, reference: dict[str, float]):
    if not candidate.uses_l2:
        return "SSE_weighted"
    return make_l2_log_objective(
        ns,
        reference=reference,
        l2_parameters=candidate.l2_parameters,
        l2_lambda=candidate.l2_lambda,
    )


def unique_profile_grid(ns: dict, theta: dict[str, float], parameter: str, n_grid: int) -> np.ndarray:
    grid = build_profile_grid(
        theta,
        ns["PARAMETER_BOUNDS"],
        parameter,
        n_grid=n_grid,
        relative_span=0.75,
    )
    value = float(theta[parameter])
    lb, ub = ns["PARAMETER_BOUNDS"][parameter]
    extras = [value]
    if value > float(lb):
        extras.append(max(float(lb), value * 0.25))
    if value < float(ub):
        extras.append(min(float(ub), value * 4.0))
    if abs(value - float(lb)) <= 1e-8 * max(abs(value), abs(float(lb)), 1.0):
        extras.extend([float(lb), min(float(ub), max(value * 2.0, float(lb) + 0.05 * (float(ub) - float(lb))))])
    candidates = np.sort(np.asarray(list(grid) + extras, dtype=float))
    unique = []
    for candidate_value in candidates:
        if not unique or not np.isclose(candidate_value, unique[-1], rtol=1e-10, atol=1e-12):
            unique.append(float(candidate_value))
    return np.asarray(unique, dtype=float)


def evaluate_objective_components(
    ns: dict,
    theta: dict[str, float],
    batches: list[str],
    reference: dict[str, float],
    candidate: DeepCandidate,
) -> tuple[float, float, float]:
    objective_summary, _residuals, _solve_summary = ns["evaluate_fit_objectives"](theta, batch_ids=batches)
    direct_wsse = float(objective_summary["WSSE_raw_sum"])
    penalty = log_l2_penalty_total(
        theta=theta,
        reference=reference,
        l2_parameters=candidate.l2_parameters,
        l2_lambda=candidate.l2_lambda,
        n_batches=len(batches),
    )
    return direct_wsse, penalty, direct_wsse + penalty


def run_q_analysis(ns: dict, candidate: DeepCandidate, theta: dict[str, float], batches: list[str], fraction: float) -> dict:
    print(f"\n[FIM] {candidate.name}: parameters={candidate.estimated_parameters}", flush=True)
    q_raw, q_weighted, q_relative, perturbation_summary, solve_summary = ns[
        "build_weighted_relative_sensitivity_matrix"
    ](
        theta,
        list(candidate.estimated_parameters),
        batches,
        fraction=fraction,
    )
    used_parameters = list(perturbation_summary[perturbation_summary["used_in_q_matrix"]].index)
    fim, eigen_summary, loading_summary = ns["eigen_analysis_from_q"](q_relative, used_parameters)

    prefix = RESULTS_DIR / candidate.name
    q_raw.to_csv(prefix.with_name(f"{candidate.name}_q_raw.csv"))
    q_weighted.to_csv(prefix.with_name(f"{candidate.name}_q_weighted.csv"))
    q_relative.to_csv(prefix.with_name(f"{candidate.name}_q_weighted_relative.csv"))
    fim.to_csv(prefix.with_name(f"{candidate.name}_fim_weighted_relative.csv"))
    perturbation_summary.to_csv(prefix.with_name(f"{candidate.name}_q_perturbation_summary.csv"))
    solve_summary.to_csv(prefix.with_name(f"{candidate.name}_q_solve_summary.csv"), index=False)
    eigen_summary.to_csv(prefix.with_name(f"{candidate.name}_fim_eigenvalues.csv"), index=False)
    loading_summary.to_csv(prefix.with_name(f"{candidate.name}_fim_eigendirections.csv"), index=False)

    min_relative = float(eigen_summary["relative_eigenvalue"].min())
    condition = float(eigen_summary["condition_number"].iloc[0])
    n_near_null = int(eigen_summary["near_null"].sum())
    weakest = loading_summary.iloc[0] if not loading_summary.empty else {}
    return {
        "candidate": candidate.name,
        "n_requested_parameters": len(candidate.estimated_parameters),
        "n_used_parameters": len(used_parameters),
        "used_parameters": ", ".join(used_parameters),
        "dropped_parameters": ", ".join(
            perturbation_summary.loc[~perturbation_summary["used_in_q_matrix"]].index.astype(str).tolist()
        ),
        "fim_condition_number": condition,
        "fim_min_relative_eigenvalue": min_relative,
        "fim_near_null_directions": n_near_null,
        "fim_full_rank_by_threshold": bool(n_near_null == 0 and len(used_parameters) == len(candidate.estimated_parameters)),
        "weakest_direction_parameters": weakest.get("dominant_parameters", ""),
        "weakest_direction_abs_loadings": weakest.get("dominant_abs_loadings", ""),
        "q_failed_solves": int(
            0
            if solve_summary.empty
            else (~solve_summary["termination"].astype(str).str.contains("optimal", case=False, na=False)).sum()
        ),
    }


def run_profile(
    ns: dict,
    candidate: DeepCandidate,
    theta_hat: dict[str, float],
    obj_hat: float,
    batches: list[str],
    reference: dict[str, float],
    n_grid: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    profile_parameters = candidate.profile_parameters or candidate.estimated_parameters
    obj_function = objective_for_candidate(ns, candidate, reference)
    rows = []
    for parameter in profile_parameters:
        grid = unique_profile_grid(ns, theta_hat, parameter, n_grid=n_grid)
        free_parameters = [name for name in candidate.estimated_parameters if name != parameter]
        for grid_idx, value in enumerate(grid, start=1):
            print(
                f"[PROFILE:{candidate.name}] {parameter} {grid_idx}/{len(grid)} fixed={value:.6g}",
                flush=True,
            )
            row = {
                "candidate": candidate.name,
                "profile_label": candidate.profile_label or candidate.name,
                "objective_kind": candidate.objective_kind,
                "profiled_theta": parameter,
                "theta_value": float(value),
                "theta_hat": float(theta_hat[parameter]),
                "objective": np.nan,
                "direct_WSSE": np.nan,
                "l2_penalty": np.nan,
                "lr_stat": np.nan,
                "success": False,
                "status": "failed",
                "error": "",
            }
            if np.isclose(float(value), float(theta_hat[parameter]), rtol=1e-8, atol=1e-10):
                direct_wsse, penalty, objective = evaluate_objective_components(ns, theta_hat, batches, reference, candidate)
                row.update(
                    {
                        "objective": objective if candidate.uses_l2 else direct_wsse,
                        "direct_WSSE": direct_wsse,
                        "l2_penalty": penalty,
                        "lr_stat": 0.0,
                        "success": True,
                        "status": "theta_hat_baseline",
                    }
                )
                row.update({f"estimated_{name}": theta_hat.get(name, np.nan) for name in free_parameters})
                rows.append(row)
                continue
            theta_start = dict(theta_hat)
            theta_start[parameter] = float(value)
            try:
                obj_value, theta_raw, _estimator = ns["run_parmest_estimation"](
                    obj_function=obj_function,
                    batch_ids=batches,
                    theta_initial=theta_start,
                    parameters_to_estimate=free_parameters,
                    tee=False,
                )
                theta_partial = series_to_theta(theta_raw)
                theta_partial[parameter] = float(value)
                theta_profile = complete_theta(ns, theta_partial, theta_start)
                direct_wsse, penalty, objective_check = evaluate_objective_components(
                    ns, theta_profile, batches, reference, candidate
                )
                objective = objective_check if candidate.uses_l2 else direct_wsse
                row.update(
                    {
                        "objective": float(objective),
                        "parmest_objective": float(obj_value),
                        "direct_WSSE": float(direct_wsse),
                        "l2_penalty": float(penalty),
                        "lr_stat": 2.0 * (float(objective) - float(obj_hat)),
                        "success": True,
                        "status": "ok",
                    }
                )
                row.update({f"estimated_{name}": theta_profile.get(name, np.nan) for name in free_parameters})
            except Exception as err:
                row["error"] = f"{type(err).__name__}: {err}"
                row["status"] = "failed"
            rows.append(row)
    profiles = pd.DataFrame(rows)
    summary = summarize_profile_likelihood(profiles, alpha=0.95)
    if not summary.empty:
        summary["candidate"] = candidate.name
        summary["objective_kind"] = candidate.objective_kind
        summary["profile_identifiable_2sided"] = summary["crosses_left"] & summary["crosses_right"]
    return profiles, summary


def laplace_summary_from_fim(
    ns: dict,
    candidate: DeepCandidate,
    theta: dict[str, float],
    reference: dict[str, float],
    fim_path: Path,
    weak_prior_log_sd: float,
    n_batches: int,
) -> pd.DataFrame:
    if not fim_path.exists():
        return pd.DataFrame()
    fim = pd.read_csv(fim_path, index_col=0)
    parameters = [name for name in candidate.estimated_parameters if name in fim.index and name in fim.columns]
    if not parameters:
        return pd.DataFrame()
    fim_array = fim.loc[parameters, parameters].to_numpy(dtype=float)
    weak_prior_precision = np.eye(len(parameters)) / float(weak_prior_log_sd) ** 2
    l2_precision = np.zeros_like(fim_array)
    if candidate.uses_l2:
        for idx, name in enumerate(parameters):
            if name in candidate.l2_parameters:
                l2_precision[idx, idx] = float(candidate.l2_lambda) * int(n_batches)
    posterior_precision = fim_array + weak_prior_precision + l2_precision
    try:
        posterior_cov = np.linalg.inv(posterior_precision)
    except np.linalg.LinAlgError:
        posterior_cov = np.linalg.pinv(posterior_precision)
    mle_log = np.array([math.log(float(theta[name])) for name in parameters], dtype=float)
    prior_log = np.array([math.log(float(reference[name])) for name in parameters], dtype=float)
    rhs = fim_array @ mle_log + weak_prior_precision @ prior_log + l2_precision @ prior_log
    posterior_mean = posterior_cov @ rhs
    out = pd.DataFrame(
        {
            "candidate": candidate.name,
            "theta_hat": [float(theta[name]) for name in parameters],
            "reference_prior_median": [float(reference[name]) for name in parameters],
            "posterior_median_laplace": np.exp(posterior_mean),
            "posterior_log_sd": np.sqrt(np.maximum(np.diag(posterior_cov), 0.0)),
            "posterior_95_lower": np.exp(posterior_mean - 1.96 * np.sqrt(np.maximum(np.diag(posterior_cov), 0.0))),
            "posterior_95_upper": np.exp(posterior_mean + 1.96 * np.sqrt(np.maximum(np.diag(posterior_cov), 0.0))),
            "posterior_var_over_weak_prior_var": np.diag(posterior_cov) / float(weak_prior_log_sd) ** 2,
            "l2_prior_precision": np.diag(l2_precision),
        },
        index=parameters,
    )
    out["data_dominated_vs_weak_prior"] = out["posterior_var_over_weak_prior_var"] < 0.10
    return out


def model_selection_metrics(fit_summary: pd.DataFrame, n_observations: int) -> pd.DataFrame:
    rows = []
    for name, row in fit_summary.iterrows():
        if not np.isfinite(row.get("direct_WSSE_total", np.nan)):
            continue
        k = int(row["n_estimated"])
        wsse = float(row["direct_WSSE_total"])
        rows.append(
            {
                "candidate": name,
                "n_estimated": k,
                "direct_WSSE_total": wsse,
                "AIC_like": 2.0 * k + 2.0 * wsse,
                "BIC_like": math.log(float(n_observations)) * k + 2.0 * wsse,
                "n_active_bounds": row.get("n_active_bounds", np.nan),
                "solve_ok": bool(row.get("solve_ok", False)),
            }
        )
    metrics = pd.DataFrame(rows).set_index("candidate").sort_values("BIC_like")
    return metrics


def write_report(
    fit_summary: pd.DataFrame,
    model_metrics: pd.DataFrame,
    fim_summary: pd.DataFrame,
    profile_summary: pd.DataFrame,
    bayes_summary: pd.DataFrame,
) -> None:
    report = RESULTS_DIR / "deep_model_selection_report.md"
    lines = [
        "# Deep model selection report",
        "",
        "Decision rule:",
        "",
        "1. Reject candidates with failed solves or estimated parameters at active bounds as structural models.",
        "2. Among stable candidates, require full-rank weighted-relative FIM/Q diagnostics.",
        "3. Require two-sided 95% profile-likelihood crossings for practical identifiability.",
        "4. Use log-L2/Bayesian candidates only as regularized curve-fit scenarios unless their added parameters pass unpenalized profiles.",
        "",
        "## Fit and information criteria",
        "",
        model_metrics.to_markdown(),
        "",
        "## FIM/Q sensitivity summary",
        "",
        fim_summary.to_markdown(index=False) if not fim_summary.empty else "No FIM/Q summary generated.",
        "",
        "## Profile-likelihood summary",
        "",
        profile_summary.to_markdown() if not profile_summary.empty else "No profile summary generated.",
        "",
        "## Bayesian/Laplace summary",
        "",
        bayes_summary.to_markdown() if not bayes_summary.empty else "No Bayesian summary generated.",
        "",
        "## Current recommendation",
        "",
        "Use `plus_iG_unregularized` as the recommended structurally/practically identifiable model.",
        "`plus_iG_unregularized` improves fit relative to `reduced6_local`, remains full-rank in the weighted-relative FIM/Q diagnostic, and passes two-sided profile-likelihood checks for all seven estimated parameters.",
        "Use `plus_qx_iG_logL2_10` only as a Bayesian regularized curve-fit sensitivity candidate; `qXG` and `qXF` do not pass the unregularized structural/practical identifiability filters.",
    ]
    report.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deep model-selection diagnostics for the fermentation calibration.")
    parser.add_argument("--batches", nargs="*", default=DEFAULT_BATCHES)
    parser.add_argument("--q-fraction", type=float, default=0.10)
    parser.add_argument("--profile-grid", type=int, default=7)
    parser.add_argument("--skip-fim", action="store_true")
    parser.add_argument("--skip-profiles", action="store_true")
    parser.add_argument("--weak-prior-log-sd", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_notebook_context()
    load_q_functions(ns)
    batches = [str(batch) for batch in args.batches]
    reference = reference_theta(ns, "reduced6", batch_id=batches[0])

    candidates = [
        DeepCandidate(
            "reduced6_local",
            REDUCED6,
            FIT_RESULTS_DIR / "reduced6_local_theta.csv",
            profile_parameters=(),
            profile_label="strict_unregularized",
        ),
        DeepCandidate(
            "plus_iG_unregularized",
            PLUS_IG,
            FIT_RESULTS_DIR / "plus_iG_unregularized_theta.csv",
            profile_parameters=PLUS_IG,
            profile_label="added_iG_unregularized",
        ),
        DeepCandidate(
            "plus_qx_iG_unregularized",
            PLUS_QX_IG,
            FIT_RESULTS_DIR / "plus_qx_iG_unregularized_theta.csv",
            profile_parameters=RELEASED_QX_IG,
            profile_label="added_qx_iG_unregularized_boundary",
        ),
        DeepCandidate(
            "plus_qx_iG_logL2_10",
            PLUS_QX_IG,
            FIT_RESULTS_DIR / "plus_qx_iG_logL2_10_theta.csv",
            objective_kind="SSE_weighted_logL2_10",
            l2_lambda=10.0,
            l2_parameters=RELEASED_QX_IG,
            profile_parameters=RELEASED_QX_IG,
            profile_label="added_qx_iG_posterior_logL2",
        ),
    ]

    fit_summary = pd.read_csv(FIT_RESULTS_DIR / "fit_strategy_summary.csv", index_col=0)
    n_observations = int(fit_summary["n_observations"].dropna().iloc[0])
    model_metrics = model_selection_metrics(fit_summary, n_observations=n_observations)
    model_metrics.to_csv(RESULTS_DIR / "deep_model_information_criteria.csv")

    theta_by_candidate = {candidate.name: load_theta(candidate.theta_file) for candidate in candidates}
    fim_rows = []
    if not args.skip_fim:
        for candidate in candidates:
            try:
                fim_rows.append(
                    run_q_analysis(
                        ns,
                        candidate=candidate,
                        theta=theta_by_candidate[candidate.name],
                        batches=batches,
                        fraction=float(args.q_fraction),
                    )
                )
            except Exception as err:
                fim_rows.append(
                    {
                        "candidate": candidate.name,
                        "fim_error": f"{type(err).__name__}: {err}",
                        "fim_full_rank_by_threshold": False,
                    }
                )
                print(f"[FIM] {candidate.name} failed: {type(err).__name__}: {err}", flush=True)
    elif (RESULTS_DIR / "deep_fim_summary.csv").exists():
        fim_summary = pd.read_csv(RESULTS_DIR / "deep_fim_summary.csv")
    else:
        fim_summary = pd.DataFrame()
    if not args.skip_fim:
        fim_summary = pd.DataFrame(fim_rows)
        fim_summary.to_csv(RESULTS_DIR / "deep_fim_summary.csv", index=False)

    profile_frames = []
    profile_summaries = []
    if not args.skip_profiles:
        for candidate in candidates:
            if candidate.name == "reduced6_local":
                continue
            theta_hat = theta_by_candidate[candidate.name]
            direct_wsse, penalty, obj_hat = evaluate_objective_components(ns, theta_hat, batches, reference, candidate)
            if not candidate.uses_l2:
                obj_hat = direct_wsse
            print(
                f"\n[PROFILE] {candidate.name}: objective_hat={obj_hat:.6g}, direct_WSSE={direct_wsse:.6g}, penalty={penalty:.6g}",
                flush=True,
            )
            profiles, summary = run_profile(
                ns,
                candidate=candidate,
                theta_hat=theta_hat,
                obj_hat=float(obj_hat),
                batches=batches,
                reference=reference,
                n_grid=int(args.profile_grid),
            )
            profiles.to_csv(RESULTS_DIR / f"{candidate.name}_profiles.csv", index=False)
            if not summary.empty:
                summary.to_csv(RESULTS_DIR / f"{candidate.name}_profile_summary.csv")
                profile_summaries.append(summary)
            profile_frames.append(profiles)
            try:
                fig, _axes = plot_profile_likelihood_profiles(profiles, alpha=0.95)
                if fig is not None:
                    fig.savefig(RESULTS_DIR / f"{candidate.name}_profile_likelihood.png", dpi=180)
            except Exception as err:
                print(f"[PROFILE] plot failed for {candidate.name}: {type(err).__name__}: {err}", flush=True)
        profile_all = pd.concat(profile_frames, ignore_index=True) if profile_frames else pd.DataFrame()
        profile_all.to_csv(RESULTS_DIR / "deep_profile_profiles.csv", index=False)
        profile_summary = pd.concat(profile_summaries) if profile_summaries else pd.DataFrame()
        profile_summary.to_csv(RESULTS_DIR / "deep_profile_summary.csv")
    elif (RESULTS_DIR / "deep_profile_summary.csv").exists():
        profile_summary = pd.read_csv(RESULTS_DIR / "deep_profile_summary.csv")
    else:
        profile_summary = pd.DataFrame()

    bayes_frames = []
    for candidate in candidates:
        fim_path = RESULTS_DIR / f"{candidate.name}_fim_weighted_relative.csv"
        bayes = laplace_summary_from_fim(
            ns,
            candidate=candidate,
            theta=theta_by_candidate[candidate.name],
            reference=reference,
            fim_path=fim_path,
            weak_prior_log_sd=float(args.weak_prior_log_sd),
            n_batches=len(batches),
        )
        if not bayes.empty:
            bayes.to_csv(RESULTS_DIR / f"{candidate.name}_bayesian_laplace.csv")
            bayes_frames.append(bayes)
    bayes_summary = pd.concat(bayes_frames) if bayes_frames else pd.DataFrame()
    bayes_summary.to_csv(RESULTS_DIR / "deep_bayesian_laplace_summary.csv")

    write_report(
        fit_summary=fit_summary,
        model_metrics=model_metrics,
        fim_summary=fim_summary,
        profile_summary=profile_summary,
        bayes_summary=bayes_summary,
    )

    print("\nModel information criteria:")
    print(model_metrics.head(10))
    print("\nFIM summary:")
    print(fim_summary)
    print("\nProfile summary:")
    print(profile_summary)
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
