"""Run profile-negative rebase analysis outside the notebook.

This script intentionally avoids the full notebook execution path. It loads the
model/data definitions from the notebook, recomputes the active profile
screening, builds full-theta candidates from negative profile rows, and
reoptimizes the base profile problem from those candidates.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd


REPO_DIR = Path(__file__).resolve().parents[1]
FERMENTATION_DIR = REPO_DIR / "fermentation_model"
NOTEBOOK = FERMENTATION_DIR / "fermentation_model_calibration.ipynb"
OUT_DIR = FERMENTATION_DIR / "results"
OUT_DIR.mkdir(exist_ok=True)


def log(message: str) -> None:
    print(message, flush=True)


def display(obj=None, *args, **kwargs) -> None:
    if obj is None:
        return
    if isinstance(obj, pd.DataFrame):
        log(obj.to_string(max_rows=20, max_cols=12))
    elif isinstance(obj, pd.Series):
        log(obj.to_string())
    else:
        log(repr(obj))


def exec_cell(nb: dict, cell_index: int, ns: dict, trim_after: str | None = None) -> None:
    src = "".join(nb["cells"][cell_index].get("source", []))
    if trim_after is not None:
        src = src.split(trim_after)[0]
    code = compile(src, f"{NOTEBOOK.name}:cell_{cell_index}", "exec")
    exec(code, ns)


def build_profile_rebase_summary(profile_df: pd.DataFrame, obj_hat: float, tol: float) -> pd.DataFrame:
    ok = profile_df[profile_df["success"]].copy()
    ok["objective"] = pd.to_numeric(ok["objective"], errors="coerce")
    ok = ok[ok["objective"].notna()].copy()
    ok["objective_improvement"] = float(obj_hat) - ok["objective"]
    ok["lr_stat_vs_obj_hat"] = 2.0 * (ok["objective"] - float(obj_hat))
    negative = ok[ok["objective"] < float(obj_hat) - float(tol)].copy()
    return negative.sort_values(["objective", "profiled_theta", "theta_value"]).reset_index(drop=True)


def main() -> None:
    os.chdir(REPO_DIR)
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))

    ns: dict = {
        "__name__": "__main__",
        "display": display,
    }

    log("Loading notebook definitions...")
    for idx in [1, 3, 5, 8, 10, 12, 21]:
        exec_cell(nb, idx, ns)

    if str(FERMENTATION_DIR) not in sys.path:
        sys.path.insert(0, str(FERMENTATION_DIR))
    from estimability_tools import build_profile_grid, summarize_profile_likelihood

    ns["build_profile_grid"] = build_profile_grid
    ns["summarize_profile_likelihood"] = summarize_profile_likelihood
    ns["plot_profile_likelihood_profiles"] = lambda *args, **kwargs: (None, None)

    # Avoid full UQ and full-data multistart here. These are only needed for
    # profile grid heuristics and warm-start selection in the notebook.
    ns["UQ_EXCLUDED_BATCHES_FOR_COVARIANCE"] = {"25085", "25171"}
    ns["UQ_COVARIANCE_BATCHES"] = [
        batch_id
        for batch_id in ns["PARAMETER_ESTIMATION_BATCHES"]
        if batch_id not in ns["UQ_EXCLUDED_BATCHES_FOR_COVARIANCE"]
    ]
    ns["parameter_uncertainty_method_comparison"] = pd.DataFrame()
    ns["multistart_best_theta"] = pd.DataFrame()

    # Last known full-data WSSE estimate, used only as a starting point for the
    # base profile-batch estimation. The base estimation is rerun below.
    ns["theta_WSSE"] = {
        "mu0": 1.0,
        "betaG0": 1.769680,
        "betaF0": 2.0,
        "Kn0": 1.083862,
        "Kg0": 0.001,
        "Kf0": 0.001,
        "Kig0": 45.255144,
        "Kie0": 6.143887,
        "Kd0": 0.001848,
        "Yxn": 7.323725,
        "Yxg": 10.0,
        "Yxf": 0.065350,
        "Yeg": 0.604107,
        "Yef": 1.669475,
    }

    log(f"Profile batches: {ns['UQ_COVARIANCE_BATCHES']}")
    log("Running base profile-batch WSSE estimation...")
    profile_obj_hat, profile_theta_raw, _ = ns["run_parmest_estimation"](
        "SSE_weighted",
        batch_ids=ns["UQ_COVARIANCE_BATCHES"],
        theta_initial=ns["theta_WSSE"],
    )
    profile_theta_hat, profile_theta_adjustments = ns["clip_theta_to_bounds"](profile_theta_raw)
    log(f"Base profile objective: {profile_obj_hat:.12g}")
    if not profile_theta_adjustments.empty:
        log("Base profile bound adjustments:")
        display(profile_theta_adjustments)

    # Define the same profile helper functions from the notebook.
    profile_src = "".join(nb["cells"][34].get("source", []))
    profile_defs = profile_src.split('if RUN_PROFILE_OPTIMIZATION or "profile_results" not in globals():')[0]
    exec(compile(profile_defs, f"{NOTEBOOK.name}:profile_defs", "exec"), ns)

    profile_parameters = [
        "mu0",
        "betaF0",
        "Kg0",
        "Kf0",
        "Kig0",
        "Yxg",
        "Yxf",
        "betaG0",
        "Yeg",
        "Yef",
    ]
    ns["PROFILE_PARAMETERS"] = profile_parameters
    ns["PROFILE_N_GRID"] = 9
    ns["PROFILE_OBJECTIVE"] = "SSE_weighted"
    ns["PROFILE_BATCHES"] = ns["UQ_COVARIANCE_BATCHES"]
    ns["PROFILE_BASELINE_TOL"] = 1e-6

    log(f"Running profile screening for {len(profile_parameters)} parameters...")
    profile_tables = []
    for parameter in profile_parameters:
        log(f"  profile {parameter}...")
        profile_tables.append(
            ns["run_profile_likelihood_parameter"](
                parameter,
                profile_theta_hat,
                profile_obj_hat,
                obj_function="SSE_weighted",
                batch_ids=ns["UQ_COVARIANCE_BATCHES"],
                n_grid=ns["PROFILE_N_GRID"],
            )
        )
    profile_profiles = pd.concat(profile_tables, ignore_index=True)
    profile_profiles.to_csv(OUT_DIR / "profile_profiles_recomputed.csv", index=False)

    negative = build_profile_rebase_summary(
        profile_profiles,
        profile_obj_hat,
        tol=ns["PROFILE_BASELINE_TOL"],
    )
    negative.to_csv(OUT_DIR / "profile_negative_rows.csv", index=False)
    log(f"Negative profile rows: {len(negative)}")
    if negative.empty:
        log("No negative profile rows found; no rebase estimation needed.")
        return
    display(
        negative[
            [
                "profiled_theta",
                "theta_value",
                "objective",
                "objective_improvement",
                "lr_stat_vs_obj_hat",
            ]
        ]
    )

    # Define candidate builders from the notebook rebase cell.
    rebase_src = "".join(nb["cells"][36].get("source", []))
    rebase_defs = rebase_src.split("profile_rebase_candidates,")[0]
    exec(compile(rebase_defs, f"{NOTEBOOK.name}:rebase_defs", "exec"), ns)
    ns["profile_results"] = {
        "profiles": profile_profiles,
        "theta_hat": profile_theta_hat,
        "obj_hat": float(profile_obj_hat),
    }
    ns["profile_profiles"] = profile_profiles

    candidates, candidate_summary = ns["collect_profile_rebase_candidates"](
        profile_profiles,
        profile_obj_hat,
        tol=ns["PROFILE_BASELINE_TOL"],
        max_candidates=None,
    )
    candidate_summary.to_csv(OUT_DIR / "profile_rebase_candidate_summary.csv", index=False)
    log(f"Rebase candidates: {len(candidates)}")

    log("Running base reoptimization from every negative-profile full-theta candidate...")
    rows = []
    for idx, candidate in enumerate(candidates, start=1):
        log(
            f"  rebase {idx}/{len(candidates)} from "
            f"{candidate['source_profiled_theta']}={candidate['source_theta_value']:.8g} "
            f"(source obj={candidate['source_objective']:.8g})"
        )
        row = {
            "candidate_id": candidate["candidate_id"],
            "source_profiled_theta": candidate["source_profiled_theta"],
            "source_theta_value": candidate["source_theta_value"],
            "source_objective": candidate["source_objective"],
            "source_lr_stat": candidate["source_lr_stat"],
            "success": False,
            "objective": np.nan,
            "error": "",
        }
        try:
            obj_value, theta_raw, _ = ns["run_parmest_estimation"](
                "SSE_weighted",
                batch_ids=ns["UQ_COVARIANCE_BATCHES"],
                theta_initial=candidate["theta_initial"],
                parameters_to_estimate=list(ns["DEFAULT_THETA"]),
            )
            theta_fit, adjustments = ns["clip_theta_to_bounds"](theta_raw)
            row["success"] = True
            row["objective"] = float(obj_value)
            row["n_bound_adjustments"] = len(adjustments)
            row.update({f"estimated_{name}": theta_fit[name] for name in ns["DEFAULT_THETA"]})
            log(f"    success objective={float(obj_value):.12g}")
        except Exception as err:  # noqa: BLE001 - diagnostic run should continue
            row["error"] = f"{type(err).__name__}: {err}"
            log(f"    failed {row['error']}")
        rows.append(row)

    rebase_results = pd.DataFrame(rows)
    rebase_results.to_csv(OUT_DIR / "profile_rebase_results.csv", index=False)
    successful = rebase_results[rebase_results["success"]].copy()
    summary = {
        "base_profile_objective": float(profile_obj_hat),
        "n_profile_rows": int(len(profile_profiles)),
        "n_negative_profile_rows": int(len(negative)),
        "n_rebase_candidates": int(len(candidates)),
        "n_successful_rebase": int(len(successful)),
    }
    if not successful.empty:
        best = successful.loc[successful["objective"].idxmin()]
        best_theta = {
            name: float(best[f"estimated_{name}"])
            for name in ns["DEFAULT_THETA"]
        }
        pd.Series(best_theta, name="profile_rebase_best_theta").to_csv(
            OUT_DIR / "profile_rebase_best_theta.csv"
        )
        summary.update(
            {
                "best_rebase_objective": float(best["objective"]),
                "best_rebase_candidate_id": best["candidate_id"],
                "best_source_profiled_theta": best["source_profiled_theta"],
                "best_source_theta_value": float(best["source_theta_value"]),
                "best_improvement_vs_base": float(profile_obj_hat - best["objective"]),
            }
        )
        log("Best rebase result:")
        display(pd.Series(summary))
        log("Best theta:")
        display(pd.Series(best_theta))
    else:
        log("No successful rebase solves.")
        display(pd.Series(summary))

    (OUT_DIR / "profile_rebase_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    log(f"Results written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
