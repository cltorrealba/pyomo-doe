"""Repeat profile likelihood using the best rebase theta as the baseline."""

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
        log(obj.to_string(max_rows=30, max_cols=14))
    elif isinstance(obj, pd.Series):
        log(obj.to_string())
    else:
        log(repr(obj))


def exec_cell(nb: dict, cell_index: int, ns: dict) -> None:
    src = "".join(nb["cells"][cell_index].get("source", []))
    code = compile(src, f"{NOTEBOOK.name}:cell_{cell_index}", "exec")
    exec(code, ns)


def load_rebase_theta() -> tuple[dict[str, float], float, dict]:
    theta_path = OUT_DIR / "profile_rebase_best_theta.csv"
    summary_path = OUT_DIR / "profile_rebase_summary.json"
    if not theta_path.exists():
        raise FileNotFoundError(f"Missing {theta_path}. Run profile rebase first.")
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing {summary_path}. Run profile rebase first.")

    theta_series = pd.read_csv(theta_path, index_col=0).iloc[:, 0]
    theta_hat = {name: float(value) for name, value in theta_series.items()}
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    obj_hat = float(summary["best_rebase_objective"])
    return theta_hat, obj_hat, summary


def main() -> None:
    os.chdir(REPO_DIR)
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    ns: dict = {"__name__": "__main__", "display": display}

    log("Loading notebook definitions...")
    for idx in [1, 3, 5, 8, 10, 12, 21]:
        exec_cell(nb, idx, ns)

    if str(FERMENTATION_DIR) not in sys.path:
        sys.path.insert(0, str(FERMENTATION_DIR))
    from estimability_tools import build_profile_grid, summarize_profile_likelihood

    ns["build_profile_grid"] = build_profile_grid
    ns["summarize_profile_likelihood"] = summarize_profile_likelihood
    ns["plot_profile_likelihood_profiles"] = lambda *args, **kwargs: (None, None)
    ns["parameter_uncertainty_method_comparison"] = pd.DataFrame()
    ns["multistart_best_theta"] = pd.DataFrame()

    ns["UQ_EXCLUDED_BATCHES_FOR_COVARIANCE"] = {"25085", "25171"}
    ns["UQ_COVARIANCE_BATCHES"] = [
        batch_id
        for batch_id in ns["PARAMETER_ESTIMATION_BATCHES"]
        if batch_id not in ns["UQ_EXCLUDED_BATCHES_FOR_COVARIANCE"]
    ]

    theta_hat, obj_hat, rebase_summary = load_rebase_theta()
    log("Rebased profile baseline:")
    display(pd.Series({"obj_hat": obj_hat, **theta_hat}))
    log(f"Rebase source: {rebase_summary.get('best_rebase_candidate_id')}")

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
    profile_batches = ns["UQ_COVARIANCE_BATCHES"]
    profile_objective = "SSE_weighted"
    n_grid = 9
    baseline_tol = 1e-6

    log(f"Profile batches: {profile_batches}")
    log(f"Running rebased profile likelihood for {len(profile_parameters)} parameters...")
    tables = []
    for parameter in profile_parameters:
        log(f"  profile {parameter}...")
        tables.append(
            ns["run_profile_likelihood_parameter"](
                parameter,
                theta_hat,
                obj_hat,
                obj_function=profile_objective,
                batch_ids=profile_batches,
                n_grid=n_grid,
            )
        )

    profiles = pd.concat(tables, ignore_index=True)
    profiles["lr_stat_rebased"] = 2.0 * (profiles["objective"] - obj_hat)
    profiles.to_csv(OUT_DIR / "profile_from_rebase_profiles.csv", index=False)

    ok = profiles[profiles["success"]].copy()
    ok["objective"] = pd.to_numeric(ok["objective"], errors="coerce")
    min_row = ok.loc[ok["objective"].idxmin()]
    negative = ok[ok["objective"] < obj_hat - baseline_tol].copy()
    negative["objective_improvement"] = obj_hat - negative["objective"]
    negative["lr_stat_vs_rebased_obj_hat"] = 2.0 * (negative["objective"] - obj_hat)
    negative = negative.sort_values(["objective", "profiled_theta", "theta_value"]).reset_index(drop=True)
    negative.to_csv(OUT_DIR / "profile_from_rebase_negative_rows.csv", index=False)

    summary = {
        "obj_hat": float(obj_hat),
        "n_profile_rows": int(len(profiles)),
        "n_successful_profile_rows": int(len(ok)),
        "n_negative_profile_rows": int(len(negative)),
        "min_profile_objective": float(min_row["objective"]),
        "min_profiled_theta": str(min_row["profiled_theta"]),
        "min_theta_value": float(min_row["theta_value"]),
        "best_improvement_vs_rebased_obj_hat": float(obj_hat - min_row["objective"]),
    }
    (OUT_DIR / "profile_from_rebase_diagnostic.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    profile_summary = ns["summarize_profile_likelihood"](profiles, alpha=0.95)
    profile_summary.reset_index().rename(columns={"index": "profiled_theta"}).to_csv(
        OUT_DIR / "profile_from_rebase_summary.csv",
        index=False,
    )

    log("Rebased profile diagnostic:")
    display(pd.Series(summary))
    log("Minimum objective by profiled parameter:")
    display(ok.groupby("profiled_theta")["objective"].min().sort_values())
    if negative.empty:
        log("No negative LR rows remain relative to the rebased baseline.")
    else:
        log("Negative LR rows remain relative to the rebased baseline:")
        display(
            negative[
                [
                    "profiled_theta",
                    "theta_value",
                    "objective",
                    "objective_improvement",
                    "lr_stat_vs_rebased_obj_hat",
                ]
            ]
        )
    log("Profile-likelihood summary:")
    display(profile_summary)
    log(f"Results written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
