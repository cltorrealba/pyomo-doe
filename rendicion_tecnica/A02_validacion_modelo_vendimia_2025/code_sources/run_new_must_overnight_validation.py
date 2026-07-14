from __future__ import annotations

import argparse
import math
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_new_must_glycerol_estimability_doe as base

OUT_DIR = SCRIPT_DIR / "results" / "new_must_glycerol_overnight_validation"
BASE_RESULTS = SCRIPT_DIR / "results" / "new_must_glycerol_estimability_doe"


def log(message: str) -> None:
    print(message, flush=True)


def load_current_theta() -> dict[str, float]:
    path = BASE_RESULTS / "theta_mixed_full17_l2.csv"
    if path.exists():
        theta = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
        return base.clip_theta({name: float(value) for name, value in theta.items()})
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    return base.complete_initial_theta(batches)


def theta_row(label: str, theta: dict[str, float], **extra) -> dict:
    row = {"label": label, **extra}
    row.update({name: float(theta[name]) for name in base.FULL17})
    return row


def random_theta_around(theta: dict[str, float], rng: np.random.Generator, sigma: float = 0.65) -> dict[str, float]:
    trial = dict(theta)
    for name in base.FULL17:
        value = max(float(theta[name]), 1e-16)
        trial[name] = value * math.exp(float(rng.normal(0.0, sigma)))
    return base.clip_theta(trial)


def run_l2_scan(batches, theta_ref, args) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    theta_rows = []
    for lam in args.l2_lambdas:
        label = f"full17_l2_lambda_{lam:g}"
        log(f"[l2] {label}")
        start = time.time()
        try:
            theta, summary = base.fit_parameters(
                label,
                batches,
                theta_ref,
                base.FULL17,
                max_nfev=args.l2_nfev,
                l2_reference=theta_ref,
                l2_parameters=tuple(name for name in base.FULL17 if name not in base.CORE_FIT),
                l2_lambda=float(lam),
            )
            summary["runtime_s"] = time.time() - start
            rows.append(summary)
            theta_rows.append(theta_row(label, theta, l2_lambda=float(lam)))
            pd.Series(theta, name=label).to_csv(OUT_DIR / f"theta_{label}.csv")
        except Exception as err:
            rows.append({"fit": label, "success": False, "error": f"{type(err).__name__}: {err}", "runtime_s": time.time() - start})
        pd.DataFrame(rows).to_csv(OUT_DIR / "l2_scan_summary.csv", index=False)
        pd.DataFrame(theta_rows).to_csv(OUT_DIR / "l2_scan_theta.csv", index=False)
    return pd.DataFrame(rows), pd.DataFrame(theta_rows)


def run_multistart(batches, theta_ref, args) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(args.seed)
    rows = []
    theta_rows = []
    for idx in range(int(args.multistarts)):
        label = f"multistart_{idx + 1:02d}"
        initial = theta_ref if idx == 0 else random_theta_around(theta_ref, rng, sigma=args.multistart_sigma)
        log(f"[multistart] {label}")
        start = time.time()
        try:
            theta, summary = base.fit_parameters(
                label,
                batches,
                initial,
                base.FULL17,
                max_nfev=args.multistart_nfev,
                l2_reference=theta_ref,
                l2_parameters=tuple(name for name in base.FULL17 if name not in base.CORE_FIT),
                l2_lambda=float(args.multistart_l2_lambda),
            )
            summary["runtime_s"] = time.time() - start
            rows.append(summary)
            theta_rows.append(theta_row(label, theta, run_id=idx + 1))
            pd.Series(theta, name=label).to_csv(OUT_DIR / f"theta_{label}.csv")
        except Exception as err:
            rows.append({"fit": label, "success": False, "error": f"{type(err).__name__}: {err}", "runtime_s": time.time() - start})
        pd.DataFrame(rows).to_csv(OUT_DIR / "multistart_summary.csv", index=False)
        pd.DataFrame(theta_rows).to_csv(OUT_DIR / "multistart_theta.csv", index=False)
    return pd.DataFrame(rows), pd.DataFrame(theta_rows)


def run_full_profiles(batches, theta_ref, args) -> tuple[pd.DataFrame, pd.DataFrame]:
    base_resid = base.residual_vector(theta_ref, batches)
    base_obj = float(np.dot(base_resid, base_resid))
    all_profiles = []
    summaries = []
    for parameter in base.FULL17:
        log(f"[profile] {parameter}")
        start = time.time()
        try:
            profile_df, summary_df = base.profile_parameters(
                theta_ref,
                batches,
                base.FULL17,
                (parameter,),
                base_obj,
                max_nfev=args.profile_nfev,
                grid_points=args.profile_grid,
            )
            profile_df["runtime_s_parameter"] = time.time() - start
            summary_df["runtime_s_parameter"] = time.time() - start
            profile_df.to_csv(OUT_DIR / f"profile_{parameter}.csv", index=False)
            all_profiles.append(profile_df)
            summaries.append(summary_df)
        except Exception as err:
            summaries.append(
                pd.DataFrame(
                    [
                        {
                            "parameter": parameter,
                            "n_success": 0,
                            "profile_identifiable_95": False,
                            "error": f"{type(err).__name__}: {err}",
                            "runtime_s_parameter": time.time() - start,
                        }
                    ]
                )
            )
        if all_profiles:
            pd.concat(all_profiles, ignore_index=True).to_csv(OUT_DIR / "profile_full_profiles.csv", index=False)
        if summaries:
            pd.concat(summaries, ignore_index=True).to_csv(OUT_DIR / "profile_full_summary.csv", index=False)
    profiles = pd.concat(all_profiles, ignore_index=True) if all_profiles else pd.DataFrame()
    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    return profiles, summary


def run_sampling_policy_benchmark(theta_ref, data, prior_fim, args) -> pd.DataFrame:
    rows = []
    designs = base.make_future_designs(data)
    design_by_name = {design.name: design for design in designs}
    for policy in ["front_loaded", "balanced", "two_per_day"]:
        log(f"[sampling] {policy}")
        candidate_fims = {}
        for design in designs:
            try:
                candidate_fims[design.name] = base.future_fim(theta_ref, design, base.FULL17, policy, args.sensitivity_step)
            except Exception as err:
                rows.append({"policy": policy, "candidate": design.name, "status": "failed", "error": f"{type(err).__name__}: {err}"})
        if not candidate_fims:
            continue
        selected, final_fim = base.greedy_campaign(
            candidate_fims,
            design_by_name,
            prior_fim,
            campaign_size=args.campaign_size,
            objective="hybrid",
        )
        selected["policy"] = policy
        selected.to_csv(OUT_DIR / f"selected_campaign_{policy}.csv", index=False)
        metrics = base.fim_metrics(final_fim, base.FULL17, prefix="final_")
        reductions = base.variance_reduction(prior_fim, final_fim, base.FULL17)
        rows.append({"policy": policy, "candidate": "__campaign__", "status": "ok", **metrics, **reductions})
    benchmark = pd.DataFrame(rows)
    benchmark.to_csv(OUT_DIR / "sampling_policy_benchmark.csv", index=False)
    return benchmark


def run_pyomo_selected(theta_ref, data, args) -> pd.DataFrame:
    selected_path = BASE_RESULTS / "selected_campaign_hybrid.csv"
    if not selected_path.exists():
        return pd.DataFrame()
    selected = pd.read_csv(selected_path)
    designs = {design.name: design for design in base.make_future_designs(data)}
    rows = []
    for name in selected["candidate"].head(int(args.pyomo_candidates)):
        if name not in designs:
            continue
        log(f"[pyomo] {name}")
        start = time.time()
        try:
            row = base.run_pyomo_doe_check(theta_ref, designs[name], base.REDUCED11, args.sample_policy, args.sensitivity_step)
            row["runtime_s"] = time.time() - start
            src = BASE_RESULTS / f"pyomo_doe_fim_{name}.csv"
            if src.exists():
                shutil.copy2(src, OUT_DIR / f"pyomo_doe_fim_{name}.csv")
        except Exception as err:
            row = {"candidate": name, "status": "failed", "error": f"{type(err).__name__}: {err}", "runtime_s": time.time() - start}
        rows.append(row)
        pd.DataFrame(rows).to_csv(OUT_DIR / "pyomo_selected_summary.csv", index=False)
    return pd.DataFrame(rows)


def write_report(l2_summary, multistart_summary, profile_summary, sampling_summary, pyomo_summary) -> None:
    lines = [
        "# Overnight validation report",
        "",
        "This run is a validation layer for the new-must glycerol calibration/DOE workflow. It does not replace the base run; it checks robustness.",
        "",
        "## L2 scan",
        "",
        l2_summary.to_markdown(index=False) if not l2_summary.empty else "_Not run._",
        "",
        "## Multistart",
        "",
        multistart_summary.sort_values("final_wsse").head(20).to_markdown(index=False)
        if not multistart_summary.empty and "final_wsse" in multistart_summary.columns
        else "_Not run._",
        "",
        "## Full profile likelihood",
        "",
        profile_summary.to_markdown(index=False) if not profile_summary.empty else "_Not run._",
        "",
        "## Sampling policy benchmark",
        "",
        sampling_summary.to_markdown(index=False) if not sampling_summary.empty else "_Not run._",
        "",
        "## Pyomo.DoE selected candidates",
        "",
        pyomo_summary.to_markdown(index=False) if not pyomo_summary.empty else "_Not run._",
    ]
    (OUT_DIR / "overnight_validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Long validation run for new-must glycerol calibration and DOE.")
    parser.add_argument("--seed", type=int, default=90609)
    parser.add_argument("--l2-lambdas", type=float, nargs="*", default=[0.1, 1.0, 10.0, 50.0])
    parser.add_argument("--l2-nfev", type=int, default=120)
    parser.add_argument("--multistarts", type=int, default=10)
    parser.add_argument("--multistart-nfev", type=int, default=120)
    parser.add_argument("--multistart-l2-lambda", type=float, default=1.0)
    parser.add_argument("--multistart-sigma", type=float, default=0.65)
    parser.add_argument("--profile-nfev", type=int, default=45)
    parser.add_argument("--profile-grid", type=int, default=7)
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--sample-policy", choices=["front_loaded", "balanced", "two_per_day"], default="front_loaded")
    parser.add_argument("--sensitivity-step", type=float, default=1e-2)
    parser.add_argument("--pyomo-candidates", type=int, default=9)
    parser.add_argument("--skip-l2", action="store_true")
    parser.add_argument("--skip-multistart", action="store_true")
    parser.add_argument("--skip-profiles", action="store_true")
    parser.add_argument("--skip-sampling", action="store_true")
    parser.add_argument("--skip-pyomo", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"[start] overnight validation results: {OUT_DIR}")
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    theta_ref = load_current_theta()
    pd.Series(theta_ref, name="theta_reference").to_csv(OUT_DIR / "theta_reference.csv")
    prior_path = BASE_RESULTS / "fim_mixed_current_full_l2.csv"
    if prior_path.exists():
        prior_fim = pd.read_csv(prior_path, index_col=0).loc[list(base.FULL17), list(base.FULL17)].to_numpy(dtype=float)
    else:
        jac, _resid = base.build_jacobian(theta_ref, base.FULL17, batches, step=args.sensitivity_step)
        prior_fim = jac.T @ jac

    l2_summary = pd.DataFrame()
    multistart_summary = pd.DataFrame()
    profile_summary = pd.DataFrame()
    sampling_summary = pd.DataFrame()
    pyomo_summary = pd.DataFrame()

    if not args.skip_l2:
        l2_summary, _theta_rows = run_l2_scan(batches, theta_ref, args)
    if not args.skip_multistart:
        multistart_summary, _theta_rows = run_multistart(batches, theta_ref, args)
    if not args.skip_profiles:
        _profile_df, profile_summary = run_full_profiles(batches, theta_ref, args)
    if not args.skip_sampling:
        sampling_summary = run_sampling_policy_benchmark(theta_ref, data, prior_fim, args)
    if not args.skip_pyomo:
        pyomo_summary = run_pyomo_selected(theta_ref, data, args)

    write_report(l2_summary, multistart_summary, profile_summary, sampling_summary, pyomo_summary)
    log("[done] overnight validation complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
