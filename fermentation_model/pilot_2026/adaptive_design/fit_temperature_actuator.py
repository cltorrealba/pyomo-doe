from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


MODEL_RUN = REPOSITORY_DIR / "fermentation_model/pilot_2026/results/adaptive_design_2026/model_dataset/20260717T125004Z_b03413"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
CONFIG = {
    "schema_version": 1,
    "model": "dT_executed_dt=(setpoint-T_executed)/tau",
    "minimum_setpoint_step_c": 0.5,
    "maximum_step_response_window_h": 12.0,
    "tau_bounds_h": [0.05, 12.0],
    "minimum_points_per_run": 30,
}


def fit_run(group: pd.DataFrame) -> dict:
    group = group.sort_values("time_h").drop_duplicates("time_h")
    time = group["time_h"].to_numpy(dtype=float)
    temperature = group["executed_temperature_c"].to_numpy(dtype=float)
    setpoint = group["commanded_setpoint_c"].to_numpy(dtype=float)
    events = np.flatnonzero(
        np.r_[False, np.abs(np.diff(setpoint)) >= float(CONFIG["minimum_setpoint_step_c"])]
    )
    responses = []
    for event_index, index in enumerate(events):
        next_change = time[events[event_index + 1]] if event_index + 1 < len(events) else np.inf
        end = min(
            time[index] + float(CONFIG["maximum_step_response_window_h"]), next_change
        )
        mask = (
            (time >= time[index])
            & (time <= end)
            & np.isfinite(temperature)
            & np.isfinite(setpoint)
        )
        if mask.sum() >= 4:
            responses.append(
                (time[mask] - time[index], temperature[mask], setpoint[index], temperature[index])
            )
    points = sum(len(response[0]) for response in responses)
    if points < int(CONFIG["minimum_points_per_run"]):
        return {
            "step_events": int(len(responses)),
            "points": int(points),
            "success": False,
            "tau_h": np.nan,
            "rmse_c": np.nan,
        }
    lower_tau, upper_tau = CONFIG["tau_bounds_h"]
    result = least_squares(
        lambda log_tau: np.concatenate(
            [
                observed
                - (new_setpoint + (initial - new_setpoint) * np.exp(-elapsed / np.exp(log_tau[0])))
                for elapsed, observed, new_setpoint, initial in responses
            ]
        ),
        np.asarray([np.log(0.75)]),
        bounds=([np.log(lower_tau)], [np.log(upper_tau)]),
        loss="soft_l1",
        f_scale=0.5,
    )
    tau = float(np.exp(result.x[0]))
    residual = np.concatenate(
        [
            observed
            - (new_setpoint + (initial - new_setpoint) * np.exp(-elapsed / tau))
            for elapsed, observed, new_setpoint, initial in responses
        ]
    )
    return {
        "step_events": int(len(responses)),
        "points": int(points),
        "success": bool(result.success),
        "tau_h": tau,
        "rmse_c": float(np.sqrt(np.mean(residual**2))),
    }


def main() -> None:
    temperature_path = MODEL_RUN / "temperature_inputs.csv"
    data = pd.read_csv(temperature_path)
    rows = []
    for run, group in data.groupby("experiment_id", sort=True):
        rows.append({"experiment_id": str(run), **fit_run(group)})
    estimates = pd.DataFrame(rows)
    valid = estimates[estimates["success"] & estimates["tau_h"].notna()]
    global_tau = float(np.median(valid["tau_h"])) if len(valid) else np.nan
    checks = {
        "all_nine_runs_evaluated": len(estimates) == 9,
        "at_least_six_run_fits_valid": len(valid) >= 6,
        "global_tau_finite": bool(np.isfinite(global_tau)),
        "global_tau_inside_bounds": bool(
            CONFIG["tau_bounds_h"][0] <= global_tau <= CONFIG["tau_bounds_h"][1]
        ),
        "median_step_response_rmse_below_1_2C": bool(valid["rmse_c"].median() <= 1.2),
    }
    gate = {
        "gate": "temperature_actuator_identification",
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "global_tau_h": global_tau,
        "interpretation": "working first-order setpoint-to-Sonda1 actuator model for MBDoE; uncertainty remains represented by run spread",
        "profiles_for_physical_execution": False,
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "temperature_actuator", CONFIG)
    estimates_path = run_dir / "actuator_estimates_by_run.csv"
    gate_path = run_dir / "actuator_gate.json"
    config_path = run_dir / "actuator_config.json"
    estimates.to_csv(estimates_path, index=False)
    write_json(gate_path, gate)
    write_json(config_path, CONFIG)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="temperature_actuator_identification",
        config=CONFIG,
        sources={"temperature_inputs": temperature_path, "model_manifest": MODEL_RUN / "run_manifest.json"},
        code_paths=[Path(__file__)],
        random_seeds=[],
        status="completed" if gate["verdict"] == "PASS" else "validation_failed",
        convergence={"valid_run_fits": int(len(valid))},
        gate=gate,
        outputs=[estimates_path, gate_path, config_path],
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    print(estimates.to_string(index=False))
    if gate["verdict"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
