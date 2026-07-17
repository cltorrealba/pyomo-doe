from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_MODEL_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from pilot_2026.adaptive_design.pilot_calibration import (  # noqa: E402
    censoring_contract,
    fit_co2,
    fit_primary,
    load_config,
    load_tables,
    validate_fit,
)
from pilot_2026.adaptive_design.resolve_historical_audit import (  # noqa: E402
    build_resolution,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


DEFAULT_CONFIG = ADAPTIVE_DIR / "calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
MODEL_RESULT_ROOT = RESULT_ROOT / "model_dataset"
INTEGRATION_RESULTS = PILOT_DIR / "results" / "data_integration_2026"


def latest_completed_model_run() -> Path:
    candidates: list[Path] = []
    for path in MODEL_RESULT_ROOT.glob("*/run_manifest.json"):
        payload = json.loads(path.read_text(encoding="utf-8"))
        adapter_path = path.parent / "adapter_config.json"
        if not adapter_path.exists():
            continue
        adapter = json.loads(adapter_path.read_text(encoding="utf-8"))
        if (
            payload.get("status") == "completed"
            and payload.get("gate", {}).get("verdict") == "PASS"
            and int(adapter.get("adapter_schema_version", 0)) >= 2
        ):
            candidates.append(path.parent)
    if not candidates:
        raise FileNotFoundError(
            "No completed active-only schema-v2 model dataset is available; run "
            "build_model_dataset.py first"
        )
    return max(candidates, key=lambda path: path.name)


def solver_capabilities() -> dict[str, Any]:
    capabilities: dict[str, Any] = {
        "calibration_backend": "scipy.optimize.least_squares_trf",
        "scipy_least_squares_available": False,
        "pyomo_doe_importable": False,
        "ipopt_available": False,
        "notes": [
            "IPOPT is reserved for the later local MBDoE refinement; the present "
            "bounded multistart calibration is solved with SciPy TRF."
        ],
        "errors": [],
    }
    try:
        from scipy.optimize import least_squares  # noqa: F401

        capabilities["scipy_least_squares_available"] = True
    except Exception as error:  # pragma: no cover
        capabilities["errors"].append(f"SciPy: {type(error).__name__}: {error}")
    try:
        import pyomo.environ as pyo
        from pyomo.contrib.doe import DesignOfExperiments  # noqa: F401

        capabilities["pyomo_doe_importable"] = True
        capabilities["ipopt_available"] = bool(
            pyo.SolverFactory("ipopt").available(exception_flag=False)
        )
    except Exception as error:
        capabilities["errors"].append(f"Pyomo/IPOPT: {type(error).__name__}: {error}")
    return capabilities


def evaluate_gate(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility evaluator for callers of the former fail-closed gate.

    The executable runner now evaluates scientific checks after the real fit.
    This helper preserves a small pure function for downstream tests and makes
    it explicit that owner review blocks physical release, not computational
    parameter estimation.
    """

    if args and len(args) >= 6:
        config, model_manifest, model_qc, _integration_qc, audit, capabilities = args[:6]
    else:
        config = kwargs["config"]
        model_manifest = kwargs["model_manifest"]
        model_qc = kwargs["model_qc"]
        audit = kwargs["audit"]
        capabilities = kwargs["capabilities"]
    checks = {
        "phase_A_model_dataset_pass": model_qc.get("verdict") == "PASS",
        "phase_A_manifest_complete": model_manifest.get("status") == "completed",
        "repository_raw_hash_audit_pass": bool(audit.get("pass")),
        "scipy_least_squares_available": bool(
            capabilities.get("scipy_least_squares_available", True)
        ),
        "existing_pilot_2025_integrated_result_not_used": not bool(
            config["pilot_2025_policy"]["use_existing_integrated_result_as_prior"]
        ),
    }
    blockers = [name for name, passed in checks.items() if not passed]
    return {
        "gate": "phase_B_computational_calibration_entry",
        "verdict": "PASS" if not blockers else "FAIL",
        "checks": checks,
        "blockers": blockers,
        "fit_executed": False,
        "profiles_for_physical_execution": False,
    }


def _write_frame(path: Path, frame, outputs: list[Path]) -> None:
    frame.to_csv(path, index=True if path.name == "primary_covariance_log_space.csv" else False)
    outputs.append(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run and validate reduced hierarchical Pilot 2026 calibration"
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--model-run", type=Path)
    parser.add_argument("--result-root", type=Path, default=RESULT_ROOT)
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = load_config(config_path)
    model_run = (args.model_run or latest_completed_model_run()).resolve()
    model_manifest_path = model_run / "run_manifest.json"
    model_qc_path = model_run / "dataset_qc.json"
    model_manifest = json.loads(model_manifest_path.read_text(encoding="utf-8"))
    model_qc = json.loads(model_qc_path.read_text(encoding="utf-8"))
    capabilities = solver_capabilities()
    audit_resolution = build_resolution()
    entry_audit = {
        "pass": audit_resolution["verdict"] == "PASS",
        "stdout": audit_resolution["current_audit"],
    }
    entry_gate = evaluate_gate(
        config,
        model_manifest,
        model_qc,
        {},
        entry_audit,
        capabilities,
    )
    if entry_gate["verdict"] != "PASS":
        raise SystemExit(
            "Calibration entry gate failed: " + ", ".join(entry_gate["blockers"])
        )

    run_dir = create_immutable_run_directory(
        Path(args.result_root), "baseline_calibration", config
    )
    outputs: list[Path] = []
    try:
        tables = load_tables(model_run)
        censoring = censoring_contract(tables, config)
        primary = fit_primary(tables, config)
        co2 = fit_co2(tables, primary, config)
        validation = validate_fit(primary, co2, audit_resolution, censoring, config)
        estimated = primary.parameter_table["parameter"].tolist() + co2.parameter_table[
            "parameter"
        ].tolist()
        gate = {
            "gate": "phase_B_hierarchical_calibration",
            "verdict": validation["release_verdict"],
            "computational_verdict": validation["computational_verdict"],
            "entry_gate": entry_gate,
            "fit_executed": True,
            "parameters_estimated": estimated,
            "residual_ESS_recomputed": True,
            "cooling_excluded_from_kinetics": True,
            "profiles_for_physical_execution": False,
            "physical_execution_authorized": False,
            "conditions": validation["conditions"],
        }

        json_outputs = {
            "calibration_config.json": config,
            "audit_resolution.json": audit_resolution,
            "solver_capabilities.json": capabilities,
            "censoring_contract.json": censoring,
            "calibration_validation.json": validation,
            "calibration_gate.json": gate,
        }
        for name, payload in json_outputs.items():
            path = run_dir / name
            write_json(path, payload)
            outputs.append(path)
        frame_outputs = {
            "primary_parameter_estimates.csv": primary.parameter_table,
            "primary_multistart_summary.csv": primary.multistart_summary,
            "primary_predictions.csv": primary.predictions,
            "primary_residual_summary.csv": primary.residual_summary,
            "primary_error_scale_estimates.csv": primary.error_scales,
            "primary_covariance_log_space.csv": primary.covariance,
            "co2_parameter_estimates.csv": co2.parameter_table,
            "co2_multistart_summary.csv": co2.multistart_summary,
            "co2_predictions.csv": co2.predictions,
            "co2_residual_ess.csv": co2.ess,
        }
        for name, frame in frame_outputs.items():
            _write_frame(run_dir / name, frame, outputs)
        ensemble = primary.multistart_summary.copy()
        ensemble_path = run_dir / "posterior_multistart_ensemble.csv"
        ensemble.to_csv(ensemble_path, index=False)
        outputs.append(ensemble_path)

        status = (
            "completed"
            if validation["computational_verdict"] == "PASS"
            else "validation_failed"
        )
        manifest = build_manifest(
            run_dir=run_dir,
            stage="baseline_calibration",
            config=config,
            sources={
                "calibration_config": config_path,
                "model_dataset_manifest": model_manifest_path,
                "model_dataset_qc": model_qc_path,
                "integration_qc": INTEGRATION_RESULTS / "qc_summary.json",
                "integration_config": PILOT_DIR / "data_integration_config.json",
                "raw_data_manifest": FERMENTATION_MODEL_DIR
                / "campaigns"
                / "raw_data_manifest.csv",
            },
            code_paths=[
                Path(__file__),
                ADAPTIVE_DIR / "pilot_calibration.py",
                ADAPTIVE_DIR / "resolve_historical_audit.py",
                ADAPTIVE_DIR / "build_model_dataset.py",
                ADAPTIVE_DIR / "run_artifacts.py",
                FERMENTATION_MODEL_DIR
                / "shared"
                / "run_new_must_glycerol_estimability_doe.py",
            ],
            random_seeds=[int(seed) for seed in config["random_seeds"]],
            status=status,
            convergence={
                "solver": capabilities["calibration_backend"],
                "primary_multistarts": primary.validation,
                "co2_multistarts": co2.validation,
            },
            gate=gate,
            outputs=outputs,
        )
        write_json(run_dir / "run_manifest.json", manifest)
    except Exception as error:
        failure = {
            "gate": "phase_B_hierarchical_calibration",
            "verdict": "FAIL",
            "fit_executed": False,
            "exception_type": type(error).__name__,
            "exception": str(error),
            "profiles_for_physical_execution": False,
        }
        write_json(run_dir / "calibration_gate.json", failure)
        manifest = build_manifest(
            run_dir=run_dir,
            stage="baseline_calibration",
            config=config,
            sources={
                "calibration_config": config_path,
                "model_dataset_manifest": model_manifest_path,
                "model_dataset_qc": model_qc_path,
            },
            code_paths=[Path(__file__), ADAPTIVE_DIR / "pilot_calibration.py"],
            random_seeds=[int(seed) for seed in config["random_seeds"]],
            status="failed",
            convergence={"exception_type": type(error).__name__, "exception": str(error)},
            gate=failure,
            outputs=[run_dir / "calibration_gate.json"],
        )
        write_json(run_dir / "run_manifest.json", manifest)
        raise

    print(
        json.dumps(
            {
                "verdict": gate["verdict"],
                "computational_verdict": gate["computational_verdict"],
                "run_directory": str(run_dir),
                "fit_executed": gate["fit_executed"],
                "parameters_estimated": len(gate["parameters_estimated"]),
                "residual_ESS_total": co2.validation["residual_ess_total"],
                "profiles_for_physical_execution": False,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    if validation["computational_verdict"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
