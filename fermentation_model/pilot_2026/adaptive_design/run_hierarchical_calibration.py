from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_MODEL_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


DEFAULT_CONFIG = ADAPTIVE_DIR / "calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
MODEL_RESULT_ROOT = RESULT_ROOT / "model_dataset"
INTEGRATION_RESULTS = PILOT_DIR / "results" / "data_integration_2026"
AUDIT_SCRIPT = FERMENTATION_MODEL_DIR / "tools" / "campaign_audit.py"


def latest_completed_model_run() -> Path:
    candidates: list[Path] = []
    for path in MODEL_RESULT_ROOT.glob("*/run_manifest.json"):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") == "completed" and payload.get("gate", {}).get(
            "verdict"
        ) == "PASS":
            candidates.append(path.parent)
    if not candidates:
        raise FileNotFoundError("No completed PASS model-dataset run is available")
    return max(candidates, key=lambda path: path.name)


def run_repository_audit() -> dict[str, Any]:
    result = subprocess.run(
        [sys.executable, str(AUDIT_SCRIPT), "--check-hashes"],
        cwd=REPOSITORY_DIR,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        payload = {"errors": ["campaign audit did not emit valid JSON"]}
    return {
        "command": "python fermentation_model/tools/campaign_audit.py --check-hashes",
        "exit_code": int(result.returncode),
        "stdout": payload,
        "stderr": result.stderr.strip(),
        "pass": result.returncode == 0 and not payload.get("errors"),
    }


def solver_capabilities() -> dict[str, Any]:
    capabilities: dict[str, Any] = {
        "pyomo_doe_importable": False,
        "ipopt_available": False,
        "errors": [],
    }
    try:
        import pyomo.environ as pyo
        from pyomo.contrib.doe import DesignOfExperiments  # noqa: F401

        capabilities["pyomo_doe_importable"] = True
        capabilities["ipopt_available"] = bool(
            pyo.SolverFactory("ipopt").available(exception_flag=False)
        )
    except Exception as error:  # capability is reported, never hidden
        capabilities["errors"].append(f"{type(error).__name__}: {error}")
    return capabilities


def evaluate_gate(
    config: dict[str, Any],
    model_manifest: dict[str, Any],
    model_qc: dict[str, Any],
    integration_qc: dict[str, Any],
    audit: dict[str, Any],
    capabilities: dict[str, Any],
) -> dict[str, Any]:
    owner = config["owner_prerequisites"]
    checks: dict[str, bool] = {
        "phase_A_model_dataset_pass": model_qc.get("verdict") == "PASS",
        "phase_A_manifest_complete": model_manifest.get("status") == "completed",
        "integration_has_no_pending_calibration_items": not bool(
            integration_qc.get("calibration_gate")
        ),
        "repository_raw_hash_audit_pass": bool(audit["pass"]),
        "pyomo_doe_importable": bool(capabilities["pyomo_doe_importable"]),
        "ipopt_available": bool(capabilities["ipopt_available"]),
        "existing_pilot_2025_integrated_result_not_used": not bool(
            config["pilot_2025_policy"]["use_existing_integrated_result_as_prior"]
        ),
    }
    checks.update({f"owner_{name}": bool(value) for name, value in owner.items()})
    blockers = [name for name, passed in checks.items() if not passed]
    return {
        "gate": "phase_B_hierarchical_calibration",
        "verdict": "PASS" if not blockers else "FAIL",
        "checks": checks,
        "blockers": blockers,
        "integration_calibration_items": integration_qc.get("calibration_gate", []),
        "raw_audit_error_count": len(audit["stdout"].get("errors", [])),
        "raw_audit_warning_count": len(audit["stdout"].get("warnings", [])),
        "fit_executed": False,
        "parameters_estimated": [],
        "profiles_for_physical_execution": False,
        "scientific_interpretation": (
            "Calibration was not attempted because all upstream scientific, data-integrity, "
            "unit-conversion, solver, and owner gates must pass first."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate and, only after all gates pass, run Pilot 2026 calibration"
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--model-run", type=Path)
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    model_run = (args.model_run or latest_completed_model_run()).resolve()
    model_manifest_path = model_run / "run_manifest.json"
    model_qc_path = model_run / "dataset_qc.json"
    integration_qc_path = INTEGRATION_RESULTS / "qc_summary.json"
    model_manifest = json.loads(model_manifest_path.read_text(encoding="utf-8"))
    model_qc = json.loads(model_qc_path.read_text(encoding="utf-8"))
    integration_qc = json.loads(integration_qc_path.read_text(encoding="utf-8"))

    run_dir = create_immutable_run_directory(RESULT_ROOT, "baseline_calibration", config)
    audit = run_repository_audit()
    capabilities = solver_capabilities()
    gate = evaluate_gate(
        config, model_manifest, model_qc, integration_qc, audit, capabilities
    )

    copied_config = run_dir / "calibration_config.json"
    audit_path = run_dir / "repository_audit.json"
    capabilities_path = run_dir / "solver_capabilities.json"
    plan_path = run_dir / "planned_calibration.json"
    gate_path = run_dir / "calibration_gate.json"
    write_json(copied_config, config)
    write_json(audit_path, audit)
    write_json(capabilities_path, capabilities)
    write_json(
        plan_path,
        {
            "campaign_hierarchy": config["campaign_hierarchy"],
            "planned_stages": config["planned_stages"],
            "nuisance_effects": config["nuisance_effects"],
            "fim_nuisance_marginalization": config["fim_nuisance_marginalization"],
            "stop_rules": config["stop_rules"],
        },
    )
    write_json(gate_path, gate)
    outputs = [copied_config, audit_path, capabilities_path, plan_path, gate_path]
    sources = {
        "calibration_config": config_path,
        "model_dataset_manifest": model_manifest_path,
        "model_dataset_qc": model_qc_path,
        "integration_qc": integration_qc_path,
        "integration_config": PILOT_DIR / "data_integration_config.json",
        "raw_data_manifest": FERMENTATION_MODEL_DIR
        / "campaigns"
        / "raw_data_manifest.csv",
    }
    status = "gate_passed" if gate["verdict"] == "PASS" else "gate_failed"
    manifest = build_manifest(
        run_dir=run_dir,
        stage="baseline_calibration",
        config=config,
        sources=sources,
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "run_artifacts.py",
            FERMENTATION_MODEL_DIR / "shared" / "run_new_must_glycerol_estimability_doe.py",
            FERMENTATION_MODEL_DIR / "shared" / "run_secondary_joint_campaign_doe.py",
            FERMENTATION_MODEL_DIR / "shared" / "aroma_partition_unifac.py",
        ],
        random_seeds=[int(seed) for seed in config["random_seeds"]],
        status=status,
        convergence={
            "solver_invoked": False,
            "multistarts_completed": 0,
            "reason": "upstream_gate" if gate["verdict"] != "PASS" else "ready",
        },
        gate=gate,
        outputs=outputs,
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(
        json.dumps(
            {
                "verdict": gate["verdict"],
                "run_directory": str(run_dir),
                "blockers": gate["blockers"],
                "solver_capabilities": capabilities,
                "fit_executed": False,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    if gate["verdict"] != "PASS":
        raise SystemExit(2)
    raise SystemExit(
        "All prerequisites passed, but calibration execution requires a reviewed "
        "configuration update and explicit rerun; no implicit continuation is allowed."
    )


if __name__ == "__main__":
    main()
