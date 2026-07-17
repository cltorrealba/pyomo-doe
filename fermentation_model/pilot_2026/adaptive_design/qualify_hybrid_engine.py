from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.hybrid_optimizer import particle_swarm  # noqa: E402
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "hybrid_engine_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def rosenbrock(values: np.ndarray) -> float:
    x, y = np.asarray(values, dtype=float)
    return float((1.0 - x) ** 2 + 100.0 * (y - x * x) ** 2)


def _ipopt_capability() -> dict:
    executable = shutil.which("ipopt")
    try:
        import pyomo  # noqa: F401
        pyomo_available = True
    except ImportError:
        pyomo_available = False
    try:
        import cyipopt  # noqa: F401
        cyipopt_available = True
    except ImportError:
        cyipopt_available = False
    return {
        "ipopt_executable": executable,
        "pyomo_available": pyomo_available,
        "cyipopt_available": cyipopt_available,
        "usable_ipopt_interface": bool(cyipopt_available or (executable and pyomo_available)),
    }


def _run_ipopt_benchmark(start: np.ndarray, bounds: np.ndarray, capability: dict) -> dict:
    if capability["cyipopt_available"]:
        from cyipopt import minimize_ipopt

        result = minimize_ipopt(
            rosenbrock,
            np.asarray(start, dtype=float),
            bounds=[tuple(row) for row in bounds],
            options={"tol": 1e-10, "max_iter": 500, "print_level": 0},
        )
        return {
            "executed": True,
            "interface": "cyipopt.minimize_ipopt",
            "success": bool(result.success),
            "best_x": np.asarray(result.x, dtype=float).tolist(),
            "best_objective": float(result.fun),
            "message": str(result.message),
        }
    if capability["pyomo_available"] and capability["ipopt_executable"]:
        import pyomo.environ as pyo

        model = pyo.ConcreteModel()
        model.x = pyo.Var(bounds=tuple(bounds[0]), initialize=float(start[0]))
        model.y = pyo.Var(bounds=tuple(bounds[1]), initialize=float(start[1]))
        model.objective = pyo.Objective(
            expr=(1.0 - model.x) ** 2 + 100.0 * (model.y - model.x**2) ** 2
        )
        solver = pyo.SolverFactory("ipopt", executable=capability["ipopt_executable"])
        result = solver.solve(model, tee=False, options={"tol": 1e-10, "max_iter": 500})
        termination = str(result.solver.termination_condition)
        return {
            "executed": True,
            "interface": "pyomo.SolverFactory(ipopt)",
            "success": termination.lower() in {"optimal", "locallyoptimal", "globallyoptimal"},
            "best_x": [float(pyo.value(model.x)), float(pyo.value(model.y))],
            "best_objective": float(pyo.value(model.objective)),
            "message": termination,
        }
    return {
        "executed": False,
        "interface": None,
        "success": False,
        "best_x": None,
        "best_objective": None,
        "message": "No usable IPOPT interface",
    }


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    bounds = np.asarray(config["qualification"]["bounds"], dtype=float)
    global_config = config["global"]
    kwargs = {
        "particles": int(global_config["particles"]),
        "iterations": int(global_config["iterations"]),
        "seed": int(config["seed"]),
        "inertia": float(global_config["inertia"]),
        "cognitive": float(global_config["cognitive"]),
        "social": float(global_config["social"]),
    }
    first = particle_swarm(rosenbrock, bounds, **kwargs)
    second = particle_swarm(rosenbrock, bounds, **kwargs)
    # Diagnostic only: proves that the PSO seed is suitable for local refinement.
    scipy_local = minimize(
        rosenbrock,
        first.x,
        method="L-BFGS-B",
        bounds=[tuple(row) for row in bounds],
        options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": 500},
    )
    capability = _ipopt_capability()
    ipopt_local = _run_ipopt_benchmark(first.x, bounds, capability)
    repeat_difference = float(np.max(np.abs(first.x - second.x)))
    checks = {
        "pso_reproducible": repeat_difference
        <= float(config["qualification"]["maximum_repeat_coordinate_difference"]),
        "pso_improves_benchmark": bool(first.history[-1] < first.history[0]),
        "diagnostic_local_refinement_reaches_benchmark": float(scipy_local.fun)
        <= float(config["qualification"]["maximum_final_objective"]),
        "ipopt_interface_available": capability["usable_ipopt_interface"],
        "ipopt_benchmark_executed": bool(ipopt_local["executed"]),
        "ipopt_refinement_reaches_benchmark": bool(
            ipopt_local["success"]
            and ipopt_local["best_objective"] is not None
            and float(ipopt_local["best_objective"])
            <= float(config["qualification"]["maximum_final_objective"])
        ),
        "scipy_not_accepted_as_final_substitute": not bool(
            config["local"]["allow_scipy_as_final_substitute"]
        ),
    }
    passed = all(checks.values())
    gate = {
        "gate": "phase_C_hybrid_engine_qualification",
        "verdict": "PASS" if passed else "FAIL",
        "checks": checks,
        "capability": capability,
        "pso": {
            "best_x": first.x.tolist(),
            "best_objective": float(first.fun),
            "evaluations": int(first.evaluations),
            "repeat_coordinate_difference": repeat_difference,
        },
        "scipy_diagnostic_local": {
            "accepted_as_final": False,
            "success": bool(scipy_local.success),
            "best_x": scipy_local.x.tolist(),
            "best_objective": float(scipy_local.fun),
        },
        "ipopt_local": ipopt_local,
        "profiles_for_physical_execution": False,
        "next_action": (
            "run qualification in the approved Pyomo/IPOPT environment"
            if not capability["usable_ipopt_interface"]
            else "implement and validate the MBDoE objective adapter"
        ),
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "hybrid_engine_qualification", config)
    gate_path = run_dir / "hybrid_engine_gate.json"
    history_path = run_dir / "pso_benchmark_history.csv"
    config_path = run_dir / "hybrid_engine_config.json"
    write_json(gate_path, gate)
    write_json(config_path, config)
    np.savetxt(history_path, first.history, delimiter=",", header="best_objective", comments="")
    manifest = build_manifest(
        run_dir=run_dir,
        stage="hybrid_engine_qualification",
        config=config,
        sources={"engine_config": CONFIG_PATH},
        code_paths=[Path(__file__), ADAPTIVE_DIR / "hybrid_optimizer.py"],
        random_seeds=[int(config["seed"])],
        status="completed" if passed else "blocked_missing_required_solver",
        convergence={
            "pso": gate["pso"],
            "scipy_diagnostic_local": gate["scipy_diagnostic_local"],
            "ipopt_local": gate["ipopt_local"],
        },
        gate=gate,
        outputs=[gate_path, history_path, config_path],
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
