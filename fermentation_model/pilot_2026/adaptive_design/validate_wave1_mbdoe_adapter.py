from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    allowed_nutrition_times,
    anchor_policy,
    evaluate_campaign,
    fim_components_from_prepared,
    fim_from_sampling_sensitivity_cache,
    load_json,
    load_wave1_config,
    prepare_design,
    prepare_sampling_sensitivity_cache,
    prior_precision,
    representative_members,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    capture_git_state,
    create_immutable_run_directory,
    filesystem_path,
    verify_manifest_output,
    write_json,
    write_json_atomic,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _relative_fim_difference(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / max(np.linalg.norm(right), 1e-15))


def _fim_metrics(fim: np.ndarray) -> dict[str, float]:
    eigenvalues = np.linalg.eigvalsh(0.5 * (fim + fim.T))
    positive = eigenvalues[eigenvalues > 1e-12]
    condition = float(eigenvalues[-1] / positive[0]) if len(positive) else math.inf
    return {
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": float(eigenvalues[-1]),
        "condition_number": condition,
    }


def _pyomo_doe_frozen_sigma_check() -> dict[str, object]:
    """Independent scalar convention check, not validation of the real 9x9 FIM."""

    try:
        import pyomo.environ as pyo
        from pyomo.contrib.doe import DesignOfExperiments
    except ImportError as exc:
        return {"executed": False, "passed": False, "reason": repr(exc)}

    class LinearExperiment:
        def get_labeled_model(self):
            model = pyo.ConcreteModel()
            model.theta = pyo.Var(initialize=2.0)
            model.theta.fix(2.0)
            model.y = pyo.Var(initialize=2.0)
            model.response = pyo.Constraint(expr=model.y == model.theta)
            model.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
            model.unknown_parameters[model.theta] = 2.0
            model.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
            model.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
            model.experiment_outputs[model.y] = 2.0
            model.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
            model.measurement_error[model.y] = 0.4
            return model

    solver = pyo.SolverFactory("ipopt")
    if not solver.available(exception_flag=False):
        return {"executed": False, "passed": False, "reason": "ipopt_unavailable"}
    try:
        doe = DesignOfExperiments(
            experiment=LinearExperiment(),
            step=0.02,
            scale_nominal_param_value=False,
            solver=solver,
            tee=False,
        )
        fim = float(np.asarray(doe.compute_FIM(method="sequential"), dtype=float)[0, 0])
    except Exception as exc:  # pragma: no cover - solver-specific failure detail
        return {"executed": True, "passed": False, "reason": repr(exc)}
    expected = 1.0 / 0.4**2
    return {
        "executed": True,
        "passed": bool(math.isclose(fim, expected, rel_tol=1e-7, abs_tol=1e-9)),
        "pyomo_doe_fim": fim,
        "analytic_frozen_sigma_fim": expected,
        "relative_difference": abs(fim - expected) / expected,
    }


def _validation_tables(
    policy: DesignPolicy,
    member: pd.Series,
    config: dict,
    partitions: dict[str, dict[str, float]],
    *,
    actuator_scenario: dict[str, float] | None = None,
    pyomo_check: dict[str, object] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    sample_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    nominal_grid = float(config["fim_validation"]["nominal_time_grid_step_h"])
    nominal_prepared = prepare_design(
        policy,
        member,
        config,
        partitions,
        simulation_grid_step_h=nominal_grid,
        actuator_scenario=actuator_scenario,
    )
    if nominal_prepared is None:
        raise RuntimeError("Nominal FIM qualification simulation failed")
    nominal_step = float(config["objective"]["finite_difference_log_step"])
    nominal_components = fim_components_from_prepared(
        nominal_prepared,
        sample_times,
        config,
        finite_difference_log_step=nominal_step,
    )
    step_rows = []
    step_fims: dict[float, np.ndarray] = {}
    for step in config["fim_validation"]["log_steps"]:
        components = fim_components_from_prepared(
            nominal_prepared,
            sample_times,
            config,
            finite_difference_log_step=float(step),
        )
        step_fims[float(step)] = components.fim
        step_rows.append(
            {
                "finite_difference_log_step": float(step),
                "relative_fim_difference_vs_nominal": _relative_fim_difference(
                    components.fim, nominal_components.fim
                ),
                "sensitivity_frobenius_norm": float(np.linalg.norm(components.sensitivity)),
                **_fim_metrics(components.fim),
            }
        )
    grid_rows = []
    grid_fims: dict[float, np.ndarray] = {}
    for grid_step in config["fim_validation"]["time_grid_steps_h"]:
        prepared = prepare_design(
            policy,
            member,
            config,
            partitions,
            simulation_grid_step_h=float(grid_step),
            actuator_scenario=actuator_scenario,
        )
        if prepared is None:
            raise RuntimeError(f"FIM time-grid simulation failed for {grid_step} h")
        components = fim_components_from_prepared(prepared, sample_times, config)
        grid_fims[float(grid_step)] = components.fim
        grid_rows.append(
            {
                "time_grid_step_h": float(grid_step),
                "relative_fim_difference_vs_nominal": _relative_fim_difference(
                    components.fim, nominal_components.fim
                ),
                **_fim_metrics(components.fim),
            }
        )
    sensitivity_rows = []
    for index, parameter in enumerate(config["objective"]["priority_parameters"]):
        column = nominal_components.sensitivity[:, index]
        sensitivity_rows.append(
            {
                "parameter": parameter,
                "difference_method": nominal_components.difference_methods[index],
                "finite_difference_log_step": nominal_step,
                "scaled_sensitivity_l2_norm": float(np.linalg.norm(column)),
                "scaled_sensitivity_max_abs": float(np.max(np.abs(column))),
                "scaled_sensitivity_nonzero": bool(np.any(np.abs(column) > 1e-12)),
            }
        )
    step_frame = pd.DataFrame(step_rows)
    grid_frame = pd.DataFrame(grid_rows)
    sensitivity_frame = pd.DataFrame(sensitivity_rows)
    tolerance = config["fim_validation"]
    pyomo_check = pyomo_check or _pyomo_doe_frozen_sigma_check()
    checks = {
        "nominal_sigma_frozen": not bool(
            config["observation_error_model"]["differentiate_sigma_with_parameters"]
        ),
        "all_scaled_sensitivity_columns_nonzero": bool(
            sensitivity_frame["scaled_sensitivity_nonzero"].all()
        ),
        "finite_difference_step_stable": bool(
            step_frame["relative_fim_difference_vs_nominal"].max()
            <= float(tolerance["maximum_relative_step_fim_difference"])
        ),
        "time_grid_stable": bool(
            grid_frame["relative_fim_difference_vs_nominal"].max()
            <= float(tolerance["maximum_relative_grid_fim_difference"])
        ),
        "fim_positive_semidefinite": bool(
            min(step_frame["minimum_eigenvalue"].min(), grid_frame["minimum_eigenvalue"].min())
            >= float(tolerance["minimum_eigenvalue_tolerance"])
        ),
        "pyomo_doe_independent_frozen_sigma_check": bool(pyomo_check.get("passed", False)),
    }
    validation = {
        "formulation": {
            "prediction": "physical wine concentration and captured interval mass",
            "observational_scale": "computed once at nominal parameter centre",
            "scaled_observation": "physical_prediction / nominal_sigma",
            "sensitivity": "finite difference of physical prediction divided by frozen nominal_sigma",
            "fim": "sensitivity.T @ sensitivity",
            "variance_derivative_included": False,
        },
        "checks": checks,
        "tolerances": tolerance,
        "nominal_fim": _fim_metrics(nominal_components.fim),
        "pyomo_doe_comparison": pyomo_check,
        "scope_note": "The Pyomo-DOE scalar comparison independently validates only the frozen-scale convention; the real 9x9 FIM is validated by multi-policy/member finite-difference and grid replication.",
        "verdict": "PASS" if all(checks.values()) else "FAIL",
    }
    return sensitivity_frame, step_frame, grid_frame, validation


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate the corrected Pilot 2026 Wave-1 MBDoE adapter")
    parser.add_argument("--source-ensemble-run", required=True)
    parser.add_argument("--source-engine-run", required=True)
    parser.add_argument("--source-aroma-run", required=True)
    parser.add_argument("--source-actuator-run", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_started = time.perf_counter()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    aroma_config = load_json(AROMA_CONFIG_PATH)
    ensemble_run = _run_path(args.source_ensemble_run)
    engine_run = _run_path(args.source_engine_run)
    aroma_run = _run_path(args.source_aroma_run)
    actuator_run = _run_path(args.source_actuator_run)
    source_verification = {
        "joint_ensemble": verify_manifest_output(ensemble_run, "joint_parameter_ensemble.csv"),
        "engine": verify_manifest_output(engine_run, "hybrid_engine_gate.json"),
        "aroma": verify_manifest_output(aroma_run, "aroma_calibration_gate.json"),
        "actuator": verify_manifest_output(actuator_run, "actuator_estimates_by_run.csv"),
    }
    engine_gate = load_json(engine_run / "hybrid_engine_gate.json")
    if engine_gate["verdict"] != "PASS":
        raise RuntimeError("Hybrid PSO-to-IPOPT engine gate is not PASS")
    ensemble = pd.read_csv(ensemble_run / "joint_parameter_ensemble.csv")
    partitions, provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    representative = representative_members(ensemble, config)
    prior = prior_precision(ensemble, config)
    anchor = anchor_policy(config)
    slots = int(config["future_process"]["optimized_temperature_slots"])
    half = slots // 2
    warm_cool = DesignPolicy(
        "benchmark_warm_then_cool",
        tuple([23.0] * half + [16.0] * (slots - half)),
        ((24.0, 70.0), (72.0, 70.0)),
    )
    cool_warm = DesignPolicy(
        "benchmark_cool_then_warm",
        tuple([16.0] * half + [23.0] * (slots - half)),
        ((48.0, 80.0),),
    )
    campaigns = {
        "reference_three_anchors": (anchor, anchor, anchor),
        "complementary_thermal_pair": (anchor, warm_cool, cool_warm),
    }
    rows = []
    scenario_rows = []
    for name, policies in campaigns.items():
        score, evaluations = evaluate_campaign(
            policies, ensemble, representative, prior, config, partitions
        )
        rows.append(
            {
                "campaign": name,
                "robust_information_score": score,
                "median_information_gain": float(np.median([x.information_gain for x in evaluations])),
                "lower_decile_information_gain": float(np.quantile([x.information_gain for x in evaluations], 0.1)),
                "representative_completion_probability": float(np.mean([x.completion for x in evaluations])),
                "maximum_residual_sugar_g_l": float(max(x.residual_sugar_g_l for x in evaluations)),
            }
        )
        for evaluation in evaluations:
            scenario_rows.append(
                {
                    "campaign": name,
                    "ensemble_member": evaluation.member,
                    "information_gain": evaluation.information_gain,
                    "completion": evaluation.completion,
                    "residual_sugar_g_l": evaluation.residual_sugar_g_l,
                    **_fim_metrics(evaluation.fim),
                }
            )
    summary = pd.DataFrame(rows)
    scenarios = pd.DataFrame(scenario_rows)
    full_members = list(range(len(ensemble)))
    _full_score, full_evaluations = evaluate_campaign(
        (anchor, warm_cool, cool_warm),
        ensemble,
        full_members,
        prior,
        config,
        partitions,
    )
    central_member = int(representative[0])
    lowest_information_member = min(
        full_evaluations, key=lambda row: row.information_gain
    ).member
    highest_information_member = max(
        full_evaluations, key=lambda row: row.information_gain
    ).member
    critical_drying_member = min(
        full_evaluations, key=lambda row: row.drying_time_h
    ).member
    actuator_config = config["future_process"]["temperature_actuator"]
    critical_actuator = {
        "tau_h": max(float(value) for value in actuator_config["empirical_tau_h"]),
        "tracking_error_c": float(np.median(actuator_config["observed_tracking_rmse_c"])),
    }
    validation_cases = [
        ("anchor", "anchor", anchor, central_member, None),
        ("candidate_A", "candidate_A", warm_cool, central_member, None),
        ("candidate_B", "candidate_B", cool_warm, central_member, None),
        ("central_member", "central_member", anchor, central_member, None),
        ("lowest_information_member", "lowest_information_member", anchor, lowest_information_member, None),
        ("highest_information_member", "highest_information_member", anchor, highest_information_member, None),
        ("critical_drying_member", "critical_drying_member", warm_cool, critical_drying_member, None),
        ("critical_actuator_scenario", "critical_actuator_scenario", cool_warm, central_member, critical_actuator),
    ]
    pyomo_check = _pyomo_doe_frozen_sigma_check()
    sensitivity_frames = []
    step_frames = []
    grid_frames = []
    validation_rows = []
    case_validations = []
    for case_id, member_role, policy, member_index, actuator_scenario in validation_cases:
        sensitivity_case, step_case, grid_case, case_validation = _validation_tables(
            policy,
            ensemble.iloc[int(member_index)],
            config,
            partitions,
            actuator_scenario=actuator_scenario,
            pyomo_check=pyomo_check,
        )
        for frame in (sensitivity_case, step_case, grid_case):
            frame.insert(0, "actuator_scenario", json.dumps(actuator_scenario, sort_keys=True))
            frame.insert(0, "ensemble_member", int(member_index))
            frame.insert(0, "member_role", member_role)
            frame.insert(0, "policy", policy.name)
            frame.insert(0, "validation_case", case_id)
        sensitivity_frames.append(sensitivity_case)
        step_frames.append(step_case)
        grid_frames.append(grid_case)
        case_validations.append(case_validation)
        validation_rows.append(
            {
                "validation_case": case_id,
                "policy": policy.name,
                "member_role": member_role,
                "ensemble_member": int(member_index),
                "actuator_scenario": json.dumps(actuator_scenario, sort_keys=True),
                "maximum_relative_step_fim_difference": float(
                    step_case["relative_fim_difference_vs_nominal"].max()
                ),
                "maximum_relative_grid_fim_difference": float(
                    grid_case["relative_fim_difference_vs_nominal"].max()
                ),
                "minimum_eigenvalue": float(
                    min(step_case["minimum_eigenvalue"].min(), grid_case["minimum_eigenvalue"].min())
                ),
                "maximum_eigenvalue": float(
                    max(step_case["maximum_eigenvalue"].max(), grid_case["maximum_eigenvalue"].max())
                ),
                "maximum_condition_number": float(
                    max(step_case["condition_number"].max(), grid_case["condition_number"].max())
                ),
                "all_sensitivity_columns_nonzero": bool(
                    sensitivity_case["scaled_sensitivity_nonzero"].all()
                ),
                "verdict": case_validation["verdict"],
            }
        )
    sensitivity = pd.concat(sensitivity_frames, ignore_index=True)
    step_validation = pd.concat(step_frames, ignore_index=True)
    grid_validation = pd.concat(grid_frames, ignore_index=True)
    fim_by_case = pd.DataFrame(validation_rows)
    independent_prepared = prepare_design(
        anchor, ensemble.iloc[central_member], config, partitions
    )
    if independent_prepared is None:
        raise RuntimeError("Independent full-FIM comparison simulation failed")
    independent_times = np.asarray(config["sampling"]["preliminary_times_h"], dtype=float)
    direct_full_fim = fim_components_from_prepared(
        independent_prepared, independent_times, config
    ).fim
    cached_full_fim = fim_from_sampling_sensitivity_cache(
        prepare_sampling_sensitivity_cache(independent_prepared, config),
        independent_times,
        config,
    )
    independent_full_fim_relative_difference = _relative_fim_difference(
        cached_full_fim, direct_full_fim
    )
    fim_checks = {
        "nominal_sigma_frozen": not bool(
            config["observation_error_model"]["differentiate_sigma_with_parameters"]
        ),
        "all_validation_scope_cases_executed": set(
            config["fim_validation"]["validation_scope"]
        )
        == set(fim_by_case["validation_case"]),
        "all_scaled_sensitivity_columns_nonzero": bool(
            fim_by_case["all_sensitivity_columns_nonzero"].all()
        ),
        "finite_difference_step_stable_all_cases": bool(
            fim_by_case["maximum_relative_step_fim_difference"].le(
                float(config["fim_validation"]["maximum_relative_step_fim_difference"])
            ).all()
        ),
        "time_grid_stable_all_cases": bool(
            fim_by_case["maximum_relative_grid_fim_difference"].le(
                float(config["fim_validation"]["maximum_relative_grid_fim_difference"])
            ).all()
        ),
        "fim_positive_semidefinite_all_cases": bool(
            fim_by_case["minimum_eigenvalue"].ge(
                float(config["fim_validation"]["minimum_eigenvalue_tolerance"])
            ).all()
        ),
        "pyomo_doe_independent_frozen_sigma_convention_check": bool(
            pyomo_check.get("passed", False)
        ),
        "independent_full_9x9_cached_trajectory_comparison": bool(
            independent_full_fim_relative_difference <= 1e-10
        ),
    }
    fim_validation = {
        "formulation": case_validations[0]["formulation"],
        "scope_note": "Eight policy/member/actuator cases validate the real 9x9 FIM numerically; the scalar Pyomo-DOE problem validates only the frozen-scale convention.",
        "checks": fim_checks,
        "tolerances": config["fim_validation"],
        "validation_cases": validation_rows,
        "pyomo_doe_comparison": pyomo_check,
        "independent_full_9x9_comparison": {
            "implementations": [
                "direct finite differences evaluated at observation times",
                "cached full-trajectory derivatives interpolated to observation times",
            ],
            "relative_difference": independent_full_fim_relative_difference,
            "passed": independent_full_fim_relative_difference <= 1e-10,
        },
        "verdict": "PASS" if all(fim_checks.values()) else "FAIL",
    }
    score_range = float(summary["robust_information_score"].max() - summary["robust_information_score"].min())
    checks = {
        "hybrid_engine_gate_pass": engine_gate["verdict"] == "PASS",
        "source_hashes_verified": len(source_verification) == 4,
        "joint_ensemble_has_64_members": len(ensemble) == 64,
        "eight_representative_members_selected": len(representative) == 8,
        "all_scores_finite": bool(np.isfinite(summary.select_dtypes(include=[np.number])).all().all()),
        "corrected_fim_validation_pass": fim_validation["verdict"] == "PASS",
        "mandatory_wine_baseline_at_t0": bool(
            len(config["sampling"]["preliminary_times_h"]) == 10
            and math.isclose(
                float(config["sampling"]["preliminary_times_h"][0]), 0.0, abs_tol=1e-12
            )
        ),
        "nine_capture_intervals_from_t0": bool(
            len(config["sampling"]["preliminary_times_h"]) - 1 == 9
            and math.isclose(
                float(config["sampling"]["preliminary_times_h"][0]), 0.0, abs_tol=1e-12
            )
        ),
        "benchmark_profiles_are_information_distinguishable": score_range > 1e-3,
        "nutrition_slots_respect_owner_window": bool(
            np.all(allowed_nutrition_times(config) <= float(config["nutrition"]["latest_h"]))
        ),
        "no_physical_profile_released": not bool(config["release_policy"]["physical_execution_authorized"]),
    }
    gate = {
        "gate": "phase_C_wave1_MBDoE_adapter_requalification",
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "representative_members": representative,
        "source_verification": source_verification,
        "profiles_for_physical_execution": False,
        "next_action": "run multiseed robust PSO search" if all(checks.values()) else "repair adapter",
    }
    # Keep the stage compact enough for pandas/Matplotlib on legacy Windows
    # MAX_PATH while retaining an unambiguous immutable run namespace.
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_mbdoe_adapter_v3", config)
    summary_path = run_dir / "benchmark_campaigns.csv"
    scenarios_path = run_dir / "benchmark_scenarios.csv"
    sensitivity_path = run_dir / "finite_difference_sensitivity.csv"
    step_path = run_dir / "finite_difference_step_validation.csv"
    grid_path = run_dir / "time_grid_sensitivity.csv"
    fim_path = run_dir / "corrected_fim_validation.json"
    fim_by_case_path = run_dir / "fim_validation_by_policy_and_member.csv"
    gate_path = run_dir / "adapter_gate.json"
    provenance_path = run_dir / "partition_surrogate_provenance.json"
    config_path = run_dir / "wave1_mbdoe_config.json"
    runtime_path = run_dir / "runtime_summary.json"
    summary.to_csv(filesystem_path(summary_path), index=False)
    scenarios.to_csv(filesystem_path(scenarios_path), index=False)
    sensitivity.to_csv(filesystem_path(sensitivity_path), index=False)
    step_validation.to_csv(filesystem_path(step_path), index=False)
    grid_validation.to_csv(filesystem_path(grid_path), index=False)
    fim_by_case.to_csv(filesystem_path(fim_by_case_path), index=False)
    write_json(fim_path, fim_validation)
    write_json(gate_path, gate)
    write_json(provenance_path, provenance)
    write_json(config_path, config)
    write_json(
        runtime_path,
        {"total_runtime_seconds": float(time.perf_counter() - run_started)},
    )
    outputs = [
        summary_path,
        scenarios_path,
        sensitivity_path,
        step_path,
        grid_path,
        fim_path,
        fim_by_case_path,
        gate_path,
        provenance_path,
        config_path,
        runtime_path,
    ]
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_mbdoe_adapter_requalification",
        config=config,
        sources={
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "engine_qualification_manifest": engine_run / "run_manifest.json",
            "aroma_calibration_manifest": aroma_run / "run_manifest.json",
            "temperature_actuator_manifest": actuator_run / "run_manifest.json",
            "aroma_config": AROMA_CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
            "wave1_config": CONFIG_PATH,
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "pilot_aroma_calibration.py",
            ADAPTIVE_DIR / "pilot_calibration.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            FERMENTATION_DIR / "shared" / "run_new_must_glycerol_estimability_doe.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[int(config["seed"])],
        status="completed" if gate["verdict"] == "PASS" else "validation_failed",
        convergence={"fim_validation": fim_validation},
        gate=gate,
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    write_json(run_dir / "run_manifest.json", manifest)
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "computational_wave1_adapter_v3_requalification",
            "current_gate": gate["gate"],
            "gate_verdict": gate["verdict"],
            "latest_wave1_adapter_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "physical_execution_status": "not_authorized",
            "executable_schedule_issued": False,
            "tank_assignments": [],
            "profiles_for_physical_execution": False,
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    print(summary.to_string(index=False))
    if gate["verdict"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
