from __future__ import annotations

import json
import sys
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
    load_json,
    prior_precision,
    representative_members,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def main() -> None:
    config = load_json(CONFIG_PATH)
    aroma_config = load_json(AROMA_CONFIG_PATH)
    ensemble_run = REPOSITORY_DIR / config["source_contract"]["joint_ensemble_run"]
    engine_run = REPOSITORY_DIR / config["source_contract"]["engine_qualification_run"]
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
            eigenvalues = np.linalg.eigvalsh(evaluation.fim)
            scenario_rows.append(
                {
                    "campaign": name,
                    "ensemble_member": evaluation.member,
                    "information_gain": evaluation.information_gain,
                    "completion": evaluation.completion,
                    "residual_sugar_g_l": evaluation.residual_sugar_g_l,
                    "minimum_fim_eigenvalue": float(eigenvalues.min()),
                    "maximum_fim_eigenvalue": float(eigenvalues.max()),
                }
            )
    summary = pd.DataFrame(rows)
    scenarios = pd.DataFrame(scenario_rows)
    score_range = float(summary["robust_information_score"].max() - summary["robust_information_score"].min())
    checks = {
        "hybrid_engine_gate_pass": engine_gate["verdict"] == "PASS",
        "joint_ensemble_has_64_members": len(ensemble) == 64,
        "eight_representative_members_selected": len(representative) == 8,
        "all_scores_finite": bool(np.isfinite(summary.select_dtypes(include=[np.number])).all().all()),
        "fim_positive_semidefinite_with_tolerance": bool(
            (scenarios["minimum_fim_eigenvalue"] >= -1e-7).all()
        ),
        "benchmark_profiles_are_information_distinguishable": score_range > 1e-3,
        "nutrition_slots_respect_owner_window": bool(
            np.all(allowed_nutrition_times(config) <= float(config["nutrition"]["latest_h"]))
        ),
        "no_physical_profile_released": not bool(
            config["release_policy"]["physical_execution_authorized"]
        ),
    }
    gate = {
        "gate": "phase_C_wave1_MBDoE_adapter_validation",
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "representative_members": representative,
        "profiles_for_physical_execution": False,
        "next_action": "run robust PSO pair search" if all(checks.values()) else "repair adapter",
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_mbdoe_adapter", config)
    summary_path = run_dir / "benchmark_campaigns.csv"
    scenarios_path = run_dir / "benchmark_scenarios.csv"
    gate_path = run_dir / "adapter_gate.json"
    provenance_path = run_dir / "partition_surrogate_provenance.json"
    config_path = run_dir / "wave1_mbdoe_config.json"
    summary.to_csv(summary_path, index=False)
    scenarios.to_csv(scenarios_path, index=False)
    write_json(gate_path, gate)
    write_json(provenance_path, provenance)
    write_json(config_path, config)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_mbdoe_adapter_validation",
        config=config,
        sources={
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "engine_qualification_manifest": engine_run / "run_manifest.json",
            "aroma_config": AROMA_CONFIG_PATH,
            "wave1_config": CONFIG_PATH,
        },
        code_paths=[Path(__file__), ADAPTIVE_DIR / "pilot_mbdoe_adapter.py"],
        random_seeds=[int(config["seed"])],
        status="completed" if gate["verdict"] == "PASS" else "validation_failed",
        convergence={},
        gate=gate,
        outputs=[summary_path, scenarios_path, gate_path, provenance_path, config_path],
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    print(summary.to_string(index=False))
    if gate["verdict"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
