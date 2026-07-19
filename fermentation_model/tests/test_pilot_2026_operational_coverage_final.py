from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.hybrid_optimizer import (  # noqa: E402
    initialize_checkpointed_swarm,
)
from pilot_2026.adaptive_design.local_refinement import (  # noqa: E402
    sequential_ipopt_refine,
)
from pilot_2026.adaptive_design.operational_coverage import (  # noqa: E402
    AUTHORIZATION_FLAGS,
    DesignPolicy,
    derive_robust_setpoint_envelope,
    feasible_information_policy_from_unit_vector,
    policy_member_drying_margin_row,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    load_wave1_config,
)
from pilot_2026.adaptive_design.run_wave1_operational_coverage_final import (  # noqa: E402
    DEPENDENCY_PATHS,
    _balanced_practical_convergence,
    _dependency_manifest,
    _latest_valid_checkpoint,
    _model_config,
    _save_checkpoint,
)


class FinalOperationalCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.final = json.loads(
            (ADAPTIVE_DIR / "coverage_final_design_config.json").read_text(
                encoding="utf-8"
            )
        )
        cls.coverage = json.loads(
            (ADAPTIVE_DIR / "coverage_design_config.json").read_text(
                encoding="utf-8"
            )
        )
        cls.model = _model_config(cls.final)
        cls.envelope = derive_robust_setpoint_envelope(cls.model, cls.coverage)

    def test_owner_decisions_and_authorization_are_exact(self) -> None:
        self.assertEqual(self.final["owner_selection"]["strategy"], "II_diagonal_1")
        self.assertEqual(self.final["owner_selection"]["dose_design"], "N76_all_conservative")
        self.assertEqual(self.final["authorization"], AUTHORIZATION_FLAGS)
        self.assertEqual(self.final["required_minimum_action_to_drying_margin_h"], 24.0)
        self.assertEqual(self.model["anchor"]["nutrition_pulses_mg_yan_l"], [[46.0, 76.0]])

    def test_drying_margin_never_crosses_policies(self) -> None:
        policy_a = DesignPolicy(
            "tank_a",
            (18.0, 18.0, 22.0, *([22.0] * 11)),
            ((2.0, 76.0),),
        )
        policy_b = DesignPolicy(
            "tank_b",
            tuple([18.0] * 14),
            ((46.0, 76.0),),
        )
        row_a = policy_member_drying_margin_row(policy_a, 0, 100.0, self.model)
        row_b = policy_member_drying_margin_row(policy_b, 0, 80.0, self.model)
        self.assertEqual(row_a["latest_active_action_h"], 24.0)
        self.assertEqual(row_a["action_to_drying_margin_h"], 76.0)
        self.assertEqual(row_b["latest_active_action_h"], 46.0)
        self.assertEqual(row_b["action_to_drying_margin_h"], 34.0)
        self.assertEqual(
            min(row_a["action_to_drying_margin_h"], row_b["action_to_drying_margin_h"]),
            34.0,
        )
        self.assertNotEqual(80.0 - 46.0, 100.0 - 46.0)

    def test_information_decoder_retains_physical_constraints_without_archetypes(self) -> None:
        vector = np.linspace(0.05, 0.95, 12)
        policy = feasible_information_policy_from_unit_vector(
            "information",
            76.0,
            vector,
            self.model,
            self.coverage,
            self.envelope,
        )
        jumps = np.abs(np.diff(np.asarray(policy.temperature_c)))
        active = jumps[jumps > 1e-12]
        self.assertLessEqual(max(active, default=0.0), 5.0 + 1e-9)
        self.assertTrue(np.all(active >= 1.0 - 1e-9))
        self.assertLessEqual(max(abs(policy.temperature_c[0] - t) for t in (17, 18, 19)), 5.0)
        self.assertEqual(policy.nutrition_mg_yan_l[0][1], 76.0)

    def test_practical_convergence_requires_explicit_point_one_percent_threshold(self) -> None:
        policies = (
            DesignPolicy("anchor", tuple([18.0] * 14), ((46.0, 76.0),)),
            DesignPolicy("cold", (17.0, 18.0, *([18.0] * 12)), ((2.0, 76.0),)),
            DesignPolicy("warm", (21.0, 24.0, *([24.0] * 12)), ((46.0, 76.0),)),
        )
        champions = [
            {
                "independent_seed": seed,
                "robust_score": 20.0 - index * 0.01,
                "feasible": True,
                "policies": policies,
            }
            for index, seed in enumerate(self.final["balanced_search"]["seeds"])
        ]
        seed_results = [
            {
                "continuation_windows": [
                    {"seed": seed, "continuation_improvement_fraction": 0.0011}
                ]
            }
            for seed in self.final["balanced_search"]["seeds"]
        ]
        failed = _balanced_practical_convergence(
            champions, seed_results, self.model, self.final
        )
        self.assertFalse(
            failed["checks"]["continuation_improvement_at_most_0_1pct"]
        )
        for result in seed_results:
            result["continuation_windows"][-1][
                "continuation_improvement_fraction"
            ] = 0.001
        passed = _balanced_practical_convergence(
            champions, seed_results, self.model, self.final
        )
        self.assertTrue(
            passed["checks"]["continuation_improvement_at_most_0_1pct"]
        )

    def test_ipopt_trust_region_executes_and_records_kkt(self) -> None:
        result = sequential_ipopt_refine(
            np.asarray([0.85]),
            np.asarray([[0.0, 1.0]]),
            lambda value: float((value[0] - 0.2) ** 2),
            np.asarray([0]),
            self.model,
            candidate_id="synthetic",
            validation_objective=lambda value: float((value[0] - 0.2) ** 2),
        )
        self.assertNotEqual(result.state, "not_executed")
        self.assertTrue(result.qualified)
        self.assertTrue(result.trace)
        ipopt_rows = [row for row in result.trace if "ipopt_termination" in row]
        self.assertTrue(ipopt_rows)
        solved = [row for row in ipopt_rows if row["ipopt_termination"] != "not_needed_stationary"]
        self.assertTrue(solved)
        self.assertTrue(all(math.isfinite(float(row["kkt_error"])) for row in solved))

    def test_checkpoint_roundtrip_preserves_full_swarm_and_hash_contract(self) -> None:
        state = initialize_checkpointed_swarm(
            lambda value: float(np.sum(value**2)),
            np.asarray([[0.0, 1.0], [0.0, 1.0]]),
            particles=3,
            seed=123,
            fidelity="four_member",
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = _save_checkpoint(
                checkpoint_dir=directory,
                kind="balanced",
                state=state,
                sequence=1,
                phase="complete",
                dependency_hash="dep",
                config_hash="cfg",
                objective_cache={"four_member|x": 1.0},
                base_vector=np.asarray([0.2, 0.3]),
                continuation_windows=[
                    {"continuation_improvement_fraction": 0.0005}
                ],
                complete=True,
            )
            latest = _latest_valid_checkpoint(directory, "balanced", 123, "dep", "cfg")
            self.assertIsNotNone(latest)
            self.assertEqual(latest[0].name, path.name)
            payload = latest[1]
            self.assertEqual(payload["swarm"]["rng_state"], state.rng_state)
            self.assertEqual(payload["objective_cache"]["four_member|x"], 1.0)
            self.assertTrue(payload["complete"])
            self.assertIn("candidate_by_seed", payload)

    def test_dependency_manifest_hashes_all_declared_transitive_files(self) -> None:
        manifest = _dependency_manifest()
        self.assertEqual(manifest["transitive_dependency_count"], len(DEPENDENCY_PATHS))
        names = {Path(row["path"]).name for row in manifest["dependencies"]}
        for required in (
            "operational_coverage.py",
            "run_wave1_operational_coverage.py",
            "pilot_mbdoe_adapter.py",
            "hybrid_optimizer.py",
            "local_refinement.py",
            "final_search_logic.py",
            "pilot_aroma_calibration.py",
            "pilot_calibration.py",
            "optimize_wave1_sampling_and_plots.py",
            "run_artifacts.py",
        ):
            self.assertIn(required, names)


if __name__ == "__main__":
    unittest.main()
