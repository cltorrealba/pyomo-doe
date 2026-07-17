from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.resolve_historical_audit import (  # noqa: E402
    build_resolution,
)
from pilot_2026.adaptive_design.run_hierarchical_calibration import (  # noqa: E402
    evaluate_gate,
)


def nested_value(payload: dict, dotted_path: str):
    value = payload
    for key in dotted_path.split("."):
        value = value[key]
    return value


class Pilot2026CalibrationGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(
            (ADAPTIVE_DIR / "calibration_config.json").read_text(encoding="utf-8")
        )

    def test_entry_gate_fails_on_current_data_or_solver_failure(self) -> None:
        gate = evaluate_gate(
            self.config,
            {"status": "completed"},
            {"verdict": "PASS"},
            {},
            {"pass": False},
            {"scipy_least_squares_available": False},
        )
        self.assertEqual(gate["verdict"], "FAIL")
        self.assertFalse(gate["fit_executed"])
        self.assertFalse(gate["profiles_for_physical_execution"])
        self.assertIn("repository_raw_hash_audit_pass", gate["blockers"])
        self.assertIn("scipy_least_squares_available", gate["blockers"])

    def test_entry_gate_passes_without_claiming_physical_release(self) -> None:
        gate = evaluate_gate(
            self.config,
            {"status": "completed"},
            {"verdict": "PASS"},
            {},
            {"pass": True},
            {"scipy_least_squares_available": True},
        )
        self.assertEqual(gate["verdict"], "PASS")
        self.assertFalse(gate["blockers"])
        self.assertFalse(gate["fit_executed"])
        self.assertFalse(gate["profiles_for_physical_execution"])

    def test_historical_53_errors_are_resolved_without_rebaseline(self) -> None:
        resolution = build_resolution()
        self.assertEqual(resolution["verdict"], "PASS")
        self.assertEqual(resolution["historical_error_count"], 53)
        self.assertEqual(resolution["classification"]["raw_worktree_size_drift"], 52)
        self.assertEqual(resolution["classification"]["ambiguous_generated_result_root"], 1)
        self.assertTrue(resolution["checks"]["raw_tree_unchanged_since_base_commit"])
        self.assertTrue(resolution["checks"]["raw_manifest_byte_identical_to_base"])

    def test_yan_oculyze_error_and_censoring_contracts_are_explicit(self) -> None:
        self.assertEqual(self.config["yan"]["mass_derived_per_pulse_mg_l"], 90.434783)
        self.assertEqual(self.config["yan"]["protocol_nominal_per_pulse_mg_l"], 80.0)
        self.assertEqual(len(self.config["yan"]["delivery_fraction_bounds"]), 2)
        self.assertIn("scale_bounds", self.config["oculyze"])
        self.assertIn("absolute_floor", self.config["observation_error"])
        self.assertIn("likelihood", self.config["censoring"])
        self.assertEqual(
            self.config["co2_fit"]["ess_method"],
            "per_run_AR1_from_final_standardized_residuals",
        )

    def test_existing_pilot_2025_integrated_result_is_never_an_official_prior(self) -> None:
        policy = self.config["pilot_2025_policy"]
        self.assertFalse(policy["use_existing_integrated_result_as_prior"])
        self.assertIn("regenerate", policy["required_action_before_use"])

    def test_design_constraints_list_every_missing_safety_bound(self) -> None:
        constraints = json.loads(
            (ADAPTIVE_DIR / "design_constraints.json").read_text(encoding="utf-8")
        )
        self.assertEqual(constraints["review_status"], "pending_owner_approval")
        self.assertFalse(constraints["approval"]["physical_execution_authorized"])
        for dotted_path in constraints["fail_closed_fields"]:
            self.assertIsNone(nested_value(constraints, dotted_path), dotted_path)


if __name__ == "__main__":
    unittest.main()
