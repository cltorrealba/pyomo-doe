from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

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

    def test_current_owner_prerequisites_fail_closed(self) -> None:
        gate = evaluate_gate(
            self.config,
            {"status": "completed"},
            {"verdict": "PASS"},
            {"calibration_gate": ["owner decision"]},
            {"pass": False, "stdout": {"errors": ["raw drift"], "warnings": []}},
            {"pyomo_doe_importable": True, "ipopt_available": False},
        )
        self.assertEqual(gate["verdict"], "FAIL")
        self.assertFalse(gate["fit_executed"])
        self.assertFalse(gate["profiles_for_physical_execution"])
        self.assertIn("repository_raw_hash_audit_pass", gate["blockers"])
        self.assertIn("ipopt_available", gate["blockers"])
        self.assertTrue(any(name.startswith("owner_") for name in gate["blockers"]))

    def test_gate_can_pass_only_when_every_check_is_true(self) -> None:
        config = copy.deepcopy(self.config)
        config["owner_prerequisites"] = {
            key: True for key in config["owner_prerequisites"]
        }
        gate = evaluate_gate(
            config,
            {"status": "completed"},
            {"verdict": "PASS"},
            {"calibration_gate": []},
            {"pass": True, "stdout": {"errors": [], "warnings": []}},
            {"pyomo_doe_importable": True, "ipopt_available": True},
        )
        self.assertEqual(gate["verdict"], "PASS")
        self.assertFalse(gate["blockers"])
        self.assertFalse(gate["fit_executed"])

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
