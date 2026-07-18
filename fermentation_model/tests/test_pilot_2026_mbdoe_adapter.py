from __future__ import annotations

import json
import sys
import unittest
from datetime import datetime
from pathlib import Path

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_DIR = FERMENTATION_DIR.parent
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    _future_design,
    allowed_nutrition_times,
    allowed_sampling_times,
    anchor_policy,
    load_wave1_config,
    policy_from_vector,
    vector_bounds,
)


class Pilot2026MBDoEAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_wave1_config(ADAPTIVE_DIR / "wave1_mbdoe_config.json")
        cls.state = json.loads(
            (ADAPTIVE_DIR / "campaign_state.json").read_text(encoding="utf-8")
        )

    def test_owner_windows_are_encoded_for_nutrition_and_sampling(self) -> None:
        for times in (allowed_nutrition_times(self.config), allowed_sampling_times(self.config)):
            for time in times:
                stamp = datetime.fromisoformat(self.config["future_process"]["start_local"])
                stamp = stamp + __import__("datetime").timedelta(hours=float(time))
                self.assertLess(stamp.weekday(), 5)
                self.assertIn(stamp.hour, {9, 13, 17})

    def test_design_vector_decodes_inside_temperature_and_yan_limits(self) -> None:
        bounds = vector_bounds(self.config)
        values = bounds[:, 1] + 100.0
        policy = policy_from_vector("test", values, self.config)
        self.assertTrue(all(15.0 <= value <= 25.0 for value in policy.temperature_c))
        self.assertLessEqual(sum(amount for _, amount in policy.nutrition_mg_yan_l), 232.0)
        self.assertLessEqual(len(policy.nutrition_mg_yan_l), 3)

    def test_temperature_actuator_is_not_an_instantaneous_setpoint(self) -> None:
        policy = anchor_policy(self.config)
        modified = type(policy)(
            "step",
            (25.0,) + policy.temperature_c[1:],
            policy.nutrition_mg_yan_l,
        )
        design = _future_design(modified, self.config)
        self.assertAlmostEqual(float(design.temperature_c[0]), 18.0)
        self.assertGreater(float(design.temperature_c[1]), 18.0)
        self.assertLess(float(design.temperature_c[1]), 25.0)

    def test_latest_adapter_and_sampling_gates_do_not_release_profiles(self) -> None:
        adapter = REPOSITORY_DIR / self.state["latest_wave1_adapter_run"] / "adapter_gate.json"
        sampling = REPOSITORY_DIR / self.state["latest_wave1_sampling_run"] / "sampling_gate.json"
        adapter_gate = json.loads(adapter.read_text(encoding="utf-8"))
        sampling_gate = json.loads(sampling.read_text(encoding="utf-8"))
        self.assertEqual(adapter_gate["verdict"], "PASS")
        self.assertEqual(sampling_gate["verdict"], "FAIL")
        self.assertFalse(adapter_gate["profiles_for_physical_execution"])
        self.assertFalse(sampling_gate["profiles_for_physical_execution"])
        self.assertFalse(self.state["executable_schedule_issued"])
        self.assertFalse(self.state["profiles_for_physical_execution"])
        self.assertEqual(self.state["tank_assignments"], [])


if __name__ == "__main__":
    unittest.main()
