from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import pandas as pd


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_DIR = FERMENTATION_DIR.parent
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))


class Pilot2026JointEnsembleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        state = json.loads((ADAPTIVE_DIR / "campaign_state.json").read_text(encoding="utf-8"))
        cls.run_dir = REPOSITORY_DIR / state["latest_joint_ensemble_run"]
        cls.gate = json.loads(
            (cls.run_dir / "ensemble_gate.json").read_text(encoding="utf-8")
        )
        cls.diagnostics = json.loads(
            (cls.run_dir / "ensemble_diagnostics.json").read_text(encoding="utf-8")
        )
        cls.ensemble = pd.read_csv(cls.run_dir / "joint_parameter_ensemble.csv")

    def test_joint_ensemble_gate_passes_without_releasing_profiles(self) -> None:
        self.assertEqual(self.gate["verdict"], "PASS")
        self.assertFalse(self.gate["profiles_for_physical_execution"])
        self.assertTrue(all(self.gate["checks"].values()))

    def test_ensemble_has_64_members_and_multiple_primary_centres(self) -> None:
        self.assertEqual(len(self.ensemble), 64)
        self.assertEqual(self.ensemble["ensemble_member"].nunique(), 64)
        self.assertGreaterEqual(self.ensemble["primary_centre_index"].nunique(), 3)

    def test_all_declared_weak_aroma_directions_remain_broad(self) -> None:
        for details in self.diagnostics["aroma"].values():
            for parameter in details["weak_parameters"]:
                self.assertGreaterEqual(details["sample_log_sd"][parameter], 1.5)


if __name__ == "__main__":
    unittest.main()
