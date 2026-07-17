from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_DIR = FERMENTATION_DIR.parent
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    evaluate_gate,
    load_partition_surrogates,
)


def _fit(species: str, weak: list[str], loss_identified: bool = True):
    identified = {
        "formation_growth_ug_per_g_sugar": True,
        "formation_stationary_ug_per_g_sugar": True,
        "effective_loss_scale": loss_identified,
    }
    for parameter in weak:
        identified[parameter] = False
    return SimpleNamespace(
        species=species,
        validation={
            "converged_multistarts": 12,
            "covariance_finite": True,
            "active_bound_fraction": 0.0,
            "objective_per_observation": 1.0,
            "loss_separately_identified": loss_identified,
            "parameter_identified": identified,
            "weak_parameters": weak,
        },
    )


class Pilot2026AromaCalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(
            (ADAPTIVE_DIR / "aroma_calibration_config.json").read_text(encoding="utf-8")
        )

    def test_authoritative_unifac_surrogate_loads_without_fallback(self) -> None:
        models, provenance = load_partition_surrogates(self.config, REPOSITORY_DIR)
        self.assertEqual(set(models), set(self.config["priority_analytes"]))
        self.assertRegex(provenance["sha256"], r"^[0-9a-f]{64}$")
        self.assertTrue(self.config["thermodynamics"]["silent_fallback_prohibited"])

    def test_missing_surrogate_table_fails_closed(self) -> None:
        config = json.loads(json.dumps(self.config))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "surrogate.json"
            source.write_text("{}", encoding="utf-8")
            config["thermodynamics"]["surrogate_source"] = "surrogate.json"
            with self.assertRaisesRegex(ValueError, "surrogate table is absent"):
                load_partition_surrogates(config, root)

    def test_observation_contract_matches_pilot_2026_tables(self) -> None:
        counts = self.config["observed_and_censored_counts"]
        self.assertEqual(counts["ethyl_acetate"]["condensate_observed"], 0)
        self.assertEqual(counts["ethyl_acetate"]["condensate_left_censored"], 33)
        self.assertEqual(counts["ethyl_octanoate"]["condensate_observed"], 19)
        self.assertEqual(counts["isoamyl_acetate"]["wine_left_censored"], 6)

    def test_weak_bilateral_direction_yields_conditional_gate(self) -> None:
        growth = "formation_growth_ug_per_g_sugar"
        loss = "effective_loss_scale"
        fits = [
            _fit("ethyl_acetate", [growth, loss], loss_identified=False),
            _fit("ethyl_octanoate", [growth]),
            _fit("isoamyl_acetate", [growth]),
        ]
        gate = evaluate_gate(fits)
        self.assertEqual(gate["verdict"], "PASS_CONDITIONAL")
        self.assertFalse(gate["all_parameter_directions_identified"])
        self.assertFalse(gate["profiles_for_physical_execution"])
        self.assertEqual(gate["weak_parameter_directions"]["ethyl_acetate"], [growth, loss])


if __name__ == "__main__":
    unittest.main()
