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
from pilot_2026.adaptive_design.pilot_calibration import (  # noqa: E402
    build_batches,
    load_tables,
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
        self.assertEqual(self.config["yan"]["fixed_historical_pulse_mg_l"], 80.0)
        self.assertNotIn("delivery_fraction_bounds", self.config["yan"])
        self.assertIn("scale_bounds", self.config["oculyze"])
        self.assertIn("absolute_floor", self.config["observation_error"])
        self.assertIn("likelihood", self.config["censoring"])
        self.assertEqual(
            self.config["co2_fit"]["ess_method"],
            "per_run_AR1_from_final_standardized_residuals",
        )

    def test_calibration_uses_fixed_80_mg_l_for_later_historical_pulses(self) -> None:
        state = json.loads((ADAPTIVE_DIR / "campaign_state.json").read_text(encoding="utf-8"))
        model_run = FERMENTATION_DIR.parent / state["latest_model_dataset_run"]
        tables = load_tables(model_run)
        batches = build_batches(
            tables,
            {"oculyze_biomass_scale_kg_m3_per_million_cells_ml": 0.03},
            self.config,
        )
        pulses = [amount for batch in batches for _, amount in batch.pulses["N"]]
        self.assertEqual(len(pulses), 9)
        self.assertTrue(all(abs(amount - 0.08) < 1e-12 for amount in pulses))

    def test_existing_pilot_2025_integrated_result_is_never_an_official_prior(self) -> None:
        policy = self.config["pilot_2025_policy"]
        self.assertFalse(policy["use_existing_integrated_result_as_prior"])
        self.assertIn("regenerate", policy["required_action_before_use"])

    def test_design_constraints_are_owner_approved_and_complete(self) -> None:
        constraints = json.loads(
            (ADAPTIVE_DIR / "design_constraints.json").read_text(encoding="utf-8")
        )
        self.assertEqual(constraints["review_status"], "owner_approved_2026-07-17")
        self.assertFalse(constraints["approval"]["physical_execution_authorized"])
        self.assertEqual(constraints["fail_closed_fields"], [])
        self.assertEqual(constraints["temperature"]["maximum_active_segments"], 42)
        self.assertEqual(constraints["nutrition"]["maximum_pulses"], 3)
        self.assertEqual(constraints["nutrition"]["maximum_total_yan_mg_l"], 232.0)
        self.assertEqual(
            constraints["sampling_and_capture"]["capture_stages_c"], [0.0, -40.0]
        )

    def test_aroma_gate_is_fail_closed_for_thermodynamic_fallbacks(self) -> None:
        aroma = json.loads(
            (ADAPTIVE_DIR / "aroma_calibration_config.json").read_text(encoding="utf-8")
        )
        self.assertTrue(aroma["thermodynamics"]["silent_fallback_prohibited"])
        self.assertEqual(
            aroma["observed_and_censored_counts"]["ethyl_acetate"]["condensate_observed"],
            0,
        )
        self.assertIn("cannot be declared", aroma["identifiability_warning"]["ethyl_acetate_loss"])


if __name__ == "__main__":
    unittest.main()
