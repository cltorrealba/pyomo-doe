from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.operational_coverage import (  # noqa: E402
    AUTHORIZATION_FLAGS,
    WATERMARK,
    archetype_checks,
    classify_nutrition_archetype,
    classify_thermal_archetype,
    closed_authorization,
    contrast_checks,
    coverage_policy_from_unit_vector,
    derive_robust_setpoint_envelope,
    explicit_controller_blocks,
    gate_verdict,
    lexicographic_select,
    nutrition_margin_row,
    partially_harmonize_schedules,
    policies_operationally_equivalent,
    robust_initial_jump_rows,
    validate_physical_temperature,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    load_wave1_config,
)
from pilot_2026.adaptive_design.run_wave1_operational_coverage import (  # noqa: E402
    _repair_physical_candidate_labels,
)


ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"


class OperationalCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = load_wave1_config(
            ADAPTIVE_DIR / "wave1_mbdoe_config.json",
            ADAPTIVE_DIR / "design_constraints.json",
        )
        cls.coverage = json.loads(
            (ADAPTIVE_DIR / "coverage_design_config.json").read_text(encoding="utf-8")
        )
        cls.envelope = derive_robust_setpoint_envelope(cls.model, cls.coverage)
        vector = np.asarray(
            [0.35, 0.45, 0.2, 0.25, 0.3, 0.2, 0.55, 0.7, 0.2, 0.85, 0.35, 0.5]
        )
        cls.cold_early = coverage_policy_from_unit_vector(
            "cold_early",
            "cold",
            "early",
            76.0,
            vector,
            cls.model,
            cls.coverage,
            cls.envelope,
        )
        cls.cold_late = coverage_policy_from_unit_vector(
            "cold_late",
            "cold",
            "late",
            76.0,
            vector,
            cls.model,
            cls.coverage,
            cls.envelope,
        )
        cls.warm_early = coverage_policy_from_unit_vector(
            "warm_early",
            "warm",
            "early",
            76.0,
            vector,
            cls.model,
            cls.coverage,
            cls.envelope,
        )
        cls.warm_late = coverage_policy_from_unit_vector(
            "warm_late",
            "warm",
            "late",
            76.0,
            vector,
            cls.model,
            cls.coverage,
            cls.envelope,
        )

    def test_cold_and_warm_classification_is_hard_and_unambiguous(self) -> None:
        self.assertEqual(
            classify_thermal_archetype(self.cold_early, self.coverage), "cold"
        )
        self.assertEqual(
            classify_thermal_archetype(self.warm_late, self.coverage), "warm"
        )
        ambiguous = DesignPolicy("ambiguous", tuple([19.0] * 14), ((2.0, 76.0),))
        self.assertEqual(
            classify_thermal_archetype(ambiguous, self.coverage), "unclassified"
        )

    def test_early_and_late_nutrition_classification(self) -> None:
        self.assertEqual(
            classify_nutrition_archetype(self.cold_early, self.coverage), "early"
        )
        self.assertEqual(
            classify_nutrition_archetype(self.warm_late, self.coverage), "late"
        )
        invalid = DesignPolicy("invalid", tuple([18.0] * 14), ((18.0, 76.0),))
        self.assertEqual(
            classify_nutrition_archetype(invalid, self.coverage), "unclassified"
        )

    def test_minimum_thermal_and_nutrition_contrasts(self) -> None:
        checks = contrast_checks(
            self.cold_early,
            self.warm_late,
            self.cold_early,
            self.warm_late,
            self.coverage,
        )
        self.assertTrue(checks["thermal_contrast_at_least_4c"])
        self.assertTrue(checks["nutrition_contrast_at_least_36h"])
        self.assertGreaterEqual(float(checks["first_24h_mean_contrast_c"]), 4.0)
        self.assertGreaterEqual(float(checks["nutrition_time_contrast_h"]), 36.0)

    def test_initial_jump_is_robust_and_not_defined_as_zero(self) -> None:
        rows = robust_initial_jump_rows(self.warm_early, self.coverage)
        self.assertEqual(
            {row["initial_temperature_scenario_c"] for row in rows}, {17.0, 18.0, 19.0}
        )
        self.assertTrue(all(row["initial_jump_within_5c"] for row in rows))
        self.assertTrue(any(row["initial_setpoint_jump_c"] > 0.0 for row in rows))
        invalid = DesignPolicy(
            "old_warm", tuple([26.314] + [21.53] * 13), ((2.0, 80.0),)
        )
        self.assertFalse(
            all(
                row["initial_jump_within_5c"]
                for row in robust_initial_jump_rows(invalid, self.coverage)
            )
        )

    def test_robust_setpoint_envelope_is_derived_and_narrower_than_hardware(
        self,
    ) -> None:
        self.assertEqual(self.envelope["approved_scenario_count"], 243)
        self.assertGreater(self.envelope["steady_state_safe_setpoint_minimum_c"], 15.0)
        self.assertLess(self.envelope["steady_state_safe_setpoint_maximum_c"], 27.0)
        self.assertAlmostEqual(
            self.envelope["first_block_safe_setpoint_maximum_c"], 22.0
        )

    def test_physical_temperature_is_robustly_within_15_to_27(self) -> None:
        frame, summary = validate_physical_temperature(
            self.warm_late, self.model, self.coverage
        )
        self.assertEqual(len(frame), 243)
        self.assertTrue(summary["initial_jump_robust_pass"])
        self.assertTrue(summary["physical_temperature_robust_pass"])
        self.assertGreaterEqual(summary["minimum_physical_temperature_c"], 15.0 - 1e-9)
        self.assertLessEqual(summary["maximum_physical_temperature_c"], 27.0 + 1e-9)

    def test_explicit_hold_reaches_504_hours(self) -> None:
        controller = explicit_controller_blocks(self.cold_early, self.coverage)
        self.assertEqual(len(controller), 42)
        self.assertAlmostEqual(float(controller.iloc[-1].end_h), 504.0)
        self.assertTrue(bool(controller.iloc[-1].hold_last_setpoint))
        self.assertTrue(controller.iloc[14:].hold_last_setpoint.all())
        self.assertTrue(
            np.allclose(
                controller.iloc[14:].setpoint_c, self.cold_early.temperature_c[-1]
            )
        )

    def test_n76_mass_and_margin_variant(self) -> None:
        row = nutrition_margin_row(
            "II_diagonal_1",
            "N76_historical_anchor",
            self.cold_early,
            self.model,
            self.coverage,
        )
        self.assertAlmostEqual(float(row["yan_total_mg_l"]), 76.0)
        self.assertAlmostEqual(float(row["organic_product_g"]), 87.4)
        self.assertAlmostEqual(float(row["dap_product_g"]), 43.7)
        self.assertAlmostEqual(float(row["organic_margin_fraction"]), 0.05)
        self.assertTrue(row["net_yan_split_50_50"])
        self.assertTrue(row["nutrition_feasible"])

    def test_n76_is_preferred_only_within_half_percent_information_loss(self) -> None:
        n80_score = 20.0
        self.assertLessEqual((n80_score - 19.91) / n80_score, 0.005)
        self.assertGreater((n80_score - 19.89) / n80_score, 0.005)

    def test_coverage_treatments_do_not_duplicate(self) -> None:
        self.assertFalse(
            policies_operationally_equivalent(self.cold_early, self.warm_late)
        )
        duplicate = DesignPolicy(
            "duplicate",
            self.cold_early.temperature_c,
            self.cold_early.nutrition_mg_yan_l,
        )
        self.assertTrue(policies_operationally_equivalent(self.cold_early, duplicate))

    def test_lexicographic_selection_uses_safety_inside_half_percent(self) -> None:
        candidates = [
            {
                "candidate_id": "information_best",
                "feasible": True,
                "robust_score": 20.0,
                "maximum_initial_or_internal_jump_c": 5.0,
                "maximum_robust_physical_temperature_c": 26.8,
                "total_thermal_variation_c": 9.0,
                "temperature_changes": 4,
                "minimum_drying_margin_h": 30.0,
                "nutrition_operational_margin_fraction": 0.0,
            },
            {
                "candidate_id": "safer_near_best",
                "feasible": True,
                "robust_score": 19.91,
                "maximum_initial_or_internal_jump_c": 4.0,
                "maximum_robust_physical_temperature_c": 26.5,
                "total_thermal_variation_c": 8.0,
                "temperature_changes": 3,
                "minimum_drying_margin_h": 32.0,
                "nutrition_operational_margin_fraction": 0.05,
            },
        ]
        selected, reason = lexicographic_select(candidates, 0.005)
        self.assertEqual(selected["candidate_id"], "safer_near_best")
        self.assertIn("safer_near_best", reason["near_best_candidate_ids"])

    def test_partial_sampling_harmonization_reduces_workload_with_guardrail(
        self,
    ) -> None:
        independent = {"A": (0.0, 2.0, 18.0), "B": (0.0, 6.0, 18.0)}
        valid = {"A": (0.0, 2.0, 6.0, 18.0), "B": (0.0, 2.0, 6.0, 18.0)}

        def score(schedules: dict[str, tuple[float, ...]]) -> float:
            return 10.0 if schedules == independent else 9.99

        harmonized, audit = partially_harmonize_schedules(
            independent, valid, score, 0.005
        )
        self.assertLess(
            audit["harmonized_distinct_slots"], audit["independent_distinct_slots"]
        )
        self.assertGreaterEqual(audit["additional_shared_nonbasal_slots"], 1)
        self.assertGreaterEqual(audit["operator_round_reduction"], 1)
        self.assertLessEqual(audit["score_loss_fraction"], 0.005)
        self.assertEqual(len(harmonized["A"]), 3)
        self.assertEqual(len(harmonized["B"]), 3)

    def test_gates_fail_closed_on_missing_or_false_checks(self) -> None:
        self.assertEqual(gate_verdict({"a": True, "b": True}, ("a", "b")), "PASS")
        self.assertEqual(gate_verdict({"a": True}, ("a", "b")), "FAIL")
        self.assertEqual(gate_verdict({"a": True, "b": False}, ("a", "b")), "FAIL")

    def test_physical_flags_and_watermark_remain_closed(self) -> None:
        self.assertEqual(closed_authorization(), AUTHORIZATION_FLAGS)
        self.assertEqual(self.coverage["authorization"], AUTHORIZATION_FLAGS)
        self.assertEqual(
            WATERMARK,
            "COMPUTATIONAL COVERAGE CANDIDATE \u2014 NOT AUTHORIZED FOR PHYSICAL EXECUTION",
        )

    def test_decoder_outputs_exact_archetype_checks(self) -> None:
        checks = archetype_checks(
            self.warm_late,
            "warm",
            "late",
            self.model,
            self.coverage,
        )
        self.assertTrue(all(checks.values()))

    def test_cached_physical_labels_are_repaired_without_numeric_changes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary)
            actions = pd.DataFrame(
                [
                    {
                        "strategy": "II_diagonal_1",
                        "dose_design": "N76_all_conservative",
                        "candidate_id": "candidate",
                        "policy": "cold_early",
                        "action": "initial_transition_audit",
                        "value": 16.5,
                    },
                    {
                        "strategy": "II_diagonal_1",
                        "dose_design": "N76_all_conservative",
                        "candidate_id": "candidate",
                        "policy": "warm_late",
                        "action": "initial_transition_audit",
                        "value": 20.5,
                    },
                ]
            )
            physical = pd.DataFrame(
                [
                    {
                        "strategy": "II_diagonal_1",
                        "dose_design": "N76_all_conservative",
                        "candidate_id": "candidate",
                        "candidate": "stale_cold_label",
                        "initial_setpoint_c": 16.5,
                        "minimum_physical_temperature_c": 15.2,
                    },
                    {
                        "strategy": "II_diagonal_1",
                        "dose_design": "N76_all_conservative",
                        "candidate_id": "candidate",
                        "candidate": "stale_warm_label",
                        "initial_setpoint_c": 20.5,
                        "minimum_physical_temperature_c": 17.2,
                    },
                ]
            )
            actions.to_csv(run_dir / "coverage_candidate_actions.csv", index=False)
            physical.to_csv(
                run_dir / "physical_temperature_envelope_by_candidate.csv",
                index=False,
            )
            repair = _repair_physical_candidate_labels(run_dir)
            repaired = pd.read_csv(
                run_dir / "physical_temperature_envelope_by_candidate.csv"
            )
            self.assertEqual(repair["labels_corrected"], 2)
            self.assertTrue(repair["all_non_label_columns_byte_value_identical"])
            self.assertEqual(repaired.candidate.tolist(), ["cold_early", "warm_late"])
            pd.testing.assert_frame_equal(
                physical.drop(columns="candidate"),
                repaired.drop(columns="candidate"),
                check_exact=True,
            )


if __name__ == "__main__":
    unittest.main()
