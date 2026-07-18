from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    decode_policy_vector,
    finite_difference_sensitivity,
    load_wave1_config,
    nominal_observation_scale,
    policy_to_canonical_vector,
    temperature_profile_metrics,
    vector_bounds,
)
from shared import run_new_must_glycerol_estimability_doe as model  # noqa: E402


class Wave1RequalificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_wave1_config(ADAPTIVE_DIR / "wave1_mbdoe_config.json")

    @staticmethod
    def _batch(pulses: tuple[tuple[float, float], ...]) -> model.BatchData:
        return model.BatchData(
            medium="test",
            batch="causal_pulse",
            time=np.asarray([0.0, 2.0]),
            temperature_c=np.asarray([18.0, 18.0]),
            pulses={"N": pulses, "G": tuple(), "F": tuple(), "E": tuple(), "X": tuple()},
            observations={},
            initials={"X": 0.0, "Xd": 0.0, "N": 0.2, "G": 80.0, "F": 80.0, "E": 0.0, "Gly": 0.0},
        )

    def test_nominal_sigma_is_frozen_in_scaled_sensitivity(self) -> None:
        centre = np.asarray([10.0])
        nominal = np.asarray([10.0])
        sigma = nominal_observation_scale(
            nominal,
            sample_count=1,
            config=self.config,
        )
        sensitivity, methods = finite_difference_sensitivity(
            lambda theta: np.asarray([theta[0]]), centre, 0.02, sigma
        )
        self.assertEqual(methods, ("central",))
        self.assertAlmostEqual(float(sensitivity[0, 0]), 1.0 / float(sigma[0]), places=12)
        self.assertNotEqual(float(sensitivity[0, 0]), 0.0)

    def test_second_order_unilateral_difference_at_bounds(self) -> None:
        sensitivity, methods = finite_difference_sensitivity(
            lambda theta: np.asarray([theta[0]]),
            np.asarray([0.0]),
            0.02,
            np.asarray([1.0]),
            lower=np.asarray([0.0]),
            upper=np.asarray([1.0]),
        )
        self.assertEqual(methods, ("forward_second_order",))
        self.assertAlmostEqual(float(sensitivity[0, 0]), 1.0, places=12)

    def test_temperature_policy_is_canonical_and_has_bounded_changepoints(self) -> None:
        bounds = vector_bounds(self.config)
        raw = 0.5 * (bounds[:, 0] + bounds[:, 1])
        layout_changes = int(self.config["future_process"]["maximum_temperature_changes"])
        raw[1 : 1 + layout_changes] = 1.0
        decoded = decode_policy_vector("temperature", raw, self.config)
        metrics = temperature_profile_metrics(decoded.policy, self.config)
        self.assertLessEqual(metrics["temperature_changes"], layout_changes)
        roundtrip = decode_policy_vector("temperature", decoded.canonical_vector, self.config)
        np.testing.assert_allclose(
            roundtrip.canonical_vector, decoded.canonical_vector, atol=0.0, rtol=0.0
        )

    def test_duplicate_pulses_are_repaired_without_losing_dose(self) -> None:
        bounds = vector_bounds(self.config)
        raw = bounds[:, 0].copy()
        raw[0] = 18.0
        changes = int(self.config["future_process"]["maximum_temperature_changes"])
        pulses = int(self.config["nutrition"]["maximum_pulses"])
        active_start = 1 + 3 * changes
        time_start = active_start + pulses
        amount_start = time_start + pulses
        raw[active_start : active_start + pulses] = 1.0
        raw[time_start : time_start + pulses] = 2.0
        raw[amount_start : amount_start + pulses] = [20.0, 30.0, 40.0]
        decoded = decode_policy_vector("pulses", raw, self.config)
        schedule = decoded.policy.nutrition_mg_yan_l
        self.assertEqual(len({time for time, _ in schedule}), 3)
        self.assertAlmostEqual(sum(amount for _, amount in schedule), 90.0)
        self.assertTrue(
            any(row["action"] == "relocate_duplicate_nutrition_time" for row in decoded.repair_log)
        )
        np.testing.assert_allclose(
            policy_to_canonical_vector(decoded.policy, self.config),
            decoded.canonical_vector,
            atol=0.0,
            rtol=0.0,
        )

    def test_pulse_is_strictly_causal_and_conserves_exact_dose(self) -> None:
        batch = self._batch(((1.0, 0.08),))
        grid = np.asarray([0.0, 0.5, 1.0, 1.0001, 2.0])
        simulated = model.simulate(batch, model.DEFAULT_THETA, grid)
        self.assertIsNotNone(simulated)
        self.assertAlmostEqual(float(simulated.loc[0.5, "N"]), 0.2, places=12)
        self.assertAlmostEqual(float(simulated.loc[1.0, "N"]), 0.2, places=12)
        self.assertAlmostEqual(float(simulated.loc[1.0001, "N"]), 0.28, places=10)

    def test_sample_action_order_and_pulse_at_zero_are_explicit(self) -> None:
        batch = self._batch(((0.0, 0.08),))
        grid = np.asarray([0.0, 0.1])
        before = model.simulate(
            batch, model.DEFAULT_THETA, grid, sample_event_order="sample_before_action"
        )
        after = model.simulate(
            batch, model.DEFAULT_THETA, grid, sample_event_order="action_before_sample"
        )
        self.assertAlmostEqual(float(before.loc[0.0, "N"]), 0.2, places=12)
        self.assertAlmostEqual(float(after.loc[0.0, "N"]), 0.28, places=12)
        self.assertAlmostEqual(float(before.loc[0.1, "N"]), 0.28, places=10)

    def test_close_pulses_and_integration_mesh_preserve_mass(self) -> None:
        batch = self._batch(((1.0, 0.03), (1.01, 0.05)))
        grid = np.asarray([0.0, 0.99, 1.005, 1.02, 2.0])
        fine = model.simulate(batch, model.DEFAULT_THETA, grid, integration_max_step_h=0.05)
        coarse = model.simulate(batch, model.DEFAULT_THETA, grid, integration_max_step_h=0.5)
        self.assertAlmostEqual(float(fine.loc[1.02, "N"]), 0.28, places=9)
        self.assertAlmostEqual(float(coarse.loc[1.02, "N"]), 0.28, places=9)
        self.assertLess(abs(float(fine.loc[2.0, "N"] - coarse.loc[2.0, "N"])), 1e-10)


if __name__ == "__main__":
    unittest.main()
