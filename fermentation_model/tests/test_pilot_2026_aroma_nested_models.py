from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    AromaForcing,
    LOSS_MODEL_EQUILIBRIUM,
)
from pilot_2026.adaptive_design.pilot_aroma_nested_models import (  # noqa: E402
    ASSUMED_COMPLETE,
    ASSUMED_COMPLETE_RESERVOIR,
    BASELINE,
    CAPTURE,
    COMBINED,
    DELAYED,
    ETHANOL_CAPTURE,
    ETHANOL_CAPTURE_RESERVOIR,
    RESERVOIR,
    _compile_observations,
    pulse_activation,
    simulate_nested,
)


def _forcing(*, uptake: float = 1.0, loss: float = 0.2) -> AromaForcing:
    time = np.arange(0.0, 11.0, 1.0)
    return AromaForcing(
        experiment_id="synthetic",
        time_h=time,
        growth_fraction=np.full(len(time), 0.5),
        sugar_uptake_g_l_h=np.full(len(time), uptake),
        loss_basis_h_inv=np.full(len(time), loss),
        initial_concentration_ug_l=10.0,
        volume_l=2.0,
        trap_efficiency=0.8,
        co2_rate_g_l_h=np.ones(len(time)),
        partition_basis_l_g=np.full(len(time), loss),
        gas_turnover_h_inv=np.ones(len(time)),
        temperature_c=np.full(len(time), 20.0),
        ethanol_g_l=np.full(len(time), 50.0),
        total_sugar_g_l=np.full(len(time), 100.0),
        loss_model=LOSS_MODEL_EQUILIBRIUM,
    )


class Pilot2026AromaNestedModelTests(unittest.TestCase):
    def test_wine_only_calibration_excludes_condensate_rows(self) -> None:
        wine = pd.DataFrame(
            {
                "experiment_id": ["r1", "r1"],
                "time_h": [0.0, 1.0],
                "model_observation_type": ["observed", "observed"],
                "observed_value": [1.0, 2.0],
                "upper_bound": [np.nan, np.nan],
            }
        )
        condensate = pd.DataFrame(
            {
                "experiment_id": ["r1"],
                "time_h": [1.0],
                "capture_interval_h": [1.0],
                "mix_model_observation_type": ["observed"],
                "observed_value": [10.0],
                "upper_bound": [np.nan],
            }
        )
        config = {
            "error_model": {
                "calibration_domains": ["wine"],
                "calibration_weighting": "domain_mean",
                "wine_relative_sigma": 0.2,
                "wine_minimum_sigma_ug_l": 1.0,
                "condensate_relative_sigma": 0.25,
                "condensate_minimum_sigma_ug": 1.0,
            }
        }
        wine_rows, condensate_rows = _compile_observations(
            wine, condensate, config
        )
        self.assertEqual(len(wine_rows), 1)
        self.assertEqual(condensate_rows, [])

    def test_pulse_activation_is_gradual_and_peaks_at_tau(self) -> None:
        time = np.array([0.0, 2.0, 3.0, 4.0, 7.0, 11.0])
        activation = pulse_activation(time, pulse_time_h=3.0, activation_peak_h=4.0)
        self.assertTrue(np.allclose(activation[:3], 0.0))
        self.assertGreater(activation[3], 0.0)
        self.assertAlmostEqual(activation[4], 1.0)
        self.assertLess(activation[5], activation[4])

    def test_line_reservoir_delays_capture_and_conserves_mass(self) -> None:
        forcing = _forcing()
        direct = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0]), BASELINE
        )
        delayed = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0, 5.0]), RESERVOIR
        )
        self.assertLess(delayed.captured_ug[3], direct.captured_ug[3])
        self.assertGreater(delayed.line_inventory_ug[3], 0.0)
        self.assertLess(
            float(np.max(delayed.relative_mass_balance_error)), 1e-12
        )
        self.assertLess(
            float(np.max(direct.relative_mass_balance_error)), 1e-12
        )

    def test_instantaneous_capture_matches_volatilized_mass_times_efficiency(self) -> None:
        forcing = _forcing()
        simulation = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0]), BASELINE
        )
        self.assertTrue(
            np.allclose(
                simulation.captured_ug,
                forcing.trap_efficiency * simulation.cumulative_volatilized_ug,
            )
        )
        self.assertTrue(np.allclose(simulation.line_inventory_ug, 0.0))

    def test_assumed_complete_capture_has_no_fitted_efficiency(self) -> None:
        forcing = _forcing()
        simulation = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0]), ASSUMED_COMPLETE
        )
        self.assertTrue(
            np.allclose(
                simulation.captured_ug,
                simulation.cumulative_volatilized_ug,
            )
        )
        self.assertTrue(
            np.allclose(simulation.capture_efficiency_fraction, 1.0)
        )

    def test_assumed_complete_capture_reservoir_drains_all_line_mass(self) -> None:
        forcing = _forcing()
        simulation = simulate_nested(
            forcing,
            3.0,
            np.log([2.0, 4.0, 5.0]),
            ASSUMED_COMPLETE_RESERVOIR,
        )
        self.assertTrue(
            np.allclose(
                simulation.captured_ug,
                simulation.cumulative_drained_ug,
            )
        )
        self.assertGreater(float(simulation.line_inventory_ug.max()), 0.0)
        self.assertLess(
            float(np.max(simulation.relative_mass_balance_error)), 1e-12
        )

    def test_delayed_response_adds_only_postpulse_production(self) -> None:
        forcing = _forcing()
        baseline = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0]), BASELINE
        )
        delayed = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0, 8.0, 2.0]), DELAYED
        )
        before = forcing.time_h <= 3.0
        after = forcing.time_h > 3.0
        self.assertTrue(
            np.allclose(
                delayed.production_rate_ug_l_h[before],
                baseline.production_rate_ug_l_h[before],
            )
        )
        self.assertTrue(
            np.all(
                delayed.production_rate_ug_l_h[after]
                > baseline.production_rate_ug_l_h[after]
            )
        )

    def test_combined_model_states_are_nonnegative(self) -> None:
        simulation = simulate_nested(
            _forcing(),
            3.0,
            np.log([2.0, 4.0, 8.0, 2.0, 5.0]),
            COMBINED,
        )
        for values in (
            simulation.liquid_ug_l,
            simulation.captured_ug,
            simulation.line_inventory_ug,
            simulation.cumulative_production_ug,
            simulation.cumulative_drained_ug,
        ):
            self.assertTrue(np.all(values >= 0.0))

    def test_fitted_capture_efficiency_changes_observation_not_liquid_balance(self) -> None:
        forcing = _forcing()
        nominal = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0]), BASELINE
        )
        fitted = simulate_nested(
            forcing, 3.0, np.log([2.0, 4.0, 0.04]), CAPTURE
        )
        self.assertTrue(np.allclose(fitted.liquid_ug_l, nominal.liquid_ug_l))
        self.assertTrue(
            np.allclose(
                fitted.captured_ug,
                nominal.captured_ug * 0.04 / forcing.trap_efficiency,
            )
        )
        self.assertLess(
            float(np.max(fitted.relative_mass_balance_error)), 1e-12
        )

    def test_ethanol_dependent_capture_is_bounded_and_increases_with_ethanol(self) -> None:
        forcing = _forcing()
        forcing = AromaForcing(
            **{
                **forcing.__dict__,
                "ethanol_g_l": np.linspace(0.0, 100.0, len(forcing.time_h)),
            }
        )
        simulation = simulate_nested(
            forcing,
            3.0,
            np.log([2.0, 4.0, 0.01, 2.0]),
            ETHANOL_CAPTURE,
        )
        self.assertTrue(np.all(simulation.capture_efficiency_fraction > 0.0))
        self.assertTrue(np.all(simulation.capture_efficiency_fraction <= 1.0))
        self.assertTrue(np.all(np.diff(simulation.capture_efficiency_fraction) > 0.0))

    def test_ethanol_capture_reservoir_combines_both_states_and_conserves_mass(self) -> None:
        forcing = _forcing()
        forcing = AromaForcing(
            **{
                **forcing.__dict__,
                "ethanol_g_l": np.linspace(0.0, 100.0, len(forcing.time_h)),
            }
        )
        simulation = simulate_nested(
            forcing,
            3.0,
            np.log([2.0, 4.0, 0.01, 2.0, 5.0]),
            ETHANOL_CAPTURE_RESERVOIR,
        )
        self.assertGreater(float(simulation.line_inventory_ug.max()), 0.0)
        self.assertTrue(np.all(np.diff(simulation.capture_efficiency_fraction) > 0.0))
        self.assertLess(
            float(np.max(simulation.relative_mass_balance_error)), 1e-12
        )


if __name__ == "__main__":
    unittest.main()
