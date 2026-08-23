from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_DIR = FERMENTATION_DIR.parent
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    AromaForcing,
    LOSS_MODEL_DYNAMIC_TRANSFER,
    LOSS_MODEL_EQUILIBRIUM,
    _co2_gas_density_g_l,
    _partition_basis,
    evaluate_gate,
    load_partition_surrogates,
    transfer_diagnostics,
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

    def test_dynamic_transfer_uses_liquid_side_ntu_and_bounded_efficiency(self) -> None:
        density = _co2_gas_density_g_l(20.0)
        self.assertAlmostEqual(density, 1.830, places=2)
        forcing = AromaForcing(
            experiment_id="test",
            time_h=np.array([0.0, 1.0]),
            growth_fraction=np.zeros(2),
            sugar_uptake_g_l_h=np.zeros(2),
            loss_basis_h_inv=np.array([0.0025, 0.0025]),
            initial_concentration_ug_l=1.0,
            volume_l=1.0,
            trap_efficiency=1.0,
            co2_rate_g_l_h=np.array([density * 0.5, density * 0.5]),
            partition_basis_l_g=np.array([0.005, 0.005]),
            gas_turnover_h_inv=np.array([0.5, 0.5]),
            ethanol_g_l=np.array([50.0, 50.0]),
            loss_model=LOSS_MODEL_DYNAMIC_TRANSFER,
        )
        log_values = np.log([1.0, 1.0, 0.0025, 1.0])
        diagnostics = transfer_diagnostics(forcing, log_values)
        expected_efficiency = 1.0 - np.exp(-1.0)
        self.assertTrue(
            np.allclose(diagnostics["transfer_efficiency"], expected_efficiency)
        )
        self.assertTrue(
            np.allclose(
                diagnostics["loss_coefficient_h_inv"],
                0.005 * 0.5 * expected_efficiency,
            )
        )

    def test_equilibrium_model_has_no_fitted_transfer_penalty(self) -> None:
        forcing = AromaForcing(
            experiment_id="test",
            time_h=np.array([0.0, 1.0]),
            growth_fraction=np.zeros(2),
            sugar_uptake_g_l_h=np.zeros(2),
            loss_basis_h_inv=np.array([0.0025, 0.0025]),
            initial_concentration_ug_l=1.0,
            volume_l=1.0,
            trap_efficiency=1.0,
            co2_rate_g_l_h=np.ones(2),
            partition_basis_l_g=np.array([0.005, 0.005]),
            gas_turnover_h_inv=np.array([0.5, 0.5]),
            ethanol_g_l=np.array([50.0, 50.0]),
            loss_model=LOSS_MODEL_EQUILIBRIUM,
        )
        diagnostics = transfer_diagnostics(forcing, np.log([1.0, 1.0]))
        self.assertTrue(np.allclose(diagnostics["transfer_efficiency"], 1.0))
        self.assertTrue(np.allclose(diagnostics["loss_coefficient_h_inv"], 0.0025))
        self.assertTrue(diagnostics["mass_transfer_kla_h_inv"].isna().all())

    def test_morakul_partition_matches_reference_temperature_expression(self) -> None:
        model = {
            "model_kind": "morakul_ethanol_temperature",
            "F1": -3.13,
            "F2_l_g": -1.35e-2,
            "F3_kj_mol": 52.0,
            "F4_kj_l_mol_g": -3.60e-3,
            "reference_temperature_k": 293.15,
        }
        expected = np.exp(-3.13 - 1.35e-2 * 50.0)
        self.assertAlmostEqual(_partition_basis(model, 20.0, 50.0, 100.0), expected)

    def test_morakul_temperature_term_uses_published_energy_units(self) -> None:
        model = {
            "model_kind": "morakul_ethanol_temperature",
            "F1": -3.13,
            "F2_l_g": -1.35e-2,
            "F3_kj_mol": 52.0,
            "F4_kj_l_mol_g": -3.60e-3,
            "reference_temperature_k": 293.15,
        }
        ethanol = 50.0
        temperature_k = 288.15
        enthalpy = 52.0 - 3.60e-3 * ethanol
        expected_log_k = (
            -3.13
            - 1.35e-2 * ethanol
            - enthalpy
            / 8.314462618
            * (1000.0 / temperature_k - 1000.0 / 293.15)
        )
        self.assertAlmostEqual(
            _partition_basis(model, 15.0, ethanol, 100.0),
            np.exp(expected_log_k),
        )


if __name__ == "__main__":
    unittest.main()
