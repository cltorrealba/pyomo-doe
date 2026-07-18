from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from unittest import mock
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    decode_policy_vector,
    finite_difference_sensitivity,
    load_wave1_config,
    nominal_observation_scale,
    physical_observation_vector,
    policy_to_canonical_vector,
    temperature_profile_metrics,
    vector_bounds,
)
from pilot_2026.adaptive_design.optimize_wave1_sampling_and_plots import (  # noqa: E402
    _capture_intervals,
    _nutrition_translation,
    _operational_conflicts,
    _tank_randomization,
    capture_constraints_approved,
    optimized_schedule_is_acceptable,
    policy_actions_before_drying,
    robust_score,
    sampling_gate_verdict,
)
from pilot_2026.adaptive_design.pilot_aroma_calibration import AromaForcing  # noqa: E402
from shared import run_new_must_glycerol_estimability_doe as model  # noqa: E402
from pilot_2026.adaptive_design import run_artifacts  # noqa: E402


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

    def test_manifest_verification_uses_windows_extended_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            run_dir = Path(temporary_directory)
            output_path = run_dir / "actuator_robustness_validation.csv"
            payload = b"member,feasible\n0,true\n"
            output_path.write_bytes(payload)
            output_key = run_artifacts.relative_or_absolute(output_path)
            manifest_path = run_dir / "run_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "run_id": "path-test",
                        "outputs": {
                            output_key: {"sha256": hashlib.sha256(payload).hexdigest()}
                        },
                    }
                ),
                encoding="utf-8",
            )
            original = run_artifacts.filesystem_path
            with mock.patch.object(
                run_artifacts, "filesystem_path", wraps=original
            ) as guarded_path:
                verification = run_artifacts.verify_manifest_output(
                    run_dir, output_path.name
                )
            guarded_inputs = {
                Path(call.args[0]).resolve() for call in guarded_path.call_args_list
            }
            self.assertIn(manifest_path.resolve(), guarded_inputs)
            self.assertIn(output_path.resolve(), guarded_inputs)
            self.assertEqual(verification["verification_mode"], "byte_exact")

    def test_campaign_state_keeps_all_physical_release_flags_closed(self) -> None:
        state = json.loads(
            (ADAPTIVE_DIR / "campaign_state.json").read_text(encoding="utf-8")
        )
        self.assertFalse(state["profiles_for_physical_execution"])
        self.assertFalse(state["physical_execution_authorized"])
        self.assertFalse(state["executable_schedule_issued"])
        self.assertEqual(state["tank_assignments"], [])
        self.assertFalse(state["tank_randomization_authorized"])

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

    def test_sampling_optimizer_rejects_score_regression_and_falls_back(self) -> None:
        accepted, reason = optimized_schedule_is_acceptable(
            np.asarray([1.0, 1.0, 1.0]),
            np.asarray([2.0, 2.0, 2.0]),
            self.config,
        )
        self.assertFalse(accepted)
        self.assertEqual(reason, "optimized_full_ensemble_score_below_preliminary")

    def test_robust_objective_is_identical_to_declared_formula(self) -> None:
        gains = np.asarray([1.0, 2.0, 3.0, 10.0])
        expected = 0.75 * float(np.median(gains)) + 0.25 * float(np.quantile(gains, 0.1))
        self.assertAlmostEqual(robust_score(gains, self.config), expected, places=14)

    def test_sampling_gate_fails_when_critical_capture_check_fails(self) -> None:
        checks = {"score": True, "capture": False, "physical_lock": True}
        self.assertEqual(
            sampling_gate_verdict(checks, ("score", "capture")),
            "FAIL",
        )
        self.assertTrue(capture_constraints_approved(self.config))

    def test_sub_one_degree_temperature_change_is_merged(self) -> None:
        bounds = vector_bounds(self.config)
        raw = bounds[:, 0].copy()
        raw[0] = 18.0
        changes = int(self.config["future_process"]["maximum_temperature_changes"])
        raw[1] = 1.0
        raw[1 + changes] = 2.0
        raw[1 + 2 * changes] = 18.5
        decoded = decode_policy_vector("subthreshold", raw, self.config)
        self.assertTrue(all(value == 18.0 for value in decoded.policy.temperature_c))
        self.assertTrue(
            any(
                row["action"] == "merge_subthreshold_temperature_segment"
                for row in decoded.repair_log
            )
        )

    def test_baseline_and_capture_intervals_are_distinct_observation_types(self) -> None:
        forcing = AromaForcing(
            experiment_id="capture_test",
            time_h=np.asarray([0.0, 1.0, 2.0, 3.0]),
            growth_fraction=np.asarray([0.5, 0.5, 0.5, 0.5]),
            sugar_uptake_g_l_h=np.asarray([1.0, 1.0, 1.0, 1.0]),
            loss_basis_h_inv=np.asarray([0.1, 0.1, 0.1, 0.1]),
            initial_concentration_ug_l=10.0,
            volume_l=230.0,
            trap_efficiency=1.0,
        )
        samples = np.linspace(0.0, 3.0, 10)
        observations = physical_observation_vector(
            forcing, np.log(np.asarray([1.0, 1.0, 1.0])), samples
        )
        self.assertEqual(len(observations[:10]), 10)
        self.assertEqual(len(observations[10:]), 9)
        self.assertAlmostEqual(float(observations[0]), 10.0)
        self.assertGreaterEqual(float(observations[10:].sum()), 0.0)

    def test_capture_intervals_begin_at_zero_and_conserve_mass(self) -> None:
        schedules = {"policy": tuple(float(value) for value in range(10))}
        time = np.arange(10.0)
        cache = {}
        for member in (0, 1):
            cumulative = {
                species: time * float(member + 1)
                for species in ("ethyl_acetate", "ethyl_octanoate", "isoamyl_acetate")
            }
            forcings = {
                species: SimpleNamespace(time_h=time) for species in cumulative
            }
            cache[(member, "policy")] = SimpleNamespace(
                prepared=SimpleNamespace(forcings=forcings),
                nominal_captured=cumulative,
            )
        intervals = _capture_intervals(schedules, self.config, cache, [0, 1])
        self.assertTrue(intervals.groupby(["policy", "species"]).size().eq(9).all())
        self.assertTrue(
            intervals.groupby(["policy", "species"])["start_h"].min().eq(0.0).all()
        )
        self.assertTrue(intervals["mass_conservation_max_abs_error_ug"].eq(0.0).all())

    def test_owner_approved_sampling_collisions_are_not_conflicts(self) -> None:
        policies = tuple(
            DesignPolicy(name, tuple([18.0] * 14), ((0.0, 20.0),))
            for name in ("anchor", "A", "B")
        )
        schedules = {
            policy.name: tuple(float(value) for value in range(10)) for policy in policies
        }
        conflicts = _operational_conflicts(schedules, policies, self.config)
        self.assertFalse(conflicts["status"].eq("unresolved").any())

    def test_nutrition_50_50_is_encoded_but_product_grams_stay_blocked(self) -> None:
        policy = DesignPolicy("A", tuple([18.0] * 14), ((24.0, 80.0),))
        translation = _nutrition_translation((policy,), self.config)
        self.assertEqual(
            translation.iloc[0]["mix_policy"], "50_50_net_YAN_contribution"
        )
        self.assertAlmostEqual(
            float(translation.iloc[0]["organic_yan_contribution_mg_l"]), 40.0
        )
        self.assertAlmostEqual(
            float(translation.iloc[0]["dap_yan_contribution_mg_l"]), 40.0
        )
        self.assertTrue(np.isnan(float(translation.iloc[0]["organic_product_g"])))
        self.assertTrue(bool(translation.iloc[0]["translation_blocker"]))

    def test_tank_randomization_is_reproducible_but_not_authorized(self) -> None:
        policies = tuple(
            DesignPolicy(name, tuple([18.0] * 14), tuple())
            for name in ("anchor", "A", "B")
        )
        frozen = datetime(2026, 7, 18, tzinfo=timezone.utc).isoformat()
        first = _tank_randomization(policies, self.config, frozen)
        second = _tank_randomization(policies, self.config, frozen)
        self.assertEqual(first["tank"].tolist(), second["tank"].tolist())
        self.assertEqual(first["tank"].nunique(), 3)
        self.assertFalse(first["authorized_for_physical_execution"].any())

    def test_actions_after_drying_are_rejected(self) -> None:
        policy = DesignPolicy(
            "late_action",
            tuple([18.0] * 13 + [20.0]),
            tuple(),
        )
        self.assertFalse(
            policy_actions_before_drying(
                policy,
                np.asarray([100.0] * 64),
                self.config,
            )
        )

    def test_action_drying_margin_is_exactly_24_hours(self) -> None:
        policy = DesignPolicy("margin", tuple([18.0] * 14), ((72.0, 20.0),))
        self.assertTrue(
            policy_actions_before_drying(policy, np.asarray([96.0] * 64), self.config)
        )
        self.assertFalse(
            policy_actions_before_drying(policy, np.asarray([95.9] * 64), self.config)
        )

    def test_minimum_search_budget_and_five_seeds_are_declared(self) -> None:
        search = self.config["search"]
        self.assertEqual(len(self.config["independent_seeds"]), 5)
        self.assertGreaterEqual(search["particles"], 24)
        self.assertGreaterEqual(search["maximum_iterations"], 20)
        self.assertGreaterEqual(search["top_k_candidates"], 5)
        self.assertGreaterEqual(search["top_k_local_refinement"], 3)

    def test_fim_validation_scope_covers_policies_members_and_actuator(self) -> None:
        self.assertEqual(
            set(self.config["fim_validation"]["validation_scope"]),
            {
                "anchor",
                "candidate_A",
                "candidate_B",
                "central_member",
                "lowest_information_member",
                "highest_information_member",
                "critical_drying_member",
                "critical_actuator_scenario",
            },
        )

    def test_hybrid_gate_declares_every_owner_critical_check(self) -> None:
        source = (ADAPTIVE_DIR / "run_wave1_hybrid_search.py").read_text(encoding="utf-8")
        for name in (
            "adapter_gate_pass",
            "source_hashes_verified",
            "corrected_fim_stability_pass",
            "pso_three_independent_seeds_completed",
            "pso_objectives_finite",
            "top_k_revalidated_on_64_members",
            "actual_objective_evaluated",
            "local_refinement_qualified",
            "full_ensemble_completion_probability_at_least_95pct",
            "canonical_policies_have_no_duplicate_pulses",
            "actions_respect_24h_drying_margin",
            "temperature_changes_respect_minimum_1C",
        ):
            self.assertIn(f'"{name}"', source)

    def test_all_required_v2_figure_names_are_declared(self) -> None:
        source = (ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py").read_text(
            encoding="utf-8"
        )
        for name in (
            "candidate_profiles_v2.png",
            "executed_temperature_ensemble_v2.png",
            "aroma_predictions_v2.png",
            "sampling_schedule_v2.png",
            "capture_intervals_v2.png",
            "information_gain_distribution_v2.png",
            "paired_information_comparison_v2.png",
            "drying_margin_validation_v2.png",
            "actuator_robustness_v2.png",
            "operational_summary_v2.png",
        ):
            self.assertIn(name, source)

    def test_actuator_uncertainty_contains_all_nine_empirical_taus(self) -> None:
        values = self.config["future_process"]["temperature_actuator"]["empirical_tau_h"]
        self.assertEqual(len(values), 9)
        self.assertAlmostEqual(min(values), 0.21148207178811254)
        self.assertAlmostEqual(max(values), 0.45812979693594796)

    def test_sampling_source_is_explicit_not_latest_glob(self) -> None:
        source = (ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('add_argument("--source-search-run", required=True)', source)
        self.assertNotIn("sorted((RESULT_ROOT", source)


if __name__ == "__main__":
    unittest.main()
