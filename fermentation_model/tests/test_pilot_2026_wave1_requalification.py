from __future__ import annotations

import hashlib
import gzip
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
    actuator_temperature_trajectory,
    allowed_sampling_times,
    decode_policy_vector,
    finite_difference_sensitivity,
    load_wave1_config,
    nominal_observation_scale,
    nutrition_policy_metrics,
    nutrition_product_masses,
    physical_observation_vector,
    policy_to_canonical_vector,
    temperature_profile_metrics,
    vector_bounds,
)
from pilot_2026.adaptive_design.final_search_logic import (  # noqa: E402
    approved_actuator_scenarios,
    practical_convergence,
    select_final_candidate,
)
from pilot_2026.adaptive_design.hybrid_optimizer import (  # noqa: E402
    CheckpointedSwarmState,
    advance_checkpointed_swarm,
    initialize_checkpointed_swarm,
)
from pilot_2026.adaptive_design.build_wave1_execution_package import (  # noqa: E402
    _controller_blocks,
)
from pilot_2026.adaptive_design.run_wave1_final_search import (  # noqa: E402
    _common_completed_continuation_windows,
    _history_from_checkpoints,
)
from pilot_2026.adaptive_design.optimize_wave1_sampling_and_plots import (  # noqa: E402
    FIGURE_NAMES,
    _audit_generated_figures,
    _capture_intervals,
    _finish_figure,
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
        self.assertAlmostEqual(sum(amount for _, amount in schedule), 80.0)
        self.assertTrue(
            any(row["action"] == "relocate_duplicate_nutrition_time" for row in decoded.repair_log)
        )
        self.assertTrue(
            any(row["action"] == "scale_total_yan_to_limit" for row in decoded.repair_log)
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

    def test_sampling_optimizer_enforces_95pct_paired_tail_guardrail(self) -> None:
        preliminary = np.ones(64)
        optimized = np.ones(64)
        optimized[:4] = 0.99
        optimized[4:] = 1.01
        accepted, reason = optimized_schedule_is_acceptable(
            optimized, preliminary, self.config
        )
        self.assertFalse(accepted)
        self.assertEqual(
            reason, "optimized_fraction_exceeding_preliminary_below_95pct"
        )

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

    def test_one_and_five_degree_temperature_jumps_are_valid(self) -> None:
        for jump in (1.0, 5.0):
            policy = DesignPolicy(
                f"jump_{jump}", tuple([18.0] * 7 + [18.0 + jump] * 7), tuple()
            )
            metrics = temperature_profile_metrics(policy, self.config)
            self.assertEqual(metrics["temperature_changes"], 1)
            self.assertAlmostEqual(float(metrics["maximum_temperature_jump_c"]), jump)

    def test_oversize_temperature_jump_is_explicitly_reconstructed(self) -> None:
        bounds = vector_bounds(self.config)
        raw = bounds[:, 0].copy()
        raw[0] = 18.0
        changes = int(self.config["future_process"]["maximum_temperature_changes"])
        raw[1] = 1.0
        raw[1 + changes] = 7.0
        raw[1 + 2 * changes] = 27.0
        decoded = decode_policy_vector("oversize", raw, self.config)
        metrics = temperature_profile_metrics(decoded.policy, self.config)
        self.assertLessEqual(float(metrics["maximum_temperature_jump_c"]), 5.0)
        self.assertTrue(
            any(
                row["action"] == "reconstruct_oversize_temperature_jump"
                for row in decoded.repair_log
            )
        )
        self.assertTrue(all(15.0 <= value <= 27.0 for value in decoded.policy.temperature_c))

    def test_product_derived_total_yan_limit_and_mass_units(self) -> None:
        limits = self.config["nutrition"]["derived_product_limits"]
        self.assertAlmostEqual(float(limits["maximum_total_yan_mg_l"]), 80.0)
        translated = nutrition_product_masses(80.0, self.config)
        self.assertAlmostEqual(float(translated["organic_product_g"]), 92.0)
        self.assertAlmostEqual(float(translated["dap_product_g"]), 46.0)
        self.assertAlmostEqual(float(translated["reconstructed_yan_mg_l"]), 80.0)
        self.assertAlmostEqual(float(translated["dimensional_error_mg_l"]), 0.0)

    def test_nutrition_feasibility_one_80_two_40_and_two_80(self) -> None:
        temperature = tuple([18.0] * 14)
        one_80 = DesignPolicy("one_80", temperature, ((2.0, 80.0),))
        two_40 = DesignPolicy("two_40", temperature, ((2.0, 40.0), (18.0, 40.0)))
        two_80 = DesignPolicy("two_80", temperature, ((2.0, 80.0), (18.0, 80.0)))
        self.assertTrue(bool(nutrition_policy_metrics(one_80, self.config)["feasible"]))
        self.assertTrue(bool(nutrition_policy_metrics(two_40, self.config)["feasible"]))
        self.assertFalse(bool(nutrition_policy_metrics(two_80, self.config)["feasible"]))

    def test_baseline_sampling_slot_is_added_at_campaign_start(self) -> None:
        allowed = allowed_sampling_times(self.config)
        self.assertEqual(float(allowed[0]), 0.0)
        self.assertIn(2.0, allowed.tolist())

    def test_probe_bias_is_feedback_error_not_physical_heat_addition(self) -> None:
        policy = DesignPolicy("controller", tuple([18.0] * 14), tuple())
        trajectory = actuator_temperature_trajectory(
            policy,
            self.config,
            {"probe_bias_c": 0.5, "tracking_error_c": 0.0},
        )
        self.assertAlmostEqual(float(trajectory["physical_temperature_c"][-1]), 17.5, places=6)
        self.assertAlmostEqual(float(trajectory["probe_reading_c"][-1]), 18.0, places=6)

    def test_complete_approved_actuator_envelope_has_243_scenarios(self) -> None:
        scenarios = approved_actuator_scenarios(self.config)
        self.assertEqual(len(scenarios), 9 * 3 * 3 * 3 * 1)
        self.assertEqual({row["command_delay_h"] for row in scenarios}, {0.0})

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

    def test_nutrition_50_50_is_encoded_and_product_grams_are_available(self) -> None:
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
        self.assertAlmostEqual(float(translation.iloc[0]["organic_product_g"]), 92.0)
        self.assertAlmostEqual(float(translation.iloc[0]["dap_fda_g"]), 46.0)
        self.assertFalse(bool(translation.iloc[0]["translation_blocker"]))

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
        self.assertEqual(
            dict(zip(first["logical_profile"], first["tank"])),
            {"anchor": "TK33", "A": "TK31", "B": "TK32"},
        )
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

    def test_campaign_start_and_anchor_pulse_use_approved_local_grid(self) -> None:
        self.assertEqual(
            self.config["future_process"]["start_local"],
            "2026-07-20T15:00:00-04:00",
        )
        self.assertEqual(self.config["anchor"]["nutrition_pulses_mg_yan_l"], [[46.0, 80.0]])

    def test_controller_translation_has_12h_blocks_and_hard_jump_check(self) -> None:
        actions = np.asarray([])  # sentinel keeps this test's construction explicit
        del actions
        import pandas as pd

        frame = pd.DataFrame(
            [
                {"policy": name, "action": "temperature_setpoint", "time_h": 0.0, "value": 18.0}
                for name in ("anchor_18C", "candidate_A", "candidate_B")
            ]
        )
        blocks = _controller_blocks(
            frame, {"anchor": "TK33", "A": "TK31", "B": "TK32"}, self.config
        )
        self.assertEqual(len(blocks), 42)
        self.assertTrue(blocks["end_h"].sub(blocks["start_h"]).eq(12.0).all())
        self.assertTrue(blocks["jump_within_5c"].all())

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

    def test_all_required_final_figure_names_and_exact_watermark_are_declared(self) -> None:
        source = (ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py").read_text(
            encoding="utf-8"
        )
        for name in (
            "candidate_profiles_final.png",
            "executed_temperature_ensemble_final.png",
            "aroma_predictions_final.png",
            "sampling_schedule_final.png",
            "capture_intervals_final.png",
            "information_gain_distribution_final.png",
            "paired_information_comparison_final.png",
            "drying_margin_validation_final.png",
            "actuator_robustness_final.png",
            "operational_summary_final.png",
        ):
            self.assertIn(name, source)
        self.assertIn(
            "COMPUTATIONAL CANDIDATE — NOT AUTHORIZED FOR PHYSICAL EXECUTION", source
        )

    def test_final_search_declares_selection_checkpoint_and_post_fim_outputs(self) -> None:
        source = (ADAPTIVE_DIR / "run_wave1_final_search.py").read_text(encoding="utf-8")
        for name in (
            "per_seed_champions_8_member.csv",
            "per_seed_champions_64_member.csv",
            "cross_seed_policy_distance.csv",
            "practical_convergence_gate.json",
            "checkpoint_manifest.json",
            "eligible_candidate_comparison.csv",
            "selected_candidate_provenance.json",
            "post_search_fim_validation.json",
            "post_search_fim_validation_cases.csv",
            "post_search_fd_stability.csv",
            "post_search_grid_stability.csv",
            "selected_candidate_local_qualified",
            "local_multistart_consistency",
        ):
            self.assertIn(name, source)

    def test_execution_package_keeps_all_release_flags_closed(self) -> None:
        source = (ADAPTIVE_DIR / "build_wave1_execution_package.py").read_text(
            encoding="utf-8"
        )
        for folder in (
            "controller",
            "calendar",
            "nutrition",
            "sampling",
            "capture",
            "tanks",
            "checklists",
            "authorization",
        ):
            self.assertIn(f'"{folder}"', source)
        for declaration in (
            '"profiles_for_physical_execution": False',
            '"physical_execution_authorized": False',
            '"executable_schedule_issued": False',
            '"tank_assignments": []',
        ):
            self.assertIn(declaration, source)

    def test_actuator_uncertainty_contains_all_nine_empirical_taus(self) -> None:
        values = self.config["future_process"]["temperature_actuator"]["empirical_tau_h"]
        self.assertEqual(len(values), 9)
        self.assertAlmostEqual(min(values), 0.21148207178811254)
        self.assertAlmostEqual(max(values), 0.45812979693594796)

    def test_final_selection_does_not_let_worse_refinement_displace_original(self) -> None:
        common = {
            "feasible": True,
            "total_thermal_variation_c": 4.0,
            "temperature_changes": 1,
            "total_yan_mg_l": 80.0,
            "minimum_action_margin_h": 30.0,
        }
        original = {
            **common,
            "candidate_id": "original",
            "source_type": "pso_original",
            "full_objective": -20.0,
            "maximum_temperature_jump_c": 4.0,
        }
        refinement = {
            **common,
            "candidate_id": "refined",
            "source_type": "accepted_refinement",
            "parent_candidate_id": "other",
            "full_objective": -19.0,
            "maximum_temperature_jump_c": 3.0,
        }
        selected, _ = select_final_candidate([original, refinement], self.config)
        self.assertEqual(selected["candidate_id"], "original")

    def test_checkpoint_resume_reproduces_next_iteration_exactly(self) -> None:
        bounds = np.asarray([[-2.0, 2.0], [-3.0, 3.0]])
        objective = lambda values: float(np.sum((values - 0.25) ** 2))
        state = initialize_checkpointed_swarm(
            objective, bounds, particles=6, seed=1234, fidelity="four_member"
        )
        state = advance_checkpointed_swarm(state, objective, bounds)
        restored = CheckpointedSwarmState.from_payload(state.to_payload())
        uninterrupted = advance_checkpointed_swarm(state, objective, bounds)
        resumed = advance_checkpointed_swarm(restored, objective, bounds)
        np.testing.assert_array_equal(resumed.positions, uninterrupted.positions)
        np.testing.assert_array_equal(resumed.velocities, uninterrupted.velocities)
        np.testing.assert_array_equal(
            resumed.personal_best_positions, uninterrupted.personal_best_positions
        )
        self.assertEqual(resumed.rng_state, uninterrupted.rng_state)

    def test_final_checkpoint_manifest_uses_windows_extended_paths(self) -> None:
        source = (ADAPTIVE_DIR / "run_wave1_final_search.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("filesystem_path(path).stat().st_size", source)
        self.assertNotIn('"bytes": path.stat().st_size', source)

    def test_resume_recognizes_only_common_completed_continuation_windows(self) -> None:
        complete = {
            seed: SimpleNamespace(iteration=29) for seed in (1, 2, 3, 4, 5)
        }
        partial = {**complete, 5: SimpleNamespace(iteration=28)}
        self.assertEqual(
            _common_completed_continuation_windows(complete, 26, 3), 1
        )
        self.assertEqual(_common_completed_continuation_windows(partial, 26, 3), 0)

    def test_resume_reconstructs_complete_multifidelity_history(self) -> None:
        checkpoints = (
            (1, "four_member", 0, -1.0),
            (2, "eight_member", 20, -2.0),
            (3, "eight_member", 26, -3.0),
            (4, "eight_member", 29, -3.1),
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            for sequence, fidelity, iteration, value in checkpoints:
                path = root / (
                    f"seed_7_{sequence:04d}_{fidelity}_iteration_{iteration:04d}.json.gz"
                )
                with gzip.open(path, "wt", encoding="utf-8") as stream:
                    json.dump(
                        {
                            "checkpoint_sequence": sequence,
                            "swarm": {
                                "seed": 7,
                                "fidelity": fidelity,
                                "iteration": iteration,
                                "global_best_value": value,
                            },
                        },
                        stream,
                    )
            rows = _history_from_checkpoints(root, 20, 6)
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            [row["stage"] for row in rows],
            [
                "stage1_four_member_exploration",
                "stage2_eight_member_rescore",
                "stage2_eight_member_continuation",
                "stage5_eight_member_extension",
            ],
        )

    def test_practical_convergence_uses_seed_champions_and_policy_family(self) -> None:
        champions = [
            {"independent_seed": seed, "full_ensemble_robust_score": score, "feasible": True}
            for seed, score in ((1, 20.0), (2, 19.95), (3, 19.92), (4, 18.0), (5, 17.0))
        ]
        distances = [
            {"left_seed": left, "right_seed": right, "policy_distance": 0.1}
            for left, right in ((1, 2), (1, 3), (2, 3))
        ]
        local_config = json.loads(json.dumps(self.config))
        local_config["independent_seeds"] = [1, 2, 3, 4, 5]
        result = practical_convergence(champions, distances, 0.0005, local_config)
        self.assertTrue(result["passed"])

    def test_sampling_source_is_explicit_not_latest_glob(self) -> None:
        source = (ADAPTIVE_DIR / "optimize_wave1_sampling_and_plots.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('add_argument("--source-search-run", required=True)', source)
        self.assertNotIn("sorted((RESULT_ROOT", source)

    def test_final_figure_audit_checks_exact_nonblank_inventory(self) -> None:
        import matplotlib.pyplot as plt

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            paths = {name: root / name for name in FIGURE_NAMES}
            x = np.linspace(0.0, 1.0, 250)
            for index, name in enumerate(FIGURE_NAMES):
                fig, axis = plt.subplots(figsize=(7.0, 5.0))
                axis.plot(x, np.sin((index + 1) * np.pi * x), linewidth=2)
                axis.set_title(name)
                _finish_figure(fig, paths[name], "Automated QA fixture")
            chart_map = [{"figure": name} for name in FIGURE_NAMES]
            audit = _audit_generated_figures(paths, chart_map)
        self.assertEqual(audit["verdict"], "PASS")
        self.assertTrue(audit["exact_expected_inventory_in_order"])
        self.assertTrue(all(row["passed"] for row in audit["figures"]))


if __name__ == "__main__":
    unittest.main()
    nutrition_policy_metrics,
    nutrition_product_masses,
