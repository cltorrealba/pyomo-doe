from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
ADAPTIVE_DIR = FERMENTATION_DIR / "pilot_2026" / "adaptive_design"
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.build_model_dataset import (  # noqa: E402
    AdapterConfig,
    EXCLUDED_CO2_RUNS,
    SOURCE_FILES,
    _co2_observations,
    load_model_dataset,
    write_model_dataset,
)


class Pilot2026ModelAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = AdapterConfig.from_json()
        cls.dataset = load_model_dataset(cls.config)

    def test_adapter_qc_passes_and_carries_nine_run_metadata_rows(self) -> None:
        self.assertEqual(self.dataset.qc["verdict"], "PASS")
        metadata = self.dataset.run_metadata
        self.assertEqual(len(metadata), 9)
        self.assertEqual(metadata["experiment_id"].nunique(), 9)
        self.assertTrue(metadata["campaign"].eq("pilot_2026").all())
        self.assertTrue(metadata["yeast"].eq("Zymaflore X5").all())
        for field in (
            "lot",
            "reactor",
            "storage_history",
            "capture_system",
            "co2_sensor",
            "temperature_sensor",
            "temperature_command",
        ):
            self.assertTrue(metadata[field].notna().all(), field)

    def test_primary_observation_units_operators_and_phases_are_explicit(self) -> None:
        primary = self.dataset.primary_observations
        self.assertTrue(primary["unit"].notna().all())
        self.assertTrue(primary["observation_operator"].notna().all())
        self.assertEqual(set(primary["process_phase"]), {"active_process"})
        self.assertTrue(primary["calibration_include"].all())

    def test_total_and_component_sugars_are_never_double_counted(self) -> None:
        sugar = self.dataset.primary_observations[
            self.dataset.primary_observations["state"].isin(
                {"glucose", "fructose", "total_sugar"}
            )
        ]
        states = sugar.groupby("sample_id")["state"].agg(set)
        double_counted = states.apply(
            lambda names: "total_sugar" in names
            and bool(names.intersection({"glucose", "fructose"}))
        )
        self.assertFalse(double_counted.any())
        self.assertFalse(sugar["state"].eq("total_sugar").any())

    def test_temperature_uses_sonda1_and_keeps_commanded_setpoint_separate(self) -> None:
        temperature = self.dataset.temperature_inputs
        self.assertTrue(
            temperature["executed_temperature_operator"].eq(
                "Sonda1 sensor_1_c"
            ).all()
        )
        self.assertTrue(
            temperature["command_operator"].eq("controller setpoint").all()
        )
        self.assertEqual(set(temperature["process_phase"]), {"active_process"})
        self.assertTrue(temperature["kinetic_include"].all())

    def test_co2_is_reproducible_weighted_and_excludes_lot1(self) -> None:
        co2 = self.dataset.co2_observations
        self.assertFalse(set(co2["experiment_id"]).intersection(EXCLUDED_CO2_RUNS))
        self.assertEqual(set(co2["experiment_id"]), {"26157", "26158", "26159", "26210", "26211", "26212"})
        self.assertTrue((co2["likelihood_weight"] > 0.0).all())
        self.assertTrue((co2["sigma_multiplier"] > 0.0).all())
        for _run, group in co2.groupby("experiment_id"):
            expected = float(group["effective_sample_size_minutes"].iloc[0])
            self.assertTrue(
                np.isclose(group["likelihood_weight"].sum(), expected, rtol=1e-12)
            )

    def test_co2_time_origin_is_the_authoritative_sampling_window_start(self) -> None:
        raw = pd.read_csv(
            SOURCE_FILES["co2_minute"], parse_dates=["minute"], low_memory=False
        )
        raw["experiment_id"] = raw["experiment_id"].astype(str)
        starts = raw.groupby("experiment_id")["minute"].min()
        for run, group in self.dataset.co2_observations.groupby("experiment_id"):
            expected = (
                pd.to_datetime(group["timestamp"]) - pd.Timestamp(starts.loc[run])
            ).dt.total_seconds() / 3600.0
            # Aggregated time is the mean minute in each bin, so it must remain
            # within one configured bin of the first timestamp in that bin.
            difference = group["time_h"].to_numpy(dtype=float) - expected.to_numpy(
                dtype=float
            )
            self.assertTrue((difference >= 0.0).all())
            self.assertTrue(
                (difference <= self.config.co2_bin_minutes / 60.0 + 1e-12).all()
            )

    def test_excluded_co2_reentry_fails_closed(self) -> None:
        raw = pd.read_csv(SOURCE_FILES["co2_minute"], nrows=20, parse_dates=["minute"])
        raw["experiment_id"] = "26134"
        raw["co2_model_include"] = True
        raw["mean_flow_ln_min_calibration"] = 1.0
        with self.assertRaisesRegex(ValueError, "re-entered"):
            _co2_observations(raw, self.config)

    def test_censored_wine_and_condensate_observations_have_upper_bounds(self) -> None:
        wine = self.dataset.wine_aroma_observations
        wine_censored = wine["model_observation_type"].eq("left_censored")
        self.assertGreater(int(wine_censored.sum()), 0)
        self.assertTrue(wine.loc[wine_censored, "observed_value"].isna().all())
        self.assertTrue(wine.loc[wine_censored, "upper_bound"].gt(0.0).all())

        condensate = self.dataset.condensate_interval_observations
        mix_censored = condensate["mix_model_observation_type"].eq("left_censored")
        self.assertGreater(int(mix_censored.sum()), 0)
        self.assertTrue(condensate.loc[mix_censored, "observed_value"].isna().all())
        self.assertTrue(condensate.loc[mix_censored, "upper_bound"].gt(0.0).all())
        self.assertTrue((condensate["capture_interval_h"] > 0.0).all())
        self.assertTrue(
            condensate["observation_operator"].str.startswith("integral(").all()
        )

    def test_density_triggered_pulses_have_data_derived_timing_intervals(self) -> None:
        events = self.dataset.operational_events
        self.assertTrue(events["calibration_include"].all())
        self.assertEqual(len(events), 44)
        uncertain = events[events["timing_uncertain"]]
        self.assertEqual(len(uncertain), 9)
        self.assertTrue(uncertain["timing_interval_start"].notna().all())
        self.assertTrue(uncertain["timing_interval_end"].notna().all())
        self.assertTrue(
            (
                pd.to_datetime(uncertain["timing_interval_end"])
                >= pd.to_datetime(uncertain["timing_interval_start"])
            ).all()
        )
        self.assertTrue(
            uncertain["timing_uncertainty_basis"].str.contains(
                "density|latent_event_time", regex=True
            ).all()
        )

    def test_carbon_balance_is_diagnostic_and_never_claimed_closed(self) -> None:
        carbon = self.dataset.carbon_balance_diagnostic
        self.assertEqual(len(carbon), 9)
        self.assertTrue(carbon["limitations"].str.contains("biomass carbon").all())
        lot1 = carbon[carbon["experiment_id"].isin(EXCLUDED_CO2_RUNS)]
        self.assertTrue(lot1["diagnostic_status"].eq("not_evaluable").all())
        later = carbon[~carbon["experiment_id"].isin(EXCLUDED_CO2_RUNS)]
        self.assertTrue(
            later["diagnostic_status"].eq("diagnostic_only_incomplete").all()
        )

    def test_written_run_has_source_code_environment_and_seed_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = write_model_dataset(
                self.dataset,
                config=self.config,
                config_path=ADAPTIVE_DIR / "model_dataset_config.json",
                result_root=Path(temporary),
            )
            manifest = json.loads(
                (run_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["gate"]["verdict"], "PASS")
            self.assertEqual(manifest["random_seeds"], [20260716])
            self.assertEqual(len(manifest["code"]["sha256"]), 64)
            self.assertTrue(manifest["sources"])
            self.assertIn("packages", manifest["environment"])
            self.assertFalse(manifest["gate"]["profiles_for_physical_execution"])


if __name__ == "__main__":
    unittest.main()
