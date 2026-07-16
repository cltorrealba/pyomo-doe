from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
RESULTS = FERMENTATION_DIR / "pilot_2026" / "results" / "data_integration_2026"
CONFIG = FERMENTATION_DIR / "pilot_2026" / "data_integration_config.json"


class Pilot2026DataIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))
        cls.qc = json.loads((RESULTS / "qc_summary.json").read_text(encoding="utf-8"))

    def test_all_expected_audit_totals_pass(self) -> None:
        checks = self.qc["expected_value_checks"]
        self.assertEqual(set(checks), set(self.config["expected_totals"]))
        self.assertTrue(all(item["pass"] for item in checks.values()), checks)

    def test_nine_confirmed_windows(self) -> None:
        with (RESULTS / "process_windows.csv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 9)
        self.assertEqual(len({row["experiment_id"] for row in rows}), 9)
        self.assertTrue(all(row["mapping_status"] == "confirmed" for row in rows))

    def test_confirmed_co2_secondary_counts(self) -> None:
        co2 = self.qc["co2"]
        self.assertEqual(co2["invalid_intermediate_zero_seconds"], 500_958)
        self.assertEqual(co2["valid_observed_edge_zero_seconds"], 72_740)
        self.assertEqual(co2["postprocess_observed_seconds"], 153_446)
        self.assertEqual(co2["artificial_initial_zero_seconds"], 86_709)
        self.assertEqual(co2["artificial_final_zero_seconds"], 9)

    def test_temperature_export_is_lossless_before_second_deduplication(self) -> None:
        temperature = self.qc["temperature"]
        self.assertEqual(temperature["raw_relevant_temperature_rows"], 15_617)
        self.assertEqual(temperature["duplicate_temperature_seconds_removed"], 3)
        self.assertEqual(temperature["normalized_temperature_rows"], 15_614)

    def test_protocol_c_is_correct_but_21c_is_not_claimed_as_executed(self) -> None:
        profile = self.qc["profile_c"]
        self.assertEqual(profile["confirmed_protocol"], "16_to_18_to_21C")
        self.assertEqual(profile["setpoints_inside_calibration_window"], [16.0, 18.0])
        self.assertFalse(profile["reached_21C_inside_calibration_window"])

    def test_gc_dilution_and_scope(self) -> None:
        gc = self.qc["gc"]
        self.assertEqual(gc["gc_mix_samples"], 33)
        self.assertEqual(gc["gc_wine_samples"], 39)
        self.assertEqual(gc["gc_initial_wine_samples_without_mix"], 6)
        self.assertEqual(gc["gc_mix_wine_pairs"], 33)
        self.assertEqual(gc["gc_mix_wine_pair_coverage"], 1.0)
        self.assertEqual(gc["gc_analyte_results"], 432)
        self.assertEqual(gc["gc_mix_analyte_results"], 198)
        self.assertEqual(gc["gc_wine_analyte_results"], 234)
        self.assertEqual(gc["gc_paired_analyte_results"], 198)
        self.assertEqual(gc["gc_numeric_mix_results"], 111)
        self.assertEqual(gc["lot1_gc_mix_samples_in_model_scope"], 0)
        self.assertFalse(self.qc["massview_second_normalization_applied"])

    def test_gc_dates_pairing_dilution_and_interval_mass(self) -> None:
        with (RESULTS / "gc_mix_wine_pairs_qc.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            pairs = list(csv.DictReader(handle))
        self.assertEqual(len(pairs), 33)
        self.assertTrue(all(row["wine_pair_available"] == "True" for row in pairs))
        lot2 = [row for row in pairs if row["experiment_id"] in {"26157", "26158", "26159"}]
        self.assertTrue(all(row["timestamp"].startswith("2026-04-") for row in lot2))
        self.assertTrue(all(float(row["capture_interval_h"]) > 0 for row in pairs))

        with (RESULTS / "gc_results_long_qc.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            results = list(csv.DictReader(handle))
        mix = [row for row in results if row["sample_type"] == "condensate_mix"]
        wine = [row for row in results if row["sample_type"] == "wine_mef"]
        self.assertTrue(all(row["dilution_factor"] == "1000" for row in mix))
        self.assertTrue(all(row["dilution_factor"] == "1" for row in wine))
        censored = [row for row in results if row["result_status"] in {"NQ", "below_loq"}]
        self.assertTrue(all(row["model_observation_type"] == "left_censored" for row in censored))
        self.assertTrue(all(float(row["censoring_upper_bound_ug_l"]) > 0 for row in censored))
        quantified_mix = next(row for row in mix if row["result_status"] == "quantified")
        expected_mass = (
            float(quantified_mix["sample_concentration_ug_l"])
            * float(quantified_mix["total_condensate_ml"])
            / 1000.0
        )
        self.assertTrue(
            math.isclose(float(quantified_mix["captured_mass_ug"]), expected_mass)
        )

    def test_owner_confirmed_ethanol_correction(self) -> None:
        with (RESULTS / "primary_results_qc.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        corrected = [row for row in rows if row["ethanol_correction"] != "none"]
        self.assertEqual(len(corrected), 8)
        self.assertTrue(
            all(float(row["ethanol_real_original_percent_vv"]) == -0.3 for row in corrected)
        )
        self.assertTrue(all(float(row["Cf Alcolyzer Real (% v/v)"]) == 0.0 for row in corrected))
        self.assertEqual(self.qc["primary_results"]["postprocess_rows_excluded"], 6)

    def test_lot1_co2_is_qc_only(self) -> None:
        with (RESULTS / "co2_run_qc_summary.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        lot1 = [row for row in rows if row["experiment_id"] in {"26134", "26135", "26136"}]
        later = [row for row in rows if row["experiment_id"] not in {"26134", "26135", "26136"}]
        self.assertTrue(all(row["co2_model_include"] == "False" for row in lot1))
        self.assertTrue(all(row["co2_model_include"] == "True" for row in later))

    def test_lot3_events_are_reconstructed_from_lot1(self) -> None:
        with (RESULTS / "operational_events_qc.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 45)
        lot3 = [row for row in rows if row["experiment_id"] in {"26210", "26211", "26212"}]
        self.assertEqual(len(lot3), 15)
        self.assertTrue(
            all(row["event_origin"] == "reconstructed_from_lot1_relative_schedule" for row in lot3)
        )
        references = {"26210": "26134", "26211": "26135", "26212": "26136"}
        self.assertTrue(
            all(row["reference_experiment_id"] == references[row["experiment_id"]] for row in lot3)
        )
        nutrient = [row for row in lot3 if row["dose_parse_status"] == "parsed_springferm_organic_plus_fda"]
        self.assertEqual(len(nutrient), 6)
        self.assertTrue(
            all(math.isclose(float(row["yan_added_mg_l"]), 90.4347826087) for row in nutrient)
        )

    def test_manifest_records_code_config_and_environment(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(len(manifest["pipeline"]["sha256"]), 64)
        self.assertEqual(len(manifest["config"]["sha256"]), 64)
        self.assertIn("git_head_at_run", manifest["execution_context"])
        self.assertIn("pandas", manifest["execution_context"]["packages"])


if __name__ == "__main__":
    unittest.main()
