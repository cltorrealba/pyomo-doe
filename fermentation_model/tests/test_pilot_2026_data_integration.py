from __future__ import annotations

import csv
import json
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
        self.assertEqual(gc["gc_analyte_results"], 198)
        self.assertEqual(gc["gc_numeric_results"], 111)
        self.assertEqual(gc["lot1_gc_mix_samples_in_model_scope"], 0)
        self.assertFalse(self.qc["massview_second_normalization_applied"])


if __name__ == "__main__":
    unittest.main()
