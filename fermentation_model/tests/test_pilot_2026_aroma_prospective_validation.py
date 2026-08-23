from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.validate_aroma_prospective_package import (  # noqa: E402
    DEFAULT_CONTRACT,
    ROOT_DIR,
    _sha256,
    validate_contract,
    validate_package,
)


SPECIES = ("ethyl_octanoate", "isoamyl_acetate")
RUNS = ("PV01", "PV02", "PV03")


def _write_valid_package(directory: Path) -> None:
    process_rows = []
    event_rows = []
    liquid_rows = []
    gas_rows = []
    condensate_rows = []
    for run_index, run in enumerate(RUNS):
        for time_h in (0.0, 1.0):
            process_rows.append(
                {
                    "experiment_id": run,
                    "timestamp_utc": f"2026-09-0{run_index + 1}T0{int(time_h)}:00:00Z",
                    "process_time_h": time_h,
                    "temperature_c": 18.0,
                    "setpoint_c": 18.0,
                    "rco2_g_l_h": 0.1,
                    "ethanol_g_l": 10.0 * time_h,
                    "total_sugar_g_l": 200.0 - time_h,
                    "liquid_volume_l": 230.0,
                }
            )
        event_rows.append(
            {
                "experiment_id": run,
                "event_id": f"{run}-N1",
                "event_type": "nutrient_pulse",
                "event_time_h": 0.0,
                "nitrogen_product": "validation_product",
                "product_dose_g_l": 0.1,
                "yan_equivalent_mg_l": 20.0,
                "new_setpoint_c": "",
            }
        )
        for species in SPECIES:
            for replicate in ("R1", "R2"):
                liquid_rows.append(
                    {
                        "experiment_id": run,
                        "sample_id": f"{run}-L01",
                        "process_time_h": 1.0,
                        "species": species,
                        "replicate_id": replicate,
                        "concentration_ug_l": 100.0,
                        "lod_ug_l": 1.0,
                        "loq_ug_l": 3.0,
                        "dilution_factor": 1.0,
                        "qa_status": "quantified",
                    }
                )
                gas_rows.append(
                    {
                        "experiment_id": run,
                        "sample_id": f"{run}-G01",
                        "interval_start_h": 0.0,
                        "interval_end_h": 1.0,
                        "species": species,
                        "replicate_id": replicate,
                        "concentration_ug_l_gas": 2.0,
                        "gas_flow_l_h": 10.0,
                        "gas_temperature_c": 20.0,
                        "lod_ug_l_gas": 0.1,
                        "loq_ug_l_gas": 0.3,
                        "qa_status": "quantified",
                    }
                )
                condensate_rows.append(
                    {
                        "experiment_id": run,
                        "mix_id": f"{run}-MIX-01",
                        "interval_start_h": 0.0,
                        "interval_end_h": 1.0,
                        "species": species,
                        "replicate_id": replicate,
                        "concentration_ug_l": 100.0,
                        "condensate_volume_l": 0.01,
                        "condensate_ethanol_fraction_v_v": 0.1,
                        "dilution_factor": 1.0,
                        "captured_mass_ug": 1.0,
                        "lod_ug_l": 1.0,
                        "loq_ug_l": 3.0,
                        "qa_status": "quantified",
                    }
                )
    standard_rows = []
    for session in ("before_campaign", "mid_campaign", "after_campaign"):
        for species in SPECIES:
            for level in ("low", "medium", "high"):
                for replicate in ("R1", "R2", "R3"):
                    standard_rows.append(
                        {
                            "session_id": session,
                            "standard_id": f"{session}-{species}-{level}",
                            "species": species,
                            "replicate_id": replicate,
                            "standard_level": level,
                            "input_mass_ug": 100.0,
                            "recovered_mass_ug": 50.0,
                            "gas_flow_l_h": 10.0,
                            "gas_temperature_c": 20.0,
                            "carrier_ethanol_fraction_v_v": 0.1,
                            "recovery_fraction": 0.5,
                            "qa_status": "accepted",
                        }
                    )
    for filename, rows in (
        ("process.csv", process_rows),
        ("events.csv", event_rows),
        ("liquid_aroma.csv", liquid_rows),
        ("outlet_gas.csv", gas_rows),
        ("condensate.csv", condensate_rows),
        ("trap_standards.csv", standard_rows),
    ):
        pd.DataFrame(rows).to_csv(directory / filename, index=False)


class ProspectiveAromaPackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(DEFAULT_CONTRACT.read_text(encoding="utf-8"))

    def test_locked_contract_matches_current_model_artifacts(self) -> None:
        self.assertEqual(validate_contract(self.contract, ROOT_DIR), [])
        self.assertFalse(self.contract["locked_model"]["refitting_permitted"])

    def test_locked_text_hash_is_independent_of_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lf = root / "lf.csv"
            crlf = root / "crlf.csv"
            lf.write_bytes(b"a,b\n1,2\n")
            crlf.write_bytes(b"a,b\r\n1,2\r\n")
            self.assertEqual(_sha256(lf), _sha256(crlf))

    def test_missing_tables_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            report = validate_package(Path(directory))
        self.assertFalse(report["valid"])
        self.assertEqual(report["error_count"], 6)

    def test_complete_minimal_package_passes_schema_and_mass_checks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            _write_valid_package(package)
            report = validate_package(package)
        self.assertTrue(report["valid"], report["errors"])
        self.assertAlmostEqual(
            report["diagnostics"][
                "maximum_condensate_mass_recalculation_relative_error"
            ],
            0.0,
        )


if __name__ == "__main__":
    unittest.main()
