from __future__ import annotations

import sys
import unittest
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from shared import aroma_partition_unifac  # noqa: E402
from shared import run_secondary_joint_campaign_doe  # noqa: E402


class AromaThermodynamicsContractTests(unittest.TestCase):
    def test_trap_efficiency_has_one_authoritative_definition(self) -> None:
        self.assertIs(
            run_secondary_joint_campaign_doe.TRAP_EFFICIENCY,
            aroma_partition_unifac.TRAP_EFFICIENCY,
        )
        self.assertEqual(
            set(run_secondary_joint_campaign_doe.AROMA_SPECIES),
            set(aroma_partition_unifac.TRAP_EFFICIENCY),
        )

    def test_unifac_dependency_failure_is_explicit(self) -> None:
        self.assertTrue(issubclass(aroma_partition_unifac.UNIFACUnavailable, RuntimeError))
        self.assertFalse(hasattr(aroma_partition_unifac, "FALLBACK_PARTITION_COEFFICIENT"))


if __name__ == "__main__":
    unittest.main()
