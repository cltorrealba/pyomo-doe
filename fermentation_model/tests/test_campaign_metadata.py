from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
AUDIT_PATH = FERMENTATION_DIR / "tools" / "campaign_audit.py"
SPEC = importlib.util.spec_from_file_location("campaign_audit", AUDIT_PATH)
assert SPEC is not None and SPEC.loader is not None
campaign_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(campaign_audit)


class CampaignMetadataTests(unittest.TestCase):
    def test_campaign_registry_is_structurally_valid(self) -> None:
        result = campaign_audit.audit(check_hashes=False)
        self.assertEqual(result["campaigns"], 4)
        self.assertGreaterEqual(result["workflows"], 1)
        self.assertFalse(result["errors"], result["errors"])

    def test_raw_manifest_matches_snapshot(self) -> None:
        result = campaign_audit.audit(check_hashes=True)
        self.assertFalse(result["errors"], result["errors"])


if __name__ == "__main__":
    unittest.main()
