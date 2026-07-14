"""Canonical repository paths for active fermentation workflows.

All active code imports paths from this module.  This prevents a script move from
silently redirecting its data inputs or generated results.
"""

from pathlib import Path


SHARED_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SHARED_DIR.parent
REPOSITORY_DIR = FERMENTATION_MODEL_DIR.parent

DATA_DIR = FERMENTATION_MODEL_DIR / "data"
SHARED_RESULTS_DIR = SHARED_DIR / "results"

LABORATORY_2025_DIR = FERMENTATION_MODEL_DIR / "laboratory_2025"
LABORATORY_2026_DIR = FERMENTATION_MODEL_DIR / "laboratory_2026"
LABORATORY_2026_RESULTS_DIR = LABORATORY_2026_DIR / "results"

PILOT_2025_DIR = FERMENTATION_MODEL_DIR / "pilot_2025"
PILOT_2025_RESULTS_DIR = PILOT_2025_DIR / "results"
PILOT_2026_DIR = FERMENTATION_MODEL_DIR / "pilot_2026"
PILOT_2026_RESULTS_DIR = PILOT_2026_DIR / "results"

LEGACY_DIR = FERMENTATION_MODEL_DIR / "legacy"
LEGACY_DEVELOPMENT_RESULTS_DIR = LEGACY_DIR / "development_2026" / "results"
