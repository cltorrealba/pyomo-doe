from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
AUDIT_SCRIPT = FERMENTATION_DIR / "tools" / "campaign_audit.py"
HISTORICAL_ROOT = PILOT_DIR / "results" / "adaptive_design_2026" / "baseline_calibration"
RAW_MANIFEST = FERMENTATION_DIR / "campaigns" / "raw_data_manifest.csv"
BASE_COMMIT = "d3795f79e1fc9830cbc5191c31f4291d6d2f56f6"


def _run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=REPOSITORY_DIR,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _historical_audit() -> tuple[Path, dict[str, Any]]:
    candidates = sorted(HISTORICAL_ROOT.glob("*/repository_audit.json"))
    for path in candidates:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if len(payload.get("stdout", {}).get("errors", [])) == 53:
            return path, payload
    raise FileNotFoundError("No historical 53-error audit artifact was found")


def build_resolution() -> dict[str, Any]:
    historical_path, historical = _historical_audit()
    historical_errors = historical["stdout"]["errors"]
    raw_size = [item for item in historical_errors if item.startswith("raw file size changed:")]
    structural = [
        item
        for item in historical_errors
        if item == (
            "fermentation_model/results is ambiguous; use shared/results or a "
            "campaign-owned results directory"
        )
    ]
    other = sorted(set(historical_errors) - set(raw_size) - set(structural))

    current = _run([sys.executable, str(AUDIT_SCRIPT), "--check-hashes"])
    try:
        current_payload = json.loads(current.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError("Current repository audit did not emit JSON") from error

    raw_diff = _run(
        [
            "git",
            "diff",
            "--name-only",
            f"{BASE_COMMIT}..HEAD",
            "--",
            "fermentation_model/data",
            "fermentation_model/campaigns/raw_data_manifest.csv",
        ]
    )
    raw_worktree = _run(
        [
            "git",
            "status",
            "--porcelain",
            "--",
            "fermentation_model/data",
            "fermentation_model/campaigns/raw_data_manifest.csv",
        ]
    )
    base_manifest = _run(
        [
            "git",
            "show",
            f"{BASE_COMMIT}:fermentation_model/campaigns/raw_data_manifest.csv",
        ]
    )
    current_manifest = RAW_MANIFEST.read_bytes()
    legacy_root_exists = (FERMENTATION_DIR / "results").exists()

    checks = {
        "historical_count_is_53": len(historical_errors) == 53,
        "historical_partition_is_52_raw_plus_1_structural": (
            len(raw_size) == 52 and len(structural) == 1 and not other
        ),
        "current_hash_audit_has_zero_errors": (
            current.returncode == 0 and not current_payload.get("errors")
        ),
        "raw_tree_unchanged_since_base_commit": not raw_diff.stdout.strip(),
        "raw_tree_clean_in_worktree": not raw_worktree.stdout.strip(),
        "raw_manifest_byte_identical_to_base": (
            base_manifest.returncode == 0
            and base_manifest.stdout.encode("utf-8") == current_manifest
        ),
        "ambiguous_legacy_result_root_absent": not legacy_root_exists,
    }
    verdict = "PASS" if all(checks.values()) else "FAIL"
    return {
        "verdict": verdict,
        "base_commit": BASE_COMMIT,
        "historical_artifact": historical_path.relative_to(REPOSITORY_DIR).as_posix(),
        "historical_error_count": len(historical_errors),
        "classification": {
            "raw_worktree_size_drift": len(raw_size),
            "ambiguous_generated_result_root": len(structural),
            "unclassified": other,
        },
        "resolution": {
            "raw_worktree_size_drift": (
                "Resolved by reproducing from the committed clean tree and checking every "
                "raw byte against the existing manifest. No raw file or manifest entry was "
                "changed or rebaselined."
            ),
            "ambiguous_generated_result_root": (
                "Resolved by using campaign-owned result directories; the legacy aggregate "
                "directory is absent."
            ),
        },
        "checks": checks,
        "current_audit": current_payload,
        "provenance": {
            "raw_manifest_sha256": _sha256_bytes(current_manifest),
            "base_raw_manifest_sha256": _sha256_bytes(
                base_manifest.stdout.encode("utf-8")
            ),
            "raw_diff_paths": [line for line in raw_diff.stdout.splitlines() if line],
            "raw_worktree_paths": [
                line for line in raw_worktree.stdout.splitlines() if line
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Classify and verify resolution of the historical 53-error audit"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = build_resolution()
    encoded = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    if payload["verdict"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
