from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPO_DIR = FERMENTATION_DIR.parent
CAMPAIGN_DIR = FERMENTATION_DIR / "campaigns"
CAMPAIGN_CSV = CAMPAIGN_DIR / "campaigns.csv"
EXPERIMENT_CSV = CAMPAIGN_DIR / "experiments.csv"
WORKFLOW_CSV = CAMPAIGN_DIR / "workflows.csv"
RAW_MANIFEST_CSV = CAMPAIGN_DIR / "raw_data_manifest.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def raw_files(campaigns: list[dict[str, str]]) -> list[tuple[str, Path]]:
    files: list[tuple[str, Path]] = []
    for campaign in campaigns:
        root = REPO_DIR / campaign["raw_data_root"]
        if root.is_file():
            files.append((campaign["campaign_id"], root))
        elif root.is_dir():
            files.extend(
                (campaign["campaign_id"], path)
                for path in sorted(root.rglob("*"))
                if path.is_file()
            )
    return files


def write_raw_manifest(campaigns: list[dict[str, str]]) -> None:
    rows = []
    for campaign_id, path in raw_files(campaigns):
        rows.append(
            {
                "campaign_id": campaign_id,
                "path": path.relative_to(REPO_DIR).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    with RAW_MANIFEST_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["campaign_id", "path", "size_bytes", "sha256"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {RAW_MANIFEST_CSV.relative_to(REPO_DIR)} ({len(rows)} files)")


def audit(check_hashes: bool = False) -> dict[str, object]:
    campaigns = read_csv(CAMPAIGN_CSV)
    experiments = read_csv(EXPERIMENT_CSV)
    workflows = read_csv(WORKFLOW_CSV)
    errors: list[str] = []
    warnings: list[str] = []

    campaign_ids = [row["campaign_id"] for row in campaigns]
    if len(campaign_ids) != len(set(campaign_ids)):
        errors.append("campaigns.csv contains duplicate campaign_id values")

    experiment_keys = [
        (row["campaign_id"], row["experiment_id"]) for row in experiments
    ]
    if len(experiment_keys) != len(set(experiment_keys)):
        errors.append("experiments.csv contains duplicate campaign/experiment keys")

    unknown_campaigns = sorted(
        {row["campaign_id"] for row in experiments} - set(campaign_ids)
    )
    if unknown_campaigns:
        errors.append(f"experiments reference unknown campaigns: {unknown_campaigns}")

    unknown_workflow_campaigns = sorted(
        {row["campaign_id"] for row in workflows} - set(campaign_ids) - {"shared"}
    )
    if unknown_workflow_campaigns:
        errors.append(
            f"workflows reference unknown campaigns: {unknown_workflow_campaigns}"
        )
    for workflow in workflows:
        path = REPO_DIR / workflow["path"]
        if not path.exists():
            errors.append(f"missing workflow: {workflow['workflow_id']}: {workflow['path']}")
        output = workflow.get("authoritative_output", "").strip()
        if output and not (REPO_DIR / output).exists():
            errors.append(
                f"missing workflow output: {workflow['workflow_id']}: {output}"
            )

    campaign_file_counts: dict[str, int] = {}
    for campaign in campaigns:
        campaign_id = campaign["campaign_id"]
        raw_root = REPO_DIR / campaign["raw_data_root"]
        workspace = REPO_DIR / campaign["workspace"]
        if not raw_root.exists():
            errors.append(f"missing raw_data_root for {campaign_id}: {raw_root}")
        if not workspace.exists():
            errors.append(f"missing workspace for {campaign_id}: {workspace}")
        result = campaign.get("authoritative_results", "").strip()
        if result and not (REPO_DIR / result).exists():
            errors.append(f"missing authoritative result for {campaign_id}: {result}")
        campaign_file_counts[campaign_id] = sum(
            1 for cid, _ in raw_files([campaign]) if cid == campaign_id
        )

    missing_experiment_paths = []
    for row in experiments:
        path = REPO_DIR / row["raw_location"]
        if not path.exists() and row["execution_status"] not in {
            "planned_or_in_progress"
        }:
            missing_experiment_paths.append(
                f"{row['campaign_id']}/{row['experiment_id']}: {row['raw_location']}"
            )
    if missing_experiment_paths:
        errors.extend(f"missing experiment path: {item}" for item in missing_experiment_paths)

    if check_hashes:
        if not RAW_MANIFEST_CSV.exists():
            errors.append("raw_data_manifest.csv does not exist")
        else:
            for row in read_csv(RAW_MANIFEST_CSV):
                path = REPO_DIR / row["path"]
                if not path.exists():
                    errors.append(f"manifested raw file is missing: {row['path']}")
                    continue
                if str(path.stat().st_size) != row["size_bytes"]:
                    errors.append(f"raw file size changed: {row['path']}")
                    continue
                if sha256(path) != row["sha256"]:
                    errors.append(f"raw file hash changed: {row['path']}")

    if not (FERMENTATION_DIR / "data" / "Laboratorio 2026" / "DOE_Lote_3").exists():
        warnings.append("Laboratory 2026 Lot 3 raw folder is not present")
    warnings.append(
        "Pilot 2026 still uses a composite placeholder; individual run mapping is pending"
    )

    return {
        "campaigns": len(campaigns),
        "experiments": len(experiments),
        "workflows": len(workflows),
        "raw_file_counts": campaign_file_counts,
        "errors": errors,
        "warnings": warnings,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit fermentation campaign metadata")
    parser.add_argument("--write-raw-manifest", action="store_true")
    parser.add_argument("--check-hashes", action="store_true")
    args = parser.parse_args()

    campaigns = read_csv(CAMPAIGN_CSV)
    if args.write_raw_manifest:
        write_raw_manifest(campaigns)

    result = audit(check_hashes=args.check_hashes)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if result["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
