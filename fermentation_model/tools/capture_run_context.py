from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
REPO_DIR = FERMENTATION_DIR.parent
RAW_MANIFEST = FERMENTATION_DIR / "campaigns" / "raw_data_manifest.csv"
PACKAGES = (
    "idaes-pse",
    "pyomo",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "openpyxl",
    "nbformat",
    "nbclient",
    "thermo",
    "chemicals",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=REPO_DIR, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def package_versions() -> dict[str, str | None]:
    versions = {}
    for package in PACKAGES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def main() -> None:
    parser = argparse.ArgumentParser(description="Capture a reproducible run context")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    config = args.config.resolve()
    if not config.exists():
        raise SystemExit(f"Config does not exist: {config}")

    status = git_value("status", "--porcelain")
    ipopt_path = shutil.which("ipopt")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git": {
            "commit": git_value("rev-parse", "HEAD"),
            "branch": git_value("branch", "--show-current"),
            "dirty": bool(status and status != "unknown"),
        },
        "platform": platform.platform(),
        "python": sys.version,
        "config": {
            "path": config.relative_to(REPO_DIR).as_posix()
            if config.is_relative_to(REPO_DIR)
            else str(config),
            "sha256": sha256(config),
        },
        "packages": package_versions(),
        "solver": {
            "ipopt_executable": ipopt_path,
            "ipopt_available_on_path": ipopt_path is not None,
        },
        "raw_data_manifest_sha256": sha256(RAW_MANIFEST)
        if RAW_MANIFEST.exists()
        else None,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
