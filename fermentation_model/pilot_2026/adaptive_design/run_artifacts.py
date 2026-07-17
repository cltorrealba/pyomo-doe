from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


FERMENTATION_MODEL_DIR = Path(__file__).resolve().parents[2]
REPOSITORY_DIR = FERMENTATION_MODEL_DIR.parent
PACKAGE_NAMES = (
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_payload(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def combined_code_hash(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted((Path(item).resolve() for item in paths), key=str):
        digest.update(path.relative_to(REPOSITORY_DIR).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args],
            cwd=REPOSITORY_DIR,
            text=True,
            encoding="utf-8",
            errors="replace",
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in PACKAGE_NAMES:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def create_immutable_run_directory(
    result_root: Path,
    stage: str,
    config: dict[str, Any],
) -> Path:
    """Create a unique run directory; never reuse or overwrite an old run."""

    # Keep this compact because campaign paths can already approach the legacy
    # Windows MAX_PATH limit.  Existing directories are never reused.
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = f"{timestamp}_{sha256_payload(config)[:6]}"
    run_dir = Path(result_root) / stage / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def relative_or_absolute(path: Path) -> str:
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPOSITORY_DIR).as_posix()
    except ValueError:
        return str(resolved)


def build_manifest(
    *,
    run_dir: Path,
    stage: str,
    config: dict[str, Any],
    sources: dict[str, Path],
    code_paths: Iterable[Path],
    random_seeds: list[int],
    status: str,
    convergence: dict[str, Any],
    gate: dict[str, Any],
    outputs: Iterable[Path],
) -> dict[str, Any]:
    status_text = git_value("status", "--porcelain")
    output_paths = [Path(path) for path in outputs]
    return {
        "schema_version": 1,
        "stage": stage,
        "run_id": Path(run_dir).name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "immutable_run_directory": relative_or_absolute(run_dir),
        "status": status,
        "gate": gate,
        "configuration": {
            "sha256": sha256_payload(config),
            "payload": config,
        },
        "sources": {
            name: {
                "path": relative_or_absolute(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for name, path in sorted(sources.items())
        },
        "code": {
            "sha256": combined_code_hash(code_paths),
            "files": [relative_or_absolute(path) for path in sorted(code_paths, key=str)],
        },
        "git": {
            "commit": git_value("rev-parse", "HEAD"),
            "branch": git_value("branch", "--show-current"),
            "dirty": bool(status_text and status_text != "unknown"),
        },
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "executable": sys.executable,
            "packages": package_versions(),
            "solver": {
                "ipopt_executable": shutil.which("ipopt"),
                "ipopt_available_on_path": shutil.which("ipopt") is not None,
            },
        },
        "random_seeds": [int(seed) for seed in random_seeds],
        "solver_status_and_convergence": convergence,
        "outputs": {
            relative_or_absolute(path): {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in sorted(output_paths, key=str)
            if path.exists() and path.is_file()
        },
    }


def write_json(path: Path, payload: Any) -> None:
    output_path = Path(path)
    if sys.platform == "win32":
        resolved = str(output_path.resolve())
        if not resolved.startswith("\\\\?\\"):
            if resolved.startswith("\\\\"):
                resolved = "\\\\?\\UNC\\" + resolved[2:]
            else:
                resolved = "\\\\?\\" + resolved
        output_path = Path(resolved)
    output_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
