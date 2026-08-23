from __future__ import annotations

"""Validate a sealed prospective aroma campaign before model evaluation."""

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_DIR = SCRIPT_DIR.parent
ROOT_DIR = FERMENTATION_DIR.parent
DEFAULT_CONTRACT = (
    SCRIPT_DIR / "adaptive_design" / "aroma_prospective_validation_contract.json"
)
TEXT_HASH_SUFFIXES = {".csv", ".json", ".md", ".py", ".txt", ".yaml", ".yml"}


def _sha256(path: Path) -> str:
    """Hash text artifacts independently of the checkout line-ending policy."""

    digest = hashlib.sha256()
    if path.suffix.lower() in TEXT_HASH_SUFFIXES:
        content = path.read_text(encoding="utf-8")
        canonical = content.replace("\r\n", "\n").replace("\r", "\n")
        digest.update(canonical.encode("utf-8"))
        return digest.hexdigest()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_contract(contract: dict[str, Any], repository_dir: Path) -> list[str]:
    errors: list[str] = []
    locked = contract["locked_model"]
    for path_key, hash_key in (
        ("parameter_file", "parameter_file_sha256"),
        ("fit_config", "fit_config_sha256"),
    ):
        path = repository_dir / locked[path_key]
        if not path.exists():
            errors.append(f"Locked source is missing: {path}")
        elif _sha256(path).lower() != str(locked[hash_key]).lower():
            errors.append(f"Locked source hash changed: {path}")
    if locked.get("refitting_permitted") is not False:
        errors.append("Prospective contract must prohibit refitting")
    required_species = set(contract["campaign"]["required_species"])
    if set(locked["parameters"]) != required_species:
        errors.append("Locked parameter species differ from campaign species")
    for filename, schema in contract["tables"].items():
        columns = set(schema["required_columns"])
        missing_key_columns = set(schema["primary_key"]) - columns
        if missing_key_columns:
            errors.append(
                f"{filename}: primary-key columns absent from schema: {sorted(missing_key_columns)}"
            )
    return errors


def _allowed_values(specification: str) -> set[str] | None:
    if "|" not in specification or specification.startswith("float"):
        return None
    if "string" in specification:
        return None
    return set(specification.split(" ", 1)[0].split("|"))


def _validate_table(
    filename: str,
    schema: dict[str, Any],
    package_dir: Path,
) -> tuple[pd.DataFrame | None, list[str]]:
    errors: list[str] = []
    path = package_dir / filename
    if not path.exists():
        return None, [f"Required table is missing: {filename}"]
    try:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception as exc:  # pragma: no cover - parser supplies diagnostic
        return None, [f"{filename}: cannot read CSV: {exc}"]
    required = list(schema["required_columns"])
    missing = [column for column in required if column not in frame.columns]
    if missing:
        return frame, [f"{filename}: missing columns {missing}"]
    if frame.empty:
        errors.append(f"{filename}: table is empty")
        return frame, errors
    duplicated = frame.duplicated(schema["primary_key"], keep=False)
    if duplicated.any():
        errors.append(
            f"{filename}: {int(duplicated.sum())} rows duplicate the primary key"
        )
    for column, specification in schema["required_columns"].items():
        values = frame[column].astype(str).str.strip()
        optional = "or blank" in specification
        if not optional and values.eq("").any():
            errors.append(f"{filename}.{column}: blank values are not allowed")
        if specification.startswith("float"):
            parsed = pd.to_numeric(values.mask(values.eq("")), errors="coerce")
            invalid = values.ne("") & parsed.isna()
            if invalid.any():
                errors.append(
                    f"{filename}.{column}: {int(invalid.sum())} non-numeric values"
                )
            frame[column] = parsed
        elif specification.startswith("ISO-8601"):
            parsed = pd.to_datetime(values, utc=True, errors="coerce")
            if parsed.isna().any():
                errors.append(
                    f"{filename}.{column}: {int(parsed.isna().sum())} invalid timestamps"
                )
        allowed = _allowed_values(specification)
        if allowed is not None:
            invalid_values = sorted(set(values[values.ne("")]) - allowed)
            if invalid_values:
                errors.append(
                    f"{filename}.{column}: invalid values {invalid_values}; allowed {sorted(allowed)}"
                )
    return frame, errors


def _relative_error(observed: pd.Series, expected: pd.Series) -> pd.Series:
    denominator = np.maximum(np.abs(expected.to_numpy(dtype=float)), 1e-12)
    return pd.Series(
        np.abs(observed.to_numpy(dtype=float) - expected.to_numpy(dtype=float))
        / denominator,
        index=observed.index,
    )


def _domain_checks(
    tables: dict[str, pd.DataFrame], contract: dict[str, Any]
) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    diagnostics: dict[str, Any] = {}
    required_species = set(contract["campaign"]["required_species"])
    minimum_reactors = int(contract["campaign"]["minimum_independent_reactors"])
    experiment_sets = {
        name: set(frame["experiment_id"].astype(str))
        for name, frame in tables.items()
        if "experiment_id" in frame
    }
    all_experiments = set.union(*experiment_sets.values()) if experiment_sets else set()
    if len(all_experiments) < minimum_reactors:
        errors.append(
            f"Campaign has {len(all_experiments)} reactors; requires {minimum_reactors}"
        )
    for name, experiments in experiment_sets.items():
        missing = all_experiments - experiments
        if missing and name in {"process.csv", "events.csv", "liquid_aroma.csv", "outlet_gas.csv", "condensate.csv"}:
            errors.append(f"{name}: missing reactors {sorted(missing)}")

    for name in ("liquid_aroma.csv", "outlet_gas.csv", "condensate.csv"):
        species = set(tables[name]["species"].astype(str))
        missing = required_species - species
        if missing:
            errors.append(f"{name}: missing species {sorted(missing)}")

    process = tables["process.csv"]
    numeric_nonnegative = [
        "process_time_h",
        "rco2_g_l_h",
        "ethanol_g_l",
        "total_sugar_g_l",
        "liquid_volume_l",
    ]
    for column in numeric_nonnegative:
        if (process[column] < 0.0).any():
            errors.append(f"process.csv.{column}: negative values found")
    for run, group in process.groupby("experiment_id"):
        if not group["process_time_h"].is_monotonic_increasing:
            errors.append(f"process.csv: time is not increasing for {run}")

    condensate = tables["condensate.csv"]
    quantified_condensate = condensate[condensate["qa_status"].eq("quantified")]
    expected_mass = (
        quantified_condensate["concentration_ug_l"]
        * quantified_condensate["condensate_volume_l"]
        * quantified_condensate["dilution_factor"]
    )
    mass_error = _relative_error(
        quantified_condensate["captured_mass_ug"], expected_mass
    )
    mass_tolerance = float(
        contract["quality_gates"]["captured_mass_recalculation_relative_tolerance"]
    )
    if (mass_error > mass_tolerance).any():
        errors.append(
            f"condensate.csv: {int((mass_error > mass_tolerance).sum())} captured masses fail recalculation"
        )
    diagnostics["maximum_condensate_mass_recalculation_relative_error"] = (
        float(mass_error.max()) if len(mass_error) else math.nan
    )
    if ((condensate["condensate_ethanol_fraction_v_v"] < 0.0) | (condensate["condensate_ethanol_fraction_v_v"] > 1.0)).any():
        errors.append("condensate.csv: ethanol fractions must be in [0,1]")
    if (condensate["interval_end_h"] <= condensate["interval_start_h"]).any():
        errors.append("condensate.csv: non-positive capture intervals found")

    standards = tables["trap_standards.csv"]
    accepted_standards = standards[standards["qa_status"].eq("accepted")]
    expected_recovery = (
        accepted_standards["recovered_mass_ug"]
        / accepted_standards["input_mass_ug"]
    )
    recovery_error = np.abs(
        accepted_standards["recovery_fraction"] - expected_recovery
    )
    recovery_tolerance = float(
        contract["quality_gates"]["trap_recovery_recalculation_absolute_tolerance"]
    )
    if (recovery_error > recovery_tolerance).any():
        errors.append("trap_standards.csv: recovery fractions fail recalculation")
    recovery_cv = (
        accepted_standards.groupby(["session_id", "species", "standard_level"])[
            "recovery_fraction"
        ]
        .agg(lambda values: float(values.std(ddof=1) / values.mean()) if len(values) > 1 and values.mean() > 0 else math.nan)
        .dropna()
    )
    maximum_cv = float(recovery_cv.max()) if len(recovery_cv) else math.nan
    diagnostics["maximum_trap_recovery_cv"] = maximum_cv
    if len(recovery_cv) and maximum_cv > float(
        contract["quality_gates"]["maximum_trap_recovery_cv_within_level"]
    ):
        errors.append(f"trap_standards.csv: maximum recovery CV is {maximum_cv:.3f}")

    minimum_replicates = int(contract["sampling"]["minimum_analytical_replicates"])
    replicate_keys = {
        "liquid_aroma.csv": ["experiment_id", "sample_id", "species"],
        "outlet_gas.csv": ["experiment_id", "sample_id", "species"],
        "condensate.csv": ["experiment_id", "mix_id", "species"],
    }
    for name, keys in replicate_keys.items():
        accepted = tables[name][~tables[name]["qa_status"].eq("excluded")]
        counts = accepted.groupby(keys)["replicate_id"].nunique()
        if (counts < minimum_replicates).any():
            errors.append(
                f"{name}: {int((counts < minimum_replicates).sum())} samples have fewer than {minimum_replicates} replicates"
            )
    return errors, diagnostics


def validate_package(
    package_dir: Path, contract_path: Path = DEFAULT_CONTRACT
) -> dict[str, Any]:
    contract = _load_json(contract_path)
    errors = validate_contract(contract, ROOT_DIR)
    tables: dict[str, pd.DataFrame] = {}
    for filename, schema in contract["tables"].items():
        frame, table_errors = _validate_table(filename, schema, package_dir)
        errors.extend(table_errors)
        if frame is not None and not table_errors:
            tables[filename] = frame
    diagnostics: dict[str, Any] = {}
    if len(tables) == len(contract["tables"]):
        domain_errors, diagnostics = _domain_checks(tables, contract)
        errors.extend(domain_errors)
    return {
        "contract_id": contract["contract_id"],
        "package_dir": str(package_dir.resolve()),
        "valid": not errors,
        "error_count": len(errors),
        "errors": errors,
        "diagnostics": diagnostics,
        "refitting_permitted": contract["locked_model"]["refitting_permitted"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package_dir", type=Path, nargs="?")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--contract-only", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    contract = _load_json(args.contract)
    if args.contract_only:
        errors = validate_contract(contract, ROOT_DIR)
        report = {
            "contract_id": contract["contract_id"],
            "valid": not errors,
            "errors": errors,
        }
    else:
        if args.package_dir is None:
            parser.error("package_dir is required unless --contract-only is used")
        report = validate_package(args.package_dir, args.contract)
    payload = json.dumps(report, indent=2, ensure_ascii=False)
    print(payload)
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload + "\n", encoding="utf-8")
    return 0 if report["valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
