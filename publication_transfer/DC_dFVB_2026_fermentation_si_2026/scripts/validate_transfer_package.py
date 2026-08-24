"""Dependency-free structural and semantic checks for the transferred CSVs."""

from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PACKAGE_ROOT / "data"


def read_csv(name: str) -> tuple[list[str], list[dict[str, str]]]:
    with (DATA_DIR / name).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def as_float(value: str) -> float | None:
    value = str(value).strip()
    if not value:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def require(condition: bool, message: str, failures: list[str]) -> None:
    if condition:
        print(f"PASS: {message}")
    else:
        print(f"FAIL: {message}")
        failures.append(message)


def main() -> int:
    failures: list[str] = []
    required = {
        "process_observations.csv",
        "aroma_observations.csv",
        "nitrogen_events.csv",
        "process_registry.csv",
        "quality_checks.csv",
        "source_manifest.csv",
        "measurement_dictionary.csv",
    }
    existing = {path.name for path in DATA_DIR.glob("*.csv")}
    require(required.issubset(existing), "all required CSV files exist", failures)
    if failures:
        return 1

    obs_columns, obs = read_csv("process_observations.csv")
    aroma_columns, aromas = read_csv("aroma_observations.csv")
    event_columns, events = read_csv("nitrogen_events.csv")
    registry_columns, registry = read_csv("process_registry.csv")

    processes = {row["process_token"] for row in obs}
    require(len(processes) == 7, "seven processes are in the submission dataset", failures)
    require(len(registry) == 8, "registry records seven included and one excluded process", failures)
    require(sum(as_bool(row["include_submission_dataset"]) for row in registry) == 7, "exactly seven registry rows are included", failures)
    require(not any(as_bool(row["independent_validation"]) for row in registry), "no retrospective process is labeled independent validation", failures)
    require(all(re.fullmatch(r"P25_\d{2}", token) for token in processes), "only publication tokens identify processes", failures)

    closure_ok = True
    closure_count = 0
    for row in obs:
        yan = as_float(row["yan_reported_n_mg_l"])
        pan = as_float(row["pan_n_mg_l"])
        ammonia = as_float(row["ammonia_compound_mg_l"])
        if yan is None or pan is None or ammonia is None:
            continue
        closure_count += 1
        closure_ok = closure_ok and math.isclose(yan, pan + 0.82 * ammonia, abs_tol=1e-9)
    require(closure_count > 0 and closure_ok, "nitrogen fields reproduce source YAN exactly", failures)

    aroma_ok = True
    aroma_count = 0
    for row in aromas:
        total = as_float(row["total_equivalent_mg_l"])
        condensate = as_float(row["condensate_equivalent_mg_l"])
        retained = as_float(row["retained_wine_mg_l"])
        if total is None or condensate is None or retained is None:
            continue
        aroma_count += 1
        aroma_ok = aroma_ok and math.isclose(total - condensate, retained, abs_tol=1e-12)
    require(aroma_count > 0 and aroma_ok, "retained aroma equals total minus condensate", failures)
    require(all(row["recorded_unit_basis"] == "unresolved_in_source_workbook" for row in events), "pulse units are not silently invented", failures)

    forbidden_columns = {"batch", "experiment_id", "sample_id", "fecha_hora", "timestamp"}
    file_columns = {
        "process_observations.csv": obs_columns,
        "aroma_observations.csv": aroma_columns,
        "nitrogen_events.csv": event_columns,
        "process_registry.csv": registry_columns,
    }
    for filename, columns in file_columns.items():
        require(not forbidden_columns.intersection(columns), f"{filename} has no private-ID or calendar-time columns", failures)
        contents = (DATA_DIR / filename).read_text(encoding="utf-8")
        require(re.search(r"\b25\d{3}\b", contents) is None, f"{filename} contains no five-digit industrial process IDs", failures)

    print(f"\nValidation completed with {len(failures)} failure(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
