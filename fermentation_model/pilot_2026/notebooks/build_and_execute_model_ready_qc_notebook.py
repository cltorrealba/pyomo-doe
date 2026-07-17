from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import nbformat
from nbclient import NotebookClient


HERE = Path(__file__).resolve().parent
PILOT_DIR = HERE.parent
REPOSITORY_DIR = PILOT_DIR.parents[1]
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026" / "model_dataset"
NOTEBOOK = HERE / "pilot_2026_model_ready_qc.ipynb"


def latest_completed_run() -> Path:
    candidates = []
    for manifest_path in RESULT_ROOT.glob("*/run_manifest.json"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        adapter_path = manifest_path.parent / "adapter_config.json"
        adapter = (
            json.loads(adapter_path.read_text(encoding="utf-8"))
            if adapter_path.exists()
            else {}
        )
        if manifest.get("status") == "completed" and int(
            adapter.get("adapter_schema_version", 0)
        ) >= 2:
            candidates.append(manifest_path.parent)
    if not candidates:
        raise FileNotFoundError(
            "No completed model-dataset run. Execute build_model_dataset.py first."
        )
    return max(candidates, key=lambda path: path.name)


def relative_run_path(run_dir: Path) -> str:
    return run_dir.resolve().relative_to(REPOSITORY_DIR.resolve()).as_posix()


def build_notebook(run_dir: Path) -> nbformat.NotebookNode:
    qc = json.loads((run_dir / "dataset_qc.json").read_text(encoding="utf-8"))
    counts = qc["row_counts"]
    run_relative = relative_run_path(run_dir)
    cells = [
        nbformat.v4.new_markdown_cell(
            f"""# Pilot 2026 — model-ready dataset QC

## tl;dr

The campaign adapter gate is **{qc['verdict']}** for run `{run_dir.name}`. It contains
{counts['primary_observations']:,} primary observations, {counts['temperature_inputs']:,}
active-process temperature rows, {counts['co2_observations']:,} CO₂ bins,
{counts['wine_aroma_observations']:,} wine-aroma observations, and
{counts['condensate_interval_observations']:,} interval-condensate observations.

This is a data-contract pass, **not a calibration pass**. YAN delivery and Oculyze
conversion are parameterized downstream; independent verification and owner review
remain release conditions. No profile in this notebook is approved for physical execution.
"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Context & Methods

### Key Assumptions

- `Sonda1` is executed temperature; controller setpoint is the commanded input.
- Total glucose+fructose is omitted whenever either component is modeled, preventing
  duplicate likelihood contributions.
- Pilot Lot 1 CO₂ is QC-only and cannot enter the model table.
- Cooling and postprocess observations remain in integration QC but are absent from
  these kinetic tables.
- CO₂ is aggregated in deterministic relative-time bins. Signal-lag weights here are
  provisional and calibration must replace them with residual-based ESS.
- NQ and below-LOQ GC rows are interval-censored observations, not zero-valued points.
- Carbon recovery is diagnostic-only because biomass, unmeasured products, sampling
  losses, masked CO₂ intervals, and verified normal reference conditions are absent.
"""
        ),
        nbformat.v4.new_code_cell(
            f"""from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-whitegrid")
pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 180)

def find_repo_root(start: Path) -> Path:
    for candidate in [start, *start.parents]:
        if (candidate / "fermentation_model" / "pilot_2026").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root")

REPO = find_repo_root(Path.cwd().resolve())
RUN_RELATIVE = {run_relative!r}
RUN = REPO / RUN_RELATIVE
print(f"Run: {{RUN_RELATIVE}}")
"""
        ),
        nbformat.v4.new_markdown_cell("## Data"),
        nbformat.v4.new_code_cell(
            """qc = json.loads((RUN / "dataset_qc.json").read_text(encoding="utf-8"))
manifest = json.loads((RUN / "run_manifest.json").read_text(encoding="utf-8"))
primary = pd.read_csv(RUN / "primary_observations.csv", parse_dates=["timestamp"])
temperature = pd.read_csv(RUN / "temperature_inputs.csv", parse_dates=["timestamp"])
co2 = pd.read_csv(RUN / "co2_observations.csv", parse_dates=["timestamp"])
events = pd.read_csv(RUN / "operational_events.csv", parse_dates=["timestamp", "timing_interval_start", "timing_interval_end"])
wine = pd.read_csv(RUN / "wine_aroma_observations.csv", parse_dates=["sample_timestamp"])
condensate = pd.read_csv(RUN / "condensate_observations.csv", parse_dates=["timestamp", "capture_interval_start", "capture_interval_end"])
metadata = pd.read_csv(RUN / "run_metadata.csv")
carbon = pd.read_csv(RUN / "carbon_balance_diagnostic.csv")

display(pd.DataFrame({"rows": qc["row_counts"]}))
assert qc["verdict"] == "PASS"
assert manifest["status"] == "completed"
assert manifest["gate"]["profiles_for_physical_execution"] is False
"""
        ),
        nbformat.v4.new_markdown_cell("## Results"),
        nbformat.v4.new_markdown_cell("### 1. Observation operators, units, and process phases"),
        nbformat.v4.new_code_cell(
            """operator_summary = (
    primary.groupby(["state", "unit", "observation_operator", "process_phase", "calibration_include"])
    .size().rename("n").reset_index()
)
display(operator_summary)

sugar_sets = primary[primary["state"].isin(["glucose", "fructose", "total_sugar"])].groupby("sample_id")["state"].agg(set)
double_counted = sugar_sets.apply(lambda states: "total_sugar" in states and bool(states.intersection({"glucose", "fructose"})))
print(f"Sugar samples double counted: {int(double_counted.sum())}")
assert not double_counted.any()
assert set(primary["process_phase"]) == {"active_process"}
assert primary["calibration_include"].astype(str).str.lower().isin(["true", "1"]).all()
assert set(temperature["process_phase"]) == {"active_process"}
assert temperature["kinetic_include"].astype(str).str.lower().isin(["true", "1"]).all()
"""
        ),
        nbformat.v4.new_markdown_cell("### 2. CO₂ exclusion, reduction, and effective sample size"),
        nbformat.v4.new_code_cell(
            """co2_summary = co2.groupby("experiment_id").agg(
    bins=("observed_flow_ln_min", "size"),
    source_valid_minutes=("source_valid_minutes", "first"),
    lag1_autocorrelation=("lag1_autocorrelation", "first"),
    effective_sample_size=("effective_sample_size_minutes", "first"),
    total_likelihood_weight=("likelihood_weight", "sum"),
).reset_index()
co2_summary["effective_fraction"] = co2_summary["effective_sample_size"] / co2_summary["source_valid_minutes"]
display(co2_summary.round(4))
assert not set(co2["experiment_id"].astype(str)).intersection({"26134", "26135", "26136"})
assert np.allclose(co2_summary["effective_sample_size"], co2_summary["total_likelihood_weight"])

fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(co2_summary["experiment_id"].astype(str), co2_summary["effective_fraction"], color="#4c78a8")
ax.set_ylabel("Effective / valid minute observations")
ax.set_xlabel("Experiment")
ax.set_title("CO₂ provisional signal-lag fraction — replaced after calibration")
ax.set_ylim(0, max(0.05, co2_summary["effective_fraction"].max() * 1.15))
fig.tight_layout()
plt.show()
"""
        ),
        nbformat.v4.new_markdown_cell("### 3. GC censoring and interval capture operators"),
        nbformat.v4.new_code_cell(
            """gc_summary = pd.DataFrame({
    "wine": wine["model_observation_type"].value_counts(),
    "condensate": condensate["mix_model_observation_type"].value_counts(),
}).fillna(0).astype(int)
display(gc_summary)

wine_censored = wine["model_observation_type"].eq("left_censored")
mix_censored = condensate["mix_model_observation_type"].eq("left_censored")
assert wine.loc[wine_censored, "observed_value"].isna().all()
assert wine.loc[wine_censored, "upper_bound"].gt(0).all()
assert condensate.loc[mix_censored, "observed_value"].isna().all()
assert condensate.loc[mix_censored, "upper_bound"].gt(0).all()
assert (condensate["capture_interval_h"] > 0).all()
display(condensate[["experiment_id", "mix_id", "analyte", "capture_interval_start", "capture_interval_end", "capture_interval_h", "observed_value", "upper_bound", "mix_model_observation_type"]].head(12))
"""
        ),
        nbformat.v4.new_markdown_cell("### 4. Density-trigger timing uncertainty"),
        nbformat.v4.new_code_cell(
            """uncertain = events[events["timing_uncertain"].astype(str).str.lower().eq("true")][
    ["experiment_id", "timestamp", "relative_time_h", "timing_interval_start", "timing_interval_end", "timing_interval_start_h", "timing_interval_end_h", "timing_uncertainty_basis"]
]
display(uncertain)
assert len(uncertain) == 9
assert uncertain[["timing_interval_start", "timing_interval_end"]].notna().all().all()
"""
        ),
        nbformat.v4.new_markdown_cell("### 5. Metadata coverage and carbon diagnostic"),
        nbformat.v4.new_code_cell(
            """display(metadata[["experiment_id", "lot", "reactor", "yeast", "storage_history", "capture_system", "co2_sensor", "temperature_sensor", "sensor_metadata_status"]])
display(carbon.round(4))

valid_carbon = carbon[carbon["diagnostic_status"].eq("diagnostic_only_incomplete")]
fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(valid_carbon["experiment_id"].astype(str), valid_carbon["partial_observed_carbon_recovery"], color="#f58518")
ax.axhline(1.0, color="black", linewidth=1, linestyle="--")
ax.set_ylabel("Partial observed carbon recovery")
ax.set_xlabel("Experiment")
ax.set_title("Diagnostic only — not a closed carbon balance")
fig.tight_layout()
plt.show()
"""
        ),
        nbformat.v4.new_markdown_cell("### 6. Automated gate table"),
        nbformat.v4.new_code_cell(
            """checks = pd.Series(qc["checks"], name="PASS").to_frame()
display(checks)
assert checks["PASS"].all()
print("PHASE A MODEL-READY ADAPTER: PASS")
print("CALIBRATION STATUS: evaluate the separate immutable calibration run")
print("PHYSICAL EXECUTION STATUS: EXPLORATORY ONLY; no schedule issued")
"""
        ),
        nbformat.v4.new_markdown_cell(
            """## Takeaways

- The adapter enforces the Pilot 2026 observation contract and fails if Lot 1 CO₂
  re-enters calibration.
- Cooling is formally absent from kinetic inputs.
- Correlated minute CO₂ data no longer dominate by raw row count; final ESS is
  calculated from fit residuals in the calibration workflow.
- GC censoring and interval-condensate accumulation remain explicit in model space.
- Density-triggered second pulses are intervals/latent timing inputs, not exact timestamps.
- The next authorized step is the conditional computational MBDoE workflow, not an
  executable thermal or nutrition schedule.
"""
        ),
    ]
    return nbformat.v4.new_notebook(
        cells=cells,
        metadata={
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
            "pilot_2026_model_dataset_run": run_relative,
        },
    )


def validate_notebook(notebook: nbformat.NotebookNode) -> None:
    nbformat.validate(notebook)
    errors = [
        output
        for cell in notebook.cells
        if cell.cell_type == "code"
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    if errors:
        raise RuntimeError(f"Executed notebook contains errors: {errors}")
    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    if not code_cells or any(cell.execution_count is None for cell in code_cells):
        raise RuntimeError("Notebook did not execute all code cells")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and execute model-ready QC notebook")
    parser.add_argument("--run-dir", type=Path)
    args = parser.parse_args()
    run_dir = (args.run_dir or latest_completed_run()).resolve()
    notebook = build_notebook(run_dir)
    NOTEBOOK.write_text(nbformat.writes(notebook), encoding="utf-8")
    executed = NotebookClient(
        notebook,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(REPOSITORY_DIR)}},
    ).execute(cwd=str(REPOSITORY_DIR))
    validate_notebook(executed)
    NOTEBOOK.write_text(nbformat.writes(executed), encoding="utf-8")
    print(f"Built and executed {NOTEBOOK}")
    print(f"Model dataset run: {relative_run_path(run_dir)}")


if __name__ == "__main__":
    main()
