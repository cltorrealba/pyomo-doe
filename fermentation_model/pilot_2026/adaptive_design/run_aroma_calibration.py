from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    build_forcings,
    evaluate_gate,
    fit_species,
    load_json,
    load_partition_surrogates,
)
from pilot_2026.adaptive_design.pilot_calibration import load_tables  # noqa: E402
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
CALIBRATION_CONFIG_PATH = ADAPTIVE_DIR / "calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def main() -> None:
    config = load_json(CONFIG_PATH)
    calibration_config = load_json(CALIBRATION_CONFIG_PATH)
    model_run = REPOSITORY_DIR / config["source_contract"]["model_dataset_run"]
    baseline_run = REPOSITORY_DIR / config["source_contract"]["fixed_yan_baseline_run"]
    partitions, partition_provenance = load_partition_surrogates(config, REPOSITORY_DIR)
    tables = load_tables(model_run)
    run_dir = create_immutable_run_directory(RESULT_ROOT, "aroma_calibration", config)
    outputs: list[Path] = []
    fits = []
    forcing_metadata = {}
    for species, analyte_label in config["priority_analytes"].items():
        forcings, wine, condensate, _theta = build_forcings(
            tables,
            calibration_config,
            baseline_run,
            species,
            analyte_label,
            partitions[species],
        )
        fit = fit_species(species, forcings, wine, condensate, config)
        fits.append(fit)
        forcing_metadata[species] = {
            "runs": sorted(forcings),
            "wine_rows": len(wine),
            "condensate_rows": len(condensate),
        }
    gate = evaluate_gate(fits)
    frames = {
        "aroma_parameter_estimates.csv": pd.concat(
            [fit.parameter_table for fit in fits], ignore_index=True
        ),
        "aroma_multistart_summary.csv": pd.concat(
            [fit.multistart_summary for fit in fits], ignore_index=True
        ),
        "aroma_wine_predictions.csv": pd.concat(
            [fit.wine_predictions for fit in fits], ignore_index=True
        ),
        "aroma_condensate_predictions.csv": pd.concat(
            [fit.condensate_predictions for fit in fits], ignore_index=True
        ),
        "aroma_profile_likelihood.csv": pd.concat(
            [fit.profiles for fit in fits], ignore_index=True
        ),
        "aroma_covariance_log_space.csv": pd.concat(
            [
                fit.covariance.rename_axis("row_parameter").reset_index().assign(
                    species=fit.species
                )
                for fit in fits
            ],
            ignore_index=True,
        ),
    }
    for name, frame in frames.items():
        path = run_dir / name
        frame.to_csv(path, index=False)
        outputs.append(path)
    json_payloads = {
        "aroma_calibration_config.json": config,
        "partition_surrogate_provenance.json": {
            "source": partition_provenance,
            "models": partitions,
            "silent_fallback_used": False,
        },
        "aroma_validation.json": {
            "species": {fit.species: fit.validation for fit in fits},
            "forcing_metadata": forcing_metadata,
        },
        "aroma_calibration_gate.json": gate,
    }
    for name, payload in json_payloads.items():
        path = run_dir / name
        write_json(path, payload)
        outputs.append(path)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="aroma_calibration",
        config=config,
        sources={
            "aroma_config": CONFIG_PATH,
            "calibration_config": CALIBRATION_CONFIG_PATH,
            "model_dataset_manifest": model_run / "run_manifest.json",
            "fixed_yan_calibration_manifest": baseline_run / "run_manifest.json",
            "partition_surrogate_source": REPOSITORY_DIR
            / config["thermodynamics"]["surrogate_source"],
        },
        code_paths=[Path(__file__), ADAPTIVE_DIR / "pilot_aroma_calibration.py"],
        random_seeds=[int(config["optimization"]["seed"])],
        status="completed" if gate["verdict"] != "FAIL" else "validation_failed",
        convergence={fit.species: fit.validation for fit in fits},
        gate=gate,
        outputs=outputs,
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(
        json.dumps(
            {
                "verdict": gate["verdict"],
                "run_directory": str(run_dir),
                "loss_identifiability": gate["loss_identifiability"],
                "profiles_for_physical_execution": False,
            },
            indent=2,
        )
    )
    if gate["verdict"] == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
