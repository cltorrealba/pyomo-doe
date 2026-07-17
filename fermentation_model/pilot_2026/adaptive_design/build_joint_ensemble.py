from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "joint_ensemble_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _nearest_psd(covariance: np.ndarray) -> np.ndarray:
    symmetric = 0.5 * (covariance + covariance.T)
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1.0)
    values = np.maximum(values, scale * 1e-10)
    return (vectors * values) @ vectors.T


def _primary_samples(config: dict, source: Path, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    estimates = pd.read_csv(source / "primary_parameter_estimates.csv")
    summary = pd.read_csv(source / "primary_multistart_summary.csv")
    covariance = pd.read_csv(source / "primary_covariance_log_space.csv", index_col=0)
    names = estimates["parameter"].tolist()
    centres = summary[summary["success"]].copy()
    estimate_columns = [f"estimate__{name}" for name in names]
    centres = centres.dropna(subset=estimate_columns).reset_index(drop=True)
    if len(centres) < int(config["primary"]["minimum_distinct_multistart_centres"]):
        raise ValueError("Too few successful primary multistart centres")
    objective = centres["robust_objective"].to_numpy(dtype=float)
    delta = objective - objective.min()
    temperature = float(config["primary"]["centre_weight_temperature"])
    weights = np.exp(-0.5 * delta / temperature)
    weights /= weights.sum()
    size = int(config["ensemble_size"])
    # Guarantee representation of each retained basin before weighted sampling.
    centre_indices = np.concatenate(
        [np.arange(len(centres)), rng.choice(len(centres), size=size - len(centres), p=weights)]
    )[:size]
    cov = covariance.loc[names, names].to_numpy(dtype=float)
    cov = _nearest_psd(cov) * float(config["primary"]["local_covariance_scale"]) ** 2
    lower = np.log(estimates.set_index("parameter").loc[names, "lower_bound"].to_numpy(float))
    upper = np.log(estimates.set_index("parameter").loc[names, "upper_bound"].to_numpy(float))
    rows = []
    for member, centre_index in enumerate(centre_indices):
        centre = np.log(centres.loc[int(centre_index), estimate_columns].to_numpy(dtype=float))
        draw = np.clip(rng.multivariate_normal(centre, cov), lower, upper)
        row = {
            "ensemble_member": member,
            "primary_centre_index": int(centre_index),
            "primary_centre_stage": centres.loc[int(centre_index), "stage"],
            "primary_centre_objective": float(centres.loc[int(centre_index), "robust_objective"]),
        }
        row.update({f"primary__{name}": math.exp(value) for name, value in zip(names, draw)})
        rows.append(row)
    metadata = {
        "parameter_names": names,
        "available_centres": int(len(centres)),
        "represented_centres": int(len(set(centre_indices.tolist()))),
        "centre_weights": [float(value) for value in weights],
    }
    return pd.DataFrame(rows), metadata


def _aroma_samples(config: dict, source: Path, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    estimates = pd.read_csv(source / "aroma_parameter_estimates.csv")
    covariance = pd.read_csv(source / "aroma_covariance_log_space.csv")
    gate = load_json(source / "aroma_calibration_gate.json")
    size = int(config["ensemble_size"])
    output = pd.DataFrame({"ensemble_member": np.arange(size)})
    diagnostics = {}
    for species, species_rows in estimates.groupby("species", sort=True):
        names = species_rows["parameter"].tolist()
        centre = np.log(species_rows["estimate"].to_numpy(dtype=float))
        lower = np.log(species_rows["lower_bound"].to_numpy(dtype=float))
        upper = np.log(species_rows["upper_bound"].to_numpy(dtype=float))
        cov_rows = covariance[covariance["species"].eq(species)].set_index("row_parameter")
        cov = _nearest_psd(cov_rows.loc[names, names].to_numpy(dtype=float))
        cov *= float(config["aroma"]["local_covariance_scale"]) ** 2
        weak = set(gate["weak_parameter_directions"][species])
        sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
        for index, name in enumerate(names):
            if name in weak:
                floor = float(config["aroma"]["weak_direction_minimum_log_sd"].get(name, 2.0))
                cov[index, index] = max(cov[index, index], floor**2)
            else:
                ceiling = float(config["aroma"]["identified_direction_maximum_log_sd"])
                cov[index, index] = min(cov[index, index], ceiling**2)
        cov = _nearest_psd(cov)
        draws = np.clip(rng.multivariate_normal(centre, cov, size=size), lower, upper)
        for index, name in enumerate(names):
            output[f"aroma__{species}__{name}"] = np.exp(draws[:, index])
        diagnostics[species] = {
            "weak_parameters": sorted(weak),
            "sample_log_sd": {
                name: float(np.std(draws[:, index], ddof=1)) for index, name in enumerate(names)
            },
        }
    return output, diagnostics


def build_gate(config: dict, ensemble: pd.DataFrame, primary: dict, aroma: dict) -> dict:
    weak_floor = config["aroma"]["weak_direction_minimum_log_sd"]
    weak_coverage = all(
        details["sample_log_sd"][parameter] >= 0.75 * float(weak_floor.get(parameter, 2.0))
        for details in aroma.values()
        for parameter in details["weak_parameters"]
    )
    checks = {
        "requested_ensemble_size": len(ensemble) == int(config["ensemble_size"]),
        "all_values_finite": bool(np.isfinite(ensemble.select_dtypes(include=[np.number])).all().all()),
        "primary_alternative_minima_represented": primary["represented_centres"]
        >= int(config["primary"]["minimum_distinct_multistart_centres"]),
        "weak_aroma_directions_remain_broad": bool(weak_coverage),
        "point_estimate_only_design_prohibited": bool(
            config["policy"]["point_estimate_only_design_prohibited"]
        ),
    }
    passed = all(checks.values())
    return {
        "gate": "phase_C_joint_uncertainty_ensemble",
        "verdict": "PASS" if passed else "FAIL",
        "checks": checks,
        "ensemble_members": len(ensemble),
        "profiles_for_physical_execution": False,
    }


def main() -> None:
    config = load_json(CONFIG_PATH)
    primary_source = REPOSITORY_DIR / config["primary_source_run"]
    aroma_source = REPOSITORY_DIR / config["aroma_source_run"]
    rng = np.random.default_rng(int(config["seed"]))
    primary, primary_meta = _primary_samples(config, primary_source, rng)
    aroma, aroma_meta = _aroma_samples(config, aroma_source, rng)
    ensemble = primary.merge(aroma, on="ensemble_member", validate="one_to_one")
    gate = build_gate(config, ensemble, primary_meta, aroma_meta)
    run_dir = create_immutable_run_directory(RESULT_ROOT, "joint_ensemble", config)
    ensemble_path = run_dir / "joint_parameter_ensemble.csv"
    ensemble.to_csv(ensemble_path, index=False)
    diagnostics_path = run_dir / "ensemble_diagnostics.json"
    write_json(diagnostics_path, {"primary": primary_meta, "aroma": aroma_meta})
    gate_path = run_dir / "ensemble_gate.json"
    write_json(gate_path, gate)
    config_path = run_dir / "joint_ensemble_config.json"
    write_json(config_path, config)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="joint_uncertainty_ensemble",
        config=config,
        sources={
            "primary_manifest": primary_source / "run_manifest.json",
            "aroma_manifest": aroma_source / "run_manifest.json",
            "ensemble_config": CONFIG_PATH,
        },
        code_paths=[Path(__file__)],
        random_seeds=[int(config["seed"])],
        status="completed" if gate["verdict"] == "PASS" else "validation_failed",
        convergence={},
        gate=gate,
        outputs=[ensemble_path, diagnostics_path, gate_path, config_path],
    )
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if gate["verdict"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
