from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import run_aroma_campaign_doe as aroma_doe


SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_joint_campaign_doe"
OUT_DIR = SCRIPT_DIR / "results" / "operational_campaign_schedule"

RECOMMENDED_BATCHES = {
    1: (
        "fructose_rich_iG_probe",
        "low_temp_growth_separation_plus_co2",
        "yan_saturation_scan",
    ),
    2: (
        "cold_synthesis_hot_stripping",
        "combined_stress_long_horizon",
        "glucose_rich_growth_yield",
    ),
    3: (
        "high_biomass_low_N_maintenance",
        "ethanol_initial_challenge",
        "glucose_pulse_after_N_depletion",
    ),
}


def load_matrix(path: Path, parameters: tuple[str, ...]) -> np.ndarray:
    frame = pd.read_csv(path, index_col=0).reindex(index=parameters, columns=parameters)
    arr = frame.to_numpy(dtype=float)
    return 0.5 * (arr + arr.T)


def campaign_metrics(names: list[str], prior: np.ndarray, fims: dict[str, np.ndarray], parameters: tuple[str, ...]) -> dict[str, object]:
    fim = prior.copy()
    for name in names:
        fim += fims[name]
    fim_metrics = aroma_doe.fim_metrics(fim, parameters)
    reductions = aroma_doe.variance_reduction(prior, fim, parameters)
    return {
        "experiments": " | ".join(names),
        "n_experiments": len(names),
        "hybrid_score": aroma_doe.campaign_score(prior, fim, parameters),
        "logdet": fim_metrics["logdet"],
        "min_relative_eigenvalue": fim_metrics["min_relative_eigenvalue"],
        "condition_number": fim_metrics["condition_number"],
        "mean_var_reduction": reductions["mean_var_reduction"],
        "worst_var_reduction": reductions["worst_var_reduction"],
        "sN_var_reduction": reductions["var_reduction_sN"],
        "qN_var_reduction": reductions["var_reduction_qN"],
        "mu0_var_reduction": reductions["var_reduction_mu0"],
        "k_EO_growth_var_reduction": reductions["var_reduction_k_EO_growth"],
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((RESULTS_DIR / "aroma_campaign_metadata.json").read_text(encoding="utf-8"))
    parameters = tuple(metadata["parameters"])
    prior = load_matrix(RESULTS_DIR / "aroma_prior_fim.csv", parameters)
    selected = pd.read_csv(RESULTS_DIR / "aroma_campaign_selected.csv")["candidate"].astype(str).tolist()
    fims = {name: load_matrix(RESULTS_DIR / f"fim_aroma_{name}.csv", parameters) for name in selected}

    combo_rows = []
    for combo in combinations(selected, 3):
        combo_rows.append(campaign_metrics(list(combo), prior, fims, parameters))
    combos = pd.DataFrame(combo_rows).sort_values("hybrid_score", ascending=False)
    combos.to_csv(OUT_DIR / "batch1_combo_ranking_selected9.csv", index=False)

    sequence_rows = []
    cumulative: list[str] = []
    for batch_id in sorted(RECOMMENDED_BATCHES):
        names = list(RECOMMENDED_BATCHES[batch_id])
        cumulative.extend(names)
        row = campaign_metrics(cumulative, prior, fims, parameters)
        row["batch"] = batch_id
        row["batch_experiments"] = " | ".join(names)
        sequence_rows.append(row)
    sequence = pd.DataFrame(sequence_rows)
    sequence.to_csv(OUT_DIR / "recommended_batch_sequence_metrics.csv", index=False)

    report = [
        "# Batch prioritization",
        "",
        "The selected nine-fermentation campaign is executed as three batches of three fermentations.",
        "The first batch is chosen for adaptive learning, not only for maximum global logdet.",
        "",
        "## Recommended batches",
        "",
        sequence[[
            "batch",
            "batch_experiments",
            "n_experiments",
            "logdet",
            "mean_var_reduction",
            "worst_var_reduction",
            "sN_var_reduction",
            "condition_number",
        ]].to_markdown(index=False),
        "",
        "## Top first-batch combinations by hybrid score",
        "",
        combos.head(10).to_markdown(index=False),
        "",
        "## Rationale",
        "",
        "The maximum-hybrid-score three-experiment batch is globally informative, but it leaves the weakest nitrogen saturation direction less protected.",
        "For an adaptive first batch, the recommended set includes `yan_saturation_scan` because `sN/qN` is the most important practical bottleneck before committing to the remaining six fermentations.",
        "",
    ]
    (OUT_DIR / "batch_prioritization_report.md").write_text("\n".join(report), encoding="utf-8")
    print(sequence.to_string(index=False))


if __name__ == "__main__":
    main()
