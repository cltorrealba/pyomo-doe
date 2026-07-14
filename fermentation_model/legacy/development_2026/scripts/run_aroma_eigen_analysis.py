from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_joint_campaign_doe"


def load_matrix(path: Path, parameters: tuple[str, ...]) -> np.ndarray:
    frame = pd.read_csv(path, index_col=0)
    frame = frame.reindex(index=parameters, columns=parameters)
    if frame.isnull().any().any():
        raise RuntimeError(f"Matrix {path} is missing entries after parameter reindexing.")
    arr = np.asarray(frame, dtype=float)
    return 0.5 * (arr + arr.T)


def load_design_fim(design: str, selected: pd.DataFrame, prior_fim: np.ndarray, parameters: tuple[str, ...]) -> np.ndarray:
    current = prior_fim.copy()
    for name in selected.loc[selected["design"].eq(design), "candidate"].astype(str):
        current += load_matrix(RESULTS_DIR / f"fim_aroma_{name}.csv", parameters)
    return 0.5 * (current + current.T)


def parameter_group(parameter: str) -> str:
    return "synthesis" if parameter.startswith("k_") else "fermentation"


def eigen_analysis(design: str, fim: np.ndarray, parameters: tuple[str, ...]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    eigvals, eigvecs = np.linalg.eigh(0.5 * (fim + fim.T))
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    max_eig = float(np.max(eigvals))
    positive_floor = max(max_eig * 1e-12, np.finfo(float).tiny)
    clipped = np.clip(eigvals, positive_floor, None)
    weights = clipped / np.sum(clipped)
    effective_rank = float(np.exp(-np.sum(weights * np.log(weights))))
    numeric_rank = int(np.sum(eigvals > positive_floor))
    spectrum_rows: list[dict[str, object]] = []
    loading_rows: list[dict[str, object]] = []
    for idx, eig in enumerate(eigvals, start=1):
        vector = eigvecs[:, idx - 1]
        # Fix sign for reproducible tables; eigenvector sign is arbitrary.
        max_abs_idx = int(np.argmax(np.abs(vector)))
        if vector[max_abs_idx] < 0:
            vector = -vector
        rel = float(eig / max_eig) if max_eig > 0 else np.nan
        spectrum_rows.append(
            {
                "design": design,
                "direction_rank_weak_to_strong": idx,
                "eigenvalue": float(eig),
                "relative_eigenvalue": rel,
                "log10_eigenvalue": float(np.log10(max(float(eig), positive_floor))),
                "numeric_rank": numeric_rank,
                "condition_number": float(np.max(clipped) / np.min(clipped)),
                "effective_rank": effective_rank,
            }
        )
        for parameter, loading in zip(parameters, vector):
            loading_rows.append(
                {
                    "design": design,
                    "direction_rank_weak_to_strong": idx,
                    "eigenvalue": float(eig),
                    "relative_eigenvalue": rel,
                    "parameter": parameter,
                    "group": parameter_group(parameter),
                    "loading": float(loading),
                    "abs_loading": float(abs(loading)),
                    "loading_sq": float(loading * loading),
                }
            )
    return spectrum_rows, loading_rows


def write_plots(spectrum: pd.DataFrame, loadings: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    for design, group in spectrum.groupby("design", sort=False):
        ax.plot(
            group["direction_rank_weak_to_strong"],
            group["relative_eigenvalue"],
            marker="o",
            linewidth=1.8,
            label=design,
        )
    ax.set_yscale("log")
    ax.set_xlabel("direction rank, weak to strong")
    ax.set_ylabel("relative eigenvalue")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "aroma_eigen_spectrum.png", dpi=180)
    plt.close(fig)

    weak = loadings[loadings["direction_rank_weak_to_strong"].le(5)].copy()
    weak["direction"] = weak["design"] + " d" + weak["direction_rank_weak_to_strong"].astype(str)
    pivot = weak.pivot_table(index="parameter", columns="direction", values="abs_loading", aggfunc="first").fillna(0.0)
    order = pivot.max(axis=1).sort_values(ascending=False).index[:14]
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    image = ax.imshow(pivot.loc[order].to_numpy(dtype=float), aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_title("Dominant loadings in weakest eigendirections")
    fig.colorbar(image, ax=ax, label="abs loading")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "aroma_weak_eigendirection_loadings.png", dpi=180)
    plt.close(fig)


def write_report(summary: pd.DataFrame, top_loadings: pd.DataFrame) -> None:
    lines = [
        "# Eigenvalue analysis for joint aroma campaign",
        "",
        "This analysis decomposes the posterior FIM in scaled parameter coordinates.",
        "Small eigenvalues correspond to weakly informed linear combinations of parameters.",
        "",
        "## Design summary",
        "",
        summary.to_markdown(index=False),
        "",
        "## Dominant weak-direction loadings",
        "",
        top_loadings.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- All tested posterior FIMs are numerically full rank under the `lambda_max * 1e-12` threshold.",
        "- The hybrid design has the better weakest relative eigenvalue and lower condition number.",
        "- The weakest directions are mixed fermentation-aroma combinations, not isolated single parameters.",
        "- The recurring weak components are `sN`, `m0`, `qXG/qXF`, and growth-phase aroma yields, especially `k_EO_growth`.",
        "- Stationary aroma yields are not among the weakest eigendirection loadings, so the remaining aroma weakness is mostly tied to growth/uptake coupling.",
        "- Because partition parameters were not included, this analysis should not be read as evidence that partition corrections are identifiable.",
        "",
    ]
    (RESULTS_DIR / "aroma_eigen_analysis_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    metadata = json.loads((RESULTS_DIR / "aroma_campaign_metadata.json").read_text(encoding="utf-8"))
    parameters = tuple(metadata["parameters"])
    prior_fim = load_matrix(RESULTS_DIR / "aroma_prior_fim.csv", parameters)
    selected = pd.read_csv(RESULTS_DIR / "aroma_dopt_benchmark_selected.csv")

    designs = ["hybrid_current", "dopt_exchange_from_hybrid"]
    spectrum_rows: list[dict[str, object]] = []
    loading_rows: list[dict[str, object]] = []
    for design in designs:
        fim = load_design_fim(design, selected, prior_fim, parameters)
        spectrum, loadings = eigen_analysis(design, fim, parameters)
        spectrum_rows.extend(spectrum)
        loading_rows.extend(loadings)

    spectrum_frame = pd.DataFrame(spectrum_rows)
    loading_frame = pd.DataFrame(loading_rows)
    top_loadings = (
        loading_frame[loading_frame["direction_rank_weak_to_strong"].le(5)]
        .sort_values(["design", "direction_rank_weak_to_strong", "abs_loading"], ascending=[True, True, False])
        .groupby(["design", "direction_rank_weak_to_strong"])
        .head(5)
        .reset_index(drop=True)
    )

    summary = (
        spectrum_frame.groupby("design", sort=False)
        .agg(
            n_directions=("direction_rank_weak_to_strong", "max"),
            min_relative_eigenvalue=("relative_eigenvalue", "min"),
            condition_number=("condition_number", "first"),
            numeric_rank=("numeric_rank", "first"),
            effective_rank=("effective_rank", "first"),
        )
        .reset_index()
    )

    spectrum_frame.to_csv(RESULTS_DIR / "aroma_eigen_spectrum.csv", index=False)
    loading_frame.to_csv(RESULTS_DIR / "aroma_eigendirection_loadings.csv", index=False)
    top_loadings.to_csv(RESULTS_DIR / "aroma_weak_eigendirection_top_loadings.csv", index=False)
    summary.to_csv(RESULTS_DIR / "aroma_eigen_summary.csv", index=False)
    write_plots(spectrum_frame, loading_frame)
    write_report(summary, top_loadings)

    print(summary.to_string(index=False))
    print("\nWeakest direction top loadings:")
    print(top_loadings[top_loadings["direction_rank_weak_to_strong"].le(2)].to_string(index=False))


if __name__ == "__main__":
    main()
