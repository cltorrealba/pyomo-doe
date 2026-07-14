from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_aroma_campaign_doe as aroma_doe


RESULTS_DIR = Path(__file__).resolve().parent / "results" / "aroma_joint_campaign_doe"


def load_matrix(path: Path, parameters: tuple[str, ...]) -> np.ndarray:
    frame = pd.read_csv(path, index_col=0)
    frame = frame.reindex(index=parameters, columns=parameters)
    if frame.isnull().any().any():
        missing = frame.columns[frame.isnull().any()].tolist()
        raise RuntimeError(f"Matrix {path} is missing parameters: {missing}")
    arr = np.asarray(frame, dtype=float)
    return 0.5 * (arr + arr.T)


def load_candidate_fims(results_dir: Path, parameters: tuple[str, ...]) -> dict[str, np.ndarray]:
    ranking = pd.read_csv(results_dir / "aroma_candidate_ranking.csv")
    fims: dict[str, np.ndarray] = {}
    for candidate in ranking.loc[ranking["status"] == "ok", "candidate"]:
        fim_path = results_dir / f"fim_aroma_{candidate}.csv"
        if fim_path.exists():
            fims[str(candidate)] = load_matrix(fim_path, parameters)
    return fims


def campaign_fim(prior_fim: np.ndarray, fims: dict[str, np.ndarray], names: list[str]) -> np.ndarray:
    current = prior_fim.copy()
    for name in names:
        current = current + fims[name]
    return 0.5 * (current + current.T)


def dopt_logdet(fim: np.ndarray, parameters: tuple[str, ...]) -> float:
    return float(aroma_doe.fim_metrics(fim, parameters)["logdet"])


def dopt_greedy(
    prior_fim: np.ndarray,
    fims: dict[str, np.ndarray],
    parameters: tuple[str, ...],
    campaign_size: int,
) -> list[str]:
    selected: list[str] = []
    remaining = set(fims)
    current = prior_fim.copy()
    for _order in range(int(campaign_size)):
        best_name = None
        best_score = -np.inf
        for name in sorted(remaining):
            score = dopt_logdet(current + fims[name], parameters)
            if score > best_score:
                best_name = name
                best_score = score
        if best_name is None:
            break
        selected.append(best_name)
        remaining.remove(best_name)
        current = current + fims[best_name]
    return selected


def dopt_exchange(
    seed: list[str],
    prior_fim: np.ndarray,
    fims: dict[str, np.ndarray],
    parameters: tuple[str, ...],
    tolerance: float = 1e-8,
) -> tuple[list[str], list[dict[str, object]]]:
    selected = list(seed)
    exchange_log: list[dict[str, object]] = []
    current_score = dopt_logdet(campaign_fim(prior_fim, fims, selected), parameters)
    improved = True
    iteration = 0
    while improved:
        improved = False
        best_payload = None
        remaining = sorted(set(fims) - set(selected))
        for out_idx, out_name in enumerate(selected):
            for in_name in remaining:
                trial = list(selected)
                trial[out_idx] = in_name
                trial_score = dopt_logdet(campaign_fim(prior_fim, fims, trial), parameters)
                gain = trial_score - current_score
                if gain > tolerance and (best_payload is None or gain > best_payload["gain_logdet"]):
                    best_payload = {
                        "iteration": iteration + 1,
                        "swap_out": out_name,
                        "swap_in": in_name,
                        "previous_logdet": current_score,
                        "new_logdet": trial_score,
                        "gain_logdet": gain,
                        "out_index": out_idx,
                    }
        if best_payload is not None:
            selected[int(best_payload["out_index"])] = str(best_payload["swap_in"])
            current_score = float(best_payload["new_logdet"])
            exchange_log.append(best_payload)
            improved = True
            iteration += 1
    return selected, exchange_log


def parameter_group(parameter: str) -> str:
    return "synthesis" if parameter.startswith("k_") else "fermentation"


def design_summary(
    design: str,
    names: list[str],
    prior_fim: np.ndarray,
    fims: dict[str, np.ndarray],
    parameters: tuple[str, ...],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    fim = campaign_fim(prior_fim, fims, names)
    metrics = aroma_doe.fim_metrics(fim, parameters)
    reductions = aroma_doe.variance_reduction(prior_fim, fim, parameters)
    per_parameter: list[dict[str, object]] = []
    min_name = None
    min_reduction = np.inf
    for parameter in parameters:
        reduction = float(reductions[f"var_reduction_{parameter}"])
        if reduction < min_reduction:
            min_reduction = reduction
            min_name = parameter
        per_parameter.append(
            {
                "design": design,
                "parameter": parameter,
                "group": parameter_group(parameter),
                "var_ratio": float(reductions[f"var_ratio_{parameter}"]),
                "var_reduction": reduction,
            }
        )
    summary = {
        "design": design,
        "n_experiments": len(names),
        "logdet": float(metrics["logdet"]),
        "min_relative_eigenvalue": float(metrics["min_relative_eigenvalue"]),
        "condition_number": float(metrics["condition_number"]),
        "trace_inv": float(metrics["trace_inv"]),
        "mean_var_reduction": float(reductions["mean_var_reduction"]),
        "worst_var_reduction": float(reductions["worst_var_reduction"]),
        "worst_parameter": min_name,
        "sN_var_reduction": float(reductions["var_reduction_sN"]),
        "qN_var_reduction": float(reductions["var_reduction_qN"]),
        "mu0_var_reduction": float(reductions["var_reduction_mu0"]),
        "k_EA_growth_var_reduction": float(reductions["var_reduction_k_EA_growth"]),
        "k_IAA_growth_var_reduction": float(reductions["var_reduction_k_IAA_growth"]),
        "k_EO_growth_var_reduction": float(reductions["var_reduction_k_EO_growth"]),
    }
    return summary, per_parameter


def selected_rows(design: str, names: list[str], ranking: pd.DataFrame) -> list[dict[str, object]]:
    info = ranking.set_index("candidate").to_dict(orient="index")
    rows = []
    for order, name in enumerate(names, start=1):
        item = info.get(name, {})
        rows.append(
            {
                "design": design,
                "campaign_order": order,
                "candidate": name,
                "family": item.get("family", ""),
                "horizon_h": item.get("horizon_h", np.nan),
                "temperature_c": item.get("temperature_c", ""),
            }
        )
    return rows


def write_report(
    path: Path,
    summary: pd.DataFrame,
    selected: pd.DataFrame,
    reductions: pd.DataFrame,
    exchange_log: pd.DataFrame,
) -> None:
    lines = [
        "# D-opt benchmark for joint fermentation-aroma DOE",
        "",
        "This benchmark reuses the successful candidate FIMs from the joint fermentation-aroma run.",
        "No new Ipopt/Pyomo DoE sensitivity solves are performed.",
        "",
        "Compared designs:",
        "",
        "- `hybrid_current`: selected with the hybrid score used in the current campaign.",
        "- `dopt_greedy`: selected from scratch with pure D-optimality, `logdet(F_post)`.",
        "- `dopt_exchange_from_hybrid`: one-swap local search initialized at `hybrid_current`, accepting only pure D-opt logdet gains.",
        "",
        "## Summary",
        "",
        summary.to_markdown(index=False),
        "",
        "## Selected campaigns",
        "",
        selected.to_markdown(index=False),
        "",
        "## Exchange log",
        "",
        exchange_log.to_markdown(index=False) if not exchange_log.empty else "No improving D-opt exchange was found.",
        "",
        "## Weakest reductions by design",
        "",
        reductions.sort_values(["design", "var_reduction"]).groupby("design").head(6).to_markdown(index=False),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Pure D-opt benchmark for the joint aroma DOE campaign.")
    parser.add_argument("--results-dir", default=str(RESULTS_DIR))
    parser.add_argument("--campaign-size", type=int, default=9)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    metadata = json.loads((results_dir / "aroma_campaign_metadata.json").read_text(encoding="utf-8"))
    parameters = tuple(str(p) for p in metadata["parameters"])
    prior_fim = load_matrix(results_dir / "aroma_prior_fim.csv", parameters)
    fims = load_candidate_fims(results_dir, parameters)
    if not fims:
        raise RuntimeError("No successful candidate FIMs found.")

    ranking = pd.read_csv(results_dir / "aroma_candidate_ranking.csv")
    current_selected = pd.read_csv(results_dir / "aroma_campaign_selected.csv")["candidate"].astype(str).tolist()
    current_selected = current_selected[: min(int(args.campaign_size), len(current_selected))]

    dopt_from_scratch = dopt_greedy(prior_fim, fims, parameters, min(int(args.campaign_size), len(fims)))
    dopt_from_seed, exchange_log = dopt_exchange(current_selected, prior_fim, fims, parameters)

    designs = {
        "hybrid_current": current_selected,
        "dopt_greedy": dopt_from_scratch,
        "dopt_exchange_from_hybrid": dopt_from_seed,
    }

    summary_rows = []
    parameter_rows = []
    selected_campaign_rows = []
    for design, names in designs.items():
        summary, per_parameter = design_summary(design, names, prior_fim, fims, parameters)
        summary_rows.append(summary)
        parameter_rows.extend(per_parameter)
        selected_campaign_rows.extend(selected_rows(design, names, ranking))

    summary_frame = pd.DataFrame(summary_rows)
    reduction_frame = pd.DataFrame(parameter_rows)
    selected_frame = pd.DataFrame(selected_campaign_rows)
    exchange_frame = pd.DataFrame(exchange_log)

    summary_frame.to_csv(results_dir / "aroma_dopt_benchmark_summary.csv", index=False)
    selected_frame.to_csv(results_dir / "aroma_dopt_benchmark_selected.csv", index=False)
    reduction_frame.to_csv(results_dir / "aroma_dopt_benchmark_parameter_reduction.csv", index=False)
    exchange_frame.to_csv(results_dir / "aroma_dopt_benchmark_exchange_log.csv", index=False)
    write_report(
        results_dir / "aroma_dopt_benchmark_report.md",
        summary_frame,
        selected_frame,
        reduction_frame,
        exchange_frame,
    )

    print(summary_frame.to_string(index=False))
    if exchange_frame.empty:
        print("\nNo improving D-opt exchange was found from the hybrid seed.")
    else:
        print("\nExchange log:")
        print(exchange_frame.to_string(index=False))


if __name__ == "__main__":
    main()
