from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_FILE = SCRIPT_DIR / "data" / "Calibration_data_vl3.xlsx"
RESULTS_DIR = SCRIPT_DIR / "results" / "medium_transfer_diagnostics"
FINAL_THETA_PATH = SCRIPT_DIR / "results" / "identifiability_reduction" / "theta_final_identifiable.csv"
ACCEPTED7 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF", "iG")

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_fit_strategy_analysis import complete_theta, load_notebook_context, reference_theta, series_to_theta


STATE_COLUMNS = {
    "X": "Viability",
    "N": "YAN",
    "G": "GLUCOSE",
    "F": "FRUCTOSE",
    "E": "ETANOL",
}

AROMA_COLUMNS = {
    "ethyl_acetate": "Ethyl_Acetate_total",
    "isoamyl_acetate": "isoamil_acetate_total",
    "ethyl_octanoate": "octanoate_de_etilo_total",
}


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def load_raw_batches(data_file: Path) -> dict[str, pd.DataFrame]:
    xls = pd.ExcelFile(data_file)
    batches = {}
    for sheet in xls.sheet_names:
        df = pd.read_excel(data_file, sheet_name=sheet)
        df.columns = [str(col).strip() for col in df.columns]
        if "t" in df.columns:
            df["t"] = safe_numeric(df["t"])
            df = df.dropna(subset=["t"]).sort_values("t")
        batches[str(sheet)] = df.reset_index(drop=True)
    return batches


def first_valid(series: pd.Series) -> float:
    values = safe_numeric(series).dropna()
    return float(values.iloc[0]) if not values.empty else math.nan


def last_valid(series: pd.Series) -> float:
    values = safe_numeric(series).dropna()
    return float(values.iloc[-1]) if not values.empty else math.nan


def batch_audit(raw_batches: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for batch, df in raw_batches.items():
        row: dict[str, object] = {
            "batch": batch,
            "n_rows": int(len(df)),
            "t_min_h": float(safe_numeric(df.get("t", pd.Series(dtype=float))).min()),
            "t_max_h": float(safe_numeric(df.get("t", pd.Series(dtype=float))).max()),
        }
        for state, column in STATE_COLUMNS.items():
            if column not in df.columns:
                row[f"{state}_initial"] = math.nan
                row[f"{state}_final"] = math.nan
                row[f"{state}_n_obs_raw"] = 0
                continue
            values = safe_numeric(df[column])
            row[f"{state}_initial"] = first_valid(values)
            row[f"{state}_final"] = last_valid(values)
            row[f"{state}_n_obs_raw"] = int(values.notna().sum())
        for column in ("temperatura", "pulso_nut", "PAN", "AMMONIA", "GLYCEROL", "PYRUVIC ACID", "ACETALDEHIDO"):
            if column in df.columns:
                values = safe_numeric(df[column])
                row[f"{column}_initial"] = first_valid(values)
                row[f"{column}_mean"] = float(values.mean()) if values.notna().any() else math.nan
                row[f"{column}_n_obs_raw"] = int(values.notna().sum())
        for species, column in AROMA_COLUMNS.items():
            values = safe_numeric(df[column]) if column in df.columns else pd.Series(dtype=float)
            row[f"{species}_initial"] = first_valid(values)
            row[f"{species}_final"] = last_valid(values)
            row[f"{species}_n_obs_raw"] = int(values.notna().sum())
        pulse_values = safe_numeric(df.get("pulso_nut", pd.Series(dtype=float))).dropna()
        pulse_values = pulse_values[pulse_values > 0.0]
        row["n_positive_nutrient_pulses"] = int(len(pulse_values))
        row["nutrient_pulse_total_mg_L"] = float(pulse_values.sum()) if not pulse_values.empty else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("batch").reset_index(drop=True)


def write_medium_template(audit: pd.DataFrame, path: Path) -> pd.DataFrame:
    template = audit[["batch", "n_rows", "t_max_h", "G_initial", "F_initial", "N_initial", "temperatura_mean"]].copy()
    template.insert(1, "medium_type", "unknown")
    template.insert(2, "medium_lot", "")
    template["notes"] = ""
    template.to_csv(path, index=False)
    return template


def load_medium_map(path: Path | None, batches: list[str]) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        frame = pd.DataFrame(
            [
                {
                    "batch": str(batch),
                    "medium_type": str(value.get("medium_type", value) if isinstance(value, dict) else value),
                    "medium_lot": str(value.get("medium_lot", "") if isinstance(value, dict) else ""),
                    "notes": str(value.get("notes", "") if isinstance(value, dict) else ""),
                }
                for batch, value in data.items()
            ]
        )
    else:
        frame = pd.read_csv(path)
    required = {"batch", "medium_type"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Medium map is missing required columns: {missing}")
    frame = frame.copy()
    frame["batch"] = frame["batch"].astype(str)
    frame["medium_type"] = frame["medium_type"].astype(str).str.strip().str.lower()
    if "medium_lot" not in frame.columns:
        frame["medium_lot"] = ""
    if "notes" not in frame.columns:
        frame["notes"] = ""
    unknown_batches = sorted(set(frame["batch"]) - set(batches))
    if unknown_batches:
        raise ValueError(f"Medium map contains batches not present in data file: {unknown_batches}")
    missing_batches = sorted(set(batches) - set(frame["batch"]))
    if missing_batches:
        filler = pd.DataFrame(
            {
                "batch": missing_batches,
                "medium_type": "unknown",
                "medium_lot": "",
                "notes": "not provided in medium map",
            }
        )
        frame = pd.concat([frame, filler], ignore_index=True)
    return frame.sort_values("batch").reset_index(drop=True)


def load_final_theta(ns: dict, batches: list[str]) -> dict[str, float]:
    reference = reference_theta(ns, "reduced6", batch_id=batches[0])
    if FINAL_THETA_PATH.exists():
        reference.update(series_to_theta(pd.read_csv(FINAL_THETA_PATH, index_col=0).iloc[:, 0]))
    clipped, _ = ns["clip_theta_to_bounds"](reference)
    return {name: float(clipped[name]) for name in ns["DEFAULT_THETA"]}


def objective_by_batch(ns: dict, theta: dict[str, float], batches: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    objective_rows = []
    residual_frames = []
    solve_frames = []
    for batch in batches:
        summary, residuals, solve_summary = ns["evaluate_fit_objectives"](theta, batch_ids=[batch])
        row = {"batch": str(batch), **summary}
        row["WSSE_per_observation"] = (
            float(row["WSSE_raw_sum"]) / float(row["n_observations"]) if float(row["n_observations"]) > 0.0 else math.nan
        )
        objective_rows.append(row)
        if residuals is not None and not residuals.empty:
            residual_frames.append(residuals.assign(batch=str(batch)))
        if solve_summary is not None and not solve_summary.empty:
            solve_frames.append(solve_summary.assign(batch=str(batch)))
    residuals_all = pd.concat(residual_frames, ignore_index=True) if residual_frames else pd.DataFrame()
    solve_all = pd.concat(solve_frames, ignore_index=True) if solve_frames else pd.DataFrame()
    return pd.DataFrame(objective_rows), residuals_all, solve_all


def summarize_residuals(residuals: pd.DataFrame, medium_map: pd.DataFrame) -> pd.DataFrame:
    if residuals.empty:
        return pd.DataFrame()
    merged = residuals.merge(medium_map[["batch", "medium_type", "medium_lot"]], on="batch", how="left")
    grouped = (
        merged.groupby(["medium_type", "state"], dropna=False)
        .agg(
            n=("weighted_squared_error", "size"),
            WSSE=("weighted_squared_error", "sum"),
            RMSE=("residual", lambda x: float(np.sqrt(np.mean(np.square(x))))),
            mean_signed_residual=("residual", "mean"),
            mean_abs_residual=("residual", lambda x: float(np.mean(np.abs(x)))),
        )
        .reset_index()
    )
    grouped["WSSE_per_observation"] = grouped["WSSE"] / grouped["n"].replace(0, np.nan)
    return grouped


def summarize_objective_by_medium(objective: pd.DataFrame, medium_map: pd.DataFrame) -> pd.DataFrame:
    merged = objective.merge(medium_map[["batch", "medium_type", "medium_lot"]], on="batch", how="left")
    grouped = (
        merged.groupby("medium_type", dropna=False)
        .agg(
            n_batches=("batch", "nunique"),
            n_observations=("n_observations", "sum"),
            WSSE_raw_sum=("WSSE_raw_sum", "sum"),
            SSE_raw_sum=("SSE_raw_sum", "sum"),
        )
        .reset_index()
    )
    grouped["WSSE_per_observation"] = grouped["WSSE_raw_sum"] / grouped["n_observations"].replace(0, np.nan)
    return grouped


def run_medium_fims(
    ns: dict,
    theta: dict[str, float],
    medium_map: pd.DataFrame,
    parameters: tuple[str, ...],
    fraction: float,
    out_dir: Path,
) -> pd.DataFrame:
    rows = []
    for label, group in [("pooled", medium_map)] + [
        (medium, medium_map[medium_map["medium_type"].eq(medium)])
        for medium in sorted(set(medium_map["medium_type"]) - {"unknown"})
    ]:
        batches = group["batch"].astype(str).tolist()
        if not batches:
            continue
        q_raw, q_weighted, q_relative, perturbation_summary, solve_summary = ns[
            "build_weighted_relative_sensitivity_matrix"
        ](theta, list(parameters), batches, fraction=fraction)
        fim, eigen_summary, loading_summary = ns["eigen_analysis_from_q"](q_relative, list(parameters))
        prefix = out_dir / f"fim_{label}"
        q_relative.to_csv(prefix.with_name(f"{prefix.name}_q_weighted_relative.csv"))
        fim.to_csv(prefix.with_name(f"{prefix.name}_matrix.csv"))
        eigen_summary.to_csv(prefix.with_name(f"{prefix.name}_eigenvalues.csv"), index=False)
        loading_summary.to_csv(prefix.with_name(f"{prefix.name}_eigendirections.csv"), index=False)
        perturbation_summary.to_csv(prefix.with_name(f"{prefix.name}_perturbation_summary.csv"))
        solve_summary.to_csv(prefix.with_name(f"{prefix.name}_solve_summary.csv"), index=False)
        rows.append(
            {
                "fim_group": label,
                "batches": ", ".join(batches),
                "n_batches": len(batches),
                "n_parameters": len(parameters),
                "condition_number": float(eigen_summary["condition_number"].iloc[0]),
                "min_relative_eigenvalue": float(eigen_summary["relative_eigenvalue"].min()),
                "near_null_directions": int(eigen_summary["near_null"].sum()),
                "weakest_direction": ""
                if loading_summary.empty
                else str(loading_summary.iloc[0].get("dominant_parameters", "")),
            }
        )
    return pd.DataFrame(rows)


def run_transfer_refits(
    ns: dict,
    medium_map: pd.DataFrame,
    parameters: tuple[str, ...],
    initial_theta: dict[str, float],
    out_dir: Path,
) -> pd.DataFrame:
    from run_fit_strategy_analysis import FitStrategy, run_parmest_strategy

    rows = []
    media = [m for m in sorted(medium_map["medium_type"].dropna().unique()) if m != "unknown"]
    groups = {"pooled": medium_map["batch"].astype(str).tolist()}
    groups.update({medium: medium_map.loc[medium_map["medium_type"].eq(medium), "batch"].astype(str).tolist() for medium in media})
    for fit_group, train_batches in groups.items():
        if not train_batches:
            continue
        strategy = FitStrategy(
            name=f"medium_{fit_group}",
            estimated_parameters=parameters,
            initial_source="reduced6",
        )
        obj_value, theta_fit, _estimator, penalty = run_parmest_strategy(
            ns,
            strategy,
            batches=train_batches,
            initial_theta=initial_theta,
            tee=False,
        )
        pd.Series(theta_fit, name=fit_group).to_csv(out_dir / f"theta_fit_{fit_group}.csv")
        for eval_group, eval_batches in groups.items():
            summary, _residuals, solve_summary = ns["evaluate_fit_objectives"](theta_fit, batch_ids=eval_batches)
            rows.append(
                {
                    "fit_group": fit_group,
                    "eval_group": eval_group,
                    "train_batches": ", ".join(train_batches),
                    "eval_batches": ", ".join(eval_batches),
                    "n_observations": int(summary["n_observations"]),
                    "WSSE_raw_sum": float(summary["WSSE_raw_sum"]),
                    "WSSE_per_observation": float(summary["WSSE_raw_sum"]) / max(float(summary["n_observations"]), 1.0),
                    "parmest_train_objective": float(obj_value),
                    "l2_penalty": float(penalty),
                    "failed_solves": int(
                        0
                        if solve_summary.empty
                        else (~solve_summary["termination"].astype(str).str.contains("optimal", case=False, na=False)).sum()
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_report(
    out_dir: Path,
    audit: pd.DataFrame,
    medium_map: pd.DataFrame | None,
    objective_medium: pd.DataFrame | None,
    residual_summary: pd.DataFrame | None,
    fim_summary: pd.DataFrame | None,
    transfer_summary: pd.DataFrame | None,
    template_path: Path,
) -> None:
    lines = [
        "# Medium transfer diagnostics",
        "",
        "This report separates the historical fermentation database by must medium when metadata are available.",
        "The goal is to decide whether synthetic-must experiments can be used directly for wine-must parameter estimation, or whether natural-must data require medium-specific correction terms.",
        "",
        "## Data audit",
        "",
        audit[
            [
                "batch",
                "n_rows",
                "t_max_h",
                "G_initial",
                "F_initial",
                "N_initial",
                "E_initial",
                "temperatura_mean",
                "n_positive_nutrient_pulses",
                "isoamyl_acetate_n_obs_raw",
                "ethyl_octanoate_n_obs_raw",
                "ethyl_acetate_n_obs_raw",
            ]
        ].to_markdown(index=False),
        "",
    ]
    if medium_map is None or medium_map["medium_type"].eq("unknown").all():
        lines.extend(
            [
                "## Medium metadata status",
                "",
                "No usable synthetic/natural medium mapping was provided or found in the Excel file.",
                f"Fill `{template_path.name}` with `medium_type = synthetic` or `natural`, then rerun this script with `--medium-map {template_path}`.",
                "",
                "Until that mapping exists, the current DOE can use the historical data as a pooled prior, but it cannot distinguish model-structure limitations from medium-transfer limitations.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## Medium map",
                "",
                medium_map.to_markdown(index=False),
                "",
            ]
        )
        if objective_medium is not None and not objective_medium.empty:
            lines.extend(["## Current-model objective by medium", "", objective_medium.to_markdown(index=False), ""])
        if residual_summary is not None and not residual_summary.empty:
            lines.extend(["## Current-model residuals by medium and state", "", residual_summary.to_markdown(index=False), ""])
        if fim_summary is not None and not fim_summary.empty:
            lines.extend(["## FIM/eigen diagnostics by medium", "", fim_summary.to_markdown(index=False), ""])
        if transfer_summary is not None and not transfer_summary.empty:
            lines.extend(["## Transfer refit diagnostics", "", transfer_summary.to_markdown(index=False), ""])
        lines.extend(
            [
                "## Deterministic interpretation criteria",
                "",
                "1. If pooled, synthetic-only, and natural-only fits have similar WSSE per observation and residual sign patterns, keep a shared kinetic parameter vector.",
                "2. If one medium has systematic residual bias in X, N, G, F, or E, introduce medium-specific initial-composition or matrix correction terms before adding new biological parameters.",
                "3. If FIM eigenvectors differ by medium, use the medium that excites the weak directions in the next DOE batch.",
                "4. If synthetic-fit to natural-validation degrades strongly while natural-fit to synthetic-validation does not, prioritize natural must in the first campaign batch.",
                "5. If both transfer directions degrade, use a mixed campaign and estimate medium offsets with shrinkage instead of treating all data as one homogeneous population.",
                "",
            ]
        )
    (out_dir / "medium_transfer_diagnostics_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-file", type=Path, default=DATA_FILE)
    parser.add_argument("--medium-map", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--parameters", nargs="*", default=list(ACCEPTED7))
    parser.add_argument("--fim-fraction", type=float, default=0.01)
    parser.add_argument("--run-fim", action="store_true")
    parser.add_argument("--run-transfer-refits", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_batches = load_raw_batches(args.data_file)
    batches = sorted(raw_batches)
    audit = batch_audit(raw_batches)
    audit.to_csv(args.out_dir / "historical_batch_audit.csv", index=False)

    template_path = args.out_dir / "medium_map_template.csv"
    write_medium_template(audit, template_path)
    medium_map = load_medium_map(args.medium_map, batches)
    if medium_map is not None:
        medium_map.to_csv(args.out_dir / "medium_map_used.csv", index=False)

    objective = residuals = solve_summary = pd.DataFrame()
    objective_medium = residual_summary = fim_summary = transfer_summary = None

    if medium_map is not None and not medium_map["medium_type"].eq("unknown").all():
        ns = load_notebook_context()
        theta = load_final_theta(ns, batches)
        objective, residuals, solve_summary = objective_by_batch(ns, theta, batches)
        objective = objective.merge(medium_map[["batch", "medium_type", "medium_lot"]], on="batch", how="left")
        objective.to_csv(args.out_dir / "current_model_objective_by_batch.csv", index=False)
        residuals.to_csv(args.out_dir / "current_model_residuals.csv", index=False)
        solve_summary.to_csv(args.out_dir / "current_model_solve_summary.csv", index=False)
        objective_medium = summarize_objective_by_medium(objective, medium_map)
        residual_summary = summarize_residuals(residuals, medium_map)
        objective_medium.to_csv(args.out_dir / "current_model_objective_by_medium.csv", index=False)
        residual_summary.to_csv(args.out_dir / "current_model_residuals_by_medium_state.csv", index=False)

        parameters = tuple(str(name) for name in args.parameters)
        if args.run_fim:
            fim_summary = run_medium_fims(ns, theta, medium_map, parameters, args.fim_fraction, args.out_dir)
            fim_summary.to_csv(args.out_dir / "medium_fim_summary.csv", index=False)
        if args.run_transfer_refits:
            transfer_summary = run_transfer_refits(ns, medium_map, parameters, theta, args.out_dir)
            transfer_summary.to_csv(args.out_dir / "medium_transfer_refit_summary.csv", index=False)

    write_report(
        args.out_dir,
        audit,
        medium_map,
        objective_medium,
        residual_summary,
        fim_summary,
        transfer_summary,
        template_path,
    )
    print(f"Wrote medium transfer diagnostics to {args.out_dir}")
    print(f"Batch audit: {args.out_dir / 'historical_batch_audit.csv'}")
    if medium_map is None:
        print(f"Fill medium map template: {template_path}")


if __name__ == "__main__":
    main()
