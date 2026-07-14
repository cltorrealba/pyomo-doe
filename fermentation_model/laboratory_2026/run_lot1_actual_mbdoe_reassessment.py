from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SCRIPT_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from laboratory_2026 import run_final_operational_doe_v2 as final
from laboratory_2026 import run_final_operational_doe_volume_constrained as vol
from shared import run_new_must_glycerol_estimability_doe as base
from shared import run_secondary_joint_campaign_doe as joint
from shared.paths import LABORATORY_2026_RESULTS_DIR


RESULTS_DIR = LABORATORY_2026_RESULTS_DIR / "lot1_actual_mbdoe_reassessment"
LOT1_DIR = LABORATORY_2026_RESULTS_DIR / "lot1_data_preview" / "processed"

POLICY = "full14_small6"
EXECUTED_NAMES = (
    "synthetic_lit_SM410_18C_highN_ester",
    "synthetic_high_biomass_low_N_maintenance",
    "synthetic_fructose_rich_glucose_pulse",
)
ORIGINAL_REMAINING = (
    "synthetic_glucose_rich_fructose_pulse",
    "synthetic_lit_SM410_24C_highN_strip",
    "synthetic_viable_biomass_step",
    "synthetic_high_sugar_reference",
    "natural_glucose_pulse_after_growth",
    "synthetic_fast_CO2_aroma_strip_highN",
)
LOT2_FIXED = (
    "synthetic_glucose_rich_fructose_pulse",
    "synthetic_lit_SM410_24C_highN_strip",
    "synthetic_viable_biomass_step",
)
LOT3_ORIGINAL = (
    "synthetic_high_sugar_reference",
    "natural_glucose_pulse_after_growth",
    "synthetic_fast_CO2_aroma_strip_highN",
)
LOT3_DEATH_REPLACEMENT = (
    "synthetic_late_ethanol_death_probe",
    "natural_glucose_pulse_after_growth",
    "synthetic_fast_CO2_aroma_strip_highN",
)
LOT3_NATURAL_REPLACEMENT = (
    "natural_aroma_matrix_cold_hot",
    "natural_glucose_pulse_after_growth",
    "synthetic_fast_CO2_aroma_strip_highN",
)


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def load_lot1_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y15 = pd.read_csv(LOT1_DIR / "lot1_y15_wide_processed.csv")
    oculyze = pd.read_csv(LOT1_DIR / "lot1_oculyze_processed.csv")
    inputs = pd.read_csv(LOT1_DIR / "lot1_actual_operational_inputs_loaded.csv")
    co2 = pd.read_csv(LOT1_DIR / "lot1_co2_filt_volume_corrected_10min.csv")
    for df in (y15, oculyze, inputs, co2):
        for column in df.columns:
            if column.endswith("_h") or column in {
                "actual_t_h",
                "actual_plot_t_h",
                "actual_amount",
                "sugar_total_g_l",
                "glucose_g_l",
                "fructose_g_l",
                "yan_for_model_physical_mg_l",
                "ethanol_g_l_for_model_manual",
                "glycerol_g_l",
                "pyruvic_acid_mg_l",
                "acetaldehyde_mg_l",
                "acetic_acid_g_l",
                "x_viable_g_l",
                "x_dead_g_l",
                "flow_filt_sccm",
            }:
                df[column] = _num(df[column])
    return y15, oculyze, inputs, co2


def first_row_at_zero(df: pd.DataFrame, process: str) -> pd.Series:
    sub = df[df["process"].eq(process)].copy()
    sub["abs_t0"] = (sub["t_h"] - 0.0).abs()
    return sub.sort_values("abs_t0").iloc[0]


def actual_design_from_process(
    process: str,
    nominal: base.FutureDesign,
    y15: pd.DataFrame,
    oculyze: pd.DataFrame,
    inputs: pd.DataFrame,
) -> base.FutureDesign:
    y0 = first_row_at_zero(y15, process)
    x0 = first_row_at_zero(oculyze, process)
    initials = dict(nominal.initials)
    initials.update(
        {
            "X": float(x0.get("x_viable_g_l", nominal.initials.get("X", 0.4))),
            "Xd": float(x0.get("x_dead_g_l", nominal.initials.get("Xd", 0.0))),
            "N": float(y0.get("yan_for_model_physical_mg_l", nominal.initials.get("N", 0.2) * 1000.0)) / 1000.0,
            "G": float(y0.get("glucose_g_l", nominal.initials.get("G", 100.0))),
            "F": float(y0.get("fructose_g_l", nominal.initials.get("F", 100.0))),
            "E": float(y0.get("ethanol_g_l_for_model_manual", nominal.initials.get("E", 0.0))),
            "Gly": float(y0.get("glycerol_g_l", nominal.initials.get("Gly", 0.0))),
            "Pyr": float(y0.get("pyruvic_acid_mg_l", nominal.initials.get("Pyr", 0.0))),
            "AcAld": float(y0.get("acetaldehyde_mg_l", nominal.initials.get("AcAld", 0.0))),
            "Acetate": float(y0.get("acetic_acid_g_l", nominal.initials.get("Acetate", 0.0))),
        }
    )
    pulses = {channel: [] for channel in base.INPUT_CHANNELS}
    sub = inputs[inputs["process"].eq(process)].copy()
    for row in sub.itertuples(index=False):
        if str(getattr(row, "executed_bool", getattr(row, "executed", ""))).lower() not in {"true", "1", "yes"}:
            continue
        channel = str(getattr(row, "channel")).strip().upper()
        if channel not in pulses:
            continue
        time_h = float(getattr(row, "actual_plot_t_h", np.nan))
        amount = float(getattr(row, "actual_amount_for_plot", getattr(row, "actual_amount", np.nan)))
        if np.isfinite(time_h) and np.isfinite(amount):
            pulses[channel].append((time_h, amount))
    return base.FutureDesign(
        name=f"actual_{process}_{nominal.name}",
        family=f"actual_lot1_{nominal.family}",
        medium=nominal.medium,
        horizon_h=float(max(y15[y15["process"].eq(process)]["t_h"].max(), nominal.horizon_h)),
        initials=initials,
        temperature_segments=nominal.temperature_segments,
        pulses={channel: tuple(sorted(rows)) for channel, rows in pulses.items()},
        rationale=f"Lot 1 executed version of {nominal.name}.",
    )


def process_times(y15: pd.DataFrame, oculyze: pd.DataFrame, co2: pd.DataFrame, process: str) -> dict[str, np.ndarray]:
    y = y15[y15["process"].eq(process)].sort_values("t_h")
    x = oculyze[oculyze["process"].eq(process)].sort_values("t_h")
    liquid = np.asarray(sorted(set(y["t_h"].dropna().round(8)).union(set(x["t_h"].dropna().round(8)))), dtype=float)
    if "ethanol_status" in y.columns:
        full = y[y["ethanol_status"].astype(str).eq("measured")]["t_h"].dropna().to_numpy(dtype=float)
    else:
        full = y["t_h"].dropna().to_numpy(dtype=float)
    co2_times = co2[co2["process"].eq(process)]["t_h"].dropna().to_numpy(dtype=float)
    if len(co2_times):
        horizon = float(np.nanmax(co2_times))
        ctimes = np.arange(0.0, horizon + 1e-9, 6.0, dtype=float)
    else:
        ctimes = np.array([], dtype=float)
    return {
        "liquid": np.asarray(sorted(set(liquid)), dtype=float),
        "full": np.asarray(sorted(set(np.round(full, 8))), dtype=float),
        "co2": ctimes,
    }


def actual_residual_vector(
    theta: dict[str, float],
    design: base.FutureDesign,
    times: dict[str, np.ndarray],
    include_aroma_if_available: bool,
) -> np.ndarray:
    all_times = set(np.arange(0.0, float(design.horizon_h) + 1e-9, 2.0).round(8))
    for key in ("liquid", "full", "co2"):
        all_times.update(float(t) for t in times.get(key, np.array([], dtype=float)))
    for schedule in design.pulses.values():
        for time_h, _amount in schedule:
            all_times.add(round(max(0.0, float(time_h) - 2.0), 8))
            all_times.add(round(float(time_h), 8))
            all_times.add(round(min(float(design.horizon_h), float(time_h) + 2.0), 8))
    grid = np.asarray(sorted(t for t in all_times if 0.0 <= t <= float(design.horizon_h)), dtype=float)

    original = final.dense_time_grid
    try:
        final.dense_time_grid = lambda _d, _p: grid
        sim = final.simulate_v2(design, theta, POLICY)
    finally:
        final.dense_time_grid = original
    if sim is None:
        return np.ones(1000, dtype=float) * 1e6

    residuals: list[np.ndarray] = []
    liquid = [float(t) for t in times["liquid"] if float(t) in sim.core.index]
    full = [float(t) for t in times["full"] if float(t) in sim.core.index]
    for state in base.STATE_NAMES:
        state_times = full if state == "E" else liquid
        if not state_times:
            continue
        center = sim.core.loc[state_times, state].to_numpy(dtype=float)
        residuals.append(center / base.sigma_for_state(state, center))

    for state in ("Pyr", "AcAld", "Acetate"):
        state_times = [float(t) for t in liquid if float(t) in sim.secondary.index]
        if state_times:
            center = sim.secondary.loc[state_times, state].to_numpy(dtype=float)
            residuals.append(center / joint.SIGMA[state])

    ctimes = [float(t) for t in times["co2"] if float(t) in sim.secondary.index]
    if ctimes:
        residuals.append(sim.secondary.loc[ctimes, "CO2"].to_numpy(dtype=float) / joint.SIGMA["CO2"])

    if include_aroma_if_available and full:
        liq_aroma = sim.aromas_liq.loc[full]
        for species in joint.AROMA_SPECIES:
            residuals.append(liq_aroma[species].to_numpy(dtype=float) / joint.SIGMA[f"{species}_liq"])
            residuals.append(np.array([sim.aromas_cond[species] / joint.SIGMA[f"{species}_cond"]], dtype=float))

    if not residuals:
        return np.ones(1000, dtype=float) * 1e6
    return np.concatenate(residuals)


def add_metrics(row: dict, prior: np.ndarray, fim: np.ndarray) -> dict:
    combined = prior + fim
    return {
        **row,
        **joint.fim_metrics(combined, final.TARGET_PARAMETERS, prefix="campaign_"),
        **joint.variance_reduction(prior, combined, final.TARGET_PARAMETERS),
        "hybrid_score": joint.score_fim(combined, final.TARGET_PARAMETERS, "hybrid"),
        "dopt_score": joint.score_fim(combined, final.TARGET_PARAMETERS, "d_opt"),
    }


def fim_sum(fims: dict[str, np.ndarray], names: tuple[str, ...], shape: tuple[int, int]) -> np.ndarray:
    out = np.zeros(shape, dtype=float)
    for name in names:
        if name not in fims:
            raise KeyError(f"Missing FIM for {name}")
        out = out + fims[name]
    return out


def scenario_frame(
    scenarios: dict[str, tuple[str, ...]],
    prior: np.ndarray,
    fims: dict[str, np.ndarray],
) -> pd.DataFrame:
    rows = []
    shape = prior.shape
    for label, names in scenarios.items():
        fim = fim_sum(fims, names, shape)
        row = add_metrics(
            {
                "scenario": label,
                "n_experiments": len(names),
                "selected_candidates": ", ".join(names),
            },
            prior,
            fim,
        )
        rows.append(row)
    return pd.DataFrame(rows)


def weakest_table(fim: np.ndarray, label: str, prior: np.ndarray, n: int = 12) -> pd.DataFrame:
    row = joint.variance_reduction(prior, fim, final.TARGET_PARAMETERS)
    rows = []
    for name in final.TARGET_PARAMETERS:
        rows.append(
            {
                "analysis": label,
                "parameter": name,
                "var_ratio": row.get(f"var_ratio_{name}", np.nan),
                "var_reduction": row.get(f"var_reduction_{name}", np.nan),
            }
        )
    return pd.DataFrame(rows).sort_values("var_reduction", ascending=True).head(n)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=float, default=0.04)
    parser.add_argument("--include-aroma-in-actual-lot1", action="store_true")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    theta = final.load_theta_final()
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    prior, prior_estimability, prior_eigen = final.current_prior_fim(theta, batches, args.step)
    prior_estimability.to_csv(RESULTS_DIR / "current_prior_estimability.csv", index=False)
    prior_eigen.to_csv(RESULTS_DIR / "current_prior_eigen_diagnostics.csv", index=False)

    designs = final.candidate_library(data)
    y15, oculyze, inputs, co2 = load_lot1_tables()

    actual_designs = {}
    process_map = {
        "F1": "synthetic_lit_SM410_18C_highN_ester",
        "F2": "synthetic_high_biomass_low_N_maintenance",
        "F3": "synthetic_fructose_rich_glucose_pulse",
    }
    actual_times = {}
    for process, name in process_map.items():
        actual_designs[process] = actual_design_from_process(process, designs[name], y15, oculyze, inputs)
        actual_times[process] = process_times(y15, oculyze, co2, process)

    rows = []
    planned_fims = {}
    for name in EXECUTED_NAMES:
        print(f"[planned_lot1] {name}", flush=True)
        fim = final.finite_difference_fim(
            theta,
            final.TARGET_PARAMETERS,
            lambda th, d=designs[name]: vol.future_residual_vector(th, d, POLICY),
            args.step,
        )
        planned_fims[name] = fim
    planned_lot1 = fim_sum(planned_fims, EXECUTED_NAMES, prior.shape)
    rows.append(add_metrics({"basis": "planned_lot1", "description": "Nominal Lot 1 F1-F3, volume-constrained residuals."}, prior, planned_lot1))

    actual_fims = {}
    for process, design in actual_designs.items():
        print(f"[actual_lot1] {process} {design.name}", flush=True)
        fim = final.finite_difference_fim(
            theta,
            final.TARGET_PARAMETERS,
            lambda th, d=design, ts=actual_times[process]: actual_residual_vector(
                th,
                d,
                ts,
                include_aroma_if_available=bool(args.include_aroma_in_actual_lot1),
            ),
            args.step,
        )
        actual_fims[process] = fim
    actual_lot1 = sum(actual_fims.values(), np.zeros_like(prior))
    rows.append(
        add_metrics(
            {
                "basis": "actual_lot1_observed_channels",
                "description": "Executed Lot 1 with actual initials, pulses, liquid times, CO2; aroma excluded unless flag is used.",
            },
            prior,
            actual_lot1,
        )
    )
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "lot1_planned_vs_actual_fim_metrics.csv", index=False)

    actual_prior = prior + actual_lot1
    candidate_fims = {}
    candidate_rows = []
    remaining_designs = {name: design for name, design in designs.items() if name not in EXECUTED_NAMES}
    for idx, (name, design) in enumerate(remaining_designs.items(), start=1):
        print(f"[candidate_after_actual_lot1] {idx}/{len(remaining_designs)} {name}", flush=True)
        fim = final.finite_difference_fim(
            theta,
            final.TARGET_PARAMETERS,
            lambda th, d=design: vol.future_residual_vector(th, d, POLICY),
            args.step,
        )
        candidate_fims[name] = fim
        combined = actual_prior + fim
        candidate_rows.append(
            {
                "candidate": name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "hybrid_score": joint.score_fim(combined, final.TARGET_PARAMETERS, "hybrid"),
                "dopt_score": joint.score_fim(combined, final.TARGET_PARAMETERS, "d_opt"),
                **joint.fim_metrics(fim, final.TARGET_PARAMETERS, prefix="new_"),
                **joint.fim_metrics(combined, final.TARGET_PARAMETERS, prefix="combined_"),
                **joint.variance_reduction(actual_prior, combined, final.TARGET_PARAMETERS),
            }
        )
    ranking = pd.DataFrame(candidate_rows).sort_values("hybrid_score", ascending=False)
    ranking.to_csv(RESULTS_DIR / "candidate_ranking_after_actual_lot1.csv", index=False)

    selected_hybrid, fim_hybrid = joint.greedy_select(
        candidate_fims,
        remaining_designs,
        actual_prior,
        final.TARGET_PARAMETERS,
        6,
        "hybrid",
    )
    selected_hybrid.insert(0, "basis", "actual_lot1_free_greedy6")
    selected_hybrid.to_csv(RESULTS_DIR / "selected_after_actual_lot1_hybrid_free6.csv", index=False)

    selected_dopt, fim_dopt = joint.greedy_select(
        candidate_fims,
        remaining_designs,
        actual_prior,
        final.TARGET_PARAMETERS,
        6,
        "d_opt",
    )
    selected_dopt.insert(0, "basis", "actual_lot1_free_dopt6")
    selected_dopt.to_csv(RESULTS_DIR / "selected_after_actual_lot1_dopt_free6.csv", index=False)

    original_seed = [name for name in ORIGINAL_REMAINING if name in candidate_fims]
    greedy_dopt_seed = selected_dopt["candidate"].astype(str).tolist()
    exchange_from_original = joint.dopt_exchange(
        original_seed,
        candidate_fims,
        actual_prior,
        final.TARGET_PARAMETERS,
    )
    exchange_from_greedy = joint.dopt_exchange(
        greedy_dopt_seed,
        candidate_fims,
        actual_prior,
        final.TARGET_PARAMETERS,
    )
    frame_exchange_original, _fim_exchange_original = joint.selection_frame(
        exchange_from_original,
        remaining_designs,
        actual_prior,
        candidate_fims,
        final.TARGET_PARAMETERS,
        "d_opt_exchange_from_current",
    )
    frame_exchange_original.insert(0, "basis", "actual_lot1_dopt_exchange_from_current")
    frame_exchange_original.to_csv(RESULTS_DIR / "selected_after_actual_lot1_dopt_exchange_from_current.csv", index=False)
    frame_exchange_greedy, _fim_exchange_greedy = joint.selection_frame(
        exchange_from_greedy,
        remaining_designs,
        actual_prior,
        candidate_fims,
        final.TARGET_PARAMETERS,
        "d_opt_exchange_from_greedy",
    )
    frame_exchange_greedy.insert(0, "basis", "actual_lot1_dopt_exchange_from_greedy")
    frame_exchange_greedy.to_csv(RESULTS_DIR / "selected_after_actual_lot1_dopt_exchange_from_greedy.csv", index=False)

    lot2_fim = fim_sum(candidate_fims, LOT2_FIXED, prior.shape)
    lot2_prior = actual_prior + lot2_fim
    lot3_pool = {
        name: fim
        for name, fim in candidate_fims.items()
        if name not in set(LOT2_FIXED)
    }
    lot3_designs = {name: remaining_designs[name] for name in lot3_pool}
    selected_lot3, fim_lot3 = joint.greedy_select(
        lot3_pool,
        lot3_designs,
        lot2_prior,
        final.TARGET_PARAMETERS,
        3,
        "hybrid",
    )
    selected_lot3.insert(0, "basis", "actual_lot1_lot2_fixed_greedy_lot3")
    selected_lot3.to_csv(RESULTS_DIR / "selected_lot3_after_actual_lot1_and_fixed_lot2.csv", index=False)

    selected_lot3_dopt, fim_lot3_dopt = joint.greedy_select(
        lot3_pool,
        lot3_designs,
        lot2_prior,
        final.TARGET_PARAMETERS,
        3,
        "d_opt",
    )
    selected_lot3_dopt.insert(0, "basis", "actual_lot1_lot2_fixed_dopt_lot3")
    selected_lot3_dopt.to_csv(RESULTS_DIR / "selected_lot3_dopt_after_actual_lot1_and_fixed_lot2.csv", index=False)

    scenarios = {
        "original_remaining_F4_F9": ORIGINAL_REMAINING,
        "lot2_fixed_plus_original_lot3": LOT2_FIXED + LOT3_ORIGINAL,
        "lot2_fixed_plus_death_replacement": LOT2_FIXED + LOT3_DEATH_REPLACEMENT,
        "lot2_fixed_plus_natural_replacement": LOT2_FIXED + LOT3_NATURAL_REPLACEMENT,
        "free_greedy6_after_actual_lot1": tuple(selected_hybrid["candidate"].astype(str).tolist()),
        "free_dopt6_after_actual_lot1": tuple(selected_dopt["candidate"].astype(str).tolist()),
        "dopt_exchange_from_current": tuple(exchange_from_original),
        "dopt_exchange_from_greedy": tuple(exchange_from_greedy),
        "lot2_fixed_plus_greedy_lot3": LOT2_FIXED + tuple(selected_lot3["candidate"].astype(str).tolist()),
        "lot2_fixed_plus_dopt_lot3": LOT2_FIXED + tuple(selected_lot3_dopt["candidate"].astype(str).tolist()),
    }
    scenario_metrics = scenario_frame(scenarios, actual_prior, candidate_fims)
    scenario_metrics.to_csv(RESULTS_DIR / "scenario_metrics_after_actual_lot1.csv", index=False)

    weakest_after_actual = weakest_table(actual_prior, "after_actual_lot1", prior)
    weakest_after_original = weakest_table(actual_prior + fim_sum(candidate_fims, ORIGINAL_REMAINING, prior.shape), "after_original_remaining", prior)
    weakest_after_lot2 = weakest_table(actual_prior + lot2_fim, "after_actual_lot1_plus_fixed_lot2", prior)
    pd.concat([weakest_after_actual, weakest_after_lot2, weakest_after_original], ignore_index=True).to_csv(
        RESULTS_DIR / "weakest_parameters_key_stages.csv",
        index=False,
    )

    summary = {
        "include_aroma_in_actual_lot1": bool(args.include_aroma_in_actual_lot1),
        "actual_lot1_logdet": float(joint.fim_metrics(actual_prior, final.TARGET_PARAMETERS)["logdet"]),
        "planned_lot1_logdet": float(joint.fim_metrics(prior + planned_lot1, final.TARGET_PARAMETERS)["logdet"]),
        "best_free_hybrid6": tuple(selected_hybrid["candidate"].astype(str).tolist()),
        "best_free_dopt6": tuple(selected_dopt["candidate"].astype(str).tolist()),
        "dopt_exchange_from_current": tuple(exchange_from_original),
        "dopt_exchange_from_greedy": tuple(exchange_from_greedy),
        "best_lot3_after_fixed_lot2": tuple(selected_lot3["candidate"].astype(str).tolist()),
        "best_lot3_dopt_after_fixed_lot2": tuple(selected_lot3_dopt["candidate"].astype(str).tolist()),
    }
    pd.Series(summary).to_json(RESULTS_DIR / "summary.json", indent=2)
    print("[done] Lot 1 actual MBDoE reassessment complete", flush=True)


if __name__ == "__main__":
    main()
