from __future__ import annotations

import argparse
import math
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_final_operational_doe_v2 as final
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint


RESULTS_DIR = SCRIPT_DIR / "results" / "final_operational_doe_volume_constrained"

TARGET_PARAMETERS = final.TARGET_PARAMETERS
SAMPLE_VOLUME_ML = {
    "full": 50.0,  # ethanol 40 mL + aroma 10 mL; other aliquots should be split from this if possible.
    "small": 8.0,  # biomass/YAN/sugars/glycerol/secondary acids without ethanol/aroma.
}
INITIAL_REACTOR_VOLUME_ML = 2000.0
MIN_FINAL_VOLUME_ML = 1100.0

VOLUME_POLICIES = {
    "full10_small10": {"max_full_216": 10, "max_full_168": 8, "max_small_216": 10, "max_small_168": 8},
    "full12_small8": {"max_full_216": 12, "max_full_168": 10, "max_small_216": 8, "max_small_168": 6},
    "full14_small6": {"max_full_216": 14, "max_full_168": 12, "max_small_216": 6, "max_small_168": 4},
}

# Working stock concentrations for feasibility calculations.
# Units are chosen to match pulse amounts used in the model.
STOCK_CONCENTRATIONS = {
    "N": {"concentration": 20.0, "unit": "mg_N_per_mL", "max_single_pulse_ml": 15.0},
    "G": {"concentration": 700.0, "unit": "g_per_L", "max_single_pulse_ml": 125.0},
    "F": {"concentration": 700.0, "unit": "g_per_L", "max_single_pulse_ml": 125.0},
    "E": {"concentration": 789.0, "unit": "g_per_L", "max_single_pulse_ml": 80.0},
    "X": {"concentration": 100.0, "unit": "g_dry_biomass_per_L", "max_single_pulse_ml": 35.0},
}


def _max_for_horizon(policy: str, horizon_h: float, kind: str) -> int:
    cfg = VOLUME_POLICIES[policy]
    key = f"max_{kind}_{216 if horizon_h > 180 else 168}"
    return int(cfg[key])


def _nearest_available(available: np.ndarray, target: float, used: set[float]) -> float | None:
    candidates = [float(t) for t in available if float(t) not in used]
    if not candidates:
        return None
    return min(candidates, key=lambda t: abs(t - float(target)))


def _select_times(available: np.ndarray, targets: list[float], max_n: int) -> np.ndarray:
    used: set[float] = set()
    selected: list[float] = []
    for target in targets:
        t = _nearest_available(available, target, used)
        if t is not None:
            used.add(t)
            selected.append(t)
        if len(selected) >= max_n:
            break
    if len(selected) < max_n:
        remaining = [float(t) for t in available if float(t) not in used]
        if remaining:
            indices = np.linspace(0, len(remaining) - 1, max_n - len(selected)).round().astype(int)
            for idx in indices:
                t = remaining[int(idx)]
                if t not in used:
                    used.add(t)
                    selected.append(t)
    return np.asarray(sorted(selected), dtype=float)


def full_sample_times(design: base.FutureDesign, policy: str) -> np.ndarray:
    available = base.operational_sample_times(design.horizon_h, policy="balanced")
    max_n = _max_for_horizon(policy, float(design.horizon_h), "full")
    pulse_targets: list[float] = []
    for schedule in design.pulses.values():
        for t, _amount in schedule:
            pulse_targets.extend([float(t), float(t) + 6.0])
    phase_targets = [2.0, 8.0, 26.0, 32.0, 50.0, 56.0, 74.0, 80.0, 98.0, 104.0, 122.0, 128.0, 170.0, 176.0]
    targets = [t for t in phase_targets + pulse_targets if 0.0 < t <= float(design.horizon_h)]
    return _select_times(available, targets, max_n)


def small_sample_times(design: base.FutureDesign, policy: str) -> np.ndarray:
    available = base.operational_sample_times(design.horizon_h, policy="balanced")
    full = set(float(t) for t in full_sample_times(design, policy))
    available = np.asarray([float(t) for t in available if float(t) not in full], dtype=float)
    max_n = _max_for_horizon(policy, float(design.horizon_h), "small")
    pulse_targets: list[float] = []
    for schedule in design.pulses.values():
        for t, _amount in schedule:
            pulse_targets.extend([max(2.0, float(t) - 2.0), float(t) + 2.0])
    phase_targets = [4.0, 6.0, 28.0, 30.0, 52.0, 54.0, 76.0, 78.0, 100.0, 102.0, 124.0, 126.0, 172.0, 174.0]
    targets = [t for t in pulse_targets + phase_targets if 0.0 < t <= float(design.horizon_h)]
    return _select_times(available, targets, max_n)


def all_liquid_times(design: base.FutureDesign, policy: str) -> np.ndarray:
    return np.asarray(sorted(set(full_sample_times(design, policy)).union(set(small_sample_times(design, policy)))), dtype=float)


def o2_times(design: base.FutureDesign, policy: str) -> np.ndarray:
    liquid = all_liquid_times(design, policy)
    early = liquid[liquid <= min(float(design.horizon_h), 72.0)]
    return np.unique(np.concatenate([np.array([0.0], dtype=float), early]))


def dense_time_grid(design: base.FutureDesign, policy: str) -> np.ndarray:
    times = set(np.arange(0.0, float(design.horizon_h) + 1e-9, 2.0).round(8))
    times.update(float(t) for t in all_liquid_times(design, policy))
    times.update(float(t) for t in o2_times(design, policy))
    times.update(float(t) for t in final.co2_times(design.horizon_h))
    times.add(float(design.horizon_h))
    for schedule in design.pulses.values():
        for time_h, _amount in schedule:
            times.add(round(max(0.0, float(time_h) - 2.0), 8))
            times.add(round(float(time_h), 8))
            times.add(round(min(float(design.horizon_h), float(time_h) + 2.0), 8))
    return np.asarray(sorted(t for t in times if 0.0 <= t <= float(design.horizon_h)), dtype=float)


def simulate(design: base.FutureDesign, theta: dict[str, float], policy: str) -> final.SimV2 | None:
    original = final.dense_time_grid
    try:
        final.dense_time_grid = lambda d, p: dense_time_grid(d, policy)
        return final.simulate_v2(design, theta, policy)
    finally:
        final.dense_time_grid = original


def future_residual_vector(theta: dict[str, float], design: base.FutureDesign, policy: str) -> np.ndarray:
    sim = simulate(design, theta, policy)
    if sim is None:
        return np.ones(1000, dtype=float) * 1e6
    residuals: list[np.ndarray] = []
    full = [t for t in full_sample_times(design, policy) if t in sim.core.index]
    small = [t for t in small_sample_times(design, policy) if t in sim.core.index]
    liquid = sorted(set(full).union(small))
    for state in base.STATE_NAMES:
        times = full if state == "E" else liquid
        if not times:
            continue
        center = sim.core.loc[times, state].to_numpy(dtype=float)
        residuals.append(center / base.sigma_for_state(state, center))
    for state in ("Pyr", "AcAld", "Acetate"):
        if liquid:
            center = sim.secondary.loc[liquid, state].to_numpy(dtype=float)
            residuals.append(center / joint.SIGMA[state])
    oxy = [t for t in o2_times(design, policy) if t in sim.secondary.index]
    if oxy:
        residuals.append(sim.secondary.loc[oxy, "O2"].to_numpy(dtype=float) / joint.SIGMA["O2"])
    ctimes = [t for t in final.co2_times(design.horizon_h) if t in sim.secondary.index]
    if ctimes:
        residuals.append(sim.secondary.loc[ctimes, "CO2"].to_numpy(dtype=float) / joint.SIGMA["CO2"])
    if full:
        liq_aroma = sim.aromas_liq.loc[full]
        for species in joint.AROMA_SPECIES:
            residuals.append(liq_aroma[species].to_numpy(dtype=float) / joint.SIGMA[f"{species}_liq"])
            residuals.append(np.array([sim.aromas_cond[species] / joint.SIGMA[f"{species}_cond"]], dtype=float))
    return np.concatenate(residuals)


def pulse_stock_volume_ml(channel: str, amount: float, current_volume_ml: float) -> float:
    stock = STOCK_CONCENTRATIONS[channel]["concentration"]
    volume_l = current_volume_ml / 1000.0
    if channel == "N":
        mass_mg = float(amount) * 1000.0 * volume_l  # kg/m3 = g/L, convert to mg/L.
        return mass_mg / stock
    mass_g = float(amount) * volume_l
    return mass_g / stock * 1000.0


def selected_protocol(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], policy: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    protocol_rows = []
    schedule_rows = []
    volume_rows = []
    hybrid = selected[selected["objective"].eq("hybrid")].sort_values("campaign_order").reset_index(drop=True)
    for _idx, record in hybrid.iterrows():
        design = designs[str(record["candidate"])]
        order = int(record["campaign_order"])
        lot = int((order - 1) // 3 + 1)
        start_dt = datetime.fromisoformat(final.START_DATES[lot - 1])
        volume_ml = INITIAL_REACTOR_VOLUME_ML
        protocol_rows.append(
            {
                "campaign_order": order,
                "lot": lot,
                "start_datetime": start_dt.isoformat(sep=" "),
                "candidate": design.name,
                "medium": design.medium,
                "family": design.family,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{v:g}" for v in design.temperature_segments),
                "initial_N_mg_l": float(design.initials.get("N", np.nan)) * 1000.0,
                "initial_G_g_l": float(design.initials.get("G", np.nan)),
                "initial_F_g_l": float(design.initials.get("F", np.nan)),
                "initial_E_g_l": float(design.initials.get("E", np.nan)),
                "initial_X_kg_m3": float(design.initials.get("X", np.nan)),
                "full_samples": int(len(full_sample_times(design, policy))),
                "small_samples": int(len(small_sample_times(design, policy))),
                "planned_sample_volume_ml": float(len(full_sample_times(design, policy)) * SAMPLE_VOLUME_ML["full"] + len(small_sample_times(design, policy)) * SAMPLE_VOLUME_ML["small"]),
                "rationale": design.rationale,
            }
        )
        event_rows = []
        for kind, times in [("full_liquid_sample", full_sample_times(design, policy)), ("small_liquid_sample", small_sample_times(design, policy))]:
            for t in times:
                dt = start_dt + timedelta(hours=float(t))
                event_rows.append(
                    {
                        "campaign_order": order,
                        "lot": lot,
                        "candidate": design.name,
                        "event": kind,
                        "relative_time_h": float(t),
                        "datetime": dt.isoformat(sep=" "),
                        "weekday": dt.strftime("%A"),
                        "details": "ethanol+aromas+full panel" if kind == "full_liquid_sample" else "small panel: biomass/YAN/sugars/glycerol/secondary acids; no ethanol/aromas",
                        "channel": "",
                        "amount": np.nan,
                        "volume_ml": SAMPLE_VOLUME_ML["full"] if kind == "full_liquid_sample" else SAMPLE_VOLUME_ML["small"],
                        "operational": bool(dt.weekday() < 5 and 9 <= dt.hour <= 16),
                    }
                )
        for t in o2_times(design, policy):
            dt = start_dt + timedelta(hours=float(t))
            event_rows.append(
                {
                    "campaign_order": order,
                    "lot": lot,
                    "candidate": design.name,
                    "event": "DO_measurement",
                    "relative_time_h": float(t),
                    "datetime": dt.isoformat(sep=" "),
                    "weekday": dt.strftime("%A"),
                    "details": "DO spot measurement; no controlled oxygen pulse",
                    "channel": "",
                    "amount": np.nan,
                    "volume_ml": 0.0,
                    "operational": bool(t == 0.0 or (dt.weekday() < 5 and 9 <= dt.hour <= 16)),
                }
            )
        for channel, rows in design.pulses.items():
            for t, amount in rows:
                dt = start_dt + timedelta(hours=float(t))
                event_rows.append(
                    {
                        "campaign_order": order,
                        "lot": lot,
                        "candidate": design.name,
                        "event": "manual_pulse",
                        "relative_time_h": float(t),
                        "datetime": dt.isoformat(sep=" "),
                        "weekday": dt.strftime("%A"),
                        "details": f"manual {channel} addition from stock",
                        "channel": channel,
                        "amount": float(amount),
                        "volume_ml": np.nan,
                        "operational": bool(dt.weekday() < 5 and 9 <= dt.hour <= 16),
                    }
                )
        end_dt = start_dt + timedelta(hours=float(design.horizon_h))
        retrieval = end_dt
        if not (retrieval.weekday() < 5 and 9 <= retrieval.hour <= 16):
            cursor = retrieval
            while not (cursor.weekday() < 5 and 10 <= cursor.hour <= 16):
                cursor += timedelta(hours=1)
            retrieval = cursor
        event_rows.append(
            {
                "campaign_order": order,
                "lot": lot,
                "candidate": design.name,
                "event": "terminal_condensate_retrieval",
                "relative_time_h": float((retrieval - start_dt).total_seconds() / 3600.0),
                "datetime": retrieval.isoformat(sep=" "),
                "weekday": retrieval.strftime("%A"),
                "details": "single terminal condensate sample",
                "channel": "",
                "amount": np.nan,
                "volume_ml": 0.0,
                "operational": True,
            }
        )
        event_rows = sorted(event_rows, key=lambda row: (row["relative_time_h"], row["event"]))
        for row in event_rows:
            before = volume_ml
            if row["event"] in {"full_liquid_sample", "small_liquid_sample"}:
                volume_ml -= float(row["volume_ml"])
            elif row["event"] == "manual_pulse":
                pulse_volume = pulse_stock_volume_ml(str(row["channel"]), float(row["amount"]), volume_ml)
                row["volume_ml"] = pulse_volume
                volume_ml += pulse_volume
            row["volume_before_ml"] = before
            row["volume_after_ml"] = volume_ml
            schedule_rows.append(row)
        volume_rows.append(
            {
                "campaign_order": order,
                "candidate": design.name,
                "initial_volume_ml": INITIAL_REACTOR_VOLUME_ML,
                "final_volume_ml": volume_ml,
                "net_volume_removed_ml": INITIAL_REACTOR_VOLUME_ML - volume_ml,
                "below_min_final_volume": bool(volume_ml < MIN_FINAL_VOLUME_ML),
                "max_single_pulse_volume_ml": max([r["volume_ml"] for r in event_rows if r["event"] == "manual_pulse"] or [0.0]),
            }
        )
    return pd.DataFrame(protocol_rows), pd.DataFrame(schedule_rows), pd.DataFrame(volume_rows)


def evaluate_policy(theta: dict[str, float], prior: np.ndarray, designs: dict[str, base.FutureDesign], policy: str, step: float, campaign_size: int):
    rows = []
    fims: dict[str, np.ndarray] = {}
    for idx, design in enumerate(designs.values(), start=1):
        print(f"[{policy}] candidate {idx}/{len(designs)} {design.name}", flush=True)
        try:
            fim = final.finite_difference_fim(theta, TARGET_PARAMETERS, lambda th, d=design: future_residual_vector(th, d, policy), step)
            fims[design.name] = fim
            combined = prior + fim
            full_n = len(full_sample_times(design, policy))
            small_n = len(small_sample_times(design, policy))
            sample_volume = full_n * SAMPLE_VOLUME_ML["full"] + small_n * SAMPLE_VOLUME_ML["small"]
            rows.append(
                {
                    "policy": policy,
                    "candidate": design.name,
                    "family": design.family,
                    "medium": design.medium,
                    "horizon_h": design.horizon_h,
                    "n_full_samples": int(full_n),
                    "n_small_samples": int(small_n),
                    "sample_volume_ml": float(sample_volume),
                    "hybrid_score": joint.score_fim(combined, TARGET_PARAMETERS, "hybrid"),
                    "dopt_score": joint.score_fim(combined, TARGET_PARAMETERS, "d_opt"),
                    **joint.fim_metrics(fim, TARGET_PARAMETERS, prefix="new_"),
                    **joint.fim_metrics(combined, TARGET_PARAMETERS, prefix="combined_"),
                    **joint.variance_reduction(prior, combined, TARGET_PARAMETERS),
                }
            )
        except Exception as err:
            rows.append({"policy": policy, "candidate": design.name, "family": design.family, "medium": design.medium, "status": "failed", "error": f"{type(err).__name__}: {err}"})
    ranking = pd.DataFrame(rows)
    selected_hybrid, fim_hybrid = joint.greedy_select(fims, designs, prior, TARGET_PARAMETERS, campaign_size, "hybrid")
    selected_hybrid.insert(0, "policy", policy)
    selected_dopt, _fim_dopt = joint.greedy_select(fims, designs, prior, TARGET_PARAMETERS, campaign_size, "d_opt")
    selected_dopt.insert(0, "policy", policy)
    selected = pd.concat([selected_hybrid, selected_dopt], ignore_index=True)
    return ranking, selected, selected_hybrid, fim_hybrid


def policy_summary(selected_by_policy: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for policy, selected in selected_by_policy.items():
        last = selected[selected["objective"].eq("hybrid")].sort_values("campaign_order").tail(1)
        if last.empty:
            continue
        row = last.iloc[0].to_dict()
        rows.append(
            {
                "policy": policy,
                "campaign_logdet": row.get("campaign_logdet", np.nan),
                "campaign_min_relative_eigenvalue": row.get("campaign_min_relative_eigenvalue", np.nan),
                "campaign_trace_inv": row.get("campaign_trace_inv", np.nan),
                "new_param_mean_var_reduction": row.get("new_param_mean_var_reduction", np.nan),
                "new_param_worst_var_reduction": row.get("new_param_worst_var_reduction", np.nan),
                "selected_candidates": ", ".join(selected[selected["objective"].eq("hybrid")]["candidate"].astype(str).tolist()),
            }
        )
    return pd.DataFrame(rows).sort_values(["new_param_worst_var_reduction", "campaign_logdet"], ascending=False)


def choose_policy(summary: pd.DataFrame) -> str:
    if summary.empty:
        return "full12_small8"
    return str(summary.iloc[0]["policy"])


def plot_volume(protocol: pd.DataFrame, schedule: pd.DataFrame) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for candidate, group in schedule.groupby("candidate"):
        fig, ax = plt.subplots(figsize=(9, 4))
        group = group.sort_values("relative_time_h")
        ax.step(group["relative_time_h"], group["volume_after_ml"], where="post")
        ax.axhline(MIN_FINAL_VOLUME_ML, color="tab:red", linestyle="--", label="minimum target")
        full = group[group["event"].eq("full_liquid_sample")]
        small = group[group["event"].eq("small_liquid_sample")]
        pulse = group[group["event"].eq("manual_pulse")]
        ax.scatter(full["relative_time_h"], full["volume_after_ml"], label="full sample", s=24)
        ax.scatter(small["relative_time_h"], small["volume_after_ml"], label="small sample", s=18)
        ax.scatter(pulse["relative_time_h"], pulse["volume_after_ml"], label="pulse", s=36)
        ax.set_title(candidate)
        ax.set_xlabel("time [h]")
        ax.set_ylabel("reactor volume [mL]")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(plot_dir / f"volume_profile_{candidate}.png", dpi=170)
        plt.close(fig)


def write_report(summary: pd.DataFrame, selected: pd.DataFrame, protocol: pd.DataFrame, volume: pd.DataFrame, stock: pd.DataFrame) -> None:
    report = f"""# Volume-constrained operational DOE

## Sampling volume constraint

The previous campaign used every operational liquid sample as a full analytical sample. With 40 mL for ethanol and 10 mL for aromas, that would remove about 1.4 L from a 2 L reactor for most 216 h fermentations. This version separates:

- `full_liquid_sample`: 50 mL, ethanol and aroma must be sampled together.
- `small_liquid_sample`: 8 mL, no ethanol/aroma; intended for biomass, sugars, YAN/PAN/NH4, glycerol, pyruvate, acetaldehyde, and acetate.

The selected policy keeps the final estimated reactor volume above {MIN_FINAL_VOLUME_ML:.0f} mL after sampling and stock additions.

## Policy benchmark

{summary.to_markdown(index=False)}

## Selected campaign

{selected[selected["objective"].eq("hybrid")].to_markdown(index=False)}

## Protocol

{protocol.to_markdown(index=False)}

## Volume audit

{volume.to_markdown(index=False)}

## Stock concentration assumptions

{stock.to_markdown(index=False)}

## Interpretation

- The design remains model-based: candidate selection is still driven by the extended FIM for core fermentation, reduced secondary metabolites, CO2, and aroma outputs.
- Full ethanol/aroma samples are now the scarce resource. They are concentrated around early phase, pulse response, mid-fermentation, and late/final windows.
- Small samples preserve core kinetic information without forcing a 50 mL withdrawal at every operational time.
- No automated injection is required by the selected schedule; all selected pulses are inside working windows.
"""
    (RESULTS_DIR / "volume_constrained_doe_report.md").write_text(report, encoding="utf-8")


def stock_table() -> pd.DataFrame:
    rows = []
    for channel, cfg in STOCK_CONCENTRATIONS.items():
        rows.append(
            {
                "channel": channel,
                "stock_concentration": cfg["concentration"],
                "unit": cfg["unit"],
                "max_single_pulse_ml_target": cfg["max_single_pulse_ml"],
                "comment": {
                    "N": "20 mg YAN-N/mL gives 5-6 mL for a 50-60 mg/L pulse in 2 L.",
                    "G": "700 g/L is concentrated and viscous but keeps 35 g/L in 2 L near 100 mL.",
                    "F": "Use same concentration as glucose because fructose is treated as equivalent sugar.",
                    "E": "Pure ethanol basis; not used in selected hybrid campaign.",
                    "X": "Dry-biomass-equivalent concentrated slurry; verify viable cell concentration before dosing.",
                }[channel],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Volume-constrained final operational DOE.")
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--step", type=float, default=0.04)
    parser.add_argument("--policies", nargs="*", default=list(VOLUME_POLICIES))
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    theta = final.load_theta_final()
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    prior, prior_estimability, prior_eigen = final.current_prior_fim(theta, batches, args.step)
    prior_estimability.to_csv(RESULTS_DIR / "current_prior_estimability.csv", index=False)
    prior_eigen.to_csv(RESULTS_DIR / "current_prior_eigen_diagnostics.csv", index=False)
    designs = final.candidate_library(data)

    rankings = []
    selected_by_policy = {}
    fims_by_policy = {}
    for policy in args.policies:
        ranking, selected, selected_hybrid, fim_hybrid = evaluate_policy(theta, prior, designs, policy, args.step, args.campaign_size)
        ranking.to_csv(RESULTS_DIR / f"candidate_ranking_{policy}.csv", index=False)
        selected.to_csv(RESULTS_DIR / f"selected_campaign_{policy}.csv", index=False)
        rankings.append(ranking)
        selected_by_policy[policy] = selected
        fims_by_policy[policy] = fim_hybrid

    all_rankings = pd.concat(rankings, ignore_index=True, sort=False)
    all_rankings.to_csv(RESULTS_DIR / "candidate_rankings_all_policies.csv", index=False)
    summary = policy_summary(selected_by_policy)
    summary.to_csv(RESULTS_DIR / "policy_summary.csv", index=False)
    chosen = choose_policy(summary)
    selected = selected_by_policy[chosen]
    selected.to_csv(RESULTS_DIR / "selected_campaign.csv", index=False)
    combined = prior + fims_by_policy[chosen]
    final.parameter_estimability(combined, theta, TARGET_PARAMETERS, f"post_campaign_{chosen}").to_csv(RESULTS_DIR / "post_campaign_estimability.csv", index=False)

    protocol, schedule, volume = selected_protocol(selected, designs, chosen)
    protocol.to_csv(RESULTS_DIR / "campaign_protocol.csv", index=False)
    schedule.to_csv(RESULTS_DIR / "operational_schedule.csv", index=False)
    volume.to_csv(RESULTS_DIR / "volume_audit.csv", index=False)
    stock = stock_table()
    stock.to_csv(RESULTS_DIR / "stock_concentration_assumptions.csv", index=False)
    plot_volume(protocol, schedule)
    write_report(summary, selected, protocol, volume, stock)
    print(f"[done] chosen policy {chosen}", flush=True)


if __name__ == "__main__":
    main()
