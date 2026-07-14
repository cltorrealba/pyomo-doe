from __future__ import annotations

import math
import os
import sys
from dataclasses import replace
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
import run_final_operational_doe_volume_constrained as vc
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "lot1_pulse_timing_mbdoe"
START_REAL = datetime(2026, 6, 15, 15, 30)
HOURS = (10, 12, 14, 16)
TARGET_PARAMETERS = final.TARGET_PARAMETERS

SCENARIOS = {
    "early_10_old_clock": 10,
    "midday_12_response": 12,
    "late_16_near_nominal": 16,
}

FOCUS_CORE = ("qN", "betaG0", "betaF0", "qEG", "qEF", "gammaG0", "gammaF0")
FOCUS_SECONDARY = final.SECONDARY_TARGETS
FOCUS_AROMA = final.AROMA_TARGETS


def actual_operational_times(horizon_h: float) -> np.ndarray:
    points = [0.0]
    current = START_REAL.date()
    end = START_REAL + timedelta(hours=float(horizon_h))
    while datetime.combine(current, datetime.min.time()) <= end + timedelta(days=1):
        if current > START_REAL.date() and datetime.combine(current, datetime.min.time()).weekday() < 5:
            for hour in HOURS:
                dt = datetime.combine(current, datetime.min.time()) + timedelta(hours=hour)
                t = (dt - START_REAL).total_seconds() / 3600.0
                if 0.0 < t <= float(horizon_h):
                    points.append(float(t))
        current = current + timedelta(days=1)
    return np.asarray(sorted(set(round(t, 8) for t in points)), dtype=float)


def move_pulse_time(nominal_h: float, strategy_hour: int) -> float:
    nominal_dt = START_REAL + timedelta(hours=float(nominal_h))
    moved_dt = datetime.combine(nominal_dt.date(), datetime.min.time()) + timedelta(hours=int(strategy_hour))
    return float((moved_dt - START_REAL).total_seconds() / 3600.0)


def scenario_design(design: base.FutureDesign, scenario: str) -> base.FutureDesign:
    strategy_hour = SCENARIOS[scenario]
    initials = dict(design.initials)
    if design.name == "synthetic_lit_SM410_18C_highN_ester":
        initials["N"] = 0.320  # user-reported first-lot deviation: 320 mg/L instead of 410 mg/L.
    pulses: dict[str, tuple[tuple[float, float], ...]] = {}
    for channel, schedule in design.pulses.items():
        moved = []
        for nominal_t, amount in schedule:
            moved_t = move_pulse_time(float(nominal_t), strategy_hour)
            if 0.0 < moved_t <= float(design.horizon_h):
                moved.append((moved_t, float(amount)))
        pulses[channel] = tuple(sorted(moved))
    return base.FutureDesign(
        design.name,
        design.family,
        design.medium,
        design.horizon_h,
        initials,
        design.temperature_segments,
        pulses,
        design.rationale,
    )


def _nearest_available(available: np.ndarray, target: float, used: set[float]) -> float | None:
    candidates = [float(t) for t in available if float(t) not in used]
    if not candidates:
        return None
    return min(candidates, key=lambda value: abs(value - float(target)))


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


def full_sample_times(design: base.FutureDesign) -> np.ndarray:
    available = actual_operational_times(design.horizon_h)
    max_n = vc._max_for_horizon("full14_small6", float(design.horizon_h), "full")
    pulse_targets: list[float] = []
    for schedule in design.pulses.values():
        for t, _amount in schedule:
            pulse_targets.extend([float(t), float(t) + 4.0, float(t) + 6.0])
    phase_targets = [0.0, 18.5, 24.5, 42.5, 48.5, 66.5, 72.5, 90.5, 96.5, 114.5, 120.5, 162.5, 168.5]
    targets = [t for t in phase_targets + pulse_targets if 0.0 <= t <= float(design.horizon_h)]
    return _select_times(available, targets, max_n)


def small_sample_times(design: base.FutureDesign) -> np.ndarray:
    available = actual_operational_times(design.horizon_h)
    full = set(float(t) for t in full_sample_times(design))
    available = np.asarray([float(t) for t in available if float(t) not in full], dtype=float)
    max_n = vc._max_for_horizon("full14_small6", float(design.horizon_h), "small")
    pulse_targets: list[float] = []
    for schedule in design.pulses.values():
        for t, _amount in schedule:
            pulse_targets.extend([max(0.0, float(t) - 2.0), float(t) + 2.0])
    phase_targets = [20.5, 22.5, 44.5, 46.5, 68.5, 70.5, 92.5, 94.5, 116.5, 118.5, 164.5, 166.5]
    targets = [t for t in pulse_targets + phase_targets if 0.0 <= t <= float(design.horizon_h)]
    return _select_times(available, targets, max_n)


def all_liquid_times(design: base.FutureDesign) -> np.ndarray:
    return np.asarray(sorted(set(full_sample_times(design)).union(set(small_sample_times(design)))), dtype=float)


def o2_times(design: base.FutureDesign) -> np.ndarray:
    liquid = all_liquid_times(design)
    early = liquid[liquid <= min(float(design.horizon_h), 72.0)]
    return np.unique(np.concatenate([np.array([0.0], dtype=float), early]))


def dense_time_grid(design: base.FutureDesign) -> np.ndarray:
    times = set(np.arange(0.0, float(design.horizon_h) + 1e-9, 2.0).round(8))
    times.update(float(t) for t in all_liquid_times(design))
    times.update(float(t) for t in o2_times(design))
    times.update(float(t) for t in final.co2_times(design.horizon_h))
    times.add(float(design.horizon_h))
    for schedule in design.pulses.values():
        for time_h, _amount in schedule:
            times.add(round(max(0.0, float(time_h) - 2.0), 8))
            times.add(round(float(time_h), 8))
            times.add(round(min(float(design.horizon_h), float(time_h) + 2.0), 8))
            times.add(round(min(float(design.horizon_h), float(time_h) + 4.0), 8))
    return np.asarray(sorted(t for t in times if 0.0 <= t <= float(design.horizon_h)), dtype=float)


def simulate_design(design: base.FutureDesign, theta: dict[str, float]) -> final.SimV2 | None:
    times = dense_time_grid(design)
    core = base.simulate(design, theta, times)
    if core is None:
        return None
    core = core.sort_index()
    sec = v2.integrate_secondary_v2(design, theta, core)
    try:
        liq, loss, cond = joint.integrate_aroma_euler(design, theta, core)
    except Exception:
        return None
    return final.SimV2(core=core, secondary=sec, aromas_liq=liq, aromas_loss=loss, aromas_cond=cond)


def residual_vector(theta: dict[str, float], design: base.FutureDesign) -> np.ndarray:
    sim = simulate_design(design, theta)
    if sim is None:
        return np.ones(1000, dtype=float) * 1e6
    residuals: list[np.ndarray] = []
    full = [t for t in full_sample_times(design) if t in sim.core.index]
    small = [t for t in small_sample_times(design) if t in sim.core.index]
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
    oxy = [t for t in o2_times(design) if t in sim.secondary.index]
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


def finite_difference_fim(theta: dict[str, float], design: base.FutureDesign, step: float) -> np.ndarray:
    return final.finite_difference_fim(theta, TARGET_PARAMETERS, lambda th: residual_vector(th, design), step)


def focus_mean(reductions: dict[str, float], names: tuple[str, ...]) -> float:
    values = [float(reductions.get(f"var_reduction_{name}", np.nan)) for name in names]
    values = [v for v in values if np.isfinite(v)]
    return float(np.mean(values)) if values else np.nan


def post_pulse_counts(design: base.FutureDesign) -> dict[str, int]:
    full = full_sample_times(design)
    small = small_sample_times(design)
    rows: dict[str, int] = {}
    idx = 0
    for channel, schedule in design.pulses.items():
        for t, _amount in schedule:
            idx += 1
            rows[f"pulse_{idx}_{channel}_time_h"] = float(t)
            rows[f"pulse_{idx}_{channel}_post_full_0_6h"] = int(np.sum((full > t) & (full <= t + 6.0)))
            rows[f"pulse_{idx}_{channel}_post_small_0_6h"] = int(np.sum((small > t) & (small <= t + 6.0)))
            rows[f"pulse_{idx}_{channel}_post_liquid_0_6h"] = int(np.sum((all_liquid_times(design) > t) & (all_liquid_times(design) <= t + 6.0)))
    return rows


def plot_schedules(scenario_rows: list[dict], designs_by_scenario: dict[str, list[base.FutureDesign]]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for scenario, designs in designs_by_scenario.items():
        fig, axes = plt.subplots(len(designs), 1, figsize=(11, 7), sharex=True)
        if len(designs) == 1:
            axes = [axes]
        for ax, design in zip(axes, designs):
            full = full_sample_times(design)
            small = small_sample_times(design)
            ax.vlines(full, 0.0, 1.0, color="tab:blue", label="full")
            ax.vlines(small, 0.0, 0.65, color="tab:orange", label="small")
            ax.vlines(o2_times(design), 0.0, 0.35, color="tab:green", label="DO")
            for channel, schedule in design.pulses.items():
                for t, amount in schedule:
                    ax.axvline(float(t), color="tab:red", lw=2.0)
                    ax.text(float(t), 1.08, f"{channel} {amount:g}", rotation=90, va="bottom", ha="center", fontsize=8)
            ax.set_title(design.name)
            ax.set_ylim(0, 1.3)
            ax.set_yticks([])
            ax.grid(True, axis="x", alpha=0.25)
        axes[0].legend(fontsize=8, ncol=4, loc="upper right")
        axes[-1].set_xlabel("process time from real inoculation [h]")
        fig.suptitle(scenario)
        fig.tight_layout()
        fig.savefig(plot_dir / f"lot1_schedule_{scenario}.png", dpi=170)
        plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    theta = final.load_theta_final()
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    prior, _, _ = final.current_prior_fim(theta, batches, step=0.04)
    designs = final.candidate_library(data)
    protocol = pd.read_csv(SCRIPT_DIR / "results" / "final_operational_doe_volume_constrained" / "campaign_protocol.csv")
    lot1_names = protocol.loc[protocol["lot"].eq(1), "candidate"].astype(str).tolist()

    rows = []
    schedule_rows = []
    designs_by_scenario: dict[str, list[base.FutureDesign]] = {}
    for scenario in SCENARIOS:
        lot_fim = np.zeros((len(TARGET_PARAMETERS), len(TARGET_PARAMETERS)), dtype=float)
        scenario_designs = []
        for name in lot1_names:
            design = scenario_design(designs[name], scenario)
            scenario_designs.append(design)
            print(f"[{scenario}] {name}", flush=True)
            fim = finite_difference_fim(theta, design, step=0.04)
            lot_fim += fim
            for t in full_sample_times(design):
                schedule_rows.append({"scenario": scenario, "candidate": name, "event": "full", "time_h": float(t)})
            for t in small_sample_times(design):
                schedule_rows.append({"scenario": scenario, "candidate": name, "event": "small", "time_h": float(t)})
            for t in o2_times(design):
                schedule_rows.append({"scenario": scenario, "candidate": name, "event": "DO", "time_h": float(t)})
            for channel, schedule in design.pulses.items():
                for t, amount in schedule:
                    schedule_rows.append({"scenario": scenario, "candidate": name, "event": f"pulse_{channel}", "time_h": float(t), "amount": float(amount)})
        designs_by_scenario[scenario] = scenario_designs
        combined = prior + lot_fim
        reductions = joint.variance_reduction(prior, combined, TARGET_PARAMETERS)
        row = {
            "scenario": scenario,
            **joint.fim_metrics(lot_fim, TARGET_PARAMETERS, prefix="lot1_new_"),
            **joint.fim_metrics(combined, TARGET_PARAMETERS, prefix="combined_"),
            **reductions,
            "focus_core_mean_var_reduction": focus_mean(reductions, FOCUS_CORE),
            "focus_secondary_mean_var_reduction": focus_mean(reductions, FOCUS_SECONDARY),
            "focus_aroma_mean_var_reduction": focus_mean(reductions, FOCUS_AROMA),
        }
        # Add post-pulse observation counts aggregated across lot.
        pulse_liquid_counts = []
        pulse_full_counts = []
        for design in scenario_designs:
            counts = post_pulse_counts(design)
            pulse_liquid_counts.extend(v for k, v in counts.items() if k.endswith("post_liquid_0_6h"))
            pulse_full_counts.extend(v for k, v in counts.items() if k.endswith("post_full_0_6h"))
        row["mean_post_pulse_liquid_0_6h"] = float(np.mean(pulse_liquid_counts)) if pulse_liquid_counts else np.nan
        row["min_post_pulse_liquid_0_6h"] = float(np.min(pulse_liquid_counts)) if pulse_liquid_counts else np.nan
        row["mean_post_pulse_full_0_6h"] = float(np.mean(pulse_full_counts)) if pulse_full_counts else np.nan
        row["min_post_pulse_full_0_6h"] = float(np.min(pulse_full_counts)) if pulse_full_counts else np.nan
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values(["combined_logdet", "focus_core_mean_var_reduction"], ascending=False)
    summary.to_csv(RESULTS_DIR / "lot1_pulse_timing_mbdoe_summary.csv", index=False)
    pd.DataFrame(schedule_rows).to_csv(RESULTS_DIR / "lot1_pulse_timing_schedules.csv", index=False)
    plot_schedules(rows, designs_by_scenario)

    report = [
        "# Lot 1 pulse-timing MBDoE benchmark",
        "",
        f"Real inoculation start assumed: `{START_REAL.isoformat(sep=' ')}`.",
        "",
        "Known actual-initial override used in this quick benchmark: `synthetic_lit_SM410_18C_highN_ester` N = 0.320 kg/m3 (320 mg/L). Other initials remain nominal until measured values are loaded.",
        "",
        "Scenarios:",
        "",
        "- `early_10_old_clock`: pulses at 10:00 on the target calendar day, equivalent to keeping the old 08:00-clock schedule after a 15:30 inoculation.",
        "- `midday_12_response`: pulses at 12:00 on the target calendar day to allow +2 h and +4 h response samples.",
        "- `late_16_near_nominal`: pulses at 16:00 on the target calendar day, closest to the original process-time target without leaving the work window.",
        "",
        "## Summary",
        "",
        summary[
            [
                "scenario",
                "lot1_new_logdet",
                "combined_logdet",
                "combined_trace_inv",
                "new_param_mean_var_reduction",
                "new_param_worst_var_reduction",
                "focus_core_mean_var_reduction",
                "focus_secondary_mean_var_reduction",
                "focus_aroma_mean_var_reduction",
                "mean_post_pulse_liquid_0_6h",
                "min_post_pulse_liquid_0_6h",
                "mean_post_pulse_full_0_6h",
                "min_post_pulse_full_0_6h",
            ]
        ].to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "Higher `combined_logdet` and lower `combined_trace_inv` are better. Post-pulse counts help detect designs that cannot observe the immediate response to a pulse.",
    ]
    (RESULTS_DIR / "lot1_pulse_timing_mbdoe_report.md").write_text("\n".join(report), encoding="utf-8")
    print(summary[["scenario", "lot1_new_logdet", "combined_logdet", "combined_trace_inv", "focus_core_mean_var_reduction", "mean_post_pulse_liquid_0_6h"]].to_string(index=False))


if __name__ == "__main__":
    main()
