from __future__ import annotations

import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_final_operational_doe_v2 as final
import run_lot1_pulse_timing_mbdoe as timing
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "lot1_express_optimal_sampling"
TARGET_PARAMETERS = final.TARGET_PARAMETERS
SCENARIOS = ("midday_12_response", "late_16_near_nominal")
FULL_BUDGET_216 = 14
SMALL_BUDGET_216 = 6
MAX_LIQUID_PER_DAY_PER_FERMENTATION = 4


@dataclass(frozen=True)
class Group:
    key: str
    scenario: str
    candidate: str
    kind: str
    time_h: float
    date: str
    fim: np.ndarray


@dataclass(frozen=True)
class SimBundle:
    core: pd.DataFrame
    secondary: pd.DataFrame
    aroma_liq: pd.DataFrame
    aroma_cond: pd.Series


def dt_for_time(time_h: float):
    return timing.START_REAL + pd.Timedelta(hours=float(time_h))


def date_for_time(time_h: float) -> str:
    return dt_for_time(time_h).date().isoformat()


def feasible_times(design: base.FutureDesign) -> np.ndarray:
    return timing.actual_operational_times(design.horizon_h)


def simulate_at_times(design: base.FutureDesign, theta: dict[str, float], times: np.ndarray) -> SimBundle | None:
    grid = set(float(t) for t in times)
    grid.update(float(t) for t in final.co2_times(design.horizon_h))
    grid.add(0.0)
    grid.add(float(design.horizon_h))
    for schedule in design.pulses.values():
        for t, _amount in schedule:
            for dt in (-2.0, 0.0, 2.0, 4.0, 6.0):
                grid.add(float(np.clip(float(t) + dt, 0.0, float(design.horizon_h))))
    grid_arr = np.asarray(sorted(grid), dtype=float)
    core = base.simulate(design, theta, grid_arr)
    if core is None:
        return None
    core = core.sort_index()
    sec = v2.integrate_secondary_v2(design, theta, core)
    try:
        liq, _loss, cond = joint.integrate_aroma_euler(design, theta, core)
    except Exception:
        return None
    return SimBundle(core=core, secondary=sec, aroma_liq=liq, aroma_cond=cond)


def obs_vector(sim: SimBundle, kind: str, time_h: float | None = None, design: base.FutureDesign | None = None) -> np.ndarray:
    pieces: list[np.ndarray] = []
    if kind in {"full", "small"}:
        assert time_h is not None
        t = float(time_h)
        core_states = list(base.STATE_NAMES)
        if kind == "small":
            core_states = [state for state in core_states if state != "E"]
        for state in core_states:
            value = np.asarray([float(sim.core.loc[t, state])], dtype=float)
            pieces.append(value / base.sigma_for_state(state, value))
        for state in ("Pyr", "AcAld", "Acetate"):
            pieces.append(np.asarray([float(sim.secondary.loc[t, state])], dtype=float) / joint.SIGMA[state])
        if t <= 72.0:
            pieces.append(np.asarray([float(sim.secondary.loc[t, "O2"])], dtype=float) / joint.SIGMA["O2"])
        if kind == "full":
            for species in joint.AROMA_SPECIES:
                pieces.append(np.asarray([float(sim.aroma_liq.loc[t, species])], dtype=float) / joint.SIGMA[f"{species}_liq"])
    elif kind == "co2_online":
        assert design is not None
        ctimes = [float(t) for t in final.co2_times(design.horizon_h) if float(t) in sim.secondary.index]
        if ctimes:
            pieces.append(sim.secondary.loc[ctimes, "CO2"].to_numpy(dtype=float) / joint.SIGMA["CO2"])
    elif kind == "terminal_condensate":
        for species in joint.AROMA_SPECIES:
            pieces.append(np.asarray([float(sim.aroma_cond[species])], dtype=float) / joint.SIGMA[f"{species}_cond"])
    else:
        raise ValueError(kind)
    if not pieces:
        return np.zeros(0, dtype=float)
    return np.concatenate(pieces)


def perturb(theta: dict[str, float], name: str, sign: float, step: float) -> dict[str, float]:
    out = dict(theta)
    low, high = final.parameter_bounds(name)
    out[name] = float(np.clip(float(theta[name]) * math.exp(sign * step), low, high))
    return out


def group_fim(sim_plus: dict[str, SimBundle], sim_minus: dict[str, SimBundle], group_kind: str, time_h: float | None, design: base.FutureDesign, step: float) -> np.ndarray:
    columns = []
    base_len = None
    for name in TARGET_PARAMETERS:
        vp = obs_vector(sim_plus[name], group_kind, time_h, design)
        vm = obs_vector(sim_minus[name], group_kind, time_h, design)
        if base_len is None:
            base_len = len(vp)
        if len(vp) != base_len or len(vm) != base_len:
            columns.append(np.zeros(base_len or 0, dtype=float))
        else:
            columns.append((vp - vm) / (2.0 * step))
    if not columns or (base_len or 0) == 0:
        return np.zeros((len(TARGET_PARAMETERS), len(TARGET_PARAMETERS)), dtype=float)
    jac = np.column_stack(columns)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def build_groups(theta: dict[str, float], design: base.FutureDesign, scenario: str, step: float) -> tuple[list[Group], list[Group]]:
    times = feasible_times(design)
    sim_plus: dict[str, SimBundle] = {}
    sim_minus: dict[str, SimBundle] = {}
    for name in TARGET_PARAMETERS:
        plus = simulate_at_times(design, perturb(theta, name, +1.0, step), times)
        minus = simulate_at_times(design, perturb(theta, name, -1.0, step), times)
        if plus is None or minus is None:
            raise RuntimeError(f"simulation failed for {design.name}, {name}")
        sim_plus[name] = plus
        sim_minus[name] = minus

    fixed: list[Group] = []
    candidates: list[Group] = []
    # Initial extended sample was already collected at real t=0.
    for kind, time_h in [("full", 0.0)]:
        fim = group_fim(sim_plus, sim_minus, kind, time_h, design, step)
        fixed.append(Group(f"{design.name}:{kind}:{time_h:g}", scenario, design.name, kind, time_h, date_for_time(time_h), fim))
    for kind in ("co2_online", "terminal_condensate"):
        fim = group_fim(sim_plus, sim_minus, kind, None, design, step)
        fixed.append(Group(f"{design.name}:{kind}", scenario, design.name, kind, np.nan, "", fim))

    for time_h in times:
        if abs(float(time_h)) < 1e-9:
            continue
        for kind in ("full", "small"):
            fim = group_fim(sim_plus, sim_minus, kind, float(time_h), design, step)
            candidates.append(Group(f"{design.name}:{kind}:{float(time_h):g}", scenario, design.name, kind, float(time_h), date_for_time(float(time_h)), fim))
    return fixed, candidates


def is_feasible(group: Group, counts: dict, selected_times: dict, day_counts: dict) -> bool:
    if group.kind not in {"full", "small"}:
        return False
    if group.time_h in selected_times[group.candidate]:
        return False
    if counts[(group.candidate, "full")] >= FULL_BUDGET_216 and group.kind == "full":
        return False
    if counts[(group.candidate, "small")] >= SMALL_BUDGET_216 and group.kind == "small":
        return False
    if day_counts[(group.candidate, group.date)] >= MAX_LIQUID_PER_DAY_PER_FERMENTATION:
        return False
    return True


def greedy_select(prior_plus_fixed: np.ndarray, candidates: list[Group], objective: str = "d_opt") -> tuple[pd.DataFrame, np.ndarray]:
    fim = prior_plus_fixed.copy()
    selected: list[dict] = []
    counts = defaultdict(int)
    selected_times = defaultdict(set)
    day_counts = defaultdict(int)
    score = joint.score_fim(fim, TARGET_PARAMETERS, objective)

    while True:
        best = None
        best_score = -np.inf
        for group in candidates:
            if not is_feasible(group, counts, selected_times, day_counts):
                continue
            candidate_score = joint.score_fim(fim + group.fim, TARGET_PARAMETERS, objective)
            if candidate_score > best_score:
                best_score = candidate_score
                best = group
        if best is None or best_score <= score + 1e-10:
            break
        fim = fim + best.fim
        score_gain = best_score - score
        score = best_score
        counts[(best.candidate, best.kind)] += 1
        selected_times[best.candidate].add(best.time_h)
        day_counts[(best.candidate, best.date)] += 1
        selected.append(
            {
                "order": len(selected) + 1,
                "candidate": best.candidate,
                "kind": best.kind,
                "time_h": best.time_h,
                "datetime": dt_for_time(best.time_h).isoformat(sep=" "),
                "date": best.date,
                "score_gain": score_gain,
                "score": score,
            }
        )
        # Stop when all per-fermentation budgets are filled.
        candidates_names = sorted(set(group.candidate for group in candidates))
        done = all(
            counts[(name, "full")] >= FULL_BUDGET_216 and counts[(name, "small")] >= SMALL_BUDGET_216
            for name in candidates_names
        )
        if done:
            break
    return pd.DataFrame(selected), fim


def fixed_budget_summary(selected: pd.DataFrame) -> pd.DataFrame:
    return selected.groupby(["candidate", "kind"]).size().unstack(fill_value=0).reset_index()


def post_pulse_counts(selected: pd.DataFrame, designs: list[base.FutureDesign]) -> dict[str, float]:
    rows = []
    for design in designs:
        chosen = selected[selected["candidate"].eq(design.name)]
        liquid_times = chosen["time_h"].astype(float).to_numpy()
        full_times = chosen[chosen["kind"].eq("full")]["time_h"].astype(float).to_numpy()
        for channel, schedule in design.pulses.items():
            for t, _amount in schedule:
                rows.append(
                    {
                        "candidate": design.name,
                        "channel": channel,
                        "pulse_time_h": float(t),
                        "post_liquid_0_6h": int(np.sum((liquid_times > t) & (liquid_times <= t + 6.0))),
                        "post_full_0_6h": int(np.sum((full_times > t) & (full_times <= t + 6.0))),
                    }
                )
    df = pd.DataFrame(rows)
    if df.empty:
        return {
            "mean_post_pulse_liquid_0_6h": np.nan,
            "min_post_pulse_liquid_0_6h": np.nan,
            "mean_post_pulse_full_0_6h": np.nan,
            "min_post_pulse_full_0_6h": np.nan,
        }
    df.to_csv(RESULTS_DIR / "post_pulse_observation_counts.csv", index=False)
    return {
        "mean_post_pulse_liquid_0_6h": float(df["post_liquid_0_6h"].mean()),
        "min_post_pulse_liquid_0_6h": float(df["post_liquid_0_6h"].min()),
        "mean_post_pulse_full_0_6h": float(df["post_full_0_6h"].mean()),
        "min_post_pulse_full_0_6h": float(df["post_full_0_6h"].min()),
    }


def plot_selection(scenario: str, selected: pd.DataFrame, designs: list[base.FutureDesign]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(designs), 1, figsize=(11, 7), sharex=True)
    if len(designs) == 1:
        axes = [axes]
    for ax, design in zip(axes, designs):
        rows = selected[selected["candidate"].eq(design.name)]
        for kind, color, height in [("full", "tab:blue", 1.0), ("small", "tab:orange", 0.65)]:
            times = rows[rows["kind"].eq(kind)]["time_h"].astype(float).to_numpy()
            ax.vlines(times, 0.0, height, color=color, lw=1.8, label=kind)
        for channel, schedule in design.pulses.items():
            for t, amount in schedule:
                ax.axvline(float(t), color="tab:red", lw=2.0)
                ax.text(float(t), 1.07, f"{channel} {amount:g}", rotation=90, ha="center", va="bottom", fontsize=8)
        ax.set_ylim(0, 1.25)
        ax.set_yticks([])
        ax.set_title(design.name)
        ax.grid(True, axis="x", alpha=0.25)
    axes[0].legend(fontsize=8, ncol=3)
    axes[-1].set_xlabel("process time from real inoculation [h]")
    fig.suptitle(f"Greedy local FIM sampling: {scenario}")
    fig.tight_layout()
    fig.savefig(plot_dir / f"express_sampling_{scenario}.png", dpi=170)
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    theta = final.load_theta_final()
    data = base.load_normalized_data()
    batches = base.make_batches(data)
    prior, _, _ = final.current_prior_fim(theta, batches, step=0.04)
    designs_all = final.candidate_library(data)
    protocol = pd.read_csv(SCRIPT_DIR / "results" / "final_operational_doe_volume_constrained" / "campaign_protocol.csv")
    lot1_names = protocol.loc[protocol["lot"].eq(1), "candidate"].astype(str).tolist()

    scenario_summaries = []
    for scenario in SCENARIOS:
        print(f"[scenario] {scenario}", flush=True)
        fixed_groups: list[Group] = []
        candidate_groups: list[Group] = []
        scenario_designs = []
        for name in lot1_names:
            design = timing.scenario_design(designs_all[name], scenario)
            scenario_designs.append(design)
            fixed, candidates = build_groups(theta, design, scenario, step=0.04)
            fixed_groups.extend(fixed)
            candidate_groups.extend(candidates)
        fixed_fim = np.sum([group.fim for group in fixed_groups], axis=0)
        selected, final_fim = greedy_select(prior + fixed_fim, candidate_groups, objective="d_opt")
        selected.insert(0, "scenario", scenario)
        selected.to_csv(RESULTS_DIR / f"selected_samples_{scenario}.csv", index=False)
        plot_selection(scenario, selected, scenario_designs)
        combined = final_fim
        reductions = joint.variance_reduction(prior, combined, TARGET_PARAMETERS)
        summary = {
            "scenario": scenario,
            "n_selected": int(len(selected)),
            **joint.fim_metrics(combined, TARGET_PARAMETERS, prefix="combined_"),
            **reductions,
            **post_pulse_counts(selected, scenario_designs),
        }
        for names, label in [
            (timing.FOCUS_CORE, "focus_core"),
            (timing.FOCUS_SECONDARY, "focus_secondary"),
            (timing.FOCUS_AROMA, "focus_aroma"),
        ]:
            values = [float(reductions.get(f"var_reduction_{name}", np.nan)) for name in names]
            summary[f"{label}_mean_var_reduction"] = float(np.nanmean(values))
        scenario_summaries.append(summary)

    summary_df = pd.DataFrame(scenario_summaries).sort_values("combined_logdet", ascending=False)
    summary_df.to_csv(RESULTS_DIR / "express_optimal_sampling_summary.csv", index=False)
    report_cols = [
        "scenario",
        "n_selected",
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
    report = [
        "# Lot 1 express optimal sampling",
        "",
        "Greedy local FIM sampling refinement for Lot 1 after real inoculation at 2026-06-15 15:30.",
        "",
        "Budgets per 216 h fermentation: 14 `full` samples and 6 `small` samples. The initial extended `full` sample at t=0, online CO2, and terminal condensate are treated as fixed. Greedy selection fills the remaining liquid sample budget on the 10:00/12:00/14:00/16:00 operational grid with at most 4 liquid samples per fermentation per day.",
        "",
        "## Summary",
        "",
        summary_df[report_cols].to_markdown(index=False),
    ]
    (RESULTS_DIR / "express_optimal_sampling_report.md").write_text("\n".join(report), encoding="utf-8")
    print(summary_df[report_cols].to_string(index=False))


if __name__ == "__main__":
    main()
