from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint
import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = SCRIPT_DIR / "results" / "final_operational_doe_v2"
NOTEBOOK_PATH = (
    SCRIPT_DIR
    / "laboratory_2026"
    / "notebooks"
    / "fermentation_final_operational_doe_v2.ipynb"
)

FERMENTATION_TARGETS = joint.FERMENTATION_TARGETS
SECONDARY_TARGETS = tuple(name for name in v2.V2_REDUCED_PARAMETERS if name != "kAldRed")
AROMA_TARGETS = joint.AROMA_TARGETS
TARGET_PARAMETERS = FERMENTATION_TARGETS + SECONDARY_TARGETS + AROMA_TARGETS
FIXED_SECONDARY = ("qO2", "kLaO2", "O2sat", "kAldRed", "kPyrS_stat", "kAldPyr", "kAldS_stat", "kAldO2", "kAcAssim")
SAMPLE_POLICIES = ("front_loaded", "balanced", "two_per_day")
START_DATES = ("2026-06-15 08:00", "2026-06-29 08:00", "2026-07-13 08:00")


@dataclass(frozen=True)
class SimV2:
    core: pd.DataFrame
    secondary: pd.DataFrame
    aromas_liq: pd.DataFrame
    aromas_loss: pd.DataFrame
    aromas_cond: pd.Series


def load_theta_final() -> dict[str, float]:
    theta = joint.load_reference_theta()
    for path in [
        joint.RESULTS_DIR / "theta_secondary_joint.csv",
        SCRIPT_DIR / "results" / "secondary_v2_model_evaluation" / "theta_secondary_v2_reduced_o2fixed.csv",
    ]:
        if path.exists():
            loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
            theta.update({str(k): float(v) for k, v in loaded.items() if np.isfinite(float(v))})
    theta.update(
        {
            "kPyrS_stat": v2.V2_BOUNDS["kPyrS_stat"][0],
            "kAldPyr": v2.V2_BOUNDS["kAldPyr"][0],
            "kAldS_stat": v2.V2_BOUNDS["kAldS_stat"][0],
            "kAldO2": v2.V2_BOUNDS["kAldO2"][0],
            "kAcAssim": v2.V2_BOUNDS["kAcAssim"][0],
        }
    )
    out = dict(theta)
    for name, (lb, ub) in v2.V2_BOUNDS.items():
        if name in out:
            out[name] = float(np.clip(out[name], lb, ub))
    return joint.clip_all(out)


def parameter_bounds(name: str) -> tuple[float, float]:
    if name in v2.V2_BOUNDS:
        return v2.V2_BOUNDS[name]
    return joint.PARAMETER_BOUNDS[name]


def sample_times(horizon_h: float, policy: str) -> np.ndarray:
    return base.operational_sample_times(float(horizon_h), policy=policy)


def o2_times(horizon_h: float, policy: str) -> np.ndarray:
    op = sample_times(horizon_h, policy="balanced")
    early = op[op <= min(float(horizon_h), 72.0)]
    return np.unique(np.concatenate([np.array([0.0], dtype=float), early]))


def co2_times(horizon_h: float) -> np.ndarray:
    return np.arange(0.0, float(horizon_h) + 1e-9, 6.0, dtype=float)


def dense_time_grid(design: base.FutureDesign, policy: str) -> np.ndarray:
    times = set(np.arange(0.0, float(design.horizon_h) + 1e-9, 2.0).round(8))
    times.update(float(t) for t in sample_times(design.horizon_h, policy))
    times.update(float(t) for t in o2_times(design.horizon_h, policy))
    times.update(float(t) for t in co2_times(design.horizon_h))
    times.add(float(design.horizon_h))
    for schedule in design.pulses.values():
        for time_h, _amount in schedule:
            times.add(round(max(0.0, float(time_h) - 2.0), 8))
            times.add(round(float(time_h), 8))
            times.add(round(min(float(design.horizon_h), float(time_h) + 2.0), 8))
    return np.asarray(sorted(t for t in times if 0.0 <= t <= float(design.horizon_h)), dtype=float)


def simulate_v2(design: base.FutureDesign, theta: dict[str, float], policy: str) -> SimV2 | None:
    times = dense_time_grid(design, policy)
    core = base.simulate(design, theta, times)
    if core is None:
        return None
    core = core.sort_index()
    sec = v2.integrate_secondary_v2(design, theta, core)
    try:
        liq, loss, cond = joint.integrate_aroma_euler(design, theta, core)
    except Exception:
        return None
    return SimV2(core=core, secondary=sec, aromas_liq=liq, aromas_loss=loss, aromas_cond=cond)


def future_residual_vector(theta: dict[str, float], design: base.FutureDesign, policy: str) -> np.ndarray:
    sim = simulate_v2(design, theta, policy)
    if sim is None:
        return np.ones(1000, dtype=float) * 1e6
    residuals: list[np.ndarray] = []
    liquid = [t for t in sample_times(design.horizon_h, policy) if t in sim.core.index]
    core = sim.core.loc[liquid]
    for state in base.STATE_NAMES:
        center = core[state].to_numpy(dtype=float)
        residuals.append(center / base.sigma_for_state(state, center))
    for state in ("Pyr", "AcAld", "Acetate"):
        center = sim.secondary.loc[liquid, state].to_numpy(dtype=float)
        residuals.append(center / joint.SIGMA[state])
    oxy = [t for t in o2_times(design.horizon_h, policy) if t in sim.secondary.index]
    if oxy:
        residuals.append(sim.secondary.loc[oxy, "O2"].to_numpy(dtype=float) / joint.SIGMA["O2"])
    ctimes = [t for t in co2_times(design.horizon_h) if t in sim.secondary.index]
    if ctimes:
        residuals.append(sim.secondary.loc[ctimes, "CO2"].to_numpy(dtype=float) / joint.SIGMA["CO2"])
    liq_aroma = sim.aromas_liq.loc[liquid]
    for species in joint.AROMA_SPECIES:
        residuals.append(liq_aroma[species].to_numpy(dtype=float) / joint.SIGMA[f"{species}_liq"])
        residuals.append(np.array([sim.aromas_cond[species] / joint.SIGMA[f"{species}_cond"]], dtype=float))
    return np.concatenate(residuals)


def finite_difference_fim(theta: dict[str, float], parameters: tuple[str, ...], residual_fun, step: float) -> np.ndarray:
    base_res = residual_fun(theta)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        low, high = parameter_bounds(name)
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), low, high))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), low, high))
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
    jac = np.column_stack(cols)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def parameter_estimability(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...], label: str) -> pd.DataFrame:
    cov = joint.stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[idx, idx], 0.0)))
        value = float(theta[name])
        low, high = parameter_bounds(name)
        active = value <= low * 1.01 or value >= high / 1.01
        if std_log <= 0.35 and not active:
            cls = "well_estimated"
        elif std_log <= 0.75 and not active:
            cls = "moderate"
        elif std_log <= 1.25 and not active:
            cls = "weak_but_actionable"
        else:
            cls = "weak_or_confounded"
        rows.append(
            {
                "analysis": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std_log,
                "approx_95_multiplier": float(math.exp(1.96 * min(std_log, 20.0))),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return pd.DataFrame(rows)


def current_prior_fim(theta: dict[str, float], batches: list[base.BatchData], step: float) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    ferment_batches = [b for b in batches if b.medium in {"natural", "synthetic"}]
    jac_f, _ = base.build_jacobian(theta, FERMENTATION_TARGETS, ferment_batches, step=step)
    fim_f = jac_f.T @ jac_f
    prior = np.eye(len(TARGET_PARAMETERS), dtype=float) * 1e-6
    for i, name_i in enumerate(FERMENTATION_TARGETS):
        ii = TARGET_PARAMETERS.index(name_i)
        for j, name_j in enumerate(FERMENTATION_TARGETS):
            jj = TARGET_PARAMETERS.index(name_j)
            prior[ii, jj] += fim_f[i, j]

    core_cache = joint.precompute_core_cache(batches, theta)
    fim_s = finite_difference_fim(
        theta,
        SECONDARY_TARGETS,
        lambda th: v2.residual_v2(th, batches, core_cache),
        step,
    )
    for i, name_i in enumerate(SECONDARY_TARGETS):
        ii = TARGET_PARAMETERS.index(name_i)
        for j, name_j in enumerate(SECONDARY_TARGETS):
            jj = TARGET_PARAMETERS.index(name_j)
            prior[ii, jj] += fim_s[i, j]

    for name in AROMA_TARGETS:
        idx = TARGET_PARAMETERS.index(name)
        cv = 0.75 if name.startswith("alpha_") else 1.25
        prior[idx, idx] += 1.0 / (cv * cv)

    prior = 0.5 * (prior + prior.T)
    estim = parameter_estimability(prior, theta, TARGET_PARAMETERS, "current_prior_v2")
    spectrum, weak = joint.eigen_diagnostics(prior, TARGET_PARAMETERS, "current_prior_v2")
    return prior, estim, pd.concat([spectrum.assign(kind="spectrum"), weak.assign(kind="weak_loading")], ignore_index=True, sort=False)


def secondary_initials(data: pd.DataFrame) -> dict[str, dict[str, float]]:
    return joint.secondary_initials_by_medium(data)


def pulse_builder(horizon: float, **kwargs) -> dict[str, tuple[tuple[float, float], ...]]:
    schedules = {channel: [] for channel in base.INPUT_CHANNELS}
    for channel, rows in kwargs.items():
        for target, amount in rows:
            schedules[channel].append((base.nearest_operational_time(float(target), float(horizon)), float(amount)))
    return {channel: tuple(sorted(rows)) for channel, rows in schedules.items()}


def add_literature_candidates(data: pd.DataFrame, designs: list[base.FutureDesign]) -> list[base.FutureDesign]:
    by_name = {design.name: design for design in designs}
    synthetic = dict(by_name["synthetic_high_sugar_reference"].initials)
    natural = dict(by_name["natural_control_18C"].initials)
    extras = [
        base.FutureDesign(
            "synthetic_lit_SM70_18C_lowN_phase",
            "literature_nitrogen_low",
            "synthetic",
            216.0,
            {**synthetic, "N": 0.070, "G": 115.0, "F": 115.0, "E": 0.0},
            (18.0, 18.0, 18.0, 18.0),
            pulse_builder(216.0),
            "Literature SM70-like low YAN level tests nitrogen-limited pyruvate/aldehyde and aroma phase behavior.",
        ),
        base.FutureDesign(
            "synthetic_lit_SM230_18C_midN_reference",
            "literature_nitrogen_mid",
            "synthetic",
            216.0,
            {**synthetic, "N": 0.230, "G": 115.0, "F": 115.0, "E": 0.0},
            (18.0, 18.0, 18.0, 18.0),
            pulse_builder(216.0),
            "Literature SM230-like mid YAN reference near reported higher-alcohol optimum.",
        ),
        base.FutureDesign(
            "synthetic_lit_SM410_18C_highN_ester",
            "literature_nitrogen_high",
            "synthetic",
            216.0,
            {**synthetic, "N": 0.410, "G": 115.0, "F": 115.0, "E": 0.0},
            (18.0, 18.0, 18.0, 18.0),
            pulse_builder(216.0),
            "Literature SM410-like high YAN tests ester capacity and nitrogen-saturated aroma phase.",
        ),
        base.FutureDesign(
            "synthetic_lit_SM410_24C_highN_strip",
            "literature_temp_strip",
            "synthetic",
            168.0,
            {**synthetic, "N": 0.410, "G": 125.0, "F": 125.0, "E": 0.0},
            (20.0, 24.0, 24.0, 20.0),
            pulse_builder(168.0),
            "High YAN plus warm high-rate fermentation excites CO2 stripping and temperature-partition effects.",
        ),
        base.FutureDesign(
            "synthetic_lit_SM70_warm_shift_Nrescue",
            "literature_nitrogen_timing",
            "synthetic",
            216.0,
            {**synthetic, "N": 0.070, "G": 115.0, "F": 115.0, "E": 0.0},
            (16.0, 18.0, 24.0, 20.0),
            pulse_builder(216.0, N=((50.0, 0.060),)),
            "Low initial YAN with operational N rescue tests phase transition timing and practical nutrient intervention.",
        ),
        base.FutureDesign(
            "natural_lit_midN_24C_matrix_transfer",
            "natural_literature_transfer",
            "natural",
            216.0,
            {**natural, "N": max(float(natural["N"]), 0.230), "E": 0.0},
            (18.0, 20.0, 24.0, 20.0),
            pulse_builder(216.0),
            "Natural must at literature-like mid/high YAN and warm phase checks synthetic-to-natural transfer.",
        ),
    ]
    return designs + [design for design in extras if design.name not in by_name]


def candidate_library(data: pd.DataFrame) -> dict[str, base.FutureDesign]:
    sec_init = secondary_initials(data)
    designs = joint.extend_candidate_set(data)
    patched = []
    for design in designs:
        patched.append(
            base.FutureDesign(
                design.name,
                design.family,
                design.medium,
                design.horizon_h,
                {**dict(design.initials), **sec_init.get(design.medium, {})},
                design.temperature_segments,
                design.pulses,
                design.rationale,
            )
        )
    patched = add_literature_candidates(data, patched)
    out = {design.name: design for design in patched}
    return out


def select_campaign(fims: dict[str, np.ndarray], designs: dict[str, base.FutureDesign], prior: np.ndarray, objective: str, size: int) -> tuple[pd.DataFrame, np.ndarray]:
    return joint.greedy_select(fims, designs, prior, TARGET_PARAMETERS, size, objective)


def dopt_exchange(seed: list[str], fims: dict[str, np.ndarray], prior: np.ndarray) -> list[str]:
    return joint.dopt_exchange(seed, fims, prior, TARGET_PARAMETERS)


def selection_frame(names: list[str], designs: dict[str, base.FutureDesign], prior: np.ndarray, fims: dict[str, np.ndarray], label: str) -> tuple[pd.DataFrame, np.ndarray]:
    return joint.selection_frame(names, designs, prior, fims, TARGET_PARAMETERS, label)


def evaluate_policy(theta: dict[str, float], prior: np.ndarray, designs: dict[str, base.FutureDesign], policy: str, step: float, campaign_size: int) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray], pd.DataFrame, np.ndarray]:
    rows = []
    fims: dict[str, np.ndarray] = {}
    for idx, design in enumerate(designs.values(), start=1):
        print(f"[{policy}] candidate {idx}/{len(designs)} {design.name}", flush=True)
        try:
            fim = finite_difference_fim(theta, TARGET_PARAMETERS, lambda th, d=design: future_residual_vector(th, d, policy), step)
            fims[design.name] = fim
            pd.DataFrame(fim, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / f"candidate_fim_{policy}_{design.name}.csv")
            combined = prior + fim
            rows.append(
                {
                    "policy": policy,
                    "candidate": design.name,
                    "family": design.family,
                    "medium": design.medium,
                    "horizon_h": design.horizon_h,
                    "n_liquid_samples": int(len(sample_times(design.horizon_h, policy))),
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
    selected_hybrid, fim_hybrid = select_campaign(fims, designs, prior, "hybrid", campaign_size)
    selected_hybrid.insert(0, "policy", policy)
    selected_dopt, fim_dopt = select_campaign(fims, designs, prior, "d_opt", campaign_size)
    selected_dopt.insert(0, "policy", policy)
    exchange_names = dopt_exchange(selected_hybrid["candidate"].astype(str).tolist(), fims, prior)
    selected_exchange, fim_exchange = selection_frame(exchange_names, designs, prior, fims, f"dopt_exchange_from_hybrid_{policy}")
    selected_exchange.insert(0, "policy", policy)
    selected = pd.concat([selected_hybrid, selected_dopt, selected_exchange], ignore_index=True)
    return ranking, selected, fims, selected_hybrid, fim_hybrid


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
        return "front_loaded"
    return str(summary.iloc[0]["policy"])


def selected_protocol(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], policy: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    sample_rows = []
    hybrid = selected[selected["objective"].eq("hybrid")].sort_values("campaign_order").reset_index(drop=True)
    for idx, record in hybrid.iterrows():
        design = designs[str(record["candidate"])]
        order = int(record["campaign_order"])
        lot = int((order - 1) // 3 + 1)
        start_dt = datetime.fromisoformat(START_DATES[lot - 1])
        rows.append(
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
                "rationale": design.rationale,
            }
        )
        for t in sample_times(design.horizon_h, policy):
            dt = start_dt + timedelta(hours=float(t))
            sample_rows.append(
                {
                    "campaign_order": order,
                    "lot": lot,
                    "candidate": design.name,
                    "event": "liquid_sample",
                    "relative_time_h": float(t),
                    "datetime": dt.isoformat(sep=" "),
                    "weekday": dt.strftime("%A"),
                    "details": "G/F/YAN/PAN/NH4/glycerol/ethanol/biomass/Pyr/AcAld/Acetate/aromas",
                    "channel": "",
                    "amount": np.nan,
                    "operational": bool(dt.weekday() < 5 and 9 <= dt.hour <= 16),
                }
            )
        for t in o2_times(design.horizon_h, policy):
            dt = start_dt + timedelta(hours=float(t))
            sample_rows.append(
                {
                    "campaign_order": order,
                    "lot": lot,
                    "candidate": design.name,
                    "event": "DO_measurement",
                    "relative_time_h": float(t),
                    "datetime": dt.isoformat(sep=" "),
                    "weekday": dt.strftime("%A"),
                    "details": "DO spot measurement; do not treat as controlled oxygen pulse",
                    "channel": "",
                    "amount": np.nan,
                    "operational": bool(t == 0.0 or (dt.weekday() < 5 and 9 <= dt.hour <= 16)),
                }
            )
        for channel, schedule in design.pulses.items():
            for t, amount in schedule:
                dt = start_dt + timedelta(hours=float(t))
                sample_rows.append(
                    {
                        "campaign_order": order,
                        "lot": lot,
                        "candidate": design.name,
                        "event": "manual_pulse",
                        "relative_time_h": float(t),
                        "datetime": dt.isoformat(sep=" "),
                        "weekday": dt.strftime("%A"),
                        "details": "manual addition inside working window",
                        "channel": channel,
                        "amount": float(amount),
                        "operational": bool(dt.weekday() < 5 and 9 <= dt.hour <= 16),
                    }
                )
        end_dt = start_dt + timedelta(hours=float(design.horizon_h))
        retrieval = end_dt
        if not (retrieval.weekday() < 5 and 9 <= retrieval.hour <= 16):
            # Condensate can be retained in the trap; retrieve at the first practical slot after nominal horizon.
            cursor = retrieval
            while not (cursor.weekday() < 5 and 10 <= cursor.hour <= 16):
                cursor += timedelta(hours=1)
            retrieval = cursor
        sample_rows.append(
            {
                "campaign_order": order,
                "lot": lot,
                "candidate": design.name,
                "event": "terminal_condensate_retrieval",
                "relative_time_h": float((retrieval - start_dt).total_seconds() / 3600.0),
                "datetime": retrieval.isoformat(sep=" "),
                "weekday": retrieval.strftime("%A"),
                "details": "single terminal condensate GC sample; CO2 online remains continuous",
                "channel": "",
                "amount": np.nan,
                "operational": True,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(sample_rows).sort_values(["lot", "campaign_order", "relative_time_h", "event"]).reset_index(drop=True)


def plot_inputs(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], policy: str, theta: dict[str, float]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    hybrid = selected[selected["objective"].eq("hybrid")].sort_values("campaign_order")
    for name in hybrid["candidate"].astype(str):
        design = designs[name]
        time = np.linspace(0.0, float(design.horizon_h), 350)
        sim = simulate_v2(design, theta, policy)
        fig, axes = plt.subplots(5, 1, figsize=(11, 12), sharex=True)
        axes[0].step(time, [base.temperature_at(design, t) for t in time], where="post", color="tab:red")
        axes[0].set_ylabel("T [C]")
        for channel in ("N", "G", "F", "E", "X"):
            schedule = design.pulses.get(channel, tuple())
            if schedule:
                axes[1].scatter([t for t, _ in schedule], [a for _, a in schedule], label=channel)
        axes[1].set_ylabel("Pulse amount")
        axes[1].legend(fontsize=8)
        axes[2].vlines(sample_times(design.horizon_h, policy), 0.0, 1.0, color="tab:blue", label="liquid")
        axes[2].vlines(o2_times(design.horizon_h, policy), 1.1, 1.45, color="tab:green", label="DO")
        axes[2].set_ylim(0, 1.6)
        axes[2].set_yticks([])
        axes[2].legend(fontsize=8)
        if sim is not None:
            axes[3].plot(sim.core.index, sim.core["G"], label="G")
            axes[3].plot(sim.core.index, sim.core["F"], label="F")
            axes[3].plot(sim.core.index, sim.core["N"] * 1000.0, label="YAN mg/L")
            axes[3].set_ylabel("Core")
            axes[3].legend(fontsize=8)
            for state in ("Pyr", "AcAld", "Acetate", "O2"):
                axes[4].plot(sim.secondary.index, sim.secondary[state], label=state)
            axes[4].legend(fontsize=8)
            axes[4].set_ylabel("Secondary")
        for ax in axes:
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("time [h]")
        fig.suptitle(name)
        fig.tight_layout()
        fig.savefig(plot_dir / f"final_operational_inputs_{name}.png", dpi=170)
        plt.close(fig)


def write_report(
    theta: dict[str, float],
    prior_estim: pd.DataFrame,
    policy_summary_df: pd.DataFrame,
    ranking: pd.DataFrame,
    selected: pd.DataFrame,
    campaign_table: pd.DataFrame,
    schedule: pd.DataFrame,
) -> None:
    weak = prior_estim[prior_estim["classification"].eq("weak_or_confounded")]
    lines = [
        "# Final operational DOE with secondary v2 reduced model",
        "",
        "## Literature-driven modeling decisions",
        "",
        "Mouret et al. show that nitrogen and temperature are the dominant levers for fermentative aroma kinetics, and that gas-liquid balances are needed to avoid confusing biological synthesis with physical evaporation. The final DOE therefore spans low, mid, and high YAN levels, includes temperature setpoint changes, uses CO2 online, and includes terminal condensate.",
        "",
        "Morakul et al. support treating gas-liquid partitioning as composition/temperature dependent rather than estimating arbitrary loss constants from liquid data alone. The DOE keeps UNIFAC/partition parameters literature-informed and only estimates low-dimensional aroma synthesis/loss scale parameters.",
        "",
        "Henriques/Scott dynamic genome-scale studies motivate a mechanistic interpretation for pyruvate, acetaldehyde, acetate, acetyl-CoA/redox, and ester formation. However, our fit-capacity diagnostic showed that several detailed redox/acetate-drain terms are not transferable or identifiable. The final DOE therefore uses a reduced empirical-mechanistic secondary layer.",
        "",
        "## Final target parameters",
        "",
        ", ".join(TARGET_PARAMETERS),
        "",
        "Fixed or strongly regularized secondary terms:",
        "",
        ", ".join(FIXED_SECONDARY),
        "",
        "## Current weak/confounded parameters before new campaign",
        "",
        weak.to_markdown(index=False) if not weak.empty else "No weak/confounded parameters under the current local FIM.",
        "",
        "## Sampling-policy benchmark",
        "",
        policy_summary_df.to_markdown(index=False),
        "",
        "## Selected campaign",
        "",
        selected[selected["objective"].eq("hybrid")].to_markdown(index=False),
        "",
        "## Campaign protocol",
        "",
        campaign_table.to_markdown(index=False),
        "",
        "## Operational schedule preview",
        "",
        schedule.head(80).to_markdown(index=False),
        "",
        "## Top candidate ranking for selected policy",
        "",
        ranking.sort_values("hybrid_score", ascending=False).head(15).to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- This is a model-based DOE over the reduced v2 secondary layer, core fermentation states, CO2, and aroma liquid/condensate outputs.",
        "- The campaign is not asking DOE to identify the terms that fit-capacity already showed as non-transferable.",
        "- Synthetic experiments dominate the information matrix because composition is controllable; natural experiments are retained to test matrix transfer before MPCC use in wine must.",
        "- All liquid samples and manual pulses are scheduled in operational windows; CO2 is online; terminal condensate retrieval may occur after the nominal horizon because the trap integrates volatile loss.",
    ]
    (RESULTS_DIR / "final_operational_doe_v2_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_notebook() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Final operational DOE with secondary v2 reduced model\n\n"
            "This notebook summarizes the literature-informed and operationally constrained model-based DOE for the next fermentation campaign."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "ROOT = Path.cwd()\n"
            "if not (ROOT / 'results').exists() and (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'results/final_operational_doe_v2'\n"
            "display(pd.read_csv(RESULTS / 'policy_summary.csv'))\n"
            "display(pd.read_csv(RESULTS / 'selected_campaign.csv'))\n"
            "display(pd.read_csv(RESULTS / 'campaign_protocol.csv'))"
        ),
        nbformat.v4.new_markdown_cell("## Current Estimability"),
        nbformat.v4.new_code_cell("display(pd.read_csv(RESULTS / 'current_prior_estimability.csv'))"),
        nbformat.v4.new_markdown_cell("## Operational Schedule"),
        nbformat.v4.new_code_cell("display(pd.read_csv(RESULTS / 'operational_schedule.csv').head(120))"),
        nbformat.v4.new_markdown_cell("## Input and Prediction Plots"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('final_operational_inputs_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Report"),
        nbformat.v4.new_code_cell("print((RESULTS / 'final_operational_doe_v2_report.md').read_text(encoding='utf-8'))"),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> int:
    parser = argparse.ArgumentParser(description="Final literature-informed operational DOE using secondary v2 reduced model.")
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--step", type=float, default=1e-2)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    theta = load_theta_final()
    pd.Series(theta, name="theta_final_operational_v2").to_csv(RESULTS_DIR / "theta_final_operational_v2.csv")
    data = joint.load_current_review_data()
    batches = joint.make_secondary_batches(data)
    designs = candidate_library(data)
    (RESULTS_DIR / "design_library.json").write_text(
        json.dumps({name: {"family": d.family, "medium": d.medium, "horizon_h": d.horizon_h, "rationale": d.rationale} for name, d in designs.items()}, indent=2),
        encoding="utf-8",
    )
    prior, prior_estim, prior_diag = current_prior_fim(theta, batches, args.step)
    pd.DataFrame(prior, index=TARGET_PARAMETERS, columns=TARGET_PARAMETERS).to_csv(RESULTS_DIR / "current_prior_fim.csv")
    prior_estim.to_csv(RESULTS_DIR / "current_prior_estimability.csv", index=False)
    prior_diag.to_csv(RESULTS_DIR / "current_prior_eigen_diagnostics.csv", index=False)

    all_rankings = []
    all_selected = []
    selected_by_policy: dict[str, pd.DataFrame] = {}
    fims_by_policy: dict[str, dict[str, np.ndarray]] = {}
    for policy in SAMPLE_POLICIES:
        ranking, selected, fims, _hybrid, _fim_hybrid = evaluate_policy(theta, prior, designs, policy, args.step, args.campaign_size)
        ranking.to_csv(RESULTS_DIR / f"candidate_ranking_{policy}.csv", index=False)
        selected.to_csv(RESULTS_DIR / f"selected_campaign_{policy}.csv", index=False)
        all_rankings.append(ranking)
        all_selected.append(selected)
        selected_by_policy[policy] = selected
        fims_by_policy[policy] = fims

    ranking_all = pd.concat(all_rankings, ignore_index=True)
    selected_all = pd.concat(all_selected, ignore_index=True)
    ranking_all.to_csv(RESULTS_DIR / "candidate_ranking_all_policies.csv", index=False)
    selected_all.to_csv(RESULTS_DIR / "selected_campaign_all_policies.csv", index=False)
    summary = policy_summary(selected_by_policy)
    summary.to_csv(RESULTS_DIR / "policy_summary.csv", index=False)
    chosen_policy = choose_policy(summary)
    selected = selected_by_policy[chosen_policy]
    selected.to_csv(RESULTS_DIR / "selected_campaign.csv", index=False)
    ranking = ranking_all[ranking_all["policy"].eq(chosen_policy)].copy()
    campaign_table, schedule = selected_protocol(selected, designs, chosen_policy)
    campaign_table.to_csv(RESULTS_DIR / "campaign_protocol.csv", index=False)
    schedule.to_csv(RESULTS_DIR / "operational_schedule.csv", index=False)
    plot_inputs(selected, designs, chosen_policy, theta)

    metadata = {
        "chosen_policy": chosen_policy,
        "campaign_size": int(args.campaign_size),
        "target_parameters": list(TARGET_PARAMETERS),
        "fixed_secondary_parameters": list(FIXED_SECONDARY),
        "start_dates": list(START_DATES),
        "literature_basis": [
            "Mouret et al. 2014 Food Research International: nitrogen-temperature aroma kinetics and gas-liquid balances",
            "Morakul et al. 2011/2013: gas-liquid partitioning and volatile losses during fermentation",
            "Henriques/Scott dynamic genome-scale models: mechanistic context for acetate, pyruvate, acetaldehyde, and ester formation",
        ],
    }
    (RESULTS_DIR / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_report(theta, prior_estim, summary, ranking, selected, campaign_table, schedule)
    write_notebook()

    print(f"Chosen sampling policy: {chosen_policy}")
    print("\nPolicy summary:")
    print(summary.to_string(index=False))
    print("\nSelected campaign:")
    print(selected[selected["objective"].eq("hybrid")][["campaign_order", "candidate", "medium", "campaign_logdet", "new_param_worst_var_reduction"]].to_string(index=False))
    print(f"\n[done] results={RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
