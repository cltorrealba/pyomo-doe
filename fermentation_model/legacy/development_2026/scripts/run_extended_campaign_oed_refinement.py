from __future__ import annotations

import argparse
import hashlib
import contextlib
import io
import json
import logging
import math
import os
import sys
import textwrap
from dataclasses import replace
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_extended_campaign_doe as base
from run_fit_strategy_analysis import load_notebook_context

RESULTS_DIR = SCRIPT_DIR / "results" / "extended_campaign_oed_refinement"
TRIAL_FIM_DIR = RESULTS_DIR / "trial_fims"


def _matrix_from_csv(path: Path, parameters: tuple[str, ...]) -> np.ndarray:
    return pd.read_csv(path, index_col=0).loc[list(parameters), list(parameters)].to_numpy(dtype=float)


def _save_matrix(path: Path, matrix: np.ndarray, parameters: tuple[str, ...]) -> None:
    pd.DataFrame(matrix, index=parameters, columns=parameters).to_csv(path)


def _has_x_pulse(candidate: base.ExtendedDesignCandidate) -> bool:
    return any(float(amount) > 0.0 for _time, amount in candidate.pulses.get("X", tuple()))


def _design_payload(candidate: base.ExtendedDesignCandidate) -> dict:
    return {
        "name": candidate.name,
        "family": candidate.family,
        "horizon_h": round(float(candidate.horizon_h), 8),
        "initials": {key: round(float(value), 8) for key, value in sorted(candidate.initials.items())},
        "temperature_c": [round(float(value), 8) for value in candidate.temperature_c],
        "pulses": {
            channel: [[round(float(time), 8), round(float(amount), 8)] for time, amount in candidate.pulses.get(channel, tuple())]
            for channel in base.CHANNELS
        },
    }


def design_hash(candidate: base.ExtendedDesignCandidate) -> str:
    payload = json.dumps(_design_payload(candidate), sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:14]


def _safe_name(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in name)


def campaign_score(prior_fim: np.ndarray, campaign_fim: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    metrics = base.fim_metrics(campaign_fim, parameters, prefix="")
    reductions = base.variance_reduction(prior_fim, campaign_fim, parameters)
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-16)
    mean_correctable = float(reductions.get("correctable_mean_var_reduction", 0.0))
    worst_correctable = float(reductions.get("correctable_worst_var_reduction", 0.0))
    sn_reduction = float(reductions.get("var_reduction_sN", 0.0))
    # Mixed D/E/A style score: keep determinant pressure, but explicitly reward the weakest
    # practically correctable direction. This avoids choosing only high-logdet campaigns.
    score = (
        float(metrics["logdet"])
        + 10.0 * mean_correctable
        + 20.0 * worst_correctable
        + 5.0 * sn_reduction
        + 2.0 * math.log(min_rel)
    )
    return {
        "targeted_oed_score": float(score),
        "logdet": float(metrics["logdet"]),
        "min_relative_eigenvalue": float(metrics["min_relative_eigenvalue"]),
        "condition_number": float(metrics["condition_number"]),
        "correctable_mean_var_reduction": mean_correctable,
        "correctable_worst_var_reduction": worst_correctable,
        "sN_var_reduction": sn_reduction,
    }


def load_prior_fim(ns: dict, theta: dict[str, float], parameters: tuple[str, ...], step: float, prior_batches: tuple[str, ...]) -> np.ndarray:
    cached = base.RESULTS_DIR / "extended_historical_prior_fim.csv"
    if cached.exists():
        return _matrix_from_csv(cached, parameters)
    solver = base.make_solver(ns)
    return base.current_prior_fim(ns, theta, parameters, prior_batches, step, solver)


def load_candidate_fims(
    ns: dict,
    theta: dict[str, float],
    candidates: dict[str, base.ExtendedDesignCandidate],
    parameters: tuple[str, ...],
    step: float,
    pulse_slots: int,
    liquid_interval_h: float,
    co2_interval_h: float,
) -> dict[str, np.ndarray]:
    fims: dict[str, np.ndarray] = {}
    solver = base.make_solver(ns)
    for name, candidate in candidates.items():
        cached = base.RESULTS_DIR / f"fim_new_{name}.csv"
        if cached.exists():
            fims[name] = _matrix_from_csv(cached, parameters)
            continue
        print(f"[base FIM] recomputing missing candidate {name}", flush=True)
        experiment = base.FutureExtendedExperiment(ns, candidate, theta, parameters, pulse_slots, liquid_interval_h, co2_interval_h)
        fims[name] = base.fim_for_experiment(experiment, parameters, step, solver)
    return fims


def campaign_fim_for_names(names: list[str], fims: dict[str, np.ndarray], prior_fim: np.ndarray) -> np.ndarray:
    current = prior_fim.copy()
    for name in names:
        current = current + fims[name]
    return current


def select_targeted_campaign(
    candidates: dict[str, base.ExtendedDesignCandidate],
    fims: dict[str, np.ndarray],
    prior_fim: np.ndarray,
    parameters: tuple[str, ...],
    campaign_size: int,
    require_x_pulse: bool,
) -> tuple[list[str], pd.DataFrame, np.ndarray]:
    selected: list[str] = []
    current = prior_fim.copy()
    remaining = set(fims)
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best_name = None
        best_score = -np.inf
        best_components = None
        for name in sorted(remaining):
            trial = current + fims[name]
            components = campaign_score(prior_fim, trial, parameters)
            if components["targeted_oed_score"] > best_score:
                best_name = name
                best_score = components["targeted_oed_score"]
                best_components = components
        if best_name is None:
            break
        selected.append(best_name)
        current = current + fims[best_name]
        remaining.remove(best_name)
        c = candidates[best_name]
        rows.append(
            {
                "campaign_order": order,
                "candidate": best_name,
                "family": c.family,
                "horizon_h": c.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                "rationale": c.rationale,
                "selection_note": "targeted greedy",
                **(best_components or {}),
            }
        )

    if require_x_pulse and selected and not any(_has_x_pulse(candidates[name]) for name in selected):
        best = None
        x_candidates = [name for name in fims if _has_x_pulse(candidates[name])]
        for replacement in x_candidates:
            if replacement in selected:
                continue
            for idx, dropped in enumerate(selected):
                trial_names = selected.copy()
                trial_names[idx] = replacement
                if len(set(trial_names)) != len(trial_names):
                    continue
                trial_fim = campaign_fim_for_names(trial_names, fims, prior_fim)
                components = campaign_score(prior_fim, trial_fim, parameters)
                if best is None or components["targeted_oed_score"] > best[0]:
                    best = (components["targeted_oed_score"], dropped, replacement, trial_names, trial_fim)
        if best is not None:
            _score, dropped, replacement, selected, current = best
            rows = []
            running = prior_fim.copy()
            for order, name in enumerate(selected, start=1):
                running = running + fims[name]
                c = candidates[name]
                note = "targeted greedy"
                if name == replacement:
                    note = f"forced X-pulse coverage; replaced {dropped}"
                rows.append(
                    {
                        "campaign_order": order,
                        "candidate": name,
                        "family": c.family,
                        "horizon_h": c.horizon_h,
                        "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                        "rationale": c.rationale,
                        "selection_note": note,
                        **campaign_score(prior_fim, running, parameters),
                    }
                )
    return selected, pd.DataFrame(rows), current


def _quantize(value: float, step: float) -> float:
    if step <= 0:
        return float(value)
    return float(round(float(value) / step) * step)


def _project_pulses(pulses: dict[str, list[tuple[float, float]]], horizon_h: float, pulse_slots: int) -> dict[str, tuple[tuple[float, float], ...]]:
    projected: dict[str, tuple[tuple[float, float], ...]] = {}
    for channel in base.CHANNELS:
        lo, hi = base.CHANNEL_AMOUNT_BOUNDS[channel]
        merged: dict[float, float] = {}
        for time_h, amount in pulses.get(channel, []):
            amount = min(max(float(amount), lo), hi)
            if amount <= 1e-9:
                continue
            time_h = min(max(float(time_h), 0.0), float(horizon_h))
            merged[float(time_h)] = merged.get(float(time_h), 0.0) + float(amount)
        rows = list(merged.items())
        rows = sorted(rows, key=lambda item: item[0])[:pulse_slots]
        total = sum(amount for _time, amount in rows)
        limit = float(base.CHANNEL_TOTAL_LIMITS[channel])
        if total > limit and total > 0.0:
            rows = [(time_h, amount * limit / total) for time_h, amount in rows]
        projected[channel] = tuple(rows)
    return projected


def project_candidate(
    candidate: base.ExtendedDesignCandidate,
    pulse_slots: int,
    temp_step_c: float,
    time_step_h: float,
) -> base.ExtendedDesignCandidate:
    initials = {
        "X": min(max(float(candidate.initials.get("X", 0.6)), 0.05), 5.0),
        "N": min(max(float(candidate.initials.get("N", 0.18)), 0.015), 0.35),
        "G": min(max(float(candidate.initials.get("G", 80.0)), 5.0), 260.0),
        "F": min(max(float(candidate.initials.get("F", 80.0)), 5.0), 260.0),
        "E": min(max(float(candidate.initials.get("E", 0.0)), 0.0), 70.0),
        "Xd": 0.0,
    }
    sugar0 = initials["G"] + initials["F"]
    if sugar0 > 260.0:
        scale = 260.0 / sugar0
        initials["G"] *= scale
        initials["F"] *= scale
    if sugar0 < 40.0:
        scale = 40.0 / max(sugar0, 1e-8)
        initials["G"] = min(initials["G"] * scale, 260.0)
        initials["F"] = min(initials["F"] * scale, 260.0)

    temperature = tuple(
        min(max(_quantize(float(value), temp_step_c), 15.0), 25.0)
        for value in candidate.temperature_c
    )
    pulse_rows = {
        channel: [(_quantize(float(time_h), time_step_h), float(amount)) for time_h, amount in candidate.pulses.get(channel, tuple())]
        for channel in base.CHANNELS
    }
    pulses = _project_pulses(pulse_rows, candidate.horizon_h, pulse_slots)

    n_total = initials["N"] + sum(amount for _time, amount in pulses.get("N", tuple()))
    if n_total > 0.42:
        excess_scale = max((0.42 - initials["N"]) / max(sum(amount for _time, amount in pulses.get("N", tuple())), 1e-8), 0.0)
        pulses = {**pulses, "N": tuple((time_h, amount * excess_scale) for time_h, amount in pulses.get("N", tuple()))}

    e_total = initials["E"] + sum(amount for _time, amount in pulses.get("E", tuple()))
    if e_total > 75.0:
        excess_scale = max((75.0 - initials["E"]) / max(sum(amount for _time, amount in pulses.get("E", tuple())), 1e-8), 0.0)
        pulses = {**pulses, "E": tuple((time_h, amount * excess_scale) for time_h, amount in pulses.get("E", tuple()))}

    sugar_total = initials["G"] + initials["F"] + sum(amount for ch in ("G", "F") for _time, amount in pulses.get(ch, tuple()))
    if sugar_total > 330.0:
        added = sum(amount for ch in ("G", "F") for _time, amount in pulses.get(ch, tuple()))
        scale = max((330.0 - initials["G"] - initials["F"]) / max(added, 1e-8), 0.0)
        pulses = {
            **pulses,
            "G": tuple((time_h, amount * scale) for time_h, amount in pulses.get("G", tuple())),
            "F": tuple((time_h, amount * scale) for time_h, amount in pulses.get("F", tuple())),
        }

    x_total = initials["X"] + sum(amount for _time, amount in pulses.get("X", tuple()))
    if x_total > 5.0:
        scale = max((5.0 - initials["X"]) / max(sum(amount for _time, amount in pulses.get("X", tuple())), 1e-8), 0.0)
        pulses = {**pulses, "X": tuple((time_h, amount * scale) for time_h, amount in pulses.get("X", tuple()))}

    return replace(candidate, initials=initials, temperature_c=temperature, pulses=pulses)


def mutate_candidate(
    candidate: base.ExtendedDesignCandidate,
    rng: np.random.Generator,
    pulse_slots: int,
    allow_new_pulses: bool,
    temp_step_c: float,
    time_step_h: float,
) -> base.ExtendedDesignCandidate:
    temperature = [float(value) + rng.normal(0.0, 1.2) for value in candidate.temperature_c]
    initials = dict(candidate.initials)
    initials["X"] = float(initials.get("X", 0.6)) * math.exp(rng.normal(0.0, 0.18))
    initials["N"] = float(initials.get("N", 0.18)) + rng.normal(0.0, 0.025)
    initials["G"] = float(initials.get("G", 80.0)) * math.exp(rng.normal(0.0, 0.16))
    initials["F"] = float(initials.get("F", 80.0)) * math.exp(rng.normal(0.0, 0.16))
    initials["E"] = float(initials.get("E", 0.0)) + rng.normal(0.0, 5.0)
    initials["Xd"] = 0.0

    pulses = {channel: list(candidate.pulses.get(channel, tuple())) for channel in base.CHANNELS}
    for channel in base.CHANNELS:
        updated = []
        for time_h, amount in pulses[channel]:
            new_time = float(time_h) + rng.normal(0.0, 12.0)
            new_amount = float(amount) * math.exp(rng.normal(0.0, 0.28))
            updated.append((new_time, new_amount))
        pulses[channel] = updated

    if allow_new_pulses and rng.random() < 0.35:
        channel = str(rng.choice(list(base.CHANNELS)))
        if len(pulses[channel]) < pulse_slots:
            lo, hi = base.CHANNEL_AMOUNT_BOUNDS[channel]
            amount_scale = {
                "N": 0.055,
                "G": 30.0,
                "F": 30.0,
                "E": 18.0,
                "X": 1.0,
            }[channel]
            amount = min(max(rng.lognormal(math.log(max(amount_scale, 1e-6)), 0.35), lo), hi)
            time_h = float(rng.uniform(0.15 * candidate.horizon_h, 0.85 * candidate.horizon_h))
            pulses[channel].append((time_h, amount))

    mutated = replace(
        candidate,
        initials=initials,
        temperature_c=tuple(temperature),
        pulses={channel: tuple(rows) for channel, rows in pulses.items()},
    )
    return project_candidate(mutated, pulse_slots, temp_step_c, time_step_h)


def compute_candidate_fim(
    ns: dict,
    theta: dict[str, float],
    candidate: base.ExtendedDesignCandidate,
    parameters: tuple[str, ...],
    step: float,
    solver,
    pulse_slots: int,
    liquid_interval_h: float,
    co2_interval_h: float,
    cache_label: str,
) -> tuple[np.ndarray, str]:
    cache_dir = TRIAL_FIM_DIR / _safe_name(cache_label)
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = design_hash(candidate)
    path = cache_dir / f"fim_{key}.csv"
    if path.exists():
        return _matrix_from_csv(path, parameters), "cache"
    experiment = base.FutureExtendedExperiment(ns, candidate, theta, parameters, pulse_slots, liquid_interval_h, co2_interval_h)
    capture = io.StringIO()
    log_messages: list[str] = []

    class _CaptureHandler(logging.Handler):
        def emit(self, record):
            log_messages.append(self.format(record))

    handler = _CaptureHandler()
    handler.setLevel(logging.WARNING)
    root_logger = logging.getLogger()
    root_logger.addHandler(handler)
    try:
        with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
            fim = base.fim_for_experiment(experiment, parameters, step, solver)
    finally:
        root_logger.removeHandler(handler)
    warning_text = "\n".join([capture.getvalue(), *log_messages])
    if "maxIterations" in warning_text or "Maximum Number of Iterations" in warning_text:
        raise RuntimeError("Solver warning during FIM evaluation: Ipopt maxIterations.")
    _save_matrix(path, fim, parameters)
    return fim, "computed"


def _candidate_protocol_row(order: int, candidate: base.ExtendedDesignCandidate, source: str) -> dict:
    row = {
        "campaign_order": order,
        "candidate": candidate.name,
        "family": candidate.family,
        "horizon_h": candidate.horizon_h,
        "temperature_c": ", ".join(f"{value:g}" for value in candidate.temperature_c),
        "source": source,
        "rationale": candidate.rationale,
    }
    for key in ("X", "N", "G", "F", "E", "Xd"):
        row[f"{key}0"] = float(candidate.initials.get(key, 0.0))
    return row


def pulse_table(candidates: list[base.ExtendedDesignCandidate]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        for channel in base.CHANNELS:
            for time_h, amount in candidate.pulses.get(channel, tuple()):
                if float(amount) <= 0.0:
                    continue
                rows.append(
                    {
                        "candidate": candidate.name,
                        "family": candidate.family,
                        "channel": channel,
                        "time_h": float(time_h),
                        "amount": float(amount),
                    }
                )
    return pd.DataFrame(rows)


def write_input_plots(candidates: list[base.ExtendedDesignCandidate], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pulse_df = pulse_table(candidates)
    for stale in output_dir.glob("oed_refined_inputs_*.png"):
        stale.unlink()
    for candidate in candidates:
        fig, axes = plt.subplots(7, 1, figsize=(10, 12), sharex=False)
        time = np.linspace(0.0, candidate.horizon_h, 400)
        axes[0].step(time, base._temperature_profile(candidate, time), where="post", color="tab:red")
        axes[0].set_ylabel("T C")
        axes[0].set_xlim(0.0, candidate.horizon_h)
        init_labels = ["X", "N", "G", "F", "E"]
        init_values = [candidate.initials.get(name, 0.0) for name in init_labels]
        axes[1].bar(init_labels, init_values, color=["tab:purple", "tab:brown", "tab:blue", "tab:orange", "tab:green"])
        axes[1].set_ylabel("initial")
        for ax, channel in zip(axes[2:], base.CHANNELS):
            rows_ch = pulse_df[(pulse_df["candidate"].eq(candidate.name)) & (pulse_df["channel"].eq(channel))]
            if rows_ch.empty:
                ax.axhline(0.0, color="0.8", lw=1)
            else:
                ax.vlines(rows_ch["time_h"], 0.0, rows_ch["amount"], lw=3)
                ax.scatter(rows_ch["time_h"], rows_ch["amount"], s=35)
            ax.set_ylabel(channel)
            ax.set_xlim(0.0, candidate.horizon_h)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("t [h]")
        fig.suptitle("\n".join(textwrap.wrap(f"{candidate.name}: {candidate.rationale}", width=95)), fontsize=10, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(output_dir / f"oed_refined_inputs_{candidate.name}.png", dpi=170)
        plt.close(fig)


def write_report(
    base_selected: pd.DataFrame,
    targeted_selected: pd.DataFrame,
    refined_selected: pd.DataFrame,
    parameter_reduction: pd.DataFrame,
    trial_log: pd.DataFrame,
    direct_run_doe_note: str,
) -> None:
    lines = [
        "# Refinamiento OED de campana extendida",
        "",
        "Este reporte complementa el screening inicial. La campana inicial ya usaba `compute_FIM()` de Pyomo DoE, pero la seleccion era principalmente greedy sobre D-optimalidad/logdet con reduccion media de varianza.",
        "",
        "Aqui se usa un criterio mixto model-based:",
        "",
        "`score = logdet(FIM) + 10 mean(R_correctable) + 20 worst(R_correctable) + 5 R_sN + 2 log(lambda_min/lambda_max)`",
        "",
        "La razon tecnica es que el diagnostico extendido dejo `sN` como direccion practica mas debil. Por eso el score conserva informacion global, pero penaliza campanas que dejan un parametro corregible con poca reduccion de varianza.",
        "",
        "## Nota sobre `run_doe()` directo",
        "",
        direct_run_doe_note,
        "",
        "## Campana candidata original",
        "",
        base_selected.to_markdown(index=False),
        "",
        "## Campana re-seleccionada con criterio dirigido",
        "",
        targeted_selected.to_markdown(index=False),
        "",
        "## Campana refinada final",
        "",
        refined_selected.to_markdown(index=False),
        "",
        "## Reduccion de varianza por parametro",
        "",
        parameter_reduction.to_markdown(index=False),
        "",
        "## Ensayos de refinamiento local",
        "",
        trial_log.to_markdown(index=False) if not trial_log.empty else "No se ejecutaron ensayos locales (`--trials-per-candidate 0`).",
        "",
        "## Interpretacion",
        "",
        "- La inclusion de `yan_saturation_scan` es una decision deliberada: mejora la direccion `sN`, que era el peor caso practico.",
        "- Se mantiene al menos un pulso de biomasa viable (`X`) para conservar capacidad de separar escalamiento por biomasa, muerte y conversion a `Xd`.",
        "- El resultado no libera automaticamente los 15 parametros; indica que la campana propuesta deberia entregar mas informacion para perfil likelihood y recalibracion posterior.",
    ]
    (RESULTS_DIR / "oed_refinement_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Targeted/refined model-based OED for the extended fermentation campaign.")
    parser.add_argument("--prior-batches", nargs="*", default=["25026", "25086", "25150", "25170"])
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--step", type=float, default=1e-2)
    parser.add_argument("--pulse-slots", type=int, default=3)
    parser.add_argument("--liquid-interval-h", type=float, default=12.0)
    parser.add_argument("--co2-interval-h", type=float, default=6.0)
    parser.add_argument("--trials-per-candidate", type=int, default=3)
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20260607)
    parser.add_argument("--require-x-pulse", action="store_true", default=True)
    parser.add_argument("--allow-new-pulses", action="store_true", default=True)
    parser.add_argument("--temp-step-c", type=float, default=0.5)
    parser.add_argument("--time-step-h", type=float, default=1.0)
    parser.add_argument("--solver-max-iter", type=int, default=8000)
    parser.add_argument("--cache-label", default="max8000")
    parser.add_argument("--skip-simulations", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_notebook_context()
    theta = base.complete_final_theta(ns)
    solver = base.make_solver(ns)
    solver.options["max_iter"] = max(int(solver.options.get("max_iter", 0) or 0), int(args.solver_max_iter))
    parameters = base.FULL15
    rng = np.random.default_rng(int(args.seed))

    candidate_map = {candidate.name: candidate for candidate in base.candidate_library()}
    prior_fim = load_prior_fim(ns, theta, parameters, float(args.step), tuple(args.prior_batches))
    base_fims = load_candidate_fims(
        ns,
        theta,
        candidate_map,
        parameters,
        float(args.step),
        int(args.pulse_slots),
        float(args.liquid_interval_h),
        float(args.co2_interval_h),
    )

    base_selected_path = base.RESULTS_DIR / "extended_campaign_selected.csv"
    if base_selected_path.exists():
        base_selected = pd.read_csv(base_selected_path)
    else:
        names, base_selected, _fim = select_targeted_campaign(candidate_map, base_fims, prior_fim, parameters, args.campaign_size, args.require_x_pulse)
    base_selected.to_csv(RESULTS_DIR / "oed_reference_original_selected.csv", index=False)

    selected_names, targeted_selected, targeted_fim = select_targeted_campaign(
        candidate_map,
        base_fims,
        prior_fim,
        parameters,
        min(int(args.campaign_size), len(base_fims)),
        bool(args.require_x_pulse),
    )
    targeted_selected.to_csv(RESULTS_DIR / "oed_targeted_campaign_selected.csv", index=False)
    _save_matrix(RESULTS_DIR / "oed_targeted_campaign_fim.csv", targeted_fim, parameters)

    current_candidates = {name: candidate_map[name] for name in selected_names}
    current_fims = {name: base_fims[name] for name in selected_names}
    trial_rows = []

    print("[targeted selection]", ", ".join(selected_names), flush=True)
    for pass_idx in range(1, int(args.passes) + 1):
        print(f"[refinement pass {pass_idx}/{args.passes}]", flush=True)
        for name in list(selected_names):
            context_fim = prior_fim.copy()
            for other in selected_names:
                if other != name:
                    context_fim = context_fim + current_fims[other]
            baseline_fim = context_fim + current_fims[name]
            baseline_components = campaign_score(prior_fim, baseline_fim, parameters)
            best_candidate = current_candidates[name]
            best_fim = current_fims[name]
            best_components = baseline_components
            print(
                f"  [{name}] baseline score={baseline_components['targeted_oed_score']:.6g} "
                f"worst={baseline_components['correctable_worst_var_reduction']:.4f}",
                flush=True,
            )
            for trial_idx in range(1, int(args.trials_per_candidate) + 1):
                trial_candidate = mutate_candidate(
                    current_candidates[name],
                    rng,
                    int(args.pulse_slots),
                    bool(args.allow_new_pulses),
                    float(args.temp_step_c),
                    float(args.time_step_h),
                )
                try:
                    trial_fim, source = compute_candidate_fim(
                        ns,
                        theta,
                        trial_candidate,
                        parameters,
                        float(args.step),
                        solver,
                        int(args.pulse_slots),
                        float(args.liquid_interval_h),
                        float(args.co2_interval_h),
                        str(args.cache_label),
                    )
                    trial_campaign_fim = context_fim + trial_fim
                    components = campaign_score(prior_fim, trial_campaign_fim, parameters)
                    status = "ok"
                    error = ""
                except Exception as err:
                    trial_fim = None
                    source = "failed"
                    components = {
                        "targeted_oed_score": np.nan,
                        "logdet": np.nan,
                        "min_relative_eigenvalue": np.nan,
                        "condition_number": np.nan,
                        "correctable_mean_var_reduction": np.nan,
                        "correctable_worst_var_reduction": np.nan,
                        "sN_var_reduction": np.nan,
                    }
                    status = "failed"
                    error = f"{type(err).__name__}: {err}"
                accepted = bool(
                    status == "ok"
                    and components["targeted_oed_score"] > best_components["targeted_oed_score"] + 1e-6
                    and trial_fim is not None
                )
                trial_rows.append(
                    {
                        "pass": pass_idx,
                        "candidate": name,
                        "trial": trial_idx,
                        "status": status,
                        "fim_source": source,
                        "accepted_best_so_far": accepted,
                        "error": error,
                        "design_hash": design_hash(trial_candidate),
                        **components,
                    }
                )
                if accepted:
                    best_candidate = trial_candidate
                    best_fim = trial_fim
                    best_components = components
                    print(
                        f"    accepted trial {trial_idx}: score={components['targeted_oed_score']:.6g} "
                        f"worst={components['correctable_worst_var_reduction']:.4f}",
                        flush=True,
                    )
            if best_candidate is not current_candidates[name]:
                current_candidates[name] = best_candidate
                current_fims[name] = best_fim

    refined_fim = campaign_fim_for_names(selected_names, current_fims, prior_fim)
    refined_components = campaign_score(prior_fim, refined_fim, parameters)
    _save_matrix(RESULTS_DIR / "oed_refined_campaign_fim.csv", refined_fim, parameters)

    refined_rows = []
    running = prior_fim.copy()
    for order, name in enumerate(selected_names, start=1):
        running = running + current_fims[name]
        candidate = current_candidates[name]
        refined_rows.append(
            {
                "campaign_order": order,
                "candidate": name,
                "family": candidate.family,
                "horizon_h": candidate.horizon_h,
                "temperature_c": ", ".join(f"{value:g}" for value in candidate.temperature_c),
                "rationale": candidate.rationale,
                "has_X_pulse": _has_x_pulse(candidate),
                **campaign_score(prior_fim, running, parameters),
            }
        )
    refined_selected = pd.DataFrame(refined_rows)
    refined_selected.to_csv(RESULTS_DIR / "oed_refined_campaign_selected.csv", index=False)

    reductions = base.variance_reduction(prior_fim, refined_fim, parameters)
    parameter_rows = []
    for name in parameters:
        parameter_rows.append(
            {
                "parameter": name,
                "group": "accepted_current" if name in base.ACCEPTED7 else "currently_fixed_target",
                "campaign_var_ratio": reductions.get(f"var_ratio_{name}", np.nan),
                "campaign_var_reduction": reductions.get(f"var_reduction_{name}", np.nan),
            }
        )
    parameter_reduction = pd.DataFrame(parameter_rows)
    parameter_reduction.to_csv(RESULTS_DIR / "oed_refined_parameter_reduction.csv", index=False)

    trial_log = pd.DataFrame(trial_rows)
    trial_log.to_csv(RESULTS_DIR / "oed_refinement_trials.csv", index=False)

    refined_candidates = [current_candidates[name] for name in selected_names]
    pd.DataFrame(
        [_candidate_protocol_row(order, candidate, "targeted+local_refinement") for order, candidate in enumerate(refined_candidates, start=1)]
    ).to_csv(RESULTS_DIR / "oed_refined_protocols.csv", index=False)
    pulse_table(refined_candidates).to_csv(RESULTS_DIR / "oed_refined_pulses.csv", index=False)
    write_input_plots(refined_candidates, RESULTS_DIR)

    if not args.skip_simulations:
        sim_rows = []
        for candidate in refined_candidates:
            try:
                sim_rows.append(
                    base.simulate_candidate(ns, theta, candidate, int(args.pulse_slots), float(args.liquid_interval_h), float(args.co2_interval_h))
                )
            except Exception as err:
                print(f"[simulation failed] {candidate.name}: {type(err).__name__}: {err}", flush=True)
        if sim_rows:
            pd.concat(sim_rows, ignore_index=True).to_csv(RESULTS_DIR / "oed_refined_simulations.csv", index=False)

    direct_run_doe_note = (
        "Se probo `DesignOfExperiments.run_doe()` directamente sobre el modelo extendido con variables de diseno "
        "libres. En esta instalacion, Ipopt termino con condicion `other` y la version local de Pyomo DoE levanto "
        "`UnboundLocalError` al intentar reportar `results_message`. Por eso este refinamiento usa una ruta mas "
        "robusta: optimizacion externa acotada, con cada propuesta evaluada por `DesignOfExperiments.compute_FIM(method='sequential')`."
    )
    write_report(base_selected, targeted_selected, refined_selected, parameter_reduction, trial_log, direct_run_doe_note)

    metadata = {
        "prior_batches": list(args.prior_batches),
        "campaign_size": int(args.campaign_size),
        "step": float(args.step),
        "pulse_slots": int(args.pulse_slots),
        "liquid_interval_h": float(args.liquid_interval_h),
        "co2_interval_h": float(args.co2_interval_h),
        "trials_per_candidate": int(args.trials_per_candidate),
        "passes": int(args.passes),
        "seed": int(args.seed),
        "require_x_pulse": bool(args.require_x_pulse),
        "allow_new_pulses": bool(args.allow_new_pulses),
        "solver_max_iter": int(args.solver_max_iter),
        "cache_label": str(args.cache_label),
        "parameters": list(parameters),
        **refined_components,
    }
    (RESULTS_DIR / "oed_refinement_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("\nRefined campaign:")
    print(
        refined_selected[
            [
                "campaign_order",
                "candidate",
                "targeted_oed_score",
                "correctable_mean_var_reduction",
                "correctable_worst_var_reduction",
                "condition_number",
            ]
        ].to_string(index=False)
    )
    print("\nParameter reductions:")
    print(parameter_reduction.to_string(index=False))
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
