from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_aroma_calibration import (  # noqa: E402
    load_partition_surrogates,
    simulate_aroma,
)
from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    DesignPolicy,
    _future_design,
    allowed_sampling_times,
    fim_from_prepared,
    load_json,
    logdet,
    prepare_design,
    prior_precision,
    representative_members,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    create_immutable_run_directory,
    write_json,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def _load_policies(path: Path) -> tuple[DesignPolicy, ...]:
    actions = pd.read_csv(path)
    policies = []
    for name, group in actions.groupby("policy", sort=False):
        temperature = tuple(
            group[group["action"].eq("temperature_setpoint")]
            .sort_values("time_h")["value"]
            .astype(float)
        )
        nutrition = tuple(
            (float(row.time_h), float(row.value))
            for row in group[group["action"].eq("nutrition")].itertuples()
        )
        policies.append(DesignPolicy(str(name), temperature, nutrition))
    if len(policies) != 3:
        raise ValueError(f"Expected three Wave-1 policies; found {len(policies)}")
    return tuple(policies)


def _score(
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    prepared: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
) -> tuple[float, np.ndarray]:
    base = logdet(prior)
    gains = []
    for member in members:
        total = prior.copy()
        for policy in policies:
            total += fim_from_prepared(
                prepared[(member, policy.name)],
                np.asarray(schedules[policy.name], dtype=float),
                config,
            )
        gains.append(logdet(total) - base)
    gains = np.asarray(gains, dtype=float)
    score = (
        float(config["objective"]["median_weight"]) * float(np.median(gains))
        + float(config["objective"]["lower_decile_weight"])
        * float(np.quantile(gains, 0.1))
    )
    return score, gains


def _individual_score(
    schedule: tuple[float, ...],
    policy: DesignPolicy,
    prepared: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
) -> float:
    base = logdet(prior)
    gains = []
    for member in members:
        total = prior + fim_from_prepared(
            prepared[(member, policy.name)], np.asarray(schedule), config
        )
        gains.append(logdet(total) - base)
    return float(np.median(gains)) + 0.25 * float(np.quantile(gains, 0.1))


def _three_anchor_reference(
    schedule: tuple[float, ...],
    anchor: DesignPolicy,
    prepared: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
) -> tuple[float, np.ndarray]:
    base = logdet(prior)
    gains = []
    for member in members:
        fim = fim_from_prepared(
            prepared[(member, anchor.name)], np.asarray(schedule, dtype=float), config
        )
        gains.append(logdet(prior + 3.0 * fim) - base)
    gains = np.asarray(gains, dtype=float)
    score = (
        float(config["objective"]["median_weight"]) * float(np.median(gains))
        + float(config["objective"]["lower_decile_weight"])
        * float(np.quantile(gains, 0.1))
    )
    return score, gains


def _greedy_schedule(
    policy: DesignPolicy,
    valid: np.ndarray,
    prepared: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
) -> tuple[float, ...]:
    count = int(config["sampling"]["samples_per_process"])
    selected = {float(valid[0]), float(valid[-1])}
    while len(selected) < count:
        best_time, best_score = None, -np.inf
        for candidate in valid:
            if float(candidate) in selected:
                continue
            schedule = tuple(sorted(selected | {float(candidate)}))
            score = _individual_score(schedule, policy, prepared, members, prior, config)
            if score > best_score:
                best_time, best_score = float(candidate), score
        if best_time is None:
            raise RuntimeError(f"Not enough valid sample slots for {policy.name}")
        selected.add(best_time)
    return tuple(sorted(selected))


def _coordinate_improve(
    schedules: dict[str, tuple[float, ...]],
    valid_by_policy: dict[str, np.ndarray],
    policies: tuple[DesignPolicy, ...],
    prepared: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
) -> dict[str, tuple[float, ...]]:
    schedules = dict(schedules)
    current, _ = _score(schedules, policies, prepared, members, prior, config)
    for policy in policies:
        fixed = {schedules[policy.name][0], schedules[policy.name][-1]}
        for old in list(schedules[policy.name][1:-1]):
            best_schedule, best_score = schedules[policy.name], current
            for candidate in valid_by_policy[policy.name]:
                if float(candidate) in schedules[policy.name]:
                    continue
                proposed = set(schedules[policy.name])
                proposed.remove(old)
                proposed.add(float(candidate))
                trial = dict(schedules)
                trial[policy.name] = tuple(sorted(proposed))
                score, _ = _score(trial, policies, prepared, members, prior, config)
                if score > best_score + 1e-8:
                    best_schedule, best_score = trial[policy.name], score
            schedules[policy.name] = best_schedule
            current = best_score
    return schedules


def _plot_profiles(
    path: Path,
    policies: tuple[DesignPolicy, ...],
    schedules: dict[str, tuple[float, ...]],
    config: dict,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    for ax, policy in zip(axes, policies):
        design = _future_design(policy, config)
        setpoint_time = np.arange(len(policy.temperature_c) + 1) * 12.0
        setpoint_value = np.r_[policy.temperature_c, policy.temperature_c[-1]]
        ax.step(setpoint_time, setpoint_value, where="post", label="Setpoint", color="tab:orange")
        mask = design.time <= 200.0
        ax.plot(design.time[mask], design.temperature_c[mask], label="T ejecutada estimada", color="tab:red")
        for time in schedules[policy.name]:
            ax.axvline(time, color="tab:blue", alpha=0.25, lw=0.9)
        for time, amount in policy.nutrition_mg_yan_l:
            ax.scatter([time], [24.5], s=30 + amount, marker="v", color="tab:green", zorder=4)
        ax.set_ylabel("°C")
        ax.set_ylim(14.0, 26.0)
        ax.set_title(policy.name)
        ax.grid(alpha=0.2)
    axes[0].legend(ncol=2, loc="lower right")
    axes[-1].set_xlabel("Tiempo desde inoculación (h)")
    fig.suptitle("Candidatos Wave 1 — no autorizados para ejecución física")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_aromas(
    path: Path,
    policies: tuple[DesignPolicy, ...],
    schedules: dict[str, tuple[float, ...]],
    prepared: dict[tuple[int, str], object],
    members: list[int],
) -> None:
    species_names = ["ethyl_acetate", "ethyl_octanoate", "isoamyl_acetate"]
    fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
    for column, policy in enumerate(policies):
        for row, species in enumerate(species_names):
            grid = np.arange(0.0, 192.01, 2.0)
            curves = []
            for member in members:
                item = prepared[(member, policy.name)]
                liquid, _ = simulate_aroma(item.forcings[species], item.aroma_log_values[species])
                curves.append(np.interp(grid, item.forcings[species].time_h, liquid))
            curves = np.asarray(curves)
            ax = axes[row, column]
            ax.fill_between(grid, np.quantile(curves, 0.1, axis=0), np.quantile(curves, 0.9, axis=0), alpha=0.25)
            ax.plot(grid, np.median(curves, axis=0), lw=1.5)
            for time in schedules[policy.name]:
                ax.axvline(time, color="black", alpha=0.18, lw=0.7)
            ax.grid(alpha=0.2)
            if row == 0:
                ax.set_title(policy.name)
            if column == 0:
                ax.set_ylabel(species.replace("_", " ") + "\nµg/L")
    for ax in axes[-1]:
        ax.set_xlabel("Tiempo (h)")
    fig.suptitle("Predicción aromática: mediana y banda 10–90 % del ensamble")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_validation(path: Path, gains: np.ndarray, drying: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].hist(gains, bins=14, color="tab:purple", alpha=0.8)
    axes[0].axvline(np.median(gains), color="black", ls="--", label="Mediana")
    axes[0].set_xlabel("Ganancia log-det")
    axes[0].set_ylabel("Miembros")
    axes[0].legend()
    data = [drying[drying.policy.eq(name)].drying_time_h for name in drying.policy.unique()]
    axes[1].boxplot(data, tick_labels=list(drying.policy.unique()), showfliers=True)
    axes[1].set_ylabel("Tiempo previsto de secado (h)")
    axes[1].tick_params(axis="x", rotation=20)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle("Validación robusta Wave 1 — 64 miembros")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    config = load_json(CONFIG_PATH)
    aroma_config = load_json(AROMA_CONFIG_PATH)
    search_runs = sorted((RESULT_ROOT / "wave1_hybrid_search").glob("*/hybrid_search_gate.json"))
    if not search_runs:
        raise RuntimeError("No Wave-1 hybrid-search result exists")
    search_gate_path = search_runs[-1]
    search_gate = load_json(search_gate_path)
    if search_gate["verdict"] not in {"PASS", "PASS_CONDITIONAL"}:
        raise RuntimeError("Latest Wave-1 search is not usable")
    policies = _load_policies(search_gate_path.parent / "candidate_policy_actions.csv")
    ensemble_run = REPOSITORY_DIR / config["source_contract"]["joint_ensemble_run"]
    ensemble = pd.read_csv(ensemble_run / "joint_parameter_ensemble.csv")
    partitions, provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    representative = representative_members(ensemble, config)[: int(config["search"]["representative_members"])]
    prepared = {}
    drying_rows = []
    for member in range(len(ensemble)):
        for policy in policies:
            item = prepare_design(policy, ensemble.iloc[member], config, partitions)
            if item is None:
                raise RuntimeError(f"Simulation failed: member {member}, {policy.name}")
            prepared[(member, policy.name)] = item
            drying_rows.append(
                {"ensemble_member": member, "policy": policy.name, "drying_time_h": item.drying_time_h}
            )
    drying = pd.DataFrame(drying_rows)
    allowed = allowed_sampling_times(config)
    valid_by_policy = {}
    active_probabilities = {}
    for policy in policies:
        times = drying[drying.policy.eq(policy.name)].drying_time_h.to_numpy(dtype=float)
        probability = {float(time): float(np.mean(times >= float(time))) for time in allowed}
        valid = np.asarray(
            [time for time in allowed if probability[float(time)] >= float(config["completion"]["minimum_probability"])],
            dtype=float,
        )
        valid_by_policy[policy.name] = valid
        active_probabilities[policy.name] = probability
    schedules = {
        policy.name: _greedy_schedule(policy, valid_by_policy[policy.name], prepared, representative, prior, config)
        for policy in policies
    }
    schedules = _coordinate_improve(
        schedules, valid_by_policy, policies, prepared, representative, prior, config
    )
    optimized_score, optimized_gains = _score(
        schedules, policies, prepared, list(range(len(ensemble))), prior, config
    )
    preliminary = {
        policy.name: tuple(float(value) for value in config["sampling"]["preliminary_times_h"])
        for policy in policies
    }
    preliminary_score, preliminary_gains = _score(
        preliminary, policies, prepared, list(range(len(ensemble))), prior, config
    )
    reference_score, reference_gains = _three_anchor_reference(
        schedules[policies[0].name],
        policies[0],
        prepared,
        list(range(len(ensemble))),
        prior,
        config,
    )
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    schedule_rows = []
    for policy in policies:
        for index, time in enumerate(schedules[policy.name], start=1):
            timestamp = start + timedelta(hours=float(time))
            schedule_rows.append(
                {
                    "policy": policy.name,
                    "sample_number": index,
                    "time_h": time,
                    "local_timestamp": timestamp.isoformat(),
                    "weekday": timestamp.strftime("%A"),
                    "active_probability": active_probabilities[policy.name][float(time)],
                }
            )
    schedule = pd.DataFrame(schedule_rows)
    checks = {
        "exactly_ten_samples_per_process": bool(
            schedule.groupby("policy").size().eq(int(config["sampling"]["samples_per_process"])).all()
        ),
        "all_samples_weekdays": bool(schedule["weekday"].isin(["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"]).all()),
        "all_samples_inside_09_17_window": bool(
            schedule["local_timestamp"].map(lambda value: 9 <= datetime.fromisoformat(value).hour <= 17).all()
        ),
        "all_samples_at_least_95pct_active": bool(
            (schedule["active_probability"] >= float(config["completion"]["minimum_probability"])).all()
        ),
        "candidate_beats_three_anchor_terminal_sampling_reference": optimized_score
        > reference_score + 1e-6,
        "full_ensemble_information_finite": bool(np.isfinite(optimized_gains).all()),
        "physical_execution_still_locked": not bool(config["release_policy"]["physical_execution_authorized"]),
    }
    gate = {
        "gate": "phase_D_wave1_sampling_and_plot_review",
        "verdict": "PASS_CONDITIONAL" if all(checks.values()) else "FAIL",
        "checks": checks,
        "information": {
            "preliminary_robust_score": preliminary_score,
            "three_anchor_terminal_sampling_reference_score": reference_score,
            "optimized_robust_score": optimized_score,
            "optimized_median_gain": float(np.median(optimized_gains)),
            "optimized_lower_decile_gain": float(np.quantile(optimized_gains, 0.1)),
        },
        "profiles_for_physical_execution": False,
        "conditions": [
            "Run IPOPT local refinement in the approved solver environment.",
            "Owner must review candidate plots.",
            "Independent Ultra audit is required before physical release."
        ],
    }
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_sampling", config)
    schedule_path = run_dir / "optimized_sampling_schedule.csv"
    drying_path = run_dir / "drying_time_ensemble.csv"
    information_path = run_dir / "full_ensemble_information.csv"
    gate_path = run_dir / "sampling_gate.json"
    profiles_plot = run_dir / "candidate_profiles.png"
    aromas_plot = run_dir / "aroma_predictions.png"
    validation_plot = run_dir / "robust_validation.png"
    config_path = run_dir / "wave1_mbdoe_config.json"
    schedule.to_csv(schedule_path, index=False)
    drying.to_csv(drying_path, index=False)
    pd.DataFrame(
        {
            "ensemble_member": np.arange(len(ensemble)),
            "preliminary_information_gain": preliminary_gains,
            "three_anchor_reference_information_gain": reference_gains,
            "optimized_information_gain": optimized_gains,
        }
    ).to_csv(information_path, index=False)
    write_json(gate_path, gate)
    write_json(config_path, config)
    _plot_profiles(profiles_plot, policies, schedules, config)
    _plot_aromas(aromas_plot, policies, schedules, prepared, list(range(len(ensemble))))
    _plot_validation(validation_plot, optimized_gains, drying)
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_sampling_and_plot_review",
        config=config,
        sources={
            "hybrid_search_manifest": search_gate_path.parent / "run_manifest.json",
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "wave1_config": CONFIG_PATH,
        },
        code_paths=[Path(__file__), ADAPTIVE_DIR / "pilot_mbdoe_adapter.py"],
        random_seeds=[int(config["seed"])],
        status="completed" if gate["verdict"] != "FAIL" else "validation_failed",
        convergence={},
        gate=gate,
        outputs=[schedule_path, drying_path, information_path, gate_path, profiles_plot, aromas_plot, validation_plot, config_path],
    )
    write_json(run_dir / "partition_surrogate_provenance.json", provenance)
    write_json(run_dir / "run_manifest.json", manifest)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    print(schedule.to_string(index=False))
    if gate["verdict"] == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
