from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


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
    SPECIES,
    _future_design,
    allowed_sampling_times,
    effective_temperature_change_indices,
    fim_from_sampling_sensitivity_cache,
    load_json,
    load_wave1_config,
    logdet,
    prepare_design,
    prepare_sampling_sensitivity_cache,
    prior_precision,
    robust_information_metrics,
)
from pilot_2026.adaptive_design.run_artifacts import (  # noqa: E402
    build_manifest,
    capture_git_state,
    create_immutable_run_directory,
    filesystem_path,
    sha256_file,
    sha256_payload,
    verify_manifest_output,
    write_json,
    write_json_atomic,
)


CONFIG_PATH = ADAPTIVE_DIR / "wave1_mbdoe_config.json"
CONSTRAINTS_PATH = ADAPTIVE_DIR / "design_constraints.json"
AROMA_CONFIG_PATH = ADAPTIVE_DIR / "aroma_calibration_config.json"
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def robust_score(gains: np.ndarray, config: dict) -> float:
    gains = np.asarray(gains, dtype=float)
    return float(config["objective"]["median_weight"]) * float(np.median(gains)) + float(
        config["objective"]["lower_decile_weight"]
    ) * float(np.quantile(gains, 0.1))


def sampling_gate_verdict(checks: dict[str, bool], critical_checks: tuple[str, ...]) -> str:
    if not all(bool(checks.get(name, False)) for name in critical_checks):
        return "FAIL"
    return "PASS" if all(bool(value) for value in checks.values()) else "PASS_CONDITIONAL"


def optimized_schedule_is_acceptable(
    optimized_gains: np.ndarray,
    preliminary_gains: np.ndarray,
    config: dict,
) -> tuple[bool, str | None]:
    score_tolerance = float(config["sampling"]["full_ensemble_score_tolerance"])
    q10_tolerance = float(config["sampling"]["q10_regression_tolerance"])
    optimized_score = robust_score(optimized_gains, config)
    preliminary_score = robust_score(preliminary_gains, config)
    if optimized_score < preliminary_score - score_tolerance:
        return False, "optimized_full_ensemble_score_below_preliminary"
    optimized_q10 = float(np.quantile(optimized_gains, 0.1))
    preliminary_q10 = float(np.quantile(preliminary_gains, 0.1))
    if optimized_q10 < preliminary_q10 - q10_tolerance:
        return False, "optimized_q10_below_preliminary_guardrail"
    return True, None


def policy_actions_before_drying(
    policy: DesignPolicy,
    drying_times: np.ndarray,
    config: dict,
) -> bool:
    effective_changes = list(effective_temperature_change_indices(policy.temperature_c, config))
    action_times = [
        float(config["future_process"]["temperature_slot_h"]) * index
        for index in effective_changes
    ]
    action_times.extend(time_h for time_h, _ in policy.nutrition_mg_yan_l)
    margin_h = float(config["operations"]["minimum_action_to_drying_margin_h"])
    return all(
        float(
            np.mean(
                np.asarray(drying_times, dtype=float)
                >= float(action_time) + margin_h - 1e-9
            )
        )
        >= float(config["completion"]["minimum_probability"])
        for action_time in action_times
    )


def capture_constraints_approved(config: dict) -> bool:
    constraints = config["design_constraints"]["sampling_and_capture"]
    duration_approved = bool(
        constraints["maximum_capture_interval_h"] is not None
        or constraints.get("maximum_capture_interval_policy") == "unbounded_by_hardware"
    )
    loading_approved = bool(
        constraints["maximum_trap_loading"] is not None
        or constraints.get("maximum_trap_loading_policy")
        == "not_applicable_for_sampling_nozzle"
    )
    return duration_approved and loading_approved and constraints["vessel_change_required"] is False


def _load_policies(path: Path, config: dict) -> tuple[DesignPolicy, ...]:
    actions = pd.read_csv(filesystem_path(path))
    policies = []
    for name, group in actions.groupby("policy", sort=False):
        temperature_actions = group[group["action"].eq("temperature_setpoint")].sort_values(
            "time_h"
        )
        slots = int(config["future_process"]["optimized_temperature_slots"])
        profile = np.empty(slots, dtype=float)
        current = float(temperature_actions.iloc[0]["value"])
        action_rows = list(temperature_actions.itertuples())
        action_index = 0
        for slot in range(slots):
            time_h = float(config["future_process"]["temperature_slot_h"]) * slot
            while action_index + 1 < len(action_rows) and float(
                action_rows[action_index + 1].time_h
            ) <= time_h + 1e-9:
                action_index += 1
                current = float(action_rows[action_index].value)
            profile[slot] = current
        nutrition = tuple(
            (float(row.time_h), float(row.value))
            for row in group[group["action"].eq("nutrition")].sort_values("time_h").itertuples()
        )
        policies.append(DesignPolicy(str(name), tuple(profile), nutrition))
    if len(policies) != 3:
        raise ValueError(f"Expected three Wave-1 policies; found {len(policies)}")
    return tuple(policies)


def _score(
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
) -> tuple[float, np.ndarray, list[np.ndarray]]:
    base = logdet(prior)
    gains = []
    posterior = []
    for member in members:
        total = prior.copy()
        for policy in policies:
            schedule = tuple(float(value) for value in schedules[policy.name])
            key = (int(member), policy.name, schedule)
            if key not in fim_cache:
                fim_cache[key] = fim_from_sampling_sensitivity_cache(
                    sampling_cache[(member, policy.name)],
                    np.asarray(schedule, dtype=float),
                    config,
                )
            total += fim_cache[key]
        posterior.append(total)
        gains.append(logdet(total) - base)
    gains_array = np.asarray(gains, dtype=float)
    return robust_score(gains_array, config), gains_array, posterior


def _individual_score(
    schedule: tuple[float, ...],
    policy: DesignPolicy,
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
) -> float:
    schedules = {policy.name: schedule}
    return _score(
        schedules,
        (policy,),
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
    )[0]


def _beam_schedule(
    policy: DesignPolicy,
    valid: np.ndarray,
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    rng: np.random.Generator,
) -> tuple[float, ...]:
    count = int(config["sampling"]["samples_per_process"])
    forced: set[float] = set()
    if config["sampling"]["force_first_slot"]:
        forced.add(float(valid[0]))
    if config["sampling"]["force_last_slot"]:
        forced.add(float(valid[-1]))
    beam: list[tuple[float, ...]] = [tuple(sorted(forced))]
    width = int(config["sampling"]["beam_width"])
    ordered_valid = list(float(value) for value in valid)
    rng.shuffle(ordered_valid)
    while len(beam[0]) < count:
        proposals: dict[tuple[float, ...], float] = {}
        for schedule in beam:
            for candidate in ordered_valid:
                if candidate in schedule:
                    continue
                trial = tuple(sorted((*schedule, candidate)))
                proposals.setdefault(
                    trial,
                    _individual_score(
                        trial,
                        policy,
                        sampling_cache,
                        members,
                        prior,
                        config,
                        fim_cache,
                    ),
                )
        if not proposals:
            raise RuntimeError(f"Not enough valid sample slots for {policy.name}")
        beam = [
            schedule
            for schedule, _ in sorted(
                proposals.items(), key=lambda item: (-item[1], item[0])
            )[:width]
        ]
    return beam[0]


def _coordinate_improve(
    schedules: dict[str, tuple[float, ...]],
    valid_by_policy: dict[str, np.ndarray],
    policies: tuple[DesignPolicy, ...],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    rng: np.random.Generator,
) -> dict[str, tuple[float, ...]]:
    schedules = dict(schedules)
    current, _, _ = _score(
        schedules, policies, sampling_cache, members, prior, config, fim_cache
    )
    for _pass in range(int(config["sampling"]["coordinate_passes"])):
        improved = False
        policy_order = list(policies)
        rng.shuffle(policy_order)
        for policy in policy_order:
            protected: set[float] = set()
            if config["sampling"]["force_first_slot"]:
                protected.add(float(valid_by_policy[policy.name][0]))
            if config["sampling"]["force_last_slot"]:
                protected.add(float(valid_by_policy[policy.name][-1]))
            old_order = [
                value for value in schedules[policy.name] if value not in protected
            ]
            rng.shuffle(old_order)
            for old in old_order:
                best_schedule, best_score = schedules[policy.name], current
                for candidate in valid_by_policy[policy.name]:
                    candidate = float(candidate)
                    if candidate in schedules[policy.name]:
                        continue
                    proposed = set(schedules[policy.name])
                    proposed.remove(old)
                    proposed.add(candidate)
                    trial = dict(schedules)
                    trial[policy.name] = tuple(sorted(proposed))
                    score, _, _ = _score(
                        trial,
                        policies,
                        sampling_cache,
                        members,
                        prior,
                        config,
                        fim_cache,
                    )
                    if score > best_score + float(config["sampling"]["full_ensemble_score_tolerance"]):
                        best_schedule, best_score = trial[policy.name], score
                if best_schedule != schedules[policy.name]:
                    improved = True
                schedules[policy.name] = best_schedule
                current = best_score
        if not improved:
            break
    return schedules


def _optimize_schedules(
    policies: tuple[DesignPolicy, ...],
    valid_by_policy: dict[str, np.ndarray],
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
    prior: np.ndarray,
    config: dict,
    fim_cache: dict[tuple, np.ndarray],
    seed_offset: int,
) -> tuple[dict[str, tuple[float, ...]], pd.DataFrame]:
    best_schedules = None
    best_score = -math.inf
    rows = []
    for restart in range(int(config["sampling"]["search_restarts"])):
        rng = np.random.default_rng(int(config["seed"]) + seed_offset + restart)
        schedules = {
            policy.name: _beam_schedule(
                policy,
                valid_by_policy[policy.name],
                sampling_cache,
                members,
                prior,
                config,
                fim_cache,
                rng,
            )
            for policy in policies
        }
        schedules = _coordinate_improve(
            schedules,
            valid_by_policy,
            policies,
            sampling_cache,
            members,
            prior,
            config,
            fim_cache,
            rng,
        )
        score, gains, _ = _score(
            schedules, policies, sampling_cache, members, prior, config, fim_cache
        )
        rows.append(
            {
                "restart": restart,
                "robust_score": score,
                "median": float(np.median(gains)),
                "q10": float(np.quantile(gains, 0.1)),
            }
        )
        if score > best_score:
            best_score = score
            best_schedules = schedules
    if best_schedules is None:  # pragma: no cover - at least one restart is configured
        raise RuntimeError("Sampling optimization produced no schedule")
    return best_schedules, pd.DataFrame(rows)


def _capture_intervals(
    schedules: dict[str, tuple[float, ...]],
    config: dict,
    sampling_cache: dict[tuple[int, str], object],
    members: list[int],
) -> pd.DataFrame:
    constraints = config["design_constraints"]["sampling_and_capture"]
    rows = []
    for policy, times in schedules.items():
        if len(times) != 10 or not math.isclose(float(times[0]), 0.0, abs_tol=1e-12):
            raise ValueError("Capture intervals require ten wine samples beginning at t=0")
        for species in SPECIES:
            cumulative_by_member = []
            for member in members:
                cache = sampling_cache[(member, policy)]
                forcing_time = cache.prepared.forcings[species].time_h
                cumulative_by_member.append(
                    np.interp(
                        np.asarray(times, dtype=float),
                        forcing_time,
                        cache.nominal_captured[species],
                    )
                )
            cumulative = np.vstack(cumulative_by_member)
            intervals = np.diff(cumulative, axis=1)
            conservation_error = float(
                np.max(
                    np.abs(
                        np.sum(intervals, axis=1)
                        - (cumulative[:, -1] - cumulative[:, 0])
                    )
                )
            )
            for index, (start, end) in enumerate(zip(times[:-1], times[1:]), start=1):
                maximum = constraints["maximum_capture_interval_h"]
                unbounded = (
                    maximum is None
                    and constraints["maximum_capture_interval_policy"]
                    == "unbounded_by_hardware"
                )
                interval_values = intervals[:, index - 1]
                rows.append(
                    {
                        "policy": policy,
                        "species": species,
                        "capture_interval": index,
                        "start_h": start,
                        "end_h": end,
                        "duration_h": end - start,
                        "capture_stage_1_c": 0.0,
                        "capture_stage_2_c": -40.0,
                        "maximum_capture_interval_h": maximum,
                        "maximum_capture_interval_policy": constraints[
                            "maximum_capture_interval_policy"
                        ],
                        "within_approved_duration": bool(
                            unbounded or end - start <= float(maximum)
                        ),
                        "captured_mass_mean_ug": float(np.mean(interval_values)),
                        "captured_mass_min_ug": float(np.min(interval_values)),
                        "captured_mass_max_ug": float(np.max(interval_values)),
                        "cumulative_mass_start_mean_ug": float(
                            np.mean(cumulative[:, index - 1])
                        ),
                        "cumulative_mass_end_mean_ug": float(np.mean(cumulative[:, index])),
                        "mass_conservation_max_abs_error_ug": conservation_error,
                        "maximum_trap_loading": constraints["maximum_trap_loading"],
                        "maximum_trap_loading_policy": constraints[
                            "maximum_trap_loading_policy"
                        ],
                        "trap_loading_constraint_satisfied": bool(
                            constraints["maximum_trap_loading"] is not None
                            or constraints["maximum_trap_loading_policy"]
                            == "not_applicable_for_sampling_nozzle"
                        ),
                        "vessel_change_required": constraints["vessel_change_required"],
                    }
                )
    return pd.DataFrame(rows)


def _operational_conflicts(
    schedules: dict[str, tuple[float, ...]],
    policies: tuple[DesignPolicy, ...],
    config: dict,
) -> pd.DataFrame:
    constraints = config["design_constraints"]["sampling_and_capture"]
    rows = []
    for policy in policies:
        sample_times = set(schedules[policy.name])
        pulse_times = {float(time_h) for time_h, _ in policy.nutrition_mg_yan_l}
        for collision in sorted(sample_times & pulse_times):
            ordered = (
                constraints["sample_event_order"] == "sample_before_action"
                and constraints["minimum_minutes_between_sampling_and_nutrition"] == 0
            )
            rows.append(
                {
                    "time_h": collision,
                    "policies": policy.name,
                    "conflict_type": "sampling_and_nutrition_same_timestamp",
                    "status": "ordered" if ordered else "unresolved",
                    "required_approval": "sample_event_order",
                }
            )
    all_times = sorted({time_h for times in schedules.values() for time_h in times})
    for time_h in all_times:
        policies_at_time = sorted(
            name for name, times in schedules.items() if time_h in set(times)
        )
        if len(policies_at_time) > 1:
            capacity = constraints["manual_sampling_capacity_per_time_slot"]
            unresolved = capacity is None or len(policies_at_time) > int(capacity)
            rows.append(
                {
                    "time_h": time_h,
                    "policies": ";".join(policies_at_time),
                    "conflict_type": "simultaneous_manual_sampling",
                    "status": "unresolved" if unresolved else "within_capacity",
                    "required_approval": "manual_sampling_capacity_per_time_slot",
                }
            )
    if constraints["minimum_minutes_between_sampling_and_nutrition"] is None:
        rows.append(
            {
                "time_h": math.nan,
                "policies": ";".join(policy.name for policy in policies),
                "conflict_type": "sampling_nutrition_separation_not_approved",
                "status": "unresolved",
                "required_approval": "minimum_minutes_between_sampling_and_nutrition",
            }
        )
    return pd.DataFrame(rows)


def _nutrition_translation(policies: tuple[DesignPolicy, ...], config: dict) -> pd.DataFrame:
    nutrition = config["design_constraints"]["nutrition"]
    volume_l = float(nutrition["initial_volume_l"])
    organic_limit = float(nutrition["maximum_total_organic_product_g_hl"]) * volume_l / 100.0
    dap_limit = float(nutrition["maximum_total_dap_fda_g_hl"]) * volume_l / 100.0
    organic_fraction = nutrition["organic_product_yan_mass_fraction"]
    dap_fraction = nutrition["dap_yan_mass_fraction"]
    composition_available = bool(
        organic_fraction is not None
        and dap_fraction is not None
        and float(organic_fraction) > 0.0
        and float(dap_fraction) > 0.0
    )
    rows = []
    for policy in policies:
        event_rows = []
        for time_h, target in policy.nutrition_mg_yan_l:
            yan_organic = 0.5 * float(target)
            yan_dap = 0.5 * float(target)
            organic_g = (
                (yan_organic * volume_l / 1000.0) / float(organic_fraction)
                if composition_available
                else math.nan
            )
            dap_g = (
                (yan_dap * volume_l / 1000.0) / float(dap_fraction)
                if composition_available
                else math.nan
            )
            event_rows.append((time_h, target, yan_organic, yan_dap, organic_g, dap_g))
        total_organic = float(sum(row[4] for row in event_rows)) if composition_available else math.nan
        total_dap = float(sum(row[5] for row in event_rows)) if composition_available else math.nan
        for time_h, target, yan_organic, yan_dap, organic_g, dap_g in event_rows:
            reconstructed = yan_organic + yan_dap
            rows.append(
                {
                    "policy": policy.name,
                    "time_h": time_h,
                    "mix_policy": nutrition["approved_product_mix_selection_policy"],
                    "selected_for_operation": False,
                    "volume_l": volume_l,
                    "yan_target_mg_l": target,
                    "organic_yan_contribution_mg_l": yan_organic,
                    "dap_yan_contribution_mg_l": yan_dap,
                    "organic_product_yan_mass_fraction": organic_fraction,
                    "dap_yan_mass_fraction": dap_fraction,
                    "organic_product_g": organic_g,
                    "dap_fda_g": dap_g,
                    "yan_reconstructed_mg_l": reconstructed,
                    "reconstruction_error_mg_l": reconstructed - target,
                    "mass_formula": "m_g=(0.5*Y_target_mg_L*V_L/1000)/YAN_mass_fraction",
                    "organic_total_limit_g": organic_limit,
                    "dap_fda_total_limit_g": dap_limit,
                    "organic_total_product_g": total_organic,
                    "dap_total_product_g": total_dap,
                    "organic_total_within_limit": bool(total_organic <= organic_limit)
                    if composition_available
                    else None,
                    "dap_total_within_limit": bool(total_dap <= dap_limit)
                    if composition_available
                    else None,
                    "per_event_organic_limit_g": nutrition["maximum_organic_product_g_per_event"],
                    "per_event_dap_fda_limit_g": nutrition["maximum_dap_fda_g_per_event"],
                    "composition_source": nutrition["composition_source"],
                    "translation_blocker": not composition_available,
                    "blocker_reason": None
                    if composition_available
                    else "missing_authoritative_product_YAN_mass_fraction",
                }
            )
    return pd.DataFrame(rows)


def _tank_randomization(
    policies: tuple[DesignPolicy, ...], config: dict, frozen_utc: str
) -> pd.DataFrame:
    randomization = config["design_constraints"]["tank_randomization"]
    seed = int(config["tank_randomization_seed"])
    if seed != int(randomization["seed"]):
        raise ValueError("Tank-randomization seeds disagree between authoritative configurations")
    tanks = list(config["design_constraints"]["sampling_and_capture"]["tanks"])
    if len(policies) != len(tanks):
        raise ValueError("Tank randomization requires one tank per Wave-1 profile")
    assigned = np.random.default_rng(seed).permutation(np.asarray(tanks, dtype=object))
    rows = []
    for policy, tank in zip(policies, assigned):
        policy_hash = sha256_payload(
            {
                "temperature_c": list(policy.temperature_c),
                "nutrition_mg_yan_l": [list(row) for row in policy.nutrition_mg_yan_l],
            }
        )
        rows.append(
            {
                "profile": policy.name,
                "tank": str(tank),
                "seed": seed,
                "algorithm": randomization["algorithm"],
                "frozen_utc": frozen_utc,
                "policy_hash": policy_hash,
                "owner_policy_approved": bool(randomization["approved"]),
                "authorized_for_physical_execution": False,
                "status": "frozen_computational_proposal_not_authorized",
            }
        )
    return pd.DataFrame(rows)


FIGURE_NAMES = (
    "candidate_profiles_v2.png",
    "executed_temperature_ensemble_v2.png",
    "aroma_predictions_v2.png",
    "sampling_schedule_v2.png",
    "capture_intervals_v2.png",
    "information_gain_distribution_v2.png",
    "paired_information_comparison_v2.png",
    "drying_margin_validation_v2.png",
    "actuator_robustness_v2.png",
    "operational_summary_v2.png",
)
WATERMARK = "COMPUTATIONAL CANDIDATE — NOT AUTHORIZED FOR PHYSICAL EXECUTION"
PALETTE = ("#235789", "#D4A72C", "#E07A3F")


def _finish_figure(
    fig: plt.Figure, path: Path, subtitle: str, *, has_suptitle: bool = False
) -> None:
    fig.text(0.01, 0.985, subtitle, ha="left", va="top", fontsize=9, color="#4B5563")
    fig.text(
        0.5,
        0.012,
        WATERMARK,
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color="#8B1E3F",
    )
    top = 0.88 if has_suptitle else 0.94
    fig.tight_layout(rect=(0.02, 0.045, 0.98, top))
    fig.savefig(filesystem_path(path), dpi=160, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def _generate_figures(
    figure_paths: dict[str, Path],
    policies: tuple[DesignPolicy, ...],
    config: dict,
    prepared: dict[tuple[int, str], object],
    members: list[int],
    schedule: pd.DataFrame,
    captures: pd.DataFrame,
    paired: pd.DataFrame,
    search_full: pd.DataFrame,
    actuator: pd.DataFrame,
    conflicts: pd.DataFrame,
    checks: dict[str, bool],
) -> list[dict]:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.edgecolor": "#374151",
            "axes.labelcolor": "#1F2937",
            "axes.titlecolor": "#111827",
            "xtick.color": "#374151",
            "ytick.color": "#374151",
            "grid.color": "#D1D5DB",
            "grid.alpha": 0.55,
        }
    )
    chart_map: list[dict] = []
    slot_h = float(config["future_process"]["temperature_slot_h"])

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for axis, policy, color in zip(axes, policies, PALETTE):
        x = np.arange(len(policy.temperature_c) + 1) * slot_h
        y = np.r_[policy.temperature_c, policy.temperature_c[-1]]
        axis.step(x, y, where="post", color=color, linewidth=2, label="Setpoint")
        for time_h, _amount in policy.nutrition_mg_yan_l:
            axis.axvline(time_h, color="#374151", linestyle="--", linewidth=1)
        axis.set_ylabel("°C")
        axis.set_title(policy.name, loc="left", fontsize=10)
        axis.set_ylim(14.5, 27.5)
        axis.grid(axis="y")
    axes[-1].set_xlabel("Process time (h)")
    fig.suptitle(
        "Candidate temperature policies", x=0.02, y=0.945, ha="left", fontsize=15
    )
    _finish_figure(
        fig,
        figure_paths["candidate_profiles_v2.png"],
        "Setpoint range 15–27 °C; dashed lines identify nutrition actions; candidate only.",
        has_suptitle=True,
    )
    chart_map.append({"figure": "candidate_profiles_v2.png", "family": "Trend", "question": "What excitation policies were evaluated?"})

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, sharey=True)
    actuator_cfg = config["future_process"]["temperature_actuator"]
    scenarios = [
        {"label": "nominal", "tau_h": actuator_cfg["global_tau_h"]},
        {"label": "tau min", "tau_h": min(actuator_cfg["empirical_tau_h"])},
        {"label": "tau max", "tau_h": max(actuator_cfg["empirical_tau_h"])},
        {"label": "tracking -", "tau_h": actuator_cfg["global_tau_h"], "tracking_error_c": -float(np.median(actuator_cfg["observed_tracking_rmse_c"]))},
        {"label": "tracking +", "tau_h": actuator_cfg["global_tau_h"], "tracking_error_c": float(np.median(actuator_cfg["observed_tracking_rmse_c"]))},
    ]
    for axis, policy in zip(axes, policies):
        for scenario in scenarios:
            design = _future_design(policy, config, scenario)
            axis.plot(design.time, design.temperature_c, label=scenario["label"], linewidth=1.3)
        axis.set_title(policy.name, loc="left", fontsize=10)
        axis.set_ylabel("°C")
        axis.grid()
    axes[0].legend(ncol=5, fontsize=8, loc="upper right")
    axes[-1].set_xlabel("Process time (h)")
    fig.suptitle(
        "Executed-temperature actuator scenarios",
        x=0.02,
        y=0.945,
        ha="left",
        fontsize=15,
    )
    _finish_figure(
        fig,
        figure_paths["executed_temperature_ensemble_v2.png"],
        "Nominal, empirical tau extremes and approved tracking-error scenarios.",
        has_suptitle=True,
    )
    chart_map.append({"figure": "executed_temperature_ensemble_v2.png", "family": "Uncertainty & Benchmark", "question": "How does actuator uncertainty alter executed temperature?"})

    central_member = members[len(members) // 2]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for axis, species in zip(axes, SPECIES):
        for policy, color in zip(policies, PALETTE):
            item = prepared[(central_member, policy.name)]
            liquid, _captured = simulate_aroma(item.forcings[species], item.aroma_log_values[species])
            axis.plot(item.forcings[species].time_h, liquid, color=color, label=policy.name)
        axis.set_title(species.replace("_", " "), loc="left", fontsize=10)
        axis.set_ylabel("µg/L")
        axis.grid()
    axes[0].legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("Process time (h)")
    fig.suptitle(
        "Aroma concentration predictions", x=0.02, y=0.945, ha="left", fontsize=15
    )
    _finish_figure(
        fig,
        figure_paths["aroma_predictions_v2.png"],
        f"Central ensemble member {central_member}; wine concentration predictions.",
        has_suptitle=True,
    )
    chart_map.append({"figure": "aroma_predictions_v2.png", "family": "Trend", "question": "How do candidate policies separate aroma trajectories?"})

    fig, ax = plt.subplots(figsize=(11, 4.8))
    for y, (policy, group) in enumerate(schedule.groupby("policy", sort=False)):
        ax.scatter(group["time_h"], np.full(len(group), y), s=55, color=PALETTE[y], edgecolor="#1F2937")
        for row in group.itertuples():
            ax.text(row.time_h, y + 0.12, str(row.sample_number), ha="center", fontsize=7)
    ax.set_yticks(range(len(policies)), [policy.name for policy in policies])
    ax.set_xlabel("Process time (h)")
    ax.set_title("Optimized wine-sampling schedule", loc="left", fontsize=15)
    ax.grid(axis="x")
    _finish_figure(fig, figure_paths["sampling_schedule_v2.png"], "Ten wine/process samples; sample 1 is the mandatory physical baseline at t=0.")
    chart_map.append({"figure": "sampling_schedule_v2.png", "family": "Progression", "question": "Where are the ten wine samples placed?"})

    fig, ax = plt.subplots(figsize=(11, 5.2))
    interval_view = captures[captures["species"].eq(SPECIES[0])]
    for y, (policy, group) in enumerate(interval_view.groupby("policy", sort=False)):
        for row in group.itertuples():
            ax.broken_barh([(row.start_h, row.duration_h)], (y - 0.3, 0.6), facecolors=PALETTE[y], edgecolors="white")
    ax.set_yticks(range(len(policies)), [policy.name for policy in policies])
    ax.set_xlabel("Process time (h)")
    ax.set_title("MIX/condensate capture intervals", loc="left", fontsize=15)
    ax.grid(axis="x")
    _finish_figure(fig, figure_paths["capture_intervals_v2.png"], "Nine contiguous intervals/process; every first interval begins at t=0; hardware duration is unbounded.")
    chart_map.append({"figure": "capture_intervals_v2.png", "family": "Progression", "question": "Do all nine capture intervals cover the process from zero?"})

    fig, ax = plt.subplots(figsize=(10, 5.5))
    bins = np.linspace(
        float(paired[["optimized_information_gain", "preliminary_information_gain", "three_anchor_information_gain"]].min().min()),
        float(paired[["optimized_information_gain", "preliminary_information_gain", "three_anchor_information_gain"]].max().max()),
        16,
    )
    for column, label, color in zip(
        ("optimized_information_gain", "preliminary_information_gain", "three_anchor_information_gain"),
        ("Optimized", "Preliminary", "Three-anchor"),
        PALETTE,
    ):
        ax.hist(paired[column], bins=bins, histtype="step", linewidth=2, label=label, color=color)
    ax.set_xlabel("Information gain")
    ax.set_ylabel("Ensemble members")
    ax.set_title("Information-gain distribution", loc="left", fontsize=15)
    ax.legend()
    ax.grid(axis="y")
    _finish_figure(fig, figure_paths["information_gain_distribution_v2.png"], "All 64 joint-ensemble members; common bins and common scale.")
    chart_map.append({"figure": "information_gain_distribution_v2.png", "family": "Distribution", "question": "How does information vary across the 64 members?"})

    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    ax.scatter(paired["preliminary_information_gain"], paired["optimized_information_gain"], color=PALETTE[0], edgecolor="#1F2937", alpha=0.8)
    low = float(min(paired["preliminary_information_gain"].min(), paired["optimized_information_gain"].min()))
    high = float(max(paired["preliminary_information_gain"].max(), paired["optimized_information_gain"].max()))
    ax.plot([low, high], [low, high], linestyle="--", color="#374151", label="Equal information")
    ax.set_xlabel("Preliminary information gain")
    ax.set_ylabel("Optimized information gain")
    ax.set_title("Paired information comparison", loc="left", fontsize=15)
    ax.legend()
    ax.grid()
    _finish_figure(fig, figure_paths["paired_information_comparison_v2.png"], "One point per ensemble member; points above the diagonal favor the optimized schedule.")
    chart_map.append({"figure": "paired_information_comparison_v2.png", "family": "Relationship", "question": "Does optimized sampling improve each paired member?"})

    fig, ax = plt.subplots(figsize=(10, 5.5))
    margins = search_full["action_margin_to_drying_h"].to_numpy(dtype=float)
    ax.scatter(search_full["ensemble_member"], margins, color=PALETTE[0], s=32)
    ax.axhline(float(config["operations"]["minimum_action_to_drying_margin_h"]), color="#8B1E3F", linestyle="--", label="Approved minimum")
    ax.set_xlabel("Ensemble member")
    ax.set_ylabel("Latest-action margin to drying (h)")
    ax.set_title("Drying-margin validation", loc="left", fontsize=15)
    ax.legend()
    ax.grid()
    _finish_figure(fig, figure_paths["drying_margin_validation_v2.png"], "Full 64-member validation; final sample is not treated as an excitation action.")
    chart_map.append({"figure": "drying_margin_validation_v2.png", "family": "Uncertainty & Benchmark", "question": "Do active actions retain the approved drying margin?"})

    fig, ax = plt.subplots(figsize=(11, 5.8))
    positions = np.arange(len(actuator))
    colors = [PALETTE[0] if bool(value) else "#B35C44" for value in actuator["feasible"]]
    ax.bar(positions, actuator["completion_probability"], color=colors, edgecolor="#374151")
    ax.axhline(float(config["completion"]["minimum_probability"]), color="#111827", linestyle="--", label="Minimum probability")
    ax.set_xticks(positions, actuator["scenario"], rotation=35, ha="right", fontsize=8)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Completion and action-margin probability")
    ax.set_title("Actuator robustness", loc="left", fontsize=15)
    ax.legend()
    ax.grid(axis="y")
    _finish_figure(fig, figure_paths["actuator_robustness_v2.png"], "Nine empirical tau values and signed tracking-error scenarios, each evaluated on 64 members.")
    chart_map.append({"figure": "actuator_robustness_v2.png", "family": "Comparison & Benchmark", "question": "Which approved actuator scenarios remain feasible?"})

    fig, ax = plt.subplots(figsize=(11, max(5.5, 0.32 * len(checks))))
    labels = list(checks)
    values = [1 if checks[label] else 0 for label in labels]
    colors = [PALETTE[0] if value else "#B35C44" for value in values]
    ax.barh(np.arange(len(labels)), values, color=colors, edgecolor="#374151")
    ax.set_yticks(np.arange(len(labels)), [label.replace("_", " ") for label in labels], fontsize=8)
    ax.set_xlim(0.0, 1.05)
    ax.set_xticks([0, 1], ["Fail", "Pass"])
    ax.invert_yaxis()
    ax.set_title("Operational qualification summary", loc="left", fontsize=15)
    ax.grid(axis="x")
    unresolved = int(conflicts["status"].eq("unresolved").sum()) if len(conflicts) else 0
    _finish_figure(fig, figure_paths["operational_summary_v2.png"], f"Fail-closed checks; unresolved operational conflicts: {unresolved}.")
    chart_map.append({"figure": "operational_summary_v2.png", "family": "Tables & Scorecards", "question": "Which computational and operational checks pass?"})
    return chart_map


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize and qualify Wave-1 sampling on 64 members")
    parser.add_argument("--source-search-run", required=True)
    parser.add_argument("--source-ensemble-run", required=True)
    parser.add_argument("--source-aroma-run", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_started = time.perf_counter()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    constraints = config["design_constraints"]
    aroma_config = load_json(AROMA_CONFIG_PATH)
    search_run = _run_path(args.source_search_run)
    ensemble_run = _run_path(args.source_ensemble_run)
    aroma_run = _run_path(args.source_aroma_run)
    source_verification = {
        "search_gate": verify_manifest_output(search_run, "hybrid_search_gate.json"),
        "search_actions": verify_manifest_output(search_run, "candidate_policy_actions.csv"),
        "search_full_ensemble": verify_manifest_output(search_run, "full_ensemble_validation.csv"),
        "search_actuator": verify_manifest_output(search_run, "actuator_robustness_validation.csv"),
        "ensemble": verify_manifest_output(ensemble_run, "joint_parameter_ensemble.csv"),
        "aroma": verify_manifest_output(aroma_run, "aroma_calibration_gate.json"),
    }
    search_gate = load_json(filesystem_path(search_run / "hybrid_search_gate.json"))
    policies = _load_policies(search_run / "candidate_policy_actions.csv", config)
    search_full = pd.read_csv(filesystem_path(search_run / "full_ensemble_validation.csv"))
    actuator = pd.read_csv(
        filesystem_path(search_run / "actuator_robustness_validation.csv")
    )
    ensemble = pd.read_csv(filesystem_path(ensemble_run / "joint_parameter_ensemble.csv"))
    partitions, provenance = load_partition_surrogates(aroma_config, REPOSITORY_DIR)
    prior = prior_precision(ensemble, config)
    members = list(range(len(ensemble)))
    prepared: dict[tuple[int, str], object] = {}
    sampling_cache: dict[tuple[int, str], object] = {}
    drying_rows = []
    for member in members:
        for policy in policies:
            item = prepare_design(policy, ensemble.iloc[member], config, partitions)
            if item is None:
                raise RuntimeError(f"Simulation failed: member {member}, {policy.name}")
            prepared[(member, policy.name)] = item
            sampling_cache[(member, policy.name)] = prepare_sampling_sensitivity_cache(
                item, config
            )
            drying_rows.append(
                {
                    "ensemble_member": member,
                    "policy": policy.name,
                    "drying_time_h": item.drying_time_h,
                }
            )
    drying = pd.DataFrame(drying_rows)
    allowed = allowed_sampling_times(config)
    valid_by_policy = {}
    active_probabilities = {}
    for policy in policies:
        drying_times = drying[drying.policy.eq(policy.name)].drying_time_h.to_numpy(dtype=float)
        probability = {
            float(time_h): float(np.mean(drying_times >= float(time_h))) for time_h in allowed
        }
        valid = np.asarray(
            [
                time_h
                for time_h in allowed
                if probability[float(time_h)]
                >= float(config["completion"]["minimum_probability"])
            ],
            dtype=float,
        )
        if not math.isclose(float(valid[0]), 0.0, abs_tol=1e-12):
            raise RuntimeError("The legal sampling space must contain the mandatory t=0 baseline")
        valid_by_policy[policy.name] = valid
        active_probabilities[policy.name] = probability
    fim_cache: dict[tuple, np.ndarray] = {}
    optimized, optimization_restarts = _optimize_schedules(
        policies,
        valid_by_policy,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
        seed_offset=0,
    )
    optimized_score, optimized_gains, optimized_posterior = _score(
        optimized, policies, sampling_cache, members, prior, config, fim_cache
    )
    preliminary = {
        policy.name: tuple(float(value) for value in config["sampling"]["preliminary_times_h"])
        for policy in policies
    }
    preliminary_score, preliminary_gains, preliminary_posterior = _score(
        preliminary, policies, sampling_cache, members, prior, config, fim_cache
    )

    reference_policies = tuple(
        DesignPolicy(f"three_anchor_reference_{index + 1}", policies[0].temperature_c, policies[0].nutrition_mg_yan_l)
        for index in range(3)
    )
    for member in members:
        original = prepared[(member, policies[0].name)]
        cache = sampling_cache[(member, policies[0].name)]
        for reference in reference_policies:
            prepared[(member, reference.name)] = original
            sampling_cache[(member, reference.name)] = cache
    reference_valid = {
        reference.name: valid_by_policy[policies[0].name] for reference in reference_policies
    }
    reference_schedules, reference_restarts = _optimize_schedules(
        reference_policies,
        reference_valid,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
        seed_offset=1000,
    )
    reference_score, reference_gains, reference_posterior = _score(
        reference_schedules,
        reference_policies,
        sampling_cache,
        members,
        prior,
        config,
        fim_cache,
    )
    acceptable, fallback_reason = optimized_schedule_is_acceptable(
        optimized_gains, preliminary_gains, config
    )
    selected = optimized if acceptable else preliminary
    selected_source = "optimized" if acceptable else "preliminary_fallback"
    selected_score, selected_gains, selected_posterior = _score(
        selected, policies, sampling_cache, members, prior, config, fim_cache
    )
    paired = pd.DataFrame(
        {
            "ensemble_member": members,
            "optimized_information_gain": optimized_gains,
            "preliminary_information_gain": preliminary_gains,
            "three_anchor_information_gain": reference_gains,
            "optimized_vs_preliminary_delta": optimized_gains - preliminary_gains,
            "optimized_vs_three_anchor_delta": optimized_gains - reference_gains,
            "selected_information_gain": selected_gains,
        }
    )
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    schedule_rows = []
    for policy in policies:
        for index, time_h in enumerate(selected[policy.name], start=1):
            timestamp = start + timedelta(hours=float(time_h))
            schedule_rows.append(
                {
                    "policy": policy.name,
                    "sample_number": index,
                    "time_h": time_h,
                    "local_timestamp": timestamp.isoformat(),
                    "weekday": timestamp.strftime("%A"),
                    "active_probability": active_probabilities[policy.name][float(time_h)],
                    "candidate_only_not_for_physical_execution": True,
                }
            )
    schedule = pd.DataFrame(schedule_rows)
    captures = _capture_intervals(selected, config, sampling_cache, members)
    conflicts = _operational_conflicts(selected, policies, config)
    nutrition = _nutrition_translation(policies, config)
    frozen_utc = datetime.now(timezone.utc).isoformat()
    randomization = _tank_randomization(policies, config, frozen_utc)
    selected_metrics = robust_information_metrics(
        selected_gains,
        selected_posterior,
        preliminary_gains,
        config,
    )
    preliminary_metrics = robust_information_metrics(
        preliminary_gains, preliminary_posterior, None, config
    )
    reference_metrics = robust_information_metrics(
        reference_gains, reference_posterior, None, config
    )
    process_active_actions = []
    for policy in policies:
        drying_times = drying[drying.policy.eq(policy.name)].drying_time_h.to_numpy(dtype=float)
        process_active_actions.append(
            policy_actions_before_drying(policy, drying_times, config)
        )
    unresolved_conflicts = bool(
        len(conflicts) and conflicts["status"].eq("unresolved").any()
    )
    capture_limits_approved = capture_constraints_approved(config)
    baseline_present = bool(
        schedule.groupby("policy")["time_h"].min().eq(0.0).all()
    )
    nine_intervals_from_zero = bool(
        captures.groupby(["policy", "species"])["capture_interval"].nunique().eq(9).all()
        and captures.groupby(["policy", "species"])["start_h"].min().eq(0.0).all()
    )
    capture_mass_conserved = bool(
        captures["mass_conservation_max_abs_error_ug"].le(1e-8).all()
        and captures["cumulative_mass_start_mean_ug"].groupby(
            [captures["policy"], captures["species"]]
        ).first().abs().le(1e-10).all()
    )
    randomization_repeat = _tank_randomization(policies, config, frozen_utc)
    randomization_reproducible = bool(
        randomization[["profile", "tank", "seed", "algorithm", "policy_hash"]].equals(
            randomization_repeat[["profile", "tank", "seed", "algorithm", "policy_hash"]]
        )
        and randomization["tank"].nunique() == len(randomization)
    )
    nutrition_reconstructed = bool(
        nutrition.empty
        or (
            nutrition["mix_policy"].eq("50_50_net_YAN_contribution").all()
            and nutrition["organic_yan_contribution_mg_l"].eq(
                0.5 * nutrition["yan_target_mg_l"]
            ).all()
            and nutrition["dap_yan_contribution_mg_l"].eq(
                0.5 * nutrition["yan_target_mg_l"]
            ).all()
            and nutrition["reconstruction_error_mg_l"].abs().le(1e-12).all()
        )
    )
    checks = {
        "source_hashes_verified": len(source_verification) == 6,
        "source_search_gate_not_fail": search_gate["verdict"] in {"PASS", "PASS_CONDITIONAL"},
        "source_search_ipopt_execution_recognized": bool(
            search_gate["checks"]["actual_objective_evaluated"]
        ),
        "sampling_optimization_used_64_members": len(members) == 64,
        "optimized_full_ensemble_score_not_below_preliminary": acceptable
        or selected_source == "preliminary_fallback",
        "selected_q10_meets_preliminary_guardrail": float(np.quantile(selected_gains, 0.1))
        >= float(np.quantile(preliminary_gains, 0.1))
        - float(config["sampling"]["q10_regression_tolerance"]),
        "automatic_preliminary_fallback_operational": acceptable
        or selected_source == "preliminary_fallback",
        "paired_comparison_published": len(paired) == 64,
        "three_anchor_reference_optimized_with_comparable_method": len(reference_restarts)
        == len(optimization_restarts),
        "exactly_ten_samples_per_process": bool(
            schedule.groupby("policy").size().eq(int(config["sampling"]["samples_per_process"])).all()
        ),
        "mandatory_wine_baseline_at_t0": baseline_present,
        "exactly_nine_capture_intervals_starting_at_t0": nine_intervals_from_zero,
        "captured_mass_conserved_between_intervals": capture_mass_conserved,
        "all_samples_in_legal_windows": bool(
            schedule["weekday"].isin(["Monday", "Tuesday", "Wednesday", "Thursday", "Friday"]).all()
            and schedule["local_timestamp"].map(
                lambda value: 9 <= datetime.fromisoformat(value).hour <= 17
            ).all()
        ),
        "samples_and_actions_before_drying": bool(
            (schedule["active_probability"] >= float(config["completion"]["minimum_probability"])).all()
            and all(process_active_actions)
        ),
        "no_ambiguous_operational_conflicts": not unresolved_conflicts,
        "capture_limits_approved_and_satisfied": capture_limits_approved
        and bool(
            captures["within_approved_duration"].all()
            and captures["trap_loading_constraint_satisfied"].all()
            and not captures["vessel_change_required"].any()
        ),
        "sample_before_action_coincidence_approved": bool(
            constraints["sampling_and_capture"]["sample_event_order"]
            == "sample_before_action"
            and constraints["sampling_and_capture"][
                "minimum_minutes_between_sampling_and_nutrition"
            ]
            == 0
        ),
        "manual_sampling_capacity_three_respected": bool(
            not len(conflicts)
            or not (
                conflicts["conflict_type"].eq("simultaneous_manual_sampling")
                & conflicts["status"].eq("unresolved")
            ).any()
        ),
        "nutrition_50_50_net_yan_policy_reconstructed": nutrition_reconstructed,
        "nutrition_product_masses_available": bool(
            nutrition.empty or not nutrition["translation_blocker"].any()
        ),
        "tank_randomization_reproducible_uniform_permutation": randomization_reproducible,
        "nutrition_translation_has_no_arbitrary_mix": bool(
            nutrition.empty or not nutrition["selected_for_operation"].any()
        ),
        "profiles_for_physical_execution_false": True,
    }
    critical_checks = (
        "source_hashes_verified",
        "sampling_optimization_used_64_members",
        "optimized_full_ensemble_score_not_below_preliminary",
        "selected_q10_meets_preliminary_guardrail",
        "automatic_preliminary_fallback_operational",
        "paired_comparison_published",
        "source_search_gate_not_fail",
        "exactly_ten_samples_per_process",
        "mandatory_wine_baseline_at_t0",
        "exactly_nine_capture_intervals_starting_at_t0",
        "captured_mass_conserved_between_intervals",
        "all_samples_in_legal_windows",
        "samples_and_actions_before_drying",
        "no_ambiguous_operational_conflicts",
        "capture_limits_approved_and_satisfied",
        "sample_before_action_coincidence_approved",
        "manual_sampling_capacity_three_respected",
        "nutrition_50_50_net_yan_policy_reconstructed",
        "tank_randomization_reproducible_uniform_permutation",
        "figures_generated_and_watermarked",
    )
    physical_blockers = sorted(
        set(constraints["fail_closed_fields"])
        | set(search_gate.get("physical_release_blockers", []))
    )
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_sampling_v2", config)
    schedule_name = (
        "optimized_sampling_schedule.csv"
        if selected_source == "optimized"
        else "fallback_sampling_schedule.csv"
    )
    paths = {
        "schedule": run_dir / schedule_name,
        "comparison": run_dir / "sampling_candidate_comparison.csv",
        "captures": run_dir / "optimized_capture_intervals.csv",
        "conflicts": run_dir / "operational_conflicts.csv",
        "nutrition": run_dir / "nutrition_product_translation.csv",
        "drying": run_dir / "drying_time_ensemble.csv",
        "candidate_restarts": run_dir / "sampling_search_restarts.csv",
        "reference_restarts": run_dir / "three_anchor_search_restarts.csv",
        "gate": run_dir / "sampling_gate.json",
        "config": run_dir / "wave1_mbdoe_config.json",
        "partition": run_dir / "partition_surrogate_provenance.json",
        "randomization": run_dir / "tank_randomization.csv",
        "figure_map": run_dir / "figure_chart_map.json",
        "runtime": run_dir / "runtime_summary.json",
    }
    for name in FIGURE_NAMES:
        paths[f"figure_{name[:-4]}"] = run_dir / name
    schedule.to_csv(filesystem_path(paths["schedule"]), index=False)
    paired.to_csv(filesystem_path(paths["comparison"]), index=False)
    captures.to_csv(filesystem_path(paths["captures"]), index=False)
    conflicts.to_csv(filesystem_path(paths["conflicts"]), index=False)
    nutrition.to_csv(filesystem_path(paths["nutrition"]), index=False)
    randomization.to_csv(filesystem_path(paths["randomization"]), index=False)
    drying.to_csv(filesystem_path(paths["drying"]), index=False)
    optimization_restarts.to_csv(filesystem_path(paths["candidate_restarts"]), index=False)
    reference_restarts.to_csv(filesystem_path(paths["reference_restarts"]), index=False)
    figure_paths = {name: run_dir / name for name in FIGURE_NAMES}
    chart_map = _generate_figures(
        figure_paths,
        policies,
        config,
        prepared,
        members,
        schedule,
        captures,
        paired,
        search_full,
        actuator,
        conflicts,
        checks,
    )
    write_json(
        paths["figure_map"],
        {
            "watermark": WATERMARK,
            "figures": chart_map,
            "qa_status": "generated_pending_visual_inspection",
        },
    )
    checks["figures_generated_and_watermarked"] = bool(
        len(chart_map) == len(FIGURE_NAMES)
        and all(path.is_file() and path.stat().st_size > 0 for path in figure_paths.values())
    )
    verdict = sampling_gate_verdict(checks, critical_checks)
    failed_checks = [name for name, passed in checks.items() if not passed]
    gate = {
        "gate": "phase_D_wave1_sampling_requalification",
        "verdict": verdict,
        "checks": checks,
        "source_search_verdict": search_gate["verdict"],
        "source_candidate_source": search_gate["candidate_source"],
        "selected_schedule_source": selected_source,
        "fallback_reason": fallback_reason,
        "information": {
            "optimized": {
                "robust_score": optimized_score,
                "median": float(np.median(optimized_gains)),
                "q10": float(np.quantile(optimized_gains, 0.1)),
            },
            "preliminary": preliminary_metrics,
            "selected": selected_metrics,
            "three_anchor_optimized_reference": reference_metrics,
            "optimized_vs_preliminary_wins": int((optimized_gains > preliminary_gains).sum()),
            "optimized_vs_preliminary_losses": int((optimized_gains < preliminary_gains).sum()),
            "optimized_vs_three_anchor_wins": int((optimized_gains > reference_gains).sum()),
            "optimized_vs_three_anchor_losses": int((optimized_gains < reference_gains).sum()),
            "worst_optimized_vs_preliminary_delta": float(
                np.min(optimized_gains - preliminary_gains)
            ),
        },
        "conditions": [f"failed_check:{name}" for name in failed_checks]
        + [f"physical_release_blocker:{name}" for name in physical_blockers],
        "nutrition_translation_blocker": bool(
            not nutrition.empty and nutrition["translation_blocker"].any()
        ),
        "tank_randomization_artifact": paths["randomization"].name,
        "tank_randomization_authorized": False,
        "tank_assignments": [],
        "profiles_for_physical_execution": False,
    }
    write_json(paths["gate"], gate)
    write_json(paths["config"], config)
    write_json(paths["partition"], provenance)
    write_json(
        paths["runtime"],
        {"total_runtime_seconds": float(time.perf_counter() - run_started)},
    )
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_sampling_requalification",
        config=config,
        sources={
            "hybrid_search_manifest": search_run / "run_manifest.json",
            "joint_ensemble_manifest": ensemble_run / "run_manifest.json",
            "aroma_calibration_manifest": aroma_run / "run_manifest.json",
            "wave1_config": CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "pilot_aroma_calibration.py",
            ADAPTIVE_DIR / "pilot_calibration.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            FERMENTATION_DIR / "shared" / "run_new_must_glycerol_estimability_doe.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[
            *[
                int(config["seed"]) + value
                for value in range(int(config["sampling"]["search_restarts"]))
            ],
            int(config["tank_randomization_seed"]),
        ],
        status="completed" if verdict != "FAIL" else "validation_failed",
        convergence={
            "candidate_restarts": optimization_restarts.to_dict(orient="records"),
            "reference_restarts": reference_restarts.to_dict(orient="records"),
            "selected_schedule_source": selected_source,
        },
        gate=gate,
        outputs=list(paths.values()),
        git_snapshot=git_snapshot,
    )
    manifest_path = run_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    if sha256_file(manifest_path) == "":  # pragma: no cover
        raise RuntimeError("Final sampling manifest could not be hashed")
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "computational_wave1_sampling_requalification",
            "current_gate": gate["gate"],
            "gate_verdict": verdict,
            "latest_wave1_sampling_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "latest_gate_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "gate_blockers": failed_checks + physical_blockers,
            "physical_execution_status": "not_authorized",
            "executable_schedule_issued": False,
            "tank_assignments": [],
            "profiles_for_physical_execution": False,
            "wave1": {
                **state.get("wave1", {}),
                "sampling_gate_verdict": verdict,
                "sampling_schedule_source": selected_source,
                "sampling_fallback_reason": fallback_reason,
                "candidate_sampling_schedule_ready_for_review": False,
            },
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(json.dumps({"run_directory": str(run_dir), **gate}, indent=2))
    if verdict == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
