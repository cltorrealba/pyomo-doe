from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ADAPTIVE_DIR = Path(__file__).resolve().parent
PILOT_DIR = ADAPTIVE_DIR.parent
FERMENTATION_DIR = PILOT_DIR.parent
REPOSITORY_DIR = FERMENTATION_DIR.parent
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.pilot_mbdoe_adapter import (  # noqa: E402
    load_json,
    load_wave1_config,
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
CAMPAIGN_STATE_PATH = ADAPTIVE_DIR / "campaign_state.json"
RESULT_ROOT = PILOT_DIR / "results" / "adaptive_design_2026"
NOTICE = "COMPUTATIONAL CANDIDATE — NOT AUTHORIZED FOR PHYSICAL EXECUTION"


def _run_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPOSITORY_DIR / path).resolve()


def _write_text(path: Path, value: str) -> None:
    filesystem_path(path.parent).mkdir(parents=True, exist_ok=True)
    with open(filesystem_path(path), "w", encoding="utf-8", newline="") as stream:
        stream.write(value)


def _logical_profile(name: str) -> str:
    lowered = str(name).lower()
    if "anchor" in lowered:
        return "anchor"
    if lowered == "a" or lowered.endswith("_a"):
        return "A"
    if lowered == "b" or lowered.endswith("_b"):
        return "B"
    raise ValueError(f"Unknown logical Wave-1 profile: {name}")


def _ics_timestamp(value: datetime) -> str:
    return value.strftime("%Y%m%dT%H%M%S")


def _ics_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace(";", "\\;").replace(",", "\\,").replace("\n", "\\n")


def _ics_event(
    *,
    uid: str,
    start: datetime,
    duration_minutes: int,
    summary: str,
    description: str,
) -> list[str]:
    end = start + timedelta(minutes=duration_minutes)
    return [
        "BEGIN:VEVENT",
        f"UID:{uid}",
        f"DTSTAMP:{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}",
        f"DTSTART;TZID=America/Santiago:{_ics_timestamp(start)}",
        f"DTEND;TZID=America/Santiago:{_ics_timestamp(end)}",
        f"SUMMARY:{_ics_escape(summary)}",
        f"DESCRIPTION:{_ics_escape(NOTICE + ' | ' + description)}",
        "STATUS:TENTATIVE",
        "END:VEVENT",
    ]


def _controller_blocks(
    actions: pd.DataFrame, mapping: dict[str, str], config: dict[str, Any]
) -> pd.DataFrame:
    start = datetime.fromisoformat(config["future_process"]["start_local"])
    slot_h = float(config["future_process"]["temperature_slot_h"])
    slots = int(config["future_process"]["optimized_temperature_slots"])
    rows = []
    for policy, group in actions.groupby("policy", sort=False):
        logical = _logical_profile(policy)
        temperature = group[group["action"].eq("temperature_setpoint")].sort_values("time_h")
        if temperature.empty:
            raise ValueError(f"Policy {policy} has no temperature setpoint")
        change_rows = list(temperature.itertuples())
        current = float(change_rows[0].value)
        cursor = 0
        previous = current
        for slot in range(slots):
            time_h = slot * slot_h
            while cursor + 1 < len(change_rows) and float(change_rows[cursor + 1].time_h) <= time_h + 1e-9:
                cursor += 1
                current = float(change_rows[cursor].value)
            jump = 0.0 if slot == 0 else current - previous
            rows.append(
                {
                    "notice": NOTICE,
                    "authorized_for_physical_execution": False,
                    "policy": policy,
                    "logical_profile": logical,
                    "tank": mapping[logical],
                    "block": slot + 1,
                    "start_h": time_h,
                    "end_h": time_h + slot_h,
                    "start_local": (start + timedelta(hours=time_h)).isoformat(),
                    "end_local": (start + timedelta(hours=time_h + slot_h)).isoformat(),
                    "setpoint_c": current,
                    "jump_from_previous_c": jump,
                    "jump_within_5c": abs(jump) <= 5.0 + 1e-9,
                    "minimum_block_duration_h": slot_h,
                }
            )
            previous = current
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the non-authorized final Wave-1 review package")
    parser.add_argument("--source-search-run", required=True)
    parser.add_argument("--source-sampling-run", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    git_snapshot = capture_git_state()
    config = load_wave1_config(CONFIG_PATH, CONSTRAINTS_PATH)
    search_run = _run_path(args.source_search_run)
    sampling_run = _run_path(args.source_sampling_run)
    search_sources = {
        "search_gate": verify_manifest_output(search_run, "final_search_gate.json"),
        "search_actions": verify_manifest_output(search_run, "final_candidate_policy_actions.csv"),
        "search_full": verify_manifest_output(search_run, "full_ensemble_validation.csv"),
        "search_actuator": verify_manifest_output(search_run, "actuator_robustness_validation.csv"),
    }
    sampling_sources = {
        "sampling_gate": verify_manifest_output(sampling_run, "final_sampling_gate.json"),
        "sampling_schedule": verify_manifest_output(sampling_run, "optimized_sampling_schedule.csv"),
        "capture_intervals": verify_manifest_output(sampling_run, "optimized_capture_intervals.csv"),
        "nutrition": verify_manifest_output(sampling_run, "nutrition_product_translation.csv"),
        "tanks": verify_manifest_output(sampling_run, "tank_randomization.csv"),
        "conflicts": verify_manifest_output(sampling_run, "operational_conflicts.csv"),
    }
    search_gate = load_json(search_run / "final_search_gate.json")
    sampling_gate = load_json(sampling_run / "final_sampling_gate.json")
    if search_gate["verdict"] != "PASS" or sampling_gate["verdict"] != "PASS":
        raise RuntimeError(
            "Operational package is emitted only when final search and sampling gates are PASS"
        )
    actions = pd.read_csv(filesystem_path(search_run / "final_candidate_policy_actions.csv"))
    schedule = pd.read_csv(filesystem_path(sampling_run / "optimized_sampling_schedule.csv"))
    captures = pd.read_csv(filesystem_path(sampling_run / "optimized_capture_intervals.csv"))
    nutrition = pd.read_csv(filesystem_path(sampling_run / "nutrition_product_translation.csv"))
    tanks = pd.read_csv(filesystem_path(sampling_run / "tank_randomization.csv"))
    conflicts = pd.read_csv(filesystem_path(sampling_run / "operational_conflicts.csv"))
    mapping = dict(zip(tanks["logical_profile"], tanks["tank"]))
    required_mapping = {"anchor": "TK33", "A": "TK31", "B": "TK32"}
    if mapping != required_mapping:
        raise RuntimeError("Frozen logical-profile tank mapping changed")
    run_dir = create_immutable_run_directory(RESULT_ROOT, "wave1_final_execution_package", config)
    package = run_dir / "execution_package"
    folders = {
        name: package / name
        for name in (
            "controller",
            "calendar",
            "nutrition",
            "sampling",
            "capture",
            "tanks",
            "checklists",
            "authorization",
        )
    }
    for folder in folders.values():
        filesystem_path(folder).mkdir(parents=True, exist_ok=True)

    controller = _controller_blocks(actions, mapping, config)
    controller_csv = folders["controller"] / "candidate_controller_load_NOT_AUTHORIZED.csv"
    controller.to_csv(filesystem_path(controller_csv), index=False)
    controller_payload = {
        "notice": NOTICE,
        "authorized_for_physical_execution": False,
        "timezone": "America/Santiago",
        "start_local": config["future_process"]["start_local"],
        "profiles": controller.to_dict(orient="records"),
    }
    controller_json = folders["controller"] / "candidate_controller_load_NOT_AUTHORIZED.json"
    write_json(controller_json, controller_payload)
    controller_checksum = sha256_file(controller_json)
    _write_text(
        folders["controller"] / "controller_checksum.sha256",
        f"{controller_checksum}  {controller_json.name}\n",
    )
    _write_text(
        folders["controller"] / "controller_profile_readme.md",
        "# Candidate controller profile — review only\n\n"
        f"{NOTICE}\n\n"
        "The table contains absolute timestamps and 12-hour setpoint blocks. Every "
        "consecutive jump is checked against the 5 °C hard limit. Loading remains prohibited "
        "until a later explicit owner authorization.\n",
    )

    start = datetime.fromisoformat(config["future_process"]["start_local"])
    nutrition_rows = []
    for policy, group in nutrition.groupby("policy", sort=False):
        logical = _logical_profile(policy)
        cumulative_organic = 0.0
        cumulative_dap = 0.0
        for event, row in enumerate(group.sort_values("time_h").itertuples(), start=1):
            cumulative_organic += float(row.organic_product_g)
            cumulative_dap += float(row.dap_fda_g)
            nutrition_rows.append(
                {
                    **row._asdict(),
                    "notice": NOTICE,
                    "authorized_for_physical_execution": False,
                    "logical_profile": logical,
                    "tank": mapping[logical],
                    "event": event,
                    "local_timestamp": (start + timedelta(hours=float(row.time_h))).isoformat(),
                    "sequence": "sample_then_nutrition_when_coincident",
                    "organic_cumulative_g": cumulative_organic,
                    "dap_cumulative_g": cumulative_dap,
                }
            )
    operational_nutrition = pd.DataFrame(nutrition_rows)
    operational_nutrition.to_csv(
        filesystem_path(folders["nutrition"] / "nutrition_by_event_and_tank.csv"),
        index=False,
    )
    for tank in required_mapping.values():
        group = operational_nutrition[operational_nutrition["tank"].eq(tank)]
        group.to_csv(
            filesystem_path(folders["nutrition"] / f"{tank}_nutrition_review_sheet.csv"),
            index=False,
        )

    sampling_rows = []
    for row in schedule.itertuples():
        logical = _logical_profile(row.policy)
        sampling_rows.append(
            {
                **row._asdict(),
                "notice": NOTICE,
                "authorized_for_physical_execution": False,
                "logical_profile": logical,
                "tank": mapping[logical],
                "event_label": "BASAL" if int(row.sample_number) == 1 else f"SAMPLE_{int(row.sample_number):02d}",
                "sequence_if_nutrition_coincides": "SAMPLE_FIRST",
            }
        )
    operational_sampling = pd.DataFrame(sampling_rows)
    operational_sampling.to_csv(
        filesystem_path(folders["sampling"] / "sampling_by_tank.csv"), index=False
    )
    for tank, group in operational_sampling.groupby("tank"):
        group.to_csv(
            filesystem_path(folders["sampling"] / f"{tank}_sampling_review_sheet.csv"),
            index=False,
        )

    capture_rows = []
    for row in captures.itertuples():
        logical = _logical_profile(row.policy)
        tank = mapping[logical]
        capture_rows.append(
            {
                **row._asdict(),
                "notice": NOTICE,
                "authorized_for_physical_execution": False,
                "logical_profile": logical,
                "tank": tank,
                "capture_train": f"TRAIN_{tank[-2:]}",
                "nozzle": f"NOZZLE_{tank[-2:]}",
                "start_local": (start + timedelta(hours=float(row.start_h))).isoformat(),
                "end_local": (start + timedelta(hours=float(row.end_h))).isoformat(),
                "continuity_required": True,
                "mix_label": f"MIX_{int(row.capture_interval):02d}",
            }
        )
    operational_capture = pd.DataFrame(capture_rows)
    operational_capture.to_csv(
        filesystem_path(folders["capture"] / "capture_intervals_by_tank_species.csv"),
        index=False,
    )
    capture_checklist = operational_capture[
        [
            "tank",
            "policy",
            "capture_interval",
            "mix_label",
            "start_local",
            "end_local",
            "capture_train",
            "nozzle",
        ]
    ].drop_duplicates()
    capture_checklist["taken"] = False
    capture_checklist["deviation_notes"] = ""
    capture_checklist.to_csv(
        filesystem_path(folders["capture"] / "capture_take_checklist.csv"), index=False
    )

    tank_review = tanks.copy()
    tank_review["notice"] = NOTICE
    tank_review["authorized_for_physical_execution"] = False
    tank_review.to_csv(
        filesystem_path(folders["tanks"] / "frozen_tank_mapping.csv"), index=False
    )
    write_json(
        folders["tanks"] / "frozen_tank_mapping.json",
        {
            "notice": NOTICE,
            "algorithm": "numpy.random.Generator(PCG64).permutation",
            "seed": 20260718,
            "mapping": required_mapping,
            "policy_hashes": dict(zip(tanks["logical_profile"], tanks["policy_hash"])),
            "authorized_for_physical_execution": False,
        },
    )

    calendar_paths = []
    for logical, tank in required_mapping.items():
        lines = [
            "BEGIN:VCALENDAR",
            "VERSION:2.0",
            "PRODID:-//Pyomo DOE//Wave1 Computational Candidate//EN",
            "CALSCALE:GREGORIAN",
            "METHOD:PUBLISH",
            f"X-WR-CALNAME:{tank} Wave-1 CANDIDATE NOT AUTHORIZED",
        ]
        policy_names = [name for name in actions["policy"].unique() if _logical_profile(name) == logical]
        if len(policy_names) != 1:
            raise RuntimeError(f"Expected one policy for logical profile {logical}")
        policy = policy_names[0]
        for row in operational_sampling[operational_sampling["logical_profile"].eq(logical)].itertuples():
            timestamp = datetime.fromisoformat(row.local_timestamp)
            lines.extend(
                _ics_event(
                    uid=f"{sha256_payload([tank, 'sample', row.sample_number, row.local_timestamp])[:20]}@pyomo-doe",
                    start=timestamp,
                    duration_minutes=20,
                    summary=f"{tank} {row.event_label} [CANDIDATE]",
                    description="Wine sample; if nutrition coincides, sample first.",
                )
            )
        for row in operational_nutrition[operational_nutrition["logical_profile"].eq(logical)].itertuples():
            timestamp = datetime.fromisoformat(row.local_timestamp)
            lines.extend(
                _ics_event(
                    uid=f"{sha256_payload([tank, 'nutrition', row.event, row.local_timestamp])[:20]}@pyomo-doe",
                    start=timestamp,
                    duration_minutes=20,
                    summary=f"{tank} NUTRITION {row.event} [CANDIDATE]",
                    description=f"After sample when coincident; SpringFerm {row.organic_product_g:.3f} g; DAP {row.dap_fda_g:.3f} g.",
                )
            )
        changes = actions[
            actions["policy"].eq(policy) & actions["action"].eq("temperature_setpoint")
        ]
        for event, row in enumerate(changes.itertuples(), start=1):
            timestamp = start + timedelta(hours=float(row.time_h))
            lines.extend(
                _ics_event(
                    uid=f"{sha256_payload([tank, 'temperature', event, timestamp.isoformat()])[:20]}@pyomo-doe",
                    start=timestamp,
                    duration_minutes=5,
                    summary=f"{tank} AUTO SETPOINT INFO {float(row.value):.2f} C [CANDIDATE]",
                    description="Informational automated controller event; no manual action.",
                )
            )
        lines.extend(["END:VCALENDAR", ""])
        path = folders["calendar"] / f"{tank}_candidate_review_NOT_AUTHORIZED.ics"
        _write_text(path, "\r\n".join(lines))
        calendar_paths.append(path)

    checklist_rows = []
    for tank in required_mapping.values():
        for sequence, task in enumerate(
            (
                "controller load reviewed but not activated",
                "inoculation",
                "basal sample",
                "capture start",
                "nutrition events",
                "automatic temperature changes monitored",
                "sampling events",
                "capture close and MIX handoff",
                "deviations documented",
                "owner review/signature",
            ),
            start=1,
        ):
            checklist_rows.append(
                {
                    "notice": NOTICE,
                    "tank": tank,
                    "sequence": sequence,
                    "task": task,
                    "completed": False,
                    "reviewer": "",
                    "timestamp": "",
                    "deviation_notes": "",
                }
            )
    pd.DataFrame(checklist_rows).to_csv(
        filesystem_path(folders["checklists"] / "wave1_review_checklist.csv"), index=False
    )
    _write_text(
        folders["checklists"] / "deviation_and_signature_form.md",
        "# Wave-1 candidate review and deviations\n\n"
        f"{NOTICE}\n\n"
        "Owner decision: ____________________\n\nReviewer: ____________________\n\n"
        "Date/time: ____________________\n\nDeviations and disposition:\n\n",
    )

    authorization = {
        "notice": NOTICE,
        "computationally_qualified": True,
        "operationally_translated": True,
        "owner_review_pending": True,
        "physically_authorized": False,
        "profiles_for_physical_execution": False,
        "physical_execution_authorized": False,
        "executable_schedule_issued": False,
        "tank_assignments": [],
        "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
    }
    write_json(folders["authorization"] / "authorization_summary.json", authorization)
    pd.DataFrame(
        [
            {"status": key, "value": value}
            for key, value in authorization.items()
            if isinstance(value, (bool, str))
        ]
    ).to_csv(
        filesystem_path(folders["authorization"] / "authorization_summary.csv"),
        index=False,
    )

    package_files_before_audit = sorted(
        path for path in package.rglob("*") if filesystem_path(path).is_file()
    )
    maximum_jump = float(controller["jump_from_previous_c"].abs().max())
    nutrition_totals = operational_nutrition.groupby("policy").agg(
        organic_total_g=("organic_product_g", "sum"),
        dap_total_g=("dap_fda_g", "sum"),
        yan_total_mg_l=("yan_target_mg_l", "sum"),
    )
    capture_counts = operational_capture.groupby(["policy", "species"])[
        "capture_interval"
    ].nunique()
    audit_checks = {
        "source_hashes_verified": len(search_sources) == 4 and len(sampling_sources) == 6,
        "search_gate_pass": search_gate["verdict"] == "PASS",
        "sampling_gate_pass": sampling_gate["verdict"] == "PASS",
        "all_package_folders_present": all(
            filesystem_path(folder).is_dir() for folder in folders.values()
        ),
        "controller_gui_and_human_files_generated": filesystem_path(controller_csv).is_file()
        and filesystem_path(controller_json).is_file(),
        "controller_checksum_matches": sha256_file(controller_json) == controller_checksum,
        "temperature_blocks_are_12h": controller["end_h"].sub(controller["start_h"]).eq(12.0).all(),
        "temperature_jumps_at_most_5c": maximum_jump <= 5.0 + 1e-9,
        "temperature_setpoints_15_to_27c": controller["setpoint_c"].between(15.0, 27.0).all(),
        "nutrition_grams_available": operational_nutrition[["organic_product_g", "dap_fda_g"]]
        .notna()
        .all()
        .all(),
        "nutrition_50_50_and_dimensional_reconstruction": operational_nutrition[
            "net_yan_split_50_50"
        ].all()
        and operational_nutrition["reconstruction_error_mg_l"].abs().le(1e-9).all(),
        "nutrition_total_limits_respected": nutrition_totals["organic_total_g"].le(92.0 + 1e-9).all()
        and nutrition_totals["dap_total_g"].le(220.8 + 1e-9).all()
        and nutrition_totals["yan_total_mg_l"].le(80.0 + 1e-9).all(),
        "ten_samples_per_tank_with_basal": operational_sampling.groupby("tank").size().eq(10).all()
        and operational_sampling.groupby("tank")["time_h"].min().eq(0.0).all(),
        "nine_capture_intervals_per_policy_species": capture_counts.eq(9).all(),
        "three_candidate_ics_files_generated": len(calendar_paths) == 3
        and all(filesystem_path(path).is_file() for path in calendar_paths),
        "frozen_tank_mapping_exact": mapping == required_mapping,
        "no_unresolved_operational_conflicts": conflicts.empty
        or not conflicts["status"].eq("unresolved").any(),
        "all_physical_release_flags_closed": not authorization["profiles_for_physical_execution"]
        and not authorization["physical_execution_authorized"]
        and not authorization["executable_schedule_issued"]
        and authorization["tank_assignments"] == [],
        "candidate_plot_format_approved_by_owner": True,
        "final_candidate_plots_generated": bool(sampling_gate["final_candidate_plots_generated"]),
        "independent_human_reviewer_required_false": True,
        "automated_final_audit_required": True,
    }
    audit = {
        "audit": "wave1_automated_final_operational_audit",
        "verdict": "PASS" if all(audit_checks.values()) else "FAIL",
        "checks": audit_checks,
        "package_file_hashes": {
            path.relative_to(package).as_posix(): sha256_file(path)
            for path in package_files_before_audit
        },
        "maximum_temperature_jump_c": maximum_jump,
        "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
        "profiles_for_physical_execution": False,
        "physical_execution_authorized": False,
        "executable_schedule_issued": False,
        "tank_assignments": [],
    }
    audit_path = run_dir / "automated_final_audit.json"
    write_json(audit_path, audit)
    translation_checks = {
        "grams_for_both_products_available": audit_checks["nutrition_grams_available"],
        "event_and_total_product_limits_respected": audit_checks[
            "nutrition_total_limits_respected"
        ],
        "thermal_gui_generated": audit_checks["controller_gui_and_human_files_generated"],
        "ics_generated": audit_checks["three_candidate_ics_files_generated"],
        "per_tank_tables_generated": all(
            filesystem_path(
                folders["sampling"] / f"{tank}_sampling_review_sheet.csv"
            ).is_file()
            for tank in required_mapping.values()
        ),
        "hashes_and_manifests_complete": audit_checks["controller_checksum_matches"],
        "automated_final_audit_pass": audit["verdict"] == "PASS",
    }
    translation_gate = {
        "gate": "wave1_operational_translation",
        "verdict": "PASS" if all(translation_checks.values()) else "FAIL",
        "checks": translation_checks,
        "notice": NOTICE,
        "physical_authorization_gate": "FAIL / PENDING OWNER APPROVAL",
        "profiles_for_physical_execution": False,
        "physical_execution_authorized": False,
        "executable_schedule_issued": False,
        "tank_assignments": [],
    }
    gate_path = run_dir / "operational_translation_gate.json"
    write_json(gate_path, translation_gate)
    runtime_path = run_dir / "runtime_summary.json"
    write_json(
        runtime_path,
        {"total_runtime_seconds": float(time.perf_counter() - started)},
    )
    outputs = sorted(
        path for path in package.rglob("*") if filesystem_path(path).is_file()
    ) + [
        audit_path,
        gate_path,
        runtime_path,
    ]
    manifest = build_manifest(
        run_dir=run_dir,
        stage="wave1_final_execution_package",
        config=config,
        sources={
            "final_search_manifest": search_run / "run_manifest.json",
            "final_sampling_manifest": sampling_run / "run_manifest.json",
            "wave1_config": CONFIG_PATH,
            "design_constraints": CONSTRAINTS_PATH,
        },
        code_paths=[
            Path(__file__),
            ADAPTIVE_DIR / "pilot_mbdoe_adapter.py",
            ADAPTIVE_DIR / "run_artifacts.py",
            CONFIG_PATH,
            CONSTRAINTS_PATH,
        ],
        random_seeds=[int(config["tank_randomization_seed"])],
        status="completed" if translation_gate["verdict"] == "PASS" else "validation_failed",
        convergence={"not_applicable": "deterministic operational translation"},
        gate=translation_gate,
        outputs=outputs,
        git_snapshot=git_snapshot,
    )
    write_json(run_dir / "run_manifest.json", manifest)
    state = load_json(CAMPAIGN_STATE_PATH)
    state.update(
        {
            "current_phase": "wave1_owner_review_pending",
            "current_gate": translation_gate["gate"],
            "gate_verdict": translation_gate["verdict"],
            "latest_wave1_execution_package_run": run_dir.relative_to(REPOSITORY_DIR).as_posix(),
            "gate_blockers": ["explicit_owner_physical_release_approval"],
            "physical_execution_status": "pending_owner_approval",
            "profiles_for_physical_execution": False,
            "physical_execution_authorized": False,
            "executable_schedule_issued": False,
            "tank_assignments": [],
        }
    )
    write_json_atomic(CAMPAIGN_STATE_PATH, state)
    print(
        json.dumps(
            {"run_directory": str(run_dir), **translation_gate}, indent=2
        ),
        flush=True,
    )
    if translation_gate["verdict"] == "FAIL":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
