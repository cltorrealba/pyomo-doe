from __future__ import annotations

import argparse
from datetime import datetime, time, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_joint_campaign_doe"
OUT_DIR = SCRIPT_DIR / "results" / "operational_campaign_schedule"

DEFAULT_START = datetime(2026, 6, 15, 8, 0)
DEFAULT_BATCH_STARTS = (
    datetime(2026, 6, 15, 8, 0),
    datetime(2026, 6, 29, 8, 0),
    datetime(2026, 7, 13, 8, 0),
)
RECOMMENDED_BATCHES = {
    1: (
        "fructose_rich_iG_probe",
        "low_temp_growth_separation_plus_co2",
        "yan_saturation_scan",
    ),
    2: (
        "cold_synthesis_hot_stripping",
        "combined_stress_long_horizon",
        "glucose_rich_growth_yield",
    ),
    3: (
        "high_biomass_low_N_maintenance",
        "ethanol_initial_challenge",
        "glucose_pulse_after_N_depletion",
    ),
}
SAMPLE_TIMES = (time(9, 0), time(11, 0), time(14, 0), time(16, 0))
ACTION_TIMES = (time(9, 0), time(11, 0), time(14, 0), time(16, 0))
LATE_SAMPLE_TIME = time(9, 0)
PRE_WEEKEND_SAMPLE_TIME = time(16, 0)
DEFAULT_AUTOMATED_CHANNEL = "E"
FRONT_LOADED_H = 96.0
MAX_MINOR_SHIFT_H = 6.0
MAX_ACCEPTABLE_SHIFT_H = 12.0


def is_weekday(dt: datetime) -> bool:
    return dt.weekday() < 5


def combine(day: datetime, clock: time) -> datetime:
    return datetime(day.year, day.month, day.day, clock.hour, clock.minute)


def relative_hours(start: datetime, dt: datetime) -> float:
    return (dt - start).total_seconds() / 3600.0


def relative_to_datetime(start: datetime, t_h: float) -> datetime:
    return start + timedelta(hours=float(t_h))


def iter_slots(start: datetime, end: datetime, clocks: tuple[time, ...]) -> list[datetime]:
    day = datetime(start.year, start.month, start.day)
    rows = []
    while day <= end + timedelta(days=1):
        if is_weekday(day):
            for clock in clocks:
                dt = combine(day, clock)
                if start <= dt <= end:
                    rows.append(dt)
        day += timedelta(days=1)
    return sorted(rows)


def first_slot_at_or_after(target: datetime, clocks: tuple[time, ...]) -> datetime:
    day = datetime(target.year, target.month, target.day)
    for _ in range(21):
        if is_weekday(day):
            for clock in clocks:
                dt = combine(day, clock)
                if dt >= target:
                    return dt
        day += timedelta(days=1)
    raise RuntimeError(f"No feasible slot found after {target}.")


def nearest_slot(target: datetime, start: datetime, end: datetime, clocks: tuple[time, ...]) -> datetime:
    search_start = start
    search_end = max(end, target) + timedelta(days=4)
    slots = iter_slots(search_start, search_end, clocks)
    if not slots:
        raise RuntimeError("No feasible operational slots available.")
    return min(slots, key=lambda dt: abs((dt - target).total_seconds()))


def parse_temperature_profile(text: str) -> list[float]:
    return [float(item.strip()) for item in str(text).split(",") if item.strip()]


def shift_class(shift_h: float) -> str:
    abs_shift = abs(float(shift_h))
    if abs_shift <= MAX_MINOR_SHIFT_H:
        return "minor"
    if abs_shift <= MAX_ACCEPTABLE_SHIFT_H:
        return "moderate"
    return "large_reconsider"


def load_design() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = pd.read_csv(RESULTS_DIR / "aroma_campaign_selected.csv")
    protocols = pd.read_csv(RESULTS_DIR / "aroma_campaign_protocols.csv")
    pulses = pd.read_csv(RESULTS_DIR / "aroma_campaign_pulses.csv")
    selected = selected[["campaign_order", "candidate", "family", "horizon_h", "temperature_c"]]
    protocols = protocols.drop(columns=[col for col in protocols.columns if col in selected.columns and col != "candidate"], errors="ignore")
    design = selected.merge(protocols, on="candidate", how="left")
    return design, pulses, selected


def parse_batch_starts(text: str | None) -> tuple[datetime, ...]:
    if not text:
        return DEFAULT_BATCH_STARTS
    return tuple(datetime.fromisoformat(part.strip()) for part in text.split(";") if part.strip())


def assign_batches(design: pd.DataFrame, batch_mode: str, common_start: datetime, batch_starts: tuple[datetime, ...]) -> pd.DataFrame:
    design = design.copy()
    if batch_mode == "all_at_once":
        design["batch"] = 1
        design["start_datetime"] = common_start
        return design
    if batch_mode != "recommended":
        raise ValueError(f"Unknown batch mode {batch_mode!r}.")
    rows = []
    lookup = design.set_index("candidate", drop=False)
    for batch_id, names in RECOMMENDED_BATCHES.items():
        if batch_id > len(batch_starts):
            raise RuntimeError("Not enough batch start dates for the recommended batch plan.")
        for name in names:
            if name not in lookup.index:
                raise RuntimeError(f"Recommended candidate {name!r} is missing from the selected campaign.")
            row = lookup.loc[name].copy()
            row["batch"] = batch_id
            row["start_datetime"] = batch_starts[batch_id - 1]
            rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True)


def build_sampling_schedule(design: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for item in design.itertuples(index=False):
        start = datetime.fromisoformat(str(item.start_datetime)) if not isinstance(item.start_datetime, datetime) else item.start_datetime
        horizon_h = float(item.horizon_h)
        nominal_end = relative_to_datetime(start, horizon_h)
        final_dt = first_slot_at_or_after(nominal_end, SAMPLE_TIMES)
        all_slots = iter_slots(start, final_dt, SAMPLE_TIMES)
        for slot in all_slots:
            t_h = relative_hours(start, slot)
            is_final = abs((slot - final_dt).total_seconds()) <= 1.0
            front_loaded = t_h <= FRONT_LOADED_H + 1e-9
            pre_weekend = slot.weekday() == 4 and slot.time() == PRE_WEEKEND_SAMPLE_TIME and t_h <= horizon_h + 1e-9
            late_morning = t_h > FRONT_LOADED_H and slot.time() == LATE_SAMPLE_TIME
            if front_loaded or late_morning or pre_weekend or is_final:
                rows.append(
                    {
                        "campaign_order": int(item.campaign_order),
                        "batch": int(item.batch),
                        "candidate": item.candidate,
                        "start_datetime": start.isoformat(sep=" "),
                        "nominal_horizon_h": horizon_h,
                        "operational_final_h": relative_hours(start, final_dt),
                        "sample_time_h": t_h,
                        "calendar_time": slot.isoformat(sep=" "),
                        "date": slot.date().isoformat(),
                        "weekday": slot.strftime("%A"),
                        "clock": slot.strftime("%H:%M"),
                        "is_final": bool(is_final),
                        "sampling_phase": "final"
                        if is_final
                        else ("front_loaded" if front_loaded else ("pre_weekend" if pre_weekend else "late_morning")),
                        "measurements": "X,Xd,N,G,F,E,aromas_liquid" + (",condensate_final" if is_final else ""),
                    }
                )
    return pd.DataFrame(rows)


def build_pulse_schedule(pulses: pd.DataFrame, design: pd.DataFrame, automated_channel: str) -> pd.DataFrame:
    design_map = design.set_index("candidate")["horizon_h"].to_dict()
    batch_map = design.set_index("candidate")["batch"].to_dict()
    start_map = design.set_index("candidate")["start_datetime"].to_dict()
    rows = []
    for pulse in pulses.itertuples(index=False):
        nominal_h = float(pulse.time_h)
        candidate = str(pulse.candidate)
        horizon_h = float(design_map[candidate])
        start_raw = start_map[candidate]
        start = datetime.fromisoformat(str(start_raw)) if not isinstance(start_raw, datetime) else start_raw
        nominal_dt = relative_to_datetime(start, nominal_h)
        channel = str(pulse.channel)
        if abs(nominal_h) <= 1e-9:
            adjusted_dt = start
            rule = "initial_formulation_or_inoculation"
        elif automated_channel and channel == automated_channel:
            adjusted_dt = nominal_dt
            rule = f"automated_{automated_channel}_nominal_time"
        else:
            adjusted_dt = nearest_slot(nominal_dt, start, relative_to_datetime(start, horizon_h), ACTION_TIMES)
            rule = "nearest_workday_09_11_14_or_16"
        adjusted_h = relative_hours(start, adjusted_dt)
        shift_h = adjusted_h - nominal_h
        rows.append(
            {
                "candidate": candidate,
                "batch": int(batch_map[candidate]),
                "channel": channel,
                "amount": float(pulse.amount),
                "start_datetime": start.isoformat(sep=" "),
                "nominal_time_h": nominal_h,
                "nominal_calendar_time": nominal_dt.isoformat(sep=" "),
                "adjusted_time_h": adjusted_h,
                "adjusted_calendar_time": adjusted_dt.isoformat(sep=" "),
                "date": adjusted_dt.date().isoformat(),
                "weekday": adjusted_dt.strftime("%A"),
                "clock": adjusted_dt.strftime("%H:%M"),
                "shift_h": shift_h,
                "shift_class": shift_class(shift_h),
                "adjustment_rule": rule,
            }
        )
    return pd.DataFrame(rows)


def build_temperature_schedule(design: pd.DataFrame, temperature_mode: str) -> pd.DataFrame:
    rows = []
    for item in design.itertuples(index=False):
        start = datetime.fromisoformat(str(item.start_datetime)) if not isinstance(item.start_datetime, datetime) else item.start_datetime
        temps = parse_temperature_profile(item.temperature_c)
        if len(temps) <= 1:
            continue
        edges = np.linspace(0.0, float(item.horizon_h), len(temps) + 1)
        for segment_idx, change_h in enumerate(edges[:-1]):
            nominal_dt = relative_to_datetime(start, change_h)
            if abs(change_h) <= 1e-9:
                adjusted_dt = start
                rule = "initial_setpoint_at_inoculation"
            else:
                adjusted_dt = nearest_slot(nominal_dt, start, relative_to_datetime(start, float(item.horizon_h)), ACTION_TIMES)
                rule = "nearest_workday_10_or_16_if_manual"
            adjusted_h = relative_hours(start, adjusted_dt)
            rows.append(
                {
                    "candidate": item.candidate,
                    "batch": int(item.batch),
                    "segment": int(segment_idx + 1),
                    "temperature_c": float(temps[segment_idx]),
                    "start_datetime": start.isoformat(sep=" "),
                    "nominal_time_h": float(change_h),
                    "nominal_calendar_time": nominal_dt.isoformat(sep=" "),
                    "adjusted_time_h": adjusted_h,
                    "adjusted_calendar_time": adjusted_dt.isoformat(sep=" "),
                    "date": adjusted_dt.date().isoformat(),
                    "weekday": adjusted_dt.strftime("%A"),
                    "clock": adjusted_dt.strftime("%H:%M"),
                    "shift_h": adjusted_h - float(change_h),
                    "shift_class": shift_class(adjusted_h - float(change_h)),
                    "adjustment_rule": "automatic_programmed_nominal_time" if temperature_mode == "automatic" else rule,
                    "note": "Temperature controller is assumed pre-programmed; nominal time is retained."
                    if temperature_mode == "automatic"
                    else "Use nominal time if controller can be pre-programmed; use adjusted time if manual.",
                }
            )
    return pd.DataFrame(rows)


def daily_load(samples: pd.DataFrame, pulses: pd.DataFrame, temperature: pd.DataFrame, temperature_mode: str, automated_channel: str) -> pd.DataFrame:
    sample_load = samples.groupby("date").size().rename("manual_samples")
    pulse_mask = pulses["adjustment_rule"].ne("initial_formulation_or_inoculation")
    if automated_channel:
        pulse_mask &= pulses["channel"].ne(automated_channel)
    pulse_load = pulses[pulse_mask].groupby("date").size().rename("manual_pulses")
    if temperature_mode == "automatic":
        temp_load = pd.Series(dtype=int, name="manual_temperature_changes")
    else:
        temp_load = temperature[temperature["adjustment_rule"].ne("initial_setpoint_at_inoculation")].groupby("date").size().rename("manual_temperature_changes")
    out = pd.concat([sample_load, pulse_load, temp_load], axis=1).fillna(0).astype(int).reset_index()
    if "index" in out.columns and "date" not in out.columns:
        out = out.rename(columns={"index": "date"})
    out["total_manual_events"] = out[["manual_samples", "manual_pulses", "manual_temperature_changes"]].sum(axis=1)
    return out.sort_values("date")


def automation_channel_comparison(design: pd.DataFrame, pulse_table: pd.DataFrame, samples: pd.DataFrame, temperature: pd.DataFrame, temperature_mode: str) -> pd.DataFrame:
    rows = []
    for channel in ("", "N", "G", "F", "E"):
        pulses = build_pulse_schedule(pulse_table, design, channel)
        load = daily_load(samples, pulses, temperature, temperature_mode, channel)
        manual_pulses = pulses[pulses["adjustment_rule"].ne("initial_formulation_or_inoculation")]
        if channel:
            manual_pulses = manual_pulses[manual_pulses["channel"].ne(channel)]
        large = manual_pulses[manual_pulses["shift_class"].eq("large_reconsider")]
        rows.append(
            {
                "automated_channel": channel or "none",
                "large_manual_pulse_shifts": int(len(large)),
                "max_abs_manual_pulse_shift_h": float(manual_pulses["shift_h"].abs().max()) if not manual_pulses.empty else 0.0,
                "max_total_manual_events_per_day": int(load["total_manual_events"].max()) if not load.empty else 0,
                "max_manual_samples_per_day": int(load["manual_samples"].max()) if not load.empty else 0,
                "max_manual_pulses_per_day": int(load["manual_pulses"].max()) if not load.empty else 0,
            }
        )
    return pd.DataFrame(rows)


def write_report(
    design: pd.DataFrame,
    samples: pd.DataFrame,
    pulses: pd.DataFrame,
    temperature: pd.DataFrame,
    load: pd.DataFrame,
    automation_comparison: pd.DataFrame,
    start: datetime,
    batch_mode: str,
    temperature_mode: str,
    automated_channel: str,
) -> None:
    large_pulses = pulses[pulses["shift_class"].eq("large_reconsider") & pulses["channel"].ne(automated_channel)]
    large_temp = temperature[temperature["shift_class"].eq("large_reconsider")] if temperature_mode != "automatic" else pd.DataFrame()
    lines = [
        "# Operational campaign schedule assessment",
        "",
        f"Batch mode: `{batch_mode}`.",
        f"Reference start: `{start.isoformat(sep=' ')}`.",
        f"Temperature mode: `{temperature_mode}`.",
        f"Automated dosing channel: `{automated_channel or 'none'}`.",
        "",
        "Operational assumptions:",
        "",
        "- Lab days: Monday-Friday.",
        "- Manual sampling slots: 09:00, 11:00, 14:00, and 16:00.",
        "- Manual action slots for non-automated pulses: 09:00, 11:00, 14:00, and 16:00.",
        "- At most one dosing channel is automated. The recommended/default automated channel is `E`.",
        "- Temperature setpoint changes are assumed automatically programmed unless `temperature_mode = manual`.",
        "- Initial formulation/inoculation actions at `t = 0 h` are allowed as setup actions.",
        "- Sampling is front-loaded with up to four samples per fermentation per weekday through approximately 96 h, then one morning sample per weekday plus a final feasible terminal sample.",
        "- Online CO2 is not constrained by staffing; this schedule covers manual sampling and manual actions.",
        "",
        "## Daily manual load",
        "",
        load.to_markdown(index=False),
        "",
        "## Automated dosing channel comparison",
        "",
        automation_comparison.to_markdown(index=False),
        "",
        "## Batch assignment",
        "",
        design[["batch", "candidate", "start_datetime", "horizon_h", "temperature_c"]].to_markdown(index=False),
        "",
        "## Pulse adjustments requiring reconsideration",
        "",
        large_pulses.to_markdown(index=False) if not large_pulses.empty else "No pulse shift exceeded the large-shift threshold.",
        "",
        "## Manual temperature changes requiring reconsideration",
        "",
        large_temp.to_markdown(index=False) if not large_temp.empty else "No manual temperature-change shift exceeded the large-shift threshold.",
        "",
        "## Design-level sample counts",
        "",
        samples.groupby("candidate").agg(
            nominal_horizon_h=("nominal_horizon_h", "first"),
            operational_final_h=("operational_final_h", "first"),
            n_manual_samples=("sample_time_h", "count"),
            n_final_samples=("is_final", "sum"),
        ).reset_index().to_markdown(index=False),
        "",
        "## Recommended DOE handling",
        "",
        "The current Pyomo DoE workflow uses uniform relative sampling intervals. Before freezing the campaign, the next technical step is to pass these explicit operational sample/action times into the FIM model and recompute the posterior FIM.",
        "",
        "The acceptance criterion should be: preserve the hybrid campaign unless the operational calendar reduces the weakest variance reduction or minimum relative eigenvalue materially. If the degradation is large, redesign start dates, pulse times, or candidate selection under the calendar constraints.",
        "",
    ]
    (OUT_DIR / "operational_campaign_schedule_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_plots(samples: pd.DataFrame, pulses: pd.DataFrame, load: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for idx, (candidate, group) in enumerate(samples.groupby("candidate", sort=False)):
        ax.scatter(group["sample_time_h"], np.full(len(group), idx), s=18, color="tab:blue", alpha=0.75)
    for idx, (candidate, group) in enumerate(samples.groupby("candidate", sort=False)):
        pulse_group = pulses[pulses["candidate"].eq(candidate)]
        for pulse in pulse_group.itertuples(index=False):
            color = "tab:red" if pulse.shift_class == "large_reconsider" else "tab:orange"
            ax.scatter(float(pulse.adjusted_time_h), idx, marker="|", s=140, color=color)
    ax.set_yticks(range(samples["candidate"].nunique()))
    ax.set_yticklabels(list(samples.groupby("candidate", sort=False).groups.keys()), fontsize=8)
    ax.set_xlabel("relative time from each batch inoculation [h]")
    ax.set_title("Operational manual sampling and adjusted pulse times")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "operational_sampling_pulse_timeline.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 3.8))
    x = np.arange(len(load))
    bottom = np.zeros(len(load))
    for column, color in [
        ("manual_samples", "tab:blue"),
        ("manual_pulses", "tab:orange"),
        ("manual_temperature_changes", "tab:green"),
    ]:
        ax.bar(x, load[column], bottom=bottom, label=column, color=color)
        bottom += load[column].to_numpy(dtype=float)
    ax.set_xticks(x)
    ax.set_xticklabels(load["date"], rotation=45, ha="right")
    ax.set_ylabel("manual events/day")
    ax.set_title("Daily operational load if all fermentations start together")
    ax.legend(loc="best", fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "operational_daily_load.png", dpi=180)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Operational calendar assessment for the selected aroma campaign.")
    parser.add_argument("--start", default=DEFAULT_START.isoformat(sep=" "), help="Common inoculation datetime, e.g. '2026-06-15 08:00'.")
    parser.add_argument("--batch-mode", choices=["recommended", "all_at_once"], default="recommended")
    parser.add_argument(
        "--batch-starts",
        default=";".join(dt.isoformat(sep=" ") for dt in DEFAULT_BATCH_STARTS),
        help="Semicolon-separated batch starts for recommended mode.",
    )
    parser.add_argument("--temperature-mode", choices=["automatic", "manual"], default="automatic")
    parser.add_argument("--automated-channel", choices=["", "N", "G", "F", "E"], default=DEFAULT_AUTOMATED_CHANNEL)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = datetime.fromisoformat(str(args.start))
    batch_starts = parse_batch_starts(args.batch_starts)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    design, pulse_table, _selected = load_design()
    design = assign_batches(design, args.batch_mode, start, batch_starts)
    samples = build_sampling_schedule(design)
    pulses = build_pulse_schedule(pulse_table, design, args.automated_channel)
    temperature = build_temperature_schedule(design, args.temperature_mode)
    load = daily_load(samples, pulses, temperature, args.temperature_mode, args.automated_channel)
    automation_comparison = automation_channel_comparison(design, pulse_table, samples, temperature, args.temperature_mode)

    design.to_csv(OUT_DIR / "operational_campaign_protocols.csv", index=False)
    samples.to_csv(OUT_DIR / "operational_sampling_schedule.csv", index=False)
    pulses.to_csv(OUT_DIR / "operational_pulse_schedule.csv", index=False)
    temperature.to_csv(OUT_DIR / "operational_temperature_schedule.csv", index=False)
    load.to_csv(OUT_DIR / "operational_daily_load.csv", index=False)
    automation_comparison.to_csv(OUT_DIR / "automation_channel_comparison.csv", index=False)
    write_plots(samples, pulses, load)
    write_report(design, samples, pulses, temperature, load, automation_comparison, start, args.batch_mode, args.temperature_mode, args.automated_channel)

    print("Operational schedule written to:", OUT_DIR)
    print(load.to_string(index=False))
    large_pulses = pulses[pulses["shift_class"].eq("large_reconsider") & pulses["channel"].ne(args.automated_channel)]
    large_temp = temperature[temperature["shift_class"].eq("large_reconsider")]
    print("\nLarge pulse shifts:")
    print(large_pulses[["candidate", "channel", "nominal_time_h", "adjusted_time_h", "shift_h"]].to_string(index=False) if not large_pulses.empty else "none")
    print("\nLarge manual temperature shifts:")
    if args.temperature_mode == "automatic":
        print("none; temperature setpoints are assumed pre-programmed at nominal times")
    else:
        print(large_temp[["candidate", "segment", "nominal_time_h", "adjusted_time_h", "shift_h"]].to_string(index=False) if not large_temp.empty else "none")


if __name__ == "__main__":
    main()
