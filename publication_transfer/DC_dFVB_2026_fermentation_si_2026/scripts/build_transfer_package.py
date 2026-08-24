"""Build a publication-safe Pilot 2025 data package for the DFVB article.

The source workbook uses industrial process identifiers.  This exporter never
writes those identifiers, sample IDs, or calendar timestamps.  Public tokens
are assigned in workbook order and are linked to private identifiers only by a
mapping maintained outside version control.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PACKAGE_ROOT.parents[1]
DATA_DIR = PACKAGE_ROOT / "data"
SOURCE_WORKBOOK = REPO_ROOT / "fermentation_model" / "data" / "Piloto 2025" / "Calibration_data_vl3.xlsx"
EXPERIMENT_REGISTRY = REPO_ROOT / "fermentation_model" / "campaigns" / "experiments.csv"

AMMONIA_TO_N_FACTOR = 0.82
BIOMASS_PROXY_G_L_PER_MILLION_CELLS_ML = 0.030

AROMA_COLUMNS = {
    "bencil_alcohol": "benzyl_alcohol",
    "benzaldehido": "benzaldehyde",
    "decanoato_de_etilo": "ethyl_decanoate",
    "hexil_acetate": "hexyl_acetate",
    "isoamil_acetate": "isoamyl_acetate",
    "octanoate_de_etilo": "ethyl_octanoate",
    "phenylethylacetate": "phenethyl_acetate",
    "Ethyl_Acetate": "ethyl_acetate",
}


def numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def has_reinoculation_annotation(frame: pd.DataFrame) -> bool:
    if "ID" not in frame:
        return False
    ids = frame["ID"].astype(str)
    return bool(ids.str.contains("Pre rein", case=False, na=False, regex=False).any())


def analytical_mask(frame: pd.DataFrame) -> pd.Series:
    candidates = [
        "Viability",
        "Peso Seco",
        "GLUCOSE",
        "FRUCTOSE",
        "PAN",
        "AMMONIA",
        "YAN",
        "GLYCEROL",
        "PYRUVIC ACID",
        "ACETALDEHIDO",
        "ETANOL",
    ] + [f"{name}_total" for name in AROMA_COLUMNS]
    present = [column for column in candidates if column in frame]
    if not present:
        return pd.Series(False, index=frame.index)
    return frame[present].apply(pd.to_numeric, errors="coerce").notna().any(axis=1)


def public_mapping(sheet_names: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    registry = pd.read_csv(EXPERIMENT_REGISTRY, dtype=str).fillna("")
    registry = registry[registry["campaign_id"].eq("pilot_2025")]
    pair_by_source = dict(zip(registry["experiment_id"], registry["replicate_group"]))
    pair_values = sorted({pair_by_source.get(sheet, "") for sheet in sheet_names})
    pair_values = [value for value in pair_values if value]
    pair_tokens = {value: f"C25_{idx:02d}" for idx, value in enumerate(pair_values, 1)}
    process_tokens = {sheet: f"P25_{idx:02d}" for idx, sheet in enumerate(sheet_names, 1)}
    campaign_tokens = {
        sheet: pair_tokens.get(pair_by_source.get(sheet, ""), "C25_UNMAPPED")
        for sheet in sheet_names
    }
    return process_tokens, campaign_tokens


def build() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    excel = pd.ExcelFile(SOURCE_WORKBOOK)
    sheet_names = list(excel.sheet_names)
    process_tokens, campaign_tokens = public_mapping(sheet_names)

    observation_frames: list[pd.DataFrame] = []
    aroma_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    registry_rows: list[dict[str, object]] = []
    all_source_ids = set(sheet_names)

    for sheet in sheet_names:
        raw = pd.read_excel(SOURCE_WORKBOOK, sheet_name=sheet)
        raw.columns = [str(column).strip() for column in raw.columns]
        token = process_tokens[sheet]
        campaign = campaign_tokens[sheet]
        incident = has_reinoculation_annotation(raw)
        include = not incident
        time_h = numeric(raw, "t")
        ammonia_compound = numeric(raw, "AMMONIA")
        pan_n = numeric(raw, "PAN")
        yan_reported = numeric(raw, "YAN")
        yan_reconstructed = pan_n + AMMONIA_TO_N_FACTOR * ammonia_compound

        observations = pd.DataFrame(
            {
                "process_token": token,
                "campaign_token": campaign,
                "time_h": time_h,
                "temperature_c": numeric(raw, "temperatura"),
                "density_kg_m3": numeric(raw, "densidad"),
                "viable_cell_concentration_million_cells_ml": numeric(raw, "Viability"),
                "viable_biomass_proxy_g_l_30pg_cell": numeric(raw, "Viability")
                * BIOMASS_PROXY_G_L_PER_MILLION_CELLS_ML,
                "dry_weight_g_l": numeric(raw, "Peso Seco"),
                "glucose_g_l": numeric(raw, "GLUCOSE"),
                "fructose_g_l": numeric(raw, "FRUCTOSE"),
                "pan_n_mg_l": pan_n,
                "ammonia_compound_mg_l": ammonia_compound,
                "ammonia_n_mg_l": AMMONIA_TO_N_FACTOR * ammonia_compound,
                "yan_reported_n_mg_l": yan_reported,
                "yan_reconstructed_n_mg_l": yan_reconstructed,
                "yan_closure_residual_mg_l": yan_reported - yan_reconstructed,
                "glycerol_g_l": numeric(raw, "GLYCEROL"),
                "pyruvic_acid_mg_l": numeric(raw, "PYRUVIC ACID"),
                "pyruvic_acid_2_mg_l": numeric(raw, "PYRUVIC ACID2"),
                "acetaldehyde_mg_l": numeric(raw, "ACETALDEHIDO"),
                "ethanol_g_l": numeric(raw, "ETANOL"),
            }
        )
        observations = observations[observations["time_h"].notna()].sort_values("time_h")
        if include:
            observation_frames.append(observations)

        for source_stem, species in AROMA_COLUMNS.items():
            total = numeric(raw, f"{source_stem}_total")
            condensate = numeric(raw, f"{source_stem}_condensado")
            for idx in raw.index[total.notna() | condensate.notna()]:
                total_value = float(total.loc[idx]) if pd.notna(total.loc[idx]) else np.nan
                condensate_value = float(condensate.loc[idx]) if pd.notna(condensate.loc[idx]) else np.nan
                retained = (
                    total_value - condensate_value
                    if np.isfinite(total_value) and np.isfinite(condensate_value)
                    else np.nan
                )
                aroma_rows.append(
                    {
                        "process_token": token,
                        "campaign_token": campaign,
                        "include_submission_dataset": include,
                        "time_h": float(time_h.loc[idx]),
                        "species": species,
                        "total_equivalent_mg_l": total_value,
                        "condensate_equivalent_mg_l": condensate_value,
                        "retained_wine_mg_l": retained,
                        "retained_qc": (
                            "negative_total_minus_condensate"
                            if np.isfinite(retained) and retained < -1e-9
                            else "ok"
                        ),
                    }
                )

        pulse = numeric(raw, "pulso_nut")
        for idx in raw.index[pulse.notna() & pulse.gt(0.0)]:
            event_rows.append(
                {
                    "process_token": token,
                    "campaign_token": campaign,
                    "include_submission_dataset": include,
                    "event_type": "initial_charge" if abs(float(time_h.loc[idx])) < 1e-9 else "process_pulse",
                    "time_h": float(time_h.loc[idx]),
                    "density_kg_m3_at_record": (
                        float(numeric(raw, "densidad").loc[idx])
                        if pd.notna(numeric(raw, "densidad").loc[idx])
                        else np.nan
                    ),
                    "recorded_amount": float(pulse.loc[idx]),
                    "recorded_field": "pulso_nut",
                    "recorded_unit_basis": "unresolved_in_source_workbook",
                    "model_use_policy": "do_not_convert_to_nitrogen_without_product_or_protocol_metadata",
                }
            )

        analytical = analytical_mask(raw) & time_h.notna()
        analytical_times = time_h[analytical]
        positive_pulses = pulse.notna() & pulse.gt(0.0)
        process_pulses = positive_pulses & time_h.gt(0.0)
        registry_rows.append(
            {
                "process_token": token,
                "campaign_token": campaign,
                "campaign": "pilot_2025_paired_vintage",
                "matrix": "natural_grape_must",
                "scale": "pilot",
                "temperature_input": "measured_sample_profile_piecewise_linear",
                "temperature_min_c": numeric(raw, "temperatura").min(),
                "temperature_mean_c": numeric(raw, "temperatura").mean(),
                "temperature_max_c": numeric(raw, "temperatura").max(),
                "analytical_start_h": analytical_times.min(),
                "analytical_end_h": analytical_times.max(),
                "n_source_rows": int(time_h.notna().sum()),
                "n_initial_charge_records": int((positive_pulses & time_h.eq(0.0)).sum()),
                "n_process_pulse_records": int(process_pulses.sum()),
                "used_in_current_repository_fit_or_model_selection": True,
                "include_submission_dataset": include,
                "submission_role": "development_calibration" if include else "excluded_process_incident",
                "exclusion_reason": (
                    "nonviable_initial_inoculum_and_reinoculation_requires_effective_time_restart"
                    if incident
                    else ""
                ),
                "independent_validation": False,
            }
        )

    observations = pd.concat(observation_frames, ignore_index=True)
    aromas = pd.DataFrame(aroma_rows)
    aromas = aromas[aromas["include_submission_dataset"].astype(bool)].reset_index(drop=True)
    events = pd.DataFrame(event_rows)
    events = events[events["include_submission_dataset"].astype(bool)].reset_index(drop=True)
    registry = pd.DataFrame(registry_rows)

    observations.to_csv(DATA_DIR / "process_observations.csv", index=False, lineterminator="\n")
    aromas.to_csv(DATA_DIR / "aroma_observations.csv", index=False, lineterminator="\n")
    events.to_csv(DATA_DIR / "nitrogen_events.csv", index=False, lineterminator="\n")
    registry.to_csv(DATA_DIR / "process_registry.csv", index=False, lineterminator="\n")

    closure = observations["yan_closure_residual_mg_l"].dropna()
    checks = [
        {
            "check": "source_process_count",
            "status": "pass" if len(sheet_names) == 8 else "fail",
            "value": len(sheet_names),
            "criterion": "8 source processes",
        },
        {
            "check": "submission_process_count",
            "status": "pass" if observations["process_token"].nunique() == 7 else "fail",
            "value": observations["process_token"].nunique(),
            "criterion": "7 homogeneous eligible processes",
        },
        {
            "check": "yan_exact_source_closure_max_abs_mg_l",
            "status": "pass" if closure.abs().max() <= 1e-9 else "fail",
            "value": closure.abs().max(),
            "criterion": "YAN = PAN + 0.82 * AMMONIA",
        },
        {
            "check": "negative_retained_aroma_count",
            "status": "pass" if (aromas["retained_qc"] != "ok").sum() == 0 else "warning",
            "value": int((aromas["retained_qc"] != "ok").sum()),
            "criterion": "review rather than clip total-condensate negatives",
        },
        {
            "check": "negative_dry_weight_count",
            "status": "warning" if observations["dry_weight_g_l"].lt(0.0).any() else "pass",
            "value": int(observations["dry_weight_g_l"].lt(0.0).sum()),
            "criterion": "preserve and treat through assay error/LOD policy; do not clip silently",
        },
        {
            "check": "pulse_unit_basis",
            "status": "warning",
            "value": "unresolved",
            "criterion": "confirm product dose and N-equivalent before model input",
        },
    ]
    pd.DataFrame(checks).to_csv(DATA_DIR / "quality_checks.csv", index=False, lineterminator="\n")

    manifest = pd.DataFrame(
        [
            {
                "logical_source": "pilot_2025_calibration_workbook",
                "source_filename": SOURCE_WORKBOOK.name,
                "source_size_bytes": SOURCE_WORKBOOK.stat().st_size,
                "source_sha256": sha256(SOURCE_WORKBOOK),
                "industrial_identifiers_exported": False,
            },
            {
                "logical_source": "experiment_campaign_registry",
                "source_filename": EXPERIMENT_REGISTRY.name,
                "source_size_bytes": EXPERIMENT_REGISTRY.stat().st_size,
                "source_sha256": sha256(EXPERIMENT_REGISTRY),
                "industrial_identifiers_exported": False,
            },
        ]
    )
    manifest.to_csv(DATA_DIR / "source_manifest.csv", index=False, lineterminator="\n")

    for output in DATA_DIR.glob("*.csv"):
        with output.open("r", encoding="utf-8", newline="") as handle:
            cells = {cell for row in csv.reader(handle) for cell in row}
        leaked = sorted(all_source_ids.intersection(cells))
        if leaked:
            raise RuntimeError(f"Industrial identifier leak in {output.name}: {leaked}")

    print(f"Built publication transfer package in {DATA_DIR}")
    print(f"Eligible processes: {observations['process_token'].nunique()}")
    print(f"Observation rows: {len(observations)}")
    print(f"Aroma rows: {len(aromas)}")
    print(f"Nitrogen event rows: {len(events)}")


if __name__ == "__main__":
    build()
