from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import tempfile
from typing import Iterable

import numpy as np
import pandas as pd


MODULE_DIR = Path(__file__).resolve().parent
PILOT_DIR = MODULE_DIR.parent if MODULE_DIR.name == "support" else MODULE_DIR
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
DATA_DIR = FERMENTATION_MODEL_DIR / "data" / "Piloto 2025"
CALIBRATION_FILE = DATA_DIR / "Calibration_data_vl3.xlsx"
CO2_DIR = DATA_DIR / "Sensores CO2"
INVALID_CO2_BATCHES = ("25150", "25151")
CO2_ACTIVATION_BATCHES = ("25171",)
CO2_ACTIVATION_TIME_OVERRIDES = {"25171": 114.0}
CO2_ACTIVATION_SAMPLE_PATTERNS = {"25171": ("Pre rein", "Pre reinoculo")}

ETHANOL_DENSITY_G_ML = 0.78924
CELL_MASS_PG_PER_CELL = 30.0
MILLION_CELLS_ML_TO_KG_M3 = CELL_MASS_PG_PER_CELL * 1e-3
MG_L_TO_KG_M3 = 1e-3

AROMA_COLUMNS = {
    "benzyl_alcohol_total": ["bencil_alcohol_total", "benzyl_alcohol_total"],
    "benzaldehyde_total": ["benzaldehido_total", "benzaldehyde_total"],
    "ethyl_decanoate_total": ["decanoato_de_etilo_total", "ethyl_decanoate_total"],
    "hexyl_acetate_total": ["hexil_acetate_total", "hexyl_acetate_total"],
    "isoamyl_acetate_total": ["isoamil_acetate_total", "isoamyl_acetate_total"],
    "ethyl_octanoate_total": ["octanoate_de_etilo_total", "ethyl_octanoate_total"],
    "phenylethyl_acetate_total": ["phenylethylacetate_total", "phenylethyl_acetate_total"],
    "ethyl_acetate_total": ["Ethyl_Acetate_total", "ethyl_acetate_total"],
}

AROMA_CONDENSATE_COLUMNS = {
    "benzyl_alcohol_condensate": ["bencil_alcohol_condensado", "benzyl_alcohol_condensate"],
    "benzaldehyde_condensate": ["benzaldehido_condensado", "benzaldehyde_condensate"],
    "ethyl_decanoate_condensate": ["decanoato_de_etilo_condensado", "ethyl_decanoate_condensate"],
    "hexyl_acetate_condensate": ["hexil_acetate_condensado", "hexyl_acetate_condensate"],
    "isoamyl_acetate_condensate": ["isoamil_acetate_condensado", "isoamyl_acetate_condensate"],
    "ethyl_octanoate_condensate": ["octanoate_de_etilo_condensado", "ethyl_octanoate_condensate"],
    "phenylethyl_acetate_condensate": ["phenylethylacetate_condensado", "phenylethyl_acetate_condensate"],
    "ethyl_acetate_condensate": ["Ethyl_Acetate_condensado", "ethyl_acetate_condensate"],
}


@dataclass(frozen=True)
class DensitySugarFit:
    intercept: float
    slope: float
    r2: float
    rmse_g_l: float
    n_points: int

    def predict(self, density: pd.Series | np.ndarray) -> np.ndarray:
        return self.intercept + self.slope * np.asarray(density, dtype=float)


def _clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(col).strip() for col in df.columns]
    return df


def _numeric(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def _readable_excel_path(path: Path) -> Path:
    """Return a readable workbook path, copying locked OneDrive files if needed."""
    path = Path(path)
    try:
        with path.open("rb"):
            return path
    except PermissionError:
        tmp_dir = Path(tempfile.gettempdir()) / "pyomo_doe_locked_excels"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_path = tmp_dir / path.name
        shutil.copy2(path, tmp_path)
        return tmp_path


def _first_available_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    values = pd.Series(np.nan, index=df.index, dtype=float)
    for column in columns:
        if column not in df.columns:
            continue
        candidate = pd.to_numeric(df[column], errors="coerce")
        values = values.where(values.notna(), candidate)
    return values


def fermentation_sheets(path: Path = CALIBRATION_FILE) -> list[str]:
    readable = _readable_excel_path(path)
    xls = pd.ExcelFile(readable)
    sheets: list[str] = []
    for sheet in xls.sheet_names:
        header = pd.read_excel(readable, sheet_name=sheet, nrows=0)
        columns = {str(col).strip() for col in header.columns}
        if {"t", "fecha_hora"}.intersection(columns):
            sheets.append(str(sheet))
    return sheets


def _ethanol_series(raw: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    """Return ethanol in g/L and %v/v.

    Pilot data range reaches about 95-98 in the ETANOL column. Treating that
    as %v/v is physically impossible for wine, while 95 g/L is plausible.
    """
    primary = _numeric(raw, "ETANOL").where(lambda x: x >= 0.0)
    e_g_l = primary.copy()
    e_percent_vv = e_g_l / (10.0 * ETHANOL_DENSITY_G_ML)
    unit = pd.Series("", index=raw.index, dtype=object).mask(primary.notna(), "g_l")
    return e_g_l, e_percent_vv, unit


def normalize_sheet(raw: pd.DataFrame, batch: str, source_file: Path = CALIBRATION_FILE) -> pd.DataFrame:
    raw = _clean_columns(raw)
    out = pd.DataFrame(index=raw.index)
    out["system"] = "pilot_2025"
    out["medium"] = "natural"
    out["scale"] = "pilot"
    out["batch"] = str(batch)
    out["source_file"] = str(source_file.name)
    out["sample_id"] = raw["ID"].astype(str) if "ID" in raw.columns else ""
    out["fecha_hora"] = pd.to_datetime(raw["fecha_hora"], errors="coerce") if "fecha_hora" in raw.columns else pd.NaT
    out["time_h"] = _first_available_numeric(raw, ["t", "time", "Horas"])
    out["temperature_c"] = _first_available_numeric(raw, ["temperatura", "temperature", "Temperatura"])
    out["density"] = _first_available_numeric(raw, ["densidad", "Densidad"])

    out["X_viable_mcells_ml"] = _first_available_numeric(raw, ["Viability", "C_viable"])
    out["X_viable_kg_m3"] = out["X_viable_mcells_ml"] * MILLION_CELLS_ML_TO_KG_M3
    out["X_dry_weight_g_l"] = _first_available_numeric(raw, ["Peso Seco", "Peso seco", "dry_weight"])

    out["G_g_l"] = _first_available_numeric(raw, ["GLUCOSE", "Glucose"])
    out["F_g_l"] = _first_available_numeric(raw, ["FRUCTOSE", "Fructose"])
    out["S_GF_g_l"] = out["G_g_l"] + out["F_g_l"]
    out.loc[out["G_g_l"].isna() | out["F_g_l"].isna(), "S_GF_g_l"] = np.nan

    out["PAN_mg_l"] = _first_available_numeric(raw, ["PAN"])
    out["NH4_mg_l"] = _first_available_numeric(raw, ["AMMONIA", "NH4", "Amonio"])
    out["YAN_mg_l"] = _first_available_numeric(raw, ["YAN"])
    out["N_kg_m3"] = out["YAN_mg_l"] * MG_L_TO_KG_M3
    out["PAN_kg_m3"] = out["PAN_mg_l"] * MG_L_TO_KG_M3
    out["NH4_kg_m3"] = out["NH4_mg_l"] * MG_L_TO_KG_M3
    out["N_pulse_mg_l"] = _first_available_numeric(raw, ["pulso_nut"]).fillna(0.0)
    out["N_pulse_kg_m3"] = out["N_pulse_mg_l"] * MG_L_TO_KG_M3

    out["glycerol_g_l"] = _first_available_numeric(raw, ["GLYCEROL", "Glycerol"])
    out["pyruvic_acid_mg_l"] = _first_available_numeric(raw, ["PYRUVIC ACID", "PYRUVIC_ACID"])
    out["pyruvic_acid_2_mg_l"] = _first_available_numeric(raw, ["PYRUVIC ACID2", "PYRUVIC_ACID2"])
    out["acetaldehyde_mg_l"] = _first_available_numeric(raw, ["ACETALDEHIDO", "Acetaldehido"])
    out["acetic_acid_g_l"] = _first_available_numeric(raw, ["ACETIC ACID", "ACETIC_ACID", "acetic_acid"])
    out["DO_mg_l"] = _first_available_numeric(raw, ["DO", "Dissolved Oxygen", "dissolved_oxygen", "O2"])

    out["E_g_l"], out["E_percent_vv"], out["E_unit_assumption"] = _ethanol_series(raw)

    for normalized, candidates in AROMA_COLUMNS.items():
        out[normalized] = _first_available_numeric(raw, candidates)
    for normalized, candidates in AROMA_CONDENSATE_COLUMNS.items():
        out[normalized] = _first_available_numeric(raw, candidates)

    out = out[out["time_h"].notna() | out["fecha_hora"].notna()].copy()
    out = out.sort_values(["time_h", "fecha_hora"], na_position="last").reset_index(drop=True)
    return out


def load_pilot_calibration_data(path: Path = CALIBRATION_FILE) -> pd.DataFrame:
    readable = _readable_excel_path(path)
    frames = []
    for sheet in fermentation_sheets(readable):
        raw = pd.read_excel(readable, sheet_name=sheet)
        frames.append(normalize_sheet(raw, batch=sheet, source_file=path))
    if not frames:
        raise RuntimeError(f"No pilot fermentation sheets found in {path}.")
    return pd.concat(frames, ignore_index=True).sort_values(["batch", "time_h"]).reset_index(drop=True)


def load_co2_sensor_file(path: Path, calibration_data: pd.DataFrame | None = None) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0, header=None)
    batch = path.stem
    table = raw.iloc[3:, :2].copy()
    table.columns = ["timestamp", "co2_raw"]
    table["timestamp"] = pd.to_datetime(table["timestamp"], format="%m/%d/%y %H:%M:%S", errors="coerce")
    table["co2_raw"] = pd.to_numeric(table["co2_raw"], errors="coerce")
    table = table.dropna(subset=["timestamp", "co2_raw"]).reset_index(drop=True)
    table["batch"] = batch
    table["source_file"] = path.name
    if calibration_data is not None and not calibration_data.empty:
        starts = calibration_data.groupby("batch")["fecha_hora"].min()
        if batch in starts.index and pd.notna(starts.loc[batch]):
            start = pd.Timestamp(starts.loc[batch])
            table["time_h"] = (table["timestamp"] - start).dt.total_seconds() / 3600.0
        else:
            table["time_h"] = np.nan
    else:
        table["time_h"] = np.nan
    return table


def load_all_co2_sensor_data(calibration_data: pd.DataFrame | None = None, co2_dir: Path = CO2_DIR) -> pd.DataFrame:
    frames = []
    for path in sorted(co2_dir.glob("*.xlsx")):
        frames.append(load_co2_sensor_file(path, calibration_data=calibration_data))
    if not frames:
        return pd.DataFrame(columns=["timestamp", "co2_raw", "batch", "source_file", "time_h"])
    return pd.concat(frames, ignore_index=True).sort_values(["batch", "timestamp"]).reset_index(drop=True)


def detect_co2_activation_time(
    group: pd.DataFrame,
    threshold: float = 0.5,
    window_h: float = 24.0,
    min_points: int = 360,
    min_fraction_above: float = 0.50,
) -> float:
    """Detect the first sustained CO2 activation time for a sensor batch."""
    clean = group.dropna(subset=["time_h", "co2_raw"]).sort_values("time_h")
    for _, row in clean[clean["co2_raw"].gt(threshold)].iterrows():
        t0 = float(row["time_h"])
        window = clean[clean["time_h"].between(t0, t0 + window_h, inclusive="both")]
        if (
            len(window) >= min_points
            and float(window["co2_raw"].median()) > threshold
            and float(window["co2_raw"].gt(threshold).mean()) >= min_fraction_above
        ):
            return t0
    return float("nan")


def activation_time_from_calibration(calibration_data: pd.DataFrame | None, batch: str) -> float:
    """Find an operational activation time from annotated calibration samples."""
    if calibration_data is None or calibration_data.empty:
        return float("nan")
    batch_rows = calibration_data[calibration_data["batch"].astype(str).eq(str(batch))].copy()
    if batch_rows.empty or "sample_id" not in batch_rows.columns:
        return float("nan")
    patterns = CO2_ACTIVATION_SAMPLE_PATTERNS.get(str(batch), tuple())
    if not patterns:
        return float("nan")
    sample_id = batch_rows["sample_id"].astype(str)
    mask = pd.Series(False, index=batch_rows.index)
    for pattern in patterns:
        mask = mask | sample_id.str.contains(pattern, case=False, na=False, regex=False)
    hits = batch_rows.loc[mask, "time_h"].dropna()
    if hits.empty:
        return float("nan")
    return float(hits.iloc[0])


def curate_co2_sensor_data(
    raw_co2: pd.DataFrame,
    calibration_data: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return CO2 data intended for calibration and a curation decision table.

    The pilot 25150/25151 CO2 files are excluded from calibration. Batch 25171
    is trimmed to the sustained activation point because the pre-activation
    segment appears to correspond to a process/sensor artifact rather than the
    real fermentation start.
    """
    rows: list[dict[str, object]] = []
    frames: list[pd.DataFrame] = []
    if raw_co2.empty:
        return raw_co2.copy(), pd.DataFrame()

    for batch, group in raw_co2.groupby("batch", sort=True):
        batch_str = str(batch)
        decision = {
            "batch": batch_str,
            "raw_rows": int(len(group)),
            "used_for_calibration": True,
            "reason": "kept",
            "activation_time_h": 0.0,
            "curated_rows": 0,
        }
        if batch_str in INVALID_CO2_BATCHES:
            decision.update(
                {
                    "used_for_calibration": False,
                    "reason": "excluded: CO2 file not usable for this batch",
                    "activation_time_h": np.nan,
                }
            )
            rows.append(decision)
            continue

        curated = group.copy()
        curated = curated[curated["time_h"].notna() & curated["time_h"].ge(0.0)].copy()
        activation_time = 0.0
        if batch_str in CO2_ACTIVATION_BATCHES:
            annotated = activation_time_from_calibration(calibration_data, batch_str)
            if np.isfinite(annotated):
                detected = annotated
                decision["reason"] = "trimmed to annotated process restart sample"
            elif batch_str in CO2_ACTIVATION_TIME_OVERRIDES:
                detected = CO2_ACTIVATION_TIME_OVERRIDES[batch_str]
                decision["reason"] = "trimmed to operational CO2/process activation override"
            else:
                detected = detect_co2_activation_time(curated)
                decision["reason"] = "trimmed to sustained CO2 activation"
            if np.isfinite(detected):
                activation_time = float(detected)
                curated = curated[curated["time_h"].ge(activation_time)].copy()
            else:
                decision["reason"] = "kept from t>=0; activation not detected"
        curated["time_h_original"] = curated["time_h"]
        curated["co2_activation_time_h"] = activation_time
        curated["time_h_effective"] = curated["time_h_original"] - activation_time
        curated["co2_curation_reason"] = str(decision["reason"])
        decision["activation_time_h"] = activation_time
        decision["curated_rows"] = int(len(curated))
        rows.append(decision)
        if not curated.empty:
            frames.append(curated)

    out = (
        pd.concat(frames, ignore_index=True).sort_values(["batch", "timestamp"]).reset_index(drop=True)
        if frames
        else pd.DataFrame(columns=list(raw_co2.columns) + ["time_h_original", "co2_activation_time_h", "time_h_effective", "co2_curation_reason"])
    )
    return out, pd.DataFrame(rows)


def fit_density_total_sugar(data: pd.DataFrame) -> DensitySugarFit:
    df = data[["density", "S_GF_g_l"]].dropna()
    if len(df) < 3:
        return DensitySugarFit(np.nan, np.nan, np.nan, np.nan, len(df))
    x = df["density"].to_numpy(dtype=float)
    y = df["S_GF_g_l"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, deg=1)
    pred = intercept + slope * x
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    return DensitySugarFit(float(intercept), float(slope), float(r2), rmse, int(len(df)))


def observation_coverage(data: pd.DataFrame) -> pd.DataFrame:
    variables = [
        "temperature_c",
        "density",
        "S_GF_g_l",
        "X_viable_kg_m3",
        "X_dry_weight_g_l",
        "YAN_mg_l",
        "PAN_mg_l",
        "NH4_mg_l",
        "glycerol_g_l",
        "pyruvic_acid_mg_l",
        "acetaldehyde_mg_l",
        "acetic_acid_g_l",
        "DO_mg_l",
        "E_g_l",
        "ethyl_acetate_total",
        "ethyl_acetate_condensate",
        "isoamyl_acetate_total",
        "isoamyl_acetate_condensate",
        "ethyl_octanoate_total",
        "ethyl_octanoate_condensate",
    ]
    rows = []
    for batch, group in data.groupby("batch"):
        row = {
            "batch": batch,
            "n_rows": int(len(group)),
            "t_min_h": float(group["time_h"].min()),
            "t_max_h": float(group["time_h"].max()),
        }
        for variable in variables:
            if variable in group.columns:
                row[f"n_{variable}"] = int(group[variable].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("batch").reset_index(drop=True)
