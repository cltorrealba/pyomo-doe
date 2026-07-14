from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import tempfile
from typing import Iterable

import numpy as np
import pandas as pd


SHARED_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SHARED_DIR.parent
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from shared.paths import DATA_DIR

NATURAL_MUST_FILE = DATA_DIR / "Laboratorio 2026" / "mosto_natural_xthiol.xlsx"
SYNTHETIC_MUST_FILE = DATA_DIR / "Laboratorio 2025" / "mosto_sintetico_vl3.xlsx"

NEW_MUST_SOURCES = {
    "natural": NATURAL_MUST_FILE,
    "synthetic": SYNTHETIC_MUST_FILE,
}

SECONDARY_METABOLITE_COLUMNS = {
    "PAN_mg_l": ["PAN", "Y15 PAN"],
    "NH4_mg_l": ["AMMONIA", "Y15 Amonio"],
    "pyruvic_acid": ["PYRUVIC ACID", "PYRUVIC_ACID", "Y15 Piruvic"],
    "pyruvic_acid_2": ["PYRUVIC ACID2"],
    "acetaldehyde": ["ACETALDEHIDO", "Y15 Acetaldehido"],
    "acetic_acid": ["Y15 Acetico", "ACETIC ACID", "ACETIC_ACID"],
}

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

FERMENTATION_SHEET_PREFIXES = {
    "natural": ("LAB",),
    "synthetic": ("MS",),
}

ETHANOL_DENSITY_G_ML = 0.78924
ETHANOL_G_L_PER_PERCENT_VV = 10.0 * ETHANOL_DENSITY_G_ML
CELL_MASS_PG_PER_CELL = 30.0
CELL_MASS_G_PER_CELL = CELL_MASS_PG_PER_CELL * 1e-12
MILLION_CELLS_ML_TO_KG_M3 = CELL_MASS_PG_PER_CELL * 1e-3
MG_L_TO_KG_M3 = 1e-3


@dataclass(frozen=True)
class DensitySugarFit:
    group: str
    source: str
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


def _first_available_numeric(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    values = pd.Series(np.nan, index=df.index, dtype=float)
    for column in columns:
        if column not in df.columns:
            continue
        candidate = pd.to_numeric(df[column], errors="coerce")
        values = values.where(values.notna(), candidate)
    return values


def _readable_excel_path(path: Path) -> Path:
    path = Path(path)
    try:
        with path.open("rb"):
            return path
    except PermissionError:
        copied = _save_open_excel_copy(path)
        if copied is not None:
            return copied
        raise


def _save_open_excel_copy(path: Path) -> Path | None:
    try:
        import win32com.client  # type: ignore
    except Exception:
        return None

    target = (
        Path(tempfile.gettempdir())
        / "pyomo_doe_excel_read_copies"
        / f"{Path(path).stem}_readcopy.xlsx"
    )
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        excel = win32com.client.GetActiveObject("Excel.Application")
        requested = Path(path).name.lower()
        for workbook in excel.Workbooks:
            if str(workbook.Name).lower() == requested:
                workbook.SaveCopyAs(str(target.resolve()))
                return target
    except Exception:
        return None
    return None


def _ethanol_observation(raw: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
    primary = _numeric(raw, "ETANOL")
    primary = primary.where(primary >= 0.0)

    fallback = _first_available_numeric(raw, ["Cf Alcolyzer Real (% v/v)", "Cf Alcolyzer (% v/v)"])
    fallback = fallback.where(fallback > 0.0)

    primary_values = primary.dropna()
    primary_is_g_l = bool(primary_values.gt(25.0).any())

    percent_vv = primary.where(primary.notna(), fallback)
    e_g_l = percent_vv * ETHANOL_G_L_PER_PERCENT_VV
    unit = pd.Series("", index=raw.index, dtype=object)
    source = pd.Series("", index=raw.index, dtype=object)

    if primary_is_g_l:
        e_g_l = primary.where(primary.notna(), fallback * ETHANOL_G_L_PER_PERCENT_VV)
        percent_vv = (primary / ETHANOL_G_L_PER_PERCENT_VV).where(
            primary.notna(),
            fallback,
        )
        unit = unit.mask(primary.notna(), "g_l")
    else:
        unit = unit.mask(primary.notna(), "percent_vv")

    unit = unit.mask(primary.isna() & fallback.notna(), "percent_vv")
    source = source.mask(primary.notna(), "ETANOL")
    source = source.mask(primary.isna() & fallback.notna(), "ALCOLYZER_FALLBACK")
    return percent_vv, e_g_l, source, unit


def _ethanol_percent_vv(raw: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    percent_vv, _e_g_l, source, _unit = _ethanol_observation(raw)
    return percent_vv, source


def fermentation_sheets(path: Path, medium: str) -> list[str]:
    path = _readable_excel_path(path)
    prefixes = FERMENTATION_SHEET_PREFIXES[medium]
    xls = pd.ExcelFile(path)
    sheets = []
    for sheet in xls.sheet_names:
        if not str(sheet).startswith(prefixes):
            continue
        header = pd.read_excel(path, sheet_name=sheet, nrows=0)
        columns = {str(col).strip() for col in header.columns}
        if "t" in columns or "time" in {col.lower() for col in columns}:
            sheets.append(str(sheet))
    return sheets


def load_new_must_data(sources: dict[str, Path] | None = None) -> pd.DataFrame:
    sources = NEW_MUST_SOURCES if sources is None else sources
    frames = []
    for medium, path in sources.items():
        path = Path(path)
        read_path = _readable_excel_path(path)
        for sheet in fermentation_sheets(read_path, medium):
            raw = _clean_columns(pd.read_excel(read_path, sheet_name=sheet))
            normalized = normalize_must_sheet(raw, medium=medium, batch=sheet, source_file=path)
            frames.append(normalized)
    if not frames:
        raise RuntimeError("No fermentation sheets were loaded from the new must data files.")
    data = pd.concat(frames, ignore_index=True)
    data = data.sort_values(["medium", "batch", "time_h"]).reset_index(drop=True)
    return data


def normalize_must_sheet(raw: pd.DataFrame, medium: str, batch: str, source_file: Path) -> pd.DataFrame:
    raw = _clean_columns(raw)
    out = pd.DataFrame(index=raw.index)
    out["medium"] = str(medium)
    out["batch"] = str(batch)
    out["source_file"] = str(source_file.name)
    out["sample_id"] = raw["ID"].astype(str) if "ID" in raw.columns else ""

    out["time_h"] = _first_available_numeric(raw, ["time", "t", "Horas"])
    out["temperature_c"] = _first_available_numeric(
        raw,
        ["temperature", "temperatura", "Temperatura_operacional_columna_Temperatura", "Temperatura"],
    )
    out["density"] = _first_available_numeric(raw, ["densidad", "Densidad"])
    out["brix"] = _first_available_numeric(raw, ["Brix"])
    out["DO_mg_l"] = _first_available_numeric(
        raw,
        ["DO", "Dissolved Oxygen", "Dissolved oxygen", "Oxigeno disuelto", "Oxígeno disuelto"],
    )
    out.loc[out["DO_mg_l"] < 0.0, "DO_mg_l"] = np.nan

    out["X_viable_mcells_ml"] = _first_available_numeric(raw, ["C_viable", "Viability"])
    out["X_total_mcells_ml"] = _first_available_numeric(raw, ["Oculyze Concentration", "C_total"])
    viable_for_dead = out["X_viable_mcells_ml"].interpolate(limit_direction="both")
    out["X_dead_direct_mcells_ml"] = (out["X_total_mcells_ml"] - out["X_viable_mcells_ml"]).clip(lower=0.0)
    out["X_dead_mcells_ml"] = (out["X_total_mcells_ml"] - viable_for_dead).clip(lower=0.0)
    out["X_dead_inferred"] = (
        out["X_total_mcells_ml"].notna()
        & out["X_viable_mcells_ml"].isna()
        & viable_for_dead.notna()
    )

    out["X_viable_kg_m3"] = out["X_viable_mcells_ml"] * MILLION_CELLS_ML_TO_KG_M3
    out["X_total_kg_m3"] = out["X_total_mcells_ml"] * MILLION_CELLS_ML_TO_KG_M3
    out["X_dead_kg_m3"] = out["X_dead_mcells_ml"] * MILLION_CELLS_ML_TO_KG_M3

    out["G_g_l"] = _first_available_numeric(raw, ["GLUCOSE", "GLUCOSE-320"])
    out["F_g_l"] = _first_available_numeric(raw, ["FRUCTOSE", "FRUCTOSE_inferida"])
    out["S_GF_g_l"] = out["G_g_l"] + out["F_g_l"]
    out.loc[out["G_g_l"].isna() | out["F_g_l"].isna(), "S_GF_g_l"] = np.nan
    out["S_Y15_g_l"] = _first_available_numeric(raw, ["Y15 Glucosa-Fructosa", "GLU-FRU-320"])
    out["S_observed_g_l"] = out["S_GF_g_l"].where(out["S_GF_g_l"].notna(), out["S_Y15_g_l"])

    out["YAN_mg_l"] = _first_available_numeric(raw, ["YAN"])
    out["N_kg_m3"] = out["YAN_mg_l"] * MG_L_TO_KG_M3
    out["N_pulse_mg_l"] = _first_available_numeric(raw, ["pulso_nut"]).fillna(0.0)
    out["N_pulse_kg_m3"] = out["N_pulse_mg_l"] * MG_L_TO_KG_M3
    out["PAN_mg_l"] = _first_available_numeric(raw, SECONDARY_METABOLITE_COLUMNS["PAN_mg_l"])
    out["NH4_mg_l"] = _first_available_numeric(raw, SECONDARY_METABOLITE_COLUMNS["NH4_mg_l"])
    out["PAN_kg_m3"] = out["PAN_mg_l"] * MG_L_TO_KG_M3
    out["NH4_kg_m3"] = out["NH4_mg_l"] * MG_L_TO_KG_M3
    out["YAN_components_mg_l"] = out["PAN_mg_l"] + out["NH4_mg_l"]
    out.loc[out["PAN_mg_l"].isna() | out["NH4_mg_l"].isna(), "YAN_components_mg_l"] = np.nan

    for normalized_column, raw_columns in SECONDARY_METABOLITE_COLUMNS.items():
        if normalized_column in {"PAN_mg_l", "NH4_mg_l"}:
            continue
        out[normalized_column] = _first_available_numeric(raw, raw_columns)

    out["glycerol_g_l"] = _first_available_numeric(raw, ["GLYCEROL"])
    (
        out["ethanol_percent_vv"],
        out["E_g_l"],
        out["ethanol_source"],
        out["ethanol_unit"],
    ) = _ethanol_observation(raw)

    for normalized_column, raw_columns in AROMA_COLUMNS.items():
        out[normalized_column] = _first_available_numeric(raw, raw_columns)

    for column in ["Condicion", "Condición", "Replica", "Réplica", "Variedad", "Levadura", "Adicion"]:
        if column in raw.columns:
            safe_column = (
                column.replace("ó", "o")
                .replace("é", "e")
                .replace("í", "i")
                .replace("á", "a")
                .replace("ú", "u")
            )
            out[safe_column] = raw[column]

    out = out.dropna(subset=["time_h"]).drop_duplicates(["medium", "batch", "time_h"], keep="first")
    out = out.sort_values("time_h").reset_index(drop=True)
    return out


def batch_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (medium, batch), group in data.groupby(["medium", "batch"], sort=True):
        row = {
            "medium": medium,
            "batch": batch,
            "n_rows": int(len(group)),
            "t_min_h": float(group["time_h"].min()),
            "t_max_h": float(group["time_h"].max()),
            "temperature_initial_c": _first_nonmissing(group["temperature_c"]),
            "temperature_min_c": float(group["temperature_c"].min()),
            "temperature_max_c": float(group["temperature_c"].max()),
            "G_initial_g_l": _first_nonmissing(group["G_g_l"]),
            "F_initial_g_l": _first_nonmissing(group["F_g_l"]),
            "S_initial_g_l": _first_nonmissing(group["S_observed_g_l"]),
            "YAN_initial_mg_l": _first_nonmissing(group["YAN_mg_l"]),
            "E_final_g_l": _last_nonmissing(group["E_g_l"]),
            "X_viable_initial_kg_m3": _first_nonmissing(group["X_viable_kg_m3"]),
            "X_dead_initial_kg_m3": _first_nonmissing(group["X_dead_kg_m3"]),
            "n_density": int(group["density"].notna().sum()),
            "n_gf": int(group["S_GF_g_l"].notna().sum()),
            "n_y15": int(group["S_Y15_g_l"].notna().sum()),
            "n_total_cells": int(group["X_total_mcells_ml"].notna().sum()),
            "n_ethanol": int(group["ethanol_percent_vv"].notna().sum()),
            "n_pulses": int((group["N_pulse_mg_l"].fillna(0.0) > 0.0).sum()),
            "total_N_pulse_mg_l": float(group["N_pulse_mg_l"].fillna(0.0).sum()),
        }
        rows.append(row)
    return pd.DataFrame(rows)


def _first_nonmissing(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.iloc[0]) if not values.empty else np.nan


def _last_nonmissing(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(values.iloc[-1]) if not values.empty else np.nan


def density_regression_points(data: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for source, column in [("G_plus_F", "S_GF_g_l"), ("Y15", "S_Y15_g_l")]:
        cols = ["medium", "batch", "time_h", "density", column]
        frame = data[cols].rename(columns={column: "sugar_g_l"}).copy()
        frame["sugar_source"] = source
        frames.append(frame)
    points = pd.concat(frames, ignore_index=True)
    points = points.dropna(subset=["density", "sugar_g_l"])
    points = points[np.isfinite(points["density"]) & np.isfinite(points["sugar_g_l"])]
    return points.reset_index(drop=True)


def fit_density_sugar_models(points: pd.DataFrame, by_medium: bool = True) -> pd.DataFrame:
    group_cols = ["sugar_source"]
    if by_medium:
        group_cols.insert(0, "medium")
    else:
        points = points.copy()
        points["medium"] = "pooled"
        group_cols.insert(0, "medium")

    rows = []
    for keys, group in points.groupby(group_cols, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        medium, source = keys
        fit = _fit_density_sugar_group(group, label=str(medium), source=str(source))
        if fit is not None:
            rows.append(fit.__dict__)
    return pd.DataFrame(rows)


def _fit_density_sugar_group(frame: pd.DataFrame, label: str, source: str) -> DensitySugarFit | None:
    clean = frame.dropna(subset=["density", "sugar_g_l"])
    if len(clean) < 3:
        return None
    x = clean["density"].to_numpy(dtype=float)
    y = clean["sugar_g_l"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, deg=1)
    pred = intercept + slope * x
    residual = y - pred
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else np.nan
    rmse = float(np.sqrt(np.mean(residual**2)))
    return DensitySugarFit(
        group=str(label),
        source=str(source),
        intercept=float(intercept),
        slope=float(slope),
        r2=float(r2),
        rmse_g_l=rmse,
        n_points=int(len(clean)),
    )


def add_density_based_sugar(data: pd.DataFrame, fit_table: pd.DataFrame, source: str = "G_plus_F") -> pd.DataFrame:
    out = data.copy()
    out["S_density_g_l"] = np.nan
    for row in fit_table.itertuples(index=False):
        if str(row.source) != source:
            continue
        mask = out["medium"].astype(str).eq(str(row.group))
        out.loc[mask, "S_density_g_l"] = float(row.intercept) + float(row.slope) * out.loc[mask, "density"]
    return out


def make_notebook_batch(ns: dict, data: pd.DataFrame, medium: str, batch: str):
    group = (
        data[(data["medium"].eq(medium)) & (data["batch"].eq(batch))]
        .sort_values("time_h")
        .reset_index(drop=True)
    )
    if group.empty:
        raise ValueError(f"No rows found for medium={medium!r}, batch={batch!r}.")
    time = group["time_h"].to_numpy(dtype=float)
    if len(np.unique(time)) != len(time):
        raise ValueError(f"Duplicate time points found for {medium}/{batch}.")

    temperature_c = group["temperature_c"].interpolate(limit_direction="both").to_numpy(dtype=float)
    nutrient_pulse_kg_m3 = group["N_pulse_kg_m3"].fillna(0.0).to_numpy(dtype=float)
    nutrient_pulse_rate_kg_m3_h = np.zeros_like(nutrient_pulse_kg_m3, dtype=float)
    for k in range(1, len(time)):
        dt = float(time[k] - time[k - 1])
        if dt <= 0.0:
            raise ValueError(f"Non-increasing time points found for {medium}/{batch}.")
        nutrient_pulse_rate_kg_m3_h[k] = float(nutrient_pulse_kg_m3[k]) / dt

    measurements = pd.DataFrame(index=time)
    measurements.index.name = "t"
    measurements["X"] = group["X_viable_kg_m3"].to_numpy(dtype=float)
    measurements["N"] = group["N_kg_m3"].to_numpy(dtype=float)
    measurements["G"] = group["G_g_l"].to_numpy(dtype=float)
    measurements["F"] = group["F_g_l"].to_numpy(dtype=float)
    measurements["E"] = group["E_g_l"].to_numpy(dtype=float)
    measurements = ns["apply_measurement_preprocessing"](measurements, time, nutrient_pulse_kg_m3)

    initial_guess = measurements.interpolate(limit_direction="both")
    initials = initial_guess.iloc[0].to_dict()
    if any(pd.isna(value) for value in initials.values()):
        missing = [name for name, value in initials.items() if pd.isna(value)]
        raise ValueError(f"{medium}/{batch} has no usable initial value for {missing}.")

    return ns["FermentationBatch"](
        batch_id=f"{medium}_{batch}",
        run_label=f"{medium}/{batch}",
        raw=group,
        time=time,
        temperature_c=temperature_c,
        nutrient_pulse_kg_m3=nutrient_pulse_kg_m3,
        nutrient_pulse_rate_kg_m3_h=nutrient_pulse_rate_kg_m3_h,
        measurements=measurements,
        initial_guess=initial_guess,
        initials=initials,
    )


def patch_notebook_measurement_error(ns: dict) -> None:
    """Allow simulation of external batches with the old notebook model.

    The calibration notebook derives measurement-error scales by reloading the
    historical Excel workbook. New natural/synthetic batches are not present in
    that file, so pure simulations need a local default scale provider.
    """

    default_error = {
        "X": 0.05,
        "N": 0.01,
        "G": 2.0,
        "F": 2.0,
        "E": 1.0,
    }
    if "MEASUREMENT_ERROR" in ns:
        default_error.update({key: float(value) for key, value in ns["MEASUREMENT_ERROR"].items()})

    def objective_measurement_error_for_batches(batch_ids=None, states=None):
        states_ = list(ns["STATE_LABELS"]) if states is None else list(states)
        return {state: float(default_error[state]) for state in states_}

    def objective_scale_lookup_for_batches(batch_ids=None, states=None):
        if batch_ids is None:
            batch_ids_ = []
        elif isinstance(batch_ids, (str, int)):
            batch_ids_ = [str(batch_ids)]
        else:
            batch_ids_ = [str(batch_id) for batch_id in batch_ids]
        states_ = list(ns["STATE_LABELS"]) if states is None else list(states)
        return {
            (batch_id, state): float(default_error[state])
            for batch_id in batch_ids_
            for state in states_
        }

    ns["objective_measurement_error_for_batches"] = objective_measurement_error_for_batches
    ns["objective_scale_lookup_for_batches"] = objective_scale_lookup_for_batches
