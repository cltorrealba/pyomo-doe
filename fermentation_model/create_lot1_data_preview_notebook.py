from pathlib import Path
import textwrap

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "fermentation_model" / "fermentation_lot1_data_preview.ipynb"


def md(text: str):
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip())


def code(text: str):
    return nbf.v4.new_code_cell(textwrap.dedent(text).strip())


nb = nbf.v4.new_notebook()
nb["metadata"] = {
    "kernelspec": {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    },
    "language_info": {
        "name": "python",
        "pygments_lexer": "ipython3",
    },
}

nb["cells"] = [
    md(
        """
        # DOE Lot 1 data preview

        This notebook inspects the first MBDoE experimental lot and builds a traceable
        preview dataset by fermentation process.

        Data sources:

        - `CO2_FILT_F*.csv`: online filtered CO2 flow. The column used here is
          `flow_filt_sccm`.
        - `Temp_F*.csv`: online temperature and setpoint.
        - `Oculyze*.csv`: total cell concentration, viability, density and budding
          index.
        - `Y15_lote1/*.txt`: chemical analyses exported from Y15 sessions.

        Time convention:

        - For each fermentation, `t = 0 h` is the process timestamp of the first
          Oculyze sample.
        - The process timestamp is taken from the Oculyze `Description` field,
          because `Date` and `Time` are the later Oculyze measurement timestamps.
        - Y15 chemistry is aligned to Oculyze sample IDs when possible. If a Y15
          sample has no matching Oculyze timestamp, the notebook first infers time
          from the lot sampling sequence using the sample number. If that is not
          possible, the session timestamp encoded in the TXT filename is used as
          fallback and flagged.
        - Online CO2 is corrected for sampling volume losses. The raw sensor signal
          is a total reactor flow; removing homogeneous broth reduces the liquid
          volume and therefore the total fermenting inventory. The notebook reports
          both a volumetric signal and a 2 L-equivalent signal.
        """
    ),
    code(
        """
        from pathlib import Path
        import re
        import warnings

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from IPython.display import display

        pd.set_option("display.max_columns", 120)
        pd.set_option("display.width", 180)

        ROOT = Path.cwd()
        if not (ROOT / "fermentation_model" / "data").exists() and (ROOT / "data").exists():
            ROOT = ROOT.parent
        DATA_DIR = ROOT / "fermentation_model" / "data" / "Laboratorio 2026" / "DOE_Lote_1"
        RESULTS_DIR = ROOT / "fermentation_model" / "results" / "lot1_data_preview"
        FIG_DIR = RESULTS_DIR / "figures"
        PROCESSED_DIR = RESULTS_DIR / "processed"
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

        CELL_DRY_MASS_PG = 30.0
        G_L_PER_MILLION_CELLS_ML = CELL_DRY_MASS_PG * 1e-3

        PROCESS_META = {
            "F1": {
                "sample_prefix": "DOE-LAB001",
                "design_name": "synthetic_lit_SM410_18C_highN_ester",
                "folder_glob": "DOE26-F01_*",
            },
            "F2": {
                "sample_prefix": "DOE-LAB002",
                "design_name": "synthetic_high_biomass_low_N_maintenance",
                "folder_glob": "DOE26-F02_*",
            },
            "F3": {
                "sample_prefix": "DOE-LAB003",
                "design_name": "synthetic_fructose_rich_glucose",
                "folder_glob": "DOE26-F03_*",
            },
        }
        PREFIX_TO_PROCESS = {v["sample_prefix"]: k for k, v in PROCESS_META.items()}

        ETHANOL_TEMPLATE_XLSX_DATA = (
            ROOT
            / "fermentation_model"
            / "data"
            / "Laboratorio 2026"
            / "DOE_Lote_1"
            / "DOE_Lote_1_ethanol_manual_entry.xlsx"
        )
        ETHANOL_TEMPLATE_XLSX_FALLBACK = (
            ROOT
            / "fermentation_model"
            / "results"
            / "lot1_data_preview"
            / "manual_entry_templates"
            / "DOE_Lote_1_ethanol_manual_entry.xlsx"
        )
        ETHANOL_TEMPLATE_XLSX = (
            ETHANOL_TEMPLATE_XLSX_DATA
            if ETHANOL_TEMPLATE_XLSX_DATA.exists()
            else ETHANOL_TEMPLATE_XLSX_FALLBACK
        )
        LOT1_PLANNED_SCHEDULE_CSV = ROOT / "db_20260609" / "rec" / "lot_1_schedule.csv"
        ACTUAL_INPUT_EVENTS_CSV = (
            ROOT
            / "fermentation_model"
            / "results"
            / "lot1_data_preview"
            / "manual_entry_templates"
            / "DOE_Lote_1_operational_inputs_actual.csv"
        )
        ETHANOL_DENSITY_G_ML = 0.789
        PERCENT_VV_TO_G_L = ETHANOL_DENSITY_G_ML * 10.0
        INITIAL_REACTOR_VOLUME_ML = 2000.0
        SMALL_SAMPLE_VOLUME_ML = 5.0
        ETHANOL_SAMPLE_VOLUME_ML = 50.0

        CAMPAIGN_ORDER_TO_PROCESS = {1: "F1", 2: "F2", 3: "F3"}
        PROCESS_TO_CAMPAIGN_ORDER = {v: k for k, v in CAMPAIGN_ORDER_TO_PROCESS.items()}
        INPUT_CHANNEL_COLORS = {
            "N": "tab:green",
            "G": "tab:orange",
            "F": "tab:purple",
            "E": "tab:red",
            "X": "tab:brown",
        }

        display(pd.DataFrame(PROCESS_META).T)
        """
    ),
    md(
        """
        ## Helper functions

        The loader keeps raw identifiers and creates corrected identifiers only in a
        processed layer. This is necessary because Oculyze exports contain a few
        apparent typing/copying issues in sample IDs.
        """
    ),
    code(
        r'''
        def parse_process_datetime(series: pd.Series) -> pd.Series:
            # Parse Oculyze process timestamps from Description.
            clean = (
                series.astype("string")
                .str.strip()
                .str.replace("-", "/", regex=False)
            )
            return pd.to_datetime(clean, dayfirst=True, errors="coerce")


        def parse_analysis_datetime(series: pd.Series) -> pd.Series:
            clean = (
                series.astype("string")
                .str.strip()
                .str.replace("/", "-", regex=False)
            )
            return pd.to_datetime(clean, dayfirst=True, errors="coerce")


        def parse_session_datetime(filename: str):
            patterns = [
                (r"(?:EXPAuto|EXP)\((\d{4}-\d{2}-\d{2}) (\d{2}-\d{2})\)", "%Y-%m-%d %H-%M"),
                (r"ONLINE\((\d{2}-\d{2}-\d{2}) (\d{2}-\d{2})\)", "%y-%m-%d %H-%M"),
            ]
            for pattern, fmt in patterns:
                m = re.search(pattern, filename)
                if m:
                    return pd.to_datetime(f"{m.group(1)} {m.group(2)}", format=fmt, errors="coerce")
            return pd.NaT


        def process_from_folder(path: Path) -> str:
            m = re.search(r"DOE26-(F\d{2})", str(path))
            if not m:
                return None
            return {"F01": "F1", "F02": "F2", "F03": "F3"}.get(m.group(1))


        def process_from_sample_id(sample_id: str):
            m = re.search(r"(DOE-LAB00[123])-\d+", str(sample_id).strip())
            if not m:
                return None
            return PREFIX_TO_PROCESS.get(m.group(1))


        def sample_number_for_process(sample_id: str, process: str):
            prefix = PROCESS_META[process]["sample_prefix"]
            m = re.search(rf"{re.escape(prefix)}-(\d+)", str(sample_id).strip())
            if not m:
                return np.nan
            return int(m.group(1))


        def assign_corrected_oculyze_ids(df: pd.DataFrame) -> pd.DataFrame:
            # Assign chronological IDs while preserving valid non-duplicated raw IDs.
            parts = []
            for process, g in df.sort_values(["process", "sample_datetime", "sample_id_raw"]).groupby("process", sort=True):
                prefix = PROCESS_META[process]["sample_prefix"]
                seen = set()
                corrected_numbers = []
                notes = []
                for _, row in g.iterrows():
                    raw_num = sample_number_for_process(row["sample_id_raw"], process)
                    max_seen = max(seen) if seen else 0
                    if pd.notna(raw_num) and int(raw_num) not in seen and int(raw_num) > max_seen:
                        number = int(raw_num)
                        note = "raw_id_used"
                    else:
                        number = max_seen + 1
                        note = "chronology_inferred"
                    seen.add(number)
                    corrected_numbers.append(number)
                    notes.append(note)
                gg = g.copy()
                gg["sample_number_clean"] = corrected_numbers
                gg["sample_id_clean"] = [f"{prefix}-{n}" for n in corrected_numbers]
                gg["sample_id_resolution"] = notes
                gg["sample_id_corrected"] = gg["sample_id_raw"].astype(str).str.strip() != gg["sample_id_clean"]
                parts.append(gg)
            return pd.concat(parts, ignore_index=True)


        def normalize_analyte(analyte: str):
            a = str(analyte).strip().upper()
            mapping = {
                "ACETIC ACID": "acetic_acid_g_l",
                "AMMONIA": "ammonia_mg_l",
                "PAN": "pan_mg_l",
                "GLYCEROL": "glycerol_g_l",
                "PYRUVIC ACID": "pyruvic_acid_mg_l",
                "GLU-FRU-320": "sugar_total_g_l",
                "GLUCOSE-FRUCTOSE": "sugar_total_g_l",
                "GLUCOSE-320": "glucose_g_l",
                "GLUCOSE": "glucose_g_l",
                "ACETALDEHID-1000": "acetaldehyde_mg_l",
                "ACETALDEHYDE": "acetaldehyde_mg_l",
                "YAN": "yan_reported_mg_l",
            }
            return mapping.get(a)


        def analyte_priority(analyte: str) -> int:
            a = str(analyte).strip().upper()
            priority = {
                "GLU-FRU-320": 0,
                "GLUCOSE-FRUCTOSE": 1,
                "GLUCOSE-320": 0,
                "GLUCOSE": 1,
                "ACETALDEHID-1000": 0,
                "ACETALDEHYDE": 1,
            }
            return priority.get(a, 0)


        def canonical_process_from_y15_sample(sample_id: str):
            m = re.match(r"^(DOE-LAB00[123])-\d+$", str(sample_id).strip())
            if not m:
                return None
            return PREFIX_TO_PROCESS.get(m.group(1))


        def sample_number_from_y15_sample(sample_id: str):
            m = re.match(r"^DOE-LAB00[123]-(\d+)$", str(sample_id).strip())
            if not m:
                return np.nan
            return int(m.group(1))


        def hours_since_t0(timestamp: pd.Series, process: pd.Series, t0_by_process: dict) -> pd.Series:
            t0 = process.map(t0_by_process)
            return (timestamp - t0).dt.total_seconds() / 3600.0


        def add_reactor_volume_from_sampling(df: pd.DataFrame, volume_balance_df: pd.DataFrame) -> pd.DataFrame:
            """Add a stepwise reactor liquid volume based on sampling events."""
            out = df.copy()
            out["reactor_volume_ml"] = INITIAL_REACTOR_VOLUME_ML
            if volume_balance_df is None or volume_balance_df.empty:
                out["volume_source"] = "initial_volume_no_sampling_balance"
                return out

            for process, idx in out.groupby("process").groups.items():
                events = (
                    volume_balance_df[volume_balance_df["process"].eq(process)]
                    .dropna(subset=["t_h", "reactor_volume_after_sample_ml"])
                    .sort_values("t_h")
                )
                if events.empty:
                    out.loc[idx, "volume_source"] = "initial_volume_no_process_events"
                    continue

                event_t = events["t_h"].to_numpy(dtype=float)
                event_v_after = events["reactor_volume_after_sample_ml"].to_numpy(dtype=float)
                t = out.loc[idx, "t_h"].to_numpy(dtype=float)
                pos = np.searchsorted(event_t, t, side="right") - 1
                volume = np.full(len(t), INITIAL_REACTOR_VOLUME_ML, dtype=float)
                mask = pos >= 0
                volume[mask] = event_v_after[pos[mask]]
                out.loc[idx, "reactor_volume_ml"] = volume
                out.loc[idx, "volume_source"] = "sampling_step_balance"
            return out


        def is_pending_marker(value) -> bool:
            if pd.isna(value):
                return False
            text = str(value).strip().lower()
            return text in {"tbd", "pending", "pendiente", "por medir", "por_medir"}


        def is_nonempty_entry(value) -> bool:
            if pd.isna(value):
                return False
            text = str(value).strip()
            return text != "" and text.lower() not in {"nan", "none", "na", "n/a"}


        def to_numeric_manual(value):
            if pd.isna(value) or is_pending_marker(value):
                return np.nan
            text = str(value).strip().replace("%", "").replace(",", ".")
            if text.lower() in {"", "nan", "none", "na", "n/a"}:
                return np.nan
            return pd.to_numeric(text, errors="coerce")


        def parse_optional_datetime(value):
            if pd.isna(value) or not is_nonempty_entry(value):
                return pd.NaT
            return pd.to_datetime(str(value).strip(), dayfirst=True, errors="coerce")


        def write_csv_safe(df: pd.DataFrame, path: Path, index: bool = False) -> Path:
            try:
                df.to_csv(path, index=index)
                return path
            except PermissionError:
                stamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
                fallback = path.with_name(f"{path.stem}.locked_{stamp}{path.suffix}")
                df.to_csv(fallback, index=index)
                warnings.warn(f"Could not overwrite locked CSV {path}. Wrote fallback {fallback}.")
                return fallback


        def resample_time_series(df, value_cols, freq="10min"):
            chunks = []
            for process, g in df.dropna(subset=["timestamp"]).groupby("process"):
                gg = (
                    g.set_index("timestamp")
                    .sort_index()[value_cols]
                    .resample(freq)
                    .median(numeric_only=True)
                    .dropna(how="all")
                    .reset_index()
                )
                gg["process"] = process
                chunks.append(gg)
            return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        '''
    ),
    md("## Load and process Oculyze biomass data"),
    code(
        """
        oculyze_rows = []
        for folder in sorted(DATA_DIR.glob("DOE26-F*")):
            folder_process = process_from_folder(folder)
            for path in sorted(folder.glob("*Oculyze*.csv")):
                df = pd.read_csv(path)
                df["source_file"] = str(path.relative_to(ROOT))
                df["folder_process"] = folder_process
                df["sample_id_raw"] = df["Name"].astype(str).str.strip()
                df["process_from_id"] = df["sample_id_raw"].map(process_from_sample_id)
                df["process"] = df["process_from_id"].fillna(df["folder_process"])
                df["sample_datetime"] = parse_process_datetime(df["Description"])
                for col in ["Concentration", "Viability", "Budding Index", "Density", "Temperature"]:
                    df[col] = pd.to_numeric(df[col], errors="coerce")
                oculyze_rows.append(df)

        oculyze_raw = pd.concat(oculyze_rows, ignore_index=True)
        oculyze = assign_corrected_oculyze_ids(oculyze_raw)
        oculyze["cell_total_million_ml"] = oculyze["Concentration"]
        oculyze["cell_viable_million_ml"] = oculyze["Concentration"] * oculyze["Viability"] / 100.0
        oculyze["cell_dead_million_ml"] = oculyze["Concentration"] - oculyze["cell_viable_million_ml"]
        oculyze["x_total_g_l"] = oculyze["cell_total_million_ml"] * G_L_PER_MILLION_CELLS_ML
        oculyze["x_viable_g_l"] = oculyze["cell_viable_million_ml"] * G_L_PER_MILLION_CELLS_ML
        oculyze["x_dead_g_l"] = oculyze["cell_dead_million_ml"] * G_L_PER_MILLION_CELLS_ML

        t0_by_process = (
            oculyze.dropna(subset=["sample_datetime"])
            .groupby("process")["sample_datetime"]
            .min()
            .to_dict()
        )
        oculyze["t_h"] = hours_since_t0(oculyze["sample_datetime"], oculyze["process"], t0_by_process)

        oculyze_out_cols = [
            "process", "sample_id_raw", "sample_id_clean", "sample_id_corrected",
            "sample_id_resolution", "sample_datetime", "t_h", "Description",
            "cell_total_million_ml", "cell_viable_million_ml", "cell_dead_million_ml",
            "x_total_g_l", "x_viable_g_l", "x_dead_g_l",
            "Viability", "Budding Index", "Density", "Temperature", "source_file",
        ]
        oculyze[oculyze_out_cols].sort_values(["process", "sample_datetime"]).to_csv(
            PROCESSED_DIR / "lot1_oculyze_processed.csv", index=False
        )

        print("t0 by process")
        display(pd.Series(t0_by_process, name="t0").to_frame())
        print("Oculyze processed preview")
        display(oculyze[oculyze_out_cols].sort_values(["process", "sample_datetime"]).head(12))
        """
    ),
    md("## Load and process Y15 chemistry"),
    code(
        """
        y15_rows = []
        y15_dir = DATA_DIR / "Y15_lote1"
        for path in sorted(y15_dir.glob("*.txt")):
            df = pd.read_csv(
                path,
                sep="\\t",
                header=None,
                names=["sample_id", "analyte", "station", "value_raw", "unit", "analysis_datetime_raw"],
                dtype=str,
                encoding="utf-8",
            )
            df["source_file"] = str(path.relative_to(ROOT))
            df["session_datetime"] = parse_session_datetime(path.name)
            y15_rows.append(df)

        y15_raw = pd.concat(y15_rows, ignore_index=True)
        for col in ["sample_id", "analyte", "station", "value_raw", "unit", "analysis_datetime_raw"]:
            y15_raw[col] = y15_raw[col].astype("string").str.strip()

        y15_raw["value"] = pd.to_numeric(
            y15_raw["value_raw"].str.replace(",", ".", regex=False),
            errors="coerce",
        )
        y15_raw["analysis_datetime"] = parse_analysis_datetime(y15_raw["analysis_datetime_raw"])
        y15_raw["source_kind"] = np.select(
            [
                y15_raw["source_file"].str.contains("EXPAuto", regex=False),
                y15_raw["source_file"].str.contains("ONLINE", regex=False),
                y15_raw["source_file"].str.contains("EXP(", regex=False),
            ],
            ["EXPAuto", "ONLINE", "EXP"],
            default="other",
        )

        y15 = y15_raw[y15_raw["sample_id"].str.match(r"^DOE-LAB00[123]-\\d+$", na=False)].copy()
        y15 = y15.drop_duplicates(
            ["sample_id", "analyte", "station", "value_raw", "unit", "analysis_datetime_raw"]
        )
        y15["process"] = y15["sample_id"].map(canonical_process_from_y15_sample)
        y15["analyte_canonical"] = y15["analyte"].map(normalize_analyte)
        y15["analyte_priority"] = y15["analyte"].map(analyte_priority)
        y15["is_numeric"] = y15["value"].notna()

        y15_relevant = y15.dropna(subset=["analyte_canonical"]).copy()
        y15_relevant = y15_relevant.sort_values(
            ["sample_id", "analyte_canonical", "is_numeric", "analyte_priority", "analysis_datetime"],
            ascending=[True, True, False, True, True],
        )
        y15_selected = (
            y15_relevant.groupby(["sample_id", "analyte_canonical"], as_index=False)
            .first()
        )

        y15_wide = y15_selected.pivot(
            index="sample_id",
            columns="analyte_canonical",
            values="value",
        ).reset_index()
        y15_wide.columns.name = None
        y15_wide["process"] = y15_wide["sample_id"].map(canonical_process_from_y15_sample)

        sample_time_map = (
            oculyze.sort_values("sample_datetime")
            .drop_duplicates("sample_id_clean")
            .set_index("sample_id_clean")["sample_datetime"]
        )
        sample_number_t_h_reference = (
            oculyze.dropna(subset=["sample_number_clean", "t_h"])
            .groupby("sample_number_clean")["t_h"]
            .median()
        )
        session_time_map = y15_selected.groupby("sample_id")["session_datetime"].min()
        y15_wide["sample_number"] = y15_wide["sample_id"].map(sample_number_from_y15_sample)
        y15_wide["sample_datetime_oculyze"] = y15_wide["sample_id"].map(sample_time_map)
        y15_wide["t_h_schedule_reference"] = y15_wide["sample_number"].map(sample_number_t_h_reference)
        y15_wide["sample_datetime_schedule"] = [
            (
                t0_by_process[row.process] + pd.to_timedelta(row.t_h_schedule_reference, unit="h")
                if pd.notna(row.t_h_schedule_reference) and row.process in t0_by_process
                else pd.NaT
            )
            for row in y15_wide.itertuples(index=False)
        ]
        y15_wide["sample_datetime_session"] = y15_wide["sample_id"].map(session_time_map)
        y15_wide["sample_datetime"] = y15_wide["sample_datetime_oculyze"].combine_first(
            y15_wide["sample_datetime_schedule"]
        ).combine_first(
            y15_wide["sample_datetime_session"]
        )
        y15_wide["time_source"] = np.select(
            [
                y15_wide["sample_datetime_oculyze"].notna(),
                y15_wide["sample_datetime_schedule"].notna(),
                y15_wide["sample_datetime_session"].notna(),
            ],
            ["oculyze_sample_id", "inferred_lot_sample_sequence", "y15_session_filename"],
            default="missing_time",
        )
        y15_wide["t_h"] = hours_since_t0(y15_wide["sample_datetime"], y15_wide["process"], t0_by_process)

        if {"sugar_total_g_l", "glucose_g_l"}.issubset(y15_wide.columns):
            y15_wide["fructose_g_l"] = y15_wide["sugar_total_g_l"] - y15_wide["glucose_g_l"]
        if {"ammonia_mg_l", "pan_mg_l"}.issubset(y15_wide.columns):
            y15_wide["yan_calc_mg_l"] = y15_wide["ammonia_mg_l"] + y15_wide["pan_mg_l"]

        for col in ["ammonia_mg_l", "pan_mg_l", "yan_calc_mg_l"]:
            if col in y15_wide.columns:
                y15_wide[f"{col}_raw"] = y15_wide[col]
                y15_wide[f"{col}_below_zero"] = y15_wide[col] < 0
                y15_wide[f"{col}_censored"] = y15_wide[col].clip(lower=0)

        if {"ammonia_mg_l_censored", "pan_mg_l_censored"}.issubset(y15_wide.columns):
            y15_wide["yan_calc_mg_l_censored"] = (
                y15_wide["ammonia_mg_l_censored"] + y15_wide["pan_mg_l_censored"]
            )
            y15_wide["yan_for_model_mg_l"] = y15_wide["yan_calc_mg_l_censored"]
            y15_wide["nitrogen_any_negative_raw"] = (
                y15_wide[["ammonia_mg_l_raw", "pan_mg_l_raw"]].lt(0).any(axis=1)
            )
            y15_wide["nitrogen_severe_negative_raw"] = (
                y15_wide[["ammonia_mg_l_raw", "pan_mg_l_raw"]].lt(-20).any(axis=1)
            )
            y15_wide["nitrogen_qc_class"] = np.select(
                [
                    y15_wide["nitrogen_severe_negative_raw"],
                    y15_wide["nitrogen_any_negative_raw"],
                ],
                ["severe_negative_raw_censored_to_zero", "negative_raw_censored_to_zero"],
                default="ok_raw_nonnegative",
            )

        y15.to_csv(PROCESSED_DIR / "lot1_y15_long_processed.csv", index=False)
        y15_wide.sort_values(["process", "sample_datetime", "sample_id"]).to_csv(
            PROCESSED_DIR / "lot1_y15_wide_processed.csv", index=False
        )

        print("Y15 relevant long rows:", len(y15_relevant))
        print("Y15 selected sample-analyte rows:", len(y15_selected))
        print("Y15 wide preview")
        display(y15_wide.sort_values(["process", "sample_datetime", "sample_id"]).head(12))
        """
    ),
    md(
        """
        ## Manual ethanol entries and sample-volume balance

        Ethanol measurements are entered manually in
        `DOE_Lote_1_ethanol_manual_entry.xlsx`. Only the experiment-specific sheets
        (`F1`, `F2`, `F3`) are used here, because those are the sheets updated during
        laboratory work.

        Rules used in this preview:

        - Numeric `% v/v` entries are converted as `ethanol_g_l = ethanol_%vv * 7.89`.
        - `TBD` is treated as a pending ethanol measurement.
        - If a sample has numeric ethanol or `TBD`, it is counted as a 50 mL sample.
        - Otherwise, it is counted as a 5 mL small sample.
        - The volume balance accounts only for sampling losses, not for stock
          additions or evaporation.
        """
    ),
    code(
        """
        def load_manual_ethanol_template(path: Path) -> pd.DataFrame:
            if not path.exists():
                warnings.warn(f"Manual ethanol workbook not found: {path}")
                return pd.DataFrame()

            frames = []
            try:
                for process in PROCESS_META:
                    df = pd.read_excel(path, sheet_name=process, dtype=object)
                    df["source_sheet"] = process
                    frames.append(df)
                ethanol = pd.concat(frames, ignore_index=True)
            except PermissionError as exc:
                fallback = PROCESSED_DIR / "lot1_ethanol_manual_processed.csv"
                if fallback.exists():
                    warnings.warn(
                        f"Could not read open/locked ethanol workbook ({path}). "
                        f"Using last processed CSV fallback: {fallback}"
                    )
                    return pd.read_csv(fallback, dtype=object)
                raise exc

            for col in ethanol.columns:
                if ethanol[col].dtype == object:
                    ethanol[col] = ethanol[col].map(lambda x: x.strip() if isinstance(x, str) else x)

            input_cols = [
                "ethanol_percent_vv_rep1",
                "ethanol_percent_vv_rep2",
                "ethanol_g_l_manual",
            ]
            for col in input_cols:
                if col not in ethanol.columns:
                    ethanol[col] = np.nan

            ethanol["ethanol_percent_vv_rep1_num"] = ethanol["ethanol_percent_vv_rep1"].map(to_numeric_manual)
            ethanol["ethanol_percent_vv_rep2_num"] = ethanol["ethanol_percent_vv_rep2"].map(to_numeric_manual)
            ethanol["ethanol_g_l_manual_num"] = ethanol["ethanol_g_l_manual"].map(to_numeric_manual)
            ethanol["ethanol_pending"] = ethanol[input_cols].applymap(is_pending_marker).any(axis=1)
            ethanol["ethanol_has_entry"] = ethanol[input_cols].applymap(is_nonempty_entry).any(axis=1)

            ethanol["ethanol_percent_vv_final_calc"] = ethanol[
                ["ethanol_percent_vv_rep1_num", "ethanol_percent_vv_rep2_num"]
            ].mean(axis=1, skipna=True)
            ethanol.loc[
                ethanol[["ethanol_percent_vv_rep1_num", "ethanol_percent_vv_rep2_num"]].notna().sum(axis=1).eq(0),
                "ethanol_percent_vv_final_calc",
            ] = np.nan
            ethanol["ethanol_g_l_from_percent_vv"] = (
                ethanol["ethanol_percent_vv_final_calc"] * PERCENT_VV_TO_G_L
            )
            ethanol["ethanol_g_l_for_model_manual"] = ethanol["ethanol_g_l_manual_num"].combine_first(
                ethanol["ethanol_g_l_from_percent_vv"]
            )
            ethanol["ethanol_status"] = np.select(
                [
                    ethanol["ethanol_g_l_for_model_manual"].notna(),
                    ethanol["ethanol_pending"],
                    ethanol["ethanol_has_entry"],
                ],
                ["measured", "pending_TBD", "non_numeric_entry"],
                default="not_measured",
            )
            ethanol["ethanol_sample_type"] = np.where(
                ethanol["ethanol_status"].isin(["measured", "pending_TBD", "non_numeric_entry"]),
                "ethanol_50mL",
                "small_5mL",
            )
            ethanol["sample_volume_ml"] = np.where(
                ethanol["ethanol_sample_type"].eq("ethanol_50mL"),
                ETHANOL_SAMPLE_VOLUME_ML,
                SMALL_SAMPLE_VOLUME_ML,
            )
            keep_cols = [
                "process",
                "sample_id",
                "source_sheet",
                "ethanol_percent_vv_rep1",
                "ethanol_percent_vv_rep2",
                "ethanol_percent_vv_rep1_num",
                "ethanol_percent_vv_rep2_num",
                "ethanol_percent_vv_final_calc",
                "ethanol_g_l_manual_num",
                "ethanol_g_l_from_percent_vv",
                "ethanol_g_l_for_model_manual",
                "ethanol_pending",
                "ethanol_status",
                "ethanol_sample_type",
                "sample_volume_ml",
                "ethanol_measurement_datetime",
                "ethanol_method",
                "operator",
                "qc_flag",
                "notes",
            ]
            return ethanol[[c for c in keep_cols if c in ethanol.columns]].copy()


        ethanol_manual = load_manual_ethanol_template(ETHANOL_TEMPLATE_XLSX)
        if not ethanol_manual.empty:
            for col in [
                "ethanol_percent_vv_rep1_num",
                "ethanol_percent_vv_rep2_num",
                "ethanol_percent_vv_final_calc",
                "ethanol_g_l_manual_num",
                "ethanol_g_l_from_percent_vv",
                "ethanol_g_l_for_model_manual",
                "sample_volume_ml",
            ]:
                if col in ethanol_manual.columns:
                    ethanol_manual[col] = pd.to_numeric(ethanol_manual[col], errors="coerce")
            if "sample_volume_ml" in ethanol_manual.columns:
                ethanol_manual["sample_volume_ml"] = ethanol_manual["sample_volume_ml"].fillna(SMALL_SAMPLE_VOLUME_ML)
            ethanol_manual.to_csv(PROCESSED_DIR / "lot1_ethanol_manual_processed.csv", index=False)

            merge_cols = [
                "process",
                "sample_id",
                "ethanol_percent_vv_rep1",
                "ethanol_percent_vv_rep2",
                "ethanol_percent_vv_rep1_num",
                "ethanol_percent_vv_rep2_num",
                "ethanol_percent_vv_final_calc",
                "ethanol_g_l_manual_num",
                "ethanol_g_l_from_percent_vv",
                "ethanol_g_l_for_model_manual",
                "ethanol_pending",
                "ethanol_status",
                "ethanol_sample_type",
                "sample_volume_ml",
                "ethanol_measurement_datetime",
                "ethanol_method",
                "operator",
                "qc_flag",
                "notes",
            ]
            y15_wide = y15_wide.merge(
                ethanol_manual[merge_cols],
                on=["process", "sample_id"],
                how="left",
                suffixes=("", "_ethanol_template"),
            )
            y15_wide["ethanol_status"] = y15_wide["ethanol_status"].fillna("not_measured")
            y15_wide["ethanol_sample_type"] = y15_wide["ethanol_sample_type"].fillna("small_5mL")
            y15_wide["sample_volume_ml"] = y15_wide["sample_volume_ml"].fillna(SMALL_SAMPLE_VOLUME_ML)

            balance_parts = []
            for process, g in y15_wide.sort_values(["process", "sample_datetime", "sample_id"]).groupby("process"):
                gg = g.copy()
                gg["reactor_volume_before_sample_ml"] = INITIAL_REACTOR_VOLUME_ML - gg["sample_volume_ml"].cumsum().shift(fill_value=0)
                gg["reactor_volume_after_sample_ml"] = gg["reactor_volume_before_sample_ml"] - gg["sample_volume_ml"]
                gg["cumulative_sample_volume_removed_ml"] = INITIAL_REACTOR_VOLUME_ML - gg["reactor_volume_after_sample_ml"]
                balance_parts.append(gg)
            y15_wide = pd.concat(balance_parts, ignore_index=True)

            volume_cols = [
                "process",
                "sample_id",
                "sample_datetime",
                "t_h",
                "ethanol_status",
                "ethanol_sample_type",
                "sample_volume_ml",
                "reactor_volume_before_sample_ml",
                "reactor_volume_after_sample_ml",
                "cumulative_sample_volume_removed_ml",
                "ethanol_percent_vv_rep1",
                "ethanol_percent_vv_rep2",
                "ethanol_percent_vv_final_calc",
                "ethanol_g_l_for_model_manual",
            ]
            volume_balance = y15_wide[volume_cols].sort_values(["process", "sample_datetime", "sample_id"]).copy()
            volume_balance.to_csv(PROCESSED_DIR / "lot1_sample_volume_balance.csv", index=False)
            y15_wide.sort_values(["process", "sample_datetime", "sample_id"]).to_csv(
                PROCESSED_DIR / "lot1_y15_wide_processed.csv", index=False
            )

            print("Manual ethanol workbook")
            print(ETHANOL_TEMPLATE_XLSX)
            print("Ethanol status by process")
            display(pd.crosstab(y15_wide["process"], y15_wide["ethanol_status"]))
            print("Final sampling volume balance by process")
            display(
                volume_balance.groupby("process").tail(1)[
                    [
                        "process",
                        "sample_id",
                        "t_h",
                        "cumulative_sample_volume_removed_ml",
                        "reactor_volume_after_sample_ml",
                    ]
                ]
            )
            print("Volume balance preview")
            display(volume_balance.head(12))
        else:
            volume_balance = pd.DataFrame()
            print("No manual ethanol workbook loaded.")
        """
    ),
    md(
        """
        ## Planned and executed operational input events

        The design bundle contains planned manual pulses for Lot 1. The online
        temperature files also include `nutricion_activa`, but in the current Lot 1
        files that signal remains zero, so it does not provide an executed-pulse log.

        This section therefore creates a separate editable execution table:

        - Planned pulses are read from `db_20260609/rec/lot_1_schedule.csv`.
        - `planned_design_t_h` is the model-based design time relative to process
          start.
        - `planned_clock_t_h_from_actual_t0` is the old scheduled clock time
          re-expressed relative to the actual `t = 0` from Oculyze.
        - If real execution differed, fill `actual_datetime` or `actual_t_h`,
          `actual_amount`, and `actual_stock_volume_ml` in
          `DOE_Lote_1_operational_inputs_actual.csv`.

        The plots show planned inputs as dashed vertical lines and actual executed
        inputs as solid vertical lines when the editable table contains actual
        times.
        """
    ),
    code(
        """
        def amount_unit_for_channel(channel: str) -> str:
            return {
                "N": "kg/m3 as N, equivalent to g/L or 1000 mg/L",
                "G": "g/L glucose equivalent",
                "F": "g/L fructose equivalent",
                "E": "g/L ethanol",
                "X": "g/L dry biomass",
            }.get(str(channel), "model units")


        def load_planned_input_events() -> pd.DataFrame:
            if not LOT1_PLANNED_SCHEDULE_CSV.exists():
                warnings.warn(f"Planned schedule not found: {LOT1_PLANNED_SCHEDULE_CSV}")
                return pd.DataFrame()

            schedule = pd.read_csv(LOT1_PLANNED_SCHEDULE_CSV)
            pulses = schedule[schedule["event"].eq("manual_pulse")].copy()
            if pulses.empty:
                return pulses

            pulses["process"] = pulses["campaign_order"].map(CAMPAIGN_ORDER_TO_PROCESS)
            pulses["planned_datetime"] = pd.to_datetime(pulses["datetime"], errors="coerce")
            pulses["planned_design_t_h"] = pd.to_numeric(pulses["relative_time_h"], errors="coerce")
            pulses["planned_clock_t_h_from_actual_t0"] = hours_since_t0(
                pulses["planned_datetime"], pulses["process"], t0_by_process
            )
            pulses["planned_amount"] = pd.to_numeric(pulses["amount"], errors="coerce")
            pulses["planned_stock_volume_ml"] = pd.to_numeric(pulses["volume_ml"], errors="coerce")
            pulses["amount_unit"] = pulses["channel"].map(amount_unit_for_channel)
            pulses["event_id"] = [
                f"{row.process}_{row.channel}_{row.planned_design_t_h:g}h"
                for row in pulses.itertuples(index=False)
            ]
            pulses["event_source"] = "planned_design"
            keep = [
                "event_id",
                "process",
                "campaign_order",
                "candidate",
                "channel",
                "planned_design_t_h",
                "planned_datetime",
                "planned_clock_t_h_from_actual_t0",
                "planned_amount",
                "amount_unit",
                "planned_stock_volume_ml",
                "details",
                "event_source",
            ]
            return pulses[keep].sort_values(["process", "planned_design_t_h", "channel"]).reset_index(drop=True)


        def ensure_actual_input_template(planned: pd.DataFrame) -> pd.DataFrame:
            ACTUAL_INPUT_EVENTS_CSV.parent.mkdir(parents=True, exist_ok=True)
            if ACTUAL_INPUT_EVENTS_CSV.exists():
                actual = pd.read_csv(ACTUAL_INPUT_EVENTS_CSV, dtype=object)
            else:
                actual = planned.copy()
                actual["executed"] = ""
                actual["actual_datetime"] = ""
                actual["actual_t_h"] = ""
                actual["actual_amount"] = ""
                actual["actual_amount_unit"] = actual["amount_unit"]
                actual["actual_stock_volume_ml"] = ""
                actual["operator"] = ""
                actual["qc_flag"] = ""
                actual["notes"] = ""
                actual.to_csv(ACTUAL_INPUT_EVENTS_CSV, index=False)
            return actual


        def parse_executed_flag(value) -> bool:
            if pd.isna(value):
                return False
            return str(value).strip().lower() in {"true", "1", "yes", "y", "si", "sí", "executed", "ok"}


        planned_input_events = load_planned_input_events()
        if not planned_input_events.empty:
            planned_input_events.to_csv(PROCESSED_DIR / "lot1_planned_operational_inputs.csv", index=False)
            actual_input_events = ensure_actual_input_template(planned_input_events)

            required_actual_cols = {
                "executed": "",
                "actual_datetime": "",
                "actual_t_h": "",
                "actual_amount": "",
                "actual_amount_unit": "",
                "actual_stock_volume_ml": "",
                "operator": "",
                "qc_flag": "",
                "notes": "",
            }
            for col, default in required_actual_cols.items():
                if col not in actual_input_events.columns:
                    actual_input_events[col] = default

            for col in ["planned_design_t_h", "planned_clock_t_h_from_actual_t0", "planned_amount", "planned_stock_volume_ml"]:
                if col in actual_input_events.columns:
                    actual_input_events[col] = pd.to_numeric(actual_input_events[col], errors="coerce")
            for col in ["actual_t_h", "actual_amount", "actual_stock_volume_ml"]:
                if col in actual_input_events.columns:
                    actual_input_events[col] = actual_input_events[col].map(to_numeric_manual)

            actual_input_events["process"] = actual_input_events["process"].astype(str).str.strip()
            actual_input_events["channel"] = actual_input_events["channel"].astype(str).str.strip()
            actual_input_events["actual_datetime_parsed"] = actual_input_events["actual_datetime"].map(parse_optional_datetime)
            actual_input_events["actual_t_h_from_datetime"] = hours_since_t0(
                actual_input_events["actual_datetime_parsed"],
                actual_input_events["process"],
                t0_by_process,
            )
            actual_input_events["actual_plot_t_h"] = actual_input_events["actual_t_h"].combine_first(
                actual_input_events["actual_t_h_from_datetime"]
            )
            actual_input_events["executed_bool"] = actual_input_events["executed"].map(parse_executed_flag)
            actual_input_events["has_actual_time"] = actual_input_events["actual_plot_t_h"].notna()
            actual_input_events["actual_amount_for_plot"] = actual_input_events["actual_amount"].combine_first(
                actual_input_events["planned_amount"]
            )
            actual_input_events["actual_stock_volume_ml_for_plot"] = actual_input_events["actual_stock_volume_ml"].combine_first(
                actual_input_events["planned_stock_volume_ml"]
            )
            actual_input_events["actual_minus_planned_amount"] = (
                actual_input_events["actual_amount"] - actual_input_events["planned_amount"]
            )
            actual_input_events["actual_over_planned_amount"] = (
                actual_input_events["actual_amount"] / actual_input_events["planned_amount"]
            )
            actual_input_events["deviation_vs_doe_nominal_h"] = (
                actual_input_events["actual_plot_t_h"] - actual_input_events["planned_design_t_h"]
            )
            actual_input_events["deviation_vs_original_calendar_h"] = (
                actual_input_events["actual_plot_t_h"]
                - actual_input_events["planned_clock_t_h_from_actual_t0"]
            )
            actual_input_events["actual_time_consistency_h"] = (
                actual_input_events["actual_t_h"] - actual_input_events["actual_t_h_from_datetime"]
            )
            actual_input_events.to_csv(PROCESSED_DIR / "lot1_actual_operational_inputs_loaded.csv", index=False)

            comparison_cols = [
                "event_id",
                "process",
                "candidate",
                "channel",
                "planned_design_t_h",
                "planned_clock_t_h_from_actual_t0",
                "actual_plot_t_h",
                "deviation_vs_doe_nominal_h",
                "deviation_vs_original_calendar_h",
                "planned_amount",
                "actual_amount",
                "actual_minus_planned_amount",
                "actual_over_planned_amount",
                "planned_stock_volume_ml",
                "actual_stock_volume_ml",
                "actual_datetime",
                "actual_time_consistency_h",
                "qc_flag",
                "notes",
            ]
            input_comparison = actual_input_events[
                [c for c in comparison_cols if c in actual_input_events.columns]
            ].copy()
            input_comparison.to_csv(PROCESSED_DIR / "lot1_operational_input_comparison.csv", index=False)

            print("Planned operational inputs")
            display(planned_input_events)
            print("Editable actual operational input log")
            print(ACTUAL_INPUT_EVENTS_CSV)
            display(actual_input_events)
            print("Operational input comparison")
            display(input_comparison)
        else:
            actual_input_events = pd.DataFrame()
            print("No planned manual input pulses found for Lot 1.")
        """
    ),
    md(
        """
        ## Nitrogen physical-consistency QC

        YAN is not expected to increase unless there is a nitrogen input. Therefore,
        after the measured/censored YAN trajectory reaches a near-zero level, any
        later positive rebound before the next nitrogen pulse is treated as an
        analytical artifact unless there is independent evidence of an unlogged N
        addition.

        This is especially relevant for F2: the trajectory reaches near-zero around
        the glucose pulse / 50 h window, but later Y15 points show a positive
        PAN/YAN rebound before the actual N pulse at 72.5 h. Based on the lab note,
        F2 nitrogen points with 50 < t < t_N,pulse are treated as analytical
        artifacts. Those points are retained as raw observations and flagged, but
        the physical-QC model series treats them as zero/censored.
        """
    ),
    code(
        """
        def add_nitrogen_physical_consistency_qc(y15_df: pd.DataFrame, input_events: pd.DataFrame) -> pd.DataFrame:
            out = y15_df.copy()
            if "yan_for_model_mg_l" not in out.columns:
                return out

            out["yan_for_model_physical_mg_l"] = out["yan_for_model_mg_l"]
            out["nitrogen_physical_qc_flag"] = "ok"
            out["nitrogen_positive_rebound_before_N_pulse"] = False
            out["nitrogen_manual_zero_between_50h_and_N_pulse"] = False
            out["nitrogen_next_N_pulse_t_h"] = np.nan
            out["nitrogen_previous_N_pulse_t_h"] = np.nan

            if input_events is None or input_events.empty:
                n_events = pd.DataFrame()
            else:
                n_events = input_events[
                    input_events["channel"].astype(str).str.upper().eq("N")
                    & input_events["actual_plot_t_h"].notna()
                ].copy()
                if n_events.empty and "planned_design_t_h" in input_events.columns:
                    n_events = input_events[
                        input_events["channel"].astype(str).str.upper().eq("N")
                    ].copy()
                    n_events["actual_plot_t_h"] = n_events["planned_design_t_h"]

            depletion_threshold_mg_l = 2.0
            rebound_threshold_mg_l = 5.0
            eps = 1e-9

            for process, idx in out.groupby("process").groups.items():
                g = out.loc[idx].sort_values("t_h").copy()
                pulse_times = sorted(
                    n_events.loc[n_events["process"].eq(process), "actual_plot_t_h"]
                    .dropna()
                    .astype(float)
                    .tolist()
                )

                prev_pulse_for_row = []
                next_pulse_for_row = []
                for t in g["t_h"].astype(float):
                    prev_pulses = [p for p in pulse_times if p <= t + eps]
                    next_pulses = [p for p in pulse_times if p > t + eps]
                    prev_pulse_for_row.append(max(prev_pulses) if prev_pulses else np.nan)
                    next_pulse_for_row.append(min(next_pulses) if next_pulses else np.nan)

                g["nitrogen_previous_N_pulse_t_h"] = prev_pulse_for_row
                g["nitrogen_next_N_pulse_t_h"] = next_pulse_for_row

                # Work in intervals before each N pulse and after each pulse.
                interval_starts = [-np.inf] + pulse_times
                interval_ends = pulse_times + [np.inf]
                for start, end in zip(interval_starts, interval_ends):
                    interval = g[(g["t_h"] >= start - eps) & (g["t_h"] < end - eps)].copy()
                    if interval.empty:
                        continue
                    running_min = np.inf
                    depleted_seen = False
                    for row in interval.itertuples():
                        value = getattr(row, "yan_for_model_mg_l")
                        if pd.isna(value):
                            continue
                        if running_min <= depletion_threshold_mg_l:
                            depleted_seen = True
                        is_before_next_pulse = np.isfinite(end) and row.t_h < end - eps
                        is_rebound = (
                            depleted_seen
                            and is_before_next_pulse
                            and value > rebound_threshold_mg_l
                        )
                        if is_rebound:
                            g.loc[row.Index, "nitrogen_positive_rebound_before_N_pulse"] = True
                            g.loc[row.Index, "nitrogen_physical_qc_flag"] = (
                                "positive_rebound_after_depletion_before_N_pulse_set_to_zero"
                            )
                            g.loc[row.Index, "yan_for_model_physical_mg_l"] = 0.0
                        running_min = min(running_min, float(value))

                if process == "F2":
                    # Lab note: keep values through the near-50 h point, but treat
                    # positive YAN/PAN between that point and the real N pulse as a
                    # method artifact. Raw/censored columns remain unchanged.
                    f2_next_n_pulses = [p for p in pulse_times if p > 50.0 + eps]
                    f2_n_pulse_t = min(f2_next_n_pulses) if f2_next_n_pulses else 72.5
                    manual_mask = (
                        g["t_h"].astype(float).gt(50.0 + eps)
                        & g["t_h"].astype(float).lt(f2_n_pulse_t - eps)
                        & g["yan_for_model_mg_l"].notna()
                    )
                    if manual_mask.any():
                        g.loc[manual_mask, "nitrogen_positive_rebound_before_N_pulse"] = True
                        g.loc[manual_mask, "nitrogen_manual_zero_between_50h_and_N_pulse"] = True
                        g.loc[manual_mask, "nitrogen_physical_qc_flag"] = (
                            "manual_rule_F2_between_50h_and_N_pulse_set_to_zero"
                        )
                        g.loc[manual_mask, "yan_for_model_physical_mg_l"] = 0.0

                out.loc[g.index, [
                    "yan_for_model_physical_mg_l",
                    "nitrogen_physical_qc_flag",
                    "nitrogen_positive_rebound_before_N_pulse",
                    "nitrogen_manual_zero_between_50h_and_N_pulse",
                    "nitrogen_next_N_pulse_t_h",
                    "nitrogen_previous_N_pulse_t_h",
                ]] = g[[
                    "yan_for_model_physical_mg_l",
                    "nitrogen_physical_qc_flag",
                    "nitrogen_positive_rebound_before_N_pulse",
                    "nitrogen_manual_zero_between_50h_and_N_pulse",
                    "nitrogen_next_N_pulse_t_h",
                    "nitrogen_previous_N_pulse_t_h",
                ]]
            return out


        y15_wide = add_nitrogen_physical_consistency_qc(y15_wide, actual_input_events)

        if "nitrogen_qc_class" in y15_wide.columns:
            nitrogen_qc_cols = [
                "process",
                "sample_id",
                "sample_datetime",
                "t_h",
                "ammonia_mg_l_raw",
                "pan_mg_l_raw",
                "yan_calc_mg_l_raw",
                "ammonia_mg_l_censored",
                "pan_mg_l_censored",
                "yan_for_model_mg_l",
                "yan_for_model_physical_mg_l",
                "nitrogen_qc_class",
                "nitrogen_physical_qc_flag",
                "nitrogen_positive_rebound_before_N_pulse",
                "nitrogen_manual_zero_between_50h_and_N_pulse",
                "nitrogen_previous_N_pulse_t_h",
                "nitrogen_next_N_pulse_t_h",
                "nitrogen_has_conflicting_y15_export",
            ]
            nitrogen_qc = y15_wide[
                [c for c in nitrogen_qc_cols if c in y15_wide.columns]
            ].sort_values(["process", "sample_datetime", "sample_id"])
            write_csv_safe(nitrogen_qc, PROCESSED_DIR / "lot1_nitrogen_qc_processed.csv", index=False)
            write_csv_safe(
                y15_wide.sort_values(["process", "sample_datetime", "sample_id"]),
                PROCESSED_DIR / "lot1_y15_wide_processed.csv",
                index=False,
            )

            print("Nitrogen physical-consistency QC")
            display(nitrogen_qc)
        """
    ),
    md(
        """
        ## Load online CO2 and temperature

        The CO2 sensor reports total gas flow from each reactor. After sampling, the
        liquid volume is smaller, so the same volumetric fermentation rate would
        generate a smaller total CO2 flow. The notebook therefore adds:

        - `reactor_volume_ml`: stepwise liquid volume after each sampling event.
        - `flow_filt_sccm_per_l`: raw filtered CO2 divided by current liquid volume.
        - `flow_filt_sccm_v0_corrected`: raw filtered CO2 scaled to the initial
          2 L volume:

        $$F_{CO2,2L}(t) = F_{CO2,meas}(t) \\frac{V_0}{V(t)}$$

        Use `flow_filt_sccm_per_l` when the model predicts a volumetric CO2 rate.
        Use `flow_filt_sccm_v0_corrected` when comparing to a total-flow prediction
        normalized to the initial reactor volume.
        """
    ),
    code(
        """
        co2_rows = []
        temp_rows = []
        for process, meta in PROCESS_META.items():
            folder = next(DATA_DIR.glob(meta["folder_glob"]))

            co2_path = next(folder.glob("CO2_FILT_*.csv"))
            co2 = pd.read_csv(co2_path)
            co2["source_file"] = str(co2_path.relative_to(ROOT))
            co2["process"] = process
            co2["timestamp"] = pd.to_datetime(co2["timestamp"], errors="coerce")
            co2["flow_filt_sccm"] = pd.to_numeric(co2["flow_filt_sccm"], errors="coerce")
            co2["t_h"] = hours_since_t0(co2["timestamp"], co2["process"], t0_by_process)
            co2_rows.append(co2)

            temp_path = next(folder.glob("Temp_*.csv"))
            temp = pd.read_csv(temp_path)
            temp["source_file"] = str(temp_path.relative_to(ROOT))
            temp["process"] = process
            temp["timestamp"] = pd.to_datetime(temp["timestamp"], errors="coerce")
            for col in ["T", "SP", "banda", "cold", "hot", "nutricion_activa", "freq_nut"]:
                temp[col] = pd.to_numeric(temp[col], errors="coerce")
            temp["t_h"] = hours_since_t0(temp["timestamp"], temp["process"], t0_by_process)
            temp_rows.append(temp)

        co2_raw = pd.concat(co2_rows, ignore_index=True)
        temp_raw = pd.concat(temp_rows, ignore_index=True)

        if "volume_balance" not in globals():
            volume_balance = pd.DataFrame()

        co2_raw = add_reactor_volume_from_sampling(co2_raw, volume_balance)
        co2_raw["reactor_volume_l"] = co2_raw["reactor_volume_ml"] / 1000.0
        co2_raw["co2_volume_correction_factor"] = INITIAL_REACTOR_VOLUME_ML / co2_raw["reactor_volume_ml"]
        co2_raw["flow_filt_sccm_per_l"] = co2_raw["flow_filt_sccm"] / co2_raw["reactor_volume_l"]
        co2_raw["flow_filt_sccm_v0_corrected"] = (
            co2_raw["flow_filt_sccm"] * co2_raw["co2_volume_correction_factor"]
        )
        co2_raw.loc[co2_raw["t_h"] < 0, [
            "reactor_volume_ml",
            "reactor_volume_l",
            "co2_volume_correction_factor",
            "flow_filt_sccm_per_l",
            "flow_filt_sccm_v0_corrected",
        ]] = np.nan

        co2_value_cols = [
            "flow_filt_sccm",
            "reactor_volume_ml",
            "reactor_volume_l",
            "co2_volume_correction_factor",
            "flow_filt_sccm_per_l",
            "flow_filt_sccm_v0_corrected",
        ]
        co2_10min = resample_time_series(co2_raw[co2_raw["t_h"] >= 0], co2_value_cols, freq="10min")
        co2_10min["t_h"] = hours_since_t0(co2_10min["timestamp"], co2_10min["process"], t0_by_process)

        temp_10min = resample_time_series(temp_raw[temp_raw["t_h"] >= 0], ["T", "SP"], freq="10min")
        temp_10min["t_h"] = hours_since_t0(temp_10min["timestamp"], temp_10min["process"], t0_by_process)

        co2_raw[co2_raw["t_h"] >= 0].to_csv(PROCESSED_DIR / "lot1_co2_filt_volume_corrected_raw.csv", index=False)
        co2_10min.to_csv(PROCESSED_DIR / "lot1_co2_filt_10min_preview.csv", index=False)
        co2_10min.to_csv(PROCESSED_DIR / "lot1_co2_filt_volume_corrected_10min.csv", index=False)
        temp_10min.to_csv(PROCESSED_DIR / "lot1_temp_10min_preview.csv", index=False)

        display(pd.DataFrame({
            "co2_raw_rows": co2_raw.groupby("process").size(),
            "co2_10min_rows": co2_10min.groupby("process").size(),
            "temp_raw_rows": temp_raw.groupby("process").size(),
            "temp_10min_rows": temp_10min.groupby("process").size(),
        }))

        print("CO2 volume correction summary")
        display(
            co2_10min.groupby("process")
            .agg(
                min_volume_ml=("reactor_volume_ml", "min"),
                max_correction_factor=("co2_volume_correction_factor", "max"),
                max_raw_sccm=("flow_filt_sccm", "max"),
                max_v0_corrected_sccm=("flow_filt_sccm_v0_corrected", "max"),
                max_sccm_per_l=("flow_filt_sccm_per_l", "max"),
            )
            .reset_index()
        )
        """
    ),
    md("## Inventory and quality-control checks"),
    code(
        """
        inventory = []
        for process, meta in PROCESS_META.items():
            inventory.append({
                "process": process,
                "design_name": meta["design_name"],
                "t0": t0_by_process.get(process),
                "oculyze_samples": int((oculyze["process"] == process).sum()),
                "y15_samples": int((y15_wide["process"] == process).sum()),
                "co2_rows_t_ge_0": int(((co2_raw["process"] == process) & (co2_raw["t_h"] >= 0)).sum()),
                "temp_rows_t_ge_0": int(((temp_raw["process"] == process) & (temp_raw["t_h"] >= 0)).sum()),
                "last_oculyze_h": oculyze.loc[oculyze["process"] == process, "t_h"].max(),
                "last_y15_h": y15_wide.loc[y15_wide["process"] == process, "t_h"].max(),
                "last_co2_h": co2_raw.loc[(co2_raw["process"] == process) & (co2_raw["t_h"] >= 0), "t_h"].max(),
            })
        inventory = pd.DataFrame(inventory)
        display(inventory)

        corrected_ids = oculyze.loc[
            oculyze["sample_id_corrected"],
            ["process", "sample_datetime", "sample_id_raw", "sample_id_clean", "sample_id_resolution", "source_file"],
        ].sort_values(["process", "sample_datetime"])

        y15_without_oculyze = y15_wide.loc[
            y15_wide["time_source"].eq("y15_session_filename"),
            ["process", "sample_id", "sample_datetime", "t_h", "time_source"],
        ].sort_values(["process", "sample_datetime"])

        oculyze_without_y15 = oculyze.loc[
            ~oculyze["sample_id_clean"].isin(set(y15_wide["sample_id"])),
            ["process", "sample_id_raw", "sample_id_clean", "sample_datetime", "t_h"],
        ].sort_values(["process", "sample_datetime"])

        negative_chem = y15_wide.loc[
            y15_wide.filter(regex="(ammonia|pan|yan|sugar|glucose|fructose)").lt(0).any(axis=1),
            [c for c in ["process", "sample_id", "sample_datetime", "t_h", "ammonia_mg_l", "pan_mg_l", "yan_calc_mg_l", "sugar_total_g_l", "glucose_g_l", "fructose_g_l"] if c in y15_wide.columns],
        ].sort_values(["process", "sample_datetime"])

        y15_value_conflicts = (
            y15_relevant.dropna(subset=["value"])
            .groupby(["sample_id", "analyte_canonical"])
            .agg(
                n_distinct_values=("value", lambda s: int(s.nunique())),
                values=("value", lambda s: "; ".join(map(str, sorted(set(s.dropna()))))),
                analytes=("analyte", lambda s: "; ".join(sorted(set(map(str, s))))),
                files=("source_file", lambda s: "; ".join(sorted(set(map(str, s))))),
            )
            .reset_index()
        )
        y15_value_conflicts = y15_value_conflicts[y15_value_conflicts["n_distinct_values"] > 1].copy()
        y15_value_conflicts["process"] = y15_value_conflicts["sample_id"].map(canonical_process_from_y15_sample)

        nitrogen_conflict_samples = set(
            y15_value_conflicts.loc[
                y15_value_conflicts["analyte_canonical"].isin(["ammonia_mg_l", "pan_mg_l"]),
                "sample_id",
            ]
        )
        if "nitrogen_qc_class" in y15_wide.columns:
            y15_wide["nitrogen_has_conflicting_y15_export"] = y15_wide["sample_id"].isin(nitrogen_conflict_samples)
            y15_wide.loc[
                y15_wide["nitrogen_has_conflicting_y15_export"]
                & y15_wide["nitrogen_qc_class"].eq("ok_raw_nonnegative"),
                "nitrogen_qc_class",
            ] = "conflicting_y15_export"
            y15_wide.loc[
                y15_wide["nitrogen_has_conflicting_y15_export"]
                & ~y15_wide["nitrogen_qc_class"].eq("conflicting_y15_export"),
                "nitrogen_qc_class",
            ] = "conflicting_y15_export_and_" + y15_wide.loc[
                y15_wide["nitrogen_has_conflicting_y15_export"]
                & ~y15_wide["nitrogen_qc_class"].eq("conflicting_y15_export"),
                "nitrogen_qc_class",
            ].astype(str)

            nitrogen_qc_cols = [
                "process",
                "sample_id",
                "sample_datetime",
                "t_h",
                "ammonia_mg_l_raw",
                "pan_mg_l_raw",
                "yan_calc_mg_l_raw",
                "ammonia_mg_l_censored",
                "pan_mg_l_censored",
                "yan_for_model_mg_l",
                "yan_for_model_physical_mg_l",
                "nitrogen_qc_class",
                "nitrogen_physical_qc_flag",
                "nitrogen_positive_rebound_before_N_pulse",
                "nitrogen_manual_zero_between_50h_and_N_pulse",
                "nitrogen_previous_N_pulse_t_h",
                "nitrogen_next_N_pulse_t_h",
                "nitrogen_has_conflicting_y15_export",
            ]
            nitrogen_qc = y15_wide[
                [c for c in nitrogen_qc_cols if c in y15_wide.columns]
            ].sort_values(["process", "sample_datetime", "sample_id"])
            write_csv_safe(nitrogen_qc, PROCESSED_DIR / "lot1_nitrogen_qc_processed.csv", index=False)
            write_csv_safe(
                y15_wide.sort_values(["process", "sample_datetime", "sample_id"]),
                PROCESSED_DIR / "lot1_y15_wide_processed.csv",
                index=False,
            )

        print("Corrected Oculyze sample IDs")
        display(corrected_ids)
        print("Y15 samples without direct Oculyze timestamp")
        display(y15_without_oculyze)
        print("Oculyze samples without Y15 chemistry")
        display(oculyze_without_y15)
        print("Negative chemistry values retained as measured and flagged")
        display(negative_chem)
        print("Y15 duplicate/conflicting values before canonical selection")
        display(y15_value_conflicts.sort_values(["process", "sample_id", "analyte_canonical"]))
        if "nitrogen_qc" in globals():
            print("Nitrogen QC table")
            display(nitrogen_qc)

        corrected_ids.to_csv(PROCESSED_DIR / "lot1_qc_corrected_oculyze_ids.csv", index=False)
        y15_without_oculyze.to_csv(PROCESSED_DIR / "lot1_qc_y15_without_oculyze_timestamp.csv", index=False)
        oculyze_without_y15.to_csv(PROCESSED_DIR / "lot1_qc_oculyze_without_y15.csv", index=False)
        negative_chem.to_csv(PROCESSED_DIR / "lot1_qc_negative_chemistry.csv", index=False)
        y15_value_conflicts.sort_values(["process", "sample_id", "analyte_canonical"]).to_csv(
            PROCESSED_DIR / "lot1_qc_y15_value_conflicts.csv", index=False
        )
        """
    ),
    md("## Process-level data previews"),
    code(
        """
        preview_cols = [
            "sample_id", "sample_datetime", "time_source", "t_h",
            "ethanol_status", "ethanol_percent_vv_final_calc", "ethanol_g_l_for_model_manual",
            "ethanol_sample_type", "sample_volume_ml", "reactor_volume_after_sample_ml",
            "sugar_total_g_l", "glucose_g_l", "fructose_g_l",
            "ammonia_mg_l", "pan_mg_l", "yan_calc_mg_l",
            "glycerol_g_l", "acetic_acid_g_l", "pyruvic_acid_mg_l", "acetaldehyde_mg_l",
        ]
        preview_cols = [c for c in preview_cols if c in y15_wide.columns]

        for process, meta in PROCESS_META.items():
            print(f"\\n{process}: {meta['design_name']}")
            chem_preview = y15_wide.loc[y15_wide["process"].eq(process), preview_cols].sort_values("sample_datetime")
            display(chem_preview)
            bio_preview = oculyze.loc[
                oculyze["process"].eq(process),
                [
                    "sample_id_clean", "sample_datetime", "t_h",
                    "cell_total_million_ml", "cell_viable_million_ml", "cell_dead_million_ml",
                    "x_total_g_l", "x_viable_g_l", "x_dead_g_l",
                    "Viability", "Budding Index", "Density",
                ],
            ].sort_values("sample_datetime")
            display(bio_preview)
        """
    ),
    md("## Plots"),
    code(
        """
        def pulse_dose_label(channel, amount, stock_volume_ml=None):
            channel = str(channel).strip().upper()
            if pd.isna(amount):
                dose = ""
            elif channel == "N":
                dose = f"+{float(amount) * 1000:g} mg/L YAN"
            elif channel == "G":
                dose = f"+{float(amount):g} g/L G"
            elif channel == "F":
                dose = f"+{float(amount):g} g/L F"
            elif channel == "E":
                dose = f"+{float(amount):g} g/L EtOH"
            elif channel == "X":
                dose = f"+{float(amount):g} g/L X"
            else:
                dose = f"+{float(amount):g}"

            if stock_volume_ml is not None and pd.notna(stock_volume_ml):
                return f"{channel}\\n{dose}\\n{float(stock_volume_ml):.1f} mL"
            return f"{channel}\\n{dose}"


        def axes_for_input_channel(axes, channel):
            channel = str(channel).strip().upper()
            mapping = {
                "X": [axes[2]],      # biomass panel
                "G": [axes[4]],      # sugars panel
                "F": [axes[4]],      # sugars panel
                "N": [axes[5]],      # nitrogen panel
                "E": [axes[6]],      # liquid products panel; ethanol plot below is also available
            }
            return mapping.get(channel, [])


        def draw_pulse_marker(ax, t_h, label, color, planned=True):
            if pd.isna(t_h):
                return
            linestyle = "--" if planned else "-"
            alpha = 0.55 if planned else 0.85
            lw = 1.2 if planned else 1.6
            ax.axvline(t_h, color=color, linestyle=linestyle, alpha=alpha, lw=lw)
            y0, y1 = ax.get_ylim()
            y = y1 - 0.03 * (y1 - y0) if planned else y0 + 0.03 * (y1 - y0)
            va = "top" if planned else "bottom"
            prefix = "P" if planned else "A"
            ax.text(
                t_h,
                y,
                f"{prefix} {label}",
                color=color,
                rotation=90,
                va=va,
                ha="right" if planned else "left",
                fontsize=7,
                alpha=0.95,
                bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 1.0},
            )


        def annotate_input_events(axes, process: str):
            if "planned_input_events" in globals() and not planned_input_events.empty:
                planned = planned_input_events[planned_input_events["process"].eq(process)]
            else:
                planned = pd.DataFrame()

            if "actual_input_events" in globals() and not actual_input_events.empty:
                actual = actual_input_events[
                    actual_input_events["process"].eq(process)
                    & actual_input_events["has_actual_time"]
                ]
            else:
                actual = pd.DataFrame()

            for _, row in planned.iterrows():
                color = INPUT_CHANNEL_COLORS.get(str(row["channel"]), "black")
                t = row["planned_design_t_h"]
                if pd.isna(t):
                    continue
                label = pulse_dose_label(
                    row["channel"],
                    row.get("planned_amount", np.nan),
                    row.get("planned_stock_volume_ml", np.nan),
                )
                for ax in axes_for_input_channel(axes, row["channel"]):
                    draw_pulse_marker(ax, t, label, color, planned=True)

            for _, row in actual.iterrows():
                color = INPUT_CHANNEL_COLORS.get(str(row["channel"]), "black")
                t = row["actual_plot_t_h"]
                if pd.isna(t):
                    continue
                label = pulse_dose_label(
                    row["channel"],
                    row.get("actual_amount_for_plot", np.nan),
                    row.get("actual_stock_volume_ml", np.nan),
                )
                for ax in axes_for_input_channel(axes, row["channel"]):
                    draw_pulse_marker(ax, t, label, color, planned=False)


        def plot_process(process: str):
            meta = PROCESS_META[process]
            fig, axes = plt.subplots(4, 2, figsize=(15, 15), sharex=False)
            axes = axes.ravel()

            co2 = co2_10min[co2_10min["process"].eq(process)].sort_values("t_h")
            temp = temp_10min[temp_10min["process"].eq(process)].sort_values("t_h")
            bio = oculyze[oculyze["process"].eq(process)].sort_values("t_h")
            chem = y15_wide[y15_wide["process"].eq(process)].sort_values("t_h")

            axes[0].plot(co2["t_h"], co2["flow_filt_sccm"], lw=1.0, color="tab:blue", label="raw")
            if "flow_filt_sccm_v0_corrected" in co2:
                axes[0].plot(
                    co2["t_h"],
                    co2["flow_filt_sccm_v0_corrected"],
                    lw=1.2,
                    color="tab:cyan",
                    label="2 L corrected",
                )
            axes[0].set_title("Filtered CO2 flow")
            axes[0].set_ylabel("sccm")
            axes[0].legend(frameon=False)

            axes[1].plot(temp["t_h"], temp["T"], lw=1.1, label="T", color="tab:red")
            axes[1].plot(temp["t_h"], temp["SP"], lw=1.1, label="SP", color="black", alpha=0.8)
            axes[1].set_title("Temperature")
            axes[1].set_ylabel("deg C")
            axes[1].legend(frameon=False)

            axes[2].plot(bio["t_h"], bio["x_total_g_l"], marker="o", label="total", color="tab:gray")
            axes[2].plot(bio["t_h"], bio["x_viable_g_l"], marker="o", label="viable", color="tab:green")
            axes[2].plot(bio["t_h"], bio["x_dead_g_l"], marker="o", label="dead", color="tab:brown")
            axes[2].set_title("Biomass from Oculyze")
            axes[2].set_ylabel("g/L dry biomass")
            axes[2].legend(frameon=False)

            axes[3].plot(bio["t_h"], bio["Density"], marker="o", color="tab:purple")
            axes[3].set_title("Density")
            axes[3].set_ylabel("density")

            if "sugar_total_g_l" in chem:
                axes[4].plot(chem["t_h"], chem["sugar_total_g_l"], marker="o", label="G+F", color="tab:blue")
            if "glucose_g_l" in chem:
                axes[4].plot(chem["t_h"], chem["glucose_g_l"], marker="o", label="glucose", color="tab:orange")
            if "fructose_g_l" in chem:
                axes[4].plot(chem["t_h"], chem["fructose_g_l"], marker="o", label="fructose", color="tab:green")
            axes[4].set_title("Sugars")
            axes[4].set_ylabel("g/L")
            axes[4].legend(frameon=False)

            for col, label, color in [
                ("ammonia_mg_l_censored", "NH4 censored", "tab:blue"),
                ("pan_mg_l_censored", "PAN censored", "tab:orange"),
                ("yan_for_model_physical_mg_l", "YAN physical QC", "tab:green"),
            ]:
                if col in chem:
                    axes[5].plot(chem["t_h"], chem[col], marker="o", label=label, color=color)
            if "yan_for_model_physical_mg_l" not in chem and "yan_for_model_mg_l" in chem:
                axes[5].plot(
                    chem["t_h"],
                    chem["yan_for_model_mg_l"],
                    marker="o",
                    label="YAN for model",
                    color="tab:green",
                )
            for raw_col, color in [("ammonia_mg_l", "tab:blue"), ("pan_mg_l", "tab:orange")]:
                if raw_col in chem:
                    neg = chem[chem[raw_col] < 0]
                    if not neg.empty:
                        axes[5].scatter(
                            neg["t_h"],
                            np.zeros(len(neg)),
                            marker="x",
                            s=50,
                            color=color,
                            label=f"{raw_col.replace('_mg_l', '')} raw < 0",
                        )
            if "nitrogen_qc_class" in chem:
                severe = chem[chem["nitrogen_qc_class"].astype(str).str.contains("severe_negative", na=False)]
                conflict = chem[chem["nitrogen_qc_class"].astype(str).str.contains("conflicting_y15_export", na=False)]
                if not severe.empty:
                    axes[5].scatter(
                        severe["t_h"],
                        np.zeros(len(severe)),
                        marker="s",
                        facecolors="none",
                        edgecolors="black",
                        s=85,
                        linewidths=1.3,
                        label="severe raw negative: model as near-zero",
                    )
                if not conflict.empty:
                    axes[5].scatter(
                        conflict["t_h"],
                        conflict["yan_for_model_mg_l"].fillna(0.0),
                        marker="D",
                        facecolors="none",
                        edgecolors="tab:red",
                        s=75,
                        linewidths=1.2,
                        label="conflicting Y15 export",
                    )
                if "nitrogen_positive_rebound_before_N_pulse" in chem:
                    rebound = chem[chem["nitrogen_positive_rebound_before_N_pulse"].fillna(False).astype(bool)]
                    if not rebound.empty:
                        axes[5].scatter(
                            rebound["t_h"],
                            rebound["yan_for_model_mg_l"].fillna(0.0),
                            marker="^",
                            facecolors="none",
                            edgecolors="tab:purple",
                            s=95,
                            linewidths=1.4,
                            label="pre-N rebound: artifact, set to zero",
                        )
            axes[5].set_title("Nitrogen, negative raw values censored at zero")
            axes[5].set_ylabel("mg/L")
            axes[5].legend(frameon=False)

            for col, label, color in [
                ("glycerol_g_l", "glycerol", "tab:green"),
                ("acetic_acid_g_l", "acetic acid", "tab:red"),
            ]:
                if col in chem:
                    axes[6].plot(chem["t_h"], chem[col], marker="o", label=label, color=color)
            axes[6].set_title("Liquid products")
            axes[6].set_ylabel("g/L")
            axes[6].legend(frameon=False)

            for col, label, color in [
                ("pyruvic_acid_mg_l", "pyruvic acid", "tab:purple"),
                ("acetaldehyde_mg_l", "acetaldehyde", "tab:brown"),
            ]:
                if col in chem:
                    axes[7].plot(chem["t_h"], chem[col], marker="o", label=label, color=color)
            axes[7].set_title("Secondary metabolites")
            axes[7].set_ylabel("mg/L")
            axes[7].legend(frameon=False)

            for ax in axes:
                ax.grid(True, alpha=0.25)
                ax.set_xlabel("time since first Oculyze sample, h")

            annotate_input_events(axes, process)

            fig.suptitle(f"{process} - {meta['design_name']}", y=1.01, fontsize=14)
            fig.tight_layout()
            out = FIG_DIR / f"lot1_{process}_preview.png"
            fig.savefig(out, dpi=180, bbox_inches="tight")
            return fig, out

        plot_outputs = []
        for process in PROCESS_META:
            fig, out = plot_process(process)
            plot_outputs.append(out)
            plt.show()

        print("Saved figures:")
        for out in plot_outputs:
            print(out.relative_to(ROOT))
        """
    ),
    md("## Operational input timeline plots"),
    code(
        """
        def plot_input_timeline(process: str):
            meta = PROCESS_META[process]
            planned = (
                planned_input_events[planned_input_events["process"].eq(process)].copy()
                if "planned_input_events" in globals() and not planned_input_events.empty
                else pd.DataFrame()
            )
            actual = (
                actual_input_events[
                    actual_input_events["process"].eq(process)
                    & actual_input_events["has_actual_time"]
                ].copy()
                if "actual_input_events" in globals() and not actual_input_events.empty
                else pd.DataFrame()
            )
            samples = (
                volume_balance[volume_balance["process"].eq(process)].copy()
                if "volume_balance" in globals() and not volume_balance.empty
                else pd.DataFrame()
            )

            fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
            y_map = {"N": 4, "G": 3, "F": 2, "E": 1, "X": 0}

            if not planned.empty:
                for _, row in planned.iterrows():
                    channel = str(row["channel"])
                    y = y_map.get(channel, -1)
                    color = INPUT_CHANNEL_COLORS.get(channel, "black")
                    axes[0].vlines(row["planned_design_t_h"], y - 0.3, y + 0.3, color=color, linestyles="--", lw=2)
                    axes[0].scatter(row["planned_design_t_h"], y, color=color, marker="o", s=60)
                    axes[0].text(
                        row["planned_design_t_h"],
                        y + 0.35,
                        f"P {channel}\\n{row['planned_amount']:g}",
                        color=color,
                        ha="center",
                        va="bottom",
                        fontsize=8,
                    )

            if not actual.empty:
                for _, row in actual.iterrows():
                    channel = str(row["channel"])
                    y = y_map.get(channel, -1)
                    color = INPUT_CHANNEL_COLORS.get(channel, "black")
                    axes[0].vlines(row["actual_plot_t_h"], y - 0.3, y + 0.3, color=color, linestyles="-", lw=2.5)
                    axes[0].scatter(row["actual_plot_t_h"], y, color=color, marker="X", s=90)
                    amount = row["actual_amount_for_plot"]
                    label_amount = "" if pd.isna(amount) else f"\\n{amount:g}"
                    axes[0].text(
                        row["actual_plot_t_h"],
                        y - 0.42,
                        f"A {channel}{label_amount}",
                        color=color,
                        ha="center",
                        va="top",
                        fontsize=8,
                    )

            axes[0].set_yticks(list(y_map.values()))
            axes[0].set_yticklabels(list(y_map.keys()))
            axes[0].set_ylabel("input channel")
            axes[0].set_title("Operational inputs: planned dashed, actual solid if logged")
            axes[0].grid(True, axis="x", alpha=0.25)

            if not samples.empty:
                full = samples[samples["ethanol_sample_type"].eq("ethanol_50mL")]
                small = samples[samples["ethanol_sample_type"].eq("small_5mL")]
                axes[1].scatter(small["t_h"], np.zeros(len(small)), color="tab:gray", s=25, label="small sample")
                axes[1].scatter(full["t_h"], np.ones(len(full)), color="tab:red", s=35, label="ethanol/full sample")
            axes[1].set_yticks([0, 1])
            axes[1].set_yticklabels(["small", "full"])
            axes[1].set_xlabel("time since first Oculyze sample, h")
            axes[1].set_title("Sampling events from executed analytical data")
            axes[1].grid(True, axis="x", alpha=0.25)
            axes[1].legend(frameon=False, loc="upper right")

            fig.suptitle(f"{process} - {meta['design_name']}", y=1.02)
            fig.tight_layout()
            out = FIG_DIR / f"lot1_{process}_operational_inputs.png"
            fig.savefig(out, dpi=180, bbox_inches="tight")
            return fig, out


        input_plot_outputs = []
        for process in PROCESS_META:
            fig, out = plot_input_timeline(process)
            input_plot_outputs.append(out)
            plt.show()

        print("Saved operational input figures:")
        for out in input_plot_outputs:
            print(out.relative_to(ROOT))
        """
    ),
    md("## Ethanol and reactor-volume plots"),
    code(
        """
        if "volume_balance" in globals() and not volume_balance.empty:
            ethanol_volume_outputs = []
            for process, meta in PROCESS_META.items():
                dfp = y15_wide[y15_wide["process"].eq(process)].sort_values("t_h")
                fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

                measured = dfp[dfp["ethanol_status"].eq("measured")]
                pending = dfp[dfp["ethanol_status"].eq("pending_TBD")]
                if not measured.empty:
                    axes[0].plot(
                        measured["t_h"],
                        measured["ethanol_g_l_for_model_manual"],
                        marker="o",
                        lw=1.4,
                        color="tab:blue",
                        label="measured ethanol",
                    )
                if not pending.empty:
                    axes[0].scatter(
                        pending["t_h"],
                        np.zeros(len(pending)),
                        marker="x",
                        s=70,
                        color="tab:orange",
                        label="TBD ethanol sample",
                    )
                axes[0].set_ylabel("ethanol, g/L")
                axes[0].set_title(f"{process} ethanol entries")
                axes[0].legend(frameon=False)
                axes[0].grid(True, alpha=0.25)

                axes[1].step(
                    dfp["t_h"],
                    dfp["reactor_volume_after_sample_ml"],
                    where="post",
                    color="tab:green",
                    lw=1.5,
                    label="volume after sample",
                )
                ethanol_samples = dfp[dfp["ethanol_sample_type"].eq("ethanol_50mL")]
                small_samples = dfp[dfp["ethanol_sample_type"].eq("small_5mL")]
                axes[1].scatter(
                    small_samples["t_h"],
                    small_samples["reactor_volume_after_sample_ml"],
                    color="tab:gray",
                    s=25,
                    label="5 mL sample",
                )
                axes[1].scatter(
                    ethanol_samples["t_h"],
                    ethanol_samples["reactor_volume_after_sample_ml"],
                    color="tab:red",
                    s=35,
                    label="50 mL sample",
                )
                axes[1].set_ylabel("reactor volume, mL")
                axes[1].set_xlabel("time since first Oculyze sample, h")
                axes[1].set_title("Sampling-loss volume balance")
                axes[1].legend(frameon=False)
                axes[1].grid(True, alpha=0.25)

                fig.suptitle(f"{process} - {meta['design_name']}", y=1.02)
                fig.tight_layout()
                out = FIG_DIR / f"lot1_{process}_ethanol_volume.png"
                fig.savefig(out, dpi=180, bbox_inches="tight")
                ethanol_volume_outputs.append(out)
                plt.show()

            print("Saved ethanol/volume figures:")
            for out in ethanol_volume_outputs:
                print(out.relative_to(ROOT))
        else:
            print("No volume balance available for plotting.")
        """
    ),
    md(
        """
        ## Notes for downstream calibration

        - Use `sample_id_clean` for Oculyze/Y15 joins, but retain `sample_id_raw`
          for auditability.
        - Use `time_source` to decide whether a chemical point should be trusted as
          an exact process sample time. `inferred_lot_sample_sequence` is based on
          the same sample number in the rest of the lot; `y15_session_filename` is
          a weaker fallback.
        - Biomass conversion uses 30 pg dry mass per cell:
          `1 million cells/mL = 0.03 g/L = 0.03 kg/m3`.
        - Negative low-range Y15 readings are retained and flagged, not clipped.
          The nitrogen preview plot censors negative raw readings at zero to avoid
          visually implying negative concentrations. The raw values remain in
          `*_raw` columns and in the QC tables.
        - CO2 correction is volume-normalization, not sensor filtering. The raw
          filtered signal remains as `flow_filt_sccm`. For concentration-based ODE
          calibration, prefer `flow_filt_sccm_per_l` or convert it to mol/L/h with
          the appropriate gas conditions. For comparison against a total-flow model
          normalized to the initial reactor volume, use
          `flow_filt_sccm_v0_corrected`.
        """
    ),
]


nbf.write(nb, NOTEBOOK_PATH)
print(f"Wrote {NOTEBOOK_PATH}")
