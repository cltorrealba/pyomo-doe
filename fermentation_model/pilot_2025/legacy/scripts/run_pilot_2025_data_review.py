from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import pilot_2025_data_loader as loader


RESULTS_DIR = SCRIPT_DIR / "results" / "data_review"
NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_data_review.ipynb"
EXECUTED_NOTEBOOK_PATH = SCRIPT_DIR / "pilot_2025_data_review.executed.ipynb"

PRIMARY_STATES = [
    ("density", "Density"),
    ("S_GF_g_l", "Glucose + fructose (g/L)"),
    ("G_g_l", "Glucose (g/L)"),
    ("F_g_l", "Fructose (g/L)"),
    ("YAN_mg_l", "YAN (mg/L)"),
    ("PAN_mg_l", "PAN (mg/L)"),
    ("NH4_mg_l", "Ammonia (mg/L)"),
    ("X_viable_kg_m3", "Viable biomass (kg/m3)"),
    ("X_dry_weight_g_l", "Dry weight (g/L)"),
    ("E_g_l", "Ethanol (g/L)"),
    ("glycerol_g_l", "Glycerol (g/L)"),
]

SECONDARY_AROMA_STATES = [
    ("pyruvic_acid_mg_l", "Pyruvic acid (mg/L)"),
    ("acetaldehyde_mg_l", "Acetaldehyde (mg/L)"),
    ("acetic_acid_g_l", "Acetic acid (g/L)"),
    ("ethyl_acetate_total", "Ethyl acetate total"),
    ("ethyl_acetate_condensate", "Ethyl acetate condensate"),
    ("isoamyl_acetate_total", "Isoamyl acetate total"),
    ("isoamyl_acetate_condensate", "Isoamyl acetate condensate"),
    ("ethyl_octanoate_total", "Ethyl octanoate total"),
    ("ethyl_octanoate_condensate", "Ethyl octanoate condensate"),
]


def _safe_name(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)


def plot_batch_overview(data: pd.DataFrame, co2: pd.DataFrame, batch: str, output_dir: Path) -> Path:
    group = data[data["batch"].eq(batch)].sort_values("time_h")
    co2_group = co2[co2["batch"].eq(batch)].sort_values("time_h")
    fig, axes = plt.subplots(5, 2, figsize=(14, 16), sharex=True)
    axes = axes.ravel()

    panels = [
        ("temperature_c", "Temperature (C)"),
        ("S_GF_g_l", "G+F (g/L)"),
        ("G_g_l", "Glucose (g/L)"),
        ("F_g_l", "Fructose (g/L)"),
        ("YAN_mg_l", "YAN (mg/L)"),
        ("X_viable_kg_m3", "Viable X (kg/m3)"),
        ("E_g_l", "Ethanol (g/L)"),
        ("glycerol_g_l", "Glycerol (g/L)"),
        ("ethyl_acetate_total", "Ethyl acetate"),
        ("isoamyl_acetate_total", "Isoamyl acetate"),
    ]
    for ax, (column, ylabel) in zip(axes, panels):
        if column in group:
            ax.plot(group["time_h"], group[column], marker="o", linewidth=1.4)
        if column == "YAN_mg_l" and "N_pulse_mg_l" in group:
            pulses = group[group["N_pulse_mg_l"].fillna(0.0).gt(0.0)]
            for _, pulse in pulses.iterrows():
                ax.axvline(float(pulse["time_h"]), color="tab:red", linestyle="--", alpha=0.35)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Time (h)")
    axes[-2].set_xlabel("Time (h)")
    fig.suptitle(f"Pilot fermentation {batch}: process and key aroma observations", y=0.995)
    fig.tight_layout()
    path = output_dir / f"batch_{_safe_name(batch)}_overview.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)

    if not co2_group.empty and co2_group["time_h"].notna().any():
        fig, ax = plt.subplots(figsize=(12, 4))
        clean = co2_group[co2_group["time_h"].notna()].copy()
        ax.plot(clean["time_h"], clean["co2_raw"], linewidth=0.8, alpha=0.8)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("CO2 sensor raw signal")
        ax.set_title(f"Pilot fermentation {batch}: online CO2")
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        co2_path = output_dir / f"batch_{_safe_name(batch)}_co2.png"
        fig.savefig(co2_path, dpi=160)
        plt.close(fig)
    return path


def plot_density_fit(data: pd.DataFrame, fit: loader.DensitySugarFit, output_dir: Path) -> Path:
    df = data[["density", "S_GF_g_l", "batch"]].dropna()
    fig, ax = plt.subplots(figsize=(8, 6))
    for batch, group in df.groupby("batch"):
        ax.scatter(group["density"], group["S_GF_g_l"], label=batch, s=24, alpha=0.75)
    if np.isfinite(fit.slope):
        x = np.linspace(float(df["density"].min()), float(df["density"].max()), 100)
        ax.plot(x, fit.predict(x), color="black", linewidth=2, label="linear fit")
    ax.set_xlabel("Density")
    ax.set_ylabel("Glucose + fructose (g/L)")
    ax.set_title(f"Density-to-total-sugar fit: R2={fit.r2:.3f}, RMSE={fit.rmse_g_l:.2f} g/L")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path = output_dir / "density_total_sugar_fit.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def write_report(
    data: pd.DataFrame,
    raw_co2: pd.DataFrame,
    curated_co2: pd.DataFrame,
    co2_decisions: pd.DataFrame,
    coverage: pd.DataFrame,
    fit: loader.DensitySugarFit,
) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / "pilot_2025_data_review_report.md"
    co2_summary = (
        curated_co2.groupby("batch")
        .agg(
            n_rows=("co2_raw", "size"),
            t_original_min_h=("time_h_original", "min"),
            t_original_max_h=("time_h_original", "max"),
            t_effective_min_h=("time_h_effective", "min"),
            t_effective_max_h=("time_h_effective", "max"),
            co2_min=("co2_raw", "min"),
            co2_median=("co2_raw", "median"),
            co2_max=("co2_raw", "max"),
        )
        .reset_index()
        if not curated_co2.empty
        else pd.DataFrame()
    )
    with path.open("w", encoding="utf-8") as f:
        f.write("# Pilot 2025 data review\n\n")
        f.write("## Scope\n\n")
        f.write(
            "This first pass loads the pilot-scale natural-must fermentations, "
            "normalizes units where assumptions are clear, aligns online CO2 files "
            "to fermentation start times, and checks which states are informative for calibration and MBDoE.\n\n"
        )
        f.write("## Unit assumptions\n\n")
        f.write("- `ETANOL` is treated as g/L. Values reach about 95-98, which is plausible as g/L and impossible as % v/v for wine.\n")
        f.write("- Viable cell counts are converted with 30 pg/cell, so 1 million cells/mL = 0.03 kg/m3.\n")
        f.write("- `PAN`, `AMMONIA`, `YAN`, and pyruvic/acetaldehyde observations are kept as mg/L.\n")
        f.write("- Aroma `*_total` columns are interpreted as wine-retained plus condenser-accumulated equivalent concentration in mg/L.\n")
        f.write("- Aroma `*_condensate` columns are interpreted as equivalent must concentration accumulated in the condenser in mg/L.\n\n")
        f.write("## Data volume\n\n")
        f.write(f"- Fermentation sheets loaded: {data['batch'].nunique()}.\n")
        f.write(f"- Total calibration rows: {len(data)}.\n")
        f.write(f"- Raw CO2 sensor batches loaded: {raw_co2['batch'].nunique() if not raw_co2.empty else 0}.\n")
        f.write(f"- Curated CO2 sensor batches retained for calibration: {curated_co2['batch'].nunique() if not curated_co2.empty else 0}.\n\n")
        f.write("## Density to sugar proxy\n\n")
        f.write(
            f"Linear fit: `S_GF_g_l = {fit.intercept:.3f} + {fit.slope:.3f} density`, "
            f"`R2 = {fit.r2:.3f}`, `RMSE = {fit.rmse_g_l:.2f} g/L`, `n = {fit.n_points}`.\n\n"
        )
        f.write("## Observation coverage\n\n")
        f.write(coverage.to_markdown(index=False))
        f.write("\n\n")
        if not co2_summary.empty:
            f.write("## CO2 curation decisions\n\n")
            f.write(co2_decisions.to_markdown(index=False))
            f.write("\n\n")
            f.write("## CO2 sensor coverage\n\n")
            f.write(co2_summary.to_markdown(index=False))
            f.write("\n\n")
        f.write("## Initial interpretation\n\n")
        f.write(
            "- The pilot data are mostly natural-must process trajectories, so design variables are less flexible than in synthetic-lab DOE.\n"
            "- The strongest immediate calibration value is the richer time series for ethanol, glycerol, pyruvate, acetaldehyde, ethyl acetate, isoamyl acetate, and condenser aroma accumulation.\n"
            "- Online CO2 is available for a subset of batches. Batches 25150 and 25151 are excluded from calibration. Batch 25171 is trimmed to the sustained CO2 activation point and carries both original and effective time columns.\n"
            "- The next step is to map these normalized batches into the extended model residual structure and run calibration/estimability with the existing reduced secondary/aroma framework.\n"
        )
    return path


def create_notebook() -> None:
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Pilot 2025 Data Review\n\n"
            "This notebook is the pilot-scale counterpart of the laboratory MBDoE workflow. "
            "The objective of this first stage is not parameter estimation yet. The objective is to verify that the data can be loaded, normalized, aligned in time, and interpreted before the calibration and estimability pipeline is re-used.\n\n"
            "The pilot system uses natural must, so the future design space is narrower than the synthetic-must laboratory campaign. In exchange, the dataset has richer observations for the aroma mass balance, especially ethyl acetate and isoamyl acetate."
        ),
        nbformat.v4.new_markdown_cell(
            "## Model-Relevant Data Mapping\n\n"
            "The extended fermentation model uses the following observation groups:\n\n"
            "- Primary fermentation states: viable biomass, glucose, fructose, assimilable nitrogen, ethanol, glycerol, and temperature input.\n"
            "- Secondary states: pyruvic acid and acetaldehyde. Acetic acid is included if available.\n"
            "- Aroma states: total-equivalent and condenser-equivalent observations for ethyl acetate, isoamyl acetate, and ethyl octanoate.\n"
            "- Online gas information: CO2 sensor trajectories for the subset of batches with sensor files.\n\n"
            "Aroma `*_total` columns are interpreted as the total equivalent concentration in wine plus condenser. Aroma `*_condensate` columns are interpreted as the equivalent must concentration accumulated in the condenser. Therefore the retained wine pool can be reconstructed as `total - condensate` when both are measured."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import sys\n"
            "import pandas as pd\n"
            "import matplotlib.pyplot as plt\n\n"
            "SCRIPT_DIR = Path.cwd()\n"
            "if SCRIPT_DIR.name != 'pilot_2025':\n"
            "    SCRIPT_DIR = Path('fermentation_model/pilot_2025').resolve()\n"
            "sys.path.insert(0, str(SCRIPT_DIR))\n"
            "import pilot_2025_data_loader as loader\n\n"
            "data = loader.load_pilot_calibration_data()\n"
            "raw_co2 = loader.load_all_co2_sensor_data(calibration_data=data)\n"
            "co2, co2_decisions = loader.curate_co2_sensor_data(raw_co2)\n"
            "coverage = loader.observation_coverage(data)\n"
            "density_fit = loader.fit_density_total_sugar(data)\n"
            "data.shape, raw_co2.shape, co2.shape"
        ),
        nbformat.v4.new_markdown_cell("## Fermentation Sheets and Observation Coverage"),
        nbformat.v4.new_code_cell("coverage"),
        nbformat.v4.new_markdown_cell(
            "## Unit Checks\n\n"
            "The `ETANOL` column is interpreted as g/L. Its final values are close to 95 g/L, which is realistic for wine. Interpreting the same values as % v/v would be physically impossible."
        ),
        nbformat.v4.new_code_cell(
            "unit_cols = ['batch', 'time_h', 'E_g_l', 'E_percent_vv', 'S_GF_g_l', 'density', 'X_viable_kg_m3', 'X_dry_weight_g_l']\n"
            "data[unit_cols].groupby('batch').agg(['count', 'min', 'median', 'max'])"
        ),
        nbformat.v4.new_markdown_cell(
            "## Density as Total Sugar Proxy\n\n"
            "Density is useful because it is measured more frequently than glucose and fructose. "
            "The linear relationship below estimates the information loss incurred when using density as a proxy for total sugar."
        ),
        nbformat.v4.new_code_cell(
            "print(density_fit)\n"
            "df = data[['batch', 'density', 'S_GF_g_l']].dropna()\n"
            "fig, ax = plt.subplots(figsize=(8, 6))\n"
            "for batch, group in df.groupby('batch'):\n"
            "    ax.scatter(group['density'], group['S_GF_g_l'], label=batch, s=24)\n"
            "x = pd.Series([df['density'].min(), df['density'].max()])\n"
            "ax.plot(x, density_fit.predict(x), color='black', linewidth=2)\n"
            "ax.set_xlabel('Density')\n"
            "ax.set_ylabel('Glucose + fructose (g/L)')\n"
            "ax.grid(True, alpha=0.25)\n"
            "ax.legend(ncol=2, fontsize=8)\n"
            "plt.show()"
        ),
        nbformat.v4.new_markdown_cell(
            "## Batch Trajectories\n\n"
            "These plots are intended to detect unit problems, missing-data patterns, and whether pilot trajectories excite the same parameter directions as the laboratory designs."
        ),
        nbformat.v4.new_code_cell(
            "states = [\n"
            "    ('temperature_c', 'Temperature (C)'), ('S_GF_g_l', 'G+F (g/L)'),\n"
            "    ('YAN_mg_l', 'YAN (mg/L)'), ('X_viable_kg_m3', 'Viable X (kg/m3)'),\n"
            "    ('E_g_l', 'Ethanol (g/L)'), ('glycerol_g_l', 'Glycerol (g/L)'),\n"
            "    ('pyruvic_acid_mg_l', 'Pyruvic acid (mg/L)'), ('acetaldehyde_mg_l', 'Acetaldehyde (mg/L)'),\n"
            "    ('ethyl_acetate_total', 'Ethyl acetate total'), ('ethyl_acetate_condensate', 'Ethyl acetate condensate'),\n"
            "    ('isoamyl_acetate_total', 'Isoamyl acetate total'), ('isoamyl_acetate_condensate', 'Isoamyl acetate condensate'),\n"
            "]\n"
            "for batch, group in data.groupby('batch'):\n"
            "    fig, axes = plt.subplots(6, 2, figsize=(14, 17), sharex=True)\n"
            "    for ax, (col, ylabel) in zip(axes.ravel(), states):\n"
            "        ax.plot(group['time_h'], group[col], marker='o', linewidth=1.2)\n"
            "        if col == 'YAN_mg_l':\n"
            "            pulses = group[group['N_pulse_mg_l'].fillna(0).gt(0)]\n"
            "            for _, pulse in pulses.iterrows():\n"
            "                ax.axvline(float(pulse['time_h']), color='tab:red', linestyle='--', alpha=0.35)\n"
            "        ax.set_ylabel(ylabel)\n"
            "        ax.grid(True, alpha=0.25)\n"
            "    axes[-1, 0].set_xlabel('Time (h)')\n"
            "    axes[-1, 1].set_xlabel('Time (h)')\n"
            "    fig.suptitle(f'Pilot fermentation {batch}', y=1.0)\n"
            "    fig.tight_layout()\n"
            "    plt.show()"
        ),
        nbformat.v4.new_markdown_cell(
            "## Online CO2 Sensor Files\n\n"
            "The raw sensor files are first aligned against the first timestamp in the corresponding fermentation sheet. "
            "For calibration, `25150` and `25151` are excluded. For `25171`, the pre-activation segment is removed and an effective CO2 time is created with the sustained activation point as `t = 0`."
        ),
        nbformat.v4.new_code_cell(
            "co2_decisions"
        ),
        nbformat.v4.new_code_cell(
            "co2_summary = co2.groupby('batch').agg(\n"
            "    n_rows=('co2_raw', 'size'), t_original_min_h=('time_h_original', 'min'), t_original_max_h=('time_h_original', 'max'),\n"
            "    t_effective_min_h=('time_h_effective', 'min'), t_effective_max_h=('time_h_effective', 'max'),\n"
            "    co2_min=('co2_raw', 'min'), co2_median=('co2_raw', 'median'), co2_max=('co2_raw', 'max')\n"
            ").reset_index()\n"
            "co2_summary"
        ),
        nbformat.v4.new_code_cell(
            "for batch, group in co2.groupby('batch'):\n"
            "    fig, ax = plt.subplots(figsize=(12, 3))\n"
            "    ax.plot(group['time_h_effective'], group['co2_raw'], linewidth=0.8)\n"
            "    ax.set_title(f'CO2 sensor curated signal: {batch}')\n"
            "    ax.set_xlabel('Effective CO2 time (h)')\n"
            "    ax.set_ylabel('Raw CO2 signal')\n"
            "    ax.grid(True, alpha=0.25)\n"
            "    plt.show()"
        ),
        nbformat.v4.new_markdown_cell(
            "## Immediate Technical Conclusions\n\n"
            "1. The pilot workbook is usable for the extended model workflow.\n"
            "2. Ethanol should be treated as g/L in this dataset.\n"
            "3. The pilot data are especially valuable for ethyl acetate and isoamyl acetate because these are measured repeatedly in several fermentations.\n"
            "4. The design problem should be reformulated for natural must: temperature policy and nutrient additions are realistic control levers, while arbitrary initial glucose/fructose composition and sugar pulses are not primary levers.\n"
            "5. The next notebook should reuse the reduced secondary/aroma model residual structure, add pilot batches as a separate dataset group, and compare estimability against the laboratory natural/synthetic datasets."
        ),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    data = loader.load_pilot_calibration_data()
    raw_co2 = loader.load_all_co2_sensor_data(calibration_data=data)
    curated_co2, co2_decisions = loader.curate_co2_sensor_data(raw_co2)
    coverage = loader.observation_coverage(data)
    density_fit = loader.fit_density_total_sugar(data)

    data.to_csv(RESULTS_DIR / "pilot_2025_normalized_data.csv", index=False)
    raw_co2.to_csv(RESULTS_DIR / "pilot_2025_co2_sensor_data_raw.csv", index=False)
    curated_co2.to_csv(RESULTS_DIR / "pilot_2025_co2_sensor_data_curated.csv", index=False)
    co2_decisions.to_csv(RESULTS_DIR / "pilot_2025_co2_curation_decisions.csv", index=False)
    coverage.to_csv(RESULTS_DIR / "pilot_2025_observation_coverage.csv", index=False)
    plot_density_fit(data, density_fit, plot_dir)
    for batch in sorted(data["batch"].dropna().unique()):
        plot_batch_overview(data, curated_co2, str(batch), plot_dir)
    report = write_report(data, raw_co2, curated_co2, co2_decisions, coverage, density_fit)
    create_notebook()
    print(f"[done] normalized data: {RESULTS_DIR / 'pilot_2025_normalized_data.csv'}")
    print(f"[done] report: {report}")
    print(f"[done] notebook: {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
