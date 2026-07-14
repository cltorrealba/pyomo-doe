from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd

SHARED_DIR = Path(__file__).resolve().parent
FERMENTATION_MODEL_DIR = SHARED_DIR.parent
SCRIPT_DIR = FERMENTATION_MODEL_DIR
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

from shared.new_must_data_loader import (
    AROMA_COLUMNS,
    NEW_MUST_SOURCES,
    SECONDARY_METABOLITE_COLUMNS,
    _clean_columns,
    _readable_excel_path,
    load_new_must_data,
    normalize_must_sheet,
)
from shared.paths import DATA_DIR, SHARED_RESULTS_DIR

RESULTS_DIR = SHARED_RESULTS_DIR / "secondary_metabolite_data_review"
PLOT_DIR = RESULTS_DIR / "plots"
NOTEBOOK_PATH = (
    SHARED_DIR
    / "notebooks"
    / "fermentation_secondary_metabolite_data_review.ipynb"
)
CALIBRATION_FILE = DATA_DIR / "Piloto 2025" / "Calibration_data_vl3.xlsx"

SECONDARY_PLOT_COLUMNS = [
    "YAN_mg_l",
    "PAN_mg_l",
    "NH4_mg_l",
    "pyruvic_acid",
    "pyruvic_acid_2",
    "acetaldehyde",
    "acetic_acid",
    "DO_mg_l",
]
AROMA_PLOT_COLUMNS = list(AROMA_COLUMNS)
CORE_COLUMNS = ["X_viable_kg_m3", "G_g_l", "F_g_l", "E_g_l", "glycerol_g_l"]


def load_historical_calibration_data() -> pd.DataFrame:
    if not CALIBRATION_FILE.exists():
        return pd.DataFrame()
    frames = []
    read_path = _readable_excel_path(CALIBRATION_FILE)
    xls = pd.ExcelFile(read_path)
    for sheet in xls.sheet_names:
        if not str(sheet).startswith("25"):
            continue
        raw = _clean_columns(pd.read_excel(read_path, sheet_name=sheet))
        if "t" not in raw.columns and "time" not in {str(c).lower() for c in raw.columns}:
            continue
        frame = normalize_must_sheet(raw, medium="historical_vl3", batch=str(sheet), source_file=CALIBRATION_FILE)
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).sort_values(["medium", "batch", "time_h"]).reset_index(drop=True)


def load_all_review_data() -> pd.DataFrame:
    new_data = load_new_must_data()
    old_data = load_historical_calibration_data()
    frames = [new_data]
    if not old_data.empty:
        frames.append(old_data)
    data = pd.concat(frames, ignore_index=True, sort=False)
    data = data.sort_values(["medium", "batch", "time_h"]).reset_index(drop=True)
    return data


def observation_summary(data: pd.DataFrame) -> pd.DataFrame:
    columns = SECONDARY_PLOT_COLUMNS + AROMA_PLOT_COLUMNS
    rows = []
    for medium, group in data.groupby("medium", sort=True):
        for column in columns:
            if column not in group.columns:
                continue
            values = pd.to_numeric(group[column], errors="coerce")
            rows.append(
                {
                    "medium": medium,
                    "variable": column,
                    "n_obs": int(values.notna().sum()),
                    "n_batches": int(group.loc[values.notna(), "batch"].nunique()),
                    "min": float(values.min()) if values.notna().any() else np.nan,
                    "median": float(values.median()) if values.notna().any() else np.nan,
                    "max": float(values.max()) if values.notna().any() else np.nan,
                }
            )
    return pd.DataFrame(rows)


def batch_variable_summary(data: pd.DataFrame) -> pd.DataFrame:
    columns = SECONDARY_PLOT_COLUMNS + AROMA_PLOT_COLUMNS
    rows = []
    for (medium, batch), group in data.groupby(["medium", "batch"], sort=True):
        row = {
            "medium": medium,
            "batch": batch,
            "n_time": int(group["time_h"].notna().sum()),
            "t_max_h": float(group["time_h"].max()),
        }
        for column in columns:
            if column not in group.columns:
                continue
            values = pd.to_numeric(group[column], errors="coerce")
            row[f"n_{column}"] = int(values.notna().sum())
            row[f"first_{column}"] = float(values.dropna().iloc[0]) if values.notna().any() else np.nan
            row[f"last_{column}"] = float(values.dropna().iloc[-1]) if values.notna().any() else np.nan
            row[f"max_{column}"] = float(values.max()) if values.notna().any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def yan_component_consistency(data: pd.DataFrame) -> pd.DataFrame:
    required = {"YAN_mg_l", "PAN_mg_l", "NH4_mg_l"}
    if not required.issubset(data.columns):
        return pd.DataFrame()
    df = data[["medium", "batch", "time_h", "YAN_mg_l", "PAN_mg_l", "NH4_mg_l"]].copy()
    df["YAN_from_components_mg_l"] = df["PAN_mg_l"] + df["NH4_mg_l"]
    df["YAN_component_residual_mg_l"] = df["YAN_mg_l"] - df["YAN_from_components_mg_l"]
    return df


def plot_secondary_by_batch(data: pd.DataFrame) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    for (medium, batch), group in data.groupby(["medium", "batch"], sort=True):
        group = group.sort_values("time_h")
        present_secondary = [c for c in SECONDARY_PLOT_COLUMNS if c in group.columns and group[c].notna().any()]
        present_aroma = [c for c in AROMA_PLOT_COLUMNS if c in group.columns and group[c].notna().any()]
        if not present_secondary and not present_aroma:
            continue

        fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
        axes = axes.ravel()
        axes[0].plot(group["time_h"], group["G_g_l"], "o-", label="G")
        axes[0].plot(group["time_h"], group["F_g_l"], "o-", label="F")
        axes[0].set_title("Sugars")
        axes[0].legend(fontsize=8)
        axes[1].plot(group["time_h"], group["X_viable_kg_m3"], "o-", label="X viable")
        axes[1].set_title("Biomass")
        axes[1].legend(fontsize=8)
        for col in ["YAN_mg_l", "PAN_mg_l", "NH4_mg_l"]:
            if col in group.columns and group[col].notna().any():
                axes[2].plot(group["time_h"], group[col], "o-", label=col)
        axes[2].set_title("Nitrogen components")
        axes[2].legend(fontsize=8)
        for col in ["pyruvic_acid", "pyruvic_acid_2", "acetaldehyde", "acetic_acid"]:
            if col in group.columns and group[col].notna().any():
                axes[3].plot(group["time_h"], group[col], "o-", label=col)
        axes[3].set_title("Secondary acids/aldehydes")
        axes[3].legend(fontsize=8)
        for col in ["E_g_l", "glycerol_g_l"]:
            if col in group.columns and group[col].notna().any():
                axes[4].plot(group["time_h"], group[col], "o-", label=col)
        if "DO_mg_l" in group.columns and group["DO_mg_l"].notna().any():
            axes[4].plot(group["time_h"], group["DO_mg_l"], "o-", label="DO_mg_l")
        axes[4].set_title("Ethanol and glycerol")
        axes[4].legend(fontsize=8)
        for col in present_aroma:
            axes[5].plot(group["time_h"], group[col], "o-", label=col)
        axes[5].set_title("Liquid aroma measurements")
        if present_aroma:
            axes[5].legend(fontsize=7, ncol=2)
        for ax in axes:
            ax.grid(True, alpha=0.25)
            ax.set_xlabel("time [h]")
        fig.suptitle(f"{medium}/{batch}: secondary metabolites and aromas", y=0.995)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"secondary_{medium}_{batch}.png", dpi=160)
        plt.close(fig)


def plot_medium_overlays(data: pd.DataFrame) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    variables = [c for c in SECONDARY_PLOT_COLUMNS + AROMA_PLOT_COLUMNS if c in data.columns and data[c].notna().any()]
    for variable in variables:
        fig, ax = plt.subplots(figsize=(9, 5))
        for (medium, batch), group in data.groupby(["medium", "batch"], sort=True):
            if variable not in group.columns or not group[variable].notna().any():
                continue
            ax.plot(group["time_h"], group[variable], marker="o", lw=1.2, alpha=0.65, label=f"{medium}/{batch}")
        ax.set_title(variable)
        ax.set_xlabel("time [h]")
        ax.grid(True, alpha=0.25)
        handles, labels = ax.get_legend_handles_labels()
        if len(labels) <= 18:
            ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(PLOT_DIR / f"overlay_{variable}.png", dpi=160)
        plt.close(fig)


def write_modeling_notes() -> None:
    notes = [
        "# Secondary-metabolite modeling notes",
        "",
        "## Data-driven interpretation",
        "",
        "The review plots should be read as a model-structure diagnostic. A state is a good candidate for direct calibration only when it has repeated time trajectories, coverage over the phase where the state changes, and enough perturbation across media or inputs. A state should remain an auxiliary or regularized output when it is sparsely measured, mostly terminal, or confounded with an unmeasured loss process.",
        "",
        "Current evidence supports adding pyruvate as a transient state. Acetaldehyde and acetic acid are also informative, but their dynamics should be coupled to redox and oxygen exposure instead of being fitted as independent Monod products. Liquid aromas should be treated as net liquid outputs unless gas or condensate measurements are available for the same experiment.",
        "",
        "PAN and ammonium should not be collapsed blindly into a single nitrogen pool. A minimal upgrade is to keep two external nitrogen states and let the growth/uptake terms use an effective nitrogen availability:",
        "",
        "$$N_{eff}=w_{NH4}N_{NH4}+w_{PAN}N_{PAN}$$",
        "",
        "with $w_{NH4} \\geq w_{PAN}$ or literature-informed priors. This preserves the fast ammonium depletion and slower PAN behavior visible in the data.",
        "",
        "## Acetaldehyde, acetate, and pyruvate",
        "",
        "The observed trajectories should not be forced into simple Monod-like accumulation. These compounds are intracellular/extracellular exchange nodes and often show transient peaks, lag-phase behavior, and redox-stress dependence.",
        "",
        "Recommended first dynamic layer:",
        "",
        "$$r_S=r_G+r_F$$",
        "",
        "$$\\frac{dPyr}{dt}=k_{Pyr,S}r_S+k_{Pyr,O}g_{O_2}(O_2)X-k_{Pyr,drain}PyrX$$",
        "",
        "$$\\frac{dAcAld}{dt}=k_{Ald,Pyr}PyrX+k_{Ald,S}r_S-k_{Ald,red}AcAldXh_{ana}(O_2)-k_{Ac,Ald}AcAldXg_{O_2}(O_2)$$",
        "",
        "$$\\frac{dAcetate}{dt}=k_{Ac,Ald}AcAldXg_{O_2}(O_2)+k_{Ac,stress}Xh_E(E)-k_{Ac,assim}AcetateXh_N(N_{eff})$$",
        "",
        "For the ODE calibration phase, use parsimonious empirical rates rather than a full biochemical network. Candidate drivers are sugar uptake rate, biomass, ethanol inhibition, nitrogen limitation, and an oxygen exposure factor.",
        "",
        "The practical rule is to estimate only the terms supported by the data. For example, if oxygen is not measured, fix or tightly regularize the oxygen-gated terms and estimate the sugar-driven pyruvate production plus one drain term first.",
        "",
        "## Net liquid aroma outputs",
        "",
        "With liquid-only aroma measurements, the identifiable quantity is the net liquid concentration:",
        "",
        "$$\\frac{dA_i^L}{dt}=r_{i,prod}(z,t)-k_{i,loss}(T,E,CO_2)A_i^L-k_{i,deg}A_i^L$$",
        "",
        "$$y_{i,liq}(t)=A_i^L(t)+\\epsilon_i(t)$$",
        "",
        "If a final condensate measurement exists, it can constrain the time-integrated loss:",
        "",
        "$$A_i^C(t_f)=\\int_0^{t_f} k_{i,loss}(T,E,CO_2)A_i^L(t)dt$$",
        "",
        "Here `r_prod` includes biological production and `k_loss` lumps volatilization, stripping, adsorption, and gas-liquid transfer. Without condensate/gas measurements, production and loss are not separately identifiable. Use literature or UNIFAC/Henry estimates to fix/regularize `k_loss`, and estimate only a small number of production parameters.",
        "",
        "For MPCC use, the robust output is therefore the predicted liquid concentration band, not a uniquely decomposed biological synthesis flux. The decomposition becomes reliable only when liquid, gas/condensate, and CO2-rate information are combined.",
        "",
        "## Oxygen",
        "",
        "The initial oxygen pulse should be represented, even if not fully controlled. A minimal state is dissolved oxygen plus an oxygen-exposure integral:",
        "",
        "$$\\frac{dO_2}{dt}=k_La(O_2^*(T,E)-O_2)-q_{O_2}Xf_S(S)f_N(N_{eff})$$",
        "",
        "$$\\frac{d\\Omega}{dt}=O_2$$",
        "",
        "$$\\mu_{eff}=\\mu_{anaer}(S,N,T,E)\\,f_{lag}(\\Omega)$$",
        "",
        "If `k_La` is not controlled, treat it as experiment-specific or fixed from a short oxygen measurement protocol. The oxygen state should mainly gate lag-phase growth, sterol/lipid limitation, acetaldehyde/acetate behavior, and aroma-side pathways. Do not allow unconstrained oxygen parameters to explain late fermentation unless oxygen is measured.",
        "",
        "A practical measurement protocol is sparse but early: dissolved oxygen at inoculation, 2 h, 6 h, 12 h, and 24 h for a subset of fermentations. That is enough to estimate or validate an initial oxygen exposure lump, but not enough to identify a detailed aeration/transfer model.",
    ]
    (RESULTS_DIR / "secondary_metabolite_modeling_notes.md").write_text("\n".join(notes), encoding="utf-8")


def write_notebook() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Secondary metabolites and liquid aroma data review\n\n"
            "This notebook visualizes the data now available for nitrogen fractions, pyruvic acid, acetaldehyde, acetic acid, and liquid aroma measurements. "
            "The purpose is model-structure selection, not parameter estimation."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "ROOT = Path.cwd()\n"
            "if (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'shared/results/secondary_metabolite_data_review'\n"
            "summary = pd.read_csv(RESULTS / 'observation_summary.csv')\n"
            "batch_summary = pd.read_csv(RESULTS / 'batch_variable_summary.csv')\n"
            "summary"
        ),
        nbformat.v4.new_markdown_cell(
            "## Observation availability\n\n"
            "Use this table to decide which variables can support dynamic calibration and which should remain diagnostic outputs or regularized states."
        ),
        nbformat.v4.new_code_cell(
            "display(summary.sort_values(['medium', 'variable']))\n"
            "availability = summary.pivot(index='variable', columns='medium', values='n_obs').fillna(0).astype(int)\n"
            "display(availability)\n"
            "display(batch_summary.head(20))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Immediate interpretation\n\n"
            "- Pyruvate has enough repeated trajectories to be considered a transient state.\n"
            "- Acetaldehyde is available in natural and historical batches, but not synthetic, so it is useful for validation and regularized fitting rather than unconstrained expansion.\n"
            "- Acetic acid is currently natural-only, so it should enter as a low-dimensional auxiliary state.\n"
            "- Liquid aromas are currently historical-only in this merged review; new natural aroma data should be appended before model-based aroma DOE.\n"
            "- PAN and ammonium have distinct behavior and should be retained as separate nitrogen pools in the next model version."
        ),
        nbformat.v4.new_markdown_cell(
            "## YAN component consistency\n\n"
            "The loader keeps total `YAN`, `PAN`, and `NH4` separately. The residual `YAN - (PAN + NH4)` is useful for detecting unit, assay, or preprocessing differences."
        ),
        nbformat.v4.new_code_cell(
            "yan = pd.read_csv(RESULTS / 'yan_component_consistency.csv')\n"
            "yan.dropna(subset=['YAN_component_residual_mg_l']).groupby('medium')['YAN_component_residual_mg_l'].describe()"
        ),
        nbformat.v4.new_markdown_cell(
            "## Per-batch plots\n\n"
            "The following figures show sugars, biomass, nitrogen fractions, secondary metabolites, ethanol/glycerol, and liquid aromas by batch."
        ),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_*.png'))[:12]:\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Variable overlays\n\n"
            "Overlay plots help identify non-Monod behavior, transient peaks, and medium-specific differences."
        ),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('overlay_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Modeling notes\n\n"
            "See `secondary_metabolite_modeling_notes.md` for the proposed ODE treatment of pyruvate, acetaldehyde, acetate, net liquid aroma outputs, and oxygen."
        ),
        nbformat.v4.new_code_cell(
            "print((RESULTS / 'secondary_metabolite_modeling_notes.md').read_text(encoding='utf-8'))"
        ),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> int:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_all_review_data()
    data.to_csv(RESULTS_DIR / "secondary_metabolite_long.csv", index=False)
    observation_summary(data).to_csv(RESULTS_DIR / "observation_summary.csv", index=False)
    batch_variable_summary(data).to_csv(RESULTS_DIR / "batch_variable_summary.csv", index=False)
    yan_component_consistency(data).to_csv(RESULTS_DIR / "yan_component_consistency.csv", index=False)
    plot_secondary_by_batch(data)
    plot_medium_overlays(data)
    write_modeling_notes()
    write_notebook()
    print(f"[done] data rows={len(data)} results={RESULTS_DIR}")
    print(f"[done] notebook={NOTEBOOK_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
