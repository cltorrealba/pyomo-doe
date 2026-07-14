from __future__ import annotations

import nbformat
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
NOTEBOOK_PATH = SCRIPT_DIR / "fermentation_final_operational_doe_volume_constrained.ipynb"


def md(text: str):
    return nbformat.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbformat.v4.new_code_cell(text.strip())


def main() -> None:
    nb = nbformat.v4.new_notebook()
    nb["cells"] = [
        md(
            r"""
# Volume-Constrained Final Operational DOE Campaign

This notebook documents the executable campaign selected after adding the reactor-volume constraint.

Use this notebook as the visual companion to the CSV protocol:

- `full_liquid_sample`: 50 mL. Ethanol and aromas are sampled together.
- `small_liquid_sample`: 8 mL. No ethanol/aromas; intended for biomass, sugars, YAN/PAN/NH4, glycerol, pyruvate, acetaldehyde, and acetate.
- DO is measured, not manipulated.
- CO2 is assumed online and therefore not represented as a destructive liquid sample.
- Terminal condensate retrieval is a final trap sample, not a repeated online aroma measurement.

The active design is the volume-constrained campaign in `results/final_operational_doe_volume_constrained`.
The older unconstrained DOE is retained only as a reference and should not be used as the lab protocol.
            """
        ),
        code(
            r"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import display, Markdown

ROOT = Path.cwd()
if (ROOT / "fermentation_model").exists():
    FM = ROOT / "fermentation_model"
else:
    FM = ROOT

if str(FM) not in sys.path:
    sys.path.insert(0, str(FM))

import run_final_operational_doe_volume_constrained as vc
import run_final_operational_doe_v2 as final
import run_new_must_glycerol_estimability_doe as base
import run_secondary_joint_campaign_doe as joint

RESULTS = FM / "results" / "final_operational_doe_volume_constrained"
PLOT_DIR = FM / "results" / "vc_doe_nb_plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

protocol = pd.read_csv(RESULTS / "campaign_protocol.csv")
schedule = pd.read_csv(RESULTS / "operational_schedule.csv", parse_dates=["datetime"])
volume = pd.read_csv(RESULTS / "volume_audit.csv")
stock = pd.read_csv(RESULTS / "stock_concentration_assumptions.csv")
policy_summary = pd.read_csv(RESULTS / "policy_summary.csv")
post_estimability = pd.read_csv(RESULTS / "post_campaign_estimability.csv")

theta = final.load_theta_final()
data = base.load_normalized_data()
designs = final.candidate_library(data)
POLICY = str(policy_summary.iloc[0]["policy"])
POLICY
            """
        ),
        md(
            r"""
## Design Selection Summary

The selected policy is compared against more conservative full-sample schedules. Scores are computed from the extended FIM under the volume-constrained observation model, so they should only be compared within this table.
            """
        ),
        code(
            r"""
cols = [
    "policy",
    "campaign_logdet",
    "campaign_trace_inv",
    "new_param_mean_var_reduction",
    "new_param_worst_var_reduction",
]
display(policy_summary[cols])
            """
        ),
        md(
            r"""
## Campaign Protocol

One row is one fermentation. Lots are groups of three parallel reactors.
            """
        ),
        code(
            r"""
display(protocol)
            """
        ),
        md(
            r"""
## Volume And Stock Audit

The design keeps final reactor volume above the 1100 mL safety target after destructive samples and stock additions.

Stock pulse volume is computed as:

```text
V_stock = target_dose * current_reactor_volume / stock_concentration
```
            """
        ),
        code(
            r"""
display(volume)
display(stock)
            """
        ),
        md(
            r"""
## Daily Workload

This table is useful for lab staffing and analytical capacity. It counts events at the lot level, so three parallel reactors can produce up to three simultaneous actions.
            """
        ),
        code(
            r"""
work = schedule.copy()
work["date"] = work["datetime"].dt.date
daily = work.groupby(["lot", "date", "event"]).size().unstack(fill_value=0).reset_index()
display(daily)
            """
        ),
        md(
            r"""
## Complete Operational Schedule

The key event types are:

- `full_liquid_sample`: 50 mL, ethanol + aromas + full panel.
- `small_liquid_sample`: 8 mL, no ethanol/aromas.
- `manual_pulse`: stock addition.
- `DO_measurement`: oxygen observation, no oxygen control.
- `terminal_condensate_retrieval`: final condensate trap sample.
            """
        ),
        code(
            r"""
display(schedule)
            """
        ),
        md(
            r"""
## Per-Fermentation DOE Plots

Each plot shows the actual DOE inputs and observation schedule for one selected fermentation:

1. Temperature setpoint.
2. Manual pulses.
3. Full/small/DO/terminal sampling events.
4. Reactor volume trajectory.
5. Simulated core states.
6. Simulated secondary states and CO2.
7. Simulated liquid aromas plus terminal condensate prediction.
            """
        ),
        code(
            r"""
def _event_times(group, event):
    return group.loc[group["event"].eq(event), "relative_time_h"].astype(float).to_numpy()


def _plot_one(record):
    name = str(record["candidate"])
    order = int(record["campaign_order"])
    design = designs[name]
    group = schedule[schedule["candidate"].eq(name)].sort_values("relative_time_h").copy()
    sim = vc.simulate(design, theta, POLICY)
    t_grid = np.linspace(0.0, float(design.horizon_h), 400)

    fig, axes = plt.subplots(7, 1, figsize=(12, 16), sharex=True)

    axes[0].step(t_grid, [base.temperature_at(design, t) for t in t_grid], where="post", color="tab:red")
    axes[0].set_ylabel("T [C]")
    axes[0].set_title(f"{order}. {name}")

    pulse_rows = group[group["event"].eq("manual_pulse")]
    if pulse_rows.empty:
        axes[1].text(0.02, 0.5, "no manual pulses", transform=axes[1].transAxes, va="center")
    else:
        for channel, rows in pulse_rows.groupby("channel"):
            axes[1].scatter(rows["relative_time_h"], rows["amount"], label=channel, s=45)
            for _, row in rows.iterrows():
                axes[1].annotate(
                    f"{row['channel']} {float(row['amount']):g}",
                    (float(row["relative_time_h"]), float(row["amount"])),
                    textcoords="offset points",
                    xytext=(5, 5),
                    fontsize=8,
                )
        axes[1].legend(fontsize=8, ncol=5)
    axes[1].set_ylabel("pulse amount")

    markers = [
        ("full_liquid_sample", "tab:blue", 1.0, "full 50 mL"),
        ("small_liquid_sample", "tab:orange", 0.7, "small 8 mL"),
        ("DO_measurement", "tab:green", 0.4, "DO"),
        ("terminal_condensate_retrieval", "tab:purple", 1.25, "condensate"),
    ]
    for event, color, ymax, label in markers:
        times = _event_times(group, event)
        if len(times):
            axes[2].vlines(times, 0, ymax, colors=color, lw=1.8, label=label)
    axes[2].set_ylim(0, 1.45)
    axes[2].set_yticks([])
    axes[2].set_ylabel("samples")
    axes[2].legend(fontsize=8, ncol=4)

    axes[3].step(group["relative_time_h"], group["volume_after_ml"], where="post", color="black")
    axes[3].axhline(vc.MIN_FINAL_VOLUME_ML, color="tab:red", ls="--", label="min volume")
    axes[3].set_ylabel("volume [mL]")
    axes[3].legend(fontsize=8)

    if sim is not None:
        axes[4].plot(sim.core.index, sim.core["G"], label="G")
        axes[4].plot(sim.core.index, sim.core["F"], label="F")
        axes[4].plot(sim.core.index, sim.core["E"], label="E")
        axes[4].plot(sim.core.index, sim.core["Gly"], label="Gly")
        axes[4].plot(sim.core.index, sim.core["N"] * 1000.0, label="YAN mg/L")
        axes[4].set_ylabel("core")
        axes[4].legend(fontsize=8, ncol=5)

        for state in ("Pyr", "AcAld", "Acetate", "O2", "CO2"):
            axes[5].plot(sim.secondary.index, sim.secondary[state], label=state)
        axes[5].set_ylabel("secondary")
        axes[5].legend(fontsize=8, ncol=5)

        for species in joint.AROMA_SPECIES:
            axes[6].plot(sim.aromas_liq.index, sim.aromas_liq[species], label=f"{species} liq")
        for idx, species in enumerate(joint.AROMA_SPECIES):
            axes[6].scatter(
                [float(design.horizon_h)],
                [float(sim.aromas_cond[species])],
                marker="x",
                s=60,
                label=f"{species} condensate",
            )
        axes[6].set_ylabel("aromas")
        axes[6].legend(fontsize=8, ncol=3)
    else:
        axes[4].text(0.02, 0.5, "simulation failed", transform=axes[4].transAxes, va="center")

    for ax in axes:
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("process time [h]")
    fig.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOT_DIR / f"doe_{order:02d}.png"
    fig.savefig(out, dpi=170)
    return fig


for _, record in protocol.sort_values("campaign_order").iterrows():
    fig = _plot_one(record)
    plt.show()
            """
        ),
        md(
            r"""
## Expected Post-Campaign Estimability

This is the local FIM-based expectation after adding the selected volume-constrained campaign to the current prior. Interpret it as design guidance, not as a guarantee: after each lot, the model should be recalibrated and the next lots should be reviewed.
            """
        ),
        code(
            r"""
display(post_estimability)
            """
        ),
        md(
            r"""
## Files Written By This Notebook

Per-fermentation overview plots were saved under:

`results/vc_doe_nb_plots/`
            """
        ),
        code(
            r"""
for path in sorted(PLOT_DIR.glob("doe_*.png")):
    print(path)
            """
        ),
    ]
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
    nbformat.write(nb, NOTEBOOK_PATH)
    print(NOTEBOOK_PATH)


if __name__ == "__main__":
    main()
