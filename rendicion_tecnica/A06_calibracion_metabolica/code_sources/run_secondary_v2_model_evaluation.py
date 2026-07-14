from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_secondary_joint_campaign_doe as v1


RESULTS_DIR = SCRIPT_DIR / "results" / "secondary_v2_model_evaluation"
NOTEBOOK_PATH = SCRIPT_DIR / "fermentation_secondary_v2_model_evaluation.ipynb"

V2_PARAMETERS = (
    "kPyrS_N",
    "kPyrS_stat",
    "kPyrO2",
    "kPyrDrain",
    "kAldPyr",
    "kAldS_N",
    "kAldS_stat",
    "kAldO2",
    "kAldRed",
    "kAcAld",
    "kAcStress",
    "kAcAssim",
    "qO2",
    "kLaO2",
    "O2sat",
)
V2_O2_PARAMETERS = ("qO2", "kLaO2", "O2sat")
V2_CHEM_PARAMETERS = tuple(name for name in V2_PARAMETERS if name not in V2_O2_PARAMETERS)
V2_REDUCED_PARAMETERS = (
    "kPyrS_N",
    "kPyrO2",
    "kPyrDrain",
    "kAldS_N",
    "kAldRed",
    "kAcAld",
    "kAcStress",
)

V2_BOUNDS = {
    "kPyrS_N": (1e-4, 10.0),
    "kPyrS_stat": (1e-4, 10.0),
    "kPyrO2": (1e-4, 10.0),
    "kPyrDrain": (1e-5, 2.0),
    "kAldPyr": (1e-5, 2.0),
    "kAldS_N": (1e-4, 10.0),
    "kAldS_stat": (1e-4, 10.0),
    "kAldO2": (1e-5, 2.0),
    "kAldRed": (1e-5, 2.0),
    "kAcAld": (1e-5, 2.0),
    "kAcStress": (1e-6, 0.5),
    "kAcAssim": (1e-6, 2.0),
    "qO2": (1e-4, 5.0),
    "kLaO2": (1e-5, 0.2),
    "O2sat": (0.5, 8.0),
}

SIGMA = dict(v1.SIGMA)
N_PHASE_HALF_KG_M3 = v1.N_PHASE_HALF_KG_M3
O2_HALF_MG_L = v1.O2_HALF_MG_L


def load_theta_v1() -> dict[str, float]:
    theta = v1.load_reference_theta()
    path = v1.RESULTS_DIR / "theta_secondary_joint.csv"
    if path.exists():
        loaded = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
        theta.update({str(k): float(val) for k, val in loaded.items() if np.isfinite(float(val))})
    return v1.clip_all(theta)


def default_theta_v2(theta_v1: dict[str, float]) -> dict[str, float]:
    return {
        "kPyrS_N": max(float(theta_v1.get("kPyrS", 0.4)), 1e-4),
        "kPyrS_stat": max(float(theta_v1.get("kPyrS", 0.4)) * 0.55, 1e-4),
        "kPyrO2": 0.20,
        "kPyrDrain": max(float(theta_v1.get("kPyrD", 0.008)), 1e-5),
        "kAldPyr": max(float(theta_v1.get("kAldPyr", 0.001)), 1e-5),
        "kAldS_N": max(float(theta_v1.get("kAldS", 0.8)) * 1.2, 1e-4),
        "kAldS_stat": max(float(theta_v1.get("kAldS", 0.8)) * 0.45, 1e-4),
        "kAldO2": 0.020,
        "kAldRed": max(float(theta_v1.get("kAldRed", 0.001)), 1e-5),
        "kAcAld": max(float(theta_v1.get("kAcAld", 0.002)), 1e-5),
        "kAcStress": max(float(theta_v1.get("kAcStress", 0.002)), 1e-6),
        "kAcAssim": max(float(theta_v1.get("kAcAssim", 1e-5)), 1e-6),
        "qO2": max(float(theta_v1.get("qO2", 0.04)), 1e-4),
        "kLaO2": max(float(theta_v1.get("kLaO2", 0.0015)), 1e-5),
        "O2sat": min(max(float(theta_v1.get("O2sat", 2.7)), 0.5), 8.0),
    }


def clip_v2(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in V2_BOUNDS.items():
        out[name] = float(np.clip(float(out[name]), lb, ub))
    return out


def log_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    return np.array([math.log(float(theta[name])) for name in parameters], dtype=float)


def theta_from_log(x: np.ndarray, base: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        theta[name] = float(math.exp(float(value)))
    return clip_v2(theta)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    low = []
    high = []
    for name in parameters:
        lb, ub = V2_BOUNDS[name]
        low.append(math.log(lb))
        high.append(math.log(ub))
    return np.asarray(low, dtype=float), np.asarray(high, dtype=float)


def secondary_rhs_v2(t: float, y: np.ndarray, theta: dict[str, float], batch, core: pd.DataFrame) -> list[float]:
    pyr, ald, acetate, o2, co2 = [max(float(v), 0.0) for v in y]
    r = v1._core_rates(theta, batch, core, t)
    x = max(float(r["X"]), 0.0)
    n_eff = max(float(r["N"]), 0.0)
    h_n = n_eff / (n_eff + N_PHASE_HALF_KG_M3)
    h_stat = 1.0 - h_n
    g_o2 = o2 / (o2 + O2_HALF_MG_L)
    h_ana = O2_HALF_MG_L / (o2 + O2_HALF_MG_L)
    h_e = r["E"] / (r["E"] + 60.0)
    sugar = max(float(r["sugar_uptake"]), 0.0)

    pyr_prod = (theta["kPyrS_N"] * h_n + theta["kPyrS_stat"] * h_stat) * sugar + theta["kPyrO2"] * g_o2 * x
    pyr_drain = theta["kPyrDrain"] * pyr * x * (0.25 + h_stat + 0.5 * h_e)
    d_pyr = pyr_prod - pyr_drain

    ald_prod = (
        theta["kAldPyr"] * pyr * x
        + (theta["kAldS_N"] * h_n + theta["kAldS_stat"] * h_stat) * sugar
        + theta["kAldO2"] * g_o2 * x
    )
    ald_sink = theta["kAldRed"] * ald * x * (h_ana + 0.25 * h_stat) + theta["kAcAld"] * ald * x * g_o2
    d_ald = ald_prod - ald_sink

    d_ac = theta["kAcAld"] * ald * x * g_o2 / 1000.0 + theta["kAcStress"] * x * h_e - theta["kAcAssim"] * acetate * x * h_n
    d_o2 = theta["kLaO2"] * (theta["O2sat"] - o2) - theta["qO2"] * x * o2 / (o2 + O2_HALF_MG_L)
    d_co2 = v1.CO2_G_PER_G_ETHANOL * r["ethanol_prod"]
    return [d_pyr, d_ald, d_ac, d_o2, d_co2]


def integrate_secondary_v2(batch, theta: dict[str, float], core: pd.DataFrame) -> pd.DataFrame:
    time = core.index.to_numpy(dtype=float)
    y = np.array(
        [
            float(batch.initials.get("Pyr", 0.0)),
            float(batch.initials.get("AcAld", 0.0)),
            float(batch.initials.get("Acetate", 0.0)),
            float(batch.initials.get("O2", 6.5)),
            float(batch.initials.get("CO2", 0.0)),
        ],
        dtype=float,
    )
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(secondary_rhs_v2(float(time[idx - 1]), y, theta, batch, core), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    return pd.DataFrame(rows, index=time, columns=("Pyr", "AcAld", "Acetate", "O2", "CO2"))


def residual_v2(theta: dict[str, float], batches: list, core_cache: dict[str, pd.DataFrame]) -> np.ndarray:
    residuals: list[np.ndarray] = []
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        sec = integrate_secondary_v2(batch, theta, core)
        for state in ("Pyr", "AcAld", "Acetate", "O2"):
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sec.loc[batch.time, state].to_numpy(dtype=float)
            residuals.append((pred[mask] - obs[mask]) / SIGMA[state])
    if not residuals:
        return np.array([], dtype=float)
    return np.concatenate(residuals)


def residual_v1(theta: dict[str, float], batches: list, core_cache: dict[str, pd.DataFrame]) -> np.ndarray:
    return v1.current_residual_vector(theta, batches, parameter_mode="secondary", core_cache=core_cache)


def fit_v2(
    theta0: dict[str, float],
    batches: list,
    core_cache: dict[str, pd.DataFrame],
    max_nfev: int,
    fit_parameters: tuple[str, ...],
    model_label: str,
) -> tuple[dict[str, float], pd.DataFrame]:
    theta0 = clip_v2(theta0)
    x0 = log_vector(theta0, fit_parameters)
    lb, ub = log_bounds(fit_parameters)
    x0 = np.clip(x0, lb + 1e-9, ub - 1e-9)

    prior_sigma = {
        "O2sat": 0.35,
        "kLaO2": 1.0,
        "kAcAssim": 1.5,
        "kAldRed": 1.5,
    }
    prior_sigma = {name: sigma for name, sigma in prior_sigma.items() if name in fit_parameters}

    def fun(x: np.ndarray) -> np.ndarray:
        theta = theta_from_log(x, theta0, fit_parameters)
        res = [residual_v2(theta, batches, core_cache)]
        for name, sigma in prior_sigma.items():
            res.append(np.array([math.log(theta[name] / theta0[name]) / sigma], dtype=float))
        return np.concatenate(res)

    start = fun(x0)
    result = least_squares(
        fun,
        x0,
        bounds=(lb, ub),
        method="trf",
        x_scale="jac",
        loss="soft_l1",
        f_scale=2.0,
        max_nfev=int(max_nfev),
        ftol=1e-5,
        xtol=1e-5,
        gtol=1e-5,
    )
    theta_hat = theta_from_log(result.x, theta0, fit_parameters)
    end = fun(result.x)
    summary = pd.DataFrame(
        [
            {
                "model": model_label,
                "success": bool(result.success),
                "status": int(result.status),
                "message": str(result.message),
                "nfev": int(result.nfev),
                "initial_wsse": float(np.dot(start, start)),
                "final_wsse": float(np.dot(end, end)),
                "n_residuals": int(len(end)),
                "wsse_per_residual": float(np.dot(end, end) / max(len(end), 1)),
            }
        ]
    )
    return theta_hat, summary


def per_state_metrics(model_name: str, theta: dict[str, float], batches: list, core_cache: dict[str, pd.DataFrame], model: str) -> pd.DataFrame:
    rows = []
    for state in ("Pyr", "AcAld", "Acetate", "O2"):
        residuals = []
        raw_errors = []
        raw_obs = []
        for batch in batches:
            core = core_cache.get(batch.label)
            if core is None:
                continue
            if model == "v1":
                sec = v1.integrate_secondary_euler(batch, theta, core)
            else:
                sec = integrate_secondary_v2(batch, theta, core)
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if not mask.any():
                continue
            pred = sec.loc[batch.time, state].to_numpy(dtype=float)
            error = pred[mask] - obs[mask]
            residuals.append(error / SIGMA[state])
            raw_errors.append(error)
            raw_obs.append(obs[mask])
        if not residuals:
            continue
        res = np.concatenate(residuals)
        err = np.concatenate(raw_errors)
        obs = np.concatenate(raw_obs)
        rows.append(
            {
                "model": model_name,
                "state": state,
                "n_obs": int(len(err)),
                "wsse": float(np.dot(res, res)),
                "wsse_per_obs": float(np.dot(res, res) / max(len(err), 1)),
                "rmse": float(np.sqrt(np.mean(err * err))),
                "mae": float(np.mean(np.abs(err))),
                "obs_median": float(np.nanmedian(obs)),
                "obs_range": float(np.nanmax(obs) - np.nanmin(obs)),
            }
        )
    return pd.DataFrame(rows)


def finite_difference_jacobian(
    theta: dict[str, float],
    batches: list,
    core_cache: dict[str, pd.DataFrame],
    parameters: tuple[str, ...],
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    base_res = residual_v2(theta, batches, core_cache)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        lb, ub = V2_BOUNDS[name]
        theta_plus[name] = float(np.clip(theta[name] * math.exp(step), lb, ub))
        theta_minus[name] = float(np.clip(theta[name] * math.exp(-step), lb, ub))
        r_plus = residual_v2(theta_plus, batches, core_cache)
        r_minus = residual_v2(theta_minus, batches, core_cache)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
    return np.column_stack(cols), base_res


def stable_inverse(fim: np.ndarray, ridge_fraction: float = 1e-9) -> np.ndarray:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    scale = max(float(np.trace(fim)) / max(fim.shape[0], 1), 1.0)
    return np.linalg.pinv(fim + ridge_fraction * scale * np.eye(fim.shape[0]))


def estimability_v2(fim: np.ndarray, theta: dict[str, float], parameters: tuple[str, ...]) -> pd.DataFrame:
    cov = stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[idx, idx], 0.0)))
        value = float(theta[name])
        lb, ub = V2_BOUNDS[name]
        active = value <= lb * 1.01 or value >= ub / 1.01
        if std_log <= 0.35 and not active:
            cls = "well_estimated"
        elif std_log <= 0.75 and not active:
            cls = "moderate"
        elif std_log <= 1.25 and not active:
            cls = "weak_but_actionable"
        else:
            cls = "weak_or_confounded"
        rows.append(
            {
                "parameter": name,
                "theta": value,
                "std_log_approx": std_log,
                "approx_95_multiplier": float(math.exp(1.96 * min(std_log, 20.0))),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return pd.DataFrame(rows)


def eigen_diagnostics(fim: np.ndarray, parameters: tuple[str, ...]) -> tuple[pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eigvals, eigvecs = np.linalg.eigh(fim)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    max_eig = float(np.max(eigvals)) if eigvals.size else np.nan
    spectrum = pd.DataFrame(
        {
            "direction": np.arange(1, len(eigvals) + 1),
            "eigenvalue": eigvals,
            "relative_eigenvalue": eigvals / max_eig if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        }
    )
    rows = []
    for idx in range(min(6, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        top = np.argsort(np.abs(vec))[::-1][:6]
        rows.append(
            {
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs(vec[i]):.3f}" for i in top),
            }
        )
    return spectrum, pd.DataFrame(rows)


def plot_fit_comparison(theta_v1: dict[str, float], theta_v2: dict[str, float], batches: list, core_cache: dict[str, pd.DataFrame]) -> None:
    plot_dir = RESULTS_DIR / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for batch in batches:
        states = [
            state
            for state in ("Pyr", "AcAld", "Acetate", "O2")
            if np.isfinite(np.asarray(batch.observations.get(state, []), dtype=float)).any()
        ]
        if not states:
            continue
        core = core_cache.get(batch.label)
        if core is None:
            continue
        sec1 = v1.integrate_secondary_euler(batch, theta_v1, core)
        sec2 = integrate_secondary_v2(batch, theta_v2, core)
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
        axes = axes.ravel()
        for ax, state in zip(axes, ("Pyr", "AcAld", "Acetate", "O2")):
            ax.plot(sec1.index, sec1[state], "--", lw=1.4, label="v1")
            ax.plot(sec2.index, sec2[state], "-", lw=1.7, label="v2")
            obs = np.asarray(batch.observations.get(state, np.full_like(batch.time, np.nan, dtype=float)), dtype=float)
            mask = np.isfinite(obs)
            if mask.any():
                ax.scatter(batch.time[mask], obs[mask], s=24, color="black", label="obs")
            ax.set_title(state)
            ax.set_xlabel("time [h]")
            ax.grid(True, alpha=0.25)
            ax.legend(fontsize=8)
        fig.suptitle(f"Secondary fit comparison: {batch.medium}/{batch.batch}")
        fig.tight_layout()
        fig.savefig(plot_dir / f"secondary_v2_fit_{batch.medium}_{batch.batch}.png", dpi=170)
        plt.close(fig)


def write_report(fit_table: pd.DataFrame, state_metrics: pd.DataFrame, estim: pd.DataFrame, weak: pd.DataFrame) -> None:
    lines = [
        "# Secondary v2 model evaluation",
        "",
        "## Purpose",
        "",
        "This run tests whether the poor current fits for pyruvate, acetaldehyde, acetate, and DO are mainly due to the secondary-state structure. The core fermentation simulation is kept fixed and only the secondary layer is changed.",
        "",
        "## v2 structure",
        "",
        "The v2 layer adds nitrogen-phase and oxygen gates to the apparent production and drainage terms:",
        "",
        "$$h_N=\\frac{N}{N+K_N},\\qquad h_{stat}=1-h_N,$$",
        "",
        "$$g_{O_2}=\\frac{O_2}{O_2+K_{O_2}},\\qquad h_{ana}=\\frac{K_{O_2}}{O_2+K_{O_2}}.$$",
        "",
        "Pyruvate production is split between N-rich and stationary phases plus an early oxygen-linked source. Acetaldehyde production is split similarly and includes pyruvate and oxygen-linked contributions. Acetate keeps the acetaldehyde oxidation, stress, and assimilation terms. DO keeps passive transfer plus cellular consumption.",
        "",
        "## Fit summary",
        "",
        fit_table.to_markdown(index=False),
        "",
        "## Per-state comparison",
        "",
        state_metrics.to_markdown(index=False),
        "",
        "## v2 weak/confounded parameters",
        "",
        weak.to_markdown(index=False) if not weak.empty else "No weak/confounded v2 parameters by the local FIM criterion.",
        "",
        "## Full v2 estimability",
        "",
        estim.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- v2 should be kept only if it improves state-wise errors enough to justify the extra parameters.",
        "- Parameters that remain weak after v2 should be fixed or strongly regularized before using the model in DOE/MPCC.",
        "- If v2 improves pyruvate/acetaldehyde but leaves acetate assimilation weak, this supports fixing acetate assimilation rather than forcing DOE to identify it.",
    ]
    (RESULTS_DIR / "secondary_v2_model_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_notebook() -> None:
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Secondary v2 model evaluation\n\n"
            "This notebook compares the original secondary layer against a phase/O2-gated v2 layer while keeping the calibrated core fermentation model fixed."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "ROOT = Path.cwd()\n"
            "if not (ROOT / 'results').exists() and (ROOT / 'fermentation_model').exists():\n"
            "    ROOT = ROOT / 'fermentation_model'\n"
            "RESULTS = ROOT / 'results/secondary_v2_model_evaluation'\n"
            "display(pd.read_csv(RESULTS / 'fit_comparison.csv'))\n"
            "display(pd.read_csv(RESULTS / 'state_fit_metrics.csv'))\n"
            "display(pd.read_csv(RESULTS / 'v2_reduced_o2fixed_estimability.csv'))"
        ),
        nbformat.v4.new_markdown_cell("## Weak FIM directions"),
        nbformat.v4.new_code_cell(
            "display(pd.read_csv(RESULTS / 'v2_reduced_o2fixed_eigen_spectrum.csv').head(8))\n"
            "display(pd.read_csv(RESULTS / 'v2_reduced_o2fixed_weak_loadings.csv'))"
        ),
        nbformat.v4.new_markdown_cell("## Fit comparison plots"),
        nbformat.v4.new_code_cell(
            "for png in sorted((RESULTS / 'plots').glob('secondary_v2_fit_*.png')):\n"
            "    print(png.name)\n"
            "    display(Image(filename=str(png)))"
        ),
        nbformat.v4.new_markdown_cell("## Report"),
        nbformat.v4.new_code_cell("print((RESULTS / 'secondary_v2_model_report.md').read_text(encoding='utf-8'))"),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate a phase/O2-gated secondary-state model.")
    parser.add_argument("--max-nfev", type=int, default=100)
    parser.add_argument("--step", type=float, default=1e-2)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)

    data = v1.load_current_review_data()
    batches = v1.make_secondary_batches(data)
    theta_v1 = load_theta_v1()
    core_cache = v1.precompute_core_cache(batches, theta_v1)
    theta0_v2 = {**theta_v1, **default_theta_v2(theta_v1)}

    res_v1 = residual_v1(theta_v1, batches, core_cache)
    theta_v2_free, fit_v2_free_summary = fit_v2(
        theta0_v2,
        batches,
        core_cache,
        args.max_nfev,
        V2_PARAMETERS,
        "secondary_v2_phase_o2free",
    )
    theta_v2_fixed, fit_v2_fixed_summary = fit_v2(
        theta0_v2,
        batches,
        core_cache,
        args.max_nfev,
        V2_CHEM_PARAMETERS,
        "secondary_v2_phase_o2fixed",
    )
    theta0_reduced = dict(theta0_v2)
    theta0_reduced.update(
        {
            "kPyrS_stat": V2_BOUNDS["kPyrS_stat"][0],
            "kAldPyr": V2_BOUNDS["kAldPyr"][0],
            "kAldS_stat": V2_BOUNDS["kAldS_stat"][0],
            "kAldO2": V2_BOUNDS["kAldO2"][0],
            "kAcAssim": V2_BOUNDS["kAcAssim"][0],
        }
    )
    theta_v2_reduced, fit_v2_reduced_summary = fit_v2(
        theta0_reduced,
        batches,
        core_cache,
        args.max_nfev,
        V2_REDUCED_PARAMETERS,
        "secondary_v2_reduced_o2fixed",
    )
    fit_v1_summary = pd.DataFrame(
        [
            {
                "model": "secondary_v1_current",
                "success": True,
                "status": np.nan,
                "message": "loaded from current secondary_joint fit",
                "nfev": np.nan,
                "initial_wsse": np.nan,
                "final_wsse": float(np.dot(res_v1, res_v1)),
                "n_residuals": int(len(res_v1)),
                "wsse_per_residual": float(np.dot(res_v1, res_v1) / max(len(res_v1), 1)),
            }
        ]
    )
    fit_table = pd.concat(
        [fit_v1_summary, fit_v2_free_summary, fit_v2_fixed_summary, fit_v2_reduced_summary],
        ignore_index=True,
    )

    pd.Series(theta_v2_free, name="theta_secondary_v2_o2free").to_csv(RESULTS_DIR / "theta_secondary_v2_o2free.csv")
    pd.Series(theta_v2_fixed, name="theta_secondary_v2_o2fixed").to_csv(RESULTS_DIR / "theta_secondary_v2_o2fixed.csv")
    pd.Series(theta_v2_reduced, name="theta_secondary_v2_reduced_o2fixed").to_csv(
        RESULTS_DIR / "theta_secondary_v2_reduced_o2fixed.csv"
    )
    fit_table.to_csv(RESULTS_DIR / "fit_comparison.csv", index=False)

    state_metrics = pd.concat(
        [
            per_state_metrics("secondary_v1_current", theta_v1, batches, core_cache, "v1"),
            per_state_metrics("secondary_v2_phase_o2free", theta_v2_free, batches, core_cache, "v2"),
            per_state_metrics("secondary_v2_phase_o2fixed", theta_v2_fixed, batches, core_cache, "v2"),
            per_state_metrics("secondary_v2_reduced_o2fixed", theta_v2_reduced, batches, core_cache, "v2"),
        ],
        ignore_index=True,
    )
    state_metrics.to_csv(RESULTS_DIR / "state_fit_metrics.csv", index=False)

    jac, base_res = finite_difference_jacobian(theta_v2_reduced, batches, core_cache, V2_REDUCED_PARAMETERS, args.step)
    fim = jac.T @ jac
    fim = 0.5 * (fim + fim.T)
    pd.DataFrame(fim, index=V2_REDUCED_PARAMETERS, columns=V2_REDUCED_PARAMETERS).to_csv(
        RESULTS_DIR / "v2_reduced_o2fixed_fim.csv"
    )
    estim = estimability_v2(fim, theta_v2_reduced, V2_REDUCED_PARAMETERS)
    estim.to_csv(RESULTS_DIR / "v2_reduced_o2fixed_estimability.csv", index=False)
    spectrum, loadings = eigen_diagnostics(fim, V2_REDUCED_PARAMETERS)
    spectrum.to_csv(RESULTS_DIR / "v2_reduced_o2fixed_eigen_spectrum.csv", index=False)
    loadings.to_csv(RESULTS_DIR / "v2_reduced_o2fixed_weak_loadings.csv", index=False)

    plot_fit_comparison(theta_v1, theta_v2_reduced, batches, core_cache)
    weak = estim[estim["classification"].eq("weak_or_confounded")].copy()
    write_report(fit_table, state_metrics, estim, weak)
    write_notebook()
    metadata = {
        "model": "secondary_v2_reduced_o2fixed",
        "parameters": list(V2_REDUCED_PARAMETERS),
        "fixed_parameters": list(V2_O2_PARAMETERS)
        + ["kPyrS_stat", "kAldPyr", "kAldS_stat", "kAldO2", "kAcAssim"],
        "n_batches": len(batches),
        "n_residuals": int(len(base_res)),
        "core_policy": "fixed core simulation loaded from secondary_joint theta",
    }
    (RESULTS_DIR / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(fit_table.to_string(index=False))
    print("\nPer-state metrics:")
    print(state_metrics.to_string(index=False))
    print("\nWeak/confounded v2 parameters:")
    print(weak.to_string(index=False) if not weak.empty else "none")
    print(f"\n[done] results={RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
