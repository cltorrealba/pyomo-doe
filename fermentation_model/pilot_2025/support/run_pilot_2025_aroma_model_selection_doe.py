from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

try:
    from pyomo.contrib.doe.utils import rescale_FIM
except Exception:  # pragma: no cover - optional reporting utility
    rescale_FIM = None

SCRIPT_DIR = Path(__file__).resolve().parent
PILOT_DIR = SCRIPT_DIR.parent if SCRIPT_DIR.name == "support" else SCRIPT_DIR
FERMENTATION_MODEL_DIR = PILOT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(PILOT_DIR) not in sys.path:
    sys.path.insert(0, str(PILOT_DIR))
if str(FERMENTATION_MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_MODEL_DIR))

import run_pilot_2025_calibration_estimability as pilot
from shared import aroma_partition_unifac as unifac_partition
from shared import run_new_must_glycerol_estimability_doe as base
from shared import run_secondary_joint_campaign_doe as joint
from shared import run_secondary_v2_model_evaluation as v2


RESULTS_DIR = PILOT_DIR / "results" / "aroma_model_selection_doe"
NOTEBOOK_PATH = PILOT_DIR / "pilot_2025_aroma_model_selection_doe.ipynb"
EXECUTED_NOTEBOOK_PATH = PILOT_DIR / "pilot_2025_aroma_model_selection_doe.executed.ipynb"

COMMON_AROMA_PARAMETERS = pilot.AROMA_TARGETS
CORE_TARGETS = pilot.CORE_TARGETS
SECONDARY_TARGETS = pilot.SECONDARY_TARGETS

EXTRA_BOUNDS = {
    "k_EA_X": (1e-5, 5.0),
    "k_EA_XE": (1e-5, 8.0),
    "k_EA_XE_Nlim": (1e-5, 8.0),
    "k_EA_AcAld": (1e-5, 8.0),
    "q10_EA": (1.0, 3.5),
}
EXTRA_DEFAULTS = {
    "k_EA_X": 0.03,
    "k_EA_XE": 0.05,
    "k_EA_XE_Nlim": 0.02,
    "k_EA_AcAld": 0.02,
    "q10_EA": 1.6,
}

PARAMETER_BOUNDS = {
    **pilot.PARAMETER_BOUNDS,
    **EXTRA_BOUNDS,
}

ETHANOL_HALF_G_L = 40.0
ACETALDEHYDE_HALF_MG_L = 120.0
N_LIMIT_HALF_KG_M3 = joint.N_PHASE_HALF_KG_M3
PARTITION_MODE = os.environ.get("PILOT_AROMA_PARTITION_MODE", "water_ethanol")
PARTITION_MIN_TEMP_C = 10.0
PARTITION_MAX_TEMP_C = 30.0
CONDENSATE_POSITIVE_THRESHOLD = 1e-12
AROMA_RESIDUAL_WEIGHTS = {
    "ethyl_acetate": 2.0,
    "isoamyl_acetate": 1.0,
    "ethyl_octanoate": 1.0,
}


@dataclass(frozen=True)
class AromaVariant:
    name: str
    parameters: tuple[str, ...]
    literature_basis: str
    equation_markdown: str
    description: str


def variant_library() -> dict[str, AromaVariant]:
    common = tuple(COMMON_AROMA_PARAMETERS)
    return {
        "baseline_phase": AromaVariant(
            name="baseline_phase",
            parameters=common,
            literature_basis=(
                "Reference empirical model: aroma synthesis is tied to sugar uptake and split between "
                "nitrogen-associated growth and stationary phases; volatilization is CO2-rate dependent."
            ),
            equation_markdown=(
                "$$r_{EA}=\\left(k_{EA,g}\\phi_N+k_{EA,s}(1-\\phi_N)\\right)q_S$$"
            ),
            description="Original phase-dependent aroma model used in the previous pilot calibration.",
        ),
        "ea_biomass_background": AromaVariant(
            name="ea_biomass_background",
            parameters=common + ("k_EA_X",),
            literature_basis=(
                "Adds a biomass-associated ethyl-acetate formation term, representing acetyltransferase "
                "capacity and natural-must precursor availability not captured by sugar uptake alone."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,X}X$$"
            ),
            description="Tests whether ethyl acetate behaves as a biomass/activity output rather than a direct sugar-rate output.",
        ),
        "ea_ethanol_biomass": AromaVariant(
            name="ea_ethanol_biomass",
            parameters=common + ("k_EA_XE",),
            literature_basis=(
                "Ethyl acetate formation requires ethanol as the alcohol precursor. This variant adds a "
                "saturating ethanol-biomass term while keeping volatilization unchanged."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,XE}X\\frac{E}{K_E+E}$$"
            ),
            description="Tests whether the missing ethyl-acetate source is ethanol/biomass dependent.",
        ),
        "ea_ethanol_temperature": AromaVariant(
            name="ea_ethanol_temperature",
            parameters=common + ("k_EA_XE", "q10_EA"),
            literature_basis=(
                "Fermentative ester synthesis and stripping are temperature-sensitive. This variant lets "
                "the ethanol-biomass production term follow a fitted Q10 factor around 20 C."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,XE}X\\frac{E}{K_E+E}"
                "q_{10,EA}^{(T-20)/10}$$"
            ),
            description="Tests whether the missing ethyl-acetate source is mainly temperature-modulated.",
        ),
        "ea_ethanol_nlimited": AromaVariant(
            name="ea_ethanol_nlimited",
            parameters=common + ("k_EA_XE", "k_EA_XE_Nlim"),
            literature_basis=(
                "Nitrogen limitation changes yeast aroma metabolism. This variant separates a baseline "
                "ethanol-biomass source from an additional nitrogen-limited source."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,XE}X\\frac{E}{K_E+E}"
                "+k_{EA,XE,Nlim}X\\frac{E}{K_E+E}\\frac{K_N}{K_N+N}$$"
            ),
            description="Tests whether ethyl acetate increases when ethanol is available and nitrogen is depleted.",
        ),
        "ea_redox_acetaldehyde": AromaVariant(
            name="ea_redox_acetaldehyde",
            parameters=common + ("k_EA_XE", "k_EA_AcAld"),
            literature_basis=(
                "Ethyl acetate is linked to acetyl-CoA/redox and acetaldehyde/acetate-side metabolism. "
                "The reduced secondary model provides acetaldehyde as a proxy state."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,XE}X\\frac{E}{K_E+E}"
                "+k_{EA,AcAld}X\\frac{AcAld}{K_{AcAld}+AcAld}\\frac{E}{K_E+E}$$"
            ),
            description="Tests whether secondary redox chemistry explains the ethyl-acetate deficit.",
        ),
        "ea_combined_parsimonious": AromaVariant(
            name="ea_combined_parsimonious",
            parameters=common + ("k_EA_X", "k_EA_XE", "k_EA_AcAld"),
            literature_basis=(
                "A compact combined model: background biomass capacity, ethanol precursor availability, "
                "and acetaldehyde/redox proxy. This is the largest variant considered here and is penalized "
                "by BIC/AICc."
            ),
            equation_markdown=(
                "$$r_{EA}=r_{EA,phase}+k_{EA,X}X+k_{EA,XE}X\\frac{E}{K_E+E}"
                "+k_{EA,AcAld}X\\frac{AcAld}{K_{AcAld}+AcAld}\\frac{E}{K_E+E}$$"
            ),
            description="Upper-complexity check for whether several mechanistic proxies are simultaneously needed.",
        ),
    }


def parameter_bounds(name: str) -> tuple[float, float]:
    return PARAMETER_BOUNDS[name]


def clip_theta(theta: dict[str, float]) -> dict[str, float]:
    out = dict(theta)
    for name, (lb, ub) in PARAMETER_BOUNDS.items():
        if name in out and np.isfinite(float(out[name])):
            out[name] = float(np.clip(float(out[name]), lb, ub))
    return out


def load_reference_theta() -> dict[str, float]:
    path = pilot.RESULTS_DIR / "theta_pilot_extended.csv"
    if path.exists():
        theta = pd.read_csv(path, index_col=0).iloc[:, 0].to_dict()
        theta = {str(k): float(v) for k, v in theta.items() if np.isfinite(float(v))}
    else:
        theta = pilot.load_initial_theta()
    for name, value in EXTRA_DEFAULTS.items():
        theta.setdefault(name, value)
    return clip_theta(theta)


def log_vector(theta: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    return np.asarray([math.log(max(float(theta[name]), 1e-16)) for name in parameters], dtype=float)


def theta_from_log(x: np.ndarray, base_theta: dict[str, float], parameters: tuple[str, ...]) -> dict[str, float]:
    theta = dict(base_theta)
    for name, value in zip(parameters, np.asarray(x, dtype=float)):
        theta[name] = float(math.exp(float(value)))
    return clip_theta(theta)


def log_bounds(parameters: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray]:
    lb = []
    ub = []
    for name in parameters:
        lo, hi = parameter_bounds(name)
        lb.append(math.log(float(lo)))
        ub.append(math.log(float(hi)))
    return np.asarray(lb, dtype=float), np.asarray(ub, dtype=float)


def interp_frame_value(frame: pd.DataFrame, column: str, t: float, default: float = 0.0) -> float:
    if frame is None or frame.empty or column not in frame.columns:
        return float(default)
    index = frame.index.to_numpy(dtype=float)
    values = frame[column].to_numpy(dtype=float)
    mask = np.isfinite(index) & np.isfinite(values)
    if not mask.any():
        return float(default)
    return float(np.interp(float(t), index[mask], values[mask]))


def secondary_cache_for(
    batches: list[base.BatchData],
    theta: dict[str, float],
    core_cache: dict[str, pd.DataFrame],
) -> dict[str, pd.DataFrame]:
    out: dict[str, pd.DataFrame] = {}
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            continue
        out[batch.label] = v2.integrate_secondary_v2(batch, theta, core)
    return out


@lru_cache(maxsize=120000)
def partition_k_antoine_unifac_cached(
    species: str,
    temp_c: float,
    ethanol_g_l: float,
    glucose_g_l: float,
    fructose_g_l: float,
    mode: str,
) -> float:
    if mode == "water_ethanol":
        glucose_g_l = 0.0
        fructose_g_l = 0.0
    result = unifac_partition.unifac_partition_K_antoine(
        species,
        float(temp_c),
        float(ethanol_g_l),
        glucose_g_l=float(glucose_g_l),
        fructose_g_l=float(fructose_g_l),
        mode=mode,
        process_min_temp_c=PARTITION_MIN_TEMP_C,
        process_max_temp_c=PARTITION_MAX_TEMP_C,
    )
    return max(float(result["K"]), 1e-12)


def partition_k_literature(
    species: str,
    temp_c: float,
    ethanol_g_l: float,
    glucose_g_l: float = 0.0,
    fructose_g_l: float = 0.0,
) -> float:
    """Gas/liquid concentration ratio from Antoine vapor pressure and UNIFAC."""
    return partition_k_antoine_unifac_cached(
        species,
        round(float(temp_c), 1),
        round(float(ethanol_g_l), 1),
        round(float(glucose_g_l), 1),
        round(float(fructose_g_l), 1),
        PARTITION_MODE,
    )


def write_partition_audit() -> None:
    rows = []
    antoine = unifac_partition.antoine_correlation_data(PARTITION_MIN_TEMP_C, PARTITION_MAX_TEMP_C)
    for species in joint.AROMA_SPECIES:
        for temp_c in (12.0, 18.0, 24.0):
            for ethanol_g_l in (0.0, 60.0, 120.0):
                result = unifac_partition.unifac_partition_K_antoine(
                    species,
                    temp_c,
                    ethanol_g_l,
                    mode=PARTITION_MODE,
                    process_min_temp_c=PARTITION_MIN_TEMP_C,
                    process_max_temp_c=PARTITION_MAX_TEMP_C,
                )
                rows.append(
                    {
                        "species": species,
                        "temp_c": temp_c,
                        "ethanol_g_l": ethanol_g_l,
                        "K_gas_liquid": result["K"],
                        "gamma_unifac": result["gamma"],
                        "Psat_Pa_Antoine": result["Psat_Pa"],
                        "liquid_mode": result["mode"],
                        "antoine_method": result["antoine_method"],
                        "antoine_extrapolated_for_process_window": result["antoine_extrapolated_for_process_window"],
                    }
                )
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "partition_antoine_unifac_audit.csv", index=False)
    pd.DataFrame.from_dict(antoine, orient="index").rename_axis("species").reset_index().to_csv(
        RESULTS_DIR / "partition_antoine_coefficients.csv",
        index=False,
    )


def ea_extra_source(
    variant: AromaVariant,
    theta: dict[str, float],
    core_rates: dict[str, float],
    secondary: pd.DataFrame | None,
    t: float,
) -> float:
    x = max(float(core_rates.get("X", 0.0)), 0.0)
    e = max(float(core_rates.get("E", 0.0)), 0.0)
    n = max(float(core_rates.get("N", 0.0)), 0.0)
    temp = float(core_rates.get("TempC", 20.0))
    e_sat = e / (ETHANOL_HALF_G_L + e + 1e-12)
    n_lim = N_LIMIT_HALF_KG_M3 / (N_LIMIT_HALF_KG_M3 + n + 1e-12)
    acald = max(interp_frame_value(secondary, "AcAld", t, 0.0), 0.0)
    acald_sat = acald / (ACETALDEHYDE_HALF_MG_L + acald + 1e-12)

    extra = 0.0
    if "k_EA_X" in variant.parameters:
        extra += theta["k_EA_X"] * x
    if "k_EA_XE" in variant.parameters:
        temp_factor = 1.0
        if "q10_EA" in variant.parameters:
            temp_factor = theta["q10_EA"] ** ((temp - 20.0) / 10.0)
        extra += theta["k_EA_XE"] * x * e_sat * temp_factor
    if "k_EA_XE_Nlim" in variant.parameters:
        extra += theta["k_EA_XE_Nlim"] * x * e_sat * n_lim
    if "k_EA_AcAld" in variant.parameters:
        extra += theta["k_EA_AcAld"] * x * acald_sat * e_sat
    return float(max(extra, 0.0))


def aroma_rhs_variant(
    t: float,
    y: np.ndarray,
    theta: dict[str, float],
    batch: base.BatchData | base.FutureDesign,
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    variant: AromaVariant,
) -> list[float]:
    liquid = {species: max(float(y[idx]), 0.0) for idx, species in enumerate(joint.AROMA_SPECIES)}
    rates = joint._core_rates(theta, batch, core, float(t))
    phi_growth = rates["N"] / (rates["N"] + joint.N_PHASE_HALF_KG_M3)
    co2_rate = max(joint.CO2_G_PER_G_ETHANOL * rates["ethanol_prod"], 0.0)
    sugar = max(rates["G"] + rates["F"], 0.0)
    out: list[float] = []
    loss_rates: list[float] = []
    for species in joint.AROMA_SPECIES:
        short = joint.AROMA_SHORT[species]
        k_prod = theta[f"k_{short}_growth"] * phi_growth + theta[f"k_{short}_stationary"] * (1.0 - phi_growth)
        prod = k_prod * rates["sugar_uptake"]
        if species == "ethyl_acetate":
            prod += ea_extra_source(variant, theta, rates, secondary, t)
        alpha = theta[f"alpha_{short}_loss"]
        k_lg = partition_k_literature(species, rates["TempC"], rates["E"], rates.get("G", 0.0), rates.get("F", 0.0))
        loss_rate = alpha * k_lg * co2_rate * liquid[species]
        out.append(float(prod - loss_rate))
        loss_rates.append(float(loss_rate))
    out.extend(loss_rates)
    return out


def integrate_aroma_variant_euler(
    batch: base.BatchData | base.FutureDesign,
    theta: dict[str, float],
    core: pd.DataFrame,
    secondary: pd.DataFrame | None,
    variant: AromaVariant,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    time = core.index.to_numpy(dtype=float)
    liquid0 = [float(batch.initials.get(f"{species}_liq", 0.0)) for species in joint.AROMA_SPECIES]
    loss0 = [float(batch.initials.get(f"{species}_cond", 0.0)) for species in joint.AROMA_SPECIES]
    y = np.asarray(liquid0 + loss0, dtype=float)
    rows = [y.copy()]
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        dy = np.asarray(aroma_rhs_variant(float(time[idx - 1]), y, theta, batch, core, secondary, variant), dtype=float)
        y = np.maximum(y + dt * dy, 0.0)
        rows.append(y.copy())
    arr = np.asarray(rows, dtype=float)
    liquid = pd.DataFrame(arr[:, : len(joint.AROMA_SPECIES)], index=time, columns=joint.AROMA_SPECIES)
    condensate = pd.DataFrame(arr[:, len(joint.AROMA_SPECIES) :], index=time, columns=joint.AROMA_SPECIES)
    return liquid, condensate


def aroma_prediction_rows(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: AromaVariant,
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for batch in batches:
        core = core_cache[batch.label]
        secondary = sec_cache.get(batch.label)
        try:
            liq, cond = integrate_aroma_variant_euler(batch, theta, core, secondary, variant)
        except Exception:
            continue
        for species, total_col in joint.AROMA_COLUMNS.items():
            cond_col = pilot.AROMA_CONDENSATE_COLUMNS[species]
            total_obs = np.asarray(batch.observations.get(total_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_obs = np.asarray(batch.observations.get(cond_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_valid = np.isfinite(cond_obs) & (cond_obs > CONDENSATE_POSITIVE_THRESHOLD)
            pred_liq = liq.loc[batch.time, species].to_numpy(dtype=float)
            pred_cond = cond.loc[batch.time, species].to_numpy(dtype=float)
            pred_total = pred_liq + pred_cond
            for idx, t in enumerate(batch.time):
                if np.isfinite(total_obs[idx]) and cond_valid[idx]:
                    retained = max(float(total_obs[idx] - cond_obs[idx]), 0.0)
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "retained", "obs": retained, "pred": float(pred_liq[idx])})
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "condensate", "obs": float(cond_obs[idx]), "pred": float(pred_cond[idx])})
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "total", "obs": float(total_obs[idx]), "pred": float(pred_total[idx])})
                elif np.isfinite(total_obs[idx]):
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "total", "obs": float(total_obs[idx]), "pred": float(pred_total[idx])})
                elif cond_valid[idx]:
                    rows.append({"batch": batch.batch, "time_h": float(t), "species": species, "pool": "condensate", "obs": float(cond_obs[idx]), "pred": float(pred_cond[idx])})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["residual"] = out["pred"] - out["obs"]
    return out


def aroma_residual_variant(
    theta: dict[str, float],
    batches: list[base.BatchData],
    variant: AromaVariant,
    core_cache: dict[str, pd.DataFrame] | None = None,
    sec_cache: dict[str, pd.DataFrame] | None = None,
) -> np.ndarray:
    if core_cache is None:
        core_cache = joint.precompute_core_cache(batches, theta)
    if sec_cache is None:
        sec_cache = secondary_cache_for(batches, theta, core_cache)
    residuals: list[np.ndarray] = []
    for batch in batches:
        core = core_cache.get(batch.label)
        if core is None:
            return np.ones(1000, dtype=float) * 1e6
        secondary = sec_cache.get(batch.label)
        try:
            liq, cond = integrate_aroma_variant_euler(batch, theta, core, secondary, variant)
        except Exception:
            return np.ones(1000, dtype=float) * 1e6
        for species, total_col in joint.AROMA_COLUMNS.items():
            species_weight = float(AROMA_RESIDUAL_WEIGHTS.get(species, 1.0))
            cond_col = pilot.AROMA_CONDENSATE_COLUMNS[species]
            total_obs = np.asarray(batch.observations.get(total_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_obs = np.asarray(batch.observations.get(cond_col, np.full_like(batch.time, np.nan)), dtype=float)
            cond_valid = np.isfinite(cond_obs) & (cond_obs > CONDENSATE_POSITIVE_THRESHOLD)
            pred_liq = liq.loc[batch.time, species].to_numpy(dtype=float)
            pred_cond = cond.loc[batch.time, species].to_numpy(dtype=float)

            mask_both = np.isfinite(total_obs) & cond_valid
            if mask_both.any():
                retained = np.maximum(total_obs[mask_both] - cond_obs[mask_both], 0.0)
                residuals.append(species_weight * (pred_liq[mask_both] - retained) / pilot.aroma_sigma(species, "liquid", retained))
                residuals.append(species_weight * (pred_cond[mask_both] - cond_obs[mask_both]) / pilot.aroma_sigma(species, "condensate", cond_obs[mask_both]))

            mask_total_only = np.isfinite(total_obs) & ~cond_valid
            if mask_total_only.any():
                pred_total = pred_liq[mask_total_only] + pred_cond[mask_total_only]
                residuals.append(species_weight * (pred_total - total_obs[mask_total_only]) / pilot.aroma_sigma(species, "liquid", total_obs[mask_total_only]))

            mask_cond_only = cond_valid & ~np.isfinite(total_obs)
            if mask_cond_only.any():
                residuals.append(species_weight * (pred_cond[mask_cond_only] - cond_obs[mask_cond_only]) / pilot.aroma_sigma(species, "condensate", cond_obs[mask_cond_only]))
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def active_bound_count(theta: dict[str, float], parameters: tuple[str, ...]) -> int:
    count = 0
    for name in parameters:
        lb, ub = parameter_bounds(name)
        value = float(theta[name])
        if value <= lb * 1.01 or value >= ub / 1.01:
            count += 1
    return count


def fit_aroma_variant(
    variant: AromaVariant,
    batches: list[base.BatchData],
    theta0: dict[str, float],
    core_cache: dict[str, pd.DataFrame],
    sec_cache: dict[str, pd.DataFrame],
    n_starts: int,
    max_nfev: int,
    seed: int,
) -> tuple[dict[str, float], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    lb, ub = log_bounds(variant.parameters)
    best_theta = dict(theta0)
    best_wsse = np.inf
    rows = []
    for idx in range(int(n_starts)):
        start = dict(theta0)
        if idx > 0:
            for name in variant.parameters:
                lo, hi = parameter_bounds(name)
                if name in EXTRA_BOUNDS:
                    start[name] = float(math.exp(rng.uniform(math.log(lo), math.log(hi))))
                else:
                    start[name] = float(np.clip(theta0[name] * math.exp(rng.normal(0.0, 0.75)), lo, hi))
        x0 = np.clip(log_vector(start, variant.parameters), lb + 1e-9, ub - 1e-9)

        def fun(x: np.ndarray) -> np.ndarray:
            theta = theta_from_log(x, theta0, variant.parameters)
            data_res = aroma_residual_variant(theta, batches, variant, core_cache, sec_cache)
            priors = []
            for name in variant.parameters:
                scale = 1.15 if name.startswith("alpha_") else 1.8
                if name.startswith("k_EA_") or name == "q10_EA":
                    scale = 4.0
                ref = max(float(theta0[name]), 1e-16)
                priors.append(math.log(max(float(theta[name]), 1e-16) / ref) / scale)
            return np.concatenate([data_res, np.asarray(priors, dtype=float)])

        start_res = fun(x0)
        result = least_squares(
            fun,
            x0,
            bounds=(lb, ub),
            method="trf",
            x_scale="jac",
            loss="soft_l1",
            f_scale=2.0,
            max_nfev=int(max_nfev),
            ftol=1e-6,
            xtol=1e-6,
            gtol=1e-6,
        )
        theta_hat = theta_from_log(result.x, theta0, variant.parameters)
        end_res = fun(result.x)
        data_res = aroma_residual_variant(theta_hat, batches, variant, core_cache, sec_cache)
        wsse_data = float(np.dot(data_res, data_res))
        row = {
            "model": variant.name,
            "start": idx,
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "initial_objective_wsse": float(np.dot(start_res, start_res)),
            "final_objective_wsse": float(np.dot(end_res, end_res)),
            "data_wsse": wsse_data,
            "n_data_residuals": int(len(data_res)),
            "n_parameters": int(len(variant.parameters)),
            "wsse_per_data_residual": wsse_data / max(int(len(data_res)), 1),
            "active_bound_count": active_bound_count(theta_hat, variant.parameters),
        }
        rows.append(row)
        if wsse_data < best_wsse:
            best_wsse = wsse_data
            best_theta = theta_hat
    return clip_theta(best_theta), pd.DataFrame(rows)


def summarize_prediction_metrics(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if pred.empty:
        return pd.DataFrame()
    for keys, group in pred.groupby(["species", "pool"], sort=True):
        species, pool = keys
        obs = group["obs"].to_numpy(dtype=float)
        residual = group["residual"].to_numpy(dtype=float)
        rmse = float(np.sqrt(np.mean(residual**2)))
        mae = float(np.mean(np.abs(residual)))
        bias = float(np.mean(residual))
        denom = max(float(np.nanmedian(np.abs(obs))), 1e-9)
        rows.append(
            {
                "species": species,
                "pool": pool,
                "n": int(len(group)),
                "rmse": rmse,
                "mae": mae,
                "bias": bias,
                "median_obs": float(np.nanmedian(obs)),
                "median_pred": float(np.nanmedian(group["pred"].to_numpy(dtype=float))),
                "relative_rmse_to_median": rmse / denom,
                "relative_bias_to_median": bias / denom,
            }
        )
    return pd.DataFrame(rows)


def information_criteria(data_wsse: float, n: int, k: int) -> tuple[float, float]:
    n = max(int(n), 1)
    k = max(int(k), 1)
    aic = float(data_wsse + 2 * k)
    if n > k + 1:
        aic += float((2 * k * (k + 1)) / (n - k - 1))
    bic = float(data_wsse + k * math.log(n))
    return aic, bic


def finite_difference_jacobian_generic(
    theta: dict[str, float],
    parameters: tuple[str, ...],
    residual_fun,
    step: float,
) -> tuple[np.ndarray, np.ndarray]:
    base_res = residual_fun(theta)
    cols = []
    for name in parameters:
        theta_plus = dict(theta)
        theta_minus = dict(theta)
        lb, ub = parameter_bounds(name)
        theta_plus[name] = float(np.clip(float(theta[name]) * math.exp(step), lb, ub))
        theta_minus[name] = float(np.clip(float(theta[name]) * math.exp(-step), lb, ub))
        r_plus = residual_fun(theta_plus)
        r_minus = residual_fun(theta_minus)
        if len(r_plus) != len(base_res) or len(r_minus) != len(base_res):
            cols.append(np.zeros_like(base_res))
        else:
            cols.append((r_plus - r_minus) / (2.0 * step))
    return np.column_stack(cols), base_res


def fim_metrics_generic(fim: np.ndarray, prefix: str = "") -> dict[str, float]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eig = np.linalg.eigvalsh(fim)
    max_eig = float(np.max(eig)) if eig.size else np.nan
    floor = max(max_eig * 1e-12, np.finfo(float).tiny) if np.isfinite(max_eig) and max_eig > 0 else np.finfo(float).tiny
    eig_pos = np.clip(eig, floor, None)
    cov = pilot.stable_inverse(fim)
    return {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}max_eigenvalue": max_eig,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if np.isfinite(max_eig) and max_eig > 0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace_inv": float(np.trace(cov)),
        f"{prefix}rank_1e-8": int(np.sum(eig > max_eig * 1e-8)) if np.isfinite(max_eig) and max_eig > 0 else 0,
    }


def fim_diagnostics_generic(
    fim: np.ndarray,
    theta: dict[str, float],
    parameters: tuple[str, ...],
    label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eigvals, eigvecs = np.linalg.eigh(fim)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    max_eig = float(np.max(eigvals)) if eigvals.size else np.nan
    spectrum = pd.DataFrame(
        {
            "analysis": label,
            "direction": np.arange(1, len(eigvals) + 1),
            "eigenvalue": eigvals,
            "relative_eigenvalue": eigvals / max_eig if np.isfinite(max_eig) and max_eig > 0.0 else np.nan,
        }
    )
    weak_rows = []
    for idx in range(min(8, eigvecs.shape[1])):
        vec = eigvecs[:, idx]
        top = np.argsort(np.abs(vec))[::-1][:8]
        weak_rows.append(
            {
                "analysis": label,
                "weak_direction": idx + 1,
                "eigenvalue": float(eigvals[idx]),
                "dominant_parameters": ", ".join(parameters[i] for i in top),
                "dominant_abs_loadings": ", ".join(f"{abs(vec[i]):.3f}" for i in top),
            }
        )
    cov = pilot.stable_inverse(fim)
    rows = []
    for idx, name in enumerate(parameters):
        std_log = float(math.sqrt(max(cov[idx, idx], 0.0)))
        lb, ub = parameter_bounds(name)
        value = float(theta[name])
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
                "analysis": label,
                "parameter": name,
                "theta": value,
                "std_log_approx": std_log,
                "approx_95_multiplier": float(math.exp(1.96 * min(std_log, 20.0))),
                "fim_diag": float(fim[idx, idx]),
                "active_bound": bool(active),
                "classification": cls,
            }
        )
    return spectrum, pd.DataFrame(weak_rows), pd.DataFrame(rows)


def add_extra_designs(model_data: pd.DataFrame, designs: dict[str, base.FutureDesign]) -> dict[str, base.FutureDesign]:
    out = dict(designs)
    initials = pilot.natural_initials(model_data)
    h = pilot.NATURAL_DESIGN_HORIZON_H
    zero = {channel: tuple() for channel in base.INPUT_CHANNELS}
    extras = [
        (
            "natural_pilot_EA_warm_early_noN",
            "EA_temperature_Nstress",
            (24.0, 23.0, 19.0, 17.0),
            (),
            "Warm early natural fermentation without nutrient pulse to excite ethanol-biomass and nitrogen-limited EA directions.",
        ),
        (
            "natural_pilot_EA_warm_early_lateN",
            "EA_lateN",
            (23.0, 22.0, 20.0, 18.0),
            ((96.0, 0.045),),
            "Warm early fermentation followed by late nitrogen to separate stationary and N-limited aroma terms.",
        ),
        (
            "natural_pilot_EA_cold_warm_Nsplit",
            "EA_temperature_Nsplit",
            (13.0, 18.0, 23.0, 19.0),
            ((36.0, 0.025), (84.0, 0.035)),
            "Cold start then warm acceleration with split nitrogen to test temperature and N response.",
        ),
        (
            "natural_pilot_EA_cold_retention_noN",
            "EA_retention_reference",
            (13.0, 13.0, 16.0, 18.0),
            (),
            "Low-temperature no-pulse comparator to decouple synthesis from CO2 stripping.",
        ),
    ]
    for name, family, temps, n_pulses, rationale in extras:
        pulses = dict(zero)
        pulses["N"] = pilot.snap_pulses(h, n_pulses)
        out[name] = base.FutureDesign(
            name=name,
            family=family,
            medium="natural",
            horizon_h=h,
            initials=dict(initials),
            temperature_segments=tuple(float(t) for t in temps),
            pulses={channel: tuple(rows) for channel, rows in pulses.items()},
            rationale=rationale,
        )
    return out


def combined_residual_selected(
    theta: dict[str, float],
    primary_batches: list[base.BatchData],
    extended_batches: list[base.BatchData],
    co2_downsampled: pd.DataFrame,
    variant: AromaVariant,
) -> np.ndarray:
    residuals = [base.residual_vector(theta, primary_batches)]
    core_cache = joint.precompute_core_cache(extended_batches, theta)
    sec_cache = secondary_cache_for(extended_batches, theta, core_cache)
    residuals.append(pilot.secondary_residual_with_co2(theta, extended_batches, co2_downsampled, core_cache))
    residuals.append(aroma_residual_variant(theta, extended_batches, variant, core_cache, sec_cache))
    return np.concatenate([r for r in residuals if len(r)])


def future_residual_selected(theta: dict[str, float], design: base.FutureDesign, variant: AromaVariant) -> np.ndarray:
    sample_times = pilot.future_sample_times(design)
    co2_times = pilot.future_co2_times(design)
    all_times = np.asarray(sorted(set(sample_times).union(set(co2_times)).union({0.0, float(design.horizon_h)})), dtype=float)
    core = base.simulate(design, theta, all_times)
    if core is None:
        return np.ones(1000, dtype=float) * 1e6
    secondary = v2.integrate_secondary_v2(design, theta, core)
    try:
        liq, cond = integrate_aroma_variant_euler(design, theta, core, secondary, variant)
    except Exception:
        return np.ones(1000, dtype=float) * 1e6

    residuals: list[np.ndarray] = []
    valid_samples = [float(t) for t in sample_times if float(t) in core.index]
    if valid_samples:
        core_sample = core.loc[valid_samples]
        for state in base.STATE_NAMES:
            center = core_sample[state].to_numpy(dtype=float)
            residuals.append(center / base.sigma_for_state(state, center))
        sec_sample = secondary.loc[valid_samples]
        for state in ("Pyr", "AcAld", "Acetate"):
            residuals.append(sec_sample[state].to_numpy(dtype=float) / joint.SIGMA[state])
        aroma_sample = liq.loc[valid_samples]
        cond_sample = cond.loc[valid_samples]
        for species in joint.AROMA_SPECIES:
            liq_center = aroma_sample[species].to_numpy(dtype=float)
            cond_center = cond_sample[species].to_numpy(dtype=float)
            residuals.append(liq_center / pilot.aroma_sigma(species, "liquid", liq_center))
            residuals.append(cond_center / pilot.aroma_sigma(species, "condensate", cond_center))
    valid_co2 = [float(t) for t in co2_times if float(t) in core.index]
    if valid_co2:
        rate = pilot.co2_rate_from_core(theta, design, core, np.asarray(valid_co2, dtype=float))
        sigma = np.maximum(0.05, 0.20 * np.maximum(rate, 0.05))
        residuals.append(rate / sigma)
    return np.concatenate(residuals) if residuals else np.array([], dtype=float)


def candidate_fim_selected(theta: dict[str, float], design: base.FutureDesign, variant: AromaVariant, parameters: tuple[str, ...], step: float) -> np.ndarray:
    jac, _ = finite_difference_jacobian_generic(theta, parameters, lambda th: future_residual_selected(th, design, variant), step)
    fim = jac.T @ jac
    return 0.5 * (fim + fim.T)


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = pilot.stable_inverse(prior)
    post_cov = pilot.stable_inverse(combined)
    rows: dict[str, float] = {}
    reductions = []
    aroma_reductions = []
    for idx, name in enumerate(parameters):
        before = float(prior_cov[idx, idx])
        after = float(post_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        reduction = 1.0 - ratio if np.isfinite(ratio) else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = reduction
        if np.isfinite(reduction):
            reductions.append(reduction)
            if name in COMMON_AROMA_PARAMETERS or name.startswith("k_EA_") or name == "q10_EA":
                aroma_reductions.append(reduction)
    rows["target_mean_var_reduction"] = float(np.nanmean(reductions)) if reductions else np.nan
    rows["target_worst_var_reduction"] = float(np.nanmin(reductions)) if reductions else np.nan
    rows["aroma_mean_var_reduction"] = float(np.nanmean(aroma_reductions)) if aroma_reductions else np.nan
    rows["aroma_worst_var_reduction"] = float(np.nanmin(aroma_reductions)) if aroma_reductions else np.nan
    return rows


def score_fim(fim: np.ndarray, objective: str = "hybrid") -> float:
    metrics = fim_metrics_generic(fim)
    if objective == "d_opt":
        return float(metrics["logdet"])
    if objective == "e_opt":
        return math.log(max(float(metrics["min_eigenvalue"]), 1e-18))
    min_rel = max(float(metrics["min_relative_eigenvalue"]), 1e-18)
    return float(metrics["logdet"]) + 2.0 * math.log(min_rel) - 0.05 * math.log(max(float(metrics["trace_inv"]), 1e-18))


def rank_candidates(candidate_fims: dict[str, np.ndarray], designs: dict[str, base.FutureDesign], prior: np.ndarray, parameters: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for name, fim in candidate_fims.items():
        combined = prior + fim
        design = designs[name]
        rows.append(
            {
                "candidate": name,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                **fim_metrics_generic(fim, prefix="candidate_"),
                **fim_metrics_generic(combined, prefix="combined_"),
                **variance_reduction(prior, combined, parameters),
                "hybrid_score": score_fim(combined, "hybrid"),
                "dopt_score": score_fim(combined, "d_opt"),
                "eopt_score": score_fim(combined, "e_opt"),
                "rationale": design.rationale,
            }
        )
    return pd.DataFrame(rows).sort_values(["hybrid_score", "combined_logdet"], ascending=False).reset_index(drop=True)


def greedy_select(
    candidate_fims: dict[str, np.ndarray],
    designs: dict[str, base.FutureDesign],
    prior: np.ndarray,
    parameters: tuple[str, ...],
    campaign_size: int,
    objective: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    selected: list[str] = []
    remaining = set(candidate_fims)
    current = prior.copy()
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + candidate_fims[name]
            score = score_fim(trial, objective=objective)
            if score > best_score:
                best = name
                best_score = score
                best_metrics = {**fim_metrics_generic(trial, prefix="campaign_"), **variance_reduction(prior, trial, parameters)}
        if best is None:
            break
        current = current + candidate_fims[best]
        remaining.remove(best)
        selected.append(best)
        design = designs[best]
        rows.append(
            {
                "objective": objective,
                "campaign_order": order,
                "candidate": best,
                "family": design.family,
                "medium": design.medium,
                "horizon_h": design.horizon_h,
                "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                "score": float(best_score),
                "rationale": design.rationale,
                **(best_metrics or {}),
            }
        )
    return pd.DataFrame(rows), current


def select_best_model(model_summary: pd.DataFrame) -> str:
    table = model_summary.copy()
    table["selection_penalty"] = 0.0
    table["selection_penalty"] += 20.0 * table["active_bound_count"].astype(float)
    table["selection_penalty"] += 80.0 * np.maximum(table["ea_retained_abs_relative_bias"].astype(float) - 0.35, 0.0)
    table["selection_penalty"] += 40.0 * np.maximum(table["ea_retained_relative_rmse"].astype(float) - 0.75, 0.0)
    table["selection_penalty"] += 80.0 * np.maximum(table["ea_total_abs_relative_bias"].astype(float) - 0.30, 0.0)
    table["selection_penalty"] += 40.0 * np.maximum(table["ea_total_relative_rmse"].astype(float) - 0.55, 0.0)
    table["selection_score"] = table["bic"] + table["selection_penalty"]
    return str(table.sort_values(["selection_score", "bic", "data_wsse"]).iloc[0]["model"])


def plot_model_comparison(model_summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ordered = model_summary.sort_values("selection_score")
    x = np.arange(len(ordered))
    axes[0].bar(x, ordered["bic"])
    axes[0].set_title("BIC")
    axes[1].bar(x, ordered["ea_retained_rmse"])
    axes[1].set_title("Ethyl acetate retained RMSE")
    axes[2].bar(x, ordered["selection_score"])
    axes[2].set_title("Selection score")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(ordered["model"], rotation=60, ha="right", fontsize=8)
        ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "model_selection_comparison.png", dpi=180)
    plt.close(fig)


def plot_selected_fit(pred: pd.DataFrame, selected_model: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for batch, group in pred.groupby("batch", sort=True):
        fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
        for ax, (species, pool) in zip(axes.ravel(), [(s, p) for s in joint.AROMA_SPECIES for p in ("retained", "condensate", "total")]):
            sub = group[group["species"].eq(species) & group["pool"].eq(pool)].sort_values("time_h")
            if sub.empty:
                ax.axis("off")
                continue
            ax.scatter(sub["time_h"], sub["obs"], s=18, color="black", label="obs")
            ax.plot(sub["time_h"], sub["pred"], color="tab:blue", lw=1.5, label="pred")
            ax.set_title(f"{species} {pool}")
            ax.grid(True, alpha=0.25)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="upper right")
        fig.suptitle(f"{selected_model}: aroma fit for pilot {batch}", y=0.995)
        fig.tight_layout()
        fig.savefig(output_dir / f"selected_aroma_fit_{batch}.png", dpi=170)
        plt.close(fig)


def plot_design_inputs(selected: pd.DataFrame, designs: dict[str, base.FutureDesign], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in selected["candidate"].drop_duplicates().astype(str):
        design = designs[name]
        time = np.linspace(0.0, float(design.horizon_h), 500)
        fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
        axes[0].step(time, [base.temperature_at(design, t) for t in time], where="post", color="tab:red")
        axes[0].set_ylabel("Temperature (C)")
        axes[0].grid(True, alpha=0.25)
        n_pulses = design.pulses.get("N", tuple())
        if n_pulses:
            axes[1].vlines([t for t, _a in n_pulses], 0.0, [a * 1000.0 for _t, a in n_pulses], color="tab:green", lw=3)
        axes[1].set_ylabel("N pulse (mg/L)")
        axes[1].grid(True, alpha=0.25)
        samples = pilot.future_sample_times(design)
        axes[2].vlines(samples, 0.0, 1.0, color="tab:blue", lw=1)
        axes[2].set_ylabel("Samples")
        axes[2].set_xlabel("time (h)")
        axes[2].grid(True, alpha=0.25)
        fig.suptitle(f"{design.name}: {design.rationale}", fontsize=10, y=0.99)
        fig.tight_layout()
        fig.savefig(output_dir / f"candidate_inputs_{design.name}.png", dpi=170)
        plt.close(fig)


def write_report(
    variants: dict[str, AromaVariant],
    selected_model: str,
    model_summary: pd.DataFrame,
    estimability: pd.DataFrame,
    weak: pd.DataFrame,
    ranking: pd.DataFrame,
    selected: pd.DataFrame,
) -> None:
    path = RESULTS_DIR / "pilot_2025_aroma_model_selection_doe_report.md"
    best_variant = variants[selected_model]
    with path.open("w", encoding="utf-8") as f:
        f.write("# Pilot 2025 Aroma Model Selection and MBDoE\n\n")
        f.write("## Selected structure\n\n")
        f.write(f"Selected model: `{selected_model}`.\n\n")
        f.write(best_variant.description + "\n\n")
        f.write(best_variant.equation_markdown + "\n\n")
        f.write("## Gas-liquid equilibrium and CO2 stripping\n\n")
        f.write(
            "Aroma loss is computed as `r_loss,i = alpha_i K_i(T,E) q_CO2 C_L,i`. "
            "`K_i` is evaluated directly from Antoine vapor pressure and original UNIFAC activity coefficients "
            "for a dilute aroma in a water-ethanol liquid mixture. The CO2 stripping proportionality follows the "
            "Mouret-style mass-balance assumption that volatile loss scales with fermentation CO2 flow.\n\n"
        )
        f.write(
            "Condensate entries equal to exactly zero are treated as missing/below-reporting-limit values, not as exact "
            "dynamic measurements of zero accumulated condensate. Positive condensate values are used to reconstruct "
            "the retained liquid pool and to constrain the accumulated loss state.\n\n"
        )
        f.write("## Model selection table\n\n")
        show = [
            "model",
            "n_parameters",
            "data_wsse",
            "wsse_per_data_residual",
            "aicc",
            "bic",
            "selection_score",
            "active_bound_count",
            "ea_retained_rmse",
            "ea_retained_relative_rmse",
            "ea_retained_relative_bias",
            "ea_total_rmse",
            "ea_total_relative_rmse",
            "ea_total_relative_bias",
        ]
        f.write(model_summary[show].sort_values("selection_score").to_markdown(index=False))
        f.write("\n\n")
        f.write("## Selected model estimability\n\n")
        f.write(estimability[["parameter", "theta", "std_log_approx", "approx_95_multiplier", "active_bound", "classification"]].to_markdown(index=False))
        f.write("\n\n")
        f.write("## Weak FIM directions\n\n")
        f.write(weak.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Candidate ranking\n\n")
        show_rank = [
            "candidate",
            "family",
            "combined_logdet",
            "combined_min_relative_eigenvalue",
            "aroma_mean_var_reduction",
            "aroma_worst_var_reduction",
            "N_pulses_kg_m3",
            "rationale",
        ]
        f.write(ranking[show_rank].head(12).to_markdown(index=False))
        f.write("\n\n")
        f.write("## Selected campaign\n\n")
        f.write(selected[["campaign_order", "candidate", "family", "campaign_logdet", "campaign_min_relative_eigenvalue", "aroma_mean_var_reduction", "aroma_worst_var_reduction", "N_pulses_kg_m3", "rationale"]].to_markdown(index=False))


def create_notebook(variants: dict[str, AromaVariant], selected_model: str) -> None:
    selected = variants[selected_model]
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Pilot 2025 Aroma Model Selection and MBDoE\n\n"
            "This notebook documents the model-structure iteration used to improve the aroma layer for the pilot-scale natural-must dataset. "
            "The goal is not only to reduce error, but to choose a structure that is interpretable, numerically identifiable, and useful for future model-based experimental design."
        ),
        nbformat.v4.new_markdown_cell(
            "## Workflow\n\n"
            "1. Load the curated pilot dataset, including the `25171` effective start correction at `t_original = 100 h`.\n"
            "2. Keep the calibrated primary fermentation, glycerol, and secondary v2 model as the current prior.\n"
            "3. Fit several aroma kinetic variants that differ only in the ethyl-acetate production structure.\n"
            "4. Compare variants using weighted residual error, AICc/BIC, active-bound diagnostics, and ethyl-acetate retained-pool bias.\n"
            "5. Recompute the local Fisher information matrix for the selected structure.\n"
            "6. Evaluate candidate natural-must designs using FIM accumulation and select an optimal campaign with a hybrid criterion."
        ),
        nbformat.v4.new_markdown_cell(
            "## Observation Model\n\n"
            "For each aroma species `i`, the workbook contains total aroma and condenser-equivalent aroma:\n\n"
            "$$C_{total,i}^{obs}=C_{wine,i}^{obs}+C_{cond,i}^{obs}$$\n\n"
            "so the retained liquid concentration fitted by the ODE is\n\n"
            "$$C_{wine,i}^{obs}=C_{total,i}^{obs}-C_{cond,i}^{obs}.$$\n\n"
            "This reconstruction is applied only when the condenser-equivalent value is positive. "
            "Exact zero entries in the condenser columns are treated as missing or below-reporting-limit records, because the condenser behaves as an accumulated trap rather than a continuously sampled online analyzer. "
            "When total aroma is available but positive condenser information is not, the fitted observation is `C_total = C_L + C_cond`.\n\n"
            "The model states are\n\n"
            "$$x_{a,i}=\\left[C_{L,i}, C_{cond,i}\\right]^T$$\n\n"
            "with dynamics\n\n"
            "$$\\frac{dC_{L,i}}{dt}=r_{prod,i}-r_{loss,i}, \\qquad \\frac{dC_{cond,i}}{dt}=r_{loss,i}.$$\n\n"
            "Volatilization is kept common across variants and follows the CO2-stripping balance used in the Mouret-style aroma partition framework:\n\n"
            "$$r_{loss,i}=\\alpha_i K_i(T,E)q_{CO2}C_{L,i}.$$\n\n"
            "In this run, `K_i` is not fitted. It is computed as\n\n"
            "$$K_i(T,E)=\\frac{\\gamma_i^{UNIFAC}(T,x_{water},x_{ethanol})P_i^{sat,Antoine}(T)}{RT C_{tot,L}}$$\n\n"
            "where `gamma_i` is the infinite-dilution activity coefficient from original UNIFAC for the water-ethanol mixture, "
            "`P_i^{sat}` is calculated from Antoine coefficients, and `C_{tot,L}` is the total liquid molar concentration. "
            "`alpha_i` remains an effective capture/loss scaling parameter estimated from the condenser data. "
            "The audit tables `partition_antoine_coefficients.csv` and `partition_antoine_unifac_audit.csv` document the thermodynamic inputs."
        ),
        nbformat.v4.new_markdown_cell(
            "## Candidate Kinetic Structures\n\n"
            + "\n\n".join(
                f"### `{variant.name}`\n\n{variant.description}\n\n{variant.equation_markdown}\n\n{variant.literature_basis}"
                for variant in variants.values()
            )
        ),
        nbformat.v4.new_markdown_cell(
            "## Statistical Criteria\n\n"
            "Let `r(theta)` be the weighted residual vector. The data-only weighted sum of squared errors is\n\n"
            "$$WSSE= r(\\theta)^T r(\\theta).$$\n\n"
            "For model comparison, the notebook reports\n\n"
            "$$AIC_c = WSSE + 2k + \\frac{2k(k+1)}{n-k-1},$$\n\n"
            "$$BIC = WSSE + k\\log(n),$$\n\n"
            "where `k` is the number of fitted aroma parameters and `n` is the number of aroma residuals. "
            "Because ethyl acetate is the priority output and the previous structure showed systematic underprediction, ethyl-acetate residuals are weighted twice as strongly as the other aroma residuals in the fitting objective. "
            "The final selection score is BIC plus penalties for active bounds and unacceptable ethyl-acetate retained-pool bias/RMSE. "
            "This prevents selection of a flexible model that only improves the global objective while still failing the main problematic output."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n"
            "cwd = Path.cwd()\n"
            "if (cwd / 'results' / 'aroma_model_selection_doe').exists():\n"
            "    RESULTS = cwd / 'results' / 'aroma_model_selection_doe'\n"
            "else:\n"
            "    RESULTS = cwd / 'fermentation_model' / 'pilot_2025' / 'results' / 'aroma_model_selection_doe'\n"
            "model_summary = pd.read_csv(RESULTS / 'model_selection_summary.csv')\n"
            "variant_metrics = pd.read_csv(RESULTS / 'variant_pool_metrics.csv')\n"
            "theta = pd.read_csv(RESULTS / 'theta_selected_aroma_model.csv', index_col=0)\n"
            "estimability = pd.read_csv(RESULTS / 'parameter_estimability_selected_current.csv')\n"
            "weak = pd.read_csv(RESULTS / 'weak_directions_selected_current.csv')\n"
            "ranking = pd.read_csv(RESULTS / 'candidate_ranking_selected_model.csv')\n"
            "campaign = pd.read_csv(RESULTS / 'selected_campaign_hybrid_selected_model.csv')\n"
            "model_summary.sort_values('selection_score')"
        ),
        nbformat.v4.new_markdown_cell("## Selected Parameter Vector"),
        nbformat.v4.new_code_cell("theta"),
        nbformat.v4.new_markdown_cell("## Pool-Level Fit Metrics"),
        nbformat.v4.new_code_cell("variant_metrics.sort_values(['model', 'species', 'pool']).head(80)"),
        nbformat.v4.new_markdown_cell("## Model-Comparison Plot"),
        nbformat.v4.new_code_cell(
            "plot = RESULTS / 'plots' / 'model_selection_comparison.png'\n"
            "display(Image(filename=str(plot)))"
        ),
        nbformat.v4.new_markdown_cell(
            f"## Selected Structure: `{selected_model}`\n\n"
            f"{selected.description}\n\n"
            f"{selected.equation_markdown}\n\n"
            "Interpretation: this is the simplest structure that best balances ethyl-acetate fit, global aroma error, parameter parsimony, and numerical diagnostics."
        ),
        nbformat.v4.new_markdown_cell(
            "## Fisher Information Matrix\n\n"
            "The local sensitivity matrix is computed by centered finite differences in log-parameter space:\n\n"
            "$$J_{m,j}=\\frac{r_m(\\theta_j e^{h})-r_m(\\theta_j e^{-h})}{2h}.$$\n\n"
            "The Fisher information matrix is then\n\n"
            "$$F = J^T J.$$\n\n"
            "Because perturbations are performed in log-parameter space, the covariance approximation is also in log-parameter units. "
            "The reported `approx_95_multiplier` is `exp(1.96 sigma_log)`."
        ),
        nbformat.v4.new_code_cell("estimability.sort_values('std_log_approx', ascending=False)"),
        nbformat.v4.new_markdown_cell("## Weak Eigen-Directions"),
        nbformat.v4.new_code_cell("weak"),
        nbformat.v4.new_markdown_cell(
            "## MBDoE Criterion\n\n"
            "Candidate experiments are evaluated by adding their expected FIM to the current-data FIM:\n\n"
            "$$F_{combined}=F_{current}+\\sum_{e \\in \\mathcal{E}}F_e.$$\n\n"
            "The ranking reports D-optimality through `logdet(F)`, E-like robustness through the minimum relative eigenvalue, and A-like behavior through `trace(inv(F))`. "
            "The selected campaign uses the hybrid score\n\n"
            "$$\\Phi_{hybrid}=\\log\\det(F)+2\\log(\\lambda_{min}/\\lambda_{max})-0.05\\log(\\mathrm{trace}(F^{-1})).$$\n\n"
            "This keeps the D-optimal volume objective but discourages designs that leave a nearly unobservable direction."
        ),
        nbformat.v4.new_code_cell(
            "cols = ['candidate', 'family', 'combined_logdet', 'combined_min_relative_eigenvalue', "
            "'aroma_mean_var_reduction', 'aroma_worst_var_reduction', 'N_pulses_kg_m3', 'rationale']\n"
            "ranking[cols].head(12)"
        ),
        nbformat.v4.new_markdown_cell("## Selected Campaign"),
        nbformat.v4.new_code_cell(
            "campaign[['campaign_order', 'candidate', 'family', 'campaign_logdet', 'campaign_min_relative_eigenvalue', "
            "'aroma_mean_var_reduction', 'aroma_worst_var_reduction', 'N_pulses_kg_m3', 'rationale']]"
        ),
        nbformat.v4.new_markdown_cell("## Input Profiles for Selected Designs"),
        nbformat.v4.new_code_cell(
            "for path in sorted((RESULTS / 'plots' / 'designs').glob('candidate_inputs_*.png')):\n"
            "    display(Image(filename=str(path)))"
        ),
        nbformat.v4.new_markdown_cell("## Selected Model Aroma Fits"),
        nbformat.v4.new_code_cell(
            "for path in sorted((RESULTS / 'plots' / 'selected_fit').glob('selected_aroma_fit_*.png')):\n"
            "    display(Image(filename=str(path)))"
        ),
    ]
    nbformat.write(nb, NOTEBOOK_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description="Pilot 2025 aroma kinetic model selection and MBDoE.")
    parser.add_argument("--variant-starts", type=int, default=5)
    parser.add_argument("--max-variant-nfev", type=int, default=220)
    parser.add_argument("--sensitivity-step", type=float, default=1e-2)
    parser.add_argument("--campaign-size", type=int, default=6)
    parser.add_argument("--skip-doe", action="store_true")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = RESULTS_DIR / "plots"
    (plot_dir / "selected_fit").mkdir(parents=True, exist_ok=True)
    (plot_dir / "designs").mkdir(parents=True, exist_ok=True)
    write_partition_audit()

    print("[load] pilot data and reference theta", flush=True)
    model_data, co2_curated, co2_decisions = pilot.load_pilot_model_data()
    co2_down = pilot.downsample_co2(co2_curated)
    primary_batches = pilot.make_primary_batches(model_data)
    extended_batches = pilot.make_extended_batches(model_data)
    theta0 = load_reference_theta()
    for name, value in EXTRA_DEFAULTS.items():
        theta0.setdefault(name, value)
    theta0 = clip_theta(theta0)
    core_cache = joint.precompute_core_cache(extended_batches, theta0)
    sec_cache = secondary_cache_for(extended_batches, theta0, core_cache)

    variants = variant_library()
    all_fit_rows = []
    all_metric_rows = []
    model_rows = []
    theta_by_model: dict[str, dict[str, float]] = {}
    pred_by_model: dict[str, pd.DataFrame] = {}
    baseline_theta: dict[str, float] | None = None

    for idx, variant in enumerate(variants.values(), start=1):
        print(f"[variant {idx}/{len(variants)}] {variant.name}", flush=True)
        theta_start = dict(baseline_theta) if baseline_theta is not None else dict(theta0)
        for name in variant.parameters:
            if name in EXTRA_BOUNDS and baseline_theta is not None:
                theta_start[name] = EXTRA_BOUNDS[name][0] * 1.05
            else:
                theta_start.setdefault(name, EXTRA_DEFAULTS.get(name, theta0.get(name, 1.0)))
        theta_hat, fit_rows = fit_aroma_variant(
            variant,
            extended_batches,
            theta_start,
            core_cache,
            sec_cache,
            n_starts=args.variant_starts,
            max_nfev=args.max_variant_nfev,
            seed=100 + idx,
        )
        theta_by_model[variant.name] = theta_hat
        pd.Series(theta_hat).to_csv(RESULTS_DIR / f"theta_{variant.name}.csv")
        pd.Series({name: theta_hat[name] for name in variant.parameters}).to_csv(RESULTS_DIR / f"theta_aroma_{variant.name}.csv")
        if variant.name == "baseline_phase":
            baseline_theta = dict(theta_hat)
        fit_rows.to_csv(RESULTS_DIR / f"fit_multistart_{variant.name}.csv", index=False)
        all_fit_rows.append(fit_rows)
        pred = aroma_prediction_rows(theta_hat, extended_batches, variant, core_cache, sec_cache)
        pred["model"] = variant.name
        pred.to_csv(RESULTS_DIR / f"predictions_{variant.name}.csv", index=False)
        pred_by_model[variant.name] = pred
        metrics = summarize_prediction_metrics(pred)
        metrics["model"] = variant.name
        all_metric_rows.append(metrics)
        best_fit = fit_rows.sort_values("data_wsse").iloc[0].to_dict()
        aicc, bic = information_criteria(float(best_fit["data_wsse"]), int(best_fit["n_data_residuals"]), int(best_fit["n_parameters"]))
        ea_ret = metrics[metrics["species"].eq("ethyl_acetate") & metrics["pool"].eq("retained")]
        ea_total = metrics[metrics["species"].eq("ethyl_acetate") & metrics["pool"].eq("total")]
        if ea_ret.empty:
            ea_rmse = np.nan
            ea_rel_rmse = np.nan
            ea_rel_bias = np.nan
        else:
            ea_rmse = float(ea_ret.iloc[0]["rmse"])
            ea_rel_rmse = float(ea_ret.iloc[0]["relative_rmse_to_median"])
            ea_rel_bias = float(ea_ret.iloc[0]["relative_bias_to_median"])
        if ea_total.empty:
            ea_total_rmse = np.nan
            ea_total_rel_rmse = np.nan
            ea_total_rel_bias = np.nan
        else:
            ea_total_rmse = float(ea_total.iloc[0]["rmse"])
            ea_total_rel_rmse = float(ea_total.iloc[0]["relative_rmse_to_median"])
            ea_total_rel_bias = float(ea_total.iloc[0]["relative_bias_to_median"])
        model_rows.append(
            {
                "model": variant.name,
                "n_parameters": int(best_fit["n_parameters"]),
                "data_wsse": float(best_fit["data_wsse"]),
                "wsse_per_data_residual": float(best_fit["wsse_per_data_residual"]),
                "n_data_residuals": int(best_fit["n_data_residuals"]),
                "aicc": aicc,
                "bic": bic,
                "active_bound_count": int(best_fit["active_bound_count"]),
                "ea_retained_rmse": ea_rmse,
                "ea_retained_relative_rmse": ea_rel_rmse,
                "ea_retained_relative_bias": ea_rel_bias,
                "ea_retained_abs_relative_bias": abs(ea_rel_bias) if np.isfinite(ea_rel_bias) else np.inf,
                "ea_total_rmse": ea_total_rmse,
                "ea_total_relative_rmse": ea_total_rel_rmse,
                "ea_total_relative_bias": ea_total_rel_bias,
                "ea_total_abs_relative_bias": abs(ea_total_rel_bias) if np.isfinite(ea_total_rel_bias) else np.inf,
                "description": variant.description,
                "literature_basis": variant.literature_basis,
            }
        )

    fit_summary = pd.concat(all_fit_rows, ignore_index=True, sort=False)
    fit_summary.to_csv(RESULTS_DIR / "variant_fit_multistart_summary.csv", index=False)
    pool_metrics = pd.concat(all_metric_rows, ignore_index=True, sort=False)
    pool_metrics.to_csv(RESULTS_DIR / "variant_pool_metrics.csv", index=False)
    model_summary = pd.DataFrame(model_rows)
    selected_model = select_best_model(model_summary)
    model_summary["selected"] = model_summary["model"].eq(selected_model)
    model_summary["selection_penalty"] = 20.0 * model_summary["active_bound_count"].astype(float)
    model_summary["selection_penalty"] += 80.0 * np.maximum(model_summary["ea_retained_abs_relative_bias"].astype(float) - 0.35, 0.0)
    model_summary["selection_penalty"] += 40.0 * np.maximum(model_summary["ea_retained_relative_rmse"].astype(float) - 0.75, 0.0)
    model_summary["selection_penalty"] += 80.0 * np.maximum(model_summary["ea_total_abs_relative_bias"].astype(float) - 0.30, 0.0)
    model_summary["selection_penalty"] += 40.0 * np.maximum(model_summary["ea_total_relative_rmse"].astype(float) - 0.55, 0.0)
    model_summary["selection_score"] = model_summary["bic"] + model_summary["selection_penalty"]
    model_summary = model_summary.sort_values("selection_score").reset_index(drop=True)
    model_summary.to_csv(RESULTS_DIR / "model_selection_summary.csv", index=False)
    (RESULTS_DIR / "selected_model.txt").write_text(selected_model, encoding="utf-8")
    plot_model_comparison(model_summary, plot_dir)

    selected_variant = variants[selected_model]
    theta_selected = theta_by_model[selected_model]
    pd.Series(theta_selected).to_csv(RESULTS_DIR / "theta_selected_aroma_model.csv")
    pd.Series({name: theta_selected[name] for name in selected_variant.parameters}).to_csv(RESULTS_DIR / "theta_selected_aroma_only.csv")
    pred_selected = pred_by_model[selected_model]
    pred_selected.to_csv(RESULTS_DIR / "predictions_selected_model.csv", index=False)
    plot_selected_fit(pred_selected, selected_model, plot_dir / "selected_fit")

    print(f"[selected] {selected_model}", flush=True)
    target_parameters = CORE_TARGETS + SECONDARY_TARGETS + selected_variant.parameters
    (RESULTS_DIR / "target_parameters_selected_model.json").write_text(json.dumps(list(target_parameters), indent=2), encoding="utf-8")

    print("[fim] selected model current data", flush=True)
    jac, resid = finite_difference_jacobian_generic(
        theta_selected,
        target_parameters,
        lambda th: combined_residual_selected(th, primary_batches, extended_batches, co2_down, selected_variant),
        args.sensitivity_step,
    )
    fim = 0.5 * ((jac.T @ jac) + (jac.T @ jac).T)
    pd.DataFrame(jac, columns=target_parameters).to_csv(RESULTS_DIR / "jacobian_selected_current.csv", index=False)
    pd.Series(resid).to_csv(RESULTS_DIR / "residual_selected_current.csv", index=False)
    pd.DataFrame(fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_selected_current.csv")
    if rescale_FIM is not None:
        try:
            scaled = rescale_FIM(fim, np.asarray([theta_selected[p] for p in target_parameters], dtype=float))
            pd.DataFrame(scaled, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_selected_current_pyomo_rescaled.csv")
        except Exception as err:
            (RESULTS_DIR / "fim_selected_current_pyomo_rescaled_error.txt").write_text(str(err), encoding="utf-8")
    spectrum, weak, estimability = fim_diagnostics_generic(fim, theta_selected, target_parameters, "selected_current")
    spectrum.to_csv(RESULTS_DIR / "eigen_spectrum_selected_current.csv", index=False)
    weak.to_csv(RESULTS_DIR / "weak_directions_selected_current.csv", index=False)
    estimability.to_csv(RESULTS_DIR / "parameter_estimability_selected_current.csv", index=False)

    ranking = pd.DataFrame()
    selected_campaign = pd.DataFrame()
    if not args.skip_doe:
        print("[doe] selected model natural candidates", flush=True)
        designs = add_extra_designs(model_data, pilot.natural_candidate_designs(model_data))
        design_rows = []
        for design in designs.values():
            design_rows.append(
                {
                    "candidate": design.name,
                    "family": design.family,
                    "medium": design.medium,
                    "horizon_h": design.horizon_h,
                    "temperature_segments": ", ".join(f"{t:g}" for t in design.temperature_segments),
                    "N_pulses_kg_m3": "; ".join(f"{t:g}h:{a:g}" for t, a in design.pulses.get("N", tuple())),
                    "initials_json": json.dumps(design.initials, sort_keys=True),
                    "rationale": design.rationale,
                }
            )
        pd.DataFrame(design_rows).to_csv(RESULTS_DIR / "candidate_design_library_selected_model.csv", index=False)
        candidate_fims = {}
        for idx, (name, design) in enumerate(designs.items(), start=1):
            print(f"[doe] candidate {idx}/{len(designs)} {name}", flush=True)
            cfim = candidate_fim_selected(theta_selected, design, selected_variant, target_parameters, args.sensitivity_step)
            candidate_fims[name] = cfim
            pd.DataFrame(cfim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / f"candidate_fim_selected_{name}.csv")
        ranking = rank_candidates(candidate_fims, designs, fim, target_parameters)
        ranking.to_csv(RESULTS_DIR / "candidate_ranking_selected_model.csv", index=False)
        selected_campaign, campaign_fim = greedy_select(candidate_fims, designs, fim, target_parameters, args.campaign_size, objective="hybrid")
        selected_campaign.to_csv(RESULTS_DIR / "selected_campaign_hybrid_selected_model.csv", index=False)
        pd.DataFrame(campaign_fim, index=target_parameters, columns=target_parameters).to_csv(RESULTS_DIR / "fim_selected_current_plus_campaign_hybrid.csv")
        spec_after, weak_after, estim_after = fim_diagnostics_generic(campaign_fim, theta_selected, target_parameters, "selected_current_plus_campaign")
        spec_after.to_csv(RESULTS_DIR / "eigen_spectrum_selected_current_plus_campaign.csv", index=False)
        weak_after.to_csv(RESULTS_DIR / "weak_directions_selected_current_plus_campaign.csv", index=False)
        estim_after.to_csv(RESULTS_DIR / "parameter_estimability_selected_current_plus_campaign.csv", index=False)
        plot_design_inputs(selected_campaign, designs, plot_dir / "designs")

    write_report(variants, selected_model, model_summary, estimability, weak, ranking, selected_campaign)
    create_notebook(variants, selected_model)
    print(f"[done] selected model: {selected_model}", flush=True)
    print(f"[done] results: {RESULTS_DIR}", flush=True)
    print(f"[done] notebook: {NOTEBOOK_PATH}", flush=True)


if __name__ == "__main__":
    main()
