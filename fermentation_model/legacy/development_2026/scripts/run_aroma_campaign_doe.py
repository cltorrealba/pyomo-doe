from __future__ import annotations

import argparse
import json
import math
import os
import sys
import textwrap
from functools import lru_cache
from dataclasses import replace
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.environ as pyo
import pyomo.dae as dae
from pyomo.contrib.doe import DesignOfExperiments

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_extended_campaign_doe as base
import aroma_partition_unifac as unifac_partition
from run_fit_strategy_analysis import load_notebook_context

RESULTS_DIR = SCRIPT_DIR / "results" / "aroma_campaign_doe"

AROMA_SPECIES = ("ethyl_acetate", "isoamyl_acetate", "ethyl_octanoate")
AROMA_LABELS = {
    "ethyl_acetate": "Ethyl acetate",
    "isoamyl_acetate": "Isoamyl acetate",
    "ethyl_octanoate": "Ethyl octanoate",
}
AROMA_SHORT = {
    "ethyl_acetate": "EA",
    "isoamyl_acetate": "IAA",
    "ethyl_octanoate": "EO",
}

AROMA_SYNTHESIS_PARAMS = (
    "k_EA_growth",
    "k_EA_stationary",
    "k_IAA_growth",
    "k_IAA_stationary",
    "k_EO_growth",
    "k_EO_stationary",
)
AROMA_PARTITION_ALPHA_PARAMS = (
    "alpha_EA_partition",
    "alpha_IAA_partition",
    "alpha_EO_partition",
)
FERMENTATION_PARAMS = base.FULL15
JOINT_FULL_PARAMS = FERMENTATION_PARAMS + AROMA_SYNTHESIS_PARAMS

# Nominal synthesis yields in mg aroma / g sugar consumed.
# They are intentionally treated as weak literature-scale seeds for local DOE.
AROMA_NOMINAL_YIELDS = {
    "ethyl_acetate": {"growth": 0.10, "stationary": 0.28},
    "isoamyl_acetate": {"growth": 0.003, "stationary": 0.012},
    "ethyl_octanoate": {"growth": 0.0008, "stationary": 0.0030},
}

AROMA_YIELD_BOUNDS = {
    "ethyl_acetate": (1e-4, 2.0),
    "isoamyl_acetate": (1e-5, 0.20),
    "ethyl_octanoate": (1e-6, 0.08),
}

PARTITION_MODE = "water_ethanol_total_sugar_as_glucose"
PARTITION_IMPLEMENTATION = "surrogate"
AROMA_PARTITION_LITERATURE = unifac_partition.fit_all_partition_surrogates(PARTITION_MODE)
AROMA_VAPOR_PRESSURE_DATA = unifac_partition.vapor_pressure_correlation_data()


def set_partition_model(mode: str) -> dict[str, dict[str, float]]:
    global PARTITION_MODE, AROMA_PARTITION_LITERATURE
    PARTITION_MODE = str(mode)
    AROMA_PARTITION_LITERATURE = unifac_partition.fit_all_partition_surrogates(PARTITION_MODE)
    return AROMA_PARTITION_LITERATURE


def set_partition_implementation(implementation: str) -> str:
    global PARTITION_IMPLEMENTATION
    if implementation not in {"surrogate", "symbolic_unifac"}:
        raise ValueError(f"Unknown partition implementation {implementation!r}.")
    PARTITION_IMPLEMENTATION = implementation
    return PARTITION_IMPLEMENTATION


def set_results_dir(path: Path) -> Path:
    global RESULTS_DIR
    RESULTS_DIR = Path(path)
    return RESULTS_DIR

AROMA_LIQUID_ERROR = {
    "ethyl_acetate": 2.0,
    "isoamyl_acetate": 0.05,
    "ethyl_octanoate": 0.02,
}
AROMA_CONDENSATE_ERROR = {
    "ethyl_acetate": 2.0,
    "isoamyl_acetate": 0.05,
    "ethyl_octanoate": 0.02,
}

CO2_L_PER_G = 22.414 / 44.01
N_PHASE_HALF_KG_M3 = 0.035


def target_parameters(parameter_mode: str, include_partition_alpha: bool = False) -> tuple[str, ...]:
    if parameter_mode == "aroma_only":
        parameters = AROMA_SYNTHESIS_PARAMS
    elif parameter_mode == "joint_full":
        parameters = JOINT_FULL_PARAMS
    else:
        raise ValueError(f"Unknown parameter mode {parameter_mode!r}.")
    if include_partition_alpha:
        parameters = parameters + AROMA_PARTITION_ALPHA_PARAMS
    return tuple(parameters)


def parameter_component(m, parameter_name: str):
    if parameter_name in FERMENTATION_PARAMS:
        return getattr(m, parameter_name)
    mapping = {
        "k_EA_growth": m.k_aroma_growth["ethyl_acetate"],
        "k_EA_stationary": m.k_aroma_stationary["ethyl_acetate"],
        "k_IAA_growth": m.k_aroma_growth["isoamyl_acetate"],
        "k_IAA_stationary": m.k_aroma_stationary["isoamyl_acetate"],
        "k_EO_growth": m.k_aroma_growth["ethyl_octanoate"],
        "k_EO_stationary": m.k_aroma_stationary["ethyl_octanoate"],
        "alpha_EA_partition": m.alpha_partition["ethyl_acetate"],
        "alpha_IAA_partition": m.alpha_partition["isoamyl_acetate"],
        "alpha_EO_partition": m.alpha_partition["ethyl_octanoate"],
    }
    return mapping[parameter_name]


def literature_prior_fim(parameters: tuple[str, ...], prior_cv: float) -> np.ndarray:
    cv = max(float(prior_cv), 1e-6)
    return np.eye(len(parameters), dtype=float) * (1.0 / (cv * cv))


def prior_fim_for_parameters(parameters: tuple[str, ...], prior_cv: float) -> np.ndarray:
    prior = np.zeros((len(parameters), len(parameters)), dtype=float)
    fermentation_params = [name for name in parameters if name in FERMENTATION_PARAMS]
    if fermentation_params:
        prior_path = base.RESULTS_DIR / "extended_historical_prior_fim.csv"
        if not prior_path.exists():
            prior_path = base.RESULTS_DIR / "historical_prior_fim_full15_step0.01.csv"
        if not prior_path.exists():
            raise FileNotFoundError(
                "Historical fermentation FIM not found. Run run_extended_campaign_doe.py once before joint aroma DOE."
            )
        prior_hist = pd.read_csv(prior_path, index_col=0).loc[fermentation_params, fermentation_params].to_numpy(dtype=float)
        for i, row_name in enumerate(fermentation_params):
            ii = parameters.index(row_name)
            for j, col_name in enumerate(fermentation_params):
                jj = parameters.index(col_name)
                prior[ii, jj] = prior_hist[i, j]
    cv = max(float(prior_cv), 1e-6)
    aroma_weight = 1.0 / (cv * cv)
    for name in parameters:
        if name not in FERMENTATION_PARAMS:
            idx = parameters.index(name)
            prior[idx, idx] = aroma_weight
    return 0.5 * (prior + prior.T)


def _partition_value(species: str, temp_c: float, ethanol_g_l: float, sugar_g_l: float = 100.0) -> float:
    if PARTITION_IMPLEMENTATION == "symbolic_unifac":
        return float(
            unifac_partition.unifac_partition_K(
                species,
                temp_c,
                ethanol_g_l,
                glucose_g_l=sugar_g_l,
                fructose_g_l=0.0,
                mode=PARTITION_MODE,
            )["K"]
        )
    pars = AROMA_PARTITION_LITERATURE[species]
    exponent = (
        pars["temp_slope"] * (float(temp_c) - 20.0)
        + pars["ethanol_slope"] * (float(ethanol_g_l) - 50.0)
        + pars.get("sugar_slope", 0.0) * (float(sugar_g_l) - 100.0)
    )
    return float(pars["K20"] * math.exp(exponent))


@lru_cache(maxsize=None)
def _symbolic_unifac_data(species: str, mode: str) -> dict[str, object]:
    if species not in AROMA_SPECIES:
        raise ValueError(f"Unknown aroma species {species!r}.")
    aroma_name = unifac_partition.AROMA_SPECIES_TO_CHEMICAL[species]
    component_names = ["water", "ethanol"]
    if mode in {"water_ethanol_glucose", "water_ethanol_total_sugar_as_glucose"}:
        component_names.append("glucose")
    component_names.append(aroma_name)
    groups = {name: unifac_partition.unifac_groups(name) for name in component_names}
    group_ids = tuple(sorted({group for comp_groups in groups.values() for group in comp_groups}))
    constants = unifac_partition.unifac_subgroup_constants(group_ids)
    interactions = {}
    for group_from in group_ids:
        for group_to in group_ids:
            interactions[(group_from, group_to)] = unifac_partition.unifac_interaction_parameter(group_from, group_to)
    rs = {
        name: sum(float(constants[group]["R"]) * count for group, count in comp_groups.items())
        for name, comp_groups in groups.items()
    }
    qs = {
        name: sum(float(constants[group]["Q"]) * count for group, count in comp_groups.items())
        for name, comp_groups in groups.items()
    }
    group_count_totals = {name: sum(comp_groups.values()) for name, comp_groups in groups.items()}
    molecular_weights = {name: unifac_partition.molecular_weight(name) for name in ("water", "ethanol", "glucose")}
    return {
        "species": species,
        "aroma_name": aroma_name,
        "component_names": tuple(component_names),
        "groups": groups,
        "group_ids": group_ids,
        "constants": constants,
        "interactions": interactions,
        "rs": rs,
        "qs": qs,
        "group_count_totals": group_count_totals,
        "molecular_weights": molecular_weights,
    }


def _pyomo_vapor_pressure(species: str, temp_k):
    pars = AROMA_VAPOR_PRESSURE_DATA[species]
    if pars["type"] == "Antoine":
        return pyo.exp(math.log(float(pars["base"])) * (float(pars["A"]) - float(pars["B"]) / (temp_k + float(pars["C"]))))
    if pars["type"] == "Wagner_original":
        tr = temp_k / float(pars["Tc"])
        tau = 1.0 - tr
        return float(pars["Pc"]) * pyo.exp(
            (
                float(pars["a"]) * tau
                + float(pars["b"]) * tau**1.5
                + float(pars["c"]) * tau**3
                + float(pars["d"]) * tau**6
            )
            / tr
        )
    if pars["type"] == "Wagner":
        tr = temp_k / float(pars["Tc"])
        tau = 1.0 - tr
        tau_rt = tau**0.5
        tau15 = tau * tau_rt
        tau25 = tau * tau15
        return float(pars["Pc"]) * pyo.exp((float(pars["a"]) + float(pars["b"]) * tau_rt + tau15 * (float(pars["c"]) + float(pars["d"]) * tau25)) * tau / tr)
    raise ValueError(f"Unsupported vapor-pressure correlation {pars['type']!r}.")


def _pyomo_symbolic_unifac_K(species: str, temp_c, ethanol_g_l, glucose_g_l, fructose_g_l):
    data = _symbolic_unifac_data(species, PARTITION_MODE)
    temp_k = temp_c + 273.15
    sugar_total = glucose_g_l + fructose_g_l
    water_g_l = 1000.0 - ethanol_g_l - sugar_total
    moles = {
        "water": water_g_l / float(data["molecular_weights"]["water"]),
        "ethanol": ethanol_g_l / float(data["molecular_weights"]["ethanol"]),
    }
    if PARTITION_MODE == "water_ethanol_glucose":
        moles["glucose"] = glucose_g_l / float(data["molecular_weights"]["glucose"])
    elif PARTITION_MODE == "water_ethanol_total_sugar_as_glucose":
        moles["glucose"] = sugar_total / float(data["molecular_weights"]["glucose"])
    elif PARTITION_MODE != "water_ethanol":
        raise ValueError(f"Unknown partition mode {PARTITION_MODE!r}.")
    aroma_name = str(data["aroma_name"])
    moles[aroma_name] = unifac_partition.TRACE_MOLES_PER_L
    total_moles_l = sum(moles.values())
    xs = {name: moles[name] / total_moles_l for name in data["component_names"]}

    groups = data["groups"]
    group_ids = data["group_ids"]
    constants = data["constants"]
    interactions = data["interactions"]
    rs = data["rs"]
    qs = data["qs"]

    eps = 1e-20

    def psi(group_from: int, group_to: int):
        a_mn = interactions[(group_from, group_to)]
        return 1.0 if a_mn is None else pyo.exp(-float(a_mn) / temp_k)

    rsxs = sum(float(rs[name]) * xs[name] for name in data["component_names"])
    qsxs = sum(float(qs[name]) * xs[name] for name in data["component_names"])
    v_i = float(rs[aroma_name]) / rsxs
    f_i = float(qs[aroma_name]) / qsxs
    log_gamma_c = 1.0 - v_i + pyo.log(v_i + eps) - 5.0 * float(qs[aroma_name]) * (1.0 - v_i / f_i + pyo.log(v_i / f_i + eps))

    group_sum = sum(float(data["group_count_totals"][name]) * xs[name] for name in data["component_names"])
    group_x = {}
    for group in group_ids:
        numerator = sum(float(groups[name].get(group, 0)) * xs[name] for name in data["component_names"])
        group_x[group] = numerator / group_sum
    q_sum = sum(float(constants[group]["Q"]) * group_x[group] for group in group_ids)
    area = {group: float(constants[group]["Q"]) * group_x[group] / q_sum for group in group_ids}

    def log_gamma_group(group_k: int, area_fractions: dict[int, object], active_groups: tuple[int, ...]):
        sum1 = sum(area_fractions[group_m] * psi(group_m, group_k) for group_m in active_groups)
        sum2 = 0.0
        for group_m in active_groups:
            denom = sum(area_fractions[group_n] * psi(group_n, group_m) for group_n in active_groups)
            sum2 += area_fractions[group_m] * psi(group_k, group_m) / (denom + eps)
        return float(constants[group_k]["Q"]) * (1.0 - pyo.log(sum1 + eps) - sum2)

    mix_log_gamma_groups = {group: log_gamma_group(group, area, tuple(group_ids)) for group in group_ids}

    aroma_groups = groups[aroma_name]
    aroma_group_ids = tuple(sorted(aroma_groups))
    aroma_group_sum = float(sum(aroma_groups.values()))
    pure_group_x = {group: float(aroma_groups[group]) / aroma_group_sum for group in aroma_group_ids}
    pure_q_sum = sum(float(constants[group]["Q"]) * pure_group_x[group] for group in aroma_group_ids)
    pure_area = {group: float(constants[group]["Q"]) * pure_group_x[group] / pure_q_sum for group in aroma_group_ids}
    pure_log_gamma_groups = {group: log_gamma_group(group, pure_area, aroma_group_ids) for group in aroma_group_ids}

    log_gamma_r = sum(
        float(count) * (mix_log_gamma_groups[group] - pure_log_gamma_groups[group])
        for group, count in aroma_groups.items()
    )
    gamma = pyo.exp(log_gamma_c + log_gamma_r)
    psat = _pyomo_vapor_pressure(species, temp_k)
    return gamma * psat / (unifac_partition.R_PA_M3_MOL_K * temp_k * 1000.0 * total_moles_l)


def _candidate_time_grid(candidate: base.ExtendedDesignCandidate, liquid_interval_h: float, co2_interval_h: float) -> np.ndarray:
    return base._sample_grid(candidate.horizon_h, liquid_interval_h, co2_interval_h)


def _positive_sugar_uptake_from_guess(df: pd.DataFrame, candidate: base.ExtendedDesignCandidate) -> np.ndarray:
    time = df.index.to_numpy(dtype=float)
    g = df["G"].to_numpy(dtype=float)
    f = df["F"].to_numpy(dtype=float)
    dg = np.gradient(g, time, edge_order=1)
    dfdt = np.gradient(f, time, edge_order=1)
    g_pulse = np.array([base._continuous_pulse_rate(t, candidate.pulses.get("G", tuple())) for t in time])
    f_pulse = np.array([base._continuous_pulse_rate(t, candidate.pulses.get("F", tuple())) for t in time])
    return np.maximum(0.0, -(dg - g_pulse) - (dfdt - f_pulse))


def aroma_dynamic_guess(ns: dict, theta: dict[str, float], candidate: base.ExtendedDesignCandidate, time: np.ndarray) -> dict[str, pd.DataFrame]:
    base_guess = base._extended_dynamic_guess(ns, candidate, theta, time)
    if base_guess is None:
        zeros = pd.DataFrame(0.0, index=time, columns=list(AROMA_SPECIES))
        return {"liquid": zeros.copy(), "loss": zeros.copy(), "condensate": zeros.iloc[[0]].copy()}

    sugar_uptake = _positive_sugar_uptake_from_guess(base_guess, candidate)
    liquid = pd.DataFrame(0.0, index=time, columns=list(AROMA_SPECIES))
    loss = pd.DataFrame(0.0, index=time, columns=list(AROMA_SPECIES))
    for idx in range(1, len(time)):
        dt = max(float(time[idx] - time[idx - 1]), 1e-9)
        n_prev = max(float(base_guess["N"].iloc[idx - 1]), 0.0)
        phi_g = n_prev / (n_prev + N_PHASE_HALF_KG_M3)
        temp_c = float(base._temperature_profile(candidate, np.array([time[idx - 1]]))[0])
        ethanol = max(float(base_guess["E"].iloc[idx - 1]), 0.0)
        sugar = max(float(base_guess["G"].iloc[idx - 1]), 0.0) + max(float(base_guess["F"].iloc[idx - 1]), 0.0)
        co2_rate = max(float(np.gradient(base_guess["CO2"].to_numpy(dtype=float), time, edge_order=1)[idx - 1]), 0.0)
        qco2 = CO2_L_PER_G * co2_rate
        for species in AROMA_SPECIES:
            yields = AROMA_NOMINAL_YIELDS[species]
            k_syn = yields["growth"] * phi_g + yields["stationary"] * (1.0 - phi_g)
            r_syn = k_syn * sugar_uptake[idx - 1]
            k_part = _partition_value(species, temp_c, ethanol, sugar)
            r_loss = qco2 * k_part * max(float(liquid[species].iloc[idx - 1]), 0.0)
            liquid.loc[time[idx], species] = max(0.0, float(liquid[species].iloc[idx - 1]) + dt * (r_syn - r_loss))
            loss.loc[time[idx], species] = max(0.0, float(loss[species].iloc[idx - 1]) + dt * r_loss)
    condensate = pd.DataFrame(
        {species: [loss[species].iloc[-1] * AROMA_PARTITION_LITERATURE[species]["trap_efficiency"]] for species in AROMA_SPECIES},
        index=[time[-1]],
    )
    return {"liquid": liquid, "loss": loss, "condensate": condensate}


def add_aroma_layer(ns: dict, m, candidate: base.ExtendedDesignCandidate, theta: dict[str, float], liquid_interval_h: float, co2_interval_h: float):
    time = np.array([float(t) for t in m.t], dtype=float)
    guesses = aroma_dynamic_guess(ns, theta, candidate, time)
    m.aroma_species = pyo.Set(initialize=list(AROMA_SPECIES), ordered=True)
    liq_init = {(species, float(t)): float(guesses["liquid"].loc[float(t), species]) for species in AROMA_SPECIES for t in time}
    loss_init = {(species, float(t)): float(guesses["loss"].loc[float(t), species]) for species in AROMA_SPECIES for t in time}
    m.A_liq = pyo.Var(m.aroma_species, m.t, bounds=(0.0, 500.0), initialize=liq_init)
    m.A_loss = pyo.Var(m.aroma_species, m.t, bounds=(0.0, 500.0), initialize=loss_init)
    m.dA_liq = dae.DerivativeVar(m.A_liq, wrt=m.t)
    m.dA_loss = dae.DerivativeVar(m.A_loss, wrt=m.t)

    m.k_aroma_growth = pyo.Var(
        m.aroma_species,
        bounds=lambda _m, species: AROMA_YIELD_BOUNDS[str(species)],
        initialize=lambda _m, species: AROMA_NOMINAL_YIELDS[str(species)]["growth"],
    )
    m.k_aroma_stationary = pyo.Var(
        m.aroma_species,
        bounds=lambda _m, species: AROMA_YIELD_BOUNDS[str(species)],
        initialize=lambda _m, species: AROMA_NOMINAL_YIELDS[str(species)]["stationary"],
    )
    m.alpha_partition = pyo.Var(m.aroma_species, bounds=(0.2, 5.0), initialize=1.0)
    for species in m.aroma_species:
        m.k_aroma_growth[species].fix(AROMA_NOMINAL_YIELDS[str(species)]["growth"])
        m.k_aroma_stationary[species].fix(AROMA_NOMINAL_YIELDS[str(species)]["stationary"])
        m.alpha_partition[species].fix(1.0)

    t0 = float(m.t.first())
    tf = float(m.t.last())
    for species in m.aroma_species:
        m.A_liq[species, t0].fix(0.0)
        m.A_loss[species, t0].fix(0.0)

    eps = 1e-8
    m.aroma_phi_growth = pyo.Expression(m.t, rule=lambda _m, t: m.N[t] / (m.N[t] + N_PHASE_HALF_KG_M3 + eps))
    m.aroma_G_uptake = pyo.Expression(
        m.t,
        rule=lambda _m, t: (
            m.qXG_growth[t]
            + m.qEG_ferm[t]
            + m.maintenance[t] * m.maintenance_sugar_availability[t] * (m.G_eff[t] / (m.G_eff[t] + m.F_eff[t] + eps))
        )
        * m.X[t],
    )
    m.aroma_F_uptake = pyo.Expression(
        m.t,
        rule=lambda _m, t: (
            m.qXF_growth[t]
            + m.qEF_ferm[t]
            + m.maintenance[t] * m.maintenance_sugar_availability[t] * (m.F_eff[t] / (m.G_eff[t] + m.F_eff[t] + eps))
        )
        * m.X[t],
    )
    m.aroma_sugar_uptake = pyo.Expression(m.t, rule=lambda _m, t: m.aroma_G_uptake[t] + m.aroma_F_uptake[t])
    m.Q_CO2_vol = pyo.Expression(m.t, rule=lambda _m, t: CO2_L_PER_G * m.CO2_rate[t])

    def partition_rule(_m, species, t):
        if PARTITION_IMPLEMENTATION == "symbolic_unifac":
            return _pyomo_symbolic_unifac_K(
                str(species),
                m.TempC[t],
                m.E[t],
                m.G_eff[t],
                m.F_eff[t],
            )
        pars = AROMA_PARTITION_LITERATURE[str(species)]
        sugar = m.G_eff[t] + m.F_eff[t]
        return pars["K20"] * pyo.exp(
            pars["temp_slope"] * (m.TempC[t] - 20.0)
            + pars["ethanol_slope"] * (m.E[t] - 50.0)
            + pars.get("sugar_slope", 0.0) * (sugar - 100.0)
        )

    m.K_lg = pyo.Expression(m.aroma_species, m.t, rule=partition_rule)
    m.aroma_synthesis_rate = pyo.Expression(
        m.aroma_species,
        m.t,
        rule=lambda _m, species, t: (
            m.k_aroma_growth[species] * m.aroma_phi_growth[t]
            + m.k_aroma_stationary[species] * (1.0 - m.aroma_phi_growth[t])
        )
        * m.aroma_sugar_uptake[t],
    )
    m.aroma_loss_rate = pyo.Expression(
        m.aroma_species,
        m.t,
        rule=lambda _m, species, t: m.alpha_partition[species] * m.K_lg[species, t] * m.Q_CO2_vol[t] * m.A_liq[species, t],
    )

    @m.Constraint(m.aroma_species, m.t)
    def aroma_liquid_balance(_m, species, t):
        return m.dA_liq[species, t] == m.aroma_synthesis_rate[species, t] - m.aroma_loss_rate[species, t]

    @m.Constraint(m.aroma_species, m.t)
    def aroma_loss_balance(_m, species, t):
        return m.dA_loss[species, t] == m.aroma_loss_rate[species, t]

    m.A_cond = pyo.Var(
        m.aroma_species,
        bounds=(0.0, 500.0),
        initialize=lambda _m, species: float(guesses["condensate"].iloc[-1][str(species)]),
    )

    @m.Constraint(m.aroma_species)
    def aroma_condensate_measurement(_m, species):
        eta = AROMA_PARTITION_LITERATURE[str(species)]["trap_efficiency"]
        return m.A_cond[species] == eta * m.A_loss[species, tf]

    return m


def build_aroma_model(
    ns: dict,
    candidate: base.ExtendedDesignCandidate,
    theta: dict[str, float],
    target_parameters: tuple[str, ...],
    pulse_slots: int,
    liquid_interval_h: float,
    co2_interval_h: float,
):
    fermentation_targets = tuple(name for name in target_parameters if name in FERMENTATION_PARAMS)
    batch = base.make_synthetic_batch(ns, candidate, liquid_interval_h, co2_interval_h, theta)
    m = base.build_extended_zenteno_model(
        ns,
        batch,
        theta,
        parameters=fermentation_targets or base.ACCEPTED7,
        candidate=candidate,
        pulse_slots=pulse_slots,
        future_outputs=True,
    )
    add_aroma_layer(ns, m, candidate, theta, liquid_interval_h, co2_interval_h)
    base.discretize_model(m)
    label_aroma_model(
        m,
        candidate,
        target_parameters,
        liquid_interval_h,
        co2_interval_h,
        include_fermentation_outputs=bool(fermentation_targets),
    )
    return m


def label_aroma_model(
    m,
    candidate: base.ExtendedDesignCandidate,
    target_parameters: tuple[str, ...],
    aroma_liquid_interval_h: float,
    co2_interval_h: float,
    include_fermentation_outputs: bool,
):
    m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for parameter in target_parameters:
        component = parameter_component(m, parameter)
        m.unknown_parameters[component] = pyo.value(component)

    m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for s in m.temperature_segments:
        m.experiment_inputs[m.T_set[s]] = None
    for name in ("X0", "N0", "G0", "F0", "E0", "Xd0"):
        m.experiment_inputs[getattr(m, name)] = None
    for channel in base.CHANNELS:
        for p in getattr(m, f"{channel}_pulse_slots"):
            m.experiment_inputs[getattr(m, f"{channel}_pulse_time")[p]] = None
            m.experiment_inputs[getattr(m, f"{channel}_pulse_amount")[p]] = None

    m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)

    if include_fermentation_outputs:
        liquid_times = np.arange(0.0, float(candidate.horizon_h) + 1e-9, float(aroma_liquid_interval_h))
        co2_times = np.arange(0.0, float(candidate.horizon_h) + 1e-9, float(co2_interval_h))
        state_vars = {"X": m.X, "Xd": m.Xd, "N": m.N, "G": m.G, "F": m.F, "E": m.E}
        for state, var in state_vars.items():
            for t in liquid_times:
                mt = base._time_in_set(m, float(t))
                if mt is None:
                    continue
                m.experiment_outputs[var[mt]] = 0.0
                m.measurement_error[var[mt]] = base.FUTURE_MEASUREMENT_ERROR[state]
        for t in co2_times:
            mt = base._time_in_set(m, float(t))
            if mt is None:
                continue
            m.experiment_outputs[m.CO2[mt]] = 0.0
            m.measurement_error[m.CO2[mt]] = base.FUTURE_MEASUREMENT_ERROR["CO2"]

    aroma_times = np.arange(0.0, float(candidate.horizon_h) + 1e-9, float(aroma_liquid_interval_h))
    for species in AROMA_SPECIES:
        for t in aroma_times:
            mt = base._time_in_set(m, float(t))
            if mt is None:
                continue
            m.experiment_outputs[m.A_liq[species, mt]] = 0.0
            m.measurement_error[m.A_liq[species, mt]] = AROMA_LIQUID_ERROR[species]
        m.experiment_outputs[m.A_cond[species]] = 0.0
        m.measurement_error[m.A_cond[species]] = AROMA_CONDENSATE_ERROR[species]
    return m


class AromaExperiment:
    def __init__(
        self,
        ns: dict,
        candidate: base.ExtendedDesignCandidate,
        theta: dict[str, float],
        target_parameters: tuple[str, ...],
        pulse_slots: int,
        liquid_interval_h: float,
        co2_interval_h: float,
    ):
        self.ns = ns
        self.candidate = candidate
        self.theta = dict(theta)
        self.target_parameters = tuple(target_parameters)
        self.pulse_slots = int(pulse_slots)
        self.liquid_interval_h = float(liquid_interval_h)
        self.co2_interval_h = float(co2_interval_h)

    def get_labeled_model(self):
        return build_aroma_model(
            self.ns,
            self.candidate,
            self.theta,
            self.target_parameters,
            self.pulse_slots,
            self.liquid_interval_h,
            self.co2_interval_h,
        )


def fim_for_aroma_experiment(experiment: AromaExperiment, parameters: tuple[str, ...], step: float, solver) -> np.ndarray:
    doe = DesignOfExperiments(
        experiment=experiment,
        step=float(step),
        scale_nominal_param_value=True,
        solver=solver,
        tee=False,
    )
    fim = np.asarray(doe.compute_FIM(method="sequential"), dtype=float)
    if fim.shape != (len(parameters), len(parameters)):
        raise RuntimeError(f"Unexpected FIM shape {fim.shape}; expected {(len(parameters), len(parameters))}.")
    return 0.5 * (fim + fim.T)


def stable_inverse(fim: np.ndarray, ridge_fraction: float = 1e-9) -> np.ndarray:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    scale = max(float(np.trace(fim)) / max(fim.shape[0], 1), 1.0)
    return np.linalg.pinv(fim + ridge_fraction * scale * np.eye(fim.shape[0]))


def fim_metrics(fim: np.ndarray, parameters: tuple[str, ...], prefix: str = "") -> dict[str, float]:
    fim = 0.5 * (np.asarray(fim, dtype=float) + np.asarray(fim, dtype=float).T)
    eig = np.linalg.eigvalsh(fim)
    max_eig = float(np.max(eig)) if eig.size else np.nan
    floor = max(max_eig * 1e-12, np.finfo(float).tiny) if np.isfinite(max_eig) and max_eig > 0 else np.finfo(float).tiny
    eig_pos = np.clip(eig, floor, None)
    cov = stable_inverse(fim)
    out = {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if max_eig > 0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace": float(np.trace(fim)),
        f"{prefix}trace_inv": float(np.trace(cov)),
    }
    for idx, name in enumerate(parameters):
        out[f"{prefix}var_{name}"] = float(cov[idx, idx])
    return out


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    combined_cov = stable_inverse(combined)
    ratios = []
    rows = {}
    for idx, name in enumerate(parameters):
        before = float(prior_cov[idx, idx])
        after = float(combined_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        ratios.append(ratio)
    rows["mean_var_reduction"] = float(1.0 - np.nanmean(ratios))
    rows["worst_var_reduction"] = float(1.0 - np.nanmax(ratios))
    return rows


def campaign_score(prior_fim: np.ndarray, campaign_fim: np.ndarray, parameters: tuple[str, ...]) -> float:
    metrics = fim_metrics(campaign_fim, parameters)
    reductions = variance_reduction(prior_fim, campaign_fim, parameters)
    return (
        float(metrics["logdet"])
        + 10.0 * float(reductions["mean_var_reduction"])
        + 20.0 * float(reductions["worst_var_reduction"])
        + 2.0 * math.log(max(float(metrics["min_relative_eigenvalue"]), 1e-16))
    )


def aroma_candidate_library() -> list[base.ExtendedDesignCandidate]:
    candidates = []
    base_by_name = {c.name: c for c in base.candidate_library()}
    for name in (
        "glucose_rich_growth_yield",
        "fructose_rich_iG_probe",
        "combined_stress_long_horizon",
        "yan_saturation_scan",
        "low_temp_growth_separation_plus_co2",
        "glucose_pulse_after_N_depletion",
        "high_biomass_low_N_maintenance",
        "ethanol_initial_challenge",
        "viable_biomass_step",
    ):
        if name in base_by_name:
            candidates.append(base_by_name[name])

    def pulses(**kwargs):
        return {channel: tuple(kwargs.get(channel, tuple())) for channel in base.CHANNELS}

    candidates.extend(
        [
            base.ExtendedDesignCandidate(
                "cold_synthesis_hot_stripping",
                "aroma_partition",
                192.0,
                {"X": 0.50, "N": 0.24, "G": 95.0, "F": 95.0, "E": 0.0, "Xd": 0.0},
                (15.0, 16.0, 24.0, 25.0),
                pulses(N=((0.0, 0.08), (48.0, 0.05)), G=((96.0, 35.0),)),
                "Low-temperature synthesis window followed by high-temperature stripping challenge.",
            ),
            base.ExtendedDesignCandidate(
                "hot_synthesis_cold_retention",
                "aroma_partition",
                192.0,
                {"X": 0.50, "N": 0.24, "G": 95.0, "F": 95.0, "E": 0.0, "Xd": 0.0},
                (24.0, 25.0, 17.0, 15.0),
                pulses(N=((0.0, 0.08), (48.0, 0.05)), F=((96.0, 35.0),)),
                "High-temperature production window followed by cold retention phase.",
            ),
            base.ExtendedDesignCandidate(
                "ethanol_partition_late_pulse",
                "aroma_partition",
                216.0,
                {"X": 0.65, "N": 0.20, "G": 85.0, "F": 85.0, "E": 10.0, "Xd": 0.0},
                (18.0, 20.0, 24.0, 22.0),
                pulses(N=((0.0, 0.06),), E=((108.0, 25.0),), G=((72.0, 30.0),)),
                "Late ethanol pulse perturbs partition without changing aroma synthesis parameters directly.",
            ),
            base.ExtendedDesignCandidate(
                "co2_stripping_without_growth",
                "aroma_partition",
                192.0,
                {"X": 1.20, "N": 0.030, "G": 55.0, "F": 55.0, "E": 20.0, "Xd": 0.0},
                (16.0, 20.0, 25.0, 25.0),
                pulses(G=((72.0, 45.0),), F=((96.0, 30.0),), N=((120.0, 0.015),)),
                "Sugar and temperature pulses generate CO2 stripping under weak growth pressure.",
            ),
        ]
    )
    return candidates


def evaluate_candidate(ns, theta, candidate, parameters, prior_fim, step, solver, pulse_slots, liquid_interval_h, co2_interval_h):
    try:
        experiment = AromaExperiment(ns, candidate, theta, parameters, pulse_slots, liquid_interval_h, co2_interval_h)
        fim = fim_for_aroma_experiment(experiment, parameters, step, solver)
        combined = prior_fim + fim
        row = {
            "candidate": candidate.name,
            "family": candidate.family,
            "horizon_h": candidate.horizon_h,
            "temperature_c": ", ".join(f"{v:g}" for v in candidate.temperature_c),
            "status": "ok",
            "error": "",
            "score": campaign_score(prior_fim, combined, parameters),
            **fim_metrics(fim, parameters, prefix="new_"),
            **fim_metrics(combined, parameters, prefix="combined_"),
            **variance_reduction(prior_fim, combined, parameters),
        }
        return row, fim
    except Exception as err:
        row = {
            "candidate": candidate.name,
            "family": candidate.family,
            "horizon_h": candidate.horizon_h,
            "temperature_c": ", ".join(f"{v:g}" for v in candidate.temperature_c),
            "status": "failed",
            "error": f"{type(err).__name__}: {err}",
        }
        return row, None


def greedy_select(fims: dict[str, np.ndarray], prior_fim: np.ndarray, candidates: dict[str, base.ExtendedDesignCandidate], parameters: tuple[str, ...], campaign_size: int):
    selected = []
    remaining = set(fims)
    current = prior_fim.copy()
    rows = []
    for order in range(1, int(campaign_size) + 1):
        best_name = None
        best_score = -np.inf
        best_payload = None
        for name in sorted(remaining):
            trial = current + fims[name]
            score = campaign_score(prior_fim, trial, parameters)
            if score > best_score:
                best_name = name
                best_score = score
                best_payload = {
                    "campaign_score": score,
                    **fim_metrics(trial, parameters, prefix="campaign_"),
                    **variance_reduction(prior_fim, trial, parameters),
                }
        if best_name is None:
            break
        selected.append(best_name)
        remaining.remove(best_name)
        current = current + fims[best_name]
        c = candidates[best_name]
        rows.append(
            {
                "campaign_order": order,
                "candidate": best_name,
                "family": c.family,
                "horizon_h": c.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                "rationale": c.rationale,
                **(best_payload or {}),
            }
        )
    return pd.DataFrame(rows), current


def simulate_candidate(ns, theta, candidate, target_parameters, pulse_slots, liquid_interval_h, co2_interval_h):
    m = build_aroma_model(ns, candidate, theta, target_parameters, pulse_slots, liquid_interval_h, co2_interval_h)
    solver = base.make_solver(ns)
    result = solver.solve(m, tee=False)
    rows = []
    for t in m.t:
        for species in AROMA_SPECIES:
            rows.append(
                {
                    "candidate": candidate.name,
                    "family": candidate.family,
                    "t": float(t),
                    "species": species,
                    "temperature_c": pyo.value(m.TempC[t]),
                    "X": pyo.value(m.X[t]),
                    "N": pyo.value(m.N[t]),
                    "G": pyo.value(m.G[t]),
                    "F": pyo.value(m.F[t]),
                    "E": pyo.value(m.E[t]),
                    "CO2": pyo.value(m.CO2[t]),
                    "CO2_rate": pyo.value(m.CO2_rate[t]),
                    "Q_CO2_vol": pyo.value(m.Q_CO2_vol[t]),
                    "A_liq": pyo.value(m.A_liq[species, t]),
                    "A_loss": pyo.value(m.A_loss[species, t]),
                    "A_synthesis_rate": pyo.value(m.aroma_synthesis_rate[species, t]),
                    "A_loss_rate": pyo.value(m.aroma_loss_rate[species, t]),
                    "K_lg": pyo.value(m.K_lg[species, t]),
                    "A_cond_final": pyo.value(m.A_cond[species]),
                    "termination": str(result.solver.termination_condition),
                }
            )
    return pd.DataFrame(rows)


def pulse_table(candidates: list[base.ExtendedDesignCandidate]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        for channel in base.CHANNELS:
            for time_h, amount in candidate.pulses.get(channel, tuple()):
                if float(amount) <= 0.0:
                    continue
                rows.append(
                    {
                        "candidate": candidate.name,
                        "family": candidate.family,
                        "channel": channel,
                        "time_h": float(time_h),
                        "amount": float(amount),
                    }
                )
    return pd.DataFrame(rows)


def write_input_plots(candidates: list[base.ExtendedDesignCandidate], selected: pd.DataFrame):
    for stale in RESULTS_DIR.glob("aroma_campaign_inputs_*.png"):
        stale.unlink()
    by_name = {candidate.name: candidate for candidate in candidates}
    rows = [by_name[name] for name in selected["candidate"] if name in by_name]
    pulse_df = pulse_table(rows)
    for candidate in rows:
        fig, axes = plt.subplots(7, 1, figsize=(10, 12), sharex=False)
        time = np.linspace(0.0, candidate.horizon_h, 400)
        axes[0].step(time, base._temperature_profile(candidate, time), where="post", color="tab:red")
        axes[0].set_ylabel("T C")
        axes[0].set_xlim(0.0, candidate.horizon_h)
        init_labels = ["X", "N", "G", "F", "E"]
        init_values = [candidate.initials.get(name, 0.0) for name in init_labels]
        axes[1].bar(init_labels, init_values, color=["tab:purple", "tab:brown", "tab:blue", "tab:orange", "tab:green"])
        axes[1].set_ylabel("initial")
        for ax, channel in zip(axes[2:], base.CHANNELS):
            rows_ch = pulse_df[(pulse_df["candidate"].eq(candidate.name)) & (pulse_df["channel"].eq(channel))]
            if rows_ch.empty:
                ax.axhline(0.0, color="0.8", lw=1)
            else:
                ax.vlines(rows_ch["time_h"], 0.0, rows_ch["amount"], lw=3)
                ax.scatter(rows_ch["time_h"], rows_ch["amount"], s=35)
            ax.set_ylabel(channel)
            ax.set_xlim(0.0, candidate.horizon_h)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("t [h]")
        fig.suptitle("\n".join(textwrap.wrap(f"{candidate.name}: {candidate.rationale}", width=95)), fontsize=10, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(RESULTS_DIR / f"aroma_campaign_inputs_{candidate.name}.png", dpi=170)
        plt.close(fig)


def write_report(selected: pd.DataFrame, ranking: pd.DataFrame, parameter_reduction: pd.DataFrame, parameters: tuple[str, ...]):
    selected_cols = [
        "campaign_order",
        "candidate",
        "family",
        "horizon_h",
        "temperature_c",
        "campaign_logdet",
        "campaign_min_relative_eigenvalue",
        "campaign_condition_number",
        "mean_var_reduction",
        "worst_var_reduction",
        "rationale",
    ]
    lines = [
        "# Aroma campaign DOE",
        "",
        "This report evaluates model-based candidate experiments for aroma and fermentation parameters.",
        "",
        f"The gas-liquid loss model uses `{PARTITION_IMPLEMENTATION}` with partition mode `{PARTITION_MODE}`.",
        "The target vector is:",
        "",
        ", ".join(f"`{name}`" for name in parameters),
        "",
        "Liquid aroma observations are included at the liquid sampling interval. A single terminal condensate observation per aroma is included as a final accumulated gas-loss constraint.",
        "",
        "## Selected Campaign",
        "",
        selected[selected_cols].to_markdown(index=False),
        "",
        "## Candidate Ranking",
        "",
        ranking.head(15)[
            [
                "candidate",
                "family",
                "horizon_h",
                "combined_logdet",
                "combined_min_relative_eigenvalue",
                "combined_condition_number",
                "mean_var_reduction",
                "worst_var_reduction",
                "status",
                "error",
            ]
        ].to_markdown(index=False),
        "",
        "## Parameter Variance Reductions",
        "",
        parameter_reduction.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- The terminal condensate improves the mass closure but does not provide time-resolved gas-loss dynamics.",
        "- Partition parameters are fixed by default to avoid confounding gas-liquid equilibrium with trap efficiency.",
        "- Temperature switches are intentionally included to separate synthesis windows from stripping windows.",
        "- Nitrogen ladders and sugar pulses are retained because the synthesis model is phase-dependent through nitrogen limitation and sugar uptake.",
    ]
    (RESULTS_DIR / "aroma_campaign_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Aroma-extended fermentation DOE with Pyomo DoE FIM screening.")
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--max-candidates", type=int, default=0)
    parser.add_argument("--step", type=float, default=1e-2)
    parser.add_argument("--pulse-slots", type=int, default=3)
    parser.add_argument("--liquid-interval-h", type=float, default=12.0)
    parser.add_argument("--co2-interval-h", type=float, default=6.0)
    parser.add_argument("--prior-cv", type=float, default=2.0)
    parser.add_argument("--parameter-mode", choices=["aroma_only", "joint_full"], default="aroma_only")
    parser.add_argument("--include-partition-alpha", action="store_true")
    parser.add_argument(
        "--partition-mode",
        choices=["water_ethanol", "water_ethanol_glucose", "water_ethanol_total_sugar_as_glucose"],
        default=PARTITION_MODE,
    )
    parser.add_argument(
        "--partition-implementation",
        choices=["surrogate", "symbolic_unifac"],
        default=PARTITION_IMPLEMENTATION,
        help="Use the fitted log-linear UNIFAC surrogate or the full symbolic original-UNIFAC equations in Pyomo.",
    )
    parser.add_argument("--results-dir", default="")
    parser.add_argument("--skip-simulations", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.results_dir:
        set_results_dir(Path(args.results_dir))
    elif args.parameter_mode == "joint_full":
        set_results_dir(SCRIPT_DIR / "results" / "aroma_joint_campaign_doe")
    else:
        set_results_dir(SCRIPT_DIR / "results" / "aroma_campaign_doe")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    set_partition_model(args.partition_mode)
    set_partition_implementation(args.partition_implementation)
    unifac_partition.write_partition_report(RESULTS_DIR / "aroma_unifac_partition_model.json", mode=args.partition_mode)
    ns = load_notebook_context()
    theta = base.complete_final_theta(ns)
    solver = base.make_solver(ns)
    parameters = target_parameters(args.parameter_mode, bool(args.include_partition_alpha))
    prior_fim = prior_fim_for_parameters(parameters, float(args.prior_cv))
    pd.DataFrame(prior_fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / "aroma_prior_fim.csv")

    candidates = aroma_candidate_library()
    if args.max_candidates and args.max_candidates > 0:
        candidates = candidates[: int(args.max_candidates)]
    candidate_by_name = {candidate.name: candidate for candidate in candidates}

    rows = []
    fims = {}
    for idx, candidate in enumerate(candidates, start=1):
        print(f"[aroma candidate {idx}/{len(candidates)}] {candidate.name}", flush=True)
        row, fim = evaluate_candidate(
            ns,
            theta,
            candidate,
            parameters,
            prior_fim,
            float(args.step),
            solver,
            int(args.pulse_slots),
            float(args.liquid_interval_h),
            float(args.co2_interval_h),
        )
        rows.append(row)
        if fim is not None:
            fims[candidate.name] = fim
            pd.DataFrame(fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / f"fim_aroma_{candidate.name}.csv")

    ranking = pd.DataFrame(rows)
    if not ranking.empty:
        ranking["solve_ok"] = ranking["status"].eq("ok")
        ranking = ranking.sort_values(["solve_ok", "score", "combined_logdet"], ascending=[False, False, False])
    ranking.to_csv(RESULTS_DIR / "aroma_candidate_ranking.csv", index=False)

    if not fims:
        print("No successful aroma candidate FIMs.")
        return 1

    selected, campaign_fim = greedy_select(fims, prior_fim, candidate_by_name, parameters, min(int(args.campaign_size), len(fims)))
    selected.to_csv(RESULTS_DIR / "aroma_campaign_selected.csv", index=False)
    pd.DataFrame(campaign_fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / "aroma_campaign_fim.csv")

    reductions = variance_reduction(prior_fim, campaign_fim, parameters)
    parameter_reduction = pd.DataFrame(
        [
            {
                "parameter": name,
                "group": "fermentation" if name in FERMENTATION_PARAMS else ("synthesis" if name in AROMA_SYNTHESIS_PARAMS else "partition_correction"),
                "campaign_var_ratio": reductions.get(f"var_ratio_{name}", np.nan),
                "campaign_var_reduction": reductions.get(f"var_reduction_{name}", np.nan),
            }
            for name in parameters
        ]
    )
    parameter_reduction.to_csv(RESULTS_DIR / "aroma_campaign_parameter_reduction.csv", index=False)

    selected_candidates = [candidate_by_name[name] for name in selected["candidate"]]
    pulse_table(selected_candidates).to_csv(RESULTS_DIR / "aroma_campaign_pulses.csv", index=False)
    pd.DataFrame(
        [
            {
                "candidate": candidate.name,
                "family": candidate.family,
                "horizon_h": candidate.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in candidate.temperature_c),
                **{f"{key}0": value for key, value in candidate.initials.items()},
                "rationale": candidate.rationale,
            }
            for candidate in selected_candidates
        ]
    ).to_csv(RESULTS_DIR / "aroma_campaign_protocols.csv", index=False)

    write_input_plots(candidates, selected)
    if not args.skip_simulations:
        sim_rows = []
        for candidate in selected_candidates:
            try:
                sim_rows.append(
                    simulate_candidate(
                        ns,
                        theta,
                        candidate,
                        parameters,
                        int(args.pulse_slots),
                        float(args.liquid_interval_h),
                        float(args.co2_interval_h),
                    )
                )
            except Exception as err:
                print(f"[aroma simulation failed] {candidate.name}: {type(err).__name__}: {err}", flush=True)
        if sim_rows:
            pd.concat(sim_rows, ignore_index=True).to_csv(RESULTS_DIR / "aroma_campaign_simulations.csv", index=False)

    metadata = {
        "campaign_size": int(args.campaign_size),
        "step": float(args.step),
        "pulse_slots": int(args.pulse_slots),
        "liquid_interval_h": float(args.liquid_interval_h),
        "co2_interval_h": float(args.co2_interval_h),
        "prior_cv": float(args.prior_cv),
        "parameter_mode": str(args.parameter_mode),
        "include_partition_alpha": bool(args.include_partition_alpha),
        "parameters": list(parameters),
        "n_candidates": len(candidates),
        "n_successful_candidates": len(fims),
        "partition_model": AROMA_PARTITION_LITERATURE,
        "partition_mode": PARTITION_MODE,
        "partition_implementation": PARTITION_IMPLEMENTATION,
        "vapor_pressure_correlations": AROMA_VAPOR_PRESSURE_DATA,
        "unifac_assignments": unifac_partition.available_unifac_assignments(),
        "nominal_yields": AROMA_NOMINAL_YIELDS,
    }
    (RESULTS_DIR / "aroma_campaign_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_report(selected, ranking, parameter_reduction, parameters)

    print("\nSelected aroma campaign:")
    print(
        selected[
            [
                "campaign_order",
                "candidate",
                "family",
                "horizon_h",
                "mean_var_reduction",
                "worst_var_reduction",
                "campaign_condition_number",
            ]
        ].to_string(index=False)
    )
    print("\nParameter reductions:")
    print(parameter_reduction.to_string(index=False))
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
