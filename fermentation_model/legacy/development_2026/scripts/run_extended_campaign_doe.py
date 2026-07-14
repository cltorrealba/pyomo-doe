from __future__ import annotations

import argparse
import json
import math
import os
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyomo.dae as dae
import pyomo.environ as pyo
from scipy.integrate import solve_ivp
from pyomo.contrib.doe import DesignOfExperiments

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_fit_strategy_analysis import load_notebook_context, series_to_theta

RESULTS_DIR = SCRIPT_DIR / "results" / "extended_campaign_doe"
FINAL_THETA_PATH = SCRIPT_DIR / "results" / "identifiability_reduction" / "theta_final_identifiable.csv"

ACCEPTED7 = ("mu0", "qN", "betaG0", "betaF0", "qEG", "qEF", "iG")
FULL15 = ("mu0", "sN", "qN", "qXG", "qXF", "betaG0", "sG", "betaF0", "sF", "qEG", "qEF", "iG", "iE", "Kd0", "m0")
CORRECTABLE = ("qXG", "qXF", "sN", "sG", "sF", "iE", "Kd0", "m0")
FUTURE_LIQUID_OUTPUTS = ("X", "Xd", "N", "G", "F", "E")
FUTURE_CO2_OUTPUTS = ("CO2",)
CO2_PER_ETHANOL_MASS = 44.01 / 46.07

CHANNELS = ("N", "G", "F", "E", "X")
CHANNEL_AMOUNT_BOUNDS = {
    "N": (0.0, 0.12),   # kg N/m3 per pulse
    "G": (0.0, 70.0),   # kg glucose/m3 per pulse
    "F": (0.0, 70.0),   # kg fructose/m3 per pulse
    "E": (0.0, 35.0),   # kg ethanol/m3 per pulse
    "X": (0.0, 3.0),    # kg viable biomass/m3 per pulse
}
CHANNEL_TOTAL_LIMITS = {
    "N": 0.24,
    "G": 90.0,
    "F": 90.0,
    "E": 45.0,
    "X": 4.0,
}

EXTENDED_STATE_BOUNDS = {
    "Xd": (0.0, 30.0),
    "CO2": (0.0, 350.0),
}
FUTURE_MEASUREMENT_ERROR = {
    "X": 0.50,
    "Xd": 0.50,
    "N": 0.01,
    "G": 2.0,
    "F": 2.0,
    "E": 1.0,
    "CO2": 1.0,
}


@dataclass(frozen=True)
class ExtendedDesignCandidate:
    name: str
    family: str
    horizon_h: float
    initials: dict[str, float]
    temperature_c: tuple[float, float, float, float]
    pulses: dict[str, tuple[tuple[float, float], ...]]
    rationale: str


class HistoricalExtendedExperiment:
    def __init__(self, ns: dict, batch_id: str, theta: dict[str, float], parameters: tuple[str, ...]):
        self.ns = ns
        self.batch_id = str(batch_id)
        self.theta = dict(theta)
        self.parameters = tuple(parameters)

    def get_labeled_model(self):
        batch = self.ns["load_batch"](self.batch_id)
        m = build_extended_zenteno_model(
            self.ns,
            batch,
            self.theta,
            parameters=self.parameters,
            candidate=None,
            future_outputs=False,
        )
        discretize_model(m)
        label_extended_model(
            self.ns,
            m,
            batch,
            self.parameters,
            candidate=None,
            future_outputs=False,
        )
        return m


class FutureExtendedExperiment:
    def __init__(
        self,
        ns: dict,
        candidate: ExtendedDesignCandidate,
        theta: dict[str, float],
        parameters: tuple[str, ...],
        pulse_slots: int,
        liquid_interval_h: float,
        co2_interval_h: float,
    ):
        self.ns = ns
        self.candidate = candidate
        self.theta = dict(theta)
        self.parameters = tuple(parameters)
        self.pulse_slots = int(pulse_slots)
        self.liquid_interval_h = float(liquid_interval_h)
        self.co2_interval_h = float(co2_interval_h)

    def get_labeled_model(self):
        batch = make_synthetic_batch(self.ns, self.candidate, self.liquid_interval_h, self.co2_interval_h, self.theta)
        m = build_extended_zenteno_model(
            self.ns,
            batch,
            self.theta,
            parameters=self.parameters,
            candidate=self.candidate,
            pulse_slots=self.pulse_slots,
            future_outputs=True,
        )
        discretize_model(m)
        label_extended_model(
            self.ns,
            m,
            batch,
            self.parameters,
            candidate=self.candidate,
            future_outputs=True,
            liquid_interval_h=self.liquid_interval_h,
            co2_interval_h=self.co2_interval_h,
        )
        return m


def complete_final_theta(ns: dict, batch_id: str = "25026") -> dict[str, float]:
    theta = ns["_theta_for_batch"](ns["load_batch"](batch_id), ns["DEFAULT_THETA"])
    if FINAL_THETA_PATH.exists():
        theta.update(series_to_theta(pd.read_csv(FINAL_THETA_PATH, index_col=0).iloc[:, 0]))
    clipped, _ = ns["clip_theta_to_bounds"](theta)
    return {name: float(clipped[name]) for name in ns["DEFAULT_THETA"]}


def make_solver(ns: dict):
    solver = pyo.SolverFactory("ipopt", executable=ns["IPOPT_EXECUTABLE"])
    for key, value in ns["PARMEST_SOLVER_OPTIONS"].items():
        solver.options[key] = value
    solver.options["max_iter"] = max(int(solver.options.get("max_iter", 0) or 0), 4000)
    solver.options["tol"] = min(float(solver.options.get("tol", 1e-6)), 1e-6)
    return solver


def _temperature_profile(candidate: ExtendedDesignCandidate, time: np.ndarray) -> np.ndarray:
    edges = np.linspace(0.0, float(candidate.horizon_h), len(candidate.temperature_c) + 1)
    values = []
    for t in time:
        idx = int(np.searchsorted(edges[1:-1], float(t), side="right"))
        values.append(float(candidate.temperature_c[min(idx, len(candidate.temperature_c) - 1)]))
    return np.asarray(values, dtype=float)


def _sample_grid(horizon_h: float, liquid_interval_h: float, co2_interval_h: float) -> np.ndarray:
    points = {0.0, float(horizon_h)}
    for interval in (liquid_interval_h, co2_interval_h, 3.0):
        n = int(math.ceil(float(horizon_h) / float(interval)))
        for k in range(n + 1):
            t = min(float(horizon_h), k * float(interval))
            points.add(round(t, 8))
    return np.array(sorted(points), dtype=float)


def _continuous_pulse_rate(t: float, schedule: tuple[tuple[float, float], ...], width_h: float = 3.0) -> float:
    width_h = max(float(width_h), 1e-6)
    rate = 0.0
    for pulse_time, amount in schedule:
        amount = float(amount)
        if amount <= 0.0:
            continue
        rate += amount * math.exp(-((float(t) - float(pulse_time)) / width_h) ** 2) / (math.sqrt(math.pi) * width_h)
    return rate


def _extended_dynamic_guess(ns: dict, candidate: ExtendedDesignCandidate, theta: dict[str, float], time: np.ndarray) -> pd.DataFrame | None:
    constants = ns["FIXED_CONSTANTS"]
    theta = {name: float(theta[name]) for name in ns["DEFAULT_THETA"]}
    y0 = np.array(
        [
            float(candidate.initials.get("X", 0.6)),
            float(candidate.initials.get("Xd", 0.0)),
            float(candidate.initials.get("N", 0.18)),
            float(candidate.initials.get("G", 80.0)),
            float(candidate.initials.get("F", 80.0)),
            float(candidate.initials.get("E", 0.0)),
            0.0,
        ],
        dtype=float,
    )

    def rhs(t, y):
        X, Xd, N, G, F, E, CO2 = [max(0.0, float(v)) for v in y]
        temp_c = float(_temperature_profile(candidate, np.array([float(t)]))[0])
        T_K = temp_c + 273.15
        R = constants["R"]
        eps = 1e-8
        A_mu = math.exp(constants["Eac"] * (T_K - 300.0) / (300.0 * R * T_K))
        A_beta = math.exp(constants["Eafe"] * (T_K - 296.15) / (296.15 * R * T_K))
        A_Kn = math.exp(constants["EaKn"] * (T_K - 293.15) / (293.15 * R * T_K))
        A_Kg = math.exp(constants["EaKg"] * (T_K - 293.15) / (293.15 * R * T_K))
        A_Kf = math.exp(constants["EaKf"] * (T_K - 293.15) / (293.15 * R * T_K))
        A_Kig = math.exp(constants["EaKig"] * (T_K - 293.15) / (293.15 * R * T_K))
        A_Kie = math.exp(constants["EaKie"] * (T_K - 293.15) / (293.15 * R * T_K))
        Kn = (theta["mu0"] / theta["sN"]) * A_Kn
        Kg = (theta["betaG0"] / theta["sG"]) * A_Kg
        Kf = (theta["betaF0"] / theta["sF"]) * A_Kf
        iG_T = theta["iG"] / (A_Kig + eps)
        iE_T = theta["iE"] / (A_Kie + eps)
        N_lim = N / (N + Kn + eps)
        G_lim = G / (G + Kg + eps)
        F_lim = F / (F + Kf + eps)
        E_inhibition = 1.0 / (1.0 + iE_T * E)
        G_inhibition_for_fructose = 1.0 / (1.0 + iG_T * G)
        growth_factor = A_mu * N_lim
        glucose_ferm_factor = A_beta * G_lim * E_inhibition
        fructose_ferm_factor = A_beta * F_lim * G_inhibition_for_fructose * E_inhibition
        mu = theta["mu0"] * growth_factor
        beta_G = theta["betaG0"] * glucose_ferm_factor
        beta_F = theta["betaF0"] * fructose_ferm_factor
        maintenance = theta["m0"] * math.exp(constants["Eam"] * (T_K - 293.3) / (293.3 * R * T_K))
        Td = -0.0001 * E**3 + 0.0049 * E**2 - 0.1279 * E + 315.89
        Kd = theta["Kd0"] * math.exp((constants["Cde"] * E) + (constants["Etd"] * (T_K - 305.65)) / (305.65 * R * T_K)) if T_K >= Td else 0.0
        sugar_total = G + F + eps
        maint_avail = sugar_total / (sugar_total + ns["MAINTENANCE_SUGAR_CUTOFF_KG_M3"])
        rates = {channel: _continuous_pulse_rate(t, candidate.pulses.get(channel, tuple())) for channel in CHANNELS}
        dX = (mu - Kd) * X + rates["X"]
        dXd = Kd * X
        dN = -theta["qN"] * growth_factor * X + rates["N"]
        dG = -(
            theta["qXG"] * growth_factor
            + theta["qEG"] * glucose_ferm_factor
            + maintenance * maint_avail * (G / sugar_total)
        ) * X + rates["G"]
        dF = -(
            theta["qXF"] * growth_factor
            + theta["qEF"] * fructose_ferm_factor
            + maintenance * maint_avail * (F / sugar_total)
        ) * X + rates["F"]
        dE = (beta_G + beta_F) * X + rates["E"]
        dCO2 = CO2_PER_ETHANOL_MASS * (beta_G + beta_F) * X
        return [dX, dXd, dN, dG, dF, dE, dCO2]

    try:
        sol = solve_ivp(
            rhs,
            (float(time[0]), float(time[-1])),
            y0,
            t_eval=time,
            method="LSODA",
            rtol=1e-6,
            atol=1e-8,
        )
    except Exception:
        return None
    if not sol.success or sol.y.shape[1] != len(time):
        return None
    values = np.asarray(sol.y.T, dtype=float)
    columns = ["X", "Xd", "N", "G", "F", "E", "CO2"]
    floors = {"X": 1e-8, "Xd": 0.0, "N": -1e-8, "G": -1e-8, "F": -1e-8, "E": 0.0, "CO2": 0.0}
    uppers = {"X": 20.0, "Xd": 30.0, "N": 5.0, "G": 300.0, "F": 300.0, "E": 200.0, "CO2": 350.0}
    for idx, col in enumerate(columns):
        values[:, idx] = np.clip(values[:, idx], floors[col], uppers[col])
    return pd.DataFrame(values, index=time, columns=columns)


def make_synthetic_batch(ns: dict, candidate: ExtendedDesignCandidate, liquid_interval_h: float, co2_interval_h: float, theta: dict[str, float] | None = None):
    time = _sample_grid(candidate.horizon_h, liquid_interval_h, co2_interval_h)
    temperature_c = _temperature_profile(candidate, time)
    nutrient_pulse = np.zeros_like(time, dtype=float)
    nutrient_pulse_rate = np.zeros_like(time, dtype=float)

    initials = {
        "X": float(candidate.initials.get("X", 0.6)),
        "N": float(candidate.initials.get("N", 0.18)),
        "G": float(candidate.initials.get("G", 80.0)),
        "F": float(candidate.initials.get("F", 80.0)),
        "E": float(candidate.initials.get("E", 0.0)),
    }
    raw = pd.DataFrame(
        {
            "t": time,
            "temperatura": temperature_c,
            "pulso_nut": np.zeros_like(time),
            "Viability": initials["X"] / float(ns["VIABILITY_TO_KG_M3"]),
            "YAN": initials["N"] * 1000.0,
            "GLUCOSE": initials["G"],
            "FRUCTOSE": initials["F"],
            "ETANOL": initials["E"],
            "ID": candidate.name,
        }
    )
    dynamic_guess = _extended_dynamic_guess(ns, candidate, theta, time) if theta is not None else None
    measurements = pd.DataFrame(index=time)
    measurements.index.name = "t"
    if dynamic_guess is not None:
        for state in ("X", "N", "G", "F", "E"):
            measurements[state] = dynamic_guess[state]
        initial_guess = dynamic_guess[["X", "N", "G", "F", "E"]].copy()
        initial_guess["Xd"] = dynamic_guess["Xd"]
        initial_guess["CO2"] = dynamic_guess["CO2"]
    else:
        progress = np.clip(time / max(float(candidate.horizon_h), 1.0), 0.0, 1.0)
        measurements["X"] = np.maximum(0.05, initials["X"] * (1.0 + 2.0 * np.exp(-((progress - 0.35) / 0.28) ** 2)))
        measurements["N"] = np.maximum(0.0, initials["N"] * (1.0 - progress * 1.6))
        measurements["G"] = np.maximum(0.0, initials["G"] * (1.0 - progress * 1.2))
        measurements["F"] = np.maximum(0.0, initials["F"] * (1.0 - progress * 0.95))
        measurements["E"] = np.minimum(160.0, initials["E"] + 100.0 * progress)
        initial_guess = measurements.interpolate(limit_direction="both")
    return ns["FermentationBatch"](
        batch_id=f"synthetic_{candidate.name}",
        run_label=candidate.name,
        raw=raw,
        time=time,
        temperature_c=temperature_c,
        nutrient_pulse_kg_m3=nutrient_pulse,
        nutrient_pulse_rate_kg_m3_h=nutrient_pulse_rate,
        measurements=measurements,
        initial_guess=initial_guess,
        initials=initials,
    )


def _pad_schedule(schedule: tuple[tuple[float, float], ...], slots: int, horizon_h: float) -> list[tuple[float, float]]:
    rows = [(float(t), float(a)) for t, a in schedule if float(a) > 0.0]
    rows = sorted(rows, key=lambda item: item[0])[:slots]
    pad_time = float(rows[-1][0]) if rows else 0.0
    pad_time = max(pad_time, float(horizon_h))
    while len(rows) < slots:
        rows.append((pad_time, 0.0))
    return rows


def _add_fixed_pulse_channel(m, channel: str, schedule: tuple[tuple[float, float], ...], horizon_h: float, pulse_slots: int, width_h: float):
    slots = pyo.RangeSet(0, pulse_slots - 1)
    setattr(m, f"{channel}_pulse_slots", slots)
    amount_bounds = CHANNEL_AMOUNT_BOUNDS[channel]
    rows = _pad_schedule(schedule, pulse_slots, horizon_h)
    time_init = {idx: min(max(float(t), 0.0), float(horizon_h)) for idx, (t, _a) in enumerate(rows)}
    amount_init = {idx: min(max(float(a), amount_bounds[0]), amount_bounds[1]) for idx, (_t, a) in enumerate(rows)}

    time_var = pyo.Var(slots, bounds=(0.0, float(horizon_h)), initialize=lambda _m, p: time_init[int(p)])
    amount_var = pyo.Var(slots, bounds=amount_bounds, initialize=lambda _m, p: amount_init[int(p)])
    setattr(m, f"{channel}_pulse_time", time_var)
    setattr(m, f"{channel}_pulse_amount", amount_var)
    for p in slots:
        time_var[p].fix(time_init[int(p)])
        amount_var[p].fix(amount_init[int(p)])

    def pulse_order_rule(_m, p):
        if int(p) == pulse_slots - 1:
            return pyo.Constraint.Skip
        return time_var[p] <= time_var[int(p) + 1]

    setattr(m, f"{channel}_pulse_order", pyo.Constraint(slots, rule=pulse_order_rule))

    def pulse_total_rule(_m):
        return sum(amount_var[p] for p in slots) <= CHANNEL_TOTAL_LIMITS[channel]

    setattr(m, f"{channel}_pulse_total_limit", pyo.Constraint(rule=pulse_total_rule))

    t0 = float(m.t.first())
    eps = 1e-8
    shape = pyo.Expression(slots, m.t, rule=lambda _m, p, t: pyo.exp(-((float(t) - time_var[p]) / width_h) ** 2))
    norm = pyo.Expression(slots, rule=lambda _m, p: sum(m.dt[t] * shape[p, t] for t in m.t if float(t) != t0) + eps)
    rate = pyo.Expression(
        m.t,
        rule=lambda _m, t: 0.0
        if float(t) == t0
        else sum(amount_var[p] * shape[p, t] / norm[p] for p in slots),
    )
    setattr(m, f"{channel}_pulse_shape", shape)
    setattr(m, f"{channel}_pulse_norm", norm)
    setattr(m, f"{channel}_input_rate", rate)


def _add_temperature_design(m, candidate: ExtendedDesignCandidate):
    for t in m.t:
        m.TempC[t].unfix()
    segments = len(candidate.temperature_c)
    m.temperature_segments = pyo.RangeSet(0, segments - 1)
    edges = np.linspace(0.0, float(candidate.horizon_h), segments + 1)
    segment_for_time = {}
    for t in m.t:
        idx = int(np.searchsorted(edges[1:-1], float(t), side="right"))
        segment_for_time[float(t)] = min(idx, segments - 1)
    m.T_set = pyo.Var(
        m.temperature_segments,
        bounds=(15.0, 25.0),
        initialize=lambda _m, s: float(candidate.temperature_c[int(s)]),
    )
    for s in m.temperature_segments:
        m.T_set[s].fix(float(candidate.temperature_c[int(s)]))

    @m.Constraint(m.t)
    def extended_temperature_profile(_m, t):
        return m.TempC[t] == m.T_set[segment_for_time[float(t)]]


def _add_initial_design_vars(m, candidate: ExtendedDesignCandidate):
    t0 = float(m.t.first())
    initial_specs = {
        "X": ("X0", (0.05, 5.0), "X"),
        "N": ("N0", (0.0, 0.35), "N"),
        "G": ("G0", (0.0, 260.0), "G"),
        "F": ("F0", (0.0, 260.0), "F"),
        "E": ("E0", (0.0, 70.0), "E"),
    }
    for state, (var_name, bounds, key) in initial_specs.items():
        getattr(m, state)[t0].unfix()
        value = float(candidate.initials.get(key, pyo.value(getattr(m, state)[t0])))
        var = pyo.Var(bounds=bounds, initialize=value)
        setattr(m, var_name, var)
        var.fix(value)

    m.Xd0 = pyo.Var(bounds=EXTENDED_STATE_BOUNDS["Xd"], initialize=float(candidate.initials.get("Xd", 0.0)))
    m.Xd0.fix(float(candidate.initials.get("Xd", 0.0)))

    @m.Constraint()
    def initial_X(_m):
        return m.X[t0] == m.X0

    @m.Constraint()
    def initial_N(_m):
        return m.N[t0] == m.N0

    @m.Constraint()
    def initial_G(_m):
        return m.G[t0] == m.G0

    @m.Constraint()
    def initial_F(_m):
        return m.F[t0] == m.F0

    @m.Constraint()
    def initial_E(_m):
        return m.E[t0] == m.E0

    @m.Constraint()
    def initial_Xd(_m):
        return m.Xd[t0] == m.Xd0


def build_extended_zenteno_model(
    ns: dict,
    batch,
    theta: dict[str, float],
    parameters: tuple[str, ...],
    candidate: ExtendedDesignCandidate | None,
    pulse_slots: int = 3,
    pulse_width_h: float = 3.0,
    future_outputs: bool = False,
):
    m = ns["build_zenteno_pyomo_model"](
        batch,
        theta_initial=theta,
        fix_parameters=True,
        label_model=False,
        input_mode="measured",
        fix_design_inputs=True,
    )
    model_time = [float(t) for t in m.t]
    if "Xd" in batch.initial_guess.columns:
        xd_values = np.interp(model_time, batch.initial_guess.index.to_numpy(dtype=float), batch.initial_guess["Xd"].to_numpy(dtype=float))
    else:
        xd_values = np.zeros(len(model_time), dtype=float)
    if "CO2" in batch.initial_guess.columns:
        co2_values = np.interp(model_time, batch.initial_guess.index.to_numpy(dtype=float), batch.initial_guess["CO2"].to_numpy(dtype=float))
    else:
        co2_values = np.zeros(len(model_time), dtype=float)
    xd_init = {float(t): float(v) for t, v in zip(model_time, xd_values)}
    co2_init = {float(t): float(v) for t, v in zip(model_time, co2_values)}
    m.Xd = pyo.Var(m.t, bounds=EXTENDED_STATE_BOUNDS["Xd"], initialize=xd_init)
    m.CO2 = pyo.Var(m.t, bounds=EXTENDED_STATE_BOUNDS["CO2"], initialize=co2_init)
    m.dXd = dae.DerivativeVar(m.Xd, wrt=m.t)
    m.dCO2 = dae.DerivativeVar(m.CO2, wrt=m.t)

    t0 = float(m.t.first())
    if candidate is None:
        m.Xd[t0].fix(0.0)
    m.CO2[t0].fix(0.0)

    if candidate is None:
        zero_rate = pyo.Expression(m.t, rule=lambda _m, _t: 0.0)
        m.N_input_rate = pyo.Expression(m.t, rule=lambda _m, t: m.N_pulse_rate[t])
        m.G_input_rate = zero_rate
        m.F_input_rate = pyo.Expression(m.t, rule=lambda _m, _t: 0.0)
        m.E_input_rate = pyo.Expression(m.t, rule=lambda _m, _t: 0.0)
        m.X_input_rate = pyo.Expression(m.t, rule=lambda _m, _t: 0.0)
    else:
        _add_temperature_design(m, candidate)
        _add_initial_design_vars(m, candidate)
        for channel in CHANNELS:
            _add_fixed_pulse_channel(
                m,
                channel,
                candidate.pulses.get(channel, tuple()),
                candidate.horizon_h,
                pulse_slots,
                pulse_width_h,
            )

    for cname in ("X_balance", "N_balance", "G_balance", "F_balance", "E_balance"):
        if hasattr(m, cname):
            getattr(m, cname).deactivate()

    eps = 1e-8
    m.CO2_rate = pyo.Expression(m.t, rule=lambda _m, t: CO2_PER_ETHANOL_MASS * (m.beta_G[t] + m.beta_F[t]) * m.X[t])

    @m.Constraint(m.t)
    def extended_X_balance(_m, t):
        return m.dX[t] == m.mu[t] * m.X[t] - m.Kd[t] * m.X[t] + m.X_input_rate[t]

    @m.Constraint(m.t)
    def extended_Xd_balance(_m, t):
        return m.dXd[t] == m.Kd[t] * m.X[t]

    @m.Constraint(m.t)
    def extended_N_balance(_m, t):
        return m.dN[t] == -m.qN_growth[t] * m.X[t] + m.N_input_rate[t]

    @m.Constraint(m.t)
    def extended_G_balance(_m, t):
        return m.dG[t] == -(
            m.qXG_growth[t]
            + m.qEG_ferm[t]
            + m.maintenance[t] * m.maintenance_sugar_availability[t] * (m.G_eff[t] / (m.G_eff[t] + m.F_eff[t] + eps))
        ) * m.X[t] + m.G_input_rate[t]

    @m.Constraint(m.t)
    def extended_F_balance(_m, t):
        return m.dF[t] == -(
            m.qXF_growth[t]
            + m.qEF_ferm[t]
            + m.maintenance[t] * m.maintenance_sugar_availability[t] * (m.F_eff[t] / (m.G_eff[t] + m.F_eff[t] + eps))
        ) * m.X[t] + m.F_input_rate[t]

    @m.Constraint(m.t)
    def extended_E_balance(_m, t):
        return m.dE[t] == (m.beta_G[t] + m.beta_F[t]) * m.X[t] + m.E_input_rate[t]

    @m.Constraint(m.t)
    def extended_CO2_balance(_m, t):
        return m.dCO2[t] == m.CO2_rate[t]

    return m


def discretize_model(m):
    if getattr(m, "_dae_discretized", False):
        return m
    nfe = len(list(m.t)) - 1
    pyo.TransformationFactory("dae.finite_difference").apply_to(m, scheme="BACKWARD", nfe=nfe, wrt=m.t)
    m._dae_discretized = True
    return m


def _time_in_set(m, time_h: float) -> float | None:
    key = round(float(time_h), 8)
    for t in m.t:
        if abs(float(t) - key) <= 1e-7:
            return float(t)
    return None


def label_extended_model(
    ns: dict,
    m,
    batch,
    parameters: tuple[str, ...],
    candidate: ExtendedDesignCandidate | None,
    future_outputs: bool,
    liquid_interval_h: float = 12.0,
    co2_interval_h: float = 6.0,
):
    m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    for name in parameters:
        var = getattr(m, name)
        m.unknown_parameters[var] = pyo.value(var)

    m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    if future_outputs and candidate is not None:
        for s in m.temperature_segments:
            m.experiment_inputs[m.T_set[s]] = None
        for name in ("X0", "N0", "G0", "F0", "E0", "Xd0"):
            m.experiment_inputs[getattr(m, name)] = None
        for channel in CHANNELS:
            for p in getattr(m, f"{channel}_pulse_slots"):
                m.experiment_inputs[getattr(m, f"{channel}_pulse_time")[p]] = None
                m.experiment_inputs[getattr(m, f"{channel}_pulse_amount")[p]] = None
    else:
        for t in m.t:
            if hasattr(m, "TempC"):
                m.experiment_inputs[m.TempC[t]] = None

    m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
    m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)

    if not future_outputs:
        objective_scale = ns["objective_measurement_error_for_batches"]([batch.batch_id])
        state_vars = {"X": m.X, "N": m.N, "G": m.G, "F": m.F, "E": m.E}
        for state, var in state_vars.items():
            for t, measurement in batch.measurements[state].dropna().items():
                mt = _time_in_set(m, float(t))
                if mt is None:
                    continue
                m.experiment_outputs[var[mt]] = float(measurement)
                m.measurement_error[var[mt]] = float(objective_scale[state])
        return m

    liquid_times = np.arange(0.0, float(candidate.horizon_h) + 1e-9, float(liquid_interval_h))
    co2_times = np.arange(0.0, float(candidate.horizon_h) + 1e-9, float(co2_interval_h))
    state_vars = {"X": m.X, "Xd": m.Xd, "N": m.N, "G": m.G, "F": m.F, "E": m.E}
    for state, var in state_vars.items():
        for t in liquid_times:
            mt = _time_in_set(m, float(t))
            if mt is None:
                continue
            m.experiment_outputs[var[mt]] = 0.0
            m.measurement_error[var[mt]] = FUTURE_MEASUREMENT_ERROR[state]
    for t in co2_times:
        mt = _time_in_set(m, float(t))
        if mt is None:
            continue
        m.experiment_outputs[m.CO2[mt]] = 0.0
        m.measurement_error[m.CO2[mt]] = FUTURE_MEASUREMENT_ERROR["CO2"]
    return m


def fim_for_experiment(experiment, parameters: tuple[str, ...], step: float, solver) -> np.ndarray:
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
    out = {
        f"{prefix}logdet": float(np.sum(np.log(eig_pos))),
        f"{prefix}min_eigenvalue": float(np.min(eig)) if eig.size else np.nan,
        f"{prefix}min_relative_eigenvalue": float(np.min(eig) / max_eig) if max_eig > 0 else np.nan,
        f"{prefix}condition_number": float(eig_pos.max() / eig_pos.min()) if eig_pos.size else np.nan,
        f"{prefix}trace": float(np.trace(fim)),
        f"{prefix}trace_inv": float(np.trace(stable_inverse(fim))),
    }
    cov = stable_inverse(fim)
    for name in parameters:
        idx = parameters.index(name)
        out[f"{prefix}var_{name}"] = float(cov[idx, idx])
    return out


def variance_reduction(prior: np.ndarray, combined: np.ndarray, parameters: tuple[str, ...]) -> dict[str, float]:
    prior_cov = stable_inverse(prior)
    combined_cov = stable_inverse(combined)
    rows = {}
    ratios = []
    for name in parameters:
        idx = parameters.index(name)
        before = float(prior_cov[idx, idx])
        after = float(combined_cov[idx, idx])
        ratio = after / before if before > 0 else np.nan
        rows[f"var_ratio_{name}"] = ratio
        rows[f"var_reduction_{name}"] = 1.0 - ratio if np.isfinite(ratio) else np.nan
        if name in CORRECTABLE:
            ratios.append(ratio)
    if ratios:
        rows["correctable_mean_var_reduction"] = float(1.0 - np.nanmean(ratios))
        rows["correctable_worst_var_reduction"] = float(1.0 - np.nanmax(ratios))
    return rows


def current_prior_fim(ns: dict, theta: dict[str, float], parameters: tuple[str, ...], batches: tuple[str, ...], step: float, solver) -> np.ndarray:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    param_key = "full15" if tuple(parameters) == FULL15 else f"p{len(parameters)}"
    cache = RESULTS_DIR / f"historical_prior_fim_{param_key}_step{step:g}.csv"
    if cache.exists():
        return pd.read_csv(cache, index_col=0).loc[list(parameters), list(parameters)].to_numpy(dtype=float)
    fim_total = np.zeros((len(parameters), len(parameters)), dtype=float)
    rows = []
    for batch_id in batches:
        print(f"[prior] batch={batch_id} parameters={len(parameters)}", flush=True)
        experiment = HistoricalExtendedExperiment(ns, batch_id, theta, parameters)
        fim = fim_for_experiment(experiment, parameters, step, solver)
        fim_total += fim
        rows.append({"batch": batch_id, **fim_metrics(fim, parameters, prefix="batch_")})
    pd.DataFrame(fim_total, index=parameters, columns=parameters).to_csv(cache)
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"historical_prior_fim_{param_key}_batch_metrics.csv", index=False)
    return fim_total


def candidate_library() -> list[ExtendedDesignCandidate]:
    def pulses(**kwargs):
        return {channel: tuple(kwargs.get(channel, tuple())) for channel in CHANNELS}

    base = []
    base.append(
        ExtendedDesignCandidate(
            "low_temp_growth_separation_plus_co2",
            "legacy_plus",
            168.0,
            {"X": 0.55, "N": 0.24, "G": 75.0, "F": 70.0, "E": 0.0, "Xd": 0.0},
            (15.0, 17.0, 20.0, 22.0),
            pulses(N=((0.0, 0.07), (32.0, 0.07), (64.0, 0.06))),
            "Previous robust N/T design, now measured with Xd and online CO2.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "glucose_rich_growth_yield",
            "sugar_initial",
            168.0,
            {"X": 0.45, "N": 0.26, "G": 150.0, "F": 25.0, "E": 0.0, "Xd": 0.0},
            (18.0, 22.0, 24.0, 20.0),
            pulses(N=((0.0, 0.08), (36.0, 0.06)), F=((72.0, 35.0),)),
            "High glucose must isolates glucose growth/fermentation terms and tests fructose response after glucose history.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "fructose_rich_iG_probe",
            "sugar_initial",
            192.0,
            {"X": 0.45, "N": 0.26, "G": 25.0, "F": 150.0, "E": 0.0, "Xd": 0.0},
            (18.0, 20.0, 24.0, 22.0),
            pulses(N=((0.0, 0.08), (48.0, 0.06)), G=((72.0, 35.0),)),
            "Fructose-rich must plus glucose pulse separates fructose capacity from glucose inhibition iG.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "low_sugar_saturation_scan",
            "saturation",
            144.0,
            {"X": 0.35, "N": 0.16, "G": 22.0, "F": 22.0, "E": 0.0, "Xd": 0.0},
            (16.0, 20.0, 24.0, 20.0),
            pulses(G=((36.0, 18.0),), F=((72.0, 18.0),), N=((0.0, 0.04), (48.0, 0.04))),
            "Low/intermediate G/F levels excite sG and sF rather than keeping saturation terms flat.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "yan_saturation_scan",
            "saturation",
            144.0,
            {"X": 0.35, "N": 0.035, "G": 90.0, "F": 90.0, "E": 0.0, "Xd": 0.0},
            (15.0, 18.0, 22.0, 22.0),
            pulses(N=((24.0, 0.025), (48.0, 0.045), (72.0, 0.065))),
            "Designed N ladder through low/intermediate YAN targets sN/qN separation.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "ethanol_initial_challenge",
            "ethanol",
            192.0,
            {"X": 0.60, "N": 0.22, "G": 90.0, "F": 90.0, "E": 35.0, "Xd": 0.0},
            (18.0, 22.0, 25.0, 22.0),
            pulses(N=((0.0, 0.07), (48.0, 0.05))),
            "Initial ethanol decouples ethanol inhibition iE from ethanol produced by fermentation.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "ethanol_pulse_death_probe",
            "death",
            216.0,
            {"X": 0.75, "N": 0.18, "G": 80.0, "F": 80.0, "E": 5.0, "Xd": 0.0},
            (20.0, 25.0, 25.0, 22.0),
            pulses(N=((0.0, 0.05),), E=((96.0, 25.0),)),
            "Late ethanol stress with Xd observation targets iE/Kd0 separately from growth.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "high_biomass_low_N_maintenance",
            "maintenance",
            168.0,
            {"X": 2.2, "N": 0.025, "G": 70.0, "F": 70.0, "E": 15.0, "Xd": 0.0},
            (16.0, 18.0, 20.0, 18.0),
            pulses(G=((48.0, 20.0),), F=((96.0, 20.0),), N=((72.0, 0.015),)),
            "High biomass and low N create low-growth sugar consumption windows for m0.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "viable_biomass_step",
            "biomass_input",
            168.0,
            {"X": 0.30, "N": 0.22, "G": 90.0, "F": 90.0, "E": 0.0, "Xd": 0.0},
            (18.0, 22.0, 22.0, 18.0),
            pulses(N=((0.0, 0.06), (48.0, 0.05)), X=((48.0, 1.2),)),
            "Known viable biomass addition tests whether rates scale with X and improves q/yield separation.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "glucose_pulse_after_N_depletion",
            "pulse_separation",
            192.0,
            {"X": 0.55, "N": 0.055, "G": 45.0, "F": 85.0, "E": 0.0, "Xd": 0.0},
            (17.0, 20.0, 24.0, 22.0),
            pulses(G=((72.0, 45.0),), N=((0.0, 0.035),)),
            "Glucose pulse after likely N limitation separates fermentation/maintenance from growth uptake.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "fructose_pulse_after_N_depletion",
            "pulse_separation",
            192.0,
            {"X": 0.55, "N": 0.055, "G": 85.0, "F": 45.0, "E": 0.0, "Xd": 0.0},
            (17.0, 20.0, 24.0, 22.0),
            pulses(F=((72.0, 45.0),), N=((0.0, 0.035),)),
            "Fructose pulse after N limitation targets qXF/qEF/sF and fructose-specific residuals.",
        )
    )
    base.append(
        ExtendedDesignCandidate(
            "combined_stress_long_horizon",
            "integrated",
            240.0,
            {"X": 0.70, "N": 0.20, "G": 110.0, "F": 110.0, "E": 15.0, "Xd": 0.0},
            (15.0, 22.0, 25.0, 18.0),
            pulses(N=((0.0, 0.05), (60.0, 0.05)), E=((120.0, 20.0),), G=((72.0, 25.0),), F=((96.0, 25.0),)),
            "Broad stress/input design for integrated all-parameter information.",
        )
    )
    return base


def evaluate_candidate(ns, theta, candidate, parameters, prior_fim, step, solver, pulse_slots, liquid_interval_h, co2_interval_h):
    try:
        experiment = FutureExtendedExperiment(ns, candidate, theta, parameters, pulse_slots, liquid_interval_h, co2_interval_h)
        fim = fim_for_experiment(experiment, parameters, step, solver)
        combined = prior_fim + fim
        row = {
            "candidate": candidate.name,
            "family": candidate.family,
            "horizon_h": candidate.horizon_h,
            "temperature_c": ", ".join(f"{v:g}" for v in candidate.temperature_c),
            "status": "ok",
            "error": "",
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


def greedy_select(fims: dict[str, np.ndarray], prior_fim: np.ndarray, candidates: dict[str, ExtendedDesignCandidate], parameters: tuple[str, ...], campaign_size: int):
    selected = []
    remaining = set(fims)
    current = prior_fim.copy()
    rows = []
    for step_idx in range(1, int(campaign_size) + 1):
        best_name = None
        best_score = -np.inf
        best_metrics = None
        for name in sorted(remaining):
            trial = current + fims[name]
            metrics = fim_metrics(trial, parameters, prefix="campaign_")
            reductions = variance_reduction(prior_fim, trial, parameters)
            score = (
                metrics["campaign_logdet"]
                + 5.0 * reductions.get("correctable_mean_var_reduction", 0.0)
                + math.log(max(metrics["campaign_min_relative_eigenvalue"], 1e-16))
            )
            if score > best_score:
                best_name = name
                best_score = score
                best_metrics = {**metrics, **reductions}
        if best_name is None:
            break
        current = current + fims[best_name]
        remaining.remove(best_name)
        selected.append(best_name)
        c = candidates[best_name]
        rows.append(
            {
                "campaign_order": step_idx,
                "candidate": best_name,
                "family": c.family,
                "horizon_h": c.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                "rationale": c.rationale,
                **best_metrics,
            }
        )
    return pd.DataFrame(rows), current


def campaign_rows_for_names(names: list[str], fims: dict[str, np.ndarray], prior_fim: np.ndarray, candidates: dict[str, ExtendedDesignCandidate], parameters: tuple[str, ...]):
    current = prior_fim.copy()
    rows = []
    for step_idx, name in enumerate(names, start=1):
        current = current + fims[name]
        c = candidates[name]
        metrics = fim_metrics(current, parameters, prefix="campaign_")
        reductions = variance_reduction(prior_fim, current, parameters)
        rows.append(
            {
                "campaign_order": step_idx,
                "candidate": name,
                "family": c.family,
                "horizon_h": c.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                "rationale": c.rationale,
                **metrics,
                **reductions,
            }
        )
    return pd.DataFrame(rows), current


def _has_positive_pulse(candidate: ExtendedDesignCandidate, channel: str) -> bool:
    return any(float(amount) > 0.0 for _time, amount in candidate.pulses.get(channel, tuple()))


def enforce_x_pulse_coverage(selected: pd.DataFrame, fims: dict[str, np.ndarray], prior_fim: np.ndarray, candidates: dict[str, ExtendedDesignCandidate], parameters: tuple[str, ...]):
    if selected.empty:
        return selected, prior_fim, False
    selected_names = list(selected["candidate"])
    if any(_has_positive_pulse(candidates[name], "X") for name in selected_names):
        selected_rows, selected_fim = campaign_rows_for_names(selected_names, fims, prior_fim, candidates, parameters)
        return selected_rows, selected_fim, False
    x_candidates = [name for name in fims if _has_positive_pulse(candidates[name], "X") and name not in selected_names]
    if not x_candidates:
        selected_rows, selected_fim = campaign_rows_for_names(selected_names, fims, prior_fim, candidates, parameters)
        return selected_rows, selected_fim, False

    best_names = selected_names
    best_score = -np.inf
    for replacement in x_candidates:
        for drop_idx, dropped in enumerate(selected_names):
            trial_names = selected_names.copy()
            trial_names[drop_idx] = replacement
            if len(set(trial_names)) != len(trial_names):
                continue
            trial_selected, trial_fim = campaign_rows_for_names(trial_names, fims, prior_fim, candidates, parameters)
            final = trial_selected.iloc[-1]
            score = (
                float(final["campaign_logdet"])
                + 5.0 * float(final.get("correctable_mean_var_reduction", 0.0))
                + math.log(max(float(final["campaign_min_relative_eigenvalue"]), 1e-16))
            )
            if score > best_score:
                best_score = score
                best_names = trial_names
    selected_rows, selected_fim = campaign_rows_for_names(best_names, fims, prior_fim, candidates, parameters)
    selected_rows["coverage_note"] = ""
    selected_rows.loc[selected_rows["candidate"].map(lambda name: _has_positive_pulse(candidates[name], "X")), "coverage_note"] = "forced X-pulse coverage"
    return selected_rows, selected_fim, True


def simulate_candidate(ns, theta, candidate, pulse_slots: int, liquid_interval_h: float, co2_interval_h: float):
    batch = make_synthetic_batch(ns, candidate, liquid_interval_h, co2_interval_h, theta)
    m = build_extended_zenteno_model(ns, batch, theta, parameters=ACCEPTED7, candidate=candidate, pulse_slots=pulse_slots, future_outputs=True)
    discretize_model(m)
    solver = make_solver(ns)
    result = solver.solve(m, tee=False)
    rows = []
    for t in m.t:
        rows.append(
            {
                "candidate": candidate.name,
                "family": candidate.family,
                "t": float(t),
                "temperature_c": pyo.value(m.TempC[t]),
                "X": pyo.value(m.X[t]),
                "Xd": pyo.value(m.Xd[t]),
                "N": pyo.value(m.N[t]),
                "G": pyo.value(m.G[t]),
                "F": pyo.value(m.F[t]),
                "E": pyo.value(m.E[t]),
                "CO2": pyo.value(m.CO2[t]),
                "CO2_rate": pyo.value(m.CO2_rate[t]),
                "N_input_rate": pyo.value(m.N_input_rate[t]),
                "G_input_rate": pyo.value(m.G_input_rate[t]),
                "F_input_rate": pyo.value(m.F_input_rate[t]),
                "E_input_rate": pyo.value(m.E_input_rate[t]),
                "X_input_rate": pyo.value(m.X_input_rate[t]),
                "termination": str(result.solver.termination_condition),
            }
        )
    return pd.DataFrame(rows)


def _pulse_table(candidates: list[ExtendedDesignCandidate]) -> pd.DataFrame:
    rows = []
    for c in candidates:
        for channel in CHANNELS:
            for time_h, amount in c.pulses.get(channel, tuple()):
                if amount <= 0.0:
                    continue
                rows.append(
                    {
                        "candidate": c.name,
                        "family": c.family,
                        "channel": channel,
                        "time_h": float(time_h),
                        "amount": float(amount),
                    }
                )
    return pd.DataFrame(rows)


def write_input_plots(candidates: list[ExtendedDesignCandidate], selected: pd.DataFrame):
    for stale in RESULTS_DIR.glob("extended_campaign_inputs_*.png"):
        stale.unlink()
    selected_names = list(selected["candidate"])
    candidates_by_name = {c.name: c for c in candidates}
    rows = [candidates_by_name[name] for name in selected_names if name in candidates_by_name]
    pulse_df = _pulse_table(rows)
    for c in rows:
        fig, axes = plt.subplots(7, 1, figsize=(10, 12), sharex=False)
        time = np.linspace(0.0, c.horizon_h, 400)
        axes[0].step(time, _temperature_profile(c, time), where="post", color="tab:red")
        axes[0].set_ylabel("T C")
        axes[0].set_xlim(0.0, c.horizon_h)
        init_labels = ["X", "N", "G", "F", "E"]
        init_values = [c.initials.get(name, 0.0) for name in init_labels]
        axes[1].bar(init_labels, init_values, color=["tab:purple", "tab:brown", "tab:blue", "tab:orange", "tab:green"])
        axes[1].set_ylabel("initial")
        axes[1].tick_params(axis="x", labelrotation=0)
        for ax, channel in zip(axes[2:], CHANNELS):
            rows_ch = pulse_df[(pulse_df["candidate"].eq(c.name)) & (pulse_df["channel"].eq(channel))]
            if rows_ch.empty:
                ax.axhline(0.0, color="0.8", lw=1)
            else:
                ax.vlines(rows_ch["time_h"], 0.0, rows_ch["amount"], lw=3)
                ax.scatter(rows_ch["time_h"], rows_ch["amount"], s=35)
            ax.set_ylabel(channel)
            ax.set_xlim(0.0, c.horizon_h)
            ax.grid(True, alpha=0.25)
        axes[-1].set_xlabel("t [h]")
        fig.suptitle("\n".join(textwrap.wrap(f"{c.name}: {c.rationale}", width=95)), fontsize=10, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(RESULTS_DIR / f"extended_campaign_inputs_{c.name}.png", dpi=170)
        plt.close(fig)


def write_report(selected: pd.DataFrame, ranking: pd.DataFrame, parameter_reduction: pd.DataFrame):
    selected_columns = [
        "campaign_order",
        "candidate",
        "family",
        "horizon_h",
        "temperature_c",
        "campaign_logdet",
        "campaign_min_relative_eigenvalue",
        "campaign_condition_number",
        "correctable_mean_var_reduction",
        "correctable_worst_var_reduction",
    ]
    if "coverage_note" in selected.columns:
        selected_columns.append("coverage_note")
    selected_columns.append("rationale")
    lines = [
        "# Extended fermentation campaign DOE",
        "",
        "This report screens synthetic-must experiments using the Pyomo DoE sequential FIM workflow.",
        "",
        "Extended design actions: temperature setpoint profile, initial synthetic-must composition, and fixed candidate pulses of YAN, glucose, fructose, ethanol, and viable biomass.",
        "",
        "Extended model states: viable biomass `X`, dead/non-viable biomass `Xd`, assimilable nitrogen `N`, glucose `G`, fructose `F`, ethanol `E`, and cumulative `CO2`.",
        "",
        "CO2 is represented here as cumulative mass-equivalent production with fixed stoichiometry `dCO2/dt = (44.01/46.07) * dE_prod/dt`. Zenteno (2010) includes a richer CO2 readout with growth-associated production and gas-liquid transfer/saturation; use that layer before final calibration if the online sensor reports gas flow instead of cumulative CO2.",
        "",
        "Manual liquid/biomass sampling is assumed every 12 h. CO2 is assumed online and downsampled every 6 h in the FIM to avoid overwhelming the manual observations.",
        "",
        "## Selected Campaign",
        "",
        selected[selected_columns].to_markdown(index=False),
        "",
        "## Top Single-Experiment Candidates",
        "",
        ranking.head(12)[
            [
                "candidate",
                "family",
                "horizon_h",
                "combined_logdet",
                "combined_min_relative_eigenvalue",
                "combined_condition_number",
                "correctable_mean_var_reduction",
                "correctable_worst_var_reduction",
                "status",
            ]
        ].to_markdown(index=False),
        "",
        "## Final Parameter Variance Reductions",
        "",
        parameter_reduction.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "- `qXG/qXF` are attacked through orthogonal glucose/fructose initial compositions and single-sugar pulses.",
        "- `sG/sF/sN` are attacked through low/intermediate substrate and YAN ladders rather than only rich musts.",
        "- `iE` is attacked by initial/pulsed ethanol so inhibition is not only generated endogenously by fermentation.",
        "- `Kd0` is attacked by adding `Xd` as an observed state; without `Xd`, death remains confounded with slow growth or low viability.",
        "- `m0` is attacked through high-biomass, low-N sugar consumption windows where growth is constrained.",
    ]
    (RESULTS_DIR / "extended_campaign_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Extended fermentation campaign DOE with Pyomo DoE FIM screening.")
    parser.add_argument("--prior-batches", nargs="*", default=["25026", "25086", "25150", "25170"])
    parser.add_argument("--campaign-size", type=int, default=9)
    parser.add_argument("--max-candidates", type=int, default=0)
    parser.add_argument("--step", type=float, default=1e-2)
    parser.add_argument("--pulse-slots", type=int, default=3)
    parser.add_argument("--liquid-interval-h", type=float, default=12.0)
    parser.add_argument("--co2-interval-h", type=float, default=6.0)
    parser.add_argument("--skip-prior", action="store_true", help="Use a small ridge prior instead of recomputing historical FIM.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ns = load_notebook_context()
    theta = complete_final_theta(ns)
    solver = make_solver(ns)
    parameters = FULL15

    candidates = candidate_library()
    if args.max_candidates and args.max_candidates > 0:
        candidates = candidates[: args.max_candidates]
    candidate_by_name = {c.name: c for c in candidates}

    if args.skip_prior:
        prior_fim = np.eye(len(parameters)) * 1e-6
    else:
        prior_fim = current_prior_fim(ns, theta, parameters, tuple(args.prior_batches), args.step, solver)
    pd.DataFrame(prior_fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / "extended_historical_prior_fim.csv")

    rows = []
    fims = {}
    for idx, candidate in enumerate(candidates, start=1):
        print(f"[candidate {idx}/{len(candidates)}] {candidate.name}", flush=True)
        row, fim = evaluate_candidate(
            ns,
            theta,
            candidate,
            parameters,
            prior_fim,
            args.step,
            solver,
            args.pulse_slots,
            args.liquid_interval_h,
            args.co2_interval_h,
        )
        rows.append(row)
        if fim is not None:
            fims[candidate.name] = fim
            pd.DataFrame(fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / f"fim_new_{candidate.name}.csv")
    ranking = pd.DataFrame(rows)
    if not ranking.empty:
        ranking["solve_ok"] = ranking["status"].eq("ok")
        ranking = ranking.sort_values(
            ["solve_ok", "correctable_mean_var_reduction", "combined_logdet"],
            ascending=[False, False, False],
        )
    ranking.to_csv(RESULTS_DIR / "extended_candidate_ranking.csv", index=False)

    if not fims:
        print("No successful candidate FIMs.")
        return 1
    selected, campaign_fim = greedy_select(fims, prior_fim, candidate_by_name, parameters, min(args.campaign_size, len(fims)))
    selected, campaign_fim, forced_x_coverage = enforce_x_pulse_coverage(selected, fims, prior_fim, candidate_by_name, parameters)
    selected.to_csv(RESULTS_DIR / "extended_campaign_selected.csv", index=False)
    pd.DataFrame(campaign_fim, index=parameters, columns=parameters).to_csv(RESULTS_DIR / "extended_campaign_fim.csv")

    reductions = variance_reduction(prior_fim, campaign_fim, parameters)
    parameter_rows = []
    for name in parameters:
        parameter_rows.append(
            {
                "parameter": name,
                "group": "accepted_current" if name in ACCEPTED7 else "currently_fixed_target",
                "campaign_var_ratio": reductions.get(f"var_ratio_{name}", np.nan),
                "campaign_var_reduction": reductions.get(f"var_reduction_{name}", np.nan),
            }
        )
    parameter_reduction = pd.DataFrame(parameter_rows)
    parameter_reduction.to_csv(RESULTS_DIR / "extended_campaign_parameter_reduction.csv", index=False)

    selected_candidates = [candidate_by_name[name] for name in selected["candidate"]]
    _pulse_table(selected_candidates).to_csv(RESULTS_DIR / "extended_campaign_pulses.csv", index=False)
    pd.DataFrame(
        [
            {
                "candidate": c.name,
                "family": c.family,
                "horizon_h": c.horizon_h,
                "temperature_c": ", ".join(f"{v:g}" for v in c.temperature_c),
                **{f"{key}0": value for key, value in c.initials.items()},
                "rationale": c.rationale,
            }
            for c in selected_candidates
        ]
    ).to_csv(RESULTS_DIR / "extended_campaign_protocols.csv", index=False)

    sim_rows = []
    for c in selected_candidates:
        try:
            sim_rows.append(simulate_candidate(ns, theta, c, args.pulse_slots, args.liquid_interval_h, args.co2_interval_h))
        except Exception as err:
            print(f"[simulation failed] {c.name}: {type(err).__name__}: {err}", flush=True)
    if sim_rows:
        pd.concat(sim_rows, ignore_index=True).to_csv(RESULTS_DIR / "extended_campaign_simulations.csv", index=False)

    write_input_plots(candidates, selected)
    write_report(selected, ranking, parameter_reduction)
    metadata = {
        "prior_batches": list(args.prior_batches),
        "campaign_size": int(args.campaign_size),
        "step": float(args.step),
        "pulse_slots": int(args.pulse_slots),
        "liquid_interval_h": float(args.liquid_interval_h),
        "co2_interval_h": float(args.co2_interval_h),
        "parameters": list(parameters),
        "n_candidates": len(candidates),
        "n_successful_candidates": len(fims),
        "forced_x_pulse_coverage": bool(forced_x_coverage),
    }
    (RESULTS_DIR / "extended_campaign_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("\nSelected campaign:")
    print(selected[["campaign_order", "candidate", "family", "horizon_h", "correctable_mean_var_reduction", "campaign_condition_number"]].to_string(index=False))
    print(f"\nResults written to: {RESULTS_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
