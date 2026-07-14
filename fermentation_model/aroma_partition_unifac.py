from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

AROMA_SPECIES_TO_CHEMICAL = {
    "ethyl_acetate": "ethyl acetate",
    "isoamyl_acetate": "isoamyl acetate",
    "ethyl_octanoate": "ethyl octanoate",
}

TRAP_EFFICIENCY = {
    "ethyl_acetate": 0.90,
    "isoamyl_acetate": 0.88,
    "ethyl_octanoate": 0.85,
}

R_PA_M3_MOL_K = 8.31446261815324
TRACE_MOLES_PER_L = 1e-12

SOLVENT_COMPONENTS = {
    "water": "water",
    "ethanol": "ethanol",
    "glucose": "glucose",
}


class UNIFACUnavailable(RuntimeError):
    pass


@lru_cache(maxsize=1)
def _thermo_modules():
    try:
        from chemicals.identifiers import search_chemical
        from thermo import Chemical
        from thermo.unifac import DDBST_UNIFAC_assignments, UFIP, UFSG, UNIFAC_gammas, load_group_assignments_DDBST
    except Exception as err:  # pragma: no cover - depends on local environment
        raise UNIFACUnavailable(
            "The packages 'thermo' and 'chemicals' are required for UNIFAC partition calculations."
        ) from err
    load_group_assignments_DDBST()
    return {
        "search_chemical": search_chemical,
        "Chemical": Chemical,
        "DDBST_UNIFAC_assignments": DDBST_UNIFAC_assignments,
        "UFIP": UFIP,
        "UFSG": UFSG,
        "UNIFAC_gammas": UNIFAC_gammas,
    }


@lru_cache(maxsize=None)
def _metadata(name: str):
    mods = _thermo_modules()
    return mods["search_chemical"](name)


@lru_cache(maxsize=None)
def _groups(name: str) -> dict[int, int]:
    mods = _thermo_modules()
    meta = _metadata(name)
    groups = mods["DDBST_UNIFAC_assignments"].get(meta.InChI_key)
    if groups is None:
        raise UNIFACUnavailable(f"No DDBST original-UNIFAC group assignment found for {name!r}.")
    return groups


@lru_cache(maxsize=None)
def _mw(name: str) -> float:
    return float(_metadata(name).MW)


def molecular_weight(name: str) -> float:
    return _mw(name)


def unifac_groups(name: str) -> dict[int, int]:
    return dict(_groups(name))


def unifac_subgroup_constants(group_ids: tuple[int, ...] | list[int] | set[int]) -> dict[int, dict[str, float | int]]:
    mods = _thermo_modules()
    constants = {}
    for group_id in group_ids:
        subgroup = mods["UFSG"][int(group_id)]
        constants[int(group_id)] = {
            "R": float(subgroup.R),
            "Q": float(subgroup.Q),
            "main_group_id": int(subgroup.main_group_id),
        }
    return constants


def unifac_interaction_parameter(subgroup_from: int, subgroup_to: int) -> float | None:
    mods = _thermo_modules()
    constants = unifac_subgroup_constants((int(subgroup_from), int(subgroup_to)))
    main_from = int(constants[int(subgroup_from)]["main_group_id"])
    main_to = int(constants[int(subgroup_to)]["main_group_id"])
    return mods["UFIP"].get(main_from, {}).get(main_to)


def vapor_pressure_correlation_data() -> dict[str, dict[str, float | str]]:
    """Return the literature vapor-pressure correlation selected by thermo.

    The symbolic Pyomo UNIFAC implementation uses these coefficients directly
    instead of fitting a temperature surrogate.
    """
    mods = _thermo_modules()
    data = {}
    for species, chemical_name in AROMA_SPECIES_TO_CHEMICAL.items():
        vp = mods["Chemical"](chemical_name, T=293.15).VaporPressure
        method = str(vp.method)
        payload: dict[str, float | str] = {"method": method, "chemical": chemical_name}
        if method in getattr(vp, "Antoine_parameters", {}):
            pars = vp.Antoine_parameters[method]
            payload.update(
                {
                    "type": "Antoine",
                    "A": float(pars["A"]),
                    "B": float(pars["B"]),
                    "C": float(pars["C"]),
                    "base": float(pars["base"]),
                    "Tmin": float(pars.get("Tmin", float("nan"))),
                    "Tmax": float(pars.get("Tmax", float("nan"))),
                }
            )
        elif method in getattr(vp, "Wagner_original_parameters", {}):
            pars = vp.Wagner_original_parameters[method]
            payload.update(
                {
                    "type": "Wagner_original",
                    "Tc": float(pars["Tc"]),
                    "Pc": float(pars["Pc"]),
                    "a": float(pars["a"]),
                    "b": float(pars["b"]),
                    "c": float(pars["c"]),
                    "d": float(pars["d"]),
                    "Tmin": float(pars.get("Tmin", float("nan"))),
                    "Tmax": float(pars.get("Tmax", float("nan"))),
                }
            )
        elif method in getattr(vp, "Wagner_parameters", {}):
            pars = vp.Wagner_parameters[method]
            payload.update(
                {
                    "type": "Wagner",
                    "Tc": float(pars["Tc"]),
                    "Pc": float(pars["Pc"]),
                    "a": float(pars["a"]),
                    "b": float(pars["b"]),
                    "c": float(pars["c"]),
                    "d": float(pars["d"]),
                    "Tmin": float(pars.get("Tmin", float("nan"))),
                    "Tmax": float(pars.get("Tmax", float("nan"))),
                }
            )
        else:
            raise UNIFACUnavailable(f"Unsupported vapor-pressure method {method!r} for {chemical_name!r}.")
        data[species] = payload
    return data


def antoine_correlation_data(
    process_min_temp_c: float = 10.0,
    process_max_temp_c: float = 30.0,
) -> dict[str, dict[str, float | str | bool]]:
    """Return Antoine coefficients selected for the fermentation range.

    The selection prefers Antoine correlations whose nominal temperature range
    covers the fermentation window. If no Antoine range covers the window, the
    nearest available Antoine correlation is returned and flagged as
    extrapolated.
    """
    mods = _thermo_modules()
    t_min = float(process_min_temp_c) + 273.15
    t_max = float(process_max_temp_c) + 273.15
    data: dict[str, dict[str, float | str | bool]] = {}
    for species, chemical_name in AROMA_SPECIES_TO_CHEMICAL.items():
        vp = mods["Chemical"](chemical_name, T=293.15).VaporPressure
        candidates = getattr(vp, "Antoine_parameters", {})
        if not candidates:
            raise UNIFACUnavailable(f"No Antoine vapor-pressure coefficients found for {chemical_name!r}.")
        ranked = []
        for method, pars in candidates.items():
            p_min = float(pars.get("Tmin", float("-inf")))
            p_max = float(pars.get("Tmax", float("inf")))
            covers = p_min <= t_min and p_max >= t_max
            gap = max(p_min - t_min, 0.0) + max(t_max - p_max, 0.0)
            ranked.append((0 if covers else 1, gap, str(method), pars))
        _cover_flag, gap, method, pars = sorted(ranked, key=lambda item: (item[0], item[1], item[2]))[0]
        data[species] = {
            "method": method,
            "chemical": chemical_name,
            "type": "Antoine",
            "A": float(pars["A"]),
            "B": float(pars["B"]),
            "C": float(pars["C"]),
            "base": float(pars["base"]),
            "Tmin": float(pars.get("Tmin", float("nan"))),
            "Tmax": float(pars.get("Tmax", float("nan"))),
            "process_min_temp_c": float(process_min_temp_c),
            "process_max_temp_c": float(process_max_temp_c),
            "extrapolated_for_process_window": bool(gap > 0.0),
            "range_gap_K": float(gap),
        }
    return data


def antoine_psat_pa(
    species: str,
    temp_c: float,
    process_min_temp_c: float = 10.0,
    process_max_temp_c: float = 30.0,
) -> tuple[float, dict[str, float | str | bool]]:
    """Return Antoine vapor pressure in Pa for an aroma species."""
    data = antoine_correlation_data(process_min_temp_c, process_max_temp_c)
    if species not in data:
        raise ValueError(f"Unknown aroma species {species!r}.")
    pars = data[species]
    temp_k = float(temp_c) + 273.15
    psat = float(pars["base"]) ** (float(pars["A"]) - float(pars["B"]) / (temp_k + float(pars["C"])))
    return float(psat), pars


def available_unifac_assignments() -> dict[str, dict[str, object]]:
    names = {
        "water": "water",
        "ethanol": "ethanol",
        "glucose": "glucose",
        "fructose": "fructose",
        **AROMA_SPECIES_TO_CHEMICAL,
    }
    rows = {}
    for key, name in names.items():
        try:
            meta = _metadata(name)
            groups = _groups(name)
            rows[key] = {
                "name": name,
                "CAS": meta.CASs,
                "InChI_key": meta.InChI_key,
                "groups": groups,
                "available": True,
            }
        except Exception as err:
            rows[key] = {
                "name": name,
                "available": False,
                "error": f"{type(err).__name__}: {err}",
            }
    return rows


def _liquid_component_moles_per_l(
    ethanol_g_l: float,
    glucose_g_l: float,
    fructose_g_l: float,
    mode: str,
) -> dict[str, float]:
    ethanol_g_l = max(float(ethanol_g_l), 0.0)
    glucose_g_l = max(float(glucose_g_l), 0.0)
    fructose_g_l = max(float(fructose_g_l), 0.0)
    sugar_g_l = glucose_g_l + fructose_g_l
    water_g_l = max(1.0, 1000.0 - ethanol_g_l - sugar_g_l)
    moles = {
        "water": water_g_l / _mw("water"),
        "ethanol": ethanol_g_l / _mw("ethanol"),
    }
    if mode == "water_ethanol":
        return moles
    if mode == "water_ethanol_glucose":
        moles["glucose"] = glucose_g_l / _mw("glucose")
        return moles
    if mode == "water_ethanol_total_sugar_as_glucose":
        moles["glucose"] = sugar_g_l / _mw("glucose")
        return moles
    raise ValueError(f"Unknown partition mode {mode!r}.")


def unifac_partition_K(
    species: str,
    temp_c: float,
    ethanol_g_l: float,
    glucose_g_l: float = 0.0,
    fructose_g_l: float = 0.0,
    mode: str = "water_ethanol_total_sugar_as_glucose",
) -> dict[str, float]:
    """Return a dilute gas/liquid concentration ratio from original UNIFAC.

    ``K`` is defined as ``C_gas / C_liquid`` with both concentrations in
    mol/L of their own phase. In the aroma model this is the coefficient used in
    ``loss_rate = Q_CO2 * K * A_liq``.
    """
    if species not in AROMA_SPECIES_TO_CHEMICAL:
        raise ValueError(f"Unknown aroma species {species!r}.")
    mods = _thermo_modules()
    chemical_name = AROMA_SPECIES_TO_CHEMICAL[species]
    temp_k = float(temp_c) + 273.15
    solvent_moles = _liquid_component_moles_per_l(ethanol_g_l, glucose_g_l, fructose_g_l, mode)
    component_names = list(solvent_moles) + [chemical_name]
    moles = [solvent_moles[name] for name in solvent_moles] + [TRACE_MOLES_PER_L]
    total_moles_l = float(sum(moles))
    xs = [n / total_moles_l for n in moles]
    chemgroups = [_groups(name) for name in component_names]
    gammas = mods["UNIFAC_gammas"](
        temp_k,
        xs,
        chemgroups,
        subgroup_data=mods["UFSG"],
        interaction_data=mods["UFIP"],
    )
    gamma = float(gammas[-1])
    vapor_pressure_pa = float(mods["Chemical"](chemical_name, T=temp_k).Psat)
    # p/(RT) is mol/m3 gas. Divide by 1000 to obtain mol/L gas.
    K = gamma * vapor_pressure_pa / (R_PA_M3_MOL_K * temp_k * 1000.0 * total_moles_l)
    return {
        "K": float(K),
        "gamma": gamma,
        "Psat_Pa": vapor_pressure_pa,
        "total_moles_per_l": total_moles_l,
        "temp_c": float(temp_c),
        "ethanol_g_l": float(ethanol_g_l),
        "glucose_g_l": float(glucose_g_l),
        "fructose_g_l": float(fructose_g_l),
        "mode": mode,
    }


def unifac_partition_K_antoine(
    species: str,
    temp_c: float,
    ethanol_g_l: float,
    glucose_g_l: float = 0.0,
    fructose_g_l: float = 0.0,
    mode: str = "water_ethanol",
    process_min_temp_c: float = 10.0,
    process_max_temp_c: float = 30.0,
) -> dict[str, float | str | bool]:
    """Return dilute gas/liquid ratio using Antoine vapor pressure and UNIFAC.

    ``K`` is defined as ``C_gas / C_liquid``. Activity coefficients are computed
    with original UNIFAC for the selected liquid mixture; vapor pressure is
    computed explicitly from Antoine coefficients.
    """
    if species not in AROMA_SPECIES_TO_CHEMICAL:
        raise ValueError(f"Unknown aroma species {species!r}.")
    mods = _thermo_modules()
    chemical_name = AROMA_SPECIES_TO_CHEMICAL[species]
    temp_k = float(temp_c) + 273.15
    solvent_moles = _liquid_component_moles_per_l(ethanol_g_l, glucose_g_l, fructose_g_l, mode)
    component_names = list(solvent_moles) + [chemical_name]
    moles = [solvent_moles[name] for name in solvent_moles] + [TRACE_MOLES_PER_L]
    total_moles_l = float(sum(moles))
    xs = [n / total_moles_l for n in moles]
    chemgroups = [_groups(name) for name in component_names]
    gammas = mods["UNIFAC_gammas"](
        temp_k,
        xs,
        chemgroups,
        subgroup_data=mods["UFSG"],
        interaction_data=mods["UFIP"],
    )
    gamma = float(gammas[-1])
    vapor_pressure_pa, antoine = antoine_psat_pa(
        species,
        temp_c,
        process_min_temp_c=process_min_temp_c,
        process_max_temp_c=process_max_temp_c,
    )
    K = gamma * vapor_pressure_pa / (R_PA_M3_MOL_K * temp_k * 1000.0 * total_moles_l)
    return {
        "K": float(K),
        "gamma": gamma,
        "Psat_Pa": float(vapor_pressure_pa),
        "total_moles_per_l": total_moles_l,
        "temp_c": float(temp_c),
        "ethanol_g_l": float(ethanol_g_l),
        "glucose_g_l": float(glucose_g_l),
        "fructose_g_l": float(fructose_g_l),
        "mode": mode,
        "vapor_pressure_model": "Antoine",
        "antoine_method": str(antoine["method"]),
        "antoine_Tmin_K": float(antoine["Tmin"]),
        "antoine_Tmax_K": float(antoine["Tmax"]),
        "antoine_extrapolated_for_process_window": bool(antoine["extrapolated_for_process_window"]),
    }


def fit_partition_surrogate(
    species: str,
    mode: str = "water_ethanol_total_sugar_as_glucose",
    temp_grid_c: tuple[float, ...] = (12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 28.0),
    ethanol_grid_g_l: tuple[float, ...] = (0.0, 20.0, 50.0, 80.0, 120.0, 150.0),
    sugar_grid_g_l: tuple[float, ...] = (0.0, 50.0, 100.0, 180.0, 260.0),
) -> dict[str, float]:
    rows = []
    for temp_c in temp_grid_c:
        for ethanol_g_l in ethanol_grid_g_l:
            sugar_values = (0.0,) if mode == "water_ethanol" else sugar_grid_g_l
            for sugar_g_l in sugar_values:
                out = unifac_partition_K(
                    species,
                    temp_c,
                    ethanol_g_l,
                    glucose_g_l=sugar_g_l,
                    fructose_g_l=0.0,
                    mode=mode,
                )
                rows.append((temp_c, ethanol_g_l, sugar_g_l, math.log(max(out["K"], 1e-30))))

    X = []
    y = []
    for temp_c, ethanol_g_l, sugar_g_l, logK in rows:
        X.append([1.0, temp_c - 20.0, ethanol_g_l - 50.0, sugar_g_l - 100.0])
        y.append(logK)
    X_arr = np.asarray(X, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    beta, *_ = np.linalg.lstsq(X_arr, y_arr, rcond=None)
    pred = X_arr @ beta
    residual = pred - y_arr
    return {
        "K20": float(math.exp(beta[0])),
        "logK_ref": float(beta[0]),
        "temp_slope": float(beta[1]),
        "ethanol_slope": float(beta[2]),
        "sugar_slope": float(beta[3]),
        "fit_rmse_log": float(np.sqrt(np.mean(residual * residual))),
        "fit_max_abs_log_error": float(np.max(np.abs(residual))),
        "mode": mode,
        "source": "thermo original UNIFAC, dilute aroma, fitted log-linear surrogate",
        "trap_efficiency": float(TRAP_EFFICIENCY[species]),
    }


def fit_all_partition_surrogates(mode: str = "water_ethanol_total_sugar_as_glucose") -> dict[str, dict[str, float]]:
    return {species: fit_partition_surrogate(species, mode=mode) for species in AROMA_SPECIES_TO_CHEMICAL}


def write_partition_report(path: Path, mode: str = "water_ethanol_total_sugar_as_glucose") -> None:
    assignments = available_unifac_assignments()
    coefficients = fit_all_partition_surrogates(mode)
    payload = {
        "mode": mode,
        "assignments": assignments,
        "coefficients": coefficients,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
