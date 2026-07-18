from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np


Objective = Callable[[np.ndarray], float]


@dataclass(frozen=True)
class LocalRefinementResult:
    x: np.ndarray
    fun: float
    state: str
    qualified: bool
    accepted_improvement: bool
    trace: tuple[dict[str, object], ...]
    rejected_candidates: tuple[dict[str, object], ...]


def _derivatives(
    objective: Objective,
    base_vector: np.ndarray,
    continuous_indices: np.ndarray,
    bounds: np.ndarray,
    z: np.ndarray,
    step: float,
) -> tuple[float, np.ndarray, np.ndarray, list[str]]:
    lower = bounds[continuous_indices, 0]
    span = bounds[continuous_indices, 1] - lower

    def full_vector(normalized: np.ndarray) -> np.ndarray:
        candidate = base_vector.copy()
        candidate[continuous_indices] = lower + normalized * span
        return candidate

    cache: dict[tuple[float, ...], float] = {}

    def evaluate(normalized: np.ndarray) -> float:
        key = tuple(np.round(normalized, 14))
        if key not in cache:
            cache[key] = float(objective(full_vector(normalized)))
        return cache[key]

    f0 = evaluate(z)
    gradient = np.zeros(len(z), dtype=float)
    curvature = np.zeros(len(z), dtype=float)
    methods: list[str] = []
    for index in range(len(z)):
        can_minus = z[index] - step >= -1e-14
        can_plus = z[index] + step <= 1.0 + 1e-14
        if can_minus and can_plus:
            plus, minus = z.copy(), z.copy()
            plus[index] += step
            minus[index] -= step
            fp, fm = evaluate(plus), evaluate(minus)
            gradient[index] = (fp - fm) / (2.0 * step)
            curvature[index] = (fp - 2.0 * f0 + fm) / step**2
            methods.append("central")
        elif can_plus:
            h = min(step, (1.0 - z[index]) / 2.0)
            if h <= 0.0:
                raise ValueError("No feasible forward derivative step")
            one, two = z.copy(), z.copy()
            one[index] += h
            two[index] += 2.0 * h
            f1, f2 = evaluate(one), evaluate(two)
            gradient[index] = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
            curvature[index] = (f0 - 2.0 * f1 + f2) / h**2
            methods.append("forward_second_order")
        elif can_minus:
            h = min(step, z[index] / 2.0)
            if h <= 0.0:
                raise ValueError("No feasible backward derivative step")
            one, two = z.copy(), z.copy()
            one[index] -= h
            two[index] -= 2.0 * h
            f1, f2 = evaluate(one), evaluate(two)
            gradient[index] = (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h)
            curvature[index] = (f0 - 2.0 * f1 + f2) / h**2
            methods.append("backward_second_order")
        else:  # pragma: no cover - normalized bounds always have one direction
            raise ValueError("No feasible derivative direction")
    return f0, gradient, curvature, methods


def sequential_ipopt_refine(
    x: np.ndarray,
    bounds: np.ndarray,
    objective: Objective,
    continuous_indices: np.ndarray,
    config: dict,
    *,
    candidate_id: str,
    validation_objective: Objective | None = None,
) -> LocalRefinementResult:
    """Sequential trust-region refinement with full-objective step acceptance."""

    full_objective = objective if validation_objective is None else validation_objective

    try:
        import pyomo.environ as pyo
    except ImportError:
        return LocalRefinementResult(
            np.asarray(x, dtype=float),
            float(full_objective(np.asarray(x, dtype=float))),
            "not_executed",
            False,
            False,
            tuple(),
            tuple(),
        )
    solver = pyo.SolverFactory("ipopt")
    if not solver.available(exception_flag=False):
        return LocalRefinementResult(
            np.asarray(x, dtype=float),
            float(full_objective(np.asarray(x, dtype=float))),
            "not_executed",
            False,
            False,
            tuple(),
            tuple(),
        )
    x = np.asarray(x, dtype=float).copy()
    bounds = np.asarray(bounds, dtype=float)
    continuous_indices = np.asarray(continuous_indices, dtype=int)
    lower = bounds[continuous_indices, 0]
    span = bounds[continuous_indices, 1] - lower
    z = np.clip((x[continuous_indices] - lower) / span, 0.0, 1.0)
    search = config["search"]
    derivative_step = float(search["local_finite_difference_step"])
    radius = float(search["local_initial_trust_radius"])
    minimum_radius = float(search["local_minimum_trust_radius"])
    maximum_radius = float(search["local_maximum_trust_radius"])
    maximum_iterations = int(search["local_maximum_iterations"])
    acceptance_ratio = float(search["local_acceptance_ratio"])
    gradient_tolerance = float(search["local_gradient_tolerance"])
    improvement_tolerance = float(search["improvement_tolerance"])
    trace: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    accepted_any = False
    previous_z: np.ndarray | None = None
    previous_gradient: np.ndarray | None = None
    previous_hessian: np.ndarray | None = None
    state = "rejected_model_mismatch"
    full_current = float(full_objective(x.copy()))
    for iteration in range(1, maximum_iterations + 1):
        base = x.copy()
        f0, gradient, diagonal_curvature, methods = _derivatives(
            objective,
            base,
            continuous_indices,
            bounds,
            z,
            derivative_step,
        )
        gradient_norm = float(np.linalg.norm(gradient, ord=np.inf))
        bfgs_status = "initial_diagonal_with_signed_curvature"
        hessian = np.diag(diagonal_curvature)
        if (
            previous_z is not None
            and previous_gradient is not None
            and previous_hessian is not None
        ):
            s = z - previous_z
            y = gradient - previous_gradient
            sty = float(s @ y)
            hs = previous_hessian @ s
            sths = float(s @ hs)
            if sty > 1e-12 and abs(sths) > 1e-12:
                hessian = previous_hessian - np.outer(hs, hs) / sths + np.outer(y, y) / sty
                bfgs_status = "bfgs_cross_terms_updated"
            else:
                hessian = previous_hessian
                bfgs_status = "bfgs_update_skipped_nonpositive_or_degenerate_curvature"
        if gradient_norm <= gradient_tolerance:
            trace.append(
                {
                    "candidate_id": candidate_id,
                    "iteration": iteration,
                    "state": "stationary_no_improving_step",
                    "actual_objective_before": f0,
                    "actual_objective_after": f0,
                    "full_ensemble_objective_before": full_current,
                    "full_ensemble_objective_after": full_current,
                    "reduced_model_acceptance_ratio": math.nan,
                    "full_ensemble_acceptance_ratio": math.nan,
                    "gradient_infinity_norm": gradient_norm,
                    "trust_radius": radius,
                    "derivative_methods": ";".join(methods),
                    "negative_coordinate_curvatures": int(np.sum(diagonal_curvature < 0.0)),
                    "hessian_update": bfgs_status,
                    "ipopt_termination": "not_needed_stationary",
                }
            )
            state = "stationary_no_improving_step"
            break
        model = pyo.ConcreteModel()
        model.I = pyo.RangeSet(0, len(z) - 1)

        def step_bounds(_model, index):
            index = int(index)
            return (
                max(-radius, -float(z[index])),
                min(radius, 1.0 - float(z[index])),
            )

        model.d = pyo.Var(model.I, bounds=step_bounds, initialize=0.0)
        model.objective = pyo.Objective(
            expr=f0
            + sum(float(gradient[i]) * model.d[i] for i in range(len(z)))
            + 0.5
            * sum(
                float(hessian[i, j]) * model.d[i] * model.d[j]
                for i in range(len(z))
                for j in range(len(z))
            )
        )
        try:
            result = solver.solve(
                model,
                tee=False,
                options={"tol": 1e-8, "max_iter": 500, "print_level": 0},
            )
            termination = str(result.solver.termination_condition)
            solver_iterations = getattr(result.solver, "iterations", None)
            step_vector = np.asarray([float(pyo.value(model.d[i])) for i in range(len(z))])
        except Exception as exc:  # pragma: no cover - solver failure is environment-specific
            trace.append(
                {
                    "candidate_id": candidate_id,
                    "iteration": iteration,
                    "state": "solver_failure",
                    "actual_objective_before": f0,
                    "gradient_infinity_norm": gradient_norm,
                    "trust_radius": radius,
                    "derivative_methods": ";".join(methods),
                    "negative_coordinate_curvatures": int(np.sum(diagonal_curvature < 0.0)),
                    "hessian_update": bfgs_status,
                    "ipopt_termination": repr(exc),
                }
            )
            state = "solver_failure"
            break
        proposed_z = np.clip(z + step_vector, 0.0, 1.0)
        candidate = x.copy()
        candidate[continuous_indices] = lower + proposed_z * span
        reduced_actual = float(objective(candidate))
        full_actual = float(full_objective(candidate))
        surrogate_after = float(
            f0 + gradient @ step_vector + 0.5 * step_vector @ hessian @ step_vector
        )
        predicted_improvement = float(f0 - surrogate_after)
        reduced_improvement = float(f0 - reduced_actual)
        full_improvement = float(full_current - full_actual)
        reduced_ratio = (
            reduced_improvement / predicted_improvement
            if predicted_improvement > improvement_tolerance
            else -math.inf
        )
        full_ratio = (
            full_improvement / predicted_improvement
            if predicted_improvement > improvement_tolerance
            else -math.inf
        )
        accepted = bool(
            full_improvement > improvement_tolerance and full_ratio >= acceptance_ratio
        )
        row = {
            "candidate_id": candidate_id,
            "iteration": iteration,
            "state": "accepted_improvement" if accepted else "rejected_model_mismatch",
            "actual_objective_before": f0,
            "surrogate_objective_after": surrogate_after,
            "actual_objective_after": reduced_actual,
            "predicted_improvement": predicted_improvement,
            "actual_improvement": reduced_improvement,
            "acceptance_ratio": full_ratio,
            "reduced_model_acceptance_ratio": reduced_ratio,
            "full_ensemble_acceptance_ratio": full_ratio,
            "full_ensemble_objective_before": full_current,
            "full_ensemble_objective_after": full_actual,
            "full_ensemble_actual_improvement": full_improvement,
            "accepted": accepted,
            "gradient_infinity_norm": gradient_norm,
            "trust_radius": radius,
            "step_infinity_norm": float(np.linalg.norm(step_vector, ord=np.inf)),
            "derivative_methods": ";".join(methods),
            "negative_coordinate_curvatures": int(np.sum(diagonal_curvature < 0.0)),
            "minimum_hessian_eigenvalue": float(np.linalg.eigvalsh(0.5 * (hessian + hessian.T))[0]),
            "hessian_update": bfgs_status,
            "ipopt_termination": termination,
            "ipopt_iterations": solver_iterations,
            "infeasibility": None,
            "kkt_error": None,
        }
        trace.append(row)
        if accepted:
            previous_z = z.copy()
            previous_gradient = gradient.copy()
            previous_hessian = hessian.copy()
            z = proposed_z
            x = candidate
            full_current = full_actual
            accepted_any = True
            state = "accepted_improvement"
            if full_ratio >= float(search["local_expansion_ratio"]) and np.linalg.norm(
                step_vector, ord=np.inf
            ) >= 0.8 * radius:
                radius = min(maximum_radius, 2.0 * radius)
        else:
            rejected.append({**row, "candidate_vector": candidate.tolist()})
            radius *= 0.5
            if radius < minimum_radius:
                state = "rejected_model_mismatch"
                break
    final_fun = float(full_objective(x))
    qualified = state in {"accepted_improvement", "stationary_no_improving_step"}
    return LocalRefinementResult(
        x=x,
        fun=final_fun,
        state=state,
        qualified=qualified,
        accepted_improvement=accepted_any,
        trace=tuple(trace),
        rejected_candidates=tuple(rejected),
    )
