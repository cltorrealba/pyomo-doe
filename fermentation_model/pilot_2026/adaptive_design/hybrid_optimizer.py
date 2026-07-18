from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.stats import qmc


Objective = Callable[[np.ndarray], float]


@dataclass(frozen=True)
class ParticleSwarmResult:
    x: np.ndarray
    fun: float
    history: np.ndarray
    evaluations: int
    initial_positions: np.ndarray
    initial_values: np.ndarray
    initial_sources: tuple[str, ...]
    final_positions: np.ndarray
    final_values: np.ndarray
    personal_best_positions: np.ndarray
    personal_best_values: np.ndarray
    diversity_history: np.ndarray
    improvement_history: np.ndarray
    position_history: np.ndarray
    value_history: np.ndarray
    converged: bool
    stop_reason: str
    iterations_completed: int
    restarts_completed: int


def _sobol_positions(count: int, dimension: int, seed: int) -> np.ndarray:
    if count <= 0:
        return np.empty((0, dimension), dtype=float)
    sampler = qmc.Sobol(d=dimension, scramble=True, seed=int(seed))
    power = int(math.ceil(math.log2(max(count, 1))))
    return sampler.random_base2(power)[:count]


def _normalized_diversity(positions: np.ndarray, lower: np.ndarray, span: np.ndarray) -> float:
    normalized = (positions - lower) / span
    return float(np.mean(np.std(normalized, axis=0, ddof=0)))


def particle_swarm(
    objective: Objective,
    bounds: np.ndarray,
    *,
    particles: int,
    iterations: int,
    seed: int,
    inertia: float = 0.72,
    cognitive: float = 1.49,
    social: float = 1.49,
    initial_positions: np.ndarray | None = None,
    stagnation_iterations: int | None = None,
    improvement_tolerance: float = 0.0,
    restarts: int = 0,
) -> ParticleSwarmResult:
    bounds = np.asarray(bounds, dtype=float)
    if bounds.ndim != 2 or bounds.shape[1] != 2 or np.any(bounds[:, 0] >= bounds[:, 1]):
        raise ValueError("bounds must be a finite n-by-2 array with lower < upper")
    if not np.isfinite(bounds).all() or particles < 2 or iterations < 1:
        raise ValueError("invalid particle-swarm configuration")
    if stagnation_iterations is not None and stagnation_iterations < 1:
        raise ValueError("stagnation_iterations must be positive when supplied")
    if improvement_tolerance < 0.0 or restarts < 0:
        raise ValueError("improvement_tolerance and restarts must be nonnegative")
    rng = np.random.default_rng(int(seed))
    lower, upper = bounds[:, 0], bounds[:, 1]
    span = upper - lower
    seeds = (
        np.empty((0, len(bounds)), dtype=float)
        if initial_positions is None
        else np.atleast_2d(np.asarray(initial_positions, dtype=float))
    )
    if seeds.shape[1] != len(bounds) or len(seeds) > particles:
        raise ValueError("initial_positions must have at most particles rows and match bounds")
    seeds = np.clip(seeds, lower, upper)
    sobol_count = particles - len(seeds)
    sobol = lower + _sobol_positions(sobol_count, len(bounds), int(seed)) * span
    positions = np.vstack([seeds, sobol])
    sources = tuple(["seed_profile"] * len(seeds) + ["sobol"] * sobol_count)
    velocities = rng.uniform(-0.1, 0.1, size=positions.shape) * span
    values = np.asarray([objective(row.copy()) for row in positions], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("objective returned non-finite values")
    initial_position_values = positions.copy()
    initial_values = values.copy()
    personal_positions = positions.copy()
    personal_values = values.copy()
    best_index = int(np.argmin(values))
    global_position = positions[best_index].copy()
    global_value = float(values[best_index])
    history = [global_value]
    improvement_history = [0.0]
    diversity_history = [_normalized_diversity(positions, lower, span)]
    position_history = [positions.copy()]
    value_history = [values.copy()]
    evaluations = particles
    stagnant = 0
    restarts_completed = 0
    stop_reason = "maximum_iterations"
    for iteration in range(1, iterations + 1):
        previous_global = global_value
        r1 = rng.random(positions.shape)
        r2 = rng.random(positions.shape)
        velocities = (
            inertia * velocities
            + cognitive * r1 * (personal_positions - positions)
            + social * r2 * (global_position - positions)
        )
        velocities = np.clip(velocities, -0.5 * span, 0.5 * span)
        positions = np.clip(positions + velocities, lower, upper)
        values = np.asarray([objective(row.copy()) for row in positions], dtype=float)
        evaluations += particles
        improved = values < personal_values
        personal_positions[improved] = positions[improved]
        personal_values[improved] = values[improved]
        best_index = int(np.argmin(personal_values))
        if personal_values[best_index] < global_value:
            global_value = float(personal_values[best_index])
            global_position = personal_positions[best_index].copy()
        improvement = max(previous_global - global_value, 0.0)
        stagnant = stagnant + 1 if improvement <= improvement_tolerance else 0
        history.append(global_value)
        improvement_history.append(improvement)
        diversity_history.append(_normalized_diversity(positions, lower, span))
        position_history.append(positions.copy())
        value_history.append(values.copy())
        if stagnation_iterations is not None and stagnant >= stagnation_iterations:
            if restarts_completed < restarts:
                restart_count = max(1, particles // 4)
                worst = np.argsort(personal_values)[-restart_count:]
                fresh = lower + _sobol_positions(
                    restart_count,
                    len(bounds),
                    int(seed) + 104729 * (restarts_completed + 1),
                ) * span
                positions[worst] = fresh
                velocities[worst] = rng.uniform(-0.1, 0.1, size=(restart_count, len(bounds))) * span
                values[worst] = np.asarray([objective(row.copy()) for row in fresh], dtype=float)
                evaluations += restart_count
                personal_positions[worst] = fresh
                personal_values[worst] = values[worst]
                best_index = int(np.argmin(personal_values))
                if personal_values[best_index] < global_value:
                    global_value = float(personal_values[best_index])
                    global_position = personal_positions[best_index].copy()
                restarts_completed += 1
                stagnant = 0
                continue
            stop_reason = "stagnation"
            break
    iterations_completed = len(history) - 1
    converged = bool(
        stop_reason == "stagnation"
        and iterations_completed > 0
        and improvement_history[-1] <= improvement_tolerance
    )
    return ParticleSwarmResult(
        x=global_position,
        fun=global_value,
        history=np.asarray(history),
        evaluations=evaluations,
        initial_positions=initial_position_values,
        initial_values=initial_values,
        initial_sources=sources,
        final_positions=positions.copy(),
        final_values=values.copy(),
        personal_best_positions=personal_positions.copy(),
        personal_best_values=personal_values.copy(),
        diversity_history=np.asarray(diversity_history),
        improvement_history=np.asarray(improvement_history),
        position_history=np.asarray(position_history),
        value_history=np.asarray(value_history),
        converged=converged,
        stop_reason=stop_reason,
        iterations_completed=iterations_completed,
        restarts_completed=restarts_completed,
    )
