from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


Objective = Callable[[np.ndarray], float]


@dataclass(frozen=True)
class ParticleSwarmResult:
    x: np.ndarray
    fun: float
    history: np.ndarray
    evaluations: int


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
) -> ParticleSwarmResult:
    bounds = np.asarray(bounds, dtype=float)
    if bounds.ndim != 2 or bounds.shape[1] != 2 or np.any(bounds[:, 0] >= bounds[:, 1]):
        raise ValueError("bounds must be a finite n-by-2 array with lower < upper")
    if not np.isfinite(bounds).all() or particles < 2 or iterations < 1:
        raise ValueError("invalid particle-swarm configuration")
    rng = np.random.default_rng(int(seed))
    lower, upper = bounds[:, 0], bounds[:, 1]
    span = upper - lower
    positions = lower + rng.random((particles, len(bounds))) * span
    velocities = rng.uniform(-0.1, 0.1, size=positions.shape) * span
    values = np.asarray([objective(row.copy()) for row in positions], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("objective returned non-finite values")
    personal_positions = positions.copy()
    personal_values = values.copy()
    best_index = int(np.argmin(values))
    global_position = positions[best_index].copy()
    global_value = float(values[best_index])
    history = [global_value]
    evaluations = particles
    for _ in range(iterations):
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
        history.append(global_value)
    return ParticleSwarmResult(global_position, global_value, np.asarray(history), evaluations)
