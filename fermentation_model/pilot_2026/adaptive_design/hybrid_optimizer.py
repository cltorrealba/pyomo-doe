from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Callable

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


@dataclass
class CheckpointedSwarmState:
    """Complete mutable PSO state required for bitwise reproducible continuation."""

    seed: int
    fidelity: str
    iteration: int
    evaluations: int
    positions: np.ndarray
    velocities: np.ndarray
    values: np.ndarray
    personal_best_positions: np.ndarray
    personal_best_values: np.ndarray
    global_best_position: np.ndarray
    global_best_value: float
    rng_state: dict[str, Any]
    history: list[float]
    improvement_history: list[float]

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "seed": int(self.seed),
            "fidelity": self.fidelity,
            "iteration": int(self.iteration),
            "evaluations": int(self.evaluations),
            "positions": self.positions.tolist(),
            "velocities": self.velocities.tolist(),
            "values": self.values.tolist(),
            "personal_best_positions": self.personal_best_positions.tolist(),
            "personal_best_values": self.personal_best_values.tolist(),
            "global_best_position": self.global_best_position.tolist(),
            "global_best_value": float(self.global_best_value),
            "rng_state": self.rng_state,
            "history": [float(value) for value in self.history],
            "improvement_history": [
                float(value) for value in self.improvement_history
            ],
        }

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "CheckpointedSwarmState":
        if int(payload.get("schema_version", -1)) != 1:
            raise ValueError("Unsupported checkpointed swarm schema")
        return cls(
            seed=int(payload["seed"]),
            fidelity=str(payload["fidelity"]),
            iteration=int(payload["iteration"]),
            evaluations=int(payload["evaluations"]),
            positions=np.asarray(payload["positions"], dtype=float),
            velocities=np.asarray(payload["velocities"], dtype=float),
            values=np.asarray(payload["values"], dtype=float),
            personal_best_positions=np.asarray(
                payload["personal_best_positions"], dtype=float
            ),
            personal_best_values=np.asarray(payload["personal_best_values"], dtype=float),
            global_best_position=np.asarray(payload["global_best_position"], dtype=float),
            global_best_value=float(payload["global_best_value"]),
            rng_state=dict(payload["rng_state"]),
            history=[float(value) for value in payload["history"]],
            improvement_history=[
                float(value) for value in payload["improvement_history"]
            ],
        )


def _sobol_positions(count: int, dimension: int, seed: int) -> np.ndarray:
    if count <= 0:
        return np.empty((0, dimension), dtype=float)
    sampler = qmc.Sobol(d=dimension, scramble=True, seed=int(seed))
    power = int(math.ceil(math.log2(max(count, 1))))
    return sampler.random_base2(power)[:count]


def _normalized_diversity(positions: np.ndarray, lower: np.ndarray, span: np.ndarray) -> float:
    normalized = (positions - lower) / span
    return float(np.mean(np.std(normalized, axis=0, ddof=0)))


def initialize_checkpointed_swarm(
    objective: Objective,
    bounds: np.ndarray,
    *,
    particles: int,
    seed: int,
    fidelity: str,
    initial_positions: np.ndarray | None = None,
) -> CheckpointedSwarmState:
    """Initialize a PSO state without hiding any state needed for resume."""

    bounds = np.asarray(bounds, dtype=float)
    if bounds.ndim != 2 or bounds.shape[1] != 2 or np.any(bounds[:, 0] >= bounds[:, 1]):
        raise ValueError("bounds must be a finite n-by-2 array with lower < upper")
    if particles < 2 or not np.isfinite(bounds).all():
        raise ValueError("invalid checkpointed particle-swarm configuration")
    lower, upper = bounds[:, 0], bounds[:, 1]
    span = upper - lower
    seeds = (
        np.empty((0, len(bounds)), dtype=float)
        if initial_positions is None
        else np.atleast_2d(np.asarray(initial_positions, dtype=float))
    )
    if seeds.shape[1] != len(bounds) or len(seeds) > particles:
        raise ValueError("initial_positions must match bounds and particle count")
    seeds = np.clip(seeds, lower, upper)
    sobol_count = particles - len(seeds)
    positions = np.vstack(
        [seeds, lower + _sobol_positions(sobol_count, len(bounds), int(seed)) * span]
    )
    rng = np.random.default_rng(int(seed))
    velocities = rng.uniform(-0.1, 0.1, size=positions.shape) * span
    values = np.asarray([objective(row.copy()) for row in positions], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("objective returned non-finite values")
    best = int(np.argmin(values))
    return CheckpointedSwarmState(
        seed=int(seed),
        fidelity=str(fidelity),
        iteration=0,
        evaluations=particles,
        positions=positions,
        velocities=velocities,
        values=values,
        personal_best_positions=positions.copy(),
        personal_best_values=values.copy(),
        global_best_position=positions[best].copy(),
        global_best_value=float(values[best]),
        rng_state=rng.bit_generator.state,
        history=[float(values[best])],
        improvement_history=[0.0],
    )


def rescore_checkpointed_swarm(
    state: CheckpointedSwarmState,
    objective: Objective,
    *,
    fidelity: str,
) -> CheckpointedSwarmState:
    """Continue the same swarm at a new fidelity after rescoring its memory."""

    positions_values = np.asarray(
        [objective(row.copy()) for row in state.positions], dtype=float
    )
    personal_values = np.asarray(
        [objective(row.copy()) for row in state.personal_best_positions], dtype=float
    )
    choose_current = positions_values < personal_values
    personal_positions = state.personal_best_positions.copy()
    personal_positions[choose_current] = state.positions[choose_current]
    personal_values[choose_current] = positions_values[choose_current]
    best = int(np.argmin(personal_values))
    return CheckpointedSwarmState(
        seed=state.seed,
        fidelity=str(fidelity),
        iteration=state.iteration,
        evaluations=state.evaluations + 2 * len(state.positions),
        positions=state.positions.copy(),
        velocities=state.velocities.copy(),
        values=positions_values,
        personal_best_positions=personal_positions,
        personal_best_values=personal_values,
        global_best_position=personal_positions[best].copy(),
        global_best_value=float(personal_values[best]),
        rng_state=state.rng_state,
        history=state.history + [float(personal_values[best])],
        improvement_history=state.improvement_history + [0.0],
    )


def advance_checkpointed_swarm(
    state: CheckpointedSwarmState,
    objective: Objective,
    bounds: np.ndarray,
    *,
    inertia: float = 0.72,
    cognitive: float = 1.49,
    social: float = 1.49,
) -> CheckpointedSwarmState:
    """Advance exactly one iteration so callers can checkpoint every iteration."""

    bounds = np.asarray(bounds, dtype=float)
    lower, upper = bounds[:, 0], bounds[:, 1]
    span = upper - lower
    if state.positions.shape != state.velocities.shape or state.positions.shape[1] != len(bounds):
        raise ValueError("Checkpointed swarm arrays do not match configured bounds")
    rng = np.random.default_rng()
    rng.bit_generator.state = state.rng_state
    r1 = rng.random(state.positions.shape)
    r2 = rng.random(state.positions.shape)
    velocities = (
        float(inertia) * state.velocities
        + float(cognitive)
        * r1
        * (state.personal_best_positions - state.positions)
        + float(social)
        * r2
        * (state.global_best_position - state.positions)
    )
    velocities = np.clip(velocities, -0.5 * span, 0.5 * span)
    positions = np.clip(state.positions + velocities, lower, upper)
    values = np.asarray([objective(row.copy()) for row in positions], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("objective returned non-finite values")
    personal_positions = state.personal_best_positions.copy()
    personal_values = state.personal_best_values.copy()
    improved = values < personal_values
    personal_positions[improved] = positions[improved]
    personal_values[improved] = values[improved]
    best = int(np.argmin(personal_values))
    global_value = float(personal_values[best])
    improvement = max(float(state.global_best_value) - global_value, 0.0)
    return CheckpointedSwarmState(
        seed=state.seed,
        fidelity=state.fidelity,
        iteration=state.iteration + 1,
        evaluations=state.evaluations + len(positions),
        positions=positions,
        velocities=velocities,
        values=values,
        personal_best_positions=personal_positions,
        personal_best_values=personal_values,
        global_best_position=personal_positions[best].copy(),
        global_best_value=global_value,
        rng_state=rng.bit_generator.state,
        history=state.history + [global_value],
        improvement_history=state.improvement_history + [improvement],
    )


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
