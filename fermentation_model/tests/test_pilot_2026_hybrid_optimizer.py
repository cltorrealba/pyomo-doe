from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


FERMENTATION_DIR = Path(__file__).resolve().parents[1]
if str(FERMENTATION_DIR) not in sys.path:
    sys.path.insert(0, str(FERMENTATION_DIR))

from pilot_2026.adaptive_design.hybrid_optimizer import particle_swarm  # noqa: E402


class HybridOptimizerTests(unittest.TestCase):
    def test_particle_swarm_is_reproducible_and_improves(self) -> None:
        objective = lambda x: float(np.dot(x - 0.25, x - 0.25))
        bounds = np.asarray([[-2.0, 2.0], [-2.0, 2.0]])
        kwargs = dict(particles=16, iterations=30, seed=17)
        first = particle_swarm(objective, bounds, **kwargs)
        second = particle_swarm(objective, bounds, **kwargs)
        np.testing.assert_allclose(first.x, second.x, atol=0.0, rtol=0.0)
        self.assertEqual(first.fun, second.fun)
        self.assertLess(first.history[-1], first.history[0])
        self.assertTrue(np.all(np.diff(first.history) <= 0.0))

    def test_invalid_bounds_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "bounds"):
            particle_swarm(lambda x: 0.0, np.asarray([[1.0, 1.0]]), particles=4, iterations=2, seed=1)

    def test_seeded_position_is_part_of_initial_population(self) -> None:
        objective = lambda x: float(np.dot(x, x))
        bounds = np.asarray([[-2.0, 2.0], [-2.0, 2.0]])
        result = particle_swarm(
            objective,
            bounds,
            particles=4,
            iterations=1,
            seed=1,
            initial_positions=np.asarray([[0.0, 0.0]]),
        )
        self.assertEqual(result.fun, 0.0)


if __name__ == "__main__":
    unittest.main()
