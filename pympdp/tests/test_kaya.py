import math

from mpdp.logger import logger
from mpdp.mpmd import run_example

# Silence the verbose INFO logging from the DP solver during unit tests.
logger.set_warning()

# Each entry captures (example name, discretization, refinement, absolute tolerance in meters).
TEST_CASES = [
    ("Kaya Example 1", 16, 4, 0.05),
    ("Kaya Example 2", 16, 4, 0.06),
    ("Kaya Example 3", 32, 4, 0.01),
    ("Kaya Example 4", 32, 4, 0.02),
    ("Omega", 60, 1, 0.12),
    ("Circuit", 8, 1, 0.05),
]


def test_examples_match_reference_lengths():
    """Ensure each translated reference scenario stays close to the published solution length."""
    for name, discr, refin, tol in TEST_CASES:
        result = run_example(name, discretizations=[discr], refinements=[refin])[0]
        assert math.isfinite(result.path_length), f"{name} produced a non-finite length"
        assert abs(result.diff) <= tol, (
            f"{name}: length differs by {result.diff:.6f} m (allowed {tol} m)"
        )
