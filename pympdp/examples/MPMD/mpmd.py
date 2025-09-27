from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

_PACKAGE_ROOT = Path(__file__).resolve().parent
if str(_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(_PACKAGE_ROOT))

from dp import DP
from dubins import dubins_shortest_path
from mpmd_data import EXAMPLE_RAW_DATA

DEFAULT_DISCRETIZATIONS: Sequence[int] = (90, )#(4, 16, 90, 360)
DEFAULT_REFINEMENTS: Sequence[int] = (4, )#(1, 2, 4, 8, 16)

_SPEC_TABLE: Tuple[Tuple[str, str, float, float], ...] = (
    ("kaya1", "Kaya Example 1", 3.0, 3.4155788580751487),
    ("kaya2", "Kaya Example 2", 3.0, 6.2780345503093136),
    ("kaya3", "Kaya Example 3", 5.0, 11.916212654285486),
    ("kaya4", "Kaya Example 4", 3.0, 7.4675621973384265),
    # ("omega", "Omega", 3.0, 41.07250164388393),
    # ("spa", "Circuit", 3.0, 6988.66098639943),
)


@dataclass
class ExampleSpec:
    key: str
    name: str
    k_max: float
    example_length: float
    raw_points: List[Tuple[float, float, Optional[float]]]

    def __post_init__(self) -> None:  # eager pre-computation for downstream speed
        self.points: List[Tuple[float, float]] = [
            (x, y) for x, y, _ in self.raw_points
        ]
        self.fixed_angles: List[bool] = [
            theta is not None for _, _, theta in self.raw_points
        ]
        self.def_thetas: List[float] = [
            float(theta) if theta is not None else 0.0
            for _, _, theta in self.raw_points
        ]
        self.theta_hints: List[Optional[float]] = [
            theta for _, _, theta in self.raw_points
        ]


@dataclass
class ExampleResult:
    example: str
    discretization: int
    refinement: int
    path_length: float
    expected_length: float
    diff: float
    runtime_ms: float
    angles: List[float]

    @property
    def diff_mm(self) -> float:
        return self.diff * 1000.0


def _build_example_specs() -> Dict[str, ExampleSpec]:
    specs: Dict[str, ExampleSpec] = {}
    for key, name, k_max, example_length in _SPEC_TABLE:
        raw = EXAMPLE_RAW_DATA[key]
        specs[name] = ExampleSpec(
            key=key,
            name=name,
            k_max=k_max,
            example_length=example_length,
            raw_points=[(float(x), float(y), theta if theta is None else float(theta)) for x, y, theta in raw],
        )
    return specs


EXAMPLES_BY_NAME = _build_example_specs()
EXAMPLES_BY_KEY = {spec.key: spec for spec in EXAMPLES_BY_NAME.values()}
EXAMPLE_ORDER: Tuple[str, ...] = tuple(name for _, name, _, _ in _SPEC_TABLE)

# Convenience exports mirroring the original C++ header
kaya1 = EXAMPLES_BY_KEY["kaya1"].points
kaya2 = EXAMPLES_BY_KEY["kaya2"].points
kaya3 = EXAMPLES_BY_KEY["kaya3"].points
kaya4 = EXAMPLES_BY_KEY["kaya4"].points
# omega = EXAMPLES_BY_KEY["omega"].points
# circuit = EXAMPLES_BY_KEY["spa"].points


def iter_example_results(
    example_names: Optional[Iterable[str]] = None,
    discretizations: Sequence[int] = DEFAULT_DISCRETIZATIONS,
    refinements: Sequence[int] = DEFAULT_REFINEMENTS,
) -> Iterator[ExampleResult]:
    """Yield DP results for the requested examples.

    The implementation mirrors the logic in the original C++ `allexamples` helper,
    but returns structured data instead of printing LaTeX rows.
    """
    if example_names is None:
        names = EXAMPLE_ORDER
    else:
        names = tuple(example_names)
    for example_name in names:
        spec = EXAMPLES_BY_NAME[example_name]
        for discr in discretizations:
            for refin in refinements:
                yield _solve_example(spec, discr, refin)


def run_example(
    example_name: str,
    discretizations: Sequence[int] = DEFAULT_DISCRETIZATIONS,
    refinements: Sequence[int] = DEFAULT_REFINEMENTS,
) -> List[ExampleResult]:
    return list(iter_example_results([example_name], discretizations, refinements))


def allexamples(
    discretizations: Sequence[int] = DEFAULT_DISCRETIZATIONS,
    refinements: Sequence[int] = DEFAULT_REFINEMENTS,
) -> List[ExampleResult]:
    return list(iter_example_results(None, discretizations, refinements))


def _solve_example(spec: ExampleSpec, discretization: int, refinement: int) -> ExampleResult:
    points = list(spec.points)
    fixed_angles = list(spec.fixed_angles)
    def_thetas = list(spec.def_thetas)

    dp_instance = DP(
        points,
        fixed_angles,
        spec.k_max,
        discretization,
        refinement,
        def_thetas=def_thetas,
    )

    start = time.perf_counter()
    dp_instance.solve_dp()
    runtime_ms = (time.perf_counter() - start) * 1000.0

    angles, path_length = _extract_solution(dp_instance, len(points))
    recomputed_length = _compute_path_length(points, angles, spec.k_max)
    if abs(recomputed_length - path_length) > 1e-9:
        path_length = recomputed_length

    diff = path_length - spec.example_length
    return ExampleResult(
        example=spec.name,
        discretization=discretization,
        refinement=refinement,
        path_length=path_length,
        expected_length=spec.example_length,
        diff=diff,
        runtime_ms=runtime_ms,
        angles=angles,
    )


def _extract_solution(dp_instance: DP, expected_length: int) -> Tuple[List[float], float]:
    final_row = dp_instance.dp_matrix[-1]
    finite_cells = [cell for cell in final_row if math.isfinite(cell.l())]
    if not finite_cells:
        raise ValueError('No finite solution found in DP matrix')
    best_cell = min(finite_cells, key=lambda cell: cell.l())
    chain = []
    current = best_cell
    while current is not None:
        chain.append(current)
        current = current.prev()
    chain.reverse()
    if len(chain) != expected_length:
        raise ValueError(
            f'Expected {expected_length} states in optimal chain, got {len(chain)}'
        )
    angles = [float(cell.th()) for cell in chain]
    length = float(chain[-1].l())
    return angles, length


def _compute_path_length(
    points: Sequence[Tuple[float, float]],
    angles: Sequence[float],
    k_max: float,
) -> float:
    total = 0.0
    for idx in range(len(points) - 1):
        x0, y0 = points[idx]
        x1, y1 = points[idx + 1]
        _, _, segments = dubins_shortest_path(
            x0,
            y0,
            angles[idx],
            x1,
            y1,
            angles[idx + 1],
            k_max,
        )
        total += float(sum(segments))
    if not math.isfinite(total):
        raise ValueError("Computed path length is not finite")
    return total


def _format_row(result: ExampleResult) -> str:
    return (
        f"{result.example:>16} | {result.discretization:4d} | {result.refinement:3d} | "
        f"{result.path_length: .6f} | {result.expected_length: .6f} | "
        f"{result.diff_mm: .3f} mm | {result.runtime_ms:6.2f} ms"
    )


def print_table(
    results: Iterable[ExampleResult],
) -> None:
    print(f"{'Example':>16} | DISCR | REF |   Length | Reference | Delta (mm) |   Time")
    print("-" * 86)
    for row in results:
        print(_format_row(row))


if __name__ == "__main__":
    print_table(allexamples())
