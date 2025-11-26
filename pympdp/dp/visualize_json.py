"""Build the DP visualization HTML starting from a JSON dump produced by the C++ solver."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, List, Sequence

from pympdp.dp._viz_mixin import _VizMixin


def _coerce_float(value):
    if value is None:
        return None
    return float(value)


class _JsonCell:
    __slots__ = ("_angle", "_length", "_next", "_next_index")

    def __init__(self, angle, length, next_index):
        self._angle = _coerce_float(angle)
        self._length = _coerce_float(length)
        self._next = None
        self._next_index = None if next_index is None else int(next_index)

    def th(self):
        return self._angle

    def l(self):
        return self._length

    def next(self):
        return self._next

    def link(self, next_row: Sequence["_JsonCell"]):
        if self._next_index is None:
            return
        if self._next_index < 0 or self._next_index >= len(next_row):
            return
        self._next = next_row[self._next_index]


class _JsonDP(_VizMixin):
    def __init__(self, *, points, matrix, best_angles, k_max, best_path):
        self.points = points
        self.dp_matrix = matrix
        self.best_path = best_path
        self.k_max = float(k_max) if k_max is not None else 0.0
        self._best_angles = best_angles or []

    def best_angles(self, _points):
        return list(self._best_angles)


def _build_matrix(rows_payload: Sequence[Sequence[dict]]) -> List[List[_JsonCell]]:
    matrix: List[List[_JsonCell]] = []
    for row in rows_payload:
        matrix_row: List[_JsonCell] = []
        for cell_payload in row:
            next_idx = cell_payload.get("next_col")
            if next_idx is not None:
                try:
                    next_idx = int(next_idx)
                except (TypeError, ValueError):
                    next_idx = None
            matrix_row.append(
                _JsonCell(
                    cell_payload.get("theta"),
                    cell_payload.get("length"),
                    next_idx,
                )
            )
        matrix.append(matrix_row)

    for row_idx, row in enumerate(matrix[:-1]):
        next_row = matrix[row_idx + 1]
        for cell in row:
            cell.link(next_row)

    return matrix


def _sanitize_points(points_payload: Iterable[Sequence[float]]) -> List[tuple[float, float]]:
    sanitized: List[tuple[float, float]] = []
    for idx, point in enumerate(points_payload):
        if len(point) < 2:
            raise ValueError(f"Point #{idx} does not contain two coordinates")
        sanitized.append((float(point[0]), float(point[1])))
    return sanitized


def _sanitize_best_path(best_path_payload: Iterable[Sequence[int]]):
    result = []
    for entry in best_path_payload or []:
        if len(entry) < 2:
            continue
        result.append((int(entry[0]), int(entry[1])))
    return result


def _load_dp_from_json(json_path: Path) -> _JsonDP:
    payload = json.loads(json_path.read_text())
    points = _sanitize_points(payload.get("points", []))
    matrix_payload = payload.get("matrix")
    if matrix_payload is None:
        raise ValueError("JSON payload does not contain a 'matrix' entry")
    matrix = _build_matrix(matrix_payload)
    if len(matrix) != len(points):
        raise ValueError("Matrix row count and points length must match")

    best_angles = [float(angle) for angle in payload.get("best_angles", []) if angle is not None]
    best_path = _sanitize_best_path(payload.get("best_path", []))
    k_max = payload.get("k_max")

    return _JsonDP(
        points=points,
        matrix=matrix,
        best_angles=best_angles,
        k_max=k_max,
        best_path=best_path,
    )


def visualize_from_json(
    json_path: Path | str,
    *,
    output_path: Path | str | None = None,
    open_in_browser: bool = True,
    show_optimal_path: bool = True,
    samples_per_segment: int = 80,
):
    """Render the DP visualization from data exported by the C++ solver."""

    json_path = Path(json_path)
    dp = _load_dp_from_json(json_path)
    target_html = Path(output_path) if output_path else json_path.with_suffix(".html")
    dp.visualize_dp_matrix(
        output_path=target_html,
        open_in_browser=open_in_browser,
        show_optimal_path=show_optimal_path,
        samples_per_segment=samples_per_segment,
    )
    return target_html


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Visualize a DP matrix exported from the C++ solver")
    parser.add_argument("json_path", type=Path, help="Path to the JSON file produced via DP::exportVisualizationData")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Where to store the HTML dashboard (defaults to <json_path>.html)",
    )
    parser.add_argument("--no-browser", action="store_true", help="Do not automatically open the HTML file")
    parser.add_argument(
        "--hide-optimal-path",
        action="store_true",
        help="Skip sampling and drawing the optimal Dubins path on top of the matrix",
    )
    parser.add_argument(
        "--samples-per-segment",
        type=int,
        default=80,
        help="Sample density used when drawing Dubins segments (default: 80)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)
    visualize_from_json(
        args.json_path,
        output_path=args.output,
        open_in_browser=not args.no_browser,
        show_optimal_path=not args.hide_optimal_path,
        samples_per_segment=args.samples_per_segment,
    )


if __name__ == "__main__":
    main()
