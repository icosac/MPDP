from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

try:
    from pympdp.logger import logger
    from pympdp.utility import circles
    from pympdp.dubins import dubins_shortest_path
    from pympdp.dp.cell import Cell
    from pympdp.dp._viz_mixin import _VizMixin
except ModuleNotFoundError as exc:  # pragma: no cover - developer convenience
    if exc.name not in {
        "pympdp",
        "pympdp.logger",
        "pympdp.utility",
        "pympdp.dubins",
        "pympdp.dp",
        "pympdp.dp.cell",
        "pympdp.dp._viz_mixin",
    }:
        raise
    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from pympdp.logger import logger
    from pympdp.utility import circles
    from pympdp.dubins import dubins_shortest_path
    from pympdp.dp.cell import Cell
    from pympdp.dp._viz_mixin import _VizMixin


class DP(_VizMixin):
    def __init__(self, points, fixed_angles, k_max, discretizations, refinements, def_thetas=None):
        if def_thetas:
            for i in range(len(def_thetas)):
                if fixed_angles[i] and def_thetas[i] is None:
                    raise AssertionError("if def_thetas is provided, all angles that are fixed must have a value in def_thetas")

        assert len(points) == len(fixed_angles), "Length of points and fixed_angles must be the same"
        assert def_thetas is None or len(def_thetas) == len(points), "Length of def_thetas must be the same as points if provided"
        assert def_thetas is None or all((not fixed) or (theta is not None) for fixed, theta in zip(fixed_angles, def_thetas)), "if def_thetas is provided, all angles that are fixed must have a value in def_thetas"
        assert def_thetas is not None or (def_thetas is None and all(not fixed_angles[i] for i in range(len(fixed_angles)))), "if def_thetas is None, all fixed_angles must be False"

        self.points = points
        self.fixed_angles = fixed_angles
        self.k_max = k_max
        self.discretizations = discretizations
        self.refinements = refinements

        self.def_thetas = def_thetas if def_thetas is not None else [0.0 for _ in points]

        self.dp_matrix = [[Cell() for _ in range(discretizations)] for _ in range(len(points))]
        self.best_path = None

        self.set_sampling_angles(hrange=2*np.pi)


    def set_first_angle(self, angle):
        """
        Set the angle of the first point in the DP matrix.
        :param angle: Angle in radians
        """
        if self.def_thetas is None:
            self.def_thetas = [0.0 for _ in self.points]
        
        self.def_thetas[0] = angle
        self.fixed_angles[0] = True
        self.dp_matrix[0] = [Cell(angle, 0, None)]

    
    def set_final_angle(self, angle):
        """
        Set the angle of the final point in the DP matrix.
        :param angle: Angle in radians
        """
        if self.def_thetas is None:
            self.def_thetas = [0.0 for _ in self.points]
        
        self.def_thetas[-1] = angle
        self.fixed_angles[-1] = True
        self.dp_matrix[-1] = [Cell(angle, 0, None)]


    def _reset_matrix(self):
        """
        Reset the dynamic programming matrix to an empty state.
        """
        for i in range(len(self.dp_matrix)):
            for j in range(len(self.dp_matrix[i])):
                self.dp_matrix[i][j] = Cell()


    def set_sampling_angles(self, hrange=2*np.pi):
        """
        Set the sampling angles for dynamic programming.
        :param hrange: Range of angles to sample. Default is 2pi for the first round, 3/2pi for subsequent rounds.
        """
        assert self.def_thetas is not None
        sampling_angles = [np.linspace(th-hrange/2, th+hrange/2, self.discretizations, endpoint=(hrange != 2*np.pi)).tolist() for th in self.def_thetas]

        for i in range(len(self.points)):
            if self.fixed_angles[i]:
                self.dp_matrix[i] = [Cell(self.def_thetas[i], 0, None)]
            else:
                self.dp_matrix[i] = [Cell(th, 0, None) for th in sampling_angles[i]]
                if i > 0:
                    th_prev, th_curr = self.guess_initial_angles(i)

                    # Adding guessed angles to the previous point if they are not already present
                    if not self.fixed_angles[i - 1]:
                        tmp_prev_angles = [cell.th() for cell in self.dp_matrix[i - 1]]
                        for th in th_prev:
                            if not any(np.isclose(th.th(), angle) for angle in tmp_prev_angles):
                                self.dp_matrix[i - 1].append(th)

                    # Adding guessed angles to the current point if they are not already present
                    tmp_curr_angles = [cell.th() for cell in self.dp_matrix[i]]
                    for th in th_curr:
                        if not any(np.isclose(th.th(), angle) for angle in tmp_curr_angles):
                            self.dp_matrix[i].append(th)


    def guess_initial_angles(self, i):
        """
        Guess initial angles for dynamic programming based on previous and current angles.
        :param i: Index of the current point
        """
        if i == 0:
            return [], []
        th = np.arctan2(self.points[i][1] - self.points[i - 1][1], self.points[i][0] - self.points[i - 1][0])
        th_prev = [Cell(th)]
        th_curr = [Cell(th)]

        XC, YX = circles(self.points[i - 1][0], self.points[i - 1][1], self.points[i][0], self.points[i][1], 1.0 / self.k_max)

        for xc, yc in zip(XC, YX):
            thp = np.arctan2(yc - self.points[i - 1][1], xc - self.points[i - 1][0])
            thc = np.arctan2(self.points[i][1] - yc, self.points[i][0] - xc)
            if np.isclose(thp, th_prev[-1].th()):
                th_prev.append(Cell(thp))
            if np.isclose(thc, th_curr[-1].th()):
                th_curr.append(Cell(thc))

        return th_prev, th_curr


    def solve_dp(self):
        # First round is on the house
        self.solve_dp_inner()

        # Extract optimal angles and best length
        opt_len = min(cell.l() for cell in self.dp_matrix[-1] if np.isfinite(cell.l()))
        opt_angles = self.best_angles(self.points)
        logger.debug(f"Initial optimal angles: {opt_angles}")
        logger.debug(f"Initial best length found: {opt_len}")
        self.def_thetas = opt_angles

        # Refinement rounds
        for r in range(self.refinements):
            logger.debug(f"Refinement round {r+1}/{self.refinements}")
            self._reset_matrix()
            self.set_sampling_angles(hrange = np.pi)
            
            self.solve_dp_inner()

            tmp_len = min(cell.l() for cell in self.dp_matrix[-1] if np.isfinite(cell.l()))
            opt_angles = self.best_angles(self.points)

            logger.debug(f"Optimal angles for refinement {r+1}: {opt_angles}")

            if tmp_len < opt_len:
                opt_len = tmp_len
                self.def_thetas = opt_angles
                logger.info(f"New best length found: {opt_len}")
            else:
                logger.info(f"No improvement in length: {tmp_len} >= {opt_len}")
                break

        logger.info(f"Optimal angles: {opt_angles}")

        return opt_angles


    def solve_dp_inner(self):
        for idx in range(0, len(self.points)-1):
            for j, cell_j in enumerate(self.dp_matrix[idx + 1]):
                cell_j._length = np.inf
                cell_j._next = None

                for i, cell_i in enumerate(self.dp_matrix[idx]):
                    # Compute Dubins path from (idx-1, cell_i) to (idx, cell_j)
                    _, _, lengths = dubins_shortest_path(
                        self.points[idx][0],     self.points[idx][1],     cell_i.th(),
                        self.points[idx + 1][0], self.points[idx + 1][1], cell_j.th(),
                        self.k_max
                    )

                    if not lengths or any(np.isnan(lengths)):
                        continue
                    
                    curr_length = sum(lengths) + (cell_i.l() if idx > 0 and not np.isinf(cell_i.l()) else 0)

                    if curr_length < cell_j.l():
                        cell_j._length = curr_length
                        cell_j._next = self.dp_matrix[idx][i]

                # logger.debug(f"Best from point {idx + 1} angle {j} is to point {idx} angle {i} with length {cell_j.l()}")
        logger.debug("Finished DP")


    def best_angles(self, points):
        # Backtrack to find the best angles
        best_angles = []
        for i in range(len(points)-1, -1, -1):
            if not self.dp_matrix[i]:
                continue
            best_cell = min(self.dp_matrix[i], key=lambda cell: cell.l())
            best_angles.append(best_cell.th())
            # Backtrack through the best previous cells
            while best_cell.prev():
                best_cell = best_cell.prev()
        return best_angles[::-1]
    

    def print_dp_matrix(self):
        import prettytable
        table = prettytable.PrettyTable()
        max_n_ths = max(len(row) for row in self.dp_matrix)
        table.field_names = ["Point Index"] + [f"Angle {i}" for i in range(max_n_ths)]
        for i, row in enumerate(self.dp_matrix):
            point_str = f"({self.points[i][0]:.2f}, {self.points[i][1]:.2f})"
            angles_str = []
            for cell in row:
                angle_str = f"{cell.th():.2f}" if cell.th() is not None else "None"
                length_str = f"{cell.l():.2e}" if cell.l() is not None else "None"
                id_prev = None
                if i > 0 and cell.prev() is not None:
                    for j, prev_cell in enumerate(self.dp_matrix[i - 1]):
                        if cell.prev() is not None and np.isclose(cell.prev().th(), prev_cell.th()):
                            id_prev = f"{j:4d}"
                            break

                angles_str.append(f"{angle_str}, {length_str}, {id_prev}")
            
            # Fill in empty cells if the row has fewer angles than max_n_ths
            while len(angles_str) < max_n_ths:
                angles_str.append("None, None, None")
            table.add_row([f"#{i} {point_str}"] + angles_str)

            # prev_angle = cell.prev().th() if cell.prev() else None
            # table.add_row([i, f"{cell.th():.4f}" if cell.th() is not None else "None", f"{cell.l():.4f}" if cell.l() is not None else "None", f"{prev_angle:.4f}" if prev_angle is not None else "None"])
        print(table)


import argparse
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Test DP class")
    args.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = args.parse_args()

    if args.debug:
        logger.set_debug()

    points = [(0, 0), (1, 1), (2, 1), (3, 0)]
    def_thetas = [-np.pi, 0.0, 0.0, np.pi]

    # points = [(0, 0), (1, 1), (2, 0)]

    dp_instance = DP(points, fixed_angles=[True, False, False, True], k_max=3, discretizations=90, refinements=1, def_thetas=def_thetas)
    now = time.time()
    dp_instance.solve_dp()
    logger.info(f"Solved in {time.time()-now:.4f} seconds")
    # dp_instance.print_dp_matrix()
    dp_instance.visualize_dp_matrix(show_optimal_path=True)


    # from dubins import plotdubins
    # dub1, _, _ = dubins_shortest_path(0, 0, -np.pi, 1, 1, 0.0, 2)
    # dub2, _, _ = dubins_shortest_path(1, 1, 0.0, 2, 0, np.pi, 2)
    # plotdubins(dub1)
    # plotdubins(dub2, color1='m', color2='c', color3='y', show=True)
