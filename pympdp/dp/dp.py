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
    from pympdp.profiler import profile_callable
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
    from pympdp.profiler import profile_callable


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
        opt_len = min(cell.l() for cell in self.dp_matrix[0] if np.isfinite(cell.l()))
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

            tmp_len = min(cell.l() for cell in self.dp_matrix[0] if np.isfinite(cell.l()))
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
        n_points = len(self.points)
        for idx in range(n_points - 2, -1, -1):
            for cell_i in self.dp_matrix[idx]:
                cell_i._length = np.inf
                cell_i._next = None

            for i, cell_i in enumerate(self.dp_matrix[idx]):
                for j, cell_j in enumerate(self.dp_matrix[idx + 1]):
                    _, _, lengths = dubins_shortest_path(
                        self.points[idx][0], self.points[idx][1], cell_i.th(),
                        self.points[idx + 1][0], self.points[idx + 1][1], cell_j.th(),
                        self.k_max,
                    )

                    if not lengths or any(np.isnan(lengths)) or not np.isfinite(cell_j.l()):
                        continue

                    curr_length = sum(lengths) + cell_j.l()

                    if curr_length < cell_i.l():
                        cell_i._length = curr_length
                        cell_i._next = cell_j

                # logger.debug(f"Best from point {idx} angle {i} is to point {idx + 1} angle {j} with length {cell_i.l()}")
        logger.debug("Finished DP")


    def best_angles(self, points):
        if not self.dp_matrix or not self.dp_matrix[0]:
            return []

        first_row = self.dp_matrix[0]
        finite_cells = [cell for cell in first_row if np.isfinite(cell.l())]
        current = min(finite_cells or first_row, key=lambda cell: cell.l())

        best_angles = []
        visited_ids = set()

        while current is not None and len(best_angles) < len(points):
            if id(current) in visited_ids:
                logger.warning("Cycle detected while extracting the optimal chain; aborting traversal")
                break
            visited_ids.add(id(current))
            best_angles.append(current.th())
            current = current.next()

        return best_angles
    

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
                id_next = None
                if i < len(self.dp_matrix) - 1 and cell.next() is not None:
                    for j, next_cell in enumerate(self.dp_matrix[i + 1]):
                        if cell.next() is not None and np.isclose(cell.next().th(), next_cell.th()):
                            id_next = f"{j:4d}"
                            break

                angles_str.append(f"{angle_str}, {length_str}, {id_next}")
            
            # Fill in empty cells if the row has fewer angles than max_n_ths
            while len(angles_str) < max_n_ths:
                angles_str.append("None, None, None")
            table.add_row([f"#{i} {point_str}"] + angles_str)

            # next_angle = cell.next().th() if cell.next() else None
            # table.add_row([i, f"{cell.th():.4f}" if cell.th() is not None else "None", f"{cell.l():.4f}" if cell.l() is not None else "None", f"{next_angle:.4f}" if next_angle is not None else "None"])
        print(table)


import argparse
def main():
    parser = argparse.ArgumentParser(description="Test DP class")
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    parser.add_argument('--profile', action='store_true', help='Run solve_dp under cProfile and print top entries')
    parser.add_argument('--profile-limit', type=int, default=30, help='Number of profile rows to display (default: 30)')
    parser.add_argument('--profile-sort', default='cumulative', help='Sort key for profile table (default: cumulative)')
    parser.add_argument('--profile-include-external', action='store_true', help='Include non-pympdp call sites in the profile output')
    parser.add_argument('--profile-keep-dirs', action='store_true', help='Keep full directory paths in profile output')
    args = parser.parse_args()

    if args.debug:
        logger.set_debug()

    # points = [(0, 0), (1, 1), (2, 1), (3, 0)]
    # def_thetas = [-np.pi, 0.0, 0.0, np.pi]

    points = [
        [0.010000025149524602, 0.0005770249507223292],
        [0.8289428479444483, -0.5732980947335886], 
        [1.3932863331700118, -1.3988381789714552], 
        [2.281286073185757, -1.8586821326470127], 
        [3.1792115718072362, -2.298829606711523], 
        [3.7196043491453383, -3.1402424935874094], 
        [4.3410319309946, -4.038320625487824], 
        [5.244159851813759, -3.608949039752501], 
        [6.23389945140574, -3.7518319463561465], 
        [6.687767640356647, -2.860763112281907], 
        [7.141635829307553, -1.9696942782076672], 
        [7.59550401825846, -1.0786254441334275], 
        [8.00309570257371, -0.021513553227158866], 
        [8.129899081617857, 0.9704143186783474], 
        [8.689678934469265, 1.0695426597270195]
    ]

    def_thetas = [-1.3413855108602083e-05] + [0.0] * (len(points)-2) + [0.5235987755982991]
    fixed_angles = [True] + [False] * (len(points)-2) + [True]

    dp_instance = DP(points, fixed_angles=fixed_angles, k_max=3, discretizations=30, refinements=2, def_thetas=def_thetas)
    start = time.time()
    if args.profile:
        angles, profile_table, _ = profile_callable(
            dp_instance.solve_dp,
            sort_by=args.profile_sort,
            limit=args.profile_limit,
            include_external=args.profile_include_external,
            strip_dirs=not args.profile_keep_dirs,
        )
        elapsed = time.time() - start
        logger.info(f"Solved in {elapsed:.4f} seconds (profiled)")
        print(profile_table)
    else:
        angles = dp_instance.solve_dp()
        elapsed = time.time() - start
        logger.info(f"Solved in {elapsed:.4f} seconds")
    # dp_instance.print_dp_matrix()
    dp_instance.visualize_dp_matrix(show_optimal_path=True)


    # from dubins import plotdubins
    # dub1, _, _ = dubins_shortest_path(0, 0, -np.pi, 1, 1, 0.0, 2)
    # dub2, _, _ = dubins_shortest_path(1, 1, 0.0, 2, 0, np.pi, 2)
    # plotdubins(dub1)
    # plotdubins(dub2, color1='m', color2='c', color3='y', show=True)


if __name__ == "__main__":
    main()