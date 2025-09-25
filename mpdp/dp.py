import numpy as np

from logger import logger
from cell import Cell
from utility import circles
from dubins import dubins_shortest_path

class DP:
    def __init__(self, points, fixed_angles, k_max, discretizations, refinements):
        self.points = points
        self.fixed_angles = fixed_angles
        self.k_max = k_max
        self.discretizations = discretizations
        self.refinements = refinements

        self.dp_matrix = [[Cell() for _ in range(discretizations)] for _ in range(len(points))]
        self.best_path = None

        self.set_sampling_angles([0.0 for _ in points], hrange=2*np.pi)


    def _reset_matrix(self):
        """
        Reset the dynamic programming matrix to an empty state.
        """
        for i in range(len(self.points)):
            for j in range(self.discretizations):
                self.dp_matrix[i][j] = Cell()


    def set_sampling_angles(self, def_thetas, hrange=2*np.pi):
        """
        Set the sampling angles for dynamic programming.
        :param hrange: Range of angles to sample. Default is 2pi for the first round, 3/2pi for subsequent rounds.
        """
        assert len(def_thetas) == len(self.points)
        sampling_angles = [np.linspace(th-hrange/2, th+hrange/2, self.discretizations, endpoint=(hrange != 2*np.pi)).tolist() for th in def_thetas]

        for i in range(len(self.points)):
            if self.fixed_angles[i]:
                self.dp_matrix[i] = [Cell(def_thetas[i], 0, None)]
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

        logger.info(f"curr {[cell.th() for cell in th_curr]}")
        logger.info(f"prev {[cell.th() for cell in th_prev]}")

        return th_prev, th_curr


    def solve_dp(self):
        # First round is on the house
        self.solve_dp_inner()

        # Extract optimal angles
        opt_angles = self.best_angles(self.points)

        # Refinement rounds
        for r in range(self.refinements):
            logger.info(f"Refinement round {r+1}/{self.refinements}")
            self._reset_matrix()
            self.set_sampling_angles(opt_angles, hrange=3*np.pi/2)
            self.solve_dp_inner()
            opt_angles = self.best_angles(self.points)

        return opt_angles


    def solve_dp_inner(self):
        idx = len(self.points) - 1
        while idx > 0:
            for i, cell_i in enumerate(self.dp_matrix[idx - 1]):
                best_length = np.inf
                best_prev = None

                for j, cell_j in enumerate(self.dp_matrix[idx]):
                    # Compute Dubins path from (idx-1, cell_i) to (idx, cell_j)
                    curve, i, lengths = dubins_shortest_path(
                        self.points[idx - 1][0], self.points[idx - 1][1], cell_i.th(),
                        self.points[idx][0], self.points[idx][1], cell_j.th(),
                        self.k_max
                    )

                    curr_length = sum(lengths) + (cell_j.l() if idx < len(self.points) - 1 else 0)

                    if curr_length < best_length:
                        best_length = curr_length
                        best_prev = self.dp_matrix[idx - 1][i]

                # Update the DP matrix with the best previous angle and length
                self.dp_matrix[idx-1][i] = Cell(cell_i.th(), best_length, best_prev)

            idx -= 1


    def best_angles(self, points):
        # Backtrack to find the best angles
        best_angles = []
        for i in range(len(points) - 1, -1, -1):
            if not self.dp_matrix[i]:
                continue
            best_cell = min(self.dp_matrix[i], key=lambda cell: cell.l())
            best_angles.append(best_cell.th())
            # Backtrack through the best previous cells
            while best_cell.prev():
                best_cell = best_cell.prev()
                best_angles.append(best_cell.th())
        return best_angles[::-1]



if __name__ == "__main__":
    points = [(0, 0), (1, 1), (2, 0)]
    fixed_angles = [False, True, False]
    k_max = 1.0
    discretizations = 4
    refinements = 0
    def_thetas = [0.0, 0, 0]

    dp_instance = DP(points, fixed_angles, k_max, discretizations, refinements)
    dp_instance.set_sampling_angles(def_thetas, hrange=2*np.pi)

    assert len(dp_instance.dp_matrix) == len(points), f"The matrix should have {len(points)} rows, but has {len(dp_instance.dp_matrix)}"
    for i in range(len(points)):
        if fixed_angles[i]:
            assert len(dp_instance.dp_matrix[i]) == 1, f"The matrix row {i} should have 1 angle, but has {len(dp_instance.dp_matrix[i])}"
            assert dp_instance.dp_matrix[i][0].th() == def_thetas[i] , f"The angle for point {i} should be {def_thetas[i]}, but is {dp_instance.dp_matrix[i][0].th()}"
        else:
            assert len(dp_instance.dp_matrix[i]) >= discretizations, f"Unexpected number of angles {len(dp_instance.dp_matrix[i])} for point {i}"
            angles = np.linspace(def_thetas[i]-np.pi, def_thetas[i]+np.pi, discretizations, endpoint=False)
            for j in range(discretizations):
                assert np.isclose(dp_instance.dp_matrix[i][j].th(), angles[j]) , f"Angle {j} for point {i} should be {angles[j]}, but is {dp_instance.dp_matrix[i][j].th()}"