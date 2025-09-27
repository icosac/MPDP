import numpy as np

from pympdp.dp import DP
from pympdp.dp.cell import Cell
from pympdp.logger import logger
# logger.set_critical()

"""
Create suite of tests for dp.py using pytest
"""

def test_initialization():
    points = [(0, 0), (1, 1), (2, 0)]
    fixed_angles = [False, False, False]
    k_max = 1.0
    discretizations = 8
    refinements = 2

    dp_instance = DP(points, fixed_angles, k_max, discretizations, refinements)

    assert dp_instance.points == points
    assert dp_instance.fixed_angles == fixed_angles
    assert dp_instance.k_max == k_max
    assert dp_instance.discretizations == discretizations
    assert dp_instance.refinements == refinements
    assert len(dp_instance.dp_matrix) == len(points)
    assert all(len(row) == discretizations for row in dp_instance.dp_matrix)
    assert dp_instance.best_path is None

def test_reset_matrix():
    points = [(0, 0), (1, 1), (2, 0)]
    fixed_angles = [False, False, False]
    k_max = 1.0
    discretizations = 8
    refinements = 2

    dp_instance = DP(points, fixed_angles, k_max, discretizations, refinements)

    # Modify the dp_matrix to have non-default values
    for i in range(len(points)):
        for j in range(discretizations):
            dp_instance.dp_matrix[i][j] = Cell(_angle=1.0, _length=1.0, _next=None)

    dp_instance._reset_matrix()

    for i in range(len(points)):
        for j in range(discretizations):
            assert np.isnan(dp_instance.dp_matrix[i][j].th())
            assert dp_instance.dp_matrix[i][j].l() == 0
            assert dp_instance.dp_matrix[i][j].prev() is None


def test_set_sampling_angles1():
    points = [(0, 0), (1, 1), (2, 0)]
    fixed_angles = [False, True, False]
    k_max = 1.0
    discretizations = 4
    refinements = 0
    def_thetas = [0.0, 0.0, 0.0]

    dp_instance = DP(points, fixed_angles, k_max, discretizations, refinements, def_thetas=def_thetas)
    # dp_instance.set_sampling_angles(hrange=2*np.pi)

    # assert len(dp_instance.dp_matrix) == len(points), f"The matrix should have {len(points)} rows, but has {len(dp_instance.dp_matrix)}"
    # for i in range(len(points)):
    #     if fixed_angles[i]:
    #         assert len(dp_instance.dp_matrix[i]) == 1, f"The matrix row {i} should have 1 angle, but has {len(dp_instance.dp_matrix[i])}"
    #         assert dp_instance.dp_matrix[i][0].th() == def_thetas[i] , f"The angle for point {i} should be {def_thetas[i]}, but is {dp_instance.dp_matrix[i][0].th()}"
    #     else:
    #         assert len(dp_instance.dp_matrix[i]) >= discretizations, f"Unexpected number of angles {len(dp_instance.dp_matrix[i])} for point {i}"
    #         angles = np.linspace(def_thetas[i]-np.pi, def_thetas[i]+np.pi, discretizations, endpoint=False)
    #         for j in range(discretizations):
    #             assert np.isclose(dp_instance.dp_matrix[i][j].th(), angles[j]) , f"Angle {j} for point {i} should be {angles[j]}, but is {dp_instance.dp_matrix[i][j].th()}"

    
# def test_set_sampling_angles2():
#     points = [(0, 0), (1, 1), (2, 0)]
#     fixed_angles = [True, False, False]
#     k_max = 1.0
#     discretizations = 4
#     refinements = 0
#     def_thetas = [0, -np.pi/2, np.pi/2]

#     dp_instance = DP(points, fixed_angles, k_max, discretizations, refinements)
#     dp_instance.set_sampling_angles(def_thetas, hrange=2*np.pi)

#     for i in range(len(points)):
#         if i == 0:
#             assert len(dp_instance.dp_matrix[i]) == 1, f"Unexpected number of angles {len(dp_instance.dp_matrix[i])} for point {i}, should be 1"
#             assert dp_instance.dp_matrix[i][0].th() == def_thetas[i], f"The angle for point {i} should be {def_thetas[i]}, but is {dp_instance.dp_matrix[i][0].th()}"
#         if i==1:
#             comp_angles = [cell.th() for cell in dp_instance.dp_matrix[i]]
#             assert comp_angles == [-4.71238898038469, -3.141592653589793, -1.5707963267948966, 0.0, 0.7853981633974483, -0.7853981633974483]
#         if i==2:
#             comp_angles = [cell.th() for cell in dp_instance.dp_matrix[i]]
#             assert comp_angles == [-1.5707963267948966, 0.0, 1.5707963267948966, 3.141592653589793, -0.7853981633974483]



