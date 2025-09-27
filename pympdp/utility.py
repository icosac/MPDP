import numpy as np 

def circles (x1, y1, x2, y2, r):
    """
    Calculate the centers of the circles of radius r passing through points (x1, y1) and (x2, y2).
    :param x1: x-coordinate of the first point
    :param y1: y-coordinate of the first point
    :param x2: x-coordinate of the second point
    :param y2: y-coordinate of the second point
    :param r: radius of the circles
    :return: Two lists containing the x and y coordinates of the circle centers
    """
    d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2
    if d2 > 4 * r * r:
        return [], []  # No solution

    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    q = np.sqrt(r * r - d2 / 4)
    dx = (y1 - y2) / np.sqrt(d2)
    dy = (x2 - x1) / np.sqrt(d2)

    xc1 = mid_x + q * dx
    yc1 = mid_y + q * dy
    xc2 = mid_x - q * dx
    yc2 = mid_y - q * dy

    return [xc1, xc2], [yc1, yc2]