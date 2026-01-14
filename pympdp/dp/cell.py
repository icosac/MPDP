import numpy as np


class Cell:
    def __init__(self, _angle=np.nan, _length=0.0, _next=None):
        self._angle = _angle
        self._length = _length
        self._next = _next

    def th(self):
        return self._angle

    def l(self):
        return self._length

    def next(self):
        return self._next
