# -*- coding: utf-8 -*-


from __future__ import division
#from __future__ import annotations

from math import sqrt, sin, cos, radians, acos, atan, isnan, pi
import numpy as np

from typing import Dict, Tuple, List

from .default_parameters import *
from .mathematics import *
from .geometry_utils import *
from .arrays import point_solution, arrays_are_close, arr2tuple


isfinite = np.isfinite
array = np.array


class Azimuth(object):
    """
    Azim class
    """

    def __init__(self, val: [int, float], unit: str="d"):

        # unit check
        if unit not in ("d", "r"):
            raise GeomInputException("Unit input must be d or r")

        # val check
        if not isinstance(val, (int, float)):
            raise GeomInputException("Input value must be integer or float")
        elif not isfinite(val):
            raise GeomInputException("Input value must be finite")

        if unit == "d":
            val = radians(val)

        self.a = val % (2*pi)

    def __repr__(self) -> str:

        return "Azim({:.2f})".format(degrees(self.a))

    def toXY(self):
        """
        Converts an azimuth to x-y components.

        :return: a tuple storing x and y values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Azim(0).toXY())
          (0.0, 1.0)
          >>> apprFTuple(Azim(90).toXY())
          (1.0, 0.0)
          >>> apprFTuple(Azim(180).toXY())
          (0.0, -1.0)
          >>> apprFTuple(Azim(270).toXY())
          (-1.0, 0.0)
          >>> apprFTuple(Azim(360).toXY())
          (0.0, 1.0)
        """

        return sin(self.a), cos(self.a)

    @classmethod
    def fromXY(cls, x: [int, float], y: [int, float]) -> 'Azim':
        """
        Calculates azimuth given cartesian components.

        :param cls: class
        :param x: x component
        :param y: y component
        :return: Azim instance

        Examples:
          >>> Azim.fromXY(1, 1)
          Azim(45.00)
          >>> Azim.fromXY(1, -1)
          Azim(135.00)
          >>> Azim.fromXY(-1, -1)
          Azim(225.00)
          >>> Azim.fromXY(-1, 1)
          Azim(315.00)
          >>> Azim.fromXY(0, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Input values cannot be both zero

        """

        # input vals check
        vals = [x, y]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            raise GeomInputException("Input values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise GeomInputException("Input values must be finite")
        elif x == 0.0 and y == 0.0:
            raise GeomInputException("Input values cannot be both zero")

        angle = atan2(x, y)

        return cls(angle, unit="r")


class Plunge(object):
    """
    Class representing a plunge
    """

    def __init__(self, val: [int, float], 