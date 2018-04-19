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
    Azimuth class
    """

    def __init__(self, val: [int, float], unit: str="dd"):

        # unit check
        if unit not in ("dd", "rad"):
            raise GeomInputException("Unit input must be dd or rad")

        # val check
        if not isinstance(val, (int, float)):
            raise GeomInputException("Input value must be integer or float")
        elif not isfinite(val):
            raise GeomInputException("Input value must be finite")

        if unit == "dd":
            val = radians(val)

        self.a = val

    def toXY(self):
        """
        Converts an azimuth to x-y components.

        :return: a tuple storing x and y values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Azimuth(0).toXY())
          (0.0, 1.0)
          >>> apprFTuple(Azimuth(90).toXY())
          (1.0, 0.0)
          >>> apprFTuple(Azimuth(180).toXY())
          (0.0, -1.0)
          >>> apprFTuple(Azimuth(270).toXY())
          (-1.0, 0.0)
          >>> apprFTuple(Azimuth(360).toXY())
          (0.0, 1.0)
        """

        return sin(self.a), cos(self.a)

"""
    @staticmethod
    fromXY

class Plunge():

    unit
    val
    toHZ(self)

    @staticmethod
    fromHZ


class SpatOrient():

    azim
    plunge

    toXYZ

    @staticmethod
    fromXYZ


class Vect():

    spat_or
    magnitude

    toXYZ
    x
    y
    z

    @staticmethod
    fromXYZ


class Segment():

    vect
    start_pt
    
    @staticmethod
    fromPoints
    
    @staticmethod
    fromXYZ
    
    @staticmethod
    fromPointDirMagn
"""

class GeomInputException(Exception):
    """
    Exception for geometric input.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()