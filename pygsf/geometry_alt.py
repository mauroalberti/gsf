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

        self.a = val % (2*pi)

    def __repr__(self) -> str:

        return "Azimuth({:.2f})".format(degrees(self.a))

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

    @classmethod
    def fromXY(cls, x: [int, float], y: [int, float]) -> 'Azimuth':
        """
        Calculates azimuth given cartesian components.

        :param cls: class
        :param x: x component
        :param y: y component
        :return: Azimuth instance

        Examples:
          >>> Azimuth.fromXY(1, 1)
          Azimuth(45.00)
          >>> Azimuth.fromXY(1, -1)
          Azimuth(135.00)
          >>> Azimuth.fromXY(-1, -1)
          Azimuth(225.00)
          >>> Azimuth.fromXY(-1, 1)
          Azimuth(315.00)
          >>> Azimuth.fromXY(0, 0)
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

        return cls(angle, unit="rad")


class Plunge(object):
    """
    Class representing a plunge
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
        elif unit == "dd" and not (-90.0 <= val <= 90.0):
            raise GeomInputException("Input value in degrees must be between -90° and 90°")
        elif unit == "rad" and not (-pi/2 <= val <= pi/2):
            raise GeomInputException("Input value in radians must be between -pi/2 and pi/2")

        if unit == "dd":
            val = radians(val)

        self.a = val

    def __repr__(self) -> str:

        return "Plunge({:.2f})".format(degrees(self.a))

    def toHZ(self):

        """
        Converts an azimuth to h-z components.

        :return: a tuple storing h (horizontal) and z values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Plunge(0).toHZ())
          (1.0, 0.0)
          >>> apprFTuple(Plunge(90).toHZ())
          (0.0, -1.0)
          >>> apprFTuple(Plunge(-90).toHZ())
          (0.0, 1.0)
          >>> apprFTuple(Plunge(-45).toHZ(), ndec=6)
          (0.707107, 0.707107)
          >>> apprFTuple(Plunge(45).toHZ(), ndec=6)
          (0.707107, -0.707107)
        """

        return cos(self.a), -sin(self.a)

    @classmethod
    def fromHZ(cls, h: [int, float], z: [int, float]) -> 'Plunge':
        """
        Calculates plunge from h and z components.

        :param cls: class
        :param h: horizontal component (always positive)
        :param z: vertical component (positive upward)
        :return: Plunge instance

        Examples:
          >>> Plunge.fromHZ(1, 1)
          Plunge(-45.00)
          >>> Plunge.fromHZ(1, -1)
          Plunge(45.00)
          >>> Plunge.fromHZ(0, 1)
          Plunge(-90.00)
          >>> Plunge.fromHZ(0, -1)
          Plunge(90.00)
          >>> Plunge.fromHZ(-1, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Horizontal component cannot be negative
          >>> Plunge.fromHZ(0, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Input values cannot be both zero
        """

        # input vals check
        vals = [h, z]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            raise GeomInputException("Input values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise GeomInputException("Input values must be finite")
        elif h == 0.0 and z == 0.0:
            raise GeomInputException("Input values cannot be both zero")
        elif h < 0.0:
            raise GeomInputException("Horizontal component cannot be negative")

        angle = atan2(-z, h)

        return cls(angle, unit="rad")

"""
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