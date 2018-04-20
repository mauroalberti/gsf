# -*- coding: utf-8 -*-

from math import isnan, isinf, sqrt
import numpy as np

from typing import Dict, Tuple, List

from .default_parameters import *


isfinite = np.isfinite
array = np.array


def is_number(s: str) -> bool:
    """
    Check if string can be converted to number.

    @param  s:  parameter to check.
    @type  s:  string

    @return:  boolean, whether string can be converted to a number (float).

    Example:
      >>> is_number("1.0")
      True
      >>> is_number("1")
      True
      >>> is_number(u"-10")
      True
      >>> is_number("one")
      False
      >>> is_number("1e-10")
      True
      >>> is_number("")
      False
    """

    try:
        float(s)
    except:
        return False
    else:
        return True


def almost_zero(an_val: float, tolerance: float = 1e-10) -> bool:
    """
    Check if a value for which abs can be used, is near zero.

    :param an_val: an abs-compatible object
    :param tolerance: the tolerance value
    :return: Boolean

      >>> almost_zero(1)
      False
      >>> almost_zero(1e-9)
      False
      >>> almost_zero(1e-11)
      True
    """

    return abs(an_val) <= tolerance


def are_close(a: float, b: float, rtol: float = 1e-012, atol: float = 1e-12, equal_nan: bool = False, equal_inf: bool = False) -> bool:
    """
    Mimics math.isclose from Python 3.5 (see: https://docs.python.org/3.5/library/math.html)

    Example:
      >>> are_close(1.0, 1.0)
      True
      >>> are_close(1.0, 1.000000000000001)
      True
      >>> are_close(1.0, 1.0000000001)
      False
      >>> are_close(0.0, 0.0)
      True
      >>> are_close(0.0, 0.000000000000001)
      True
      >>> are_close(0.0, 0.0000000001)
      False
      >>> are_close(100000.0, 100000.0)
      True
      >>> are_close(100000.0, 100000.0000000001)
      True
      >>> are_close(float('nan'), float('nan'))
      False
      >>> are_close(float('nan'), 1000000)
      False
      >>> are_close(1.000000000001e300, 1.0e300)
      False
      >>> are_close(1.0000000000001e300, 1.0e300)
      True
      >>> are_close(float('nan'), float('nan'), equal_nan=True)
      True
      >>> are_close(float('inf'), float('inf'))
      False
      >>> are_close(float('inf'), 1.0e300)
      False
      >>> are_close(float('inf'), float('inf'), equal_inf=True)
      True
    """

    # nan cases
    if equal_nan and isnan(a) and isnan(b):
        return True
    elif isnan(a) or isnan(b):
        return False

    # inf cases
    if equal_inf and isinf(a) and a > 0 and isinf(b) and b > 0:
        return True
    elif equal_inf and isinf(a) and a < 0 and isinf(b) and b < 0:
        return True
    elif isinf(a) or isinf(b):
        return False

    # regular case
    return abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)


def apprFloat(val: [int, float], ndec: int=1):
    """
    Rounds a numeric value to ndec.

    :param val: value to round
    :param ndec: number of decimals used
    :return: rounded float value

    Examples:
      >>> apprFloat(0.00001)
      0.0
      >>> apprFloat(1.425324e-7)
      0.0
    """

    rval = round(val, ndec)
    if rval == 0.0:
        rval = round(0.0, ndec)

    return rval


def apprFTuple(tup: tuple, ndec=1):
    """
    Rounds numeric values inside a tuple to ndec decimals

    :param tup: tuple of int/float/exp values
    :param ndec: number of decimals used
    :return: tuple with rounded numbers

    Examples:
      >>> apprFTuple((-2.4492935982947064e-16, 1.0))
      (0.0, 1.0)
      >>> apprFTuple(((-1.0, -1.8369701987210297e-16)))
      (-1.0, 0.0)
    """

    return tuple(map(lambda val: apprFloat(val, ndec), tup))


def normXYZ(x: [int, float], y: [int, float], z: [int, float]):
    """
    Normalize numeric values.

    :param x: x numeric value
    :param y: y numeric value
    :param z: z numeric value
    :return: a tuple of three float values
    """

    # input vals checks
    vals = [x, y, z]
    if not all(map(lambda val: isinstance(val, (int, float)), vals)):
        raise InputValuesException("Input values must be integer or float")
    elif not all(map(isfinite, vals)):
        raise InputValuesException("Input values must be finite")

    mag = sqrt(x*x + y*y + z*z)

    if mag <= MIN_VECTOR_MAGNITUDE:
        norm_xyz = None
    else:
        norm_xyz = x/mag, y/mag, z/mag

    return mag, norm_xyz


class InputValuesException(Exception):
    """
    Exception for values input.
    """

    pass


if __name__ == '__main__':

    import doctest
    doctest.testmod()
