
from math import sqrt
from numpy import isfinite

from pygsf.defaults.mathematics import MIN_VECTOR_MAGNITUDE
from pygsf.exceptions.vectors import InputValuesException


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