# -*- coding: utf-8 -*-

import numpy as np


def is_number(s):
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
      >>> is_number(None)
      False
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


def almost_zero(an_val, tolerance=1e-10):
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


def isclose(a, b, rtol=1e-012, atol=1e-12, equal_nan=False, equal_inf=False):
    """
    Mimics math.isclose from Python 3.5 (see: https://docs.python.org/3.5/library/math.html)

    Example:
      >>> isclose(1.0, 1.0)
      True
      >>> isclose(1.0, 1.000000000000001)
      True
      >>> isclose(1.0, 1.0000000001)
      False
      >>> isclose(0.0, 0.0)
      True
      >>> isclose(0.0, 0.000000000000001)
      True
      >>> isclose(0.0, 0.0000000001)
      False
      >>> isclose(100000.0, 100000.0)
      True
      >>> isclose(100000.0, 100000.0000000001)
      True
      >>> isclose(np.nan, np.nan)
      False
      >>> isclose(np.nan, 1000000)
      False
      >>> isclose(1.000000000001e300, 1.0e300)
      False
      >>> isclose(1.0000000000001e300, 1.0e300)
      True
      >>> isclose(np.nan, np.nan, equal_nan=True)
      True
      >>> isclose(np.inf, np.inf)
      False
      >>> isclose(np.inf, 1.0e300)
      False
      >>> isclose(np.inf, np.inf, equal_inf=True)
      True
    """

    if equal_nan and a is np.nan and b is np.nan:
        return True
    elif equal_inf and a is np.inf and b is np.inf:
        return True
    elif a is np.inf or b is np.inf:
        return False
    else:
        return abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)



if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()
