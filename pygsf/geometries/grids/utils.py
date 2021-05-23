import numbers
from typing import Tuple


def ij_array_to_ij_pixel(
        i_arr: numbers.Real,
        j_arr: numbers.Real
) -> Tuple[numbers.Real, numbers.Real]:
    """
    Converts from array indices to geotransform-related pixel indices.

    :param i_arr: the array i value.
    :type i_arr: numbers.Real.
    :param j_arr: the array j value.
    :type j_arr: numbers.Real.
    :return: the geotransform-equivalent i and j indices.
    :rtype: a tuple of two numbers.

    Examples:
      >>> ij_array_to_ij_pixel(0, 0)
      (0.5, 0.5)
      >>> ij_array_to_ij_pixel(0.5, 0.5)
      (1.0, 1.0)
      >>> ij_array_to_ij_pixel(1.5, 0.5)
      (2.0, 1.0)
    """

    return i_arr + 0.5, j_arr + 0.5


def ij_pixel_to_ij_array(
        i_pix: numbers.Real,
        j_pix: numbers.Real
) -> Tuple[numbers.Real, numbers.Real]:
    """
    Converts from pixel (geotransform-derived) to array indices.

    :param i_pix: the geotransform i value.
    :type i_pix: numbers.Real.
    :param j_pix: the geotransform j value.
    :type j_pix: numbers.Real.
    :return: the array-equivalent i and j indices.
    :rtype: a tuple of two numbers.

    Examples:
      >>> ij_pixel_to_ij_array(0, 0)
      (-0.5, -0.5)
      >>> ij_pixel_to_ij_array(0.5, 0.5)
      (0.0, 0.0)
      >>> ij_pixel_to_ij_array(0.5, 1.5)
      (0.0, 1.0)
    """

    return i_pix - 0.5, j_pix - 0.5