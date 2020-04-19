
from typing import Tuple

import numbers
from affine import Affine


def gdal_to_affine(
    geotransform: Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Real]
) -> Affine:
    """
    Create an affine transformation
    from a GDAL geotransform tuple.

    """

    return Affine.from_gdal(*geotransform)


def forward_transformation(
    trans: Affine,
    row: numbers.Real,
    col: numbers.Real
) -> Tuple[numbers.Real, numbers.Real]:
    """
    Calculate the x, y coordinates given an affine transformation
    and the row, col values.

    """

    return trans * (col, row)


def backward_transformation(
    trans: Affine,
    x: numbers.Real,
    y: numbers.Real
) -> Tuple[numbers.Real, numbers.Real]:
    """
    Calculate the row, column values give an affine transformation
    and the x, y values.

    """

    rev = ~trans
    col, row = rev * (x, y)
    return row, col