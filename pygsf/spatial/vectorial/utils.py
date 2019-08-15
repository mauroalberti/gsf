# -*- coding: utf-8 -*-


from .geometries import *

from ..projections.crs import *

from ...utils.types import *
from ...mathematics.defaults import *


def point_or_segment(
        point1: Point,
        point2: Point,
        tol: numbers.Real = MIN_POINT_POS_DIFF
) -> Union[Point, Segment]:
    """
    Creates a point or segment based on the points distance.

    :param point1: first input point.
    :type point1: Point.
    :param point2: second input point.
    :type point2: Point.
    :param tol: distance tolerance between the two points.
    :type tol: numbers.Real.
    :return: point or segment based on their distance.
    :rtype: PointOrSegment.
    :raise: Exception.
    """

    check_type(point1, "First point", Point)
    check_type(point2, "Second point", Point)

    check_crs(point1, point2)

    if point1.dist3DWith(point2) <= tol:
        return point1.clone()
    else:
        return Segment(
            start_pt=point1,
            end_pt=point2
        )


# TODO class Path
"""
The trajectory of a Point with time
"""


if __name__ == "__main__":

    import doctest
    doctest.testmod()
