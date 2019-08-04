

from .geometries import *
from ..projections.crs import *


class SegmentPair:
    """
    A pair of segments.

    """

    def __init__(self,
        segment1: Segment,
        segment2: Segment
    ):
        """

        :param segment1:
        :param segment2:
        """

        check_type(segment1, "First segment", Segment)
        check_type(segment2, "Second segment", Segment)

        check_crs(segment1, segment2)

        self._s1 = segment1
        self._s2 = segment2

    def same_start(self, tol: numbers.Real = 1e-12):
        """
        Check whether the two segments have the same start point.

        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

          Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(0,0,0), Point(0,1,0))
          >>>
        """

        return self._s1.start_pt().isCoinc(
            another=self._s2.start_pt(),
            tolerance=tol)

    def same_end(self, tol: numbers.Real = 1e-12):
        """
        Check whether the two segments have the same end point.

        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.
        """

        return self._s1.end_pt().isCoinc(
            another=self._s2.end_pt(),
            tolerance=tol)

    def first_conn_to_second(self, tol: numbers.Real = 1e-12):
        """
        Check whether the first segment is sequentially connected to the second one.

        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.
        """

        return self._s1.end_pt().isCoinc(
            another=self._s2.start_pt(),
            tolerance=tol)

    def second_conn_to_first(self, tol: numbers.Real = 1e-12):
        """
        Check whether the second segment is sequentially connected to the first one.

        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.
        """

        return self._s2.end_pt().isCoinc(
            another=self._s1.start_pt(),
            tolerance=tol)


class Segments(list):
    """
    Collection of segments, inheriting from list.

    """

    def __init__(self, segments: List[Segment]):

        check_type(segments, "Segments", List)
        for el in segments:
            check_type(el, "Segment", Segment)

        super(Segments, self).__init__(segments)


class Lines(list):
    """
    Collection of lines, inheriting from list.

    """

    def __init__(self, lines: List[Line]):

        check_type(lines, "Lines", List)
        for el in lines:
            check_type(el, "Line", Line)

        super(Lines, self).__init__(lines)


class MultiLines(list):
    """
    Collection of multilines, inheriting from list.

    """

    def __init__(self, multilines: List[MultiLine]):

        check_type(multilines, "MultiLines", List)
        for el in multilines:
            check_type(el, "MultiLine", MultiLine)

        super(MultiLines, self).__init__(multilines)

