

from math import fabs


from .geometries import *
from ..projections.crs import *


def shortest_segment_or_point(
    segment1: Segment,
    segment2: Segment
) -> Optional[Union[Segment, Point]]:

    """
    Calculates the optional shortest segment - or the intesection point - between two segments.

    Adapted from:
        http://paulbourke.net/geometry/pointlineplane/

    C code from:
        http://paulbourke.net/geometry/pointlineplane/lineline.c
[
    typedef struct {
    double x,y,z;
    } XYZ;

    /*
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. Calculate also the values of mua and mub where
      Pa = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - P3)
    Return FALSE if no solution exists.
    */
    int LineLineIntersect(
    XYZ p1,XYZ p2,XYZ p3,XYZ p4,XYZ *pa,XYZ *pb,
    double *mua, double *mub)
    {
    XYZ p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;

    p13.x = p1.x - p3.x;
    p13.y = p1.y - p3.y;
    p13.z = p1.z - p3.z;
    p43.x = p4.x - p3.x;
    p43.y = p4.y - p3.y;
    p43.z = p4.z - p3.z;
    if (ABS(p43.x) < EPS && ABS(p43.y) < EPS && ABS(p43.z) < EPS)
      return(FALSE);
    p21.x = p2.x - p1.x;
    p21.y = p2.y - p1.y;
    p21.z = p2.z - p1.z;
    if (ABS(p21.x) < EPS && ABS(p21.y) < EPS && ABS(p21.z) < EPS)
      return(FALSE);

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

    denom = d2121 * d4343 - d4321 * d4321;
    if (ABS(denom) < EPS)
      return(FALSE);
    numer = d1343 * d4321 - d1321 * d4343;

    *mua = numer / denom;
    *mub = (d1343 + d4321 * (*mua)) / d4343;

    pa->x = p1.x + *mua * p21.x;
    pa->y = p1.y + *mua * p21.y;
    pa->z = p1.z + *mua * p21.z;
    pb->x = p3.x + *mub * p43.x;
    pb->y = p3.y + *mub * p43.y;
    pb->z = p3.z + *mub * p43.z;

    return(TRUE);
    }

    :param segment1: the first segment.
    :type segment1: Segment.
    :param segment2: the second segment.
    :type segment2: Segment.
    :return: the optional shortest segment or an intersection point.
    :rtype: Optional[Union[Segment, Point]]
    """

    check_type(segment1, "First segment", Segment)
    check_type(segment2, "Second segment", Segment)

    check_crs(segment1, segment2)

    epsg_cd = segment1.epsg()

    p1 = segment1.start_pt
    p2 = segment1.end_pt

    p3 = segment2.start_pt
    p4 = segment2.end_pt

    p13 = Point(
        x=p1.x - p3.x,
        y=p1.y - p3.y,
        z=p1.z - p3.z,
        epsg_cd=epsg_cd
    )

    p43 = Point(
        x=p4.x - p3.x,
        y=p4.y - p3.y,
        z=p4.z - p3.z,
        epsg_cd=epsg_cd
    )

    if p43.asVect().isAlmostZero:
        return None

    p21 = Point(
        x=p2.x - p1.x,
        y=p2.y - p1.y,
        z=p2.z - p1.z,
        epsg_cd=epsg_cd
    )

    if p21.asVect().isAlmostZero:
        return None

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z

    denom = d2121 * d4343 - d4321 * d4321

    if fabs(denom) < MIN_SCALAR_VALUE:
        return None

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa = Point(
        x=p1.x + mua * p21.x,
        y=p1.y + mua * p21.y,
        z=p1.z + mua * p21.z,
        epsg_cd=epsg_cd
    )

    pb = Point(
        x=p3.x + mub * p43.x,
        y=p3.y + mub * p43.y,
        z=p3.z + mub * p43.z,
        epsg_cd=epsg_cd
    )

    return point_or_segment(
        point1=pa,
        point2=pb
    )



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

    def intersect(self) -> Optional[Union[Point, Segment]]:
        """
        Determines the optional point or segment intersection between the segment pair.

        :return: the optional point or segment intersection between the segment pair.
        :rtype: Optional[Union[Point, Segment]]

        Examples:
          >>> s2 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(-2,0,0), Point(-1,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect() is None
          True
          >>> s1 = Segment(Point(-2,0,0), Point(0,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Point(0.0000, 0.0000, 0.0000, 0.0000, -1)
          >>> s1 = Segment(Point(-2,0,0), Point(0.5,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(0.5000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(-2,0,0), Point(1,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(-2,0,0), Point(2,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0,0,0), Point(0.5,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(0.5000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0.25,0,0), Point(0.75,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.2500, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(0.7500, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0.25,0,0), Point(1,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.2500, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0.25,0,0), Point(1.25,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.2500, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0,0,0), Point(1.25,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.0000, 0.0000, 0.0000, 0.0000, -1), end_pt=Point(1.0000, 0.0000, 0.0000, 0.0000, -1))
          >>> s1 = Segment(Point(1,0,0), Point(1.25,0,0))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Point(1.0000, 0.0000, 0.0000, 0.0000, -1)
          >>> s2 = Segment(Point(0,0,0), Point(1,1,1))
          >>> s1 = Segment(Point(0.25,0.25,0.25), Point(0.75,0.75,0.75))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.2500, 0.2500, 0.2500, 0.0000, -1), end_pt=Point(0.7500, 0.7500, 0.7500, 0.0000, -1))
          >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,1.75,1.75))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Segment(start_pt=Point(0.2500, 0.2500, 0.2500, 0.0000, -1), end_pt=Point(1.0000, 1.0000, 1.0000, 0.0000, -1))
          >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,0,1.75))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Point(0.2500, 0.2500, 0.2500, 0.0000, -1)
          >>> s1 = Segment(Point(0.25,1,0.25), Point(0.75,0.75,0.75))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Point(0.7500, 0.7500, 0.7500, 0.0000, -1)
          >>> s2 = Segment(Point(-1,-1,-1), Point(1,1,1))
          >>> s1 = Segment(Point(-1,1,1), Point(1,-1,-1))
          >>> sp = SegmentPair(s1, s2)
          >>> sp.intersect()
          Point(0.0000, 0.0000, 0.0000, 0.0000, -1)

        """

        s1_startpt_inside = self._s1.segment_start_in(self._s2)
        s2_startpt_inside = self._s2.segment_start_in(self._s1)

        s1_endpt_inside = self._s1.segment_end_in(self._s2)
        s2_endpt_inside = self._s2.segment_end_in(self._s1)

        elements = [s1_startpt_inside, s2_startpt_inside, s1_endpt_inside, s2_endpt_inside]

        if all(elements):
            return self._s1.clone()

        if s1_startpt_inside and s1_endpt_inside:
            return self._s1.clone()

        if s2_startpt_inside and s2_endpt_inside:
            return self._s2.clone()

        if s1_startpt_inside and s2_startpt_inside:
            return point_or_segment(
                self._s1.start_pt,
                self._s2.start_pt
            )

        if s1_startpt_inside and s2_endpt_inside:
            return point_or_segment(
                self._s1.start_pt,
                self._s2.end_pt
            )

        if s1_endpt_inside and s2_startpt_inside:
            return point_or_segment(
                self._s2.start_pt,
                self._s1.end_pt

            )

        if s1_endpt_inside and s2_endpt_inside:
            return point_or_segment(
                self._s1.end_pt,
                self._s2.end_pt
            )

        if s1_startpt_inside:
            return self._s1.start_pt.clone()

        if s1_endpt_inside:
            return self._s1.end_pt.clone()

        if s2_startpt_inside:
            return self._s2.start_pt.clone()

        if s2_endpt_inside:
            return self._s2.end_pt.clone()

        shortest_segm_or_pt = shortest_segment_or_point(
            segment1=self._s1,
            segment2=self._s2
        )

        if not shortest_segm_or_pt:
            return None

        if isinstance(shortest_segm_or_pt, Point):
            return shortest_segm_or_pt
        else:
            return None


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

