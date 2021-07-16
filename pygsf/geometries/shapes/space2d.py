
import functools
import itertools

from typing import Optional, Union

import numbers

from math import fabs
import random
from array import array

import numpy as np

from .abstract import *
from ...mathematics.vectors2d import *
from ...utils.types import *


class Point2D(Point):
    """
    Cartesian point.
    Dimensions: 2D
    """

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real
                 ):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :param y: point y coordinate.
        """

        super(Point2D, self).__init__(x, y)

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y = Point2D(1,1)
          >>> x == 1
          True
          >>> y == 1
          True

        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return "Point2D({:.4f}, {:.4f})".format(self.x, self.y)

    def __eq__(self,
        another: 'Point2D'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point2D(1., 1.) == Point2D(1, 1)
          True
          >>> Point2D(1., 1.) == Point2D(1, 1)
          True
          >>> Point2D(1., 1.) == Point2D(1, -1)
          False
        """

        if not isinstance(another, Point2D):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            ]
        )

    def __ne__(self,
        another: 'Point2D'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point2D(1., 1.) != Point2D(0., 0.)
          True
          >>> Point2D(1., 1.) != Point2D(1, 1)
          False
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Return the individual values of the point.

        :return: double array of x, y values

        Examples:
          >>> Point2D(4, 3).a()
          (4.0, 3.0)
        """

        return self.x, self.y

    def __add__(self, another: 'Point2D') -> 'Point2D':
        """
        Sum of two points.

        :param another: the point to add
        :type another: Point2D
        :return: the sum of the two points
        :rtype: Point2D
        :raise: Exception

        Example:
          >>> Point2D(1, 0) + Point2D(0, 1)
          Point2D(1.0000, 1.0000)
          >>> Point2D(1, 1) + Point2D(-1, -1)
          Point2D(0.0000, 0.0000)
        """

        check_type(another, "Second point", Point2D)

        x0, y0 = self
        x1, y1 = another

        return Point2D(
            x=x0+x1,
            y=y0+y1
        )

    def __sub__(self,
        another: 'Point2D'
    ) -> 'Point2D':
        """Subtract two points.

        :param another: the point to subtract
        :type another: Point2D
        :return: the difference between the two points
        :rtype: Point2D
        :raise: Exception

        Example:
          >>> Point2D(1., 1.) - Point2D(1., 1.)
          Point2D(0.0000, 0.0000)
          >>> Point2D(1., 1.) - Point2D(1., 1.)
          Point2D(0.0000, 0.0000)
        """

        check_type(another, "Second point", Point2D)

        x0, y0 = self
        x1, y1 = another

        return Point2D(
            x=x0 - x1,
            y=y0 - y1
        )

    def clone(self) -> 'Point2D':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return Point2D(*self.a())

    def as_point2d(self) -> 'Point2D':
        """
        Convert a point to a 2D point.
        """

        return Point2D(
            x=self.x,
            y=self.y
        )

    def bounding_box(self) -> Optional[Tuple['Point2D', 'Point2D']]:

        return self.clone(), self.clone()

    def toXY(self) -> Tuple[numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of two values.

        :return: the spatial components (x, y).
        :rtype: a tuple of two floats.

        Examples:
          >>> Point2D(1, 0).toXY()
          (1.0, 0.0)
        """

        return self.x, self.y

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Point2D(1, 2).toArray(), np.array([ 1., 2.]))
          True
        """

        return np.asarray(self.toXY())

    def delta_x(self,
                another: 'Point2D'
                ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> Point2D(1, 2).delta_x(Point2D(4, 7))
          3.0
        """

        return another.x - self.x

    def delta_y(self,
                another: 'Point2D'
                ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point2D(1, 2).delta_y(Point2D(4, 7))
          5.0
        """

        return another.y - self.y

    def distance_2d(self,
                 another: 'Point2D'
                 ) -> numbers.Real:
        """
        Calculate horizontal (2D) distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :return: the 2D distance (when the two points have the same CRS).

        Examples:
          >>> Point2D(1., 1.).distance_2d(Point2D(4., 5.))
          5.0
        """

        return math.sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def distance(self,
                 another: 'Point2D'
                 ) -> numbers.Real:
        """
        Calculate horizontal (2D) distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the 2D distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point2D(1., 1.).distance(Point2D(4., 5.))
          5.0
        """

        check_type(another, "Second point", Point2D)

        return self.distance_2d(another)

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'Point2D':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point2D(1, 0).scale(2.5)
          Point2D(2.5000, 0.0000)
          >>> Point2D(1, 0).scale(2.5)
          Point2D(2.5000, 0.0000)
        """

        x, y = self.x * scale_factor, self.y * scale_factor
        return Point2D(x, y)

    def invert(self) -> 'Point2D':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point2D(1, 1).invert()
          Point2D(-1.0000, -1.0000)
          >>> Point2D(2, -1).invert()
          Point2D(-2.0000, 1.0000)
        """

        return self.scale(-1)

    def is_coincident(self,
                      other: 'Point2D',
                      tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                      ) -> bool:
        """
        Check spatial coincidence of two points, limiting to the horizontal (XY) plane.

        :param other: the point to compare.
        :type other: Point.
        :param tolerance: the maximum allowed distance between the two points.
        :type tolerance: numbers.Real.
        :return: whether the two points are coincident.
        :rtype: bool.
        :raise: Exception.

        Example:

        """

        check_type(other, "Second point", Point2D)

        return self.distance(other) <= tolerance

    def already_present(self,
                        pt_list: List['Point2D'],
                        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                        ) -> Optional[bool]:
        """
        Determines if a point is already in a given point list, using an optional distance separation,

        :param pt_list: list of points. May be empty.
        :type pt_list: List of Points.
        :param tolerance: optional maximum distance between near-coincident point pair.
        :type tolerance: numbers.Real.
        :return: True if already present, False otherwise.
        :rtype: optional boolean.
        """

        for pt in pt_list:
            if self.is_coincident(
                    pt,
                    tolerance=tolerance):
                return True
        return False

    def shift(self,
        *s
    ) -> Optional['Point2D']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point2D(1, 1).shift(0.5, 1.)
          Point2D(1.5000, 2.0000)
          >>> Point2D(1, 2).shift(0.5, 1.)
          Point2D(1.5000, 3.0000)
       """

        return Point2D(self.x + s[0], self.y + s[1])

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float =  MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: Point2D
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(2)]
        return cls(*vals)


class Segment2D(Segment):
    """
    Segment2D is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self,
                 start_pt: Point2D,
                 end_pt: Point2D
                 ):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :param end_pt: the end point.
        :return: the new segment instance if both points have the same georeferenced.
        :raises: CRSCodeException.
        """

        check_type(start_pt, "Start point", Point2D)
        check_type(end_pt, "End point", Point2D)
        if start_pt.distance(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        super(Segment2D, self).__init__()

        self._start_pt = start_pt
        self._end_pt = end_pt

    def __repr__(self) -> str:
        """
        Represents a Segment2D instance.

        :return: the Segment2D representation.
        :rtype: str.
        """

        return "Segment2D(start_pt={}, end_pt={})".format(
            self.start_pt,
            self.end_pt
        )

    def asPoints(self) -> List[Point2D]:
        """
        Return the segments as points.
        """

        return [self.start_pt, self.end_pt]

    def as_segment2d(self) -> 'Segment2D':
        """
        Convert a segment to a segment 2D.
        """

        return Segment2D(
            start_pt=self.start_pt.as_point2d(),
            end_pt=self.end_pt.as_point2d(),
        )

    def length_2d(self) -> numbers.Real:
        """
        Returns the horizontal length of the segment.

        :return: the horizontal length of the segment.
        :rtype: numbers.Real.
        """

        return self.start_pt.distance_2d(self.end_pt)

    @property
    def length(self) -> numbers.Real:
        """
        Returns the horizontal length of the segment.

        :return: the horizontal length of the segment.
        :rtype: numbers.Real.
        """

        return self.length_2d()

    def __iter__(self):
        """
        Return the elements of a Segment2D, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'Segment2D':

        return Segment2D(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment2D':

        if self.end_pt.x < self.start_pt.x:
            return Segment2D(self.end_pt, self.start_pt)
        else:
            return self.clone()

    def x_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.x < self.end_pt.x:
            return self.start_pt.x, self.end_pt.x
        else:
            return self.end_pt.x, self.start_pt.x

    def y_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.y < self.end_pt.y:
            return self.start_pt.y, self.end_pt.y
        else:
            return self.end_pt.y, self.start_pt.y

    def delta_x(self) -> numbers.Real:

        return self.end_pt.x - self.start_pt.x

    def delta_y(self) -> numbers.Real:

        return self.end_pt.y - self.start_pt.y

    def segment_2d_m(self) -> Optional[numbers.Real]:

        denom = self.end_pt.x - self.start_pt.x

        if denom == 0.0:
            return None

        return (self.end_pt.y - self.start_pt.y) / denom

    def segment_2d_p(self) -> Optional[numbers.Real]:

        s2d_m = self.segment_2d_m()

        if s2d_m is None:
            return None

        return self.start_pt.y - s2d_m * self.start_pt.x

    def shift(self,
              dx: numbers.Real,
              dy: numbers.Real
    ) -> 'Segment2D':
        """
        Shift a segment by dx and dy
        """

        return Segment2D(
            self.start_pt.shift(dx, dy),
            self.end_pt.shift(dx, dy)
        )

    def intersection_2d_pt(self,
                           another: 'Segment2D'
                           ) -> Optional[Point2D]:
        """

        :param another:
        :return:
        """

        check_type(another, "Second segment", Segment2D)

        s_len2d = self.length
        a_len2d = another.length

        if s_len2d == 0.0 or a_len2d == 0.0:
            return None

        if self.start_pt.x == self.end_pt.x:  # self segment parallel to y axis
            x0 = self.start_pt.x
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            if m1 is None:
                return None
            y0 = m1 * x0 + p1
        elif another.start_pt.x == another.end_pt.x:  # another segment parallel to y axis
            x0 = another.start_pt.x
            m1, p1 = self.segment_2d_m(), self.segment_2d_p()
            if m1 is None:
                return None
            y0 = m1 * x0 + p1
        else:  # no segment parallel to y axis
            m0, p0 = self.segment_2d_m(), self.segment_2d_p()
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            if m0 is None or m1 is None:
                return None
            x0 = (p1 - p0) / (m0 - m1)
            y0 = m0 * x0 + p0

        return Point2D(x0, y0)

    def contains_pt(self,
                    pt: Point2D
                    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment2D(Point2D(0, 0), Point2D(1, 0))
          >>> segment.contains_pt(Point2D(0, 0))
          True
          >>> segment.contains_pt(Point2D(1, 0))
          True
          >>> segment.contains_pt(Point2D(0.5, 0))
          True
          >>> segment.contains_pt(Point2D(0.5, 0.00001))
          False
          >>> segment.contains_pt(Point2D(1.00001, 0))
          False
          >>> segment.contains_pt(Point2D(0.000001, 0))
          True
          >>> segment.contains_pt(Point2D(-0.000001, 0))
          False
          >>> segment.contains_pt(Point2D(0.5, 1000))
          False
          >>> segment = Segment2D(Point2D(0, 0), Point2D(0, 1))
          >>> segment.contains_pt(Point2D(0, 0))
          True
          >>> segment.contains_pt(Point2D(0, 0.5))
          True
          >>> segment.contains_pt(Point2D(0, 1))
          True
          >>> segment.contains_pt(Point2D(0, 1.5))
          False
          >>> segment = Segment2D(Point2D(0, 0), Point2D(1, 1))
          >>> segment.contains_pt(Point2D(0.5, 0.5))
          True
          >>> segment.contains_pt(Point2D(1, 1))
          True
          >>> segment = Segment2D(Point2D(1, 2), Point2D(9, 8))
          >>> segment.contains_pt(segment.point_at_factor(0.745))
          True
          >>> segment.contains_pt(segment.point_at_factor(1.745))
          False
          >>> segment.contains_pt(segment.point_at_factor(-0.745))
          False
          >>> segment.contains_pt(segment.point_at_factor(0))
          True
        """

        check_type(pt, "Point", Point2D)

        segment_length = self.length
        length_startpt_pt = self.start_pt.distance(pt)
        length_endpt_pt = self.end_pt.distance(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def fast_2d_contains_pt(self,
                            pt2d
                            ) -> bool:
        """
        Deprecated. Use 'contains_pt'.

        to work properly, this function requires that the pt lies on the line defined by the segment
        """

        range_x = self.x_range
        range_y = self.y_range

        if range_x()[0] <= pt2d.x <= range_x()[1] or \
                range_y()[0] <= pt2d.y <= range_y()[1]:
            return True
        else:
            return False

    def point_at_factor(self,
                        scale_factor: numbers.Real
                        ) -> Point2D:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        and 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :return: Point at scale factor

        Examples:
          >>> s = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s.point_at_factor(0)
          Point2D(0.0000, 0.0000)
          >>> s.point_at_factor(0.5)
          Point2D(0.5000, 0.0000)
          >>> s.point_at_factor(1)
          Point2D(1.0000, 0.0000)
          >>> s.point_at_factor(-1)
          Point2D(-1.0000, 0.0000)
          >>> s.point_at_factor(-2)
          Point2D(-2.0000, 0.0000)
          >>> s.point_at_factor(2)
          Point2D(2.0000, 0.0000)
          >>> s = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s.point_at_factor(0.5)
          Point2D(0.5000, 0.5000)
          >>> s = Segment2D(Point2D(0,0), Point2D(4,0))
          >>> s.point_at_factor(7.5)
          Point2D(30.0000, 0.0000)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor

        return Point2D(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy
        )

    def point_at_distance(self,
                        distance: numbers.Real
                        ) -> Point2D:
        """
        Returns a point aligned with the segment
        and lying at given distance from the segment start.

        :param distance: the distance from segment start
        :return: point at provided distance from segment start

        Examples:
          >>> s = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s.point_at_distance(0)
          Point2D(0.0000, 0.0000)
          >>> s.point_at_distance(0.5)
          Point2D(0.5000, 0.0000)
          >>> s.point_at_distance(1)
          Point2D(1.0000, 0.0000)
          >>> s.point_at_distance(2)
          Point2D(2.0000, 0.0000)
          >>> s = Segment2D(Point2D(0,0), Point2D(3, 4))
          >>> s.point_at_distance(0)
          Point2D(0.0000, 0.0000)
          >>> s.point_at_distance(2.5)
          Point2D(1.5000, 2.0000)
          >>> s.point_at_distance(5)
          Point2D(3.0000, 4.0000)
          >>> s.point_at_distance(10)
          Point2D(6.0000, 8.0000)
        """

        scale_factor = distance / self.length
        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor

        return Point2D(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy
        )

    @property
    def start_pt(self) -> Point2D:

        return self._start_pt

    @property
    def end_pt(self) -> Point2D:

        return self._end_pt

    '''
    def pointProjection(self,
                        point: Point
                        ) -> Point:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = Segment2D(start_pt=Point2D(0,0), end_pt=Point2D(1,0))
          >>> p = Point2D(0.5, 1)
          >>> s.pointProjection(p)
          Point2D(0.5000, 0.0000)
          >>> s = Segment2D(start_pt=Point2D(0,0), end_pt=Point2D(4,0))
          >>> p = Point2D(7.5, 19.2)
          >>> s.pointProjection(p)
          Point2D(7.5000, 0.0000)
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        other_segment = Segment2D(
            self.start_pt,
            point
        )

        scale_factor = self.vector().scalarProjection(other_segment.vector()) / self.length2D()
        return self.pointAt(scale_factor)
    '''

    '''
    def pointDistance(self,
                      point: Point
                      ) -> numbers.Real:
        """
        Returns the point distance to the segment.

        :param point: the point to calculate the distance with
        :type point: Point
        :return: the distance of the point to the segment
        :rtype: numbers.Real

        Examples:
          >>> s = Segment2D(Point2D(0, 0), Point2D(0, 0))
          >>> s.pointDistance(Point2D(-17.2, 0.0))
          17.2
          >>> s.pointDistance(Point2D(-17.2, 1.22))
          17.24321315764553
        """

        check_type(point, "Input point", Point)

        # check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.dist2DWith(point_projection)
    '''

    def point_distance_from_start(self,
                                  point: Point2D
                                  ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point2D)

        if not self.contains_pt(point):
            return None

        return self.start_pt.distance(point)

    def point_signed_s(self,
               point: Point2D
               ) -> Optional[numbers.Real]:
        """
        Calculates the signed distance of the point along the segment.
        A zero value is for a point coinciding with the start point.

        :param point: the point to calculate the optional distance in the segment.
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point2D)

        if not self.contains_pt(point):
            return None

        return self.start_pt.distance(point)

    def scale(self,
              scale_factor
              ) -> 'Segment2D':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point2D
        """

        end_pt = self.point_at_factor(scale_factor)

        return Segment2D(
            self.start_pt,
            end_pt)

    def as_vector(self) -> Vect2D:
        """
        Convert a segment to a vector.
        """

        return Vect2D(
            x=self.delta_x(),
            y=self.delta_y()
        )

    def as_versor(self) -> Vect2D:
        """
        Convert a segment to a versor.
        """

        return self.as_vector().versor()

    def densify_as_line2d(self,
                          densify_distance: numbers.Real
                          ) -> 'Line2D':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: the densify distance
        :return: a line
        """

        return Line2D.fromPoints(self.densify_as_points2d(densify_distance=densify_distance))

    def densify_as_points2d(self,
                            densify_distance: numbers.Real,
                            start_offset: numbers.Real = 0.0
                            ) -> List[Point2D]:
        """
        Densify a segment as a list of points, by using the provided densify distance.

        :param densify_distance: the densify distance
        :return: a list of points
        """

        length2d = self.length

        vers_2d = self.as_versor()

        step_vector = vers_2d.scale(densify_distance)

        if start_offset == 0.0:
            start_point = self.start_pt
        else:
            start_point = self.point_at_distance(start_offset)

        pts = [start_point]

        n = 0
        while True:
            n += 1
            shift_vector = step_vector.scale(n)
            new_pt = start_point.shift(shift_vector.x, shift_vector.y)
            distance = self.start_pt.distance(new_pt)
            if distance >= length2d:
                break
            pts.append(new_pt)

        pts.append(self.end_pt)

        return pts

    def densify_as_array1d(self,
                           densify_distance: numbers.Real
                           ) -> array:
        """
        Defines the array storing the incremental lengths according to the provided densify distance.

        :param densify_distance: the densify distance.
        :return: array storing incremental steps, with the last step being equal to the segment length.
        """

        if not isinstance(densify_distance, numbers.Real):
            raise Exception("Densify distance must be float or int")

        if not math.isfinite(densify_distance):
            raise Exception("Densify distance must be finite")

        if densify_distance <= 0.0:
            raise Exception("Densify distance must be positive")

        segment_length = self.length

        s_list = []
        n = 0
        length = n * densify_distance

        while length < segment_length:
            s_list.append(length)
            n += 1
            length = n * densify_distance

        s_list.append(segment_length)

        return array('d', s_list)

    '''
    def vertical_plane(self) -> Optional[CPlane]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length2D() == 0.0:
            return None

        # arbitrary point on the same vertical as end point

        section_final_pt_up = self.end_pt.shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane.fromPoints(
            pt1=self.start_pt,
            pt2=self.end_pt,
            pt3=section_final_pt_up)
    '''

    def same_start(self,
                   another: 'Segment2D',
                   tol: numbers.Real = 1e-12
                   ) -> bool:
        """
        Check whether the two segments have the same start point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(0,0), Point2D(0,1))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt.is_coincident(
            other=another.start_pt,
            tolerance=tol
        )

    def same_end(self,
                 another: 'Segment2D',
                 tol: numbers.Real = 1e-12
                 ) -> bool:
        """
        Check whether the two segments have the same end point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(2,0), Point2D(1,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt.is_coincident(
            other=another.end_pt,
            tolerance=tol)

    def conn_to_other(self,
                      another: 'Segment2D',
                      tol: numbers.Real = 1e-12
                      ) -> bool:
        """
        Check whether the first segment is sequentially connected to the second one.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(1,0), Point2D(2,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt.is_coincident(
            other=another.start_pt,
            tolerance=tol)

    def other_connected(self,
                        another: 'Segment2D',
                        tol: numbers.Real = 1e-12
                        ) -> bool:
        """
        Check whether the second segment is sequentially connected to the first one.

        :param another: a segment to check for.
        :type another: Segment2D.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-1,0), Point2D(0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt.is_coincident(
            other=self.start_pt,
            tolerance=tol)

    def segment_start_in(self,
                         another: 'Segment2D'
                         ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-0.5,0), Point2D(0.5,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment2D(Point2D(0,1), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = Segment2D(Point2D(-1,-1), Point2D(1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt)

    def segment_end_in(self,
                       another: 'Segment2D'
                       ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment2D.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
          >>> s2 = Segment2D(Point2D(-0.5,0), Point2D(0.5,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment2D(Point2D(0,0), Point2D(1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment2D(Point2D(0,1), Point2D(1,1))
          >>> s2 = Segment2D(Point2D(1,1), Point2D(0.5,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = Segment2D(Point2D(-1,-1), Point2D(1,1))
          >>> s2 = Segment2D(Point2D(0,2), Point2D(2,0))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt)

    @classmethod
    def random(cls,
               lower_boundary: float = -MAX_SCALAR_VALUE,
               upper_boundary: float = MAX_SCALAR_VALUE):
        """
        Creates a random segment.

        :return: random segment
        :rtype: Segment2D
        """

        return cls(
            start_pt=Point2D.random(lower_boundary, upper_boundary),
            end_pt=Point2D.random(lower_boundary, upper_boundary)
        )

    '''
    def vertical_plane(self) -> Optional[CPlane3D]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length() == 0.0:  # collapsed segment
            return None

        # arbitrary point on the same vertical as end point

        p1 = Point3D(
            x=self.start_pt.x,
            y=self.start_pt.y,
            z=0.0
        )
        p2 = Point3D(
            x=self.end_pt.x,
            y=self.end_pt.y,
            z=0.0
        )
        p3 = Point3D(
            x=self.end_pt.x,
            y=self.end_pt.y,
            z=1000.0
        )
        return CPlane3D.fromPoints(
            pt1=p1,
            pt2=p2,
            pt3=p3)

    def vector(self) -> Vect:

        return Vect3D(self.delta_x(),
                    self.delta_y(),
                    0
        )
    '''

    def vector(self) -> Vect2D:

        return Vect2D(self.delta_x(),
                      self.delta_y()
                      )


class PointSegmentCollection2D(list):
    """
    Collection of point or segment elements.

    """

    def __init__(
            self,
            geoms: Optional[List[Union[Point2D, Segment2D]]] = None,
            # epsg_code: Optional[numbers.Integral] = None
    ):

        if geoms is not None:

            for geom in geoms:
                check_type(geom, "Spatial element", (Point2D, Segment2D))

        """
        if epsg_code is not None:
            check_type(
                var=epsg_code,
                name="EPSG code",
                expected_types=numbers.Integral
            )

        if geoms is not None and epsg_code is not None:

            for geom in geoms:
                check_epsg(
                    spatial_element=geom,
                    epsg_code=epsg_code
                )

        elif geoms is not None and len(geoms) > 0:

            epsg_code = geoms[0].epsg_code()
        """

        if geoms is not None and len(geoms) > 0:

            super(PointSegmentCollection2D, self).__init__(geoms)

        else:

            super(PointSegmentCollection2D, self).__init__()

        # self.epsg_code = epsg_code

    def append(self,
               spatial_element: Union[Point2D, Segment2D]
               ) -> None:

        check_type(
            var=spatial_element,
            name="Spatial element",
            expected_types=(Point2D, Segment2D)
        )

        """
        if self.epsg_code is not None:

            check_epsg(
                spatial_element=spatial_element,
                epsg_code=self.epsg_code
            )

        else:

            self.epsg_code = spatial_element.epsg_code()
        """

        self.append(spatial_element)


def point_or_segment2d(
        point1: Point2D,
        point2: Point2D,
        tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Union[Point2D, Segment2D]:
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

    check_type(point1, "First point", Point2D)
    check_type(point2, "Second point", Point2D)

    if point1.distance(point2) <= tol:
        return Point2D(
            x=(point1.x + point2.x) / 2,
            y=(point1.y + point2.y) / 2
        )
    else:
        return Segment2D(
            start_pt=point1,
            end_pt=point2
        )


def shortest_segment_or_point2d(
    first_segment: Segment2D,
    second_segment: Segment2D,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Segment2D, Point2D]]:
    """
    TODO: check correct implementation for 2D case, since it derives from 3D implementation.

    Calculates the optional shortest segment - or the intersection point - between two lines represented by two segments.

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

    :param first_segment: the first segment
    :param second_segment: the second segment
    :param tol: tolerance value for collapsing a segment into the midpoint.
    :return: the optional shortest segment or an intersection point.
    """

    check_type(second_segment, "Second Cartesian line", Segment2D)

    p1 = first_segment.start_pt
    p2 = first_segment.end_pt

    p3 = second_segment.start_pt
    p4 = second_segment.end_pt

    p13 = Point2D(
        x=p1.x - p3.x,
        y=p1.y - p3.y
    )

    p43 = Point2D(
        x=p4.x - p3.x,
        y=p4.y - p3.y
    )

    if p43.is_coincident(Point2D(0, 0)):
        return None

    p21 = Point2D(
        x=p2.x - p1.x,
        y=p2.y - p1.y
    )

    if p21.is_coincident(Point2D(0, 0)):
        return None

    d1343 = p13.x * p43.x + p13.y * p43.y
    d4321 = p43.x * p21.x + p43.y * p21.y
    d1321 = p13.x * p21.x + p13.y * p21.y
    d4343 = p43.x * p43.x + p43.y * p43.y
    d2121 = p21.x * p21.x + p21.y * p21.y

    denom = d2121 * d4343 - d4321 * d4321

    if fabs(denom) < MIN_SCALAR_VALUE:
        return None

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa = Point2D(
        x=p1.x + mua * p21.x,
        y=p1.y + mua * p21.y
    )

    pb = Point2D(
        x=p3.x + mub * p43.x,
        y=p3.y + mub * p43.y
    )

    intersection = point_or_segment2d(
        point1=pa,
        point2=pb,
        tol=tol
    )

    return intersection


def intersect_segments2d(
    segment1: Segment2D,
    segment2: Segment2D,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Point2D, Segment2D]]:
    """
    Determines the optional point or segment intersection between the segment pair.

    :param segment1: the first segment
    :param segment2: the second segment
    :param tol: the distance tolerance for collapsing a intersection segment into a point
    :return: the optional point or segment intersection between the segment pair.

    Examples:
      >>> s2 = Segment2D(Point2D(0,0), Point2D(1,0))
      >>> s1 = Segment2D(Point2D(0,0), Point2D(1,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(-2,0), Point2D(-1,0))
      >>> intersect_segments2d(s1, s2) is None
      True
      >>> s1 = Segment2D(Point2D(-2,0), Point2D(0,0))
      >>> intersect_segments2d(s1, s2)
      Point2D(0.0000, 0.0000)
      >>> s1 = Segment2D(Point2D(-2,0), Point2D(0.5,0.0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(0.5000, 0.0000))
      >>> s1 = Segment2D(Point2D(-2,0), Point2D(1,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(-2,0), Point2D(2,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(0,0), Point2D(0.5,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(0.5000, 0.0000))
      >>> s1 = Segment2D(Point2D(0.25,0), Point2D(0.75,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.2500, 0.0000), end_pt=Point2D(0.7500, 0.0000))
      >>> s1 = Segment2D(Point2D(0.25,0), Point2D(1,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.2500, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(0.25,0), Point2D(1.25,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.2500, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(0,0), Point2D(1.25,0))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.0000, 0.0000), end_pt=Point2D(1.0000, 0.0000))
      >>> s1 = Segment2D(Point2D(1,0), Point2D(1.25,0))
      >>> intersect_segments2d(s1, s2)
      Point2D(1.0000, 0.0000)
      >>> s2 = Segment2D(Point2D(0,0), Point2D(1,1))
      >>> s1 = Segment2D(Point2D(0.25,0.25), Point2D(0.75,0.75))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.2500, 0.2500), end_pt=Point2D(0.7500, 0.7500))
      >>> s1 = Segment2D(Point2D(0.25,0.25), Point2D(1.75,1.75))
      >>> intersect_segments2d(s1, s2)
      Segment2D(start_pt=Point2D(0.2500, 0.2500), end_pt=Point2D(1.0000, 1.0000))
      >>> s1 = Segment2D(Point2D(0.25,0.25), Point2D(1.75,0))
      >>> intersect_segments2d(s1, s2)
      Point2D(0.2500, 0.2500)
      >>> s1 = Segment2D(Point2D(0.25,1), Point2D(0.75,0.75))
      >>> intersect_segments2d(s1, s2)
      Point2D(0.7500, 0.7500)
      >>> s2 = Segment2D(Point2D(-1,-1), Point2D(1,1))
      >>> s1 = Segment2D(Point2D(-1,1), Point2D(1,-1))
      >>> intersect_segments2d(s1, s2)
      Point2D(0.0000, 0.0000)
    """

    check_type(segment1, "First segment", Segment2D)
    check_type(segment2, "Second segment", Segment2D)

    #check_crs(segment1, segment2)

    s1_startpt_inside = segment1.segment_start_in(segment2)
    s2_startpt_inside = segment2.segment_start_in(segment1)

    s1_endpt_inside = segment1.segment_end_in(segment2)
    s2_endpt_inside = segment2.segment_end_in(segment1)

    elements = [s1_startpt_inside, s2_startpt_inside, s1_endpt_inside, s2_endpt_inside]

    if all(elements):
        return segment1.clone()

    if s1_startpt_inside and s1_endpt_inside:
        return segment1.clone()

    if s2_startpt_inside and s2_endpt_inside:
        return segment2.clone()

    if s1_startpt_inside and s2_startpt_inside:
        return point_or_segment2d(
            segment1.start_pt,
            segment2.start_pt,
            tol=tol
        )

    if s1_startpt_inside and s2_endpt_inside:
        return point_or_segment2d(
            segment1.start_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_startpt_inside:
        return point_or_segment2d(
            segment2.start_pt,
            segment1.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_endpt_inside:
        return point_or_segment2d(
            segment1.end_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_startpt_inside:
        return segment1.start_pt.clone()

    if s1_endpt_inside:
        return segment1.end_pt.clone()

    if s2_startpt_inside:
        return segment2.start_pt.clone()

    if s2_endpt_inside:
        return segment2.end_pt.clone()

    shortest_segm_or_pt = shortest_segment_or_point2d(
        segment1,
        segment2,
        tol=tol
    )

    if not shortest_segm_or_pt:
        return None

    if not isinstance(shortest_segm_or_pt, Point2D):
        return None

    inters_pt = shortest_segm_or_pt

    if not segment1.contains_pt(inters_pt):
        return None

    if not segment2.contains_pt(inters_pt):
        return None

    return inters_pt


class Line2D(Line):
    """
    A list of Point objects.
    """

    def __init__(self,
        x_seq: Optional[Union[Sequence[float], np.ndarray]] = None,
        y_seq: Optional[Union[Sequence[float], np.ndarray]] = None
    ):
        """
        Creates the Line2D instance.
        """

        if x_seq is None and y_seq is not None:
            raise Exception(f"x_seq is None while y_seq is {type(y_seq)}")

        if x_seq is not None and y_seq is None:
            raise Exception(f"x_seq is {type(x_seq)} while y_seq is None")

        if x_seq is None:
            x_array = None
            y_array = None
        else:
            x_array = np.asarray(x_seq)
            y_array = np.asarray(y_seq)
            if x_array.ndim != 1:
                raise Exception(f"X input must have a dimension of 1, not {x_array.ndim}")
            if y_array.ndim != 1:
                raise Exception(f"Y input must have a dimension of 1, not {y_array.ndim}")
            if len(x_array) != len(y_array):
                raise Exception(f"X input has length {len(x_array)}, while y input has length {len(y_array)}")

        super(Line2D, self).__init__()

        self._x_array = x_array
        self._y_array = y_array

    def x_array(self):

        return self._x_array

    def y_array(self):

        return self._y_array

    def xy_zipped(self) -> List[Tuple[numbers.Real, numbers.Real]]:

        return [(x, y) for x, y in zip(self.x_list(), self.y_list())]

    @classmethod
    def fromPoints(cls,
        pts: Optional[List[Point2D]] = None
    ):
        """
        Creates the Line2D instance.

        :param pts: a list of points
        :return: a Line2D instance.
        """

        if pts is None:

            x_list = None
            y_list = None

        else:

            x_list = []
            y_list = []

            check_type(pts, "List", list)
            for pt in pts:
                check_type(pt, "Point", Point2D)
                x_list.append(pt.x)
                y_list.append(pt.y)

        return cls(
            x_list,
            y_list
        )

    @classmethod
    def fromCoordinates(cls,
                        coordinates: List[List[numbers.Real]]
                        ) -> 'Line2D':
        """
        Create a Line2D instance from a list of x and y values.

        Example:
          >>> Line2D.fromCoordinates([[0, 0], [1, 0], [0, 1]])
          Line2D with 3 points: (0.0000, 0.0000) ... (0.0000, 1.0000)
        """

        x_vals = []
        y_vals = []

        for x, y in coordinates:
            x_vals.append(x)
            y_vals.append(y)

        return cls(
            x_vals,
            y_vals
        )

    def start_pt(self) -> Optional[Point2D]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        """

        return self.pt(0) if len(self.x_array()) > 0 else None

    def end_pt(self) -> Optional[Point2D]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        """

        return self.pt(-1) if len(self.x_array()) > 0 else None

    def values_at(self,
        ndx: numbers.Integral
    ) -> Tuple[float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x and y values
        """

        return (
            self.x_array()[ndx],
            self.y_array()[ndx]
        )

    def num_pts(self) -> numbers.Integral:
        """
        Returns the number of line points.
        """

        if self.x_array() is None:
            return 0
        else:
            return len(self.x_array())

    def pt(self,
           ndx: numbers.Integral) -> Optional[Point2D]:
        """
        Returns the point at given index.
        """

        if self.num_pts() == 0:
            return None
        elif ndx < self.num_pts():
            return Point2D(
                x=self.x_array()[ndx],
                y=self.y_array()[ndx]
            )
        else:
            return None

    def pts(self) -> List[Point2D]:

        return [Point2D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def length(self):

        return self.length_2d()

    def segment(self,
        ndx: numbers.Integral
    ) -> Optional[Segment2D]:
        """
        Returns the optional segment at index ndx.

        :param ndx: the segment index.
        :type ndx: numbers.Integral
        :return: the optional segment
        :rtype: Optional[Segment2D]
        """

        start_pt = self.pt(ndx)
        end_pt = self.pt(ndx + 1)

        if start_pt.is_coincident(end_pt):
            return None
        else:
            return Segment2D(
                start_pt=self.pt(ndx),
                end_pt=self.pt(ndx + 1)
            )

    def __repr__(self) -> str:
        """
        Represents a Line instance as a shortened text.

        :return: a textual shortened representation of a Line instance.
        :rtype: str.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty Line2D"
        else:
            x1, y1 = self.start_pt()
            if num_points == 1:
                txt = "Line2D with unique point: {.4f}, {.4f}".format(x1, y1)
            else:
                x2, y2 = self.end_pt()
                txt = "Line2D with {} points: ({:.4f}, {:.4f}) ... ({:.4f}, {:.4f})".format(num_points, x1, y1, x2, y2)

        return txt

    def add_pt(self, pt: Point2D):
        """
        In-place transformation of the original Line2D instance
        by adding a new point at the end.

        :param pt: the point to add
        """

        check_type(pt, "Point", Point2D)

        if self._x_array is None:
            self._x_array = np.array([pt.x])
            self._y_array = np.array([pt.y])
        else:
            self._x_array = np.append(self._x_array, pt.x)
            self._y_array = np.append(self._y_array, pt.y)

    def add_pts(self, pt_list: List[Point2D]):
        """
        In-place transformation of the original Line instance
        by adding a new set of points at the end.

        :param pt_list: list of Points.
        """

        for pt in pt_list:
            check_type(pt, "Point", Point2D)

        x_vals = [pt.x for pt in pt_list]
        y_vals = [pt.y for pt in pt_list]

        if self._x_array is None:
            self._x_array = np.array(x_vals)
            self._y_array = np.array(y_vals)
        else:
            self._x_array = np.append(self._x_array, x_vals)
            self._y_array = np.append(self._y_array, y_vals)

    def remove_coincident_points(self) -> Optional['Line2D']:
        """
        Remove coincident successive points.
        """

        if self.num_pts() == 0:
            return

        new_line = Line2D.fromPoints(
            pts=[self.pt(0)]
        )

        for ndx in range(1, self.num_pts()):
            if not self.pt(ndx).is_coincident(new_line.pt(-1)):
                new_line.add_pt(self.pt(ndx))

        return new_line

    def as_segment(self) -> Segment2D:
        """Return the segment defined by line start and end points"""

        return Segment2D(self.start_pt(), self.end_pt())

    def segments(self) -> List[Segment2D]:
        """
        Convert to a list of segments 2d.

        :return: list of Segment2D objects
        """

        pts_pairs = zip(self.pts()[:-1], self.pts()[1:])

        segments = [Segment2D(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    def densify_as_equally_spaced_points2d(self,
       sample_distance: numbers.Real
       ) -> Tuple[List[Point2D], List[numbers.Real], numbers.Real]:
        """
        Densify a line into a set of Point2D instances, each equally spaced along the line
        (so that at corners 2D distance between points is less than 'sample_distance').

        """

        points = []
        break_distances = [0.0]
        segment_start_offset = 0.0
        for segment in self.segments():
            break_distances.append(segment.length)
            segment_points = segment.densify_as_points2d(start_offset=segment_start_offset)
            points.extend(segment_points[:-1])
            segment_start_offset = sample_distance - segment_points[-2].distance(segment_points[-1])

        return points, break_distances, sample_distance

    def densify_as_line2d(self, sample_distance) -> 'Line2D':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line2D instance
        """

        if sample_distance <= 0.0:
            raise Exception(f"Sample distance must be positive. {sample_distance} received")

        segments = self.segments()

        densified_line_list = [segment.densify_as_line2d(sample_distance) for segment in segments]

        densifyied_points = []
        for line in densified_line_list:
            densifyied_points += line.pts()

        densified_line = Line2D(densifyied_points)
        densified_line_wo_coinc_pts = densified_line.remove_coincident_points()

        return densified_line_wo_coinc_pts

    def join(self, another) -> 'Line2D':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line2D(self.pts() + another.pts())

    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).distance_2d(self.pt(ndx + 1))
        return length

    def step_lengths_2d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 2D distances.
        It is the distance between a point and its previous one.
        The list has the same length as the source point list.

        :return: the individual 2D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self.pt(ndx).distance(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list

    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))

    def reversed(self) -> 'Line2D':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        pts = [pt.clone() for pt in self.pts()]
        pts.reverse()

        return Line2D.fromPoints(
            pts=pts
        )

    def extremes_distance_2d(self) -> numbers.Real:
        """
        Calculate the 2D distance between start and end points.

        :return: the 2D distance between start and end points
        """

        return self.end_pt().distance(self.start_pt())

    def is_closed_2d(self,
                     tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                     ) -> bool:
        """
        Determine if the line is 2D-closed.

        :param tolerance: the tolerance for considering the line closed
        :return: whether the line is to be considered 2D-closed
        """

        return self.end_pt().is_coincident(self.start_pt(), tolerance=tolerance)

    def walk_backward(self) -> 'Line2D':
        """
        Create a new line by walking the line backward from the last point up to the first and thus closing it.

        :return: a closed line with zero area
        :rtype: 'Line'
        """

        return Line2D(
            np.concatenate(self.x_array(), np.flip(self.x_array())[1:]),
            np.concatenate(self.y_array(), np.flip(self.y_array())[1:])
        )

    def clone(self) -> 'Line2D':
        """
        Clone a line.

        :return: the cloned line
        :rtype: Line2D
        """

        return Line2D(
            np.copy(self.x_array()),
            np.copy(self.y_array())
        )

    def x_list(self) -> List[numbers.Real]:

        return list(self.x_array())

    def x_min(self):
        return np.nanmin(self.x_array())

    def x_max(self):
        return np.nanmax(self.x_array())

    def x_minmax(self):

        return self.x_min(), self.x_max()

    def x_mean(self):
        return np.nanmean(self.x_array())

    def x_var(self):

        return np.nanvar(self.x_array())

    def x_std(self):

        return np.nanstd(self.x_array())

    def y_list(self) -> List[numbers.Real]:

        return list(self.y_array())

    def y_min(self):
        return np.nanmin(self.y_array())

    def y_max(self):
        return np.nanmax(self.y_array())

    def y_minmax(self):

        return self.y_min(), self.y_max()

    def y_mean(self):

        return np.nanmean(self.y_array())

    def y_var(self):

        return np.nanvar(self.y_array())

    def y_std(self):

        return np.nanstd(self.y_array())

    def bounding_box(self) -> Optional[Tuple[Point2D, Point2D]]:
        """ The shape bounding box"""

        if self.length == 0:
            return None

        lower_left_point = Point2D(
            x=self.x_min(),
            y=self.y_min()
        )

        upper_right_point = Point2D(
            x=self.x_max(),
            y=self.y_max()
        )

        return lower_left_point, upper_right_point

    def close_2d(self) -> 'Line2D':
        """
        Return a line that is 2D-closed.

        :return: a 2D-closed line
        :rtype: Line2D
        """

        line = self.clone()

        if not line.is_closed_2d():

            line.add_pt(line.start_pt())

        return line

    def intersect_segment(self,
                          segment: Segment2D
                          ) -> Optional[PointSegmentCollection2D]:
        """
        Calculates the possible intersection between the line and a provided segment.

        :param segment: the input segment
        :return: the optional intersections, points or segments
        :raise: Exception
        """

        check_type(segment, "Input segment", Segment2D)

        intersections = [intersect_segments2d(curr_segment, segment) for curr_segment in self.segments() if curr_segment is not None]
        intersections = list(filter(lambda val: val is not None, intersections))
        intersections = PointSegmentCollection2D(intersections)

        return intersections


class MultiLine2D(object):
    """
    MultiLine is a list of Line objects
    """

    def __init__(self, lines_list=None):

        if lines_list is None:
            lines_list = []
        self._lines = lines_list

    @property
    def lines(self):

        return self._lines

    def add(self, line):

        return MultiLine2D(self.lines + [line])

    def clone(self):

        return MultiLine2D(self.lines)

    @property
    def num_parts(self):

        return len(self.lines)

    @property
    def num_points(self):

        num_points = 0
        for line in self.lines:
            num_points += line.num_pts

        return num_points

    @property
    def x_min(self):

        return np.nanmin([line.x_min for line in self.lines])

    @property
    def x_max(self):

        return np.nanmax([line.x_max for line in self.lines])

    @property
    def y_min(self):

        return np.nanmin([line.y_min for line in self.lines])

    @property
    def y_max(self):

        return np.nanmax([line.y_max for line in self.lines])

    def is_continuous(self):

        for line_ndx in range(len(self._lines) - 1):
            if not self.lines[line_ndx].pts[-1].is_coincident(self.lines[line_ndx + 1].pts[0]) or \
               not self.lines[line_ndx].pts[-1].is_coincident(self.lines[line_ndx + 1].pts[-1]):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines) - 1):
            if not self.lines[line_ndx].pts[-1].is_coincident(self.lines[line_ndx + 1].pts[0]):
                return False

        return True

    def to_line(self):

        return Line2D([point for line in self.lines for point in line.pts()])

    '''
    def crs_project(self, srcCrs, destCrs):

        lines = []
        for line in self.lines:
            lines.append(line.crs_project(srcCrs, destCrs))

        return MultiLine4D(lines)
    '''

    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines:
            lDensifiedLines.append(line.densify_as_line2d(sample_distance))

        return MultiLine2D(lDensifiedLines)

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines:
            cleaned_lines.append(line.remove_coincident_points())

        return MultiLine2D(cleaned_lines)


class Circle2D(Shape):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 r: numbers.Real
                 ):

        self._x = float(x)
        self._y = float(y)
        self._r = float(r)

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current circle.

        :return: x coordinate.

        Examples:
          >>> Circle2D(4, 3, 2).x
          4.0
          >>> Circle2D(-0.39, 3, 7).x
          -0.39
        """

        return self._x

    @property
    def y(self) -> numbers.Real:
        """
        Return the y coordinate of the current circle.

        :return: y coordinate.

        Examples:
          >>> Point2D(4, 3).y
          3.0
          >>> Point2D(-0.39, 17.42).y
          17.42
        """

        return self._y

    @property
    def radius(self):
        return self._r

    def area(self):
        return math.pi * self._r * self._r

    def length(self):
        return 2.0 * math.pi * self._r

    def clone(self):
        return Circle2D(self._x, self._y, self._r)


class Triangle2D(Polygon):

    def area(self):
        pass

    def length(self):
        pass

    @property
    def num_side(self):
        """Return numer of sides"""

        return 3


class Quadrilateral2D(Polygon, metaclass=abc.ABCMeta):

    def area(self):
        pass

    def length(self):
        pass

    @property
    def num_side(self):
        """Return numer of sides"""

        return 4


class Square2D(Quadrilateral2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 side: numbers.Real,
                 rotation: numbers.Real
                 ):
        """

        :param x: x coordinate of square center
        :param y: y coordinate of square center
        :param side: square side
        :param rotation: square rotation, counter-clockwise, decimal degrees
        """

        self._x = x
        self._y = y
        self._side = side
        self._cc_rotat = rotation % 360.0

    def area(self):
        return self._side * self._side

    def length(self):
        return 4.0 * self._side


@functools.singledispatch
def center(
        shape: Shape
) -> Point2D:
    """
    The 2D shape center as a point (2D)
    """

    raise NotImplementedError(f"Center method is not implemented for {type(shape)}")


@center.register(Point2D)
def _(
        shape: Point2D
) -> Point2D:

    return Point2D(
        x=shape.x,
        y=shape.y
    )


@center.register(Segment2D)
def _(
        shape: Segment2D
) -> Point2D:
    """
    Segment mean point.

    Examples:
      >>> center(Segment2D(Point2D(3, 2), Point2D(1, 1)))
      Point(2.0, 1.5)
    """

    x0, y0 = shape.start_pt.a()
    x1, y1 = shape.end_pt.a()

    return Point2D(
        x=(x0 + x1)/2.0,
        y=(y0 + y1)/2.0
    )


@center.register(Circle2D)
def _(
        shape: Circle2D
) -> Point2D:

    return Point2D(
        x=shape.x,
        y=shape.y
    )


def xytuple_list_to_line2d(
        xy_list: Tuple[numbers.Real, numbers.Real]
) -> Line2D:

    return Line2D(
        np.array([x for (x, _) in xy_list]),
        np.array([y for (_, y) in xy_list])
    )


def xytuple_l2_to_multiline2d(
        xytuple_list2
) -> MultiLine2D:

    # input is a list of list of (x,y) values

    #assert len(xytuple_list2) > 0
    lines_list = []
    for xy_list in xytuple_list2:
        #assert len(xy_list) > 0
        lines_list.append(xytuple_list_to_line2d(xy_list))

    return MultiLine2D(lines_list)


def merge_line2d(
        line
) -> Line2D:
    """
    line: a list of (x,y) tuples for line
    """

    line_type, line_geometry = line

    if line_type == 'multiline':
        path_line = xytuple_l2_to_multiline2d(line_geometry).to_line()
    elif line_type == 'line':
        path_line = xytuple_list_to_line2d(line_geometry)
    else:
        raise Exception("unknown line type")

    # transformed into a single Line

    result = path_line.remove_coincident_points()

    return result


def merge_lines(
        lines: List[Line2D],
        progress_ids
):
    """
    lines: a list of list of (x,y,z) tuples for multilines
    """

    sorted_line_list = [line for (_, line) in sorted(zip(progress_ids, lines))]

    line_list = []
    for line in sorted_line_list:

        line_type, line_geometry = line

        if line_type == 'multiline':
            path_line = xytuple_l2_to_multiline2d(line_geometry).to_line()
        elif line_type == 'line':
            path_line = xytuple_list_to_line2d(line_geometry)
        else:
            continue
        line_list.append(path_line)  # now a list of Lines

    # now the list of Lines is transformed into a single Line with coincident points removed

    line = MultiLine2D(line_list).to_line().remove_coincident_points()

    return line


def merge_lines2d(
        lines: List[Line2D],
        progress_ids
):
    """
    lines: a list of list of (x,y,z) tuples for multilines
    """

    sorted_line_list = [line for (_, line) in sorted(zip(progress_ids, lines))]

    line_list = []
    for line in sorted_line_list:

        line_type, line_geometry = line

        if line_type == 'multiline':
            path_line = xytuple_l2_to_multiline2d(line_geometry).to_line()
        elif line_type == 'line':
            path_line = xytuple_list_to_line2d(line_geometry)
        else:
            continue
        line_list.append(path_line)  # now a list of Lines

    # now the list of Lines is transformed into a single Line with coincident points removed

    line = MultiLine2D(line_list).to_line().remove_coincident_points()

    return line


class XYArrayPair:
    """
    Represent an x-y array pair (i.e., a single y value for each x value).
    X values should be sorted and should start from zero.
    """

    def __init__(self,
                 x_array: Union[array, np.ndarray],
                 y_array: Union[array, np.ndarray],
                 breaks: Optional[Union[array, np.ndarray]] = None
                 ):
        """
        Constructor.

        :param x_array: the x values array.
        :param y_array: the y values array.
        :param breaks: the internal x breaks.
        """

        check_type(x_array, "X array", (array, np.ndarray))

        if not np.all(np.isfinite(x_array)):
            raise Exception("X array values must all be finite")

        check_type(y_array, "Y array", (array, np.ndarray))
        if breaks is not None:
            check_type(breaks, "X breaks", (array, np.ndarray))

        if len(y_array) != len(x_array):
            raise Exception("Y array must have the same length as x array")

        if x_array[0] != 0.0:
            raise Exception("First value of X array should be zero")

        # from: https://stackoverflow.com/questions/47004506/check-if-a-numpy-array-is-sorted
        # answer by: luca
        is_sorted = lambda a: np.all(a[:-1] <= a[1:])

        if not is_sorted(x_array):
            raise Exception("X array must be already sorted")

        if breaks is not None:
            if breaks[0] != 0.0:
                raise Exception("First element of X breaks must always be zero")
            if not areClose(breaks[-1], x_array[-1]):
                raise Exception("Last element of X breaks must always be (almost) equal to last element of X array")

        self._x = np.array(x_array, dtype=np.float64)
        self._y = np.array(y_array, dtype=np.float64)
        self._x_breaks = np.array(breaks, dtype=np.float64) if breaks is not None else np.array([x_array[0], x_array[-1]])

    def clone(self) -> 'XYArrayPair':
        """
        Clone the XYArrayPair instance.
        """

        return XYArrayPair(
            x_array=np.copy(self._x),
            y_array=np.copy(self._y),
            breaks=None if self._x_breaks is None else np.copy(self._x_breaks)
        )

    def x_arr(self) -> np.ndarray:
        """
        Return the x array.

        :return: the x array.
        :rtype: array.
        """

        return self._x

    def y_arr(self) -> np.ndarray:
        """
        Return the y array.

        :return: the y array.
        """

        return self._y

    def x_breaks(self) -> Optional[np.ndarray]:

        return self._x_breaks

    def __repr__(self) -> str:
        """
        Representation of a topographic profile instance.

        :return: the textual representation of the instance.
        :rtype: str.
        """

        return f"XYArrayPair with {len(self.x_arr())} nodes\nx = {self._x},\ny = {self._y}"

    def num_steps(self) -> numbers.Integral:
        """
        Return the number of steps in the array pair.

        :return: number of steps in the array pair.
        """

        return len(self._x)

    def x(self,
          ndx: numbers.Integral
          ) -> numbers.Real:
        """
        Returns the x value with the index ndx.

        :param ndx: the index in the x array
        :return: the s value corresponding to the ndx index
        """

        return self._x[ndx]

    def y(self,
          ndx: numbers.Integral
          ) -> numbers.Real:
        """
        Returns the y value with the index ndx.

        :param ndx: the index in the y array
        :return: the y value corresponding to the ndx index
        """

        return self._y[ndx]

    def x_min(self) -> numbers.Real:
        """
        Returns the minimum x value.

        :return: the minimum x value.
        """

        return np.nanmin(self._x)

    def x_max(self) -> numbers.Real:
        """
        Returns the maximum x value.

        :return: the maximum x value.
        """

        return np.nanmax(self._x)

    def y_min(self) -> numbers.Real:
        """
        Returns the minimum y value.

        :return: the minimum y value.
        """

        return np.nanmin(self._y)

    def y_max(self) -> numbers.Real:
        """
        Returns the maximum y value.

        :return: the maximum y value.
        :rtype: numbers.Real.
        """

        return np.nanmax(self._y)

    def x_length(self) -> numbers.Real:
        """
        Returns the geographic length of the profile.

        :return: length of profile.
        :rtype: numbers.Real.
        """

        return self._x[-1]

    def find_index_ge(self,
      x_val: numbers.Real):
        """

        Examples:
          >>> p = XYArrayPair(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.find_index_ge(-1)
          0
          >>> p.find_index_ge(0.0)
          0
          >>> p.find_index_ge(0.5)
          1
          >>> p.find_index_ge(0.75)
          1
          >>> p.find_index_ge(1.0)
          1
          >>> p.find_index_ge(2.0)
          2
          >>> p.find_index_ge(2.5)
          3
          >>> p.find_index_ge(3.08)
          4
          >>> p.find_index_ge(3.14)
          4
          >>> p.find_index_ge(5) is None
          True
        """

        check_type(x_val, "X value", numbers.Real)
        if not np.isfinite(x_val):
            raise Exception(f"X value must be finite but {x_val} got")

        if x_val <= self.x_min():
            return 0
        elif x_val > self.x_max():
            return None
        else:
            return np.argmax(self._x >= x_val)

    def x_upper_ndx(self,
                    x_val: numbers.Real
                    ) -> Optional[numbers.Integral]:
        """
        To be possibly deprecated.
        Returns the optional index in the x array of the provided value.

        :param x_val: the value to search the index for in the x array
        :return: the optional index in the s array of the provided value

        Examples:
          >>> p = XYArrayPair(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.x_upper_ndx(-1) is None
          True
          >>> p.x_upper_ndx(5) is None
          True
          >>> p.x_upper_ndx(0.5)
          1
          >>> p.x_upper_ndx(0.75)
          1
          >>> p.x_upper_ndx(1.0)
          1
          >>> p.x_upper_ndx(2.0)
          2
          >>> p.x_upper_ndx(2.5)
          3
          >>> p.x_upper_ndx(3.11)
          4
          >>> p.x_upper_ndx(3.14)
          4
          >>> p.x_upper_ndx(0.0)
          0
        """

        check_type(x_val, "Input value", numbers.Real)

        if x_val < self.x_min() or x_val > self.x_max():
            return None

        return np.searchsorted(self._x, x_val)

    def y_linear_interpol(self,
                          x_val: numbers.Real
                          ) -> Optional[numbers.Real]:
        """
        Returns the optional interpolated z value in the z array of the provided s value.

        :param x_val: the value to search the index for in the s array
        :return: the optional interpolated z value

        Examples:
          >>> p = XYArrayPair(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.y_linear_interpol(-1) is None
          True
          >>> p.y_linear_interpol(5) is None
          True
          >>> p.y_linear_interpol(0.5)
          15.0
          >>> p.y_linear_interpol(0.75)
          17.5
          >>> p.y_linear_interpol(2.5)
          7.25
          >>> p.y_linear_interpol(3.14)
          17.9
          >>> p.y_linear_interpol(0.0)
          10.0
          >>> p.y_linear_interpol(1.0)
          20.0
        """

        check_type(x_val, "Input value", numbers.Real)

        ndx = self.x_upper_ndx(x_val)

        if ndx is not None:

            if ndx == 0:
                return self.y(0)

            val_y_i = self.y(ndx - 1)
            val_y_i_next = self.y(ndx)
            delta_val_y = val_y_i_next - val_y_i

            if delta_val_y == 0.0:
                return val_y_i

            val_x_i = self.x(ndx - 1)
            val_x_i_next = self.x(ndx)
            delta_val_x = val_x_i_next - val_x_i

            if delta_val_x == 0.0:
                return val_y_i

            d_x = x_val - val_x_i

            return val_y_i + d_x * delta_val_y / delta_val_x

        else:

            return None

    def x_subset(self,
                 x_start: numbers.Real,
                 x_end: Optional[numbers.Real] = None
                 ) -> np.ndarray:
        """
        Return the x array values defined by the provided x range (extremes included).
        When the end value is not provided, a single-valued array is returned.

        :param x_start: the start x value (distance along the profile)
        :param x_end: the optional end x value (distance along the profile)
        :return: the s array subset, possibly with just a value

        Examples:
          >>> p = XYArrayPair(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.x_subset(1.0)
          array([1.])
          >>> p.x_subset(0.0)
          array([0.])
          >>> p.x_subset(0.75)
          array([0.75])
          >>> p.x_subset(3.14)
          array([3.14])
          >>> p.x_subset(1.0, 2.0)
          array([1., 2.])
          >>> p.x_subset(0.75, 2.0)
          array([0.75, 1.  , 2.  ])
          >>> p.x_subset(0.75, 2.5)
          array([0.75, 1.  , 2.  , 2.5 ])
          >>> p.x_subset(0.75, 3.0)
          array([0.75, 1.  , 2.  , 3.  ])
          >>> p.x_subset(0.75, 0.5)
          NotImplemented
          >>> p.x_subset(-1, 1)
          NotImplemented
          >>> p.x_subset(-1)
          NotImplemented
          >>> p.x_subset(0.0, 10)
          NotImplemented
          >>> p.x_subset(0.0, 3.14)
          array([0.  ,  1.  ,  2.  ,  3.  , 3.14])
        """

        if not x_end and not self.x_min() <= x_start <= self.x_max():

            return NotImplemented

        if x_end and not self.x_min() <= x_start <= x_end <= self.x_max():

            return NotImplemented

        if x_end is None or x_end == x_start:

            return np.array([x_start], dtype=np.float64)

        values = []

        s_start_upper_index_value = self.x_upper_ndx(x_start)

        if x_start < self.x(s_start_upper_index_value):
            values.append(x_start)

        s_end_upper_index_value = self.x_upper_ndx(x_end)

        for ndx in range(s_start_upper_index_value, s_end_upper_index_value):
            values.append(self.x(ndx))

        if x_end > self.x(s_end_upper_index_value - 1):
            values.append(x_end)

        return np.array(values, dtype=np.float64)

    def ys_from_x_range(self,
                        x_start: numbers.Real,
                        x_end: Optional[numbers.Real] = None
                        ) -> np.ndarray:
        """
        Return the y array values defined by the provided x range (extremes included).
        When the end value is not provided, a single-valued array is returned.

        :param x_start: the start x value (distance along the profile)
        :param x_end: the optional end x value (distance along the profile)
        :return: the y array, possibly with just a value

        Examples:
          >>> p = XYArrayPair(array('d', [ 0.0,  1.0,  2.0,  3.0, 3.14]), array('d', [10.0, 20.0, 0.0, 14.5, 17.9]))
          >>> p.ys_from_x_range(1.0)
          array([20.])
          >>> p.ys_from_x_range(0.0)
          array([10.])
          >>> p.ys_from_x_range(0.75)
          array([17.5])
          >>> p.ys_from_x_range(3.14)
          array([17.9])
          >>> p.ys_from_x_range(1.0, 2.0)
          array([20., 0.])
          >>> p.ys_from_x_range(0.75, 2.0)
          array([17.5, 20. , 0. ])
          >>> p.ys_from_x_range(0.75, 2.5)
          array([17.5 , 20.  , 0.  , 7.25])
          >>> p.ys_from_x_range(0.75, 3.0)
          array([17.5, 20. , 0. , 14.5])
          >>> p.ys_from_x_range(0.75, 0.5)
          NotImplemented
          >>> p.ys_from_x_range(-1, 1)
          NotImplemented
          >>> p.ys_from_x_range(-1)
          NotImplemented
          >>> p.ys_from_x_range(0.0, 10)
          NotImplemented
          >>> p.ys_from_x_range(0.0, 3.14)
          array([10. , 20. , 0. , 14.5, 17.9])
        """

        s_subset = self.x_subset(x_start, x_end)

        if s_subset is NotImplemented:
            return NotImplemented

        return np.array(list(map(self.y_linear_interpol, s_subset)), dtype=np.float64)

    def extend_in_place(self,
               another: 'XYArrayPair'
               ) -> 'XYArrayPair':
        """
        Extend an XYArrayPair instance in-place.
        Note that the last element of the first x-y array pair will be dropped
        and substituted by the first element of the second x-y array pair
        (no attempt to check values equality or to calculate a mean).
        Moverover, all the x values of the second x-y array pair will be incremented
        by the last x value of the first element.

        Examples:
          >>> f = XYArrayPair(np.array([ 0.0,  14.2,  20.0]), np.array([149.17, 132.4, 159.2]))
          >>> f.x_breaks()
          array([ 0., 20.])
          >>> f.x_max()
          20.0
          >>> s = XYArrayPair(np.array([ 0.0,  7.0,  11.0]), np.array([159.17, 180.1, 199.5]))
          >>> s.x_breaks()
          array([ 0., 11.])
          >>> s._x_breaks
          array([ 0., 11.])
          >>> s._x_breaks + f.x_max()
          array([20., 31.])
          >>> f.extend_in_place(s)
          >>> f.x_arr()
          array([ 0. ,  14.2,  20. , 27. , 31. ])
          >>> f.y_arr()
          array([149.17, 132.4 , 159.17, 180.1 , 199.5 ])
          >>> f.x_breaks()
          array([ 0., 20., 31.])
        """

        check_type(another, "Second XYArrayPair instance", XYArrayPair)

        offset = self.x_max()
        self._x = np.append(self._x[:-1], another._x + offset)
        self._y = np.append(self._y[:-1], another._y)
        self._x_breaks = np.append(self._x_breaks[:-1], another._x_breaks + offset)

    def dir_slopes_deltas(self) -> Tuple[np.ndarray, np.ndarray]:

        delta_x = np.ediff1d(self._x, to_end=np.nan)
        delta_y = np.ediff1d(self._y, to_end=np.nan)

        return delta_x, delta_y

    def dir_slopes_ratios(self) -> np.ndarray:

        dx, dy = self.dir_slopes_deltas()

        return dy / dx

    def dir_slopes_percent(self) -> np.ndarray:

        return 100.0 * self.dir_slopes_ratios()

    def dir_slopes_radians(self) -> np.ndarray:

        dx, dy = self.dir_slopes_deltas()

        return np.arctan2(dy, dx)

    def dir_slopes_degrees(self) -> np.ndarray:

        return self.dir_slopes_radians() * 180.0 / np.pi

    def as_dir_slopes_degrees(self) -> 'XYArrayPair':

        return XYArrayPair(
            x_array=self._x,
            y_array=self.dir_slopes_degrees(),
            breaks=self._x_breaks
        )

    def as_abs_slopes_degrees(self) -> 'XYArrayPair':

        return XYArrayPair(
            x_array=self._x,
            y_array=np.fabs(self.dir_slopes_degrees()),
            breaks=self._x_breaks
        )


def combine_xy_arrays(
        *arrays: Tuple[XYArrayPair]
) -> Optional[XYArrayPair]:
    """
    Combines a tuple of XYArrayPair instances into a single XYArrayPair instance.

    Examples:
      >>> f = XYArrayPair(np.array([ 0.0,  14.2,  20.0]), np.array([149.17, 132.4, 159.2]))
      >>> s = XYArrayPair(np.array([ 0.0,  7.0,  11.0]), np.array([159.17, 180.1, 199.5]))
      >>> t = XYArrayPair(np.array([ 0.0,  22.0,  30.0]), np.array([199.5, 200.1, 179.1]))
      >>> combined = combine_xy_arrays(f, s, t)
      >>> combined.x_arr()
      array([ 0. ,  14.2, 20. , 27. , 31. , 53. , 61. ])
      >>> combined.y_arr()
      array([149.17, 132.4 , 159.17, 180.1 , 199.5 , 200.1 , 179.1 ])
      >>> combined.x_breaks()
      array([ 0., 20., 31., 61.])
    """

    if len(arrays) == 0:
        return None

    combined_xy_arrays = arrays[0].clone()

    for xy_array_pair in arrays[1:]:
        combined_xy_arrays.extend_in_place(xy_array_pair)

    return combined_xy_arrays


class SimpleVerticalTrapezion(Polygon):
    """
    It's a trapezion with a horizontal base at height y0,
    two vertical side at x0 and x1,
    a corner at one side at y1
    and the extreme vertex at the other side at y2

      .   y2
     /|
    / |   y1
    | |
    | |
    ___   y0


    """

    def __init__(self,
                 x0: numbers.Real,
                 x1: numbers.Real,
                 y0: numbers.Real,
                 y1: numbers.Real,
                 y2: numbers.Real
                 ):

        for val in (x0, x1, y0, y1, y2):
            check_type(val, f"Input {val} is not a real", numbers.Real)

        if not (y0 <= y1 <= y2 or y0 >= y1 >= y2):
            raise Exception("Y values not fully increasing or decreasing")

        super(SimpleVerticalTrapezion, self).__init__()

        self._x0 = x0
        self._x1 = x1
        self._y0 = y0
        self._y1 = y1
        self._y2 = y2

    def num_side(self):
        """Return number of sides"""

        return 4

    def area(self) -> numbers.Real:
        """
        The area.

        Examples:
          >>> SimpleVerticalTrapezion(0, 2, 0, 3, 5).area()
          8.0
          >>> SimpleVerticalTrapezion(2, 0, 0, 3, 5).area()
          8.0
          >>> SimpleVerticalTrapezion(2, 0, 0, -3, -5).area()
          8.0
          >>> SimpleVerticalTrapezion(-2, -4, 2, -1, -3).area()
          8.0
          >>> SimpleVerticalTrapezion(-2, -4, -3, -6, -8).area()
          8.0
          >>> SimpleVerticalTrapezion(-2, -1, -3, -6, -8).area()
          4.0
        """

        base = abs(self._x1 - self._x0)
        rectangle_height = abs(self._y1 - self._y0)
        triangle_height = abs(self._y2 - self._y1)

        rectangle_area = base * rectangle_height
        triangle_area = (base * triangle_height) / 2.0

        return rectangle_area + triangle_area

    def length(self):

        return NotImplemented

    def clone(self):

        return SimpleVerticalTrapezion(
            self._x0,
            self._x1,
            self._y0,
            self._y1,
            self._y2
        )


if __name__ == "__main__":

    import doctest
    doctest.testmod()

