# -*- coding: utf-8 -*-


from typing import Dict, Union

from math import sqrt
from enum import Enum
from array import array
import itertools

from ...types.utils import *
from ...mathematics.defaults import MIN_SEPARATION_THRESHOLD
from ...mathematics.statistics import get_statistics

from ..vectors import *
from ..projections.crs import Crs

from ...utils.lists import find_val


class Point(object):
    """
    Cartesian point.
    Dimensions: 4D (space-time)
    """

    def __init__(
        self,
        x: numbers.Real,
        y: numbers.Real,
        z: numbers.Real = 0.0,
        t: numbers.Real = 0.0,
        epsg_cd: numbers.Integral = -1):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :type x: numbers.Real.
        :param y: point y coordinate.
        :type y: numbers.Real.
        :param z: point z coordinate.
        :type z: numbers.Real.
        :param t: point time coordinate.
        :type t: numbers.Real.
        :param epsg_cd: CRS EPSG code.
        :type epsg_cd: numbers.Integral.
        """

        vals = [x, y, z, t]
        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise VectorInputException("Input values must be integer or float type")
        if not all(map(math.isfinite, vals)):
            raise VectorInputException("Input values must be finite")

        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
        self._t = float(t)
        self._crs = Crs(epsg_cd)

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7, epsg_cd=4326).x
          4.0
          >>> Point(-0.39, 3, 7).x
          -0.39
        """

        return self._x

    @property
    def y(self) -> numbers.Real:
        """
        Return the y coordinate of the current point.

        :return: y coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7, epsg_cd=4326).y
          3.0
          >>> Point(-0.39, 17.42, 7).y
          17.42
        """

        return self._y

    @property
    def z(self) -> numbers.Real:
        """
        Return the z coordinate of the current point.

        :return: z coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7, epsg_cd=4326).z
          7.0
          >>> Point(-0.39, 17.42, 8.9).z
          8.9
        """

        return self._z

    @property
    def t(self) -> numbers.Real:
        """
        Return the time coordinate of the current point.

        :return: time coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7, epsg_cd=4326).t
          0.0
          >>> Point(-0.39, 17.42, 8.9, 4112).t
          4112.0
        """

        return self._t

    def crs(self) -> Crs:
        """
        Return a copy of the current point CRS.

        :return: the CRS code as a Csr instance.
        :rtype: Crs.
        """

        return Crs(self._crs.epsg())

    def epsg(self) -> numbers.Integral:
        """
        Returns the EPSG code of the point.

        :return: the CRS code.
        :rtype: numbers.Integral.
        """

        return self.crs().epsg()

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f}, {:.4f}, {})".format(self.x, self.y, self.z, self.t, self.epsg())

    def __eq__(self, another: 'Point') -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1., epsg_cd=4326) == Point(1, 1, 1, epsg_cd=4326)
          True
          >>> Point(1., 1., 1., epsg_cd=4326) == Point(1, 1, -1, epsg_cd=4326)
          False
        """

        if not isinstance(another, Point):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            self.z == another.z,
            self.t == another.t,
            self.crs() == another.crs()])

    def __ne__(self, another: 'Point') -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
          >>> Point(1., 1., 1., epsg_cd=4326) != Point(1, 1, 1)
          True
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Integral]:
        """
        Return the individual values of the point.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7, epsg_cd=4326).a()
          (4.0, 3.0, 7.0, 0.0, 4326)
        """

        return self.x, self.y, self.z, self.t, self.epsg()

    def clone(self) -> 'Point':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return Point(*self.a())

    def toXYZ(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toXYZT(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial and time components as a tuple of four values.

        :return: the spatial components (x, y, z) and the time component.
        :rtype: a tuple of four floats.

        Examples:
          >>> Point(1, 0, 3).toXYZT()
          (1.0, 0.0, 3.0, 0.0)
        """

        return self.x, self.y, self.z, self.t

    def toArray(self) -> 'np.array':
        """
        Return a Numpy array representing the point values (without the crs code).

        :return: Numpy array

        Examples:
          >>> np.allclose(Point(1, 2, 3).toArray(), np.array([ 1., 2., 3., 0.]))
          True
        """

        return np.asarray(self.toXYZT())

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000, 0.0000, -1)
        """

        return Point(self.x, self.y, 0.0, self.t, self.epsg())

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000, 0.0000, -1)
        """

        return Point(self.x, 0.0, self.z, self.t, self.epsg())

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000, 0.0000, -1)
        """

        return Point(0.0, self.y, self.z, self.t, self.epsg())

    def deltaX(self, another: 'Point') -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3, epsg_cd=32632).deltaX(Point(4, 7, 1, epsg_cd=32632))
          3.0
          >>> Point(1, 2, 3, epsg_cd=4326).deltaX(Point(4, 7, 1)) is None
          True
        """

        if self.crs() != another.crs():
            return None
        else:
            return another.x - self.x

    def deltaY(self, another: 'Point') -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3, epsg_cd=32632).deltaY(Point(4, 7, 1, epsg_cd=32632))
          5.0
          >>> Point(1, 2, 3, epsg_cd=4326).deltaY(Point(4, 7, 1)) is None
          True
        """

        if self.crs() != another.crs():
            return None
        else:
            return another.y - self.y

    def deltaZ(self, another: 'Point') -> Optional[numbers.Real]:
        """
        Delta between z components of two Point Instances.

        :return: z coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3, epsg_cd=32632).deltaZ(Point(4, 7, 1, epsg_cd=32632))
          -2.0
          >>> Point(1, 2, 3, epsg_cd=4326).deltaZ(Point(4, 7, 1)) is None
          True
        """

        if self.crs() != another.crs():
            return None
        else:
            return another.z - self.z

    def deltaT(self, another: 'Point') -> numbers.Real:
        """
        Delta between t components of two Point Instances.

        :return: difference value
        :rtype: numbers.Real

        Examples:
          >>> Point(1, 2, 3, 17.3).deltaT(Point(4, 7, 1, 42.9))
          25.599999999999998
        """

        return another.t - self.t

    def dist3DWith(self, another: 'Point') -> numbers.Real:
        """
        Calculate Euclidean spatial distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point(1., 1., 1., epsg_cd=32632).dist3DWith(Point(4., 5., 1, epsg_cd=32632))
          5.0
          >>> Point(1, 1, 1, epsg_cd=32632).dist3DWith(Point(4, 5, 1, epsg_cd=32632))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        check_type(another, "Point", Point)

        if self.crs() != another.crs():
            raise Exception("Input point should have {} EPSG code but has {}".format(self.crs(), another.crs()))

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist2DWith(self, another: 'Point') -> numbers.Real:
        """
        Calculate horizontal (2D) distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the 2D distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point(1., 1., 1., epsg_cd=32632).dist2DWith(Point(4., 5., 7., epsg_cd=32632))
          5.0
          >>> Point(1., 1., 1., epsg_cd=32632).dist2DWith(Point(4., 5., 7.))
          Traceback (most recent call last):
          ...
          Exception: Second point should have 32632 EPSG code but has -1
        """

        check_type(another, "Second point", Point)

        if self.crs() != another.crs():
            raise Exception("Second point should have {} EPSG code but has {}".format(self.epsg(), another.epsg()))

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self, scale_factor: numbers.Real) -> 'Point':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000, 0.0000, -1)
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000, 0.0000, -1)
        """

        x, y, z = self.x * scale_factor, self.y * scale_factor, self.z * scale_factor
        return Point(x, y, z, self.t, self.epsg())

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000, 0.0000, -1)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000, 0.0000, -1)
        """

        return self.scale(-1)

    def isCoinc(self,
            another: 'Point',
            tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Check spatial coincidence of two points

        :param another: the point to compare.
        :type another: Point.
        :param tolerance: the maximum allowed distance between the two points.
        :type tolerance: numbers.Real.
        :return: whether the two points are coincident.
        :rtype: bool.
        :raise: Exception.

        Example:
          >>> Point(1., 0., -1.).isCoinc(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0., epsg_cd=32632).isCoinc(Point(1., 0., 0., epsg_cd=32632))
          True
          >>> Point(1.2, 7.4, 1.4, epsg_cd=32632).isCoinc(Point(1.2, 7.4, 1.4))
          Traceback (most recent call last):
          ...
          Exception: The first point has 32632 EPSG code while the second has -1
          >>> Point(1.2, 7.4, 1.4, epsg_cd=4326).isCoinc(Point(1.2, 7.4, 1.4))
          Traceback (most recent call last):
          ...
          Exception: The first point has 4326 EPSG code while the second has -1
        """

        check_type(another, "Second point", Point)

        if self.crs() != another.crs():
            raise Exception("The first point has {} EPSG code while the second has {}".format(self.epsg(), another.epsg()))

        return self.dist3DWith(another) <= tolerance

    def already_present(self, pt_list: List['Point'], tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD) -> Optional[bool]:
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
            if self.isCoinc(pt, tolerance=tolerance):
                return True
        return False

    def shift(self, sx: numbers.Real, sy: numbers.Real, sz: numbers.Real) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1, epsg_cd=32632).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000, 0.0000, 32632)
          >>> Point(1, 2, -1, epsg_cd=32632).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000, 0.0000, 32632)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz, self.t, self.epsg())

    def shiftByVect(self, v: Vect) -> 'Point':
        """
        Create a new point shifted from the self instance by given vector.

        :param v: the shift vector.
        :type v: Vect.
        :return: the shifted point.
        :rtype: Point.

        Example:
          >>> Point(1, 1, 1, epsg_cd=32632).shiftByVect(Vect(0.5, 1., 1.5, epsg_cd=32632))
          Point(1.5000, 2.0000, 2.5000, 0.0000, 32632)
          >>> Point(1, 2, -1, epsg_cd=32632).shiftByVect(Vect(0.5, 1., 1.5, epsg_cd=32632))
          Point(1.5000, 3.0000, 0.5000, 0.0000, 32632)
       """

        if self.crs() != v.crs():
            return None

        sx, sy, sz = v.toXYZ()

        return Point(self.x + sx, self.y + sy, self.z + sz, self.t, self.epsg())

    def asVect(self) -> 'Vect':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0).asVect()
          Vect(1.0000, 1.0000, 0.0000, EPSG: -1)
          >>> Point(0.2, 1, 6).asVect()
          Vect(0.2000, 1.0000, 6.0000, EPSG: -1)
        """

        return Vect(self.x, self.y, self.z, self.epsg())


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: CPlane is locational - its position in space is defined.
    This contrast with Plane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a: numbers.Real, b: numbers.Real, c: numbers.Real, d: numbers.Real, epsg_cd: numbers.Integral = -1):

        if not isinstance(a, numbers.Real):
            raise Exception("Input value a must be float or int but is {}".format(type(a)))
        if not isinstance(b, numbers.Real):
            raise Exception("Input value b must be float or int but is {}".format(type(b)))
        if not isinstance(c, numbers.Real):
            raise Exception("Input value c must be float or int but is {}".format(type(c)))
        if not isinstance(d, numbers.Real):
            raise Exception("Input value d must be float or int but is {}".format(type(d)))
        if not isinstance(epsg_cd, numbers.Integral):
            raise Exception("Input value epsg_cd must be int but is {}".format(type(epsg_cd)))

        norm = sqrt(a*a + b*b + c*c)
        self._a = float(a) / norm
        self._b = float(b) / norm
        self._c = float(c) / norm
        self._d = float(d) / norm
        self._crs = Crs(epsg_cd)

    def a(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).a()
          1.0
        """

        return self._a

    def b(self) -> numbers.Real:
        """
        Return b coefficient of a CPlane instance.

        Example:
          >>> CPlane(1, 4, 0, 2).b()
          0.9701425001453319
        """

        return self._b

    def c(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane instance.

        Example:
          >>> CPlane(1, 0, 5.4, 2).c()
          0.9832820049844602
        """

        return self._c

    def d(self) -> numbers.Real:
        """
        Return a coefficient of a CPlane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).d()
          2.0
        """

        return self._d

    def crs(self) -> Crs:
        """
        Returns the EPSG code as a Crs instance.

        :return: EPSG code.
        :rtype: pygsf.projections.crs.Crs

        Example:
        """

        return self._crs

    def epsg(self) -> numbers.Integral:
        """
        Returns the EPSG code.

        :return: EPSG code.
        :rtype: numbers.Integral

        Example:
        """

        return self._crs.epsg()

    def v(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real, numbers.Integral]:
        """
        Return coefficients of a CPlane instance.

        Example:
          >>> CPlane(1, 1, 7, -4).v()
          (0.14002800840280097, 0.14002800840280097, 0.9801960588196068, -0.5601120336112039, -1)
        """

        return self.a(), self.b(), self.c(), self.d(), self.epsg()

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3) -> 'CPlane':
        """
        Create a CPlane from three given Point instances.

        Example:
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000, -1)
          >>> CPlane.fromPoints(Point(0, 0, 0, epsg_cd=4326), Point(1, 0, 0, epsg_cd=4326), Point(0, 1, 0, epsg_cd=4326))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000, 4326)
          >>> CPlane.fromPoints(Point(0, 0, 0, epsg_cd=4326), Point(0, 1, 0, epsg_cd=4326), Point(0, 0, 1, epsg_cd=4326))
          CPlane(1.0000, 0.0000, 0.0000, 0.0000, 4326)
          >>> CPlane.fromPoints(Point(1,2,3), Point(2,3,4), Point(-1,7,-2))
          CPlane(-0.7956, 0.2387, 0.5569, -1.3524, -1)
        """

        if not (isinstance(pt1, Point)):
            raise Exception("First input point should be Point but is {}".format(type(pt1)))

        if not (isinstance(pt2, Point)):
            raise Exception("Second input point should be Point but is {}".format(type(pt2)))

        if not (isinstance(pt3, Point)):
            raise Exception("Third input point should be Point but is {}".format(type(pt3)))

        if pt2.crs() != pt1.crs():
            raise Exception("Second input point should have {} EPSG code like first point but has {}".format(pt1.epsg(), pt2.epsg()))

        if pt3.crs() != pt1.crs():
            raise Exception("Third input point should have {} EPSG code like first point but has {}".format(pt1.epsg(), pt3.epsg()))

        matr_a = np.array(
            [[pt1.y, pt1.z, 1],
             [pt2.y, pt2.z, 1],
             [pt3.y, pt3.z, 1]])

        matr_b = - np.array(
            [[pt1.x, pt1.z, 1],
             [pt2.x, pt2.z, 1],
             [pt3.x, pt3.z, 1]])

        matr_c = np.array(
            [[pt1.x, pt1.y, 1],
             [pt2.x, pt2.y, 1],
             [pt3.x, pt3.y, 1]])

        matr_d = - np.array(
            [[pt1.x, pt1.y, pt1.z],
             [pt2.x, pt2.y, pt2.z],
             [pt3.x, pt3.y, pt3.z]])

        return cls(
            np.linalg.det(matr_a),
            np.linalg.det(matr_b),
            np.linalg.det(matr_c),
            np.linalg.det(matr_d),
            epsg_cd=pt1.epsg())

    def __repr__(self):

        return "CPlane({:.4f}, {:.4f}, {:.4f}, {:.4f}, {:d})".format(*self.v(), self.crs())

    def normVersor(self) -> Vect:
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane(0, 0, 5, -2).normVersor()
          Vect(0.0000, 0.0000, 1.0000, EPSG: -1)
          >>> CPlane(0, 7, 0, 5, epsg_cd=32632).normVersor()
          Vect(0.0000, 1.0000, 0.0000, EPSG: 32632)
        """

        return Vect(self.a(), self.b(), self.c(), epsg_cd=self.epsg()).versor()

    def toPoint(self) -> Point:
        """
        Returns a point lying in the plane (non-unique solution).

        Examples:
          >>> CPlane(0, 0, 1, -1).toPoint()
          Point(0.0000, 0.0000, 1.0000, 0.0000, -1)
        """

        point = Point(
            *pointSolution(
                np.array([[self.a(), self.b(), self.c()]]),
                np.array([-self.d()])),
            epsg_cd=self.epsg())

        return point

    def intersVersor(self, another) -> Optional[Vect]:
        """
        Return intersection versor for two intersecting planes.

        :param another: another Cartesian plane.
        :type another: CPlane.
        :return: the intersection line as a vector.
        :rtype: Vect.
        :raise: Exception.

        Examples:
          >>> a = CPlane(1, 0, 0, 0, epsg_cd=2000)
          >>> b = CPlane(0, 0, 1, 0, epsg_cd=2000)
          >>> a.intersVersor(b)
          Vect(0.0000, -1.0000, 0.0000, EPSG: 2000)
        """

        if not isinstance(another, CPlane):
            raise Exception("The argument must be a Cartesian plane")

        if self.crs() != another.crs():
            return None

        return self.normVersor().vCross(another.normVersor()).versor()

    def intersPoint(self, another) -> Optional[Point]:
        """
        Return point on intersection line (non-unique solution)
        for two planes.

        Examples:
          >>> a = CPlane(1, 0, 0, 0, epsg_cd=32632)
          >>> b = CPlane(0, 0, 1, 0, epsg_cd=32632)
          >>> a.intersPoint(b)
          Point(0.0000, 0.0000, 0.0000, 0.0000, 32632)
        """

        if self.crs() != another.crs():
            return None

        # find a point lying on the intersection line (this is a non-unique solution)

        a = np.array([[self.a(), self.b(), self.c()], [another.a(), another.b(), another.c()]])
        b = np.array([-self.d(), -another.d()])
        x, y, z = pointSolution(a, b)

        return Point(x, y, z, epsg_cd=self.epsg())

    def pointDistance(self, pt: Point) -> numbers.Real:
        """
        Calculate the distance between a point and the cartesian plane.
        Distance expression:
        distance = a * x1 + b * y1 + c * z1 + d
        where a, b, c and d are plane parameters of the plane equation:
         a * x + b * y + c * z + d = 0
        and x1, y1, and z1 are the point coordinates.

        :param pt: the point to calculate distance with.
        :type pt: Point.
        :return: the distance value.
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> cpl = CPlane(0, 0, 1, 0, epsg_cd=32632)
          >>> pt = Point(0, 0, 1, epsg_cd=32632)
          >>> cpl.pointDistance(pt)
          1.0
          >>> pt = Point(0, 0, 0.5, epsg_cd=32632)
          >>> cpl.pointDistance(pt)
          0.5
          >>> pt = Point(0, 0, -0.5, epsg_cd=32632)
          >>> cpl.pointDistance(pt)
          -0.5
          >>> pt = Point(10, 20, 0.0, epsg_cd=32632)
          >>> cpl.pointDistance(pt)
          0.0
          >>> cpl = CPlane(0, 0, 1, 0, epsg_cd=32632)
          >>> pt = Point(10, 20, 0.0)
          >>> cpl.pointDistance(pt)
          Traceback (most recent call last):
          ...
          Exception: Input point should have 32632 EPSG code but has -1
        """

        if not (isinstance(pt, Point)):
            raise Exception("Input point should be Point but is {}".format(type(pt)))

        if self.crs() != pt.crs():
            raise Exception("Input point should have {} EPSG code but has {}".format(self.epsg(), pt.epsg()))

        return self.a() * pt.x + self.b() * pt.y + self.c() * pt.z + self.d()

    def isPointInPlane(self, pt) -> bool:
        """
        Check whether a point lies in the current plane.

        :param pt: the point to check.
        :type pt: Point.
        :return: whether the point lies in the current plane.
        :rtype: bool.
        :raise: Exception.

        Examples:
          >>> pl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
          >>> pl = CPlane(0, 0, 1, 0, epsg_cd=32632)
          >>> pt = Point(0, 1, 0, epsg_cd=32632)
          >>> pl.isPointInPlane(pt)
          True
        """

        if not (isinstance(pt, Point)):
            raise Exception("Input point should be Point but is {}".format(type(pt)))

        if self.crs() != pt.crs():
            raise Exception("Input point should have {} EPSG code like the plane but has {}".format(self.epsg(), pt.epsg()))

        if abs(self.pointDistance(pt)) < MIN_SEPARATION_THRESHOLD:
            return True
        else:
            return False

    def angle(self, another) -> numbers.Real:
        """
        Calculate angle (in degrees) between two planes.

        :param another: the CPlane instance to calculate angle with.
        :type another: CPlane.
        :return: the angle (in degrees) between the two planes.
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> CPlane(1,0,0,0).angle(CPlane(0,1,0,0))
          90.0
          >>> CPlane(1,0,0,0, epsg_cd=32632).angle(CPlane(0,1,0,0, epsg_cd=32632))
          90.0
          >>> CPlane(1,0,0,0, epsg_cd=32632).angle(CPlane(1,0,1,0, epsg_cd=32632))
          45.0
          >>> CPlane(1,0,0,0, epsg_cd=32632).angle(CPlane(1,0,0,0, epsg_cd=32632))
          0.0
          >>> CPlane(1,0,0,0, epsg_cd=32632).angle(CPlane(1,0,0,0))
          Traceback (most recent call last):
          ...
          Exception: Input plane should have 32632 EPSG code like current point but has -1
        """

        if not (isinstance(another, CPlane)):
            raise Exception("Input plane should be CPlane but is {}".format(type(another)))

        if self.crs() != another.crs():
            raise Exception("Input plane should have {} EPSG code like current point but has {}".format(self.epsg(), another.epsg()))

        angle_degr = self.normVersor().angle(another.normVersor())

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


class JoinTypes(Enum):
    """
    Enumeration for Line and Segment type.
    """

    START_START = 1  # start point coincident with start point
    START_END   = 2  # start point coincident with end point
    END_START   = 3  # end point coincident with start point
    END_END     = 4  # end point coincident with end point


class Segment(object):
    """
    Segment is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self, start_pt: Point, end_pt: Point):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :type: Point.
        :param end_pt: the end point.
        :type end_pt: Point.
        :return: the new segment instance if both points have the same crs.
        :raises: CRSCodeException.
        """

        if not isinstance(start_pt, Point):
            raise Exception("Start point must be a Point instance")

        if not isinstance(end_pt, Point):
            raise Exception("Start point must be a Point instance")

        if start_pt.dist3DWith(end_pt) == 0.0:
            raise Exception("Segment point distance must be greater than zero")

        if start_pt.crs() != end_pt.crs():
            raise Exception("Start and end point must have the same CRS code")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()
        self._crs = Crs(start_pt.epsg())

    def __repr__(self) -> str:
        """
        Represents a Segment instance.

        :return: the Segment representation.
        :rtype: str.
        """

        return "Segment(start_pt={}, end_pt={})".format(
            self.start_pt(),
            self.end_pt()
        )

    def crs(self) -> Crs:

        return Crs(self._crs.epsg())

    def epsg(self) -> numbers.Integral:

        return self.crs().epsg()

    def extract_start_pt(self) -> Point:

        return self._start_pt

    def extract_end_pt(self) -> Point:

        return self._end_pt

    def start_pt(self) -> Point:

        return self.extract_start_pt().clone()

    def end_pt(self) -> Point:

        return self.extract_end_pt().clone()

    def clone(self) -> 'Segment':

        return Segment(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment':

        if self.end_pt().x < self.start_pt().x:
            return Segment(self.end_pt(), self.start_pt())
        else:
            return self.clone()

    def x_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt().x < self.end_pt().x:
            return self.start_pt().x, self.end_pt().x
        else:
            return self.end_pt().x, self.start_pt().x

    def y_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt().y < self.end_pt().y:
            return self.start_pt().y, self.end_pt().y
        else:
            return self.end_pt().y, self.start_pt().y

    def z_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt().z < self.end_pt().z:
            return self.start_pt().z, self.end_pt().z
        else:
            return self.end_pt().z, self.start_pt().z

    def delta_x(self) -> numbers.Real:

        return self.end_pt().x - self.start_pt().x

    def delta_y(self) -> numbers.Real:

        return self.end_pt().y - self.start_pt().y

    def delta_z(self) -> numbers.Real:
        """
        Z delta between segment end point and start point.

        :return: numbers.Real.
        """

        return self.end_pt().z - self.start_pt().z

    def length2D(self) -> numbers.Real:
        """
        Returns the horizontal length of the segment.

        :return: the horizontal length of the segment.
        :rtype: numbers.Real.
        """

        return self.start_pt().dist2DWith(self.end_pt())

    def length3D(self) -> numbers.Real:

        return self.start_pt().dist3DWith(self.end_pt())

    def deltaZS(self) -> Optional[numbers.Real]:
        """
        Calculates the delta z - delta s ratio of a segment.

        :return: optional numbers.Real.
        """

        len2d = self.length2D()

        if len2d == 0.0:
            return None

        return self.delta_z() / len2d

    def slope_rad(self) -> Optional[numbers.Real]:
        """
        Calculates the slope in radians of the segment.
        Positive is downward point, negative upward pointing.

        :return: optional numbers.Real.
        """

        delta_zs = self.deltaZS()

        if delta_zs is None:
            return None
        else:
            return - math.atan(delta_zs)

    def vector(self) -> Vect:

        return Vect(self.delta_x(),
                    self.delta_y(),
                    self.delta_z(),
                    epsg_cd=self.epsg())

    def antivector(self) -> Vect:
        """
        Returns the vector pointing from the segment end to the segment start.

        :return: the vector pointing from the segment end to the segment start.
        :rtype: Vect.
        """

        return self.vector().invert()

    def segment_2d_m(self) -> Optional[numbers.Real]:

        denom = self.end_pt().x - self.start_pt().x

        if denom == 0.0:
            return None

        return (self.end_pt().y - self.start_pt().y) / denom

    def segment_2d_p(self) -> Optional[numbers.Real]:

        s2d_m = self.segment_2d_m()

        if s2d_m is None:
            return None

        return self.start_pt().y - s2d_m * self.start_pt().x

    def intersection_2d_pt(self, another) -> Optional[Point]:

        if self.crs() != another.crs():
            return None

        s_len2d = self.length2D()
        a_len2d = another.length2D()

        if s_len2d == 0.0 or a_len2d == 0.0:
            return None

        if self.start_pt().x == self.end_pt().x:  # self segment parallel to y axis
            x0 = self.start_pt().x
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            if m1 is None:
                return None
            y0 = m1 * x0 + p1
        elif another.start_pt().x == another.end_pt().x:  # another segment parallel to y axis
            x0 = another.start_pt().x
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

        return Point(x0, y0, epsg_cd=self.epsg())

    def contains_pt(self,
        pt: Point) -> bool:
        """
        Checks whether a point is contained in a segment.
        It is contained, the angle between
        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment(Point(0, 0, 0), Point(1, 0, 0),)
          >>> segment.contains_pt(Point(0, 0, 0))
          True
          >>> segment.contains_pt(Point(1, 0, 0))
          True
          >>> segment.contains_pt(Point(0.5, 0, 0))
          True
          >>> segment.contains_pt(Point(0.5, 0.00001, 0))
          False
          >>> segment.contains_pt(Point(0.5, 0, 0.00001))
          False
          >>> segment.contains_pt(Point(1.00001, 0, 0))
          False
          >>> segment.contains_pt(Point(0.000001, 0, 0))
          True
          >>> segment.contains_pt(Point(-0.000001, 0, 0))
          False
          >>> segment.contains_pt(Point(0.5, 1000, 1000))
          False
        """

        check_type(pt, "Point", Point)

        segment_length = self.length3D()
        length_startpt_pt = self.start_pt().dist3DWith(pt)
        length_endpt_pt = self.end_pt().dist3DWith(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def fast_2d_contains_pt(self, pt2d) -> bool:
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

    def scale(self, scale_factor) -> 'Segment':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: numbers.Real
        :return: Segment instance
        """

        delta_x = self.delta_x() * scale_factor
        delta_y = self.delta_y() * scale_factor
        delta_z = self.delta_z() * scale_factor

        end_pt = Point(
            x=self.start_pt().x + delta_x,
            y=self.start_pt().y + delta_y,
            z=self.start_pt().z + delta_z,
            t=self.end_pt().t,
            epsg_cd=self.epsg())

        return Segment(
            self.start_pt(),
            end_pt)

    def densify2d_asSteps(self, densify_distance: numbers.Real) -> array:
        """
        Defines the array storing the incremental lengths according to the provided densify distance.

        :param densify_distance: the step distance.
        :type densify_distance: numbers.Real.
        :return: array storing incremental steps, with the last step being equal to the segment length.
        :rtype: array.
        """

        if not isinstance(densify_distance, numbers.Real):
            raise Exception("Densify distance must be float or int")

        if not math.isfinite(densify_distance):
            raise Exception("Densify distance must be finite")

        if not densify_distance > 0.0:
            raise Exception("Densify distance must be positive")

        segment_length = self.length2D()

        s_list = []
        n = 0
        length = n * densify_distance

        while length < segment_length:
            s_list.append(length)
            n += 1
            length = n * densify_distance

        s_list.append(segment_length)

        return array('d', s_list)

    def densify2d_asPts(self, densify_distance) -> List[Point]:
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: the distance with which to densify the segment.
        :type densify_distance: numbers.Real.
        :return: the set of densified points.
        :rtype: List[Point].
        """

        if not isinstance(densify_distance, numbers.Real):
            raise Exception("Input densify distance must be float or integer")

        if not math.isfinite(densify_distance):
            raise Exception("Input densify distance must be finite")

        if densify_distance <= 0.0:
            raise Exception("Input densify distance must be positive")

        length2d = self.length2D()

        vect = self.vector()
        vers_2d = vect.versor2D()
        generator_vector = vers_2d.scale(densify_distance)

        pts = [self.start_pt()]

        n = 0
        while True:
            n += 1
            new_pt = self.start_pt().shiftByVect(generator_vector.scale(n))
            distance = self.start_pt().dist2DWith(new_pt)
            if distance >= length2d:
                break
            pts.append(new_pt)

        pts.append(self.end_pt())

        return pts

    def densify2d_asLine(self, densify_distance) -> 'Line':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: numbers.Real
        :return: Line
        """

        pts = self.densify2d_asPts(densify_distance=densify_distance)

        return Line(
            pts=pts)

    def vertical_plane(self) -> Optional[CPlane]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length2D() == 0.0:
            return None

        # arbitrary point on the same vertical as end point

        section_final_pt_up = self.end_pt().shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane.fromPoints(
            pt1=self.start_pt(),
            pt2=self.end_pt(),
            pt3=section_final_pt_up)

    def same_start(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the two segments have the same start point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(0,0,0), Point(0,1,0))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt().isCoinc(
            another=another.start_pt(),
            tolerance=tol
        )

    def same_end(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the two segments have the same end point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(2,0,0), Point(1,0,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt().isCoinc(
            another=another.end_pt(),
            tolerance=tol)

    def conn_to_other(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the first segment is sequentially connected to the second one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(1,0,0), Point(2,0,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt().isCoinc(
            another=another.start_pt(),
            tolerance=tol)

    def other_connected(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the second segment is sequentially connected to the first one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-1,0,0), Point(0,0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt().isCoinc(
            another=self.start_pt(),
            tolerance=tol)

    def segment_start_in(self,
        another: 'Segment'
    ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-0.5,0,0), Point(0.5,0,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment(Point(0,0,0), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment(Point(0,1,0), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = Segment(Point(-1,-1,-1), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt())

    def segment_end_in(self,
        another: 'Segment'
    ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-0.5,0,0), Point(0.5,0,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment(Point(0,0,0), Point(1,1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment(Point(0,1,0), Point(1,1,1))
          >>> s2 = Segment(Point(1,1,1), Point(0.5,0,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = Segment(Point(-1,-1,3), Point(1,1,3))
          >>> s2 = Segment(Point(0,2,3), Point(2,0,3))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt())


class Line(object):
    """
    A list of Point objects, all with the same CRS code.
    """

    def __init__(self, pts: Optional[List[Point]] = None, epsg_cd: numbers.Integral = -1):
        """
        Creates the Line instance, when all the provided points have the same CRS codes.

        :param pts: a list of points
        :type pts: List of Point instances.
        :param epsg_cd: the CRS code of the points.
        :type epsg_cd: numbers.Integral.
        :return: a Line instance.
        :rtype: Line.
        :raises: CRSCodeException.

        """

        if pts is None:
            pts = []

        for pt in pts:
            if not isinstance(pt, Point):
                raise Exception("All input data must be point")

        # when implicit (-1) EPSG line code, initialize it to that of the first point

        if pts and epsg_cd == -1:
            epsg_cd = pts[0].epsg()

        # check all points have the same CRS

        for ndx in range(len(pts)):
            pt = pts[ndx]
            if pt.epsg() != epsg_cd:
                raise Exception("All points must have the same '{}' EPSG code".format(epsg_cd))

        self._pts = [pt.clone() for pt in pts]
        self._crs = Crs(epsg_cd)

    @classmethod
    def fromPointList(cls, pt_list: List[List[numbers.Real]], epsg_cd: numbers.Integral = -1) -> 'Line':
        """
        Create a Line instance from a list of x, y and optional z values.

        Example:
          >>> Line.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          Line with 3 points: (0.0000, 0.0000, 0.0000) ... (0.0000, 1.0000, 0.0000) - EPSG: -1
        """

        pts = []
        for vals in pt_list:
            if len(vals) == 2:
                pt = Point(
                    x=vals[0],
                    y=vals[1],
                    epsg_cd=epsg_cd)
            elif len(vals) == 3:
                pt = Point(
                    x=vals[0],
                    y=vals[1],
                    z=vals[2],
                    epsg_cd=epsg_cd)
            elif len(vals) == 3:
                pt = Point(
                    x=vals[0],
                    y=vals[1],
                    z=vals[2],
                    t=vals[3],
                    epsg_cd=epsg_cd)
            else:
                raise Exception("Point input values should be 2, 3 or 4. {} got ({}).".format(len(vals), vals))

            pts.append(pt)

        return cls(pts, epsg_cd=epsg_cd)

    def extract_pts(self):

        return self._pts

    def extract_pt(self, pt_ndx: numbers.Integral) -> Optional[Point]:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :type pt_ndx: numbers.Integral.
        :return: the extracted Point instance or None when index out-of-range.
        :rtype: Optional[Point].

        Examples:
        """

        num_pts = self.num_pts()

        if num_pts == 0:
            return None
        elif pt_ndx not in range(num_pts):
            return None
        else:
            return self._pts[pt_ndx]

    def pts(self):

        return [pt.clone() for pt in self._pts]

    def crs(self) -> Crs:

        return self._crs

    def epsg(self) -> numbers.Integral:

        return self.crs().epsg()

    def num_pts(self):

        return len(self._pts)

    def start_pt(self) -> Optional[Point]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        :rtype: optional Point instance.
        """

        if self.num_pts() >= 1:
            return self._pts[0].clone()
        else:
            return None

    def end_pt(self) -> Optional[Point]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        :rtype: optional Point instance.
        """

        if self.num_pts() >= 1:
            return self._pts[-1].clone()
        else:
            return None

    def __repr__(self) -> str:
        """
        Represents a Line instance as a shortened text.

        :return: a textual shortened representation of a Line instance.
        :rtype: str.
        """

        num_points = self.num_pts()
        epsg = self.epsg()

        if num_points == 0:
            txt = "Empty Line - EPSG: {}".format(epsg)
        else:
            first_pt = self.start_pt()
            x1, y1, z1 = first_pt.x, first_pt.y, first_pt.z
            if num_points == 1:
                txt = "Line with unique point: {.4f}.{.4f},{.4f} - EPSG: {}".format(x1, y1, z1, epsg)
            else:
                last_pt = self.end_pt()
                x2, y2, z2 = last_pt.x, last_pt.y, last_pt.z
                txt = "Line with {} points: ({:.4f}, {:.4f}, {:.4f}) ... ({:.4f}, {:.4f}, {:.4f}) - EPSG: {}".format(num_points, x1, y1, z1, x2, y2, z2, epsg)

        return txt

    def clone(self):

        return Line(
            pts=self._pts,
            epsg_cd=self.epsg()
        )

    def add_pt(self, pt) -> bool:
        """
        In-place transformation of the original Line instance
        by adding a new point at the end.

        :param pt: the point to add
        :type pt: Point.
        :return: status of addition. True when added, False otherwise.
        :rtype: bool.
        """

        if self.num_pts() == 0 and not self.crs().valid():
            self._crs = Crs(pt.epsg())

        if self.num_pts() > 0 and pt.crs() != self.crs():
            return False

        self._pts.append(pt.clone())
        return True

    def add_pts(self, pt_list) -> numbers.Integral:
        """
        In-place transformation of the original Line instance
        by adding a new set of points at the end.

        :param pt_list: list of Points.
        :type pt_list: List of Point instances.
        :return: number of added points
        :rtype: numbers.Integral.
        """

        num_added = 0
        for pt in pt_list:
            success = self.add_pt(pt)
            if success:
                num_added += 1

        return num_added

    def x_list(self) -> List[numbers.Real]:

        return [pt.x for pt in self._pts]

    def y_list(self) -> List[numbers.Real]:

        return [pt.y for pt in self._pts]

    def z_list(self) -> List[numbers.Real]:

        return [pt.z for pt in self._pts]

    def t_list(self) -> List[numbers.Real]:

        return [pt.t for pt in self._pts]

    def z_array(self) -> np.array:

        return np.array(self.z_list())

    def xy_lists(self) -> Tuple[List[numbers.Real], List[numbers.Real]]:

        return self.x_list(), self.y_list()

    def x_min(self) -> Optional[numbers.Real]:

        return find_val(
            func=min,
            lst=self.x_list())

    def x_max(self) -> Optional[numbers.Real]:

        return find_val(
            func=max,
            lst=self.x_list())

    def y_min(self) -> Optional[numbers.Real]:

        return find_val(
            func=min,
            lst=self.y_list())

    def y_max(self) -> Optional[numbers.Real]:

        return find_val(
            func=max,
            lst=self.y_list())

    def z_stats(self) -> Dict:
        """
        Returns the line elevation statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Dictionary of numbers.Real values.
        """

        return get_statistics(self.z_array())

    def z_min(self) -> Optional[numbers.Real]:

        return find_val(
            func=min,
            lst=self.z_list())

    def z_max(self) -> Optional[numbers.Real]:

        return find_val(
            func=max,
            lst=self.z_list())

    def z_mean(self) -> Optional[numbers.Real]:

        zs = self.z_list()
        return float(np.mean(zs)) if zs else None

    def z_var(self) -> Optional[numbers.Real]:

        zs = self.z_list()
        return float(np.var(zs)) if zs else None

    def z_std(self) -> Optional[numbers.Real]:

        zs = self.z_list()
        return float(np.std(zs)) if zs else None

    def remove_coincident_points(self) -> 'Line':
        """
        Remove coincident successive points

        :return: Line instance
        """

        new_line = Line(
            pts=self._pts[:1])

        for ndx in range(1, self.num_pts()):
            if not self._pts[ndx].isCoinc(new_line._pts[-1]):
                new_line.add_pt(self._pts[ndx])

        return new_line

    def as_segments(self):
        """
        Convert to a list of segments.

        :return: list of Segment objects
        """

        pts_pairs = zip(self._pts[:-1], self._pts[1:])

        segments = [Segment(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    def densify_2d_line(self, sample_distance) -> 'Line':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line instance
        """

        if sample_distance <= 0.0:
            raise Exception("Sample distance must be positive. {} received".format(sample_distance))

        segments = self.as_segments()

        densified_line_list = [segment.densify2d_asLine(sample_distance) for segment in segments]

        densifyied_multiline = MultiLine(densified_line_list, epsg_cd=self.epsg())

        densifyied_line = densifyied_multiline.to_line()

        densifyied_line_wo_coinc_pts = densifyied_line.remove_coincident_points()

        return densifyied_line_wo_coinc_pts

    def join(self, another) -> 'Line':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line(self._pts + another._pts)

    def length_3d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self._pts[ndx].dist3DWith(self._pts[ndx + 1])
        return length

    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self._pts[ndx].dist2DWith(self._pts[ndx + 1])
        return length

    def step_delta_z(self) -> List[numbers.Real]:
        """
        Return the difference in elevation between consecutive points:
        z[ndx+1] - z[ndx]

        :return: a list of height differences.
        :rtype: list of floats.
        """

        delta_z = [0.0]

        for ndx in range(1, self.num_pts()):
            delta_z.append(self._pts[ndx].z - self._pts[ndx - 1].z)

        return delta_z

    def step_lengths_3d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 3D distances.
        It is the distance between a point and its previous one.
        The list has the same lenght as the source point list.

        :return: the individual 3D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self._pts[ndx].dist3DWith(self._pts[ndx - 1])
            step_length_list.append(length)

        return step_length_list

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
            length = self._pts[ndx].dist2DWith(self._pts[ndx - 1])
            step_length_list.append(length)

        return step_length_list

    def incremental_length_3d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 3D segment lengths.

        :return: accumulated 3D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_3d()))

    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))

    def reversed(self) -> 'Line':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        new_line = self.clone()
        new_line._pts.reverse()  # in-place operation on new_line

        return new_line

    def slopes_degr(self) -> List[Optional[numbers.Real]]:
        """
        Calculates the slopes (in degrees) of each Line segment.
        The first value is the slope of the first segment.
        The last value, always None, is the slope of the segment starting at the last point.
        The number of elements is equal to the number of points in the Line.

        :return: list of slopes (degrees).
        :rtype: List[Optional[numbers.Real]].
        """

        lSlopes = []

        segments = self.as_segments()
        for segment in segments:
            vector = segment.vector()
            lSlopes.append(-vector.slope_degr())  # minus because vector convention is positive downward

        lSlopes.append(None)  # None refers to the slope of the Segment starting with the last point

        return lSlopes

    def slopes_stats(self) -> Dict:
        """
        Returns the line directional slope statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Dictionary.
        """

        return get_statistics(self.slopes_degr())

    def abs_slopes_degr(self) -> List[Optional[numbers.Real]]:

        return [abs(val) for val in self.slopes_degr()]

    def abs_slopes_stats(self) -> Dict:
        """
        Returns the line absolute slopes statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Dictionary.
        """

        return get_statistics(self.abs_slopes_degr())


def analizeJoins(first: Union[Line, Segment], second: Union[Line, Segment]) -> List[Optional[JoinTypes]]:
    """
    Analyze join types between two lines/segments.

    :param first: a line or segment.
    :type first: Line or Segment.
    :param second: a line or segment.
    :param second: Line or Segment.
    :return: a list of join types.
    :rtype: List[Optional[JoinTypes]].

    Examples:
      >>> first = Segment(Point(x=0,y=0, epsg_cd=32632), Point(x=1,y=0, epsg_cd=32632))
      >>> second = Segment(Point(x=1,y=0, epsg_cd=32632), Point(x=0,y=0, epsg_cd=32632))
      >>> analizeJoins(first, second)
      [<JoinTypes.START_END: 2>, <JoinTypes.END_START: 3>]
      >>> first = Segment(Point(x=0,y=0, epsg_cd=32632), Point(x=1,y=0, epsg_cd=32632))
      >>> second = Segment(Point(x=2,y=0, epsg_cd=32632), Point(x=3,y=0, epsg_cd=32632))
      >>> analizeJoins(first, second)
      []
    """

    join_types = []

    if first.start_pt().isCoinc(second.start_pt()):
        join_types.append(JoinTypes.START_START)

    if first.start_pt().isCoinc(second.end_pt()):
        join_types.append(JoinTypes.START_END)

    if first.end_pt().isCoinc(second.start_pt()):
        join_types.append(JoinTypes.END_START)

    if first.end_pt().isCoinc(second.end_pt()):
        join_types.append(JoinTypes.END_END)

    return join_types


class MultiLine(object):
    """
    MultiLine is a list of Line objects, each one with the same CRS code
    """

    def __init__(self, lines: Optional[List[Line]] = None, epsg_cd: numbers.Integral = -1):

        if lines is None:
            lines = []

        if lines and epsg_cd == -1:
            epsg_cd = lines[0].epsg()

        for ndx in range(len(lines)):
            if lines[ndx].epsg() != epsg_cd:
                raise Exception("Input line with index {} should have EPSG code {} but has {}".format(
                    ndx,
                    epsg_cd,
                    lines[ndx].epsg()
                ))

        self._lines = lines
        self._crs = Crs(epsg_cd)

    def lines(self):

        return self._lines

    def crs(self):

        return self._crs

    def epsg(self) -> numbers.Integral:

        return self._crs.epsg()

    def num_lines(self):

        return len(self.lines())

    def num_tot_pts(self) -> numbers.Integral:

        num_points = 0
        for line in self._lines:
            num_points += line.num_pts()

        return num_points

    def extract_line(self, ln_ndx: numbers.Integral = 0) -> Optional[Line]:
        """
        Extracts a line from the multiline instance, based on the provided index.

        :return: Line instance or None when ln_ndx is out-of-range.
        :rtype: Optional[Line].
        """

        num_lines = self.num_lines()
        if num_lines == 0:
            return None

        if ln_ndx not in range(num_lines):
            return None

        return self.lines()[ln_ndx]

    def __repr__(self) -> str:
        """
        Represents a MultiLine instance as a shortened text.

        :return: a textual shortened representation of a MultiLine instance.
        :rtype: basestring.
        """

        num_lines = self.num_lines()
        num_tot_pts = self.num_tot_pts()
        epsg = self.epsg()

        txt = "MultiLine with {} line(s) and {} total point(s) - EPSG: {}".format(num_lines, num_tot_pts, epsg)

        return txt

    def add_line(self, line) -> bool:
        """
        In-place addition of a Line instance (that is not cloned).

        :param line: the line to add.
        :type line: Line.
        :return: status of addition. True when added, False otherwise.
        :rtype: bool.
        """

        if self.num_lines() == 0 and not self.crs().valid():
            self._crs = line.crs()

        if self.num_lines() > 0 and line.crs() != self.crs():
            return False

        self._lines += [line]
        return True

    def clone(self) -> 'MultiLine':

        return MultiLine(
            lines=[line.clone() for line in self._lines],
            epsg_cd=self.epsg()
        )

    def x_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.x_min() for line in self.lines()]))

    def x_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.x_max() for line in self.lines()]))

    def y_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.y_min() for line in self.lines()]))

    def y_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.y_max() for line in self.lines()]))

    def z_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.z_min() for line in self.lines()]))

    def z_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.z_max() for line in self.lines()]))

    def is_continuous(self) -> bool:
        """
        Checks whether all lines in a multiline are connected.

        :return: whether all lines are connected.
        :rtype: bool.
        """

        if len(self._lines) <= 1:
            return False

        for line_ndx in range(len(self._lines) - 1):
            first = self._lines[line_ndx]
            second = self._lines[line_ndx + 1]
            if not analizeJoins(first, second):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines()) - 1):
            if not self.lines()[line_ndx].extract_pts()[-1].isCoinc(self.lines()[line_ndx + 1].extract_pts()[0]):
                return False

        return True

    def to_line(self):

        return Line([point for line in self._lines for point in line.extract_pts()], epsg_cd=self.epsg())

    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines():
            lDensifiedLines.append(line.densify_2d_line(sample_distance))

        return MultiLine(lDensifiedLines, self.epsg())

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines():
            cleaned_lines.append(line.remove_coincident_points())

        return MultiLine(cleaned_lines, self.epsg())


class ParamLine3D(object):
    """
    parametric line
    srcPt: source Point
    l, m, n: line coefficients
    """

    def __init__(self, srcPt, l, m, n):

        for v in (l, m, n):
            if not (-1.0 <= v <= 1.0):
                raise Exception("Parametric line values must be in -1 to 1 range")

        self._srcPt = srcPt.clone()
        self._l = l
        self._m = m
        self._n = n

    def intersect_cartes_plane(self, cartes_plane) -> Optional[Point]:
        """
        Return intersection point between parametric line and Cartesian plane.

        :param cartes_plane: a Cartesian plane:
        :type cartes_plane: CPlane.
        :return: the intersection point between parametric line and Cartesian plane.
        :rtype: Point.
        :raise: Exception.
        """

        if not isinstance(cartes_plane, CPlane):
            raise Exception("Method argument should be a Cartesian plane but is {}".format(type(cartes_plane)))

        # line parameters
        x1, y1, z1 = self._srcPt.x, self._srcPt.y, self._srcPt.z
        l, m, n = self._l, self._m, self._n

        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a, cartes_plane.b, cartes_plane.c, cartes_plane.d

        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return Point(x1 - l * k,
                     y1 - m * k,
                     z1 - n * k)


# TODO class Path
"""
The trajectory of a Point with time
"""


if __name__ == "__main__":

    import doctest
    doctest.testmod()
