import numbers
from math import sqrt
from typing import Tuple, Optional

import numpy as np

from pygsf.geometries.geom3d.shapes import Point
from pygsf.mathematics.arrays import pointSolution
from pygsf.mathematics.defaults import MIN_SEPARATION_THRESHOLD
from pygsf.mathematics.vectors import Vect
from pygsf.utils.types import check_type


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: CPlane is locational - its position in space is defined.
    This contrast with Plane, defined just by its attitude, but with undefined position

    """

    def __init__(self,
                 a: numbers.Real,
                 b: numbers.Real,
                 c: numbers.Real,
                 d: numbers.Real
                 ):

        if not isinstance(a, numbers.Real):
            raise Exception("Input value a must be float or int but is {}".format(type(a)))
        if not isinstance(b, numbers.Real):
            raise Exception("Input value b must be float or int but is {}".format(type(b)))
        if not isinstance(c, numbers.Real):
            raise Exception("Input value c must be float or int but is {}".format(type(c)))
        if not isinstance(d, numbers.Real):
            raise Exception("Input value d must be float or int but is {}".format(type(d)))

        norm = sqrt(a*a + b*b + c*c)
        self._a = float(a) / norm
        self._b = float(b) / norm
        self._c = float(c) / norm
        self._d = float(d) / norm

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

    def v(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, numbers.Real]:
        """
        Return coefficients of a CPlane instance.

        Example:
          >>> CPlane(1, 1, 7, -4).v()
          (0.14002800840280097, 0.14002800840280097, 0.9801960588196068, -0.5601120336112039)
        """

        return self.a(), self.b(), self.c(), self.d()

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3) -> 'CPlane':
        """
        Create a CPlane from three given Point instances.

        Example:
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          CPlane(1.0000, 0.0000, 0.0000, 0.0000)
          >>> CPlane.fromPoints(Point(1,2,3), Point(2,3,4), Point(-1,7,-2))
          CPlane(-0.7956, 0.2387, 0.5569, -1.3524)
        """

        if not (isinstance(pt1, Point)):
            raise Exception("First input point should be Point but is {}".format(type(pt1)))

        if not (isinstance(pt2, Point)):
            raise Exception("Second input point should be Point but is {}".format(type(pt2)))

        if not (isinstance(pt3, Point)):
            raise Exception("Third input point should be Point but is {}".format(type(pt3)))

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
            np.linalg.det(matr_d)
        )

    def __repr__(self):

        return "CPlane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v())

    def normVersor(self) -> Vect:
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane(0, 0, 5, -2).normVersor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> CPlane(0, 7, 0, 5).normVersor()
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a(), self.b(), self.c()).versor()

    def toPoint(self) -> Point:
        """
        Returns a point lying in the plane (non-unique solution).

        Examples:
          >>> CPlane(0, 0, 1, -1).toPoint()
          Point(0.0000, 0.0000, 1.0000)
        """

        point = Point(
            *pointSolution(
                np.array([[self.a(), self.b(), self.c()]]),
                np.array([-self.d()]))
        )

        return point

    def intersVersor(self, another) -> Optional[Vect]:
        """
        Return intersection versor for two intersecting planes.
        Return None for not intersecting planes.

        :param another: another Cartesian plane.
        :type another: CPlane.
        :return: the intersection line as a vector.
        :rtype: Optional[Vect].
        :raise: Exception.

        Examples:
          >>> a = CPlane(1, 0, 0, 0)
          >>> b = CPlane(0, 0, 1, 0)
          >>> a.intersVersor(b)
          Vect(0.0000, -1.0000, 0.0000)
          >>> b = CPlane(-1, 0, 0, 0)  # parallel plane, no intersection
          >>> a.intersVersor(b) is None
          True
        """

        check_type(another, "Input Cartesian plane", CPlane)

        return self.normVersor().vCross(another.normVersor()).versor()

    def intersPoint(self,
            another) -> Optional[Point]:
        """
        Return point on intersection line (non-unique solution)
        for two planes.

        :param another: the second cartesian plane
        :type another: CPlane
        :return: the optional instersection point
        :rtype: Optional[Point]
        :raise: Exception

        Examples:
          >>> p_a = CPlane(1, 0, 0, 0)
          >>> p_b = CPlane(0, 0, 1, 0)
          >>> p_a.intersPoint(p_b)
          Point(0.0000, 0.0000, 0.0000, 0.0000)
          >>> p_b = CPlane(-1, 0, 0, 0)  # parallel plane, no intersection
          >>> p_a.intersPoint(p_b) is None
        """

        check_type(another, "Second plane", CPlane)

        # find a point lying on the intersection line (this is a non-unique solution)

        a = np.array([[self.a(), self.b(), self.c()], [another.a(), another.b(), another.c()]])
        b = np.array([-self.d(), -another.d()])
        x, y, z = pointSolution(a, b)

        if x is not None and y is not None and z is not None:
            return Point(x, y, z)
        else:
            return None

    def pointDistance(self,
        pt: Point
    ) -> numbers.Real:
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
          >>> cpl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 0, 1)
          >>> cpl.pointDistance(pt)
          1.0
          >>> pt = Point(0, 0, 0.5)
          >>> cpl.pointDistance(pt)
          0.5
          >>> pt = Point(0, 0, -0.5)
          >>> cpl.pointDistance(pt)
          -0.5
          >>> pt = Point(10, 20, 0.0)
          >>> cpl.pointDistance(pt)
          0.0
        """

        check_type(pt, "Input point", Point)

        return self.a() * pt.x + self.b() * pt.y + self.c() * pt.z + self.d()

    def isPointInPlane(self,
        pt: Point
    ) -> bool:
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
          >>> pl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
        """

        check_type(pt, "Input point", Point)

        if abs(self.pointDistance(pt)) < MIN_SEPARATION_THRESHOLD:
            return True
        else:
            return False

    def angle(self,
        another: 'CPlane'
    ) -> numbers.Real:
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
          >>> CPlane(1,0,0,0).angle(CPlane(0,1,0,0))
          90.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,1,0))
          45.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,0,0))
          0.0
        """

        check_type(another, "Second Cartesian plane", CPlane)

        angle_degr = self.normVersor().angle(another.normVersor())

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


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

    '''
    def epsg(self) -> numbers.Integral:
        """
        Return the EPSG code of the parametric line.
        """

        return self._srcPt.epsg_code()
    '''

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

        '''
        if cartes_plane.epsg_code() != self.epsg():
            raise Exception("Parametric line has EPSG {} while Cartesian plane has {}".format(self.epsg(), cartes_plane.epsg_code()))
        '''

        # line parameters
        x1, y1, z1 = self._srcPt.x, self._srcPt.y, self._srcPt.z
        l, m, n = self._l, self._m, self._n
        # Cartesian plane parameters
        a, b, c, d = cartes_plane.a(), cartes_plane.b(), cartes_plane.c(), cartes_plane.d()
        try:
            k = (a * x1 + b * y1 + c * z1 + d) / (a * l + b * m + c * n)
        except ZeroDivisionError:
            return None

        return Point(
            x=x1 - l * k,
            y=y1 - m * k,
            z=z1 - n * k
        )