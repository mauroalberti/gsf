
from typing import List, Optional, Tuple
from enum import Enum

import itertools
import numbers
from math import *
from array import array

from shapely.geometry import LineString

from pygsf.crs.geoshapes import Points, PointSegmentCollection, MultiLine
from pygsf.orientations.orientations import rotVectByAxis

from pygsf.mathematics.statistics import *
from pygsf.mathematics.quaternions import *
from pygsf.geometries.geom2d.shapes import *


class Point:
    """
    Cartesian point.
    Dimensions: 3D
    """

    def __init__(
        self,
        x: numbers.Real,
        y: numbers.Real,
        z: numbers.Real = 0.0
    ):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :type x: numbers.Real.
        :param y: point y coordinate.
        :type y: numbers.Real.
        :param z: point z coordinate.
        :type z: numbers.Real.
        """

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer or float type")
        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite")

        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    @classmethod
    def fromVect(cls,
        vect: Vect) -> 'Point':
        """

        :param vect:
        :return:
        """

        return cls(
            x=vect.x,
            y=vect.y,
            z=vect.z
        )

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7).x
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
          >>> Point(4, 3, 7).y
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
          >>> Point(4, 3, 7).z
          7.0
          >>> Point(-0.39, 17.42, 8.9).z
          8.9
        """

        return self._z

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y, z = Point(1,1)
          >>> x == 1
          True
          >>> y == 1
          True

        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self,
        another: 'Point'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, -1)
          False
        """

        if not isinstance(another, Point):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            self.z == another.z
            ]
        )

    def __ne__(self,
        another: 'Point'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
          >>> Point(1., 1., 1.) != Point(1, 1, 1)
          True
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Return the individual values of the point.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7).a()
          (4.0, 3.0, 7.0)
        """

        return self.x, self.y, self.z

    def __add__(self, another: 'Point') -> 'Point':
        """
        Sum of two points.

        :param another: the point to add
        :type another: Point
        :return: the sum of the two points
        :rtype: Point
        :raise: Exception

        Example:
          >>> Point(1, 0, 0) + Point(0, 1, 1)
          Point(1.0000, 1.0000, 1.0000)
          >>> Point(1, 1, 1) + Point(-1, -1, -1)
          Point(0.0000, 0.0000, 0.0000)
        """

        check_type(another, "Second point", Point)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point(
            x=x0+x1,
            y=y0+y1,
            z=z0+z1
        )

    def __sub__(self,
        another: 'Point'
    ) -> 'Point':
        """Subtract two points.

        :param another: the point to subtract
        :type another: Point
        :return: the difference between the two points
        :rtype: Point
        :raise: Exception

        Example:
          >>> Point(1., 1., 1.) - Point(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000)
          >>> Point(1., 1., 3.) - Point(1., 1., 2.2)
          Point(0.0000, 0.0000, 0.8000)
        """

        check_type(another, "Second point", Point)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point(
            x=x0 - x1,
            y=y0 - y1,
            z=z0 - z1
        )

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

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Point(1, 2, 3).toArray(), np.array([ 1., 2., 3., 0.]))
          True
        """

        return np.asarray(self.toXYZ())

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return Point(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return Point(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return Point(0.0, self.y, self.z)

    def deltaX(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        check_crs(self, another)

        return another.x - self.x

    def deltaY(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        check_crs(self, another)

        return another.y - self.y

    def deltaZ(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between z components of two Point Instances.

        :return: z coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        check_crs(self, another)

        return another.z - self.z

    def dist3DWith(self,
        another: 'Point'
    ) -> numbers.Real:
        """
        Calculate Euclidean spatial distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point(1., 1., 1.).dist3DWith(Point(4., 5., 1))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        check_type(another, "Point", Point)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    '''
    def dist2DWith(self,
        another: 'Point'
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
          >>> Point(1., 1., 1.).dist2DWith(Point(4., 5., 7.))
          5.0
        """

        check_type(another, "Second point", Point)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)
    '''

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'Point':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
        """

        x, y, z = self.x * scale_factor, self.y * scale_factor, self.z * scale_factor
        return Point(x, y, z)

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def reflect_vertical(self) -> 'Point':
        """
        Reflect a point along a vertical axis.

        :return: reflected point.
        :rtype: Point

        Examples:
          >>> Point(1,1,1).reflect_vertical()
          Point(-1.0000, -1.0000, 1.0000)
        """

        x, y, z = self

        return Point(
            x=-x,
            y=-y,
            z=z
        )

    def isCoinc3D(self,
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
          >>> Point(1., 0., -1.).isCoinc3D(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc3D(Point(1., 0., 0.))
          True
        """

        check_type(another, "Second point", Point)

        return self.dist3DWith(another) <= tolerance

    def already_present(self,
        pt_list: List['Point'],
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
            if self.isCoinc3D(pt, tolerance=tolerance):
                return True
        return False

    def shift(self,
        sx: numbers.Real,
        sy: numbers.Real,
        sz: numbers.Real
    ) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz)

    def shiftByVect(self,
        v: Vect
    ) -> 'Point':
        """
        Create a new point shifted from the self instance by given vector.

        :param v: the shift vector.
        :type v: Vect.
        :return: the shifted point.
        :rtype: Point.
        :raise: Exception

        Example:
          >>> Point(1, 1, 1).shiftByVect(Vect(0.5, 1., 1.5))
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shiftByVect(Vect(0.5, 1., 1.5))
          Point(1.5000, 3.0000, 0.5000)
       """

        x, y, z = self

        sx, sy, sz = v.toXYZ()

        return Point(x + sx, y + sy, z + sz)

    def asVect(self) -> 'Vect':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0).asVect()
          Vect(1.0000, 1.0000, 0.0000)
          >>> Point(0.2, 1, 6).asVect()
          Vect(0.2000, 1.0000, 6.0000)
        """

        return Vect(self.x, self.y, self.z)

    def rotate(self,
        rotation_axis: 'RotationAxis',
        center_point: 'Point' = None
        ) -> 'Point':
        """
        Rotates a point.
        :param rotation_axis:
        :param center_point:
        :return: the rotated point
        :rtype: Point

        Examples:
          >>> pt = Point(0,0,1)
          >>> rot_axis = RotationAxis(0,0,90)
          >>> center_pt = Point(0,0,0.5)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(0.5000, 0.0000, 0.5000)
          >>> center_pt = Point(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(0.0000, 0.0000, 1.0000)
          >>> center_pt = Point(0, 0, 2)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-1.0000, 0.0000, 2.0000)
          >>> rot_axis = RotationAxis(0,0,180)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-0.0000, 0.0000, 3.0000)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(0.0000, 0.0000, -1.0000)
          >>> pt = Point(1,1,1,5)
          >>> rot_axis = RotationAxis(0,90,90)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(1.0000, -1.0000, 1.0000)
          >>> rot_axis = RotationAxis(0,90,180)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(-1.0000, -1.0000, 1.0000)
          >>> center_pt = Point(1,1,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(1.0000, 1.0000, 1.0000)
          >>> center_pt = Point(2,2,10)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(3.0000, 3.0000, 1.0000)
          >>> pt = Point(1, 1, 2)
          >>> rot_axis = RotationAxis(135, 0, 180)
          >>> center_pt = Point(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-1.0000, -1.0000, 0.0000)
        """

        if not center_point:

            center_point = Point(
                x=0.0,
                y=0.0,
                z=0.0
            )

        check_type(center_point, "Center point", Point)

        p_diff = self - center_point

        p_vect = p_diff.asVect()

        rot_vect = rotVectByAxis(
            v=p_vect,
            rot_axis=rotation_axis
        )

        x, y, z, epsg_cd = rot_vect

        rot_pt = Point(
            x=x,
            y=y,
            z=z
        )

        transl_pt = center_point + rot_pt

        return transl_pt

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: Point
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(3)]
        return cls(*vals)


def pack_to_points(
    xs: array,
    ys: array,
    zs: Optional[array] = None,
) -> List[Point]:
    # Side effects: None
    """
    Create a list of points given a set
    of input arrays.

    :param xs: array of x values
    :type xs: array
    :param ys: array of y values
    :type ys: array
    :param zs: optional array of z values
    :type zs: Optional[array]
    :return: a list of Point instances
    :rtype: List[Point]
    """

    if zs is None:
        zs = [0.0] * len(xs)
    pts = []
    for x, y, z, t in zip(xs, ys, zs):
        pts.append(
            Point(
                x,
                y,
                z
            )
        )

    return pts


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


class Segment:
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

        check_type(start_pt, "Start point", Point)

        check_type(end_pt, "End point", Point)

        if start_pt.dist3DWith(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()

    def __repr__(self) -> str:
        """
        Represents a Segment instance.

        :return: the Segment representation.
        :rtype: str.
        """

        return "Segment(start_pt={}, end_pt={})".format(
            self.start_pt,
            self.end_pt
        )

    @property
    def start_pt(self) -> Point:

        return self._start_pt

    @property
    def end_pt(self) -> Point:

        return self._end_pt

    def __iter__(self):
        """
        Return the elements of a Segment, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'Segment':

        return Segment(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment':

        if self.end_pt.x < self.start_pt.x:
            return Segment(self.end_pt, self.start_pt)
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

    def z_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.z < self.end_pt.z:
            return self.start_pt.z, self.end_pt.z
        else:
            return self.end_pt.z, self.start_pt.z

    def delta_x(self) -> numbers.Real:

        return self.end_pt.x - self.start_pt.x

    def delta_y(self) -> numbers.Real:

        return self.end_pt.y - self.start_pt.y

    def delta_z(self) -> numbers.Real:
        """
        Z delta between segment end point and start point.

        :return: numbers.Real.
        """

        return self.end_pt.z - self.start_pt.z

    def length2D(self) -> numbers.Real:
        """
        Returns the horizontal length of the segment.

        :return: the horizontal length of the segment.
        :rtype: numbers.Real.
        """

        return self.start_pt.dist2DWith(self.end_pt)

    def length3D(self) -> numbers.Real:

        return self.start_pt.dist3DWith(self.end_pt)

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
                    self.delta_z()
        )

    def antivector(self) -> Vect:
        """
        Returns the vector pointing from the segment end to the segment start.

        :return: the vector pointing from the segment end to the segment start.
        :rtype: Vect.
        """

        return self.vector().invert()

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

    def intersection_2d_pt(self,
        another: 'Segment'
    ) -> Optional[Point2D]:
        """

        :param another:
        :return:
        """

        check_type(another, "Second segment", Segment)

        check_crs(self, another)

        s_len2d = self.length2D()
        a_len2d = another.length2D()

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
        pt: Point
    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment(Point(0, 0, 0), Point(1, 0, 0))
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
          >>> segment = Segment(Point(0, 0, 0), Point(0, 1, 0))
          >>> segment.contains_pt(Point(0, 0, 0))
          True
          >>> segment.contains_pt(Point(0, 0.5, 0))
          True
          >>> segment.contains_pt(Point(0, 1, 0))
          True
          >>> segment.contains_pt(Point(0, 1.5, 0))
          False
          >>> segment = Segment(Point(0, 0, 0), Point(1, 1, 1))
          >>> segment.contains_pt(Point(0.5, 0.5, 0.5))
          True
          >>> segment.contains_pt(Point(1, 1, 1))
          True
          >>> segment = Segment(Point(1,2,3), Point(9,8,2))
          >>> segment.contains_pt(segment.pointAt(0.745))
          True
          >>> segment.contains_pt(segment.pointAt(1.745))
          False
          >>> segment.contains_pt(segment.pointAt(-0.745))
          False
          >>> segment.contains_pt(segment.pointAt(0))
          True
        """

        check_type(pt, "Point", Point)

        segment_length = self.length3D()
        length_startpt_pt = self.start_pt.dist3DWith(pt)
        length_endpt_pt = self.end_pt.dist3DWith(pt)

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

    def pointAt(self,
        scale_factor: numbers.Real
    ) -> Point:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        ans 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point

        Examples:
          >>> s = Segment(Point(0,0,0), Point(1,0,0))
          >>> s.pointAt(0)
          Point(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point(0.5000, 0.0000, 0.0000)
          >>> s.pointAt(1)
          Point(1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-1)
          Point(-1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-2)
          Point(-2.0000, 0.0000, 0.0000)
          >>> s.pointAt(2)
          Point(2.0000, 0.0000, 0.0000)
          >>> s = Segment(Point(0,0,0), Point(0,0,1))
          >>> s.pointAt(0)
          Point(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point(0.0000, 0.0000, 0.5000)
          >>> s.pointAt(1)
          Point(0.0000, 0.0000, 1.0000)
          >>> s.pointAt(-1)
          Point(0.0000.0000, 0.0000, -1)
          >>> s.pointAt(-2)
          Point(0.0000, 0.0000, -2.0000)
          >>> s.pointAt(2)
          Point(0.0000, 0.0000, 2.0000)
          >>> s = Segment(Point(0,0,0), Point(1,1,1))
          >>> s.pointAt(0.5)
          Point(0.5000, 0.5000, 0.5000)
          >>> s = Segment(Point(0,0,0), Point(4,0,0))
          >>> s.pointAt(7.5)
          Point(30.0000, 0.0000, 0.0000)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor
        dz = self.delta_z() * scale_factor

        return Point(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy,
            z=self.start_pt.z + dz
        )

    def pointProjection(self,
        point: Point
    ) -> Point:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = Segment(start_pt=Point(0,0,0), end_pt=Point(1,0,0))
          >>> p = Point(0.5, 1, 4)
          >>> s.pointProjection(p)
          Point(0.5000, 0.0000, 0.0000)
          >>> s = Segment(start_pt=Point(0,0,0), end_pt=Point(4,0,0))
          >>> p = Point(7.5, 19.2, -14.72)
          >>> s.pointProjection(p)
          Point(7.5000, 0.0000, 0.0000)
        """

        check_type(point, "Input point", Point)

        check_crs(self, point)

        other_segment = Segment(
            self.start_pt,
            point
        )
        
        scale_factor = self.vector().scalarProjection(other_segment.vector()) / self.length3D()
        return self.pointAt(scale_factor)

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
          >>> s = Segment(Point(0,0,0), Point(0,0,4))
          >>> s.pointDistance(Point(-17.2, 0.0, -49,3))
          17.2
          >>> s.pointDistance(Point(-17.2, 1.22, -49,3))
          17.24321315764553
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.dist3DWith(point_projection)

    def pointS(self,
        point: Point
    ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :type point: Point
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        if not self.contains_pt(point):
            return None

        return self.start_pt.dist3DWith(point)

    def scale(self,
        scale_factor
    ) -> 'Segment':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point
        """

        end_pt = self.pointAt(scale_factor)

        return Segment(
            self.start_pt,
            end_pt)

    def densify2d_asSteps(self,
        densify_distance: numbers.Real
    ) -> array:
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

    def densify2d_asPts(self,
        densify_distance
    ) -> List[Point]:
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

        pts = [self.start_pt]

        n = 0
        while True:
            n += 1
            new_pt = self.start_pt.shiftByVect(generator_vector.scale(n))
            distance = self.start_pt.dist2DWith(new_pt)
            if distance >= length2d:
                break
            pts.append(new_pt)

        pts.append(self.end_pt)

        return pts

    def densify2d_asLine(self,
        densify_distance
    ) -> 'Points':
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: numbers.Real
        :return: Line
        """

        pts = self.densify2d_asPts(densify_distance=densify_distance)

        return Points(
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

        section_final_pt_up = self.end_pt.shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane.fromPoints(
            pt1=self.start_pt,
            pt2=self.end_pt,
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

        return self.start_pt.isCoinc3D(
            another=another.start_pt,
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

        return self.end_pt.isCoinc3D(
            another=another.end_pt,
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

        return self.end_pt.isCoinc3D(
            another=another.start_pt,
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

        return another.end_pt.isCoinc3D(
            another=self.start_pt,
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

        return another.contains_pt(self.start_pt)

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

        return another.contains_pt(self.end_pt)

    def rotate(self,
        rotation_axis: 'RotationAxis',
        center_point: 'Point' = None
        ) -> 'Segment':
        """
        Rotates a segment.
        :param rotation_axis:
        :param center_point:
        :return: the rotated segment
        :rtype: Segment

        Examples:
        >>> seg = Segment(Point(0,0,0), Point(0,0,1))
        >>> rot_ax = RotationAxis(0, 0, 90)
        >>> seg.rotate(rot_ax)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
        >>> rot_ax = RotationAxis(0, 0, 180)
        >>> seg.rotate(rot_ax)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.0000.0000))
        >>> centr_pt = Point(0,0,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(-0.0000, 0.0000, 1.0000), end_pt=Point(0.0000, 0.0000, 0.0000))
        >>> seg = Segment(Point(0,0,0), Point(1,1,0))
        >>> centr_pt = Point(1,0,0)
        >>> rot_ax = RotationAxis(0, 90, 90)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(1.0000, 1.0000, 0.0000), end_pt=Point(2.0000, 0.0000, -0.0000))
        >>> seg = Segment(Point(1,1,1), Point(0,0,0))
        >>> rot_ax = RotationAxis(135, 0, 180)
        >>> centr_pt = Point(0.5,0.5,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 1.0000, 1.0000))
        """

        start_pt, end_pt = self

        rotated_start_pt = start_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        rotated_end_pt = end_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        return Segment(
            start_pt=rotated_start_pt,
            end_pt=rotated_end_pt
        )

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE):
        """
        Creates a random segment.

        :return: random segment
        :rtype: Segment
        """

        return cls(
            start_pt=Point.random(lower_boundary, upper_boundary),
            end_pt=Point.random(lower_boundary, upper_boundary)
        )


def point_or_segment(
        point1: Point,
        point2: Point,
        tol: numbers.Real = PRACTICAL_MIN_DIST
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

    #check_crs(point1, point2)

    if point1.dist3DWith(point2) <= tol:
        return Points.fromPoints([point1, point2]).nanmean_point()
    else:
        return Segment(
            start_pt=point1,
            end_pt=point2
        )


def intersect_segments(
    segment1: Segment,
    segment2: Segment,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Point, Segment]]:
    """
    Determines the optional point or segment intersection between the segment pair.

    :param segment1: the first segment
    :type segment1: Segment
    :param segment2: the second segment
    :type segment2: Segment
    :param tol: the distance tolerance for collapsing a intersection segment into a point
    :type tol: numbers.Real
    :return: the optional point or segment intersection between the segment pair.
    :rtype: Optional[Union[Point, Segment]]

    Examples:
      >>> s2 = Segment(Point(0,0,0), Point(1,0,0))
      >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(-1,0,0))
      >>> intersect_segments(s1, s2) is None
      True
      >>> s1 = Segment(Point(-2,0,0), Point(0,0,0))
      >>> intersect_segments(s1, s2)
      Point(0.0000, 0.0000, 0.0000)
      >>> s1 = Segment(Point(-2,0,0), Point(0.5,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(2,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0,0,0), Point(0.5,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(0.75,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(0.7500, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(1,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Point(1.0000, 0.0000, 0.0000)
      >>> s2 = Segment(Point(0,0,0), Point(1,1,1))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(0.75,0.75,0.75))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.2500, 0.2500), end_pt=Point(0.7500, 0.7500, 0.7500))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,1.75,1.75))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.2500, 0.2500), end_pt=Point(1.0000, 1.0000, 1.0000))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,0,1.75))
      >>> intersect_segments(s1, s2)
      Point(0.2500, 0.2500, 0.2500)
      >>> s1 = Segment(Point(0.25,1,0.25), Point(0.75,0.75,0.75))
      >>> intersect_segments(s1, s2)
      Point(0.7500, 0.7500, 0.7500)
      >>> s2 = Segment(Point(-1,-1,-1), Point(1,1,1))
      >>> s1 = Segment(Point(-1,1,1), Point(1,-1,-1))
      >>> intersect_segments(s1, s2)
      Point(-0.0000, 0.0000, 0.0000)
    """

    check_type(segment1, "First segment", Segment)
    check_type(segment2, "Second segment", Segment)

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
        return point_or_segment(
            segment1.start_pt,
            segment2.start_pt,
            tol=tol
        )

    if s1_startpt_inside and s2_endpt_inside:
        return point_or_segment(
            segment1.start_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_startpt_inside:
        return point_or_segment(
            segment2.start_pt,
            segment1.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_endpt_inside:
        return point_or_segment(
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

    cline1 = CLine.fromSegment(segment1)
    cline2 = CLine.fromSegment(segment2)

    shortest_segm_or_pt = cline1.shortest_segment_or_point(
        cline2,
        tol=tol
    )

    if not shortest_segm_or_pt:
        return None

    if not isinstance(shortest_segm_or_pt, Point):
        return None

    inters_pt = shortest_segm_or_pt

    if not segment1.contains_pt(inters_pt):
        return None

    if not segment2.contains_pt(inters_pt):
        return None

    return inters_pt


class Segments(list):
    """
    Collection of segments, inheriting from list.

    """

    def __init__(self, segments: List[Segment]):

        check_type(segments, "Segments", List)
        for el in segments:
            check_type(el, "Segment", Segment)

        super(Segments, self).__init__(segments)


class CLine:
    """
    Cartesian line.

    Defined by a point and a unit vector.
    """

    def __init__(self,
        point: Point,
        dir_vector: Vect):

        check_type(point, "Input point", Point)
        check_type(dir_vector, "Directional vector", Vect)

        self._start_pt = point
        self._dir_vect = dir_vector.versor().upward()

    @property
    def start_pt(self) -> Point:
        """
        Returns the Cartesian line point.

        :return: the Cartesian line point
        :rtype: Point
        """

        return self._start_pt

    @property
    def end_pt(self) -> Point:
        """
        Returns the Cartesian line point.

        :return: the Cartesian line point
        :rtype: Point
        """

        return self.start_pt.shiftByVect(self.versor)

    @property
    def versor(self) -> Vect:
        """
        Returns .

        :return: the unit vector
        :rtype: Vect
        """

        return self._dir_vect

    def points(self) -> Tuple[Point, Point]:
        """
        Returns the CLine as a tuple of two points.

        :return: the CLine as a tuple of two points
        :rtype: Tuple[Point, Point]
        """

        return self.start_pt, self.end_pt

    def segment(self) -> Segment:
        """
        Returns the CLine as a segment.

        :return: the CLine as a segment
        :rtype: Segment
        """

        return Segment(
            start_pt=self.start_pt,
            end_pt=self.end_pt
        )

    @classmethod
    def fromPoints(cls,
        first_pt: Point,
        second_pt: Point,
        tolerance = PRACTICAL_MIN_DIST
    ) -> 'CLine':
        """
        Creates a CLine instance from two distinct points.

        :param tolerance:
        :param first_pt: the first input point
        :type first_pt: Point
        :param second_pt: the second input point
        :type second_pt: Point
        :return: a new CLine instance
        :rtype: CLine
        """

        check_type(first_pt, "First point", Point)
        check_type(second_pt, "Second point", Point)

        if first_pt.isCoinc3D(second_pt, tolerance=tolerance):
            raise Exception("The two input points are practically coincident")

        segment = Segment(
            start_pt=first_pt,
            end_pt=second_pt
        )

        return cls(
            point=first_pt,
            dir_vector=segment.vector()
        )

    @classmethod
    def fromSegment(cls,
        segment: Segment):
        """
        Creates a Cartesian line from a segment instance.

        :param segment: the segment to convert to Cartesian line
        :type segment: Segment
        :return: a new CLine
        :rtype: CLine
        """

        return cls.fromPoints(
            first_pt=segment.start_pt,
            second_pt=segment.end_pt
        )

    def shortest_segment_or_point(self,
        another: 'CLine',
        tol: numbers.Real = PRACTICAL_MIN_DIST
    ) -> Optional[Union[Segment, Point]]:

        """
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

        :param another: the second Cartesian line.
        :type another: Cartesian line.
        :param tol: tolerance value for collapsing a segment into the midpoint.
        :type tol: numbers.Real
        :return: the optional shortest segment or an intersection point.
        :rtype: Optional[Union[Segment, Point]]
        """

        check_type(another, "Second Cartesian line", CLine)

        p1 = self.start_pt
        p2 = self.end_pt

        p3 = another.start_pt
        p4 = another.end_pt

        p13 = Point(
            x=p1.x - p3.x,
            y=p1.y - p3.y,
            z=p1.z - p3.z
        )

        p43 = Point(
            x=p4.x - p3.x,
            y=p4.y - p3.y,
            z=p4.z - p3.z
        )

        if p43.asVect().isAlmostZero:
            return None

        p21 = Point(
            x=p2.x - p1.x,
            y=p2.y - p1.y,
            z=p2.z - p1.z,
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
            z=p1.z + mua * p21.z
        )

        pb = Point(
            x=p3.x + mub * p43.x,
            y=p3.y + mub * p43.y,
            z=p3.z + mub * p43.z
        )

        intersection = point_or_segment(
            point1=pa,
            point2=pb,
            tol=tol
        )

        return intersection

    def pointDistance(self,
        point: Point
    ) -> numbers.Real:
        """
        Returns the distance between a line and a point.

        Algorithm from Wolfram MathWorld: Point-Line Distance -- 3-Dimensional

        :param point: input point
        :type point: Point
        :return: the distance
        :rtype: numbers.Real

        Examples:
        """

        v2 = self.end_pt.asVect()
        v1 = self.start_pt.asVect()

        v0 = point.asVect()

        d = abs((v0 - v1).vCross(v0 - v2)) / abs(v2 - v1)

        return d


class Points:
    """
    A list of Point objects.
    """

    def __init__(self,
        pts: Optional[List[Point]] = None
    ):
        """
        Creates the Line instance.

        :param pts: a list of points
        :type pts: List of Point instances.
        :return: a Line instance.
        :rtype: Line.
        """

        if pts is None:
            pts = []

        for pt in pts:
            if not isinstance(pt, Point):
                raise Exception("All input data must be point")

        self._x = array('d', [pt.x for pt in pts])
        self._y = array('d', [pt.y for pt in pts])
        self._z = array('d', [pt.z for pt in pts])

    @classmethod
    def fromArrays(cls,
        xs: array,
        ys: array,
        zs: array = None
    ) -> 'Points':
        """
        Create a Line instance from a list of x, y and optional z values.

        Example:
          >>> Points.fromArrays(xs=array('d',[1,2,3]), ys=array('d', [3,4,5]), zs=array('d',[1,2,3]))
          Line with 3 points: (1.0000, 3.0000, 1.0000) ... (3.0000, 5.0000, 3.0000)
          >>> Points.fromArrays(xs=array('d',[1,2,3]), ys=array('d', [3,4,5]))
          Line with 3 points: (1.0000, 3.0000, 0.0000) ... (3.0000, 5.0000, 0.0000)
        """

        if not isinstance(xs, array):
            raise Exception("X values have type {} instead of array".format(type(xs)))

        if not isinstance(ys, array):
            raise Exception("Y values have type {} instead of array".format(type(ys)))

        if zs is not None and not isinstance(zs, array):
            raise Exception("Z values have type {} instead of array or None".format(type(zs)))

        num_vals = len(xs)
        if len(ys) != num_vals:
            raise Exception("Y array has length {} while x array has length {}".format(len(ys), num_vals))

        if zs is not None and len(zs) != num_vals:
            raise Exception("Z array has length {} while x array has length {}".format(len(zs), num_vals))

        if zs is None:
            zs = array('d', [0.0]*num_vals)

        self = cls()

        self._x = xs
        self._y = ys
        self._z = zs

        return self

    @classmethod
    def fromPointList(cls,
        pt_list: List[List[numbers.Real]]
    ) -> 'Points':
        """
        Create a Line instance from a list of x, y and optional z values.

        Example:
          >>> Points.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          Line with 3 points: (0.0000, 0.0000, 0.0000) ... (0.0000, 1.0000, 0.0000)
        """

        pts = []
        for vals in pt_list:
            if len(vals) == 2:
                pt = Point(
                    x=vals[0],
                    y=vals[1]
                )
            elif len(vals) == 3:
                pt = Point(
                    x=vals[0],
                    y=vals[1],
                    z=vals[2]
                )
            else:
                raise Exception(f"Point input values should be 2 or 3, {len(vals)} got")

            pts.append(pt)

        return cls(pts)

    def pt(self, pt_ndx: numbers.Integral) -> Point:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :type pt_ndx: numbers.Integral.
        :return: the extracted Point instance.
        :rtype: Point.

        Examples:
        """

        return Point(
            x=self._x[pt_ndx],
            y=self._y[pt_ndx],
            z=self._z[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
                  ) -> Tuple[float, float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x, y and z values
        """

        return (
            self._x[ndx],
            self._y[ndx],
            self._z[ndx]
        )

    def pts(self):

        return [Point(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def segment(self,
        ndx: numbers.Integral
    ) -> Optional[Segment]:
        """
        Returns the optional segment at index ndx.

        :param ndx: the segment index.
        :type ndx: numbers.Integral
        :return: the optional segment
        :rtype: Optional[Segment]
        """

        start_pt = self.pt(ndx)
        end_pt = self.pt(ndx + 1)

        if start_pt.isCoinc3D(end_pt):
            return None
        else:
            return Segment(
                start_pt=self.pt(ndx),
                end_pt=self.pt(ndx + 1)
            )

    def intersectSegment(self,
        segment: Segment
    ) -> Optional[PointSegmentCollection]:
        """
        Calculates the possible intersection between the line and a provided segment.

        :param segment: the input segment
        :type segment: Segment
        :return: the optional intersections, points or segments
        :rtype: Optional[List[Optional[Union[Point, Segment]]]]
        :raise: Exception
        """

        if self.num_pts() <= 1:
            return

        check_type(segment, "Input segment", Segment)

        intersections = [intersect_segments(curr_segment, segment) for curr_segment in self if curr_segment is not None]
        intersections = list(filter(lambda val: val is not None, intersections))
        intersections = PointSegmentCollection(intersections)

        return intersections

    def num_pts(self):

        return len(self._x)

    def start_pt(self) -> Optional[Point]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        :rtype: optional Point instance.
        """

        return self.pt(0) if self.num_pts() > 0 else None

    def end_pt(self) -> Optional[Point]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        :rtype: optional Point instance.
        """

        return self.pt(-1) if self.num_pts() > 0 else None

    def __iter__(self):
        """
        Return each element of a Line, i.e., its segments.
        """

        return (self.segment(i) for i in range(self.num_pts()-1))

    def __repr__(self) -> str:
        """
        Represents a Line instance as a shortened text.

        :return: a textual shortened representation of a Line instance.
        :rtype: str.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty Line"
        else:
            x1, y1, z1 = self.start_pt()
            if num_points == 1:
                txt = "Line with unique point: {.4f}.{.4f},{.4f}".format(x1, y1, z1)
            else:
                x2, y2, z2 = self.end_pt()
                txt = "Line with {} points: ({:.4f}, {:.4f}, {:.4f}) ... ({:.4f}, {:.4f}, {:.4f})".format(num_points, x1, y1, z1, x2, y2, z2)

        return txt

    def add_pt(self, pt) -> bool:
        """
        In-place transformation of the original Line instance
        by adding a new point at the end.

        :param pt: the point to add
        :type pt: Point.
        :return: status of addition. True when added, False otherwise.
        :rtype: bool.
        """

        self._x.append(pt.x)
        self._y.append(pt.y)
        self._z.append(pt.z)
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

        return list(self._x)

    def y_list(self) -> List[numbers.Real]:

        return list(self._y)

    def z_list(self) -> List[numbers.Real]:

        return list(self._z)

    def z_array(self) -> np.ndarray:

        return np.array(self._z)

    def xy_lists(self) -> Tuple[List[numbers.Real], List[numbers.Real]]:

        return self.x_list(), self.y_list()

    def xy_zipped(self) -> List[Tuple[numbers.Real, numbers.Real]]:

        return [(x, y) for x, y in zip(self.x_list(), self.y_list())]

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          >>> l.x_min()
          0.0
        """

        return min(self._x) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum of x values.

        :return: the optional maximum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          >>> l.x_max()
          1.0
        """

        return max(self._x) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of y values.

        :return: the optional minimum of y values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          >>> l.y_min()
          0.0
        """

        return min(self._y) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum of y values.

        :return: the optional maximum of y values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          >>> l.y_max()
          1.0
        """

        return max(self._y) if self.num_pts() > 0 else None

    def z_stats(self) -> Optional[Dict]:
        """
        Returns the line elevation statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Optional dictionary of numbers.Real values.
        """

        return get_statistics(self.z_array()) if self.num_pts() > 0 else None

    def z_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of z values.

        :return: the optional minimum of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 0], [-2.72, 1, -0.7]])
          >>> l.z_min()
          -0.7
        """

        return min(self._z) if self.num_pts() > 0 else None

    def z_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum of z values.

        :return: the optional maximum of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 0], [1, 0, 4.4], [0, 1, 0]])
          >>> l.z_max()
          4.4
        """

        return max(self._z) if self.num_pts() > 0 else None

    def z_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean of z values.

        :return: the optional maximum of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 2], [1, 0, 4], [0, 1, 6]])
          >>> l.z_mean()
          4.0
        """

        return np.mean(self._z) if self.num_pts() > 0 else None

    def z_var(self) -> Optional[numbers.Real]:
        """
        Optional variance of z values.

        :return: the optional variance of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 2], [1, 0, 2], [0, 1, 2]])
          >>> l.z_var()
          0.0
        """

        return np.var(self._z) if self.num_pts() > 0 else None

    def z_std(self) -> Optional[numbers.Real]:
        """
        Optional standard deviation of z values.

        :return: the optional standard deviation of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = Points.fromPointList([[0, 0, 2], [1, 0, 2], [0, 1, 2]])
          >>> l.z_std()
          0.0
        """

        return np.std(self._z) if self.num_pts() > 0 else None

    def remove_coincident_points(self) -> Optional['Points']:
        """
        Remove coincident successive points

        :return: Line instance
        :rtype: Optional[Points]
        """

        if self.num_pts() == 0:
            return

        new_line = Points(
            pts=[self.pt(0)]
        )

        for ndx in range(1, self.num_pts()):
            if not self.pt(ndx).isCoinc3D(new_line.pt(-1)):
                new_line.add_pt(self.pt(ndx))

        return new_line

    def as_segments(self):
        """
        Convert to a list of segments.

        :return: list of Segment objects
        """

        pts_pairs = zip(self.pts()[:-1], self.pts()[1:])

        segments = [Segment(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    '''
    def densify_2d_line(self, sample_distance) -> 'Points':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line instance
        """

        if sample_distance <= 0.0:
            raise Exception(f"Sample distance must be positive. {sample_distance} received")

        segments = self.as_segments()

        densified_line_list = [segment.densify2d_asLine(sample_distance) for segment in segments]

        densifyied_multiline = MultiLine(densified_line_list)

        densifyied_line = densifyied_multiline.to_line()

        densifyied_line_wo_coinc_pts = densifyied_line.remove_coincident_points()

        return densifyied_line_wo_coinc_pts
    '''

    def join(self, another) -> 'Points':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Points(self.pts() + another.pts())

    def length_3d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).dist3DWith(self.pt(ndx + 1))
        return length

    '''
    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).dist2DWith(self.pt(ndx + 1))
        return length
    '''

    def step_delta_z(self) -> List[numbers.Real]:
        """
        Return the difference in elevation between consecutive points:
        z[ndx+1] - z[ndx]

        :return: a list of height differences.
        :rtype: list of floats.
        """

        delta_z = [0.0]

        for ndx in range(1, self.num_pts()):
            delta_z.append(self.pt(ndx).z - self.pt(ndx - 1).z)

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
            length = self.pt(ndx).dist3DWith(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list

    '''
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
            length = self.pt(ndx).dist2DWith(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list
    '''

    def incremental_length_3d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 3D segment lengths.

        :return: accumulated 3D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_3d()))

    '''
    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))
    '''

    def reversed(self) -> 'Points':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        pts = [pt.clone() for pt in self.pts()]
        pts.reverse()

        return Points(
            pts=pts
        )

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

    def extremes_distance_3d(self) -> numbers.Real:
        """
        Calculate the 3D distance between start and end points.

        :return: the 3D distance between start and end points
        :rtype: numbers.Real
        """

        return self.end_pt().dist3DWith(self.start_pt())

    '''
    def extremes_distance_2d(self) -> numbers.Real:
        """
        Calculate the 2D distance between start and end points.

        :return: the 2D distance between start and end points
        """

        return self.end_pt().dist2DWith(self.start_pt())
    '''

    def isClosed_3d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 3D-closed.

        :param tolerance: the tolerance for considering the line closed
        :type tolerance: numbers.Real
        :return: whether the line is to be considered 3D-closed
        :rtype: bool
        """

        return self.end_pt().isCoinc3D(self.start_pt(), tolerance=tolerance)

    '''
    def isClosed_2d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 2D-closed.

        :param tolerance: the tolerance for considering the line closed
        :return: whether the line is to be considered 2D-closed
        """

        return self.end_pt().isCoinc2D(self.start_pt(), tolerance=tolerance)
    '''

    def walk_backward(self) -> 'Points':
        """
        Create a new line by walking the line backward from the last point up to the first and thus closing it.

        :return: a closed line with zero area
        :rtype: 'Line'
        """

        return Points(self.pts() + self.reversed().pts()[1:])

    def clone(self) -> 'Points':
        """
        Clone a line.

        :return: the cloned line
        :rtype: Points
        """

        return Points(self.pts())

    '''
    def close_2d(self) -> 'Points':
        """
        Return a line that is 2D-closed.

        :return: a 2D-closed line
        :rtype: Points
        """

        line = self.clone()

        if not line.isClosed_2d():

            line.add_pt(line.start_pt())

        return line
    '''

    def close_3d(self) -> 'Points':
        """
        Return a line that is 3D-closed.

        :return: a 3D-closed line
        :rtype: Points
        """

        line = self.clone()

        if not line.isClosed_3d():

            line.add_pt(line.start_pt())

        return line


def line_from_shapely(
        shapely_geom: LineString,
        epsg_code: numbers.Integral
) -> Points:
    # Side effects: none
    """
    Create a Line instance from a shapely Linestring instance.

    :param shapely_geom: the shapely input LineString instance
    :type shapely_geom: shapely.geometry.linestring.LineString
    :param epsg_code: the EPSG code of the LineString instance
    :type epsg_code: numbers.Integral
    :return: the converted Line instance
    :rtype: Points
    """

    x_array, y_array = shapely_geom.xy

    return Points.fromArrays(
        x_array,
        y_array,
        epsg_cd=epsg_code
    )


def line_to_shapely(
        src_line: Points
) -> LineString:
    """
    Create a shapely.LineString instance from a Line one.

    :param src_line: the source line to convert to the shapely format
    :type src_line: Points
    :return: the shapely LineString instance and the EPSG code
    :rtype: Tuple[LineString, numbers.Integral]
    """

    return LineString(src_line.xy_zipped()), src_line.epsg_code()


class Lines(list):
    """
    Collection of lines.

    """

    def __init__(self,
                 lines: Optional[List[Points]] = None
                 ):

        if lines:

            check_type(lines, "Lines", List)
            for line in lines:
                check_type(line, "Line", Points)
            first_line = lines[0]
            for line in lines[1:]:
                check_crs(first_line, line)

            super(Lines, self).__init__(lines)

        else:

            super(Lines, self).__init__()

    def append(self,
               item: Points
               ) -> None:

        check_type(item, "Line", Points)
        if len(self) > 0:
            check_crs(self[0], item)

        super(Lines, self).append(item)


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


class JoinTypes(Enum):
    """
    Enumeration for Line and Segment type.
    """

    START_START = 1  # start point coincident with start point
    START_END   = 2  # start point coincident with end point
    END_START   = 3  # end point coincident with start point
    END_END     = 4  # end point coincident with end point


def analizeJoins(first: Union[Points, Segment], second: Union[Points, Segment]) -> List[Optional[JoinTypes]]:
    """
    Analyze join types between two lines/segments.

    :param first: a line or segment.
    :type first: Line or Segment.
    :param second: a line or segment.
    :param second: Line or Segment.
    :return: a list of join types.
    :rtype: List[Optional[JoinTypes]].

    Examples:
      >>> first = Segment(Point(x=0,y=0), Point(x=1,y=0))
      >>> second = Segment(Point(x=1,y=0), Point(x=0,y=0))
      >>> analizeJoins(first, second)
      [<JoinTypes.START_END: 2>, <JoinTypes.END_START: 3>]
      >>> first = Segment(Point(x=0,y=0), Point(x=1,y=0))
      >>> second = Segment(Point(x=2,y=0), Point(x=3,y=0))
      >>> analizeJoins(first, second)
      []
    """

    join_types = []

    if first.start_pt.isCoinc3D(second.start_pt):
        join_types.append(JoinTypes.START_START)

    if first.start_pt.isCoinc3D(second.end_pt):
        join_types.append(JoinTypes.START_END)

    if first.end_pt.isCoinc3D(second.start_pt):
        join_types.append(JoinTypes.END_START)

    if first.end_pt.isCoinc3D(second.end_pt):
        join_types.append(JoinTypes.END_END)

    return join_types


