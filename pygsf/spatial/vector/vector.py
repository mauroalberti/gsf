# -*- coding: utf-8 -*-


from typing import Optional

from math import degrees, cos, acos

from typing import Tuple

from pygsf.mathematics.arrays import *


isfinite = np.isfinite
array = np.array


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D (space)
    """

    def __init__(self, x: [int, float], y: [int, float], z: [int, float]):
        """
        Construct a Point instance.
        """

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, (int, float)), vals)):
            raise VectorInputException("Input values must be integer of float")
        elif not all(map(isfinite, vals)):
            raise VectorInputException("Input values must be finite")
        else:
            self._a = array(vals, dtype=np.float64)

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self, another: 'Point') -> Optional[bool]:
        """
        Return True if objects are equal.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, -1)
          False
        """

        if not isinstance(another, Point):
            raise VectorInputException("Variables must be of the same type")
        else:
            return all([
                self.x == another.x,
                self.y == another.y,
                self.z == another.z])

    def __ne__(self, another: 'Point') -> Optional[bool]:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
        """

        if not isinstance(another, Point):
            return None
        else:
            return not (self == another)

    @property
    def a(self) -> 'numpy.array':
        """
        Return a copy of the object inner array.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7).a
          array([4., 3., 7.])
        """

        return np.copy(self._a)

    @property
    def x(self) -> float:
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self) -> float:
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.a[1]

    @property
    def z(self) -> float:
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.a[2]

    def toXYZ(self) -> Tuple[float, float, float]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> 'numpy.array':
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> Point(1, 2, 3).toArray()
          array([1., 2., 3.])
        """

        return self.a

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return self.__class__(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return self.__class__(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return self.__class__(0.0, self.y, self.z)

    @property
    def len3D(self) -> float:
        """
        Spatial distance of the point from the axis origin.

        :return: distance
        :rtype: float

        Examples:
          >>> Point(4.0, 3.0, 0.0).len3D
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def len2D(self) -> float:
        """
        2D distance of the point from the axis origin.

        Example:
          >>> Point(3, 4, 0).len2D
          5.0
          >>> Point(12, 5, 3).len2D
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    def deltaX(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.x - self.x

    def deltaY(self, another: 'Point') -> Optional[float]:
        """
        Delta between y components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.y - self.y

    def deltaZ(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.z - self.z

    def dist3DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate Euclidean spatial distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist3DWith(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist2DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate horizontal (2D) distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist2DWith(Point(4., 5., 7.))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self, scale_factor: [int, float]) -> Optional['Point']:
        """
        Create a scaled object.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(np.nan) is None
          True
          >>> Point(1, 0, 1).scale(np.inf) is None
          True
        """

        if not isinstance(scale_factor, (int, float)):
            return None
        elif not isfinite(scale_factor):
            return None
        else:
            x, y, z = arr2tuple(self.a * scale_factor)
            return self.__class__(x, y, z)

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def isCoinc(self, another: 'Point', tolerance: [int, float] = MIN_SEPARATION_THRESHOLD) -> Optional[bool]:
        """
        Check spatial coincidence of two points

        Example:
          >>> Point(1., 0., -1.).isCoinc(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc(Point(1., 0., 0.))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4), tolerance=np.nan) is None
          True
        """

        if not isinstance(another, Point):
            return None
        elif not isinstance(tolerance, (int, float)):
            return None
        elif not isfinite(tolerance):
            return None
        else:
            distance_2d = self.dist2DWith(another)
            if np.isnan(distance_2d) or distance_2d > tolerance:
                return False
            else:
                distance_3d = self.dist3DWith(another)
                if np.isnan(distance_3d) or distance_3d > tolerance:
                    return False
                else:
                    return True

    def shift(self, sx: float, sy: float, sz: float) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
          >>> Point(1, 2, -1).shift(0.5, np.nan, 1.5) is None
          True
       """

        vals = [sx, sy, sz]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            return None
        elif not all(map(isfinite, vals)):
            return None
        else:
            return self.__class__(self.x + sx, self.y + sy, self.z + sz)

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


class Vect(Point):
    """
    Cartesian 3D vector.
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    z axis -> Up
    """

    def __init__(self, x: [int, float], y: [int, float], z: [int, float]):
        """
        Vect constructor.

        Example;
          >>> Vect(1, 0, 1)
          Vect(1.0000, 0.0000, 1.0000)
          >>> Vect(1, np.nan, 1)
          Traceback (most recent call last):
          ...
          VectorInputException: Input values must be finite
          >>> Vect(1, 0, np.inf)
          Traceback (most recent call last):
          ...
          VectorInputException: Input values must be finite
          >>> Vect(0, 0, 0)
          Vect(0.0000, 0.0000, 0.0000)
        """

        super().__init__(x, y, z)


    def __repr__(self) -> str:

        return "Vect({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __add__(self, another: 'Vect') -> 'Vect':
        """
        Sum of two vectors.

        Example:
          >>> Vect(1, 0, 0) + Vect(0, 1, 1)
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(1, 1, 1) + Vect(-1, -1, -1)
          Vect(0.0000, 0.0000, 0.0000)
        """

        x, y, z = arr2tuple(self.a + another.a)
        return self.__class__(x, y, z)

    def __sub__(self, another: 'Vect') -> 'Vect':
        """Return object difference

        Example:
          >>> Vect(1., 1., 1.) - Vect(1., 1., 1.)
          Vect(0.0000, 0.0000, 0.0000)
          >>> Vect(1., 1., 3.) - Vect(1., 1., 2.2)
          Vect(0.0000, 0.0000, 0.8000)
        """

        x, y, z = arr2tuple(self.a - another.a)
        return self.__class__(x, y, z)

    @property
    def isAlmostZero(self) -> bool:
        """
        Check if the Vect instance length is near zero.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isAlmostZero
          False
          >>> Vect(0.0, 0.0, 0.0).isAlmostZero
          True
        """

        return are_close(self.len3D, 0)

    @property
    def isAlmostUnit(self) -> bool:
        """
        Check if the Vect instance length is near unit.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isAlmostUnit
          False
          >>> Vect(0.0, 1.0, 0.0).isAlmostUnit
          True
        """

        return are_close(self.len3D, 1)

    @property
    def isValid(self) -> bool:
        """
        Check if the Vect instance components are not all valid and the xyz not all zero-valued.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isValid
          True
          >>> Vect(0.0, 0.0, 0.0).isValid
          False
        """

        return not self.isAlmostZero

    def versor(self) -> Optional['Vect']:
        """
        Calculate versor in xyz space.

        Example:
          >>> Vect(5, 0, 0).versor()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> Vect(0, 0, 0).versor() is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.scale(1.0 / self.len3D)

    def versor2D(self) -> Optional['Vect']:
        """
        Create 2D versor version of the current vector

        :return: unit vector

        Example:
          >>> Vect(7, 0, 10).versor2D()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, 10).versor2D() is None
          True
        """

        vXY = self.pXY()
        if vXY.isValid:
            return self.pXY().versor()
        else:
            return None

    @property
    def isUpward(self) -> Optional[bool]:
        """
        Check that a vector is upward-directed.

        :return: boolean

        Example:
          >>> Vect(0,0,1).isUpward
          True
          >>> Vect(0,0,-0.5).isUpward
          False
          >>> Vect(1, 3, 0).isUpward
          False
          >>> Vect(0, 0, 0).isUpward is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.z > 0.0

    @property
    def isDownward(self) -> Optional[bool]:
        """
        Check that a vector is downward-directed.

        :return: boolean

        Example:
          >>> Vect(0,0,1).isDownward
          False
          >>> Vect(0,0,-0.5).isDownward
          True
          >>> Vect(1, 3, 0).isDownward
          False
          >>> Vect(0, 0, 0).isDownward is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.z < 0.0

    def upward(self) -> Optional['Vect']:
        """
        Calculate a new upward-pointing vector.

        Example:
          >>> Vect(1, 1, 1).upward()
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(-1, -1, -1).upward()
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(0, 0, 0).upward() is None
          True
        """

        if not self.isValid:
            return None
        elif self.z < 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def downward(self) -> Optional['Vect']:
        """
        Calculate a new vector downward-pointing.

        Example:
          >>> Vect(1, 1, 1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(-1, -1, -1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(0, 0, 0).downward() is None
          True
        """

        if not self.isValid:
            return None
        elif self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def vDot(self, another: 'Vect') -> float:
        """
        Vector scalar multiplication.

        Examples:
          >>> Vect(1, 0, 0).vDot(Vect(1, 0, 0))
          1.0
          >>> Vect(1, 0, 0).vDot(Vect(0, 1, 0))
          0.0
          >>> Vect(1, 0, 0).vDot(Vect(-1, 0, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y + self.z * another.z

    def angleCos(self, another: 'Vect') -> Optional[float]:
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect(1,0,0).angleCos(Vect(0,0,1))
          0.0
          >>> Vect(1,0,0).angleCos(Vect(-1,0,0))
          -1.0
          >>> Vect(1,0,0).angleCos(Vect(1,0,0))
          1.0
          >>> Vect(0, 0, 0).angleCos(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).angleCos(Vect(0,0,0)) is None
          True
        """

        if not (self.isValid and another.isValid):
            return None
        else:
            val = self.vDot(another) / (self.len3D * another.len3D)
            if val > 1.0:
                return 1.0
            elif val < -1.0:
                return -1.0
            else:
                return val

    def angle(self, another: 'Vect') -> Optional[float]:
        """
        Calculate angle between two vectors, as degrees
        in 0° - 180° range.

        Example:
          >>> Vect(1, 0, 0).angle(Vect(0, 0, 1))
          90.0
          >>> Vect(1, 0, 0).angle(Vect(-1, 0, 0))
          180.0
          >>> Vect(0, 0, 1).angle(Vect(0, 0, -1))
          180.0
          >>> Vect(1, 1, 1).angle(Vect(1, 1,1 ))
          0.0
          >>> Vect(0, 0, 0).angle(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).angle(Vect(0,0,0)) is None
          True
        """

        if not (self.isValid and another.isValid):
            return None
        else:
            return degrees(acos(self.angleCos(another)))

    def isAlmostParallel(self, another: 'Vect', angle_tolerance: [int, float] = VECTOR_ANGLE_THRESHOLD) -> Optional[
        bool]:
        """
        Check that two Vect are sub-parallel,

        :param another: a Vect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(1, 0, 0))
          True
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(0, 0, 1))
          False
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(-1, 0, 0))
          False
          >>> Vect(0, 0, 0).isAlmostParallel(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(0,0,0)) is None
          True
        """

        if not isinstance(another, Vect):
            return None
        elif not (self.isValid and another.isValid):
            return None
        elif not isinstance(angle_tolerance, (int, float)):
            return None
        elif not isfinite(angle_tolerance):
            return None
        else:
            return self.angle(another) <= angle_tolerance

    def isSubOrthogonal(self, another: 'Vect', angle_tolerance: [int, float] = VECTOR_ANGLE_THRESHOLD) -> Optional[
        bool]:
        """
        Check whether two vectors are sub-orhogonal.

        :param another: a second Vect instance
        :param angle_tolerance: the tolerance angle, in decimal degrees
        :return: Boolean

        Example:
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 0))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 1))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 0.9999999999999, 0))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal((Vect(0, 0, 0))) is None
          True
        """

        if not isinstance(another, Vect):
            return None
        elif not (self.isValid and another.isValid):
            return None
        elif not isinstance(angle_tolerance, (int, float)):
            return None
        elif not isfinite(angle_tolerance):
            return None
        else:
            return are_close(0, self.angleCos(another), atol=cos(angle_tolerance))

    def vCross(self, another: 'Vect') -> 'Vect':
        """
        Vector product (cross product).

        Examples:
          >>> Vect(1, 0, 0).vCross(Vect(0, 1, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> Vect(1, 0, 0).vCross(Vect(1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
          >>> (Vect(1, 0, 0).vCross(Vect(-1, 0, 0))).isAlmostZero
          True
        """

        x, y, z = arr2tuple(np.cross(self.a[:3], another.a[:3]))
        return Vect(x, y, z)

    def byMatrix(self, array3x3: 'np.array') -> 'Vect':
        """
        Matrix multiplication of a vector.

        """

        x, y, z = arr2tuple(array3x3.dot(self.a))
        return Vect(x, y, z)


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with PPlane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a, b, c, d):

        self._a = float(a)
        self._b = float(b)
        self._c = float(c)
        self._d = float(d)

    @property
    def a(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).a
          1.0
        """

        return self._a

    @property
    def b(self):
        """
        Return b coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 4, 0, 2).b
          4.0
        """

        return self._b

    @property
    def c(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 5.4, 2).c
          5.4
        """

        return self._c

    @property
    def d(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).d
          2.0
        """

        return self._d

    @property
    def v(self):
        """
        Return coefficients of a Plane instance.

        Example:
          >>> CPlane(1, 1, 7, -4).v
          (1.0, 1.0, 7.0, -4.0)
        """
        return self.a, self.b, self.c, self.d

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3):
        """
        Create a Plane from three given Point instances.

        Example:
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          CPlane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = array([[pt1.y, pt1.z, 1],
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1]])

        matr_b = - array([[pt1.x, pt1.z, 1],
                             [pt2.x, pt2.z, 1],
                             [pt3.x, pt3.z, 1]])

        matr_c = array([[pt1.x, pt1.y, 1],
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1]])

        matr_d = - array([[pt1.x, pt1.y, pt1.z],
                             [pt2.x, pt2.y, pt2.z],
                             [pt3.x, pt3.y, pt3.z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def __repr__(self):

        return "CPlane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

    def normVersor(self):
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane(0, 0, 5, -2).normVersor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> CPlane(0, 7, 0, 5).normVersor()
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a, self.b, self.c).versor()

    def toPoint(self):
        """
        Returns a point lying in the plane (non-unique solution).

        Examples:
          >>> CPlane(0, 0, 1, -1).toPoint()
          Point(0.0000, 0.0000, 1.0000)
        """

        point = Point(*point_solution(array([[self.a, self.b, self.c]]),
                                      array([-self.d])))
        return point

    def intersVersor(self, another):
        """
        Return intersection versor for two intersecting planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersVersor(b)
        Vect(0.0000, -1.0000, 0.0000)
        """

        return self.normVersor().vCross(another.normVersor()).versor()

    def intersPoint(self, another):
        """
        Return point on intersection line (obviously non-unique solution)
        for two planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersPoint(b)
        Point(0.0000, 0.0000, 0.0000)
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point(x, y, z)

    def isPointInPlane(self, pt):
        """
          Check whether a point lie in a plane.

          >>> pl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
        """

        if abs(self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d) < MIN_SCALAR_VALUE:
            return True
        else:
            return False

    def angle(self, another):
        """
        Calculate angle (in degrees) between two planes.

        Examples:
          >>> CPlane(1,0,0,0).angle(CPlane(0,1,0,0))
          90.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,1,0))
          45.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,0,0))
          0.0
        """

        angle_degr = self.normVersor().angle(another.normVersor())
        if abs(angle_degr) < MIN_ANGLE_DEGR_VALUE:
            angle_degr = 0.0
        elif angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr

    def isAlmostParallel(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two Plane are sub-parallel

        :param another: a Plane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> CPlane(1,0,0,0).isAlmostParallel(CPlane(1,0,0,0))
          True
          >>> CPlane(1,0,0,0).isAlmostParallel(CPlane(1,0,1,0))
          False
        """

        return self.angle(another) < angle_tolerance


class VectorInputException(Exception):
    """
    Exception for geometric input.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
