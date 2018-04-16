# -*- coding: utf-8 -*-


from __future__ import division

from math import sqrt, sin, cos, radians, acos, atan, isnan
import numpy as np

from typing import Dict, Tuple, List

from .default_parameters import *
from .mathematics import are_close
from .geometry_utils import *
from .arrays import point_solution, arrays_are_close, arr2tuple


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D (space)
    """

    def __init__(self, x, y, z):
        """
        Construct a Point instance.
        """

        if any(map(lambda val: not np.isfinite(val), [x, y, z])):
            raise GeomInputException("All spatial coordinates must be finite")

        self._a = np.array([x, y, z], dtype=np.float64)

    def __repr__(self):

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __sub__(self, another):
        """Return point difference

        Example:
          >>> Point(1., 1., 1.) - Point(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000)
          >>> Point(1., 1., 3.) - Point(1., 1., 2.2)
          Point(0.0000, 0.0000, 0.8000)
        """

        x, y, z = arr2tuple(self.a - another.a)
        return self.__class__(x, y, z)

    def __eq__(self, another):
        """
        Return True if points are equal.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
        """

        return all([
            self.x == another.x,
            self.y == another.y,
            self.z == another.z])

    def __ne__(self, another):
        """
        Return False if vectors are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
        """

        return not (self == another)

    @property
    def a(self):
        """
        Return values as array

        Example:
          >>> arrays_are_close(Point(1, 0, 0).a, np.array([ 1.,  0.,  0.]), equal_nan=True)
          True
        """

        return self._a

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.a[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.a[2]

    def xyz(self) -> Tuple[float, float, float]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).xyz()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def pXY(self):
        """
        Projection of point on the x-y plane

        :return: projected Point instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return Point(self.x, self.y, 0.0)

    def pXZ(self):
        """
        Projection of point on the x-z plane

        :return: projected Point instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return Point(self.x, 0.0, self.z)

    def pYZ(self):
        """
        Projection of point on the y-z plane

        :return: projected Point instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return Point(0.0, self.y, self.z)

    def len(self) -> float:
        """
        Spatial distance of the point form the axis origin.

        :return: distance
        :rtype: float

        Examples:
          >>> Point(4.0, 3.0, 0.0).len()
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def deltaX(self, another) -> float:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        return another.x - self.x

    def deltaY(self, another) -> float:
        """
        Delta between y components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        return another.y - self.y

    def deltaZ(self, another) -> float:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        return another.z - self.z

    def dist3D(self, another):
        """
        Calculate Euclidean spatial distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist3D(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1).dist3D(Point(4, 5, 1))
          5.0
        """

        return (self - another).len()

    def dist2D(self, another):
        """
        Calculate horizontal (2D) distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist2D(Point(4., 5., 7.))
          5.0
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self, scale_factor: float):
        """
        Create a scaled point.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Scale factor for vector must be finite
          >>> Point(1, 0, 1).scale(np.inf)
          Traceback (most recent call last):
          ...
          GeomInputException: Scale factor for vector must be finite
        """

        if not np.isfinite(scale_factor):
            raise GeomInputException("Scale factor for vector must be finite")

        x, y, z = arr2tuple(self.a * scale_factor)
        return self.__class__(x, y, z)

    def invert(self):
        """
        Create a new point with inverted direction.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def isCoinc(self, another, tolerance=MIN_SEPARATION_THRESHOLD):
        """
        Check spatial coincidence of two points

        Example:
          >>> Point(1., 0., -1.).isCoinc(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc(Point(1., 0., 0.))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4))
          True
        """

        distance_2d = self.dist2D(another)
        if np.isnan(distance_2d) or distance_2d > tolerance:
            return False

        distance_3d = self.dist3D(another)
        if np.isnan(distance_3d) or distance_3d > tolerance:
            return False

        return True

    def shift(self, sx=0.0, sy=0.0, sz=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz)

    def shiftWithVect(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector.

        Example:
          >>> Point(1, 2, 0).shiftWithVect(Vect(10, 5, 0))
          Point(11.0000, 7.0000, 0.0000)
          >>> Point(3, 1, 0.2).shiftWithVect(Vect(7, 5, 0))
          Point(10.0000, 6.0000, 0.2000)
        """

        return Point(self.x + displ_vect.x,
                     self.y + displ_vect.y,
                     self.z + displ_vect.z)

    def toArray(self):
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> Point(1, 2, 3).toArray()
          array([1., 2., 3.])
        """

        return self.a

    def asVect(self):
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

    def __init__(self, x=0, y=0, z=0):
        """
        Vect constructor.

        Example;
          >>> Vect(1, 0, 1)
          Vect(1.0000, 0.0000, 1.0000)
          >>> Vect(1, np.nan, 1)
          Traceback (most recent call last):
          ...
          GeomInputException: All spatial coordinates must be finite
          >>> Vect(1, 0, np.inf)
          Traceback (most recent call last):
          ...
          GeomInputException: All spatial coordinates must be finite
          >>> Vect(0, 0, 0)
          Vect(0.0000, 0.0000, 0.0000)
        """

        super().__init__(x, y, z)

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

    @property
    def len2D(self) -> float:
        """
        2D Vector magnitude.

        Example:
          >>> Vect(3, 4).len2D
          5.0
          >>> Vect(12, 5, 3).len2D
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    @property
    def len3D(self) -> float:
        """
        3D Vector magnitude.

        Example:
          >>> Vect(12, 0, 5).len3D
          13.0
          >>> Vect(3, 0, 4).len3D
          5.0
        """

        return self.len()

    @property
    def isAlmostZero(self) -> bool:
        """
        Check if the Vect instance lenght is near zero.

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
        Check if the Vect instance lenght is near unit.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isAlmostUnit
          False
          >>> Vect(0.0, 1.0, 0.0).isAlmostUnit
          True
        """

        return are_close(self.len3D, 1)

    def __repr__(self):

        return "Vect({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __add__(self, another):
        """
        Sum of two vectors.

        Example:
          >>> Vect(1, 0, 0) + Vect(0, 1, 1)
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(1, 1, 1) + Vect(-1, -1, -1)
          Vect(0.0000, 0.0000, 0.0000)
        """

        x, y, z = arr2tuple(self.a + another.a)
        return Vect(x, y, z)

    def versor(self):
        """
        Calculate versor in xyz space.

        Example:
          >>> Vect(5, 0, 0).versor()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> Vect(0, 0, 0).versor()
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not self.isValid:
            raise VectInvalidException("Vector must be isValid (not zero-valued)")
        return self.scale(1.0 / self.len3D)

    def flatten2D(self):
        """
        Create equivalent vector in xy space.

        Example:
          >>> Vect(5, 16, 43).flatten2D()
          Vect(5.0000, 16.0000, 0.0000)
        """

        return Vect(self.x, self.y, 0.0)

    def versor2D(self):
        """
        Create 2D versor equivalent of the current vector

        :return: Vector

        Example:
          >>> Vect(7, 0, 10).versor2D()
          Vect(1.0000, 0.0000, 0.0000)
        """

        return self.flatten2D().versor()

    def trend(self) -> Optional[float]:
        """
        Trend of a vector
        (degrees, clockwise from North, range 0°-360°)

        :return: an optional float value, representing trend, in decimal degrees.

        Examples:
          >>> Vect(1, 0, 0).trend()
          90.0
          >>> Vect(0, 1, 0).trend()
          0.0
          >>> Vect(1, 1, 0).trend()
          45.0
          >>> Vect(1, -1, 0).trend()
          135.0
          >>> Vect(0, -1, 0).trend()
          180.0
          >>> Vect(-1, -1, 0).trend()
          225.0
          >>> Vect(-1, 0, 0).trend()
          270.0
          >>> Vect(-1, 1, 0).trend()
          315.0
          >>> Vect(1, 1, 10).trend()
          45.0
          >>> Vect(0, 0, 0).trend() is None
          True
         """

        return angle_north_clock(self.x, self.y)

    def slope(self) -> Optional[float]:
        """
        Slope of a vector.
        Degrees, positive: downward-directed, negative: upward-dir., range -90°/90°

        :return: an optional float, representing the vector slope, in decimal degrees.

        Examples:
          >>> Vect(1, 0, -1).slope()
          45.0
          >>> Vect(1, 0, 1).slope()
          -45.0
          >>> Vect(0, 1, 0).slope()
          0.0
          >>> Vect(0, 0, 1).slope()
          -90.0
          >>> Vect(0, 0, -1).slope()
          90.0
          >>> Vect(0, 0, 0).slope() is None
          True
         """
        h = self.len2D()
        zv = self.z
        sl = slope(h, abs(zv))

        if sl is None:
            return None
        else:
            if zv <= 0.0:
                return sl
            else:
                return -sl

    @property
    def isUpward(self):
        """
        Check that a vector is upward-directed.
         
        Example:
          >>> Vect(0,0,1).isUpward
          True
          >>> Vect(0,0,-0.5).isUpward
          False
        """

        return self.z > 0.0

    @property
    def isDownward(self):
        """
        Check that a vector is downward-directed.

        Example:
          >>> Vect(0,0,1).isDownward
          False
          >>> Vect(0,0,-0.5).isDownward
          True
        """

        return self.z < 0.0

    def upward(self):
        """
        Calculate a new upward-pointing vector.

        Example:
          >>> Vect(1, 1, 1).upward()
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(-1, -1, -1).upward()
          Vect(1.0000, 1.0000, 1.0000)
        """

        if self.z < 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def downward(self):
        """
        Calculate a new vector downward-pointing.

        Example:
          >>> Vect(1, 1, 1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(-1, -1, -1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
        """

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    @property
    def slope(self):
        """
        Slope of a vector expressed as degrees.
        Positive when vector is downward pointing or horizontal,
        negative when upward pointing.

        Example:
          >>> Vect(1, 0, -1).slope
          45.0
          >>> Vect(1, 0, 1).slope
          -45.0
          >>> Vect(0, 1, 0).slope
          0.0
          >>> Vect(0, 0, 1).slope
          -90.0
          >>> Vect(0, 0, -1).slope
          90.0
          >>> Vect(0, 0, 0).slope
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not self.isValid:
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        hlen = self.len2D
        if hlen == 0.0:
            if self.z > 0.:
                return -90.
            elif self.z < 0.:
                return 90.
            else:
                raise Exception("Zero-valued vector")
        else:
            slope = - degrees(atan(self.z / self.len2D))
            if abs(slope) > MIN_SCALAR_VALUE:
                return slope
            else:
                return 0.

    def asGVect(self):
        """
        Calculate the geological vector parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(1, 1, 1).asGVect()
          GVect(045.00, -35.26)
          >>> Vect(0, 1, 1).asGVect()
          GVect(000.00, -45.00)
          >>> Vect(1, 0, 1).asGVect()
          GVect(090.00, -45.00)
          >>> Vect(0, 0, 1).asGVect()
          GVect(000.00, -90.00)
          >>> Vect(0, 0, -1).asGVect()
          GVect(000.00, +90.00)
          >>> Vect(-1, 0, 0).asGVect()
          GVect(270.00, +00.00)
          >>> Vect(0, -1, 0).asGVect()
          GVect(180.00, +00.00)
          >>> Vect(-1, -1, 0).asGVect()
          GVect(225.00, +00.00)
          >>> Vect(0, 0, 0).asGVect()
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not self.isValid:
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        plunge = self.slope  # upward pointing -> negative value, downward -> positive

        unit_vect = self.versor()
        if unit_vect.y == 0. and unit_vect.x == 0:
            trend = 0.
        else:
            trend = (90. - degrees(atan2(unit_vect.y, unit_vect.x))) % 360.

        return GVect(trend, plunge)

    def asGAxis(self):
        """
        Calculate the geological axis parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(0, 1, 1).asGAxis()
          GAxis(000.00, -45.00)
          >>> Vect(1, 0, 1).asGAxis()
          GAxis(090.00, -45.00)
          >>> Vect(0, 0, 1).asGAxis()
          GAxis(000.00, -90.00)
          >>> Vect(0, 0, -1).asGAxis()
          GAxis(000.00, +90.00)
          >>> Vect(-1, 0, 0).asGAxis()
          GAxis(270.00, +00.00)
          >>> Vect(0, -1, 0).asGAxis()
          GAxis(180.00, +00.00)
          >>> Vect(-1, -1, 0).asGAxis()
          GAxis(225.00, +00.00)
          >>> Vect(0, 0, 0).asGAxis()
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not self.isValid:
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        return self.asGVect().asGAxis()

    def vDot(self, another):
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

    def angleCos(self, another):
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect(1,0,0).angleCos(Vect(0,0,1))
          0.0
          >>> Vect(1,0,0).angleCos(Vect(-1,0,0))
          -1.0
          >>> Vect(1,0,0).angleCos(Vect(1,0,0))
          1.0
          >>> Vect(0, 0, 0).angleCos(Vect(1,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
          >>> Vect(1, 0, 0).angleCos(Vect(0,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not (self.isValid and another.isValid):
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        val = self.vDot(another) / (self.len3D * another.len3D)
        if val > 1.0:
            return 1.0
        elif val < -1.0:
            return -1.0
        else:
            return val

    def angle(self, another):
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
          >>> Vect(0, 0, 0).angle(Vect(1,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
          >>> Vect(1, 0, 0).angle(Vect(0,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not (self.isValid and another.isValid):
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        return degrees(acos(self.angleCos(another)))

    def isAlmostParallel(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
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
          >>> Vect(0, 0, 0).isAlmostParallel(Vect(1,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(0,0,0))
          Traceback (most recent call last):
          ...
          VectInvalidException: Vector must be isValid (not zero-valued)
        """

        if not (self.isValid and another.isValid):
            raise VectInvalidException("Vector must be isValid (not zero-valued)")

        return self.angle(another) <= angle_tolerance

    def isSubOrthogonal(self, another):
        """
        Check whether two vectors are sub-orhogonal.

        :param another:
        :return: Boolean

        Example:
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 0))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 1))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 0.9999999999999, 0))
          True
        """

        return are_close(0, self.angleCos(another))

    def vCross(self, another):
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

    def byMatrix(self, array3x3):
        """
        Matrix multiplication of a vector.

        """

        x, y, z = arr2tuple(array3x3.dot(self.a))
        return Vect(x, y, z)


class GVect(object):
    """
    Geological vector.
    Defined by trend and plunge (both in degrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
    """

    def __init__(self, trend, plunge, is_axis=False):
        """
        Geological vector constructor.
        trend: Trend range: [0.0, 360.0[ clockwise, from 0 (North)
        plunge: Plunge: [-90.0, 90.0],
        negative value: upward pointing is_axis, positive values: downward is_axis;
            
        Example:
          >>> GVect(120.2, -27.4)
          GVect(120.20, -27.40)
          >>> GVect(315.0, -0.0)
          GVect(315.00, -00.00)
          >>> GVect(23, 40.0)
          GVect(023.00, +40.00)
       """

        self._trend = float(trend) % 360.0
        self._plunge = float(plunge)
        self._axis = is_axis

    @property
    def tr(self):
        """
        Return trend of the geological direction.
        Range is [0, 360[

        Example:
          >>> GVect(420, -17).tr
          60.0
          >>> GVect(-20, 49).tr
          340.0
        """

        return self._trend

    @property
    def pl(self):
        """
        Return plunge of the geological direction.
        Range is [-90, 90]

        Example:
          >>> GVect(420, -17).pl
          -17.0
        """

        return self._plunge

    @property
    def ax(self):
        """
        Return whether is a geological axis.

        Example:
          >>> GVect(420, -17).ax
          False
          >>> GVect(420, -17, is_axis=True).ax
          True
        """

        return self._axis

    @property
    def tp(self):
        """
        Return trend and plunge of the geological direction.

        Example:
          >>> GVect(-90, -45).tp
          (270.0, -45.0)
        """

        return self.tr, self.pl

    @property
    def tpa(self):
        """
        Return trend and plunge of the geological direction, as well as if known movement sense..

        Example:
          >>> GVect(-90, -45).tpa
          (270.0, -45.0, False)
          >>> GVect(-90, -45, is_axis=True).tpa
          (270.0, -45.0, True)
        """

        return self.tr, self.pl, self.ax

    @property
    def tpoa(self):
        """
        Return trend and plunge of the geological direction, as well as the opposite of the known/unknown movement sense..

        Example:
          >>> GVect(-90, -45).tpoa
          (270.0, -45.0, True)
          >>> GVect(-90, -45, is_axis=True).tpoa
          (270.0, -45.0, False)
        """

        return self.tr, self.pl, not self.ax

    @property
    def pt(self):
        """
        Return plunge and trend of the geological direction.

        Example:
          >>> GVect(-90, -45).pt
          (-45.0, 270.0)
        """

        return self.pl, self.tr

    def hasDirection(self):
        """
        GVect has direction

        """

        return self._axis

    def __repr__(self):

        if self.ax:
            return "GAxis({:06.2f}, {:+06.2f})".format(*self.tp)
        else:
            return "GVect({:06.2f}, {:+06.2f})".format(*self.tp)

    @property
    def colatNorth(self) -> float:
        """
        Calculates the colatitude from the North (top).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: float

        Examples:
          >>> GVect(320, 90).colatNorth
          180.0
          >>> GVect(320, 45).colatNorth
          135.0
          >>> GVect(320, 0).colatNorth
          90.0
          >>> GVect(320, -45).colatNorth
          45.0
          >>> GVect(320, -90).colatNorth
          0.0
        """

        return plng2colatTop(self.pl)

    @property
    def colatSouth(self) -> float:
        """
        Calculates the colatitude from the South (bottom).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: float

        Examples:
          >>> GVect(320, 90).colatSouth
          0.0
          >>> GVect(320, 45).colatSouth
          45.0
          >>> GVect(320, 0).colatSouth
          90.0
          >>> GVect(320, -45).colatSouth
          135.0
          >>> GVect(320, -90).colatSouth
          180.0
        """

        return plng2colatBottom(self.pl)

    def copy(self):
        """
        Return a copy of the GVect instance.
        
        Example:
          >>> GVect(10, 20).copy()
          GVect(010.00, +20.00)
        """

        return GVect(*(self.tpa))

    def opposite(self):
        """
        Return the opposite GVect.
        
        Example:
          >>> GVect(0, 30).opposite()
          GVect(180.00, -30.00)
          >>> GVect(315, 10).opposite()
          GVect(135.00, -10.00)
          >>> GVect(135, 0).opposite()
          GVect(315.00, -00.00)
          >>> GVect(0, 30, True).opposite()
          GAxis(180.00, -30.00)
          >>> GVect(315, 10, True).opposite()
          GAxis(135.00, -10.00)
          >>> GVect(135, 0, True).opposite()
          GAxis(315.00, -00.00)
        """

        trend = (self.tr + 180.) % 360.
        plunge = -self.pl

        return GVect(trend, plunge, self.ax)

    def mirrorHoriz(self):
        """
        Return the mirror GVect using a horizontal plane.

        Example:
          >>> GVect(0, 30).mirrorHoriz()
          GVect(000.00, -30.00)
          >>> GVect(315, 10).mirrorHoriz()
          GVect(315.00, -10.00)
          >>> GVect(135, 0).mirrorHoriz()
          GVect(135.00, -00.00)
        """

        trend = self.tr
        plunge = -self.pl

        return GVect(trend, plunge, self.ax)

    def versor(self):
        """
        Return the unit vector corresponding to the geological vector.

        Examples:
          >>> GVect(0, 90).versor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> GVect(0, -90).versor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> GVect(90, 90).versor()
          Vect(0.0000, 0.0000, -1.0000)
          
        """

        north_coord = cos(radians(self.pl)) * cos(radians(self.tr))
        east_coord = cos(radians(self.pl)) * sin(radians(self.tr))
        down_coord = sin(radians(self.pl))

        return Vect(east_coord, north_coord, -down_coord)

    @property
    def isUpward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> GVect(10, 15).isUpward
          False
          >>> GVect(257.4, 0.0).isUpward
          False
          >>> GVect(90, -45).isUpward
          True
        """

        return self.pl < 0.0

    @property
    def isDownward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> GVect(10, 15).isDownward
          True
          >>> GVect(257.4, 0.0).isDownward
          False
          >>> GVect(90, -45).isDownward
          False
        """

        return self.pl > 0.0

    def upward(self):
        """
        Return upward-point geological vector.

        Examples:
          >>> GVect(90, -45).upward().isAlmostParallel(GVect(90.0, -45.0))
          True
          >>> GVect(180, 45).upward().isAlmostParallel(GVect(0.0, -45.0))
          True
          >>> GVect(0, 0).upward().isAlmostParallel(GVect(0.0, 0.0))
          True
          >>> GVect(0, 90).upward().isAlmostParallel(GVect(180.0, -90.0))
          True
          >>> GVect(90, -45).upward().isAlmostParallel(GVect(90.0, -35.0))
          False
          >>> GVect(180, 45).upward().isAlmostParallel(GVect(10.0, -45.0))
          False
          >>> GVect(0, 0).upward().isAlmostParallel(GVect(170.0, 0.0))
          False
          >>> GVect(0, 90).upward().isAlmostParallel(GVect(180.0, -80.0))
          False
        """

        if not self.isDownward:
            return self.copy()
        else:
            return self.opposite()

    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> GVect(90, -45).downward().isAlmostParallel(GVect(270.0, 45.0))
          True
          >>> GVect(180, 45).downward().isAlmostParallel(GVect(180.0, 45.0))
          True
          >>> GVect(0, 0).downward().isAlmostParallel(GVect(180.0, 0.0))
          False
          >>> GVect(0, 90).downward().isAlmostParallel(GVect(0.0, 90.0))
          True
          >>> GVect(90, -45).downward().isAlmostParallel(GVect(270.0, 35.0))
          False
          >>> GVect(180, 45).downward().isAlmostParallel(GVect(170.0, 45.0))
          False
          >>> GVect(0, 0).downward().isAlmostParallel(GVect(180.0, 10.0))
          False
          >>> GVect(0, 90).downward().isAlmostParallel(GVect(0.0, 80.0))
          False
        """

        if not self.isUpward:
            return self.copy()
        else:
            return self.opposite()

    def isAbsDipInRange(self, min_val, max_val, min_val_incl=False, max_value_incl=True):
        """
        Check whether the absolute value of the dip angle of a GVect instances is within a given range
        (default: minimum value is not included, maximum value is included).

        :param min_val: the minimum dip angle, positive, domain: 0-90°.
        :param max_val: the maximum dip angle, positive, domain: 0-90°.
        :param min_val_incl: is minimum value included, boolean.
        :param max_value_incl: is maximum value included, boolean.
        :return: Boolean

        Examples:
          >>> GVect(90, -45).isAbsDipInRange(30, 60)
          True
          >>> GVect(120, 0).isAbsDipInRange(0, 60)
          False
          >>> GVect(120, 0).isAbsDipInRange(0, 60, min_val_incl=True)
          True
          >>> GVect(120, 60).isAbsDipInRange(0, 60)
          True
          >>> GVect(120, 60).isAbsDipInRange(0, 60, max_value_incl=False)
          False
        """

        abs_dip = abs(self.pl)

        if abs_dip < min_val or abs_dip > max_val:
            return False
        elif abs_dip == min_val:
            if min_val_incl:
                return True
            else:
                return False
        elif abs_dip == max_val:
            if max_value_incl:
                return True
            else:
                return False
        else:
            return True

    def isSubHorizontal(self, max_dip_angle=DIP_ANGLE_THRESHOLD):
        """
        Check whether the instance is almost horizontal.

        Examples:
          >>> GVect(10, 15).isSubHorizontal()
          False
          >>> GVect(257, 2).isSubHorizontal()
          True
          >>> GVect(90, -5).isSubHorizontal()
          False
        """

        return abs(self.pl) < max_dip_angle

    def isSubVertical(self, min_dip_angle=90.0 - DIP_ANGLE_THRESHOLD):
        """
        Check whether the instance is almost vertical.

        Examples:
          >>> GVect(10, 15).isSubVertical()
          False
          >>> GVect(257, 89).isSubVertical()
          True
        """

        return abs(self.pl) > min_dip_angle

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two GVect instances or a GVect and a GAXis instances.
        Range is 0°-180°.

        Examples:
          >>> are_close(GVect(0, 90).angle(GVect(90, 0)), 90)
          True
          >>> are_close(GVect(0, 0).angle(GVect(270, 0)), 90)
          True
          >>> are_close(GVect(0, 0).angle(GVect(0, 0)), 0)
          True
          >>> are_close(GVect(0, 0).angle(GVect(180, 0)), 180)
          True
          >>> are_close(GVect(90, 0).angle(GVect(270, 0)), 180)
          True
          >>> are_close(GVect(0, 90, True).angle(GVect(90, 0, True)), 90)
          True
          >>> are_close(GVect(0, 0, True).angle(GVect(270, 0, True)), 90)
          True
          >>> are_close(GVect(0, 0, True).angle(GVect(0, 0, True)), 0)
          True
          >>> are_close(GVect(0, 0, True).angle(GVect(180, 0, True)), 0)
          True
          >>> are_close(GVect(0, 0, True).angle(GVect(179, 0, True)), 1)
          True
          >>> are_close(GVect(0, -90, True).angle(GVect(0, 90, True)), 0)
          True
          >>> are_close(GVect(90, 0, True).angle(GVect(315, 0, True)), 45)
          True
        """

        angle_vers = self.versor().angle(another.versor())

        if self.ax or another.ax:
            return min(angle_vers, 180. - angle_vers)
        else:
            return angle_vers

    def isAlmostParallel(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two GVect instances, or a GVect and a GPlane instances, are sub-parallel,

        :param another: a GVect or a GPlane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> GVect(0, 90).isAlmostParallel(GVect(90, 0))
          False
          >>> GVect(0, 0).isAlmostParallel(GVect(0, 1e-6))
          True
          >>> GVect(0, 90).isAlmostParallel(GVect(180, 0))
          False
          >>> GVect(0, 90).isAlmostParallel(GVect(0, -90))
          False
          >>> GVect(0, 90, True).isAlmostParallel(GVect(90, 0, True))
          False
          >>> GVect(0, 0, True).isAlmostParallel(GVect(0, 1.e-6, True))
          True
          >>> GVect(0, 0, True).isAlmostParallel(GVect(180, 0, True))
          True
          >>> GVect(90, 0, True).isAlmostParallel(GVect(270, 0, True))
          True
          >>> GVect(0, 90, True).isAlmostParallel(GVect(0, -90, True))
          True
          >>> GVect(90, 90).isAlmostParallel(GPlane(315, 90))
          True
          >>> GVect(90, 45).isAlmostParallel(GPlane(312, 27))
          False
          >>> GVect(0, 10).isAlmostParallel(GPlane(359.5, 10))
          True
          >>> GVect(0, 10).isAlmostParallel(GPlane(0, 12))
          False
        """

        fst_gvect = self

        if isinstance(another, GVect):
            snd_geoelem = another
        elif isinstance(another, GPlane):
            snd_geoelem = another.normGAxis()
        else:
            raise GeomInputException("Argument must be GPlane or GVect")

        angle = fst_gvect.angle(snd_geoelem)

        if isinstance(another, GPlane):
            return angle > (90.0 - angle_tolerance)
        else:
            return angle <= angle_tolerance

    def isAlmostAntiParallel(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two GVect instances are almost anti-parallel,

        :param another: a GVect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> GVect(0, 90).isAlmostAntiParallel(GVect(90, -89.5))
          True
          >>> GVect(0, 0).isAlmostAntiParallel(GVect(180, 1e-6))
          True
          >>> GVect(90, 45).isAlmostAntiParallel(GVect(270, -45.5))
          True
          >>> GVect(45, 90).isAlmostAntiParallel(GVect(0, -90))
          True
          >>> GVect(45, 72).isAlmostAntiParallel(GVect(140, -38))
          False
        """

        if self.ax or another.ax:
            raise GeomInputException("Both elements must be GVect instances")

        return self.angle(another) > (180.0 - angle_tolerance)

    def isSubOrthogonal(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two GVect instance are sub-orthogonal

        :param another: a GVect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees) from orthogonality
        :return: Boolean

         Examples:
          >>> GVect(0, 90).isSubOrthogonal(GVect(90, 0))
          True
          >>> GVect(0, 0).isSubOrthogonal(GVect(0, 1.e-6))
          False
          >>> GVect(0, 0).isSubOrthogonal(GVect(180, 0))
          False
          >>> GVect(90, 0).isSubOrthogonal(GVect(270, 89.5))
          True
          >>> GVect(0, 90).isSubOrthogonal(GVect(0, 0.5))
          True
        """

        return 90.0 - angle_tolerance <= self.angle(another) <= 90.0 + angle_tolerance

    def normVersor(self, another):
        """
        Calculate the versor (GVect) defined by the vector product of two GVect instances.

        Examples:
          >>> GVect(0, 0).normVersor(GVect(90, 0))
          Vect(0.0000, 0.0000, -1.0000)
          >>> GVect(45, 0).normVersor(GVect(310, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> GVect(0, 0).normVersor(GVect(90, 90))
          Vect(-1.0000, 0.0000, -0.0000)
          >>> GVect(315, 45).normVersor(GVect(315, 44.5)) is None
          True
        """

        if not isinstance(self, GAxis) and not isinstance(self, GAxis) and self.isAlmostAntiParallel(another):
            return None
        elif self.isAlmostParallel(another):
            return None
        else:
            return self.versor().vCross(another.versor()).versor()

    def normGPlane(self):
        """
        Return the geological plane that is normal to the geological vector.

        Examples:
          >>> GVect(0, 45).normGPlane()
          GPlane(180.00, +45.00)
          >>> GVect(0, -45).normGPlane()
          GPlane(000.00, +45.00)
          >>> GVect(0, 90).normGPlane()
          GPlane(180.00, +00.00)
        """

        down_axis = self.downward()
        dipdir = (down_axis.tr + 180.0) % 360.0
        dipangle = 90.0 - down_axis.pl

        return GPlane(dipdir, dipangle)

    def commonGPlane(self, another):
        """
        Calculate GPlane instance defined by the two GVect instances.

        Examples:
          >>> GVect(0, 0).commonGPlane(GVect(90, 0)).isAlmostParallel(GPlane(180.0, 0.0))
          True
          >>> GVect(0, 0).commonGPlane(GVect(90, 90)).isAlmostParallel(GPlane(90.0, 90.0))
          True
          >>> GVect(45, 0).commonGPlane(GVect(135, 45)).isAlmostParallel(GPlane(135.0, 45.0))
          True
          >>> GVect(315, 45).commonGPlane(GVect(135, 45)).isAlmostParallel(GPlane(225.0, 90.0))
          True
          >>> GVect(0, 0).commonGPlane(GVect(90, 0)).isAlmostParallel(GPlane(180.0, 10.0))
          False
          >>> GVect(0, 0).commonGPlane(GVect(90, 90)).isAlmostParallel(GPlane(90.0, 80.0))
          False
          >>> GVect(45, 0).commonGPlane(GVect(135, 45)).isAlmostParallel(GPlane(125.0, 45.0))
          False
          >>> GVect(315, 45).commonGPlane(GVect(135, 45)).isAlmostParallel(GPlane(225.0, 80.0))
          False
          >>> GVect(315, 45).commonGPlane(GVect(315, 44.5)) is None
          True
        """

        normal_versor = self.normVersor(another)
        if normal_versor is None:
            return None
        else:
            return normal_versor.asGVect().normGPlane()

    def asGAxis(self):
        """
        Create GAxis instance with the same attitude as the self instance.

        Example:
          >>> GVect(220, 32).asGAxis()
          GAxis(220.00, +32.00)
        """

        return GAxis(*self.tp)

    def normGVect(self, another):
        """
        Calculate the GVect instance that is normal to the two provided sources.
        Angle between sources must be larger than MIN_ANGLE_DEGR_DISORIENTATION,
        otherwise a SubparallelLineationException will be raised.
        
        Example:
          >>> GVect(0, 0).normGVect(GVect(0.5, 0)) is None
          True
          >>> GVect(0, 0).normGVect(GVect(179.5, 0)) is None
          True
          >>> GVect(0, 0).normGVect(GVect(5.1, 0))
          GVect(000.00, +90.00)
          >>> GVect(90, 45).normGVect(GVect(90, 0))
          GVect(180.00, +00.00)
        """

        if not isinstance(self, GAxis) and not isinstance(self, GAxis) and self.isAlmostAntiParallel(another):
            return None
        elif self.isAlmostParallel(another):
            return None
        else:
            return self.normVersor(another).asGVect()


class GAxis(GVect):
    """
    Convenience class wrapping geological axis with unknown movement sense.
    While GAxis is non-directional, the geological vector (GVect) is directional.
    Defined by trend and plunge (both in degrees):
     - trend: [0.0, 360.0[ clockwise, from 0 (North):
     - plunge: [-90.0, 90.0].
    """

    def __init__(self, trend, plunge):

        super().__init__(trend, plunge, is_axis=True)

    def __repr__(self):

        return "GAxis({:06.2f}, {:+06.2f})".format(*self.tp)

    def asGVect(self):
        """
        Create GVect instance with the same attitude as the self instance.
        
        Example:
          >>> GAxis(220, 32).asGVect()
          GVect(220.00, +32.00)
        """

        return GVect(self.tr, self.pl, is_axis=False)

    def normGAxis(self, another):
        """
        Calculate the GAxis instance that is perpendicular to the two provided.
        The two source GAxis must not be subparallel (threshold is MIN_ANGLE_DEGR_DISORIENTATION),
        otherwise a SubparallelLineationException will be raised.
        
        Example:
          >>> GAxis(0, 0).normGAxis(GAxis(0.5, 0)) is None
          True
          >>> GAxis(0, 0).normGAxis(GAxis(180, 0)) is None
          True
          >>> GAxis(90, 0).normGAxis(GAxis(180, 0))
          GAxis(000.00, +90.00)
          >>> GAxis(90, 45).normGAxis(GAxis(180, 0))
          GAxis(270.00, +45.00)
          >>> GAxis(270, 45).normGAxis(GAxis(180, 90)).isAlmostParallel(GAxis(180, 0))
          True
        """

        norm_gvect = self.normGVect(another)
        if norm_gvect is None:
            return None
        else:
            return norm_gvect.asGAxis()


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with GPlane, defined just by its attitude, but with undefined position

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
          Plane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          Plane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = np.array([[pt1.y, pt1.z, 1],
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1]])

        matr_b = - np.array([[pt1.x, pt1.z, 1],
                             [pt2.x, pt2.z, 1],
                             [pt3.x, pt3.z, 1]])

        matr_c = np.array([[pt1.x, pt1.y, 1],
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1]])

        matr_d = - np.array([[pt1.x, pt1.y, pt1.z],
                             [pt2.x, pt2.y, pt2.z],
                             [pt3.x, pt3.y, pt3.z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def __repr__(self):

        return "Plane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

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

    def toGPlanePoint(self):
        """
        Converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution).

        Examples:
          >>> gpl, pt = CPlane(0, 0, 1, -1).toGPlanePoint()
          >>> gpl
          GPlane(000.00, +00.00)
          >>> pt
          Point(0.0000, 0.0000, 1.0000)
        """

        geol_plane = self.normVersor().asGVect().normGPlane()
        point = Point(*point_solution(np.array([[self.a, self.b, self.c]]),
                                      np.array([-self.d])))
        return geol_plane, point

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
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
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


class GPlane(object):
    """
    Geological plane.
    Defined by dip direction and dip angle (both in degrees):
     - dip direction: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
    """

    def __init__(self, azim: float, dip_ang: float, is_rhr_strike=False):
        """
        Geological plane constructor.

        :param  azim:  azimuth of the plane (RHR strike or dip direction).
        :type  azim:  number or string convertible to float.
        :param  dip_ang:  Dip angle of the plane (0-90°).
        :type  dip_ang:  number or string convertible to float.
        :param is_rhr_strike: if the source azimuth is RHR strike (default is False, i.e. it is dip direction)
        :return: the instantiated geological plane.
        :rtype: GPlane.

        Example:
          >>> GPlane(0, 90)
          GPlane(000.00, +90.00)
          >>> GPlane(0, 90, is_rhr_strike=True)
          GPlane(090.00, +90.00)
          >>> GPlane(0, 90, True)
          GPlane(090.00, +90.00)
          >>> GPlane(0, "90", True)
          Traceback (most recent call last):
          ...
          GeomInputException: Source dip angle must be number
          >>> GPlane(0, 900)
          Traceback (most recent call last):
          ...
          GeomInputException: Dip angle must be between 0° and 90°
        """

        def rhrstrk2dd(rhr_strk):
            """Converts RHR strike value to dip direction value.

            Example:
                >>> rhrstrk2dd(285.5)
                15.5
            """

            return (rhr_strk + 90.0) % 360.0

        if not isinstance(azim, (int, float)):
            raise GeomInputException("Source azimuth must be number")
        if not isinstance(dip_ang, (int, float)):
            raise GeomInputException("Source dip angle must be number")
        if not isinstance(is_rhr_strike, bool):
            raise GeomInputException("Source azimuth type must be boolean")

        if not (0.0 <= dip_ang <= 90.0):
            raise GeomInputException("Dip angle must be between 0° and 90°")

        if is_rhr_strike:
            self._dipdir = rhrstrk2dd(azim)
        else:
            self._dipdir = azim % 360.0
        self._dipangle = float(dip_ang)

    @property
    def dd(self):
        """
        Return the dip direction of the geological plane.

        Example:
          >>> GPlane(34.2, 89.7).dd
          34.2
        """

        return self._dipdir

    @property
    def da(self):
        """
        Return the dip angle of the geological plane.

        Example:
          >>> GPlane(183, 77).da
          77.0

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Return a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = GPlane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)
        """
        
        return self.dd, self.da

    @property
    def rhrStrike(self):
        """
        Return the strike according to the right-hand-rule.

        Examples:
          >>> GPlane(90, 45).rhrStrike
          0.0
          >>> GPlane(45, 89).rhrStrike
          315.0
          >>> GPlane(275, 38).rhrStrike
          185.0
          >>> GPlane(0, 38).rhrStrike
          270.0
        """

        return (self.dd - 90.0) % 360.0

    @property
    def srda(self):
        """
        Return a tuple storing the right-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> GPlane(100, 17.2).srda
          (10.0, 17.2)
          >>> GPlane(10, 87).srda
          (280.0, 87.0)
        """

        return self.rhrStrike, self.da

    @property
    def lhrStrike(self):
        """
        Return the strike according to the left-hand-rule.

        Examples:
          >>> GPlane(90, 45).lhrStrike
          180.0
          >>> GPlane(45, 89).lhrStrike
          135.0
          >>> GPlane(275, 38).lhrStrike
          5.0
          >>> GPlane(0, 38).lhrStrike
          90.0
        """

        return (self.dd + 90.0) % 360.0

    @property
    def slda(self):
        """
        Return a tuple storing the left-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> GPlane(100, 17.2).slda
          (190.0, 17.2)
          >>> GPlane(10, 87).slda
          (100.0, 87.0)
        """

        return self.lhrStrike, self.da

    def __repr__(self):

        return "GPlane({:06.2f}, {:+06.2f})".format(*self.dda)

    def rhrStrikeGVect(self):
        """
        Creates a GVect instance that is parallel to the right-hand rule strike.

        :return: GVect instance,

        Examples:
          >>> GPlane(90, 45).rhrStrikeGVect()
          GVect(000.00, +00.00)
          >>> GPlane(45, 17).rhrStrikeGVect()
          GVect(315.00, +00.00)
          >>> GPlane(90, 0).rhrStrikeGVect()
          GVect(000.00, +00.00)
        """

        return GVect(
            trend=self.rhrStrike,
            plunge=0.0)

    def lhrStrikeGVect(self):
        """
        Creates a GVect instance that is parallel to the left-hand rule strike.

        :return: GVect instance.

        Examples:
          >>> GPlane(90, 45).lhrStrikeGVect()
          GVect(180.00, +00.00)
          >>> GPlane(45, 17).lhrStrikeGVect()
          GVect(135.00, +00.00)
        """

        return GVect(
            trend=self.lhrStrike,
            plunge=0.0)

    def dipDirGVect(self):
        """
        Creates a GVect instance that is parallel to the dip direction.

        :return: GVect instance.

        Examples:
          >>> GPlane(90, 45).dipDirGVect()
          GVect(090.00, +45.00)
          >>> GPlane(45, 17).dipDirGVect()
          GVect(045.00, +17.00)
        """

        return GVect(
            trend=self.dd,
            plunge=self.da)

    def dipDirOppGVect(self):
        """
        Creates a GVect instance that is anti-parallel to the dip direction.

        :return: GVect instance.

        Examples:
          >>> GPlane(90, 45).dipDirOppGVect()
          GVect(270.00, -45.00)
          >>> GPlane(45, 17).dipDirOppGVect()
          GVect(225.00, -17.00)
        """

        return self.dipDirGVect().opposite()

    def mirrorVertGPlane(self):
        """
        Mirror a geological plane around a vertical plane
        creating a new one that has a dip direction opposite
        to the original one but with downward plunge.

        :return: geological plane
        :rtype: GPlane

        Examples:
          >>> GPlane(0, 45).mirrorVertGPlane()
          GPlane(180.00, +45.00)
          >>> GPlane(225, 80).mirrorVertGPlane()
          GPlane(045.00, +80.00)
          >>> GPlane(90, 90).mirrorVertGPlane()
          GPlane(270.00, +90.00)
          >>> GPlane(270, 0).mirrorVertGPlane()
          GPlane(090.00, +00.00)
        """

        return GPlane(
            azim=opposite_trend(self.dd),
            dip_ang=self.da)

    def _normal_gv_frwrd(self):
        """
        Return the geological vector normGVect to the geological plane,
        pointing in the same direction as the geological plane.

        Example:
            >>> GPlane(90, 55)._normal_gv_frwrd()
            GVect(090.00, -35.00)
            >>> GPlane(90, 90)._normal_gv_frwrd()
            GVect(090.00, +00.00)
            >>> GPlane(90, 0)._normal_gv_frwrd()
            GVect(090.00, -90.00)
        """
        
        trend = self.dd % 360.0
        plunge = self.da - 90.0

        return GVect(
            trend=trend,
            plunge=plunge)

    def _normal_gv_anti(self):
        """
        Return the geological vector normal to the geological plane,
        pointing in the opposite direction to the geological plane.

        Example:
            >>> GPlane(90, 55)._normal_gv_anti()
            GVect(270.00, +35.00)
            >>> GPlane(90, 90)._normal_gv_anti()
            GVect(270.00, -00.00)
            >>> GPlane(90, 0)._normal_gv_anti()
            GVect(270.00, +90.00)
        """

        return self._normal_gv_frwrd().opposite()

    def downNormGVect(self):
        """
        Return the geological vector normGVect to the geological plane,
        pointing downward.

        Example:
            >>> GPlane(90, 55).downNormGVect()
            GVect(270.00, +35.00)
            >>> GPlane(90, 90).downNormGVect()
            GVect(090.00, +00.00)
            >>> GPlane(90, 0).downNormGVect()
            GVect(270.00, +90.00)
        """

        return self._normal_gv_frwrd().downward()

    def upNormGVect(self):
        """
        Return the geological vector normGVect to the geological plane,
        pointing upward.

        Example:
            >>> GPlane(90, 55).upNormGVect()
            GVect(090.00, -35.00)
            >>> GPlane(90, 90).upNormGVect()
            GVect(090.00, +00.00)
            >>> GPlane(90, 0).upNormGVect()
            GVect(090.00, -90.00)
        """

        return self._normal_gv_frwrd().upward()


    def normGVect(self):
        """
        Wrapper to down_normal_gv.

        :return: GVect normGVect to the GPlane self instance
        """

        return self.downNormGVect()

    def normGAxis(self):
        """
        Normal GAxis.

        :return: GAxis normal to the GPlane self instance
        """

        return self.downNormGVect().asGAxis()

    def plane(self, point):
        """
        Given a GPlane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> GPlane(0, 0).plane(Point(0, 0, 0))
          Plane(0.0000, 0.0000, 1.0000, -0.0000)
          >>> GPlane(90, 45).plane(Point(0, 0, 0))
          Plane(0.7071, 0.0000, 0.7071, -0.0000)
          >>> GPlane(0, 90).plane(Point(0, 0, 0))
          Plane(0.0000, 1.0000, -0.0000, -0.0000)
        """

        normal_versor = self._normal_gv_frwrd().versor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return CPlane(a, b, c, d)

    def angle(self, another):
        """
        Calculate angle (in degrees) between two geoplanes, or between a geoplane and a GAxis/GVect instance.
        Range is 0°-90°.

        >>> GPlane(100.0, 50.0).angle(GPlane(100.0, 50.0))
        0.0
        >>> GPlane(300.0, 10.0).angle(GPlane(300.0, 90.0))
        80.0
        >>> GPlane(90.0, 90.0).angle(GPlane(270.0, 90.0))
        0.0
        >>> are_close(GPlane(90.0, 90.0).angle(GPlane(130.0, 90.0)), 40)
        True
        >>> are_close(GPlane(90, 70).angle(GPlane(270, 70)), 40)
        True
        >>> are_close(GPlane(90.0, 10.0).angle(GPlane(270.0, 10.0)), 20.0)
        True
        >>> are_close(GPlane(90.0, 10.0).angle(GPlane(270.0, 30.0)), 40.0)
        True
        >>> are_close(GPlane(0, 0).angle(GVect(90, 70)), 70.0)
        True
        >>> GPlane(0, 0).angle(GAxis(90, -60))
        60.0
        """

        gpl_axis = self._normal_gv_frwrd().asGAxis()
        if isinstance(another, GPlane):
            an_axis = another._normal_gv_frwrd().asGAxis()
        elif isinstance(another, GAxis):
            an_axis = another
        elif isinstance(another, GVect):
            an_axis = another.asGAxis()
        else:
            raise GeomInputException("Provided another instance for angle is of {} type".format(type(another)))

        angle = gpl_axis.angle(an_axis)

        if isinstance(another, GPlane):
            return angle
        else:
            return 90.0 - angle

    def isAlmostParallel(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two GPlanes are sub-parallel

        :param another: a GPlane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> GPlane(0, 90).isAlmostParallel(GPlane(270, 90))
          False
          >>> GPlane(0, 90).isAlmostParallel(GPlane(180, 90))
          True
          >>> GPlane(0, 90).isAlmostParallel(GPlane(0, 0))
          False
          >>> GPlane(0, 0).isAlmostParallel(GPlane(0, 1e-6))
          True
          >>> GPlane(0, 0).isAlmostParallel(GPlane(0, 1.1))
          False
        """

        return self.angle(another) < angle_tolerance

    def isSubOrthogonal(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two GPlanes, or the Gplane instances and a GVect/GAxis, are sub-orthogonal.

        :param another: a GPlane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> GPlane(0, 90).isSubOrthogonal(GPlane(270, 90))
          True
          >>> GPlane(0, 90).isSubOrthogonal(GPlane(180, 90))
          False
          >>> GPlane(0, 90).isSubOrthogonal(GPlane(0, 0))
          True
          >>> GPlane(0, 0).isSubOrthogonal(GPlane(0, 88))
          False
          >>> GPlane(0, 0).isSubOrthogonal(GPlane(0, 45))
          False
          >>> GPlane(0, 0).isSubOrthogonal(GVect(0, 90))
          True
          >>> GPlane(0, 0).isSubOrthogonal(GVect(0, 45))
          False
          >>> GPlane(0, 0).isSubOrthogonal(GVect(0, 0))
          False
          >>> GPlane(0, 0).isSubOrthogonal(GAxis(0, -90))
          True
          >>> GPlane(0, 0).isSubOrthogonal(GAxis(0, 45))
          False
          >>> GPlane(90, 90).isSubOrthogonal(GAxis(90, 0.5))
          True
        """

        fst_gaxis = self.normGVect().asGAxis()

        if isinstance(another, GPlane):
            snd_gaxis = another.normGVect().asGAxis()
        elif isinstance(another, GAxis):
            snd_gaxis = another
        elif isinstance(another, GVect):
            snd_gaxis = another.asGAxis()
        else:
            raise GeomInputException("Not accepted argument type for isSubOrthogonal method")

        angle = fst_gaxis.angle(snd_gaxis)

        if isinstance(another, GPlane):
            return angle > 90.0 - angle_tolerance
        else:
            return angle < angle_tolerance

    def rakeToGVect(self, rake):
        """
        Calculate GVect given a GPlane instance and a rake value.
        The rake is defined according to the Aki and Richards, 1980 conventions:
        rake = 0° -> left-lateral
        rake = 90° -> reverse
        rake = +/- 180° -> right-lateral
        rake = -90° -> normal

        Examples:
          >>> GPlane(180, 45).rakeToGVect(0.0)
          GVect(090.00, +00.00)
          >>> GPlane(180, 45).rakeToGVect(90.0)
          GVect(000.00, -45.00)
          >>> GPlane(180, 45).rakeToGVect(-90.0)
          GVect(180.00, +45.00)
          >>> GPlane(180, 45).rakeToGVect(180.0).isAlmostParallel(GVect(270.00, -00.00))
          True
          >>> GPlane(180, 45).rakeToGVect(-180.0)
          GVect(270.00, +00.00)
        """

        rk = radians(rake)
        strk = radians(self.rhrStrike)
        dip = radians(self.da)

        x = cos(rk)*sin(strk)-sin(rk)*cos(dip)*cos(strk)
        y = cos(rk)*cos(strk)+sin(rk)*cos(dip)*sin(strk)
        z = sin(rk) * sin(dip)

        return Vect(x, y, z).asGVect()

    def isVLowAngle(self, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a geological plane is very low angle.

        :param threshold: the limit for the plane angle, in degrees
        :type threshold: float
        :return: bool flag indicating if it is very low angle

        Examples:
          >>> GPlane(38.9, 1.2).isVLowAngle()
          True
          >>> GPlane(38.9, 7.4).isVLowAngle()
          False
        """

        return self.da < dip_angle_threshold

    def isVHighAngle(self, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a geological plane is very high angle.

        :param threshold: the limit for the plane angle, in degrees
        :type threshold: float
        :return: bool flag indicating if it is very high angle

        Examples:
          >>> GPlane(38.9, 11.2).isVHighAngle()
          False
          >>> GPlane(38.9, 88.4).isVHighAngle()
          True
        """

        return self.da > (90.0 - dip_angle_threshold)


class GeomInputException(Exception):
    """
    Exception for geometric input.
    """

    pass


class VectInvalidException(Exception):
    """
    Exception for invalid Vect.
    """

    pass


class SubparallelLineationException(Exception):
    """
    Exception for subparallel GAxis/GVect instances.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
