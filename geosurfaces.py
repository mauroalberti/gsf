# -*- coding: utf-8 -*-

from __future__ import division

from math import sqrt, sin, cos, tan, radians, asin, acos, atan, atan2, degrees
import numpy as np

from array_utils import point_solution
from geosurf_utils import rhrstrk2dd


MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10


class CPoint(object):
    """
    Cartesian point.
    Dimensions: 3D + time
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan, t=None):

        self._x = x
        self._y = y
        self._z = z
        self._t = t

    @property
    def x(self):
        """
        Return x value
        """

        return self._x

    @property
    def y(self):
        """
        Return y value
        """
        return self._y

    @property
    def z(self):
        """
        Return z value
        """
        return self._z

    @property
    def t(self):
        """
        Return time value
        """
        return self._t

    def clone(self):
        """
        Clone the point
        """

        return CPoint(self.x, self.y, self.z, self.t)

    def dist_3d(self, another):
        """
        Calculate Euclidean spatial distance between two points.
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist_2d(self, another):
        """
        Calculate horizontal (2D) distance between two points.
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def coincident(self, another, tolerance=MINIMUM_SEPARATION_THRESHOLD):
        """
        Check spatial coincidence of two points
        """

        if self.dist_3d(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0):
        """
        Create a new point shifted by given amount from the self instance.
       """

        return CPoint(self.x + sx, self.y + sy, self.z + sz, self.t)

    def translate_with_vector(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector
        """

        return CPoint(self.x + displ_vect.x, self.y + displ_vect.y,
                      self.z + displ_vect.z, self.t)

    @property
    def vector(self):
        """
        Create a vector based on the point coordinates
        """

        return CVect(self.x, self.y, self.z)

    def delta_time(self, another):
        """
        Calculate the time difference between two points
        """

        return another.t - self.t

    def speed(self, another):
        """
        Calculate the speed required to displace self to another
        """
        try:
            return self.dist_3d(another) / self.delta_time(another)
        except:
            return np.Infinity


class CVect(object):
    """
    Cartesian vector, 3D
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan):
        """
        Constructor
        """
        self._x = x
        self._y = y
        self._z = z

    def __repr__(self):

        return "CVect({:.3f}, {:.3f}, {:.3f})".format(self._x, self._y, self._z)

    @property
    def x(self):
        """
        Return x value of vector
        """

        return self._x

    @property
    def y(self):
        """
        Return y value of vector
        """

        return self._y

    @property
    def z(self):
        """
        Return z value of vector
        """

        return self._z

    def __add__(self, another):
        """
        Sum of two vectors
        """

        return CVect(self.x + another.x,
                     self.y + another.y,
                     self.z + another.z)

    def clone(self):
        """
        Clone the vector
        """
        return CVect(self.x, self.y, self.z)

    def __abs__(self):
        """
        Vector magnitude
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def length_horiz(self):
        """
        Vector length on the horizontal (xy) plane
        """
        return sqrt(self.x * self.x + self.y * self.y)

    def scale(self, scale_factor):
        """
        Create a scaled vector
        """

        return CVect(self.x * scale_factor,
                     self.y * scale_factor,
                     self.z * scale_factor)

    @property
    def versor(self):
        """
        Calculate versor

        Example:
          >>> print CVect(5, 0, 0).versor
          CVect(1.000, 0.000, 0.000)
        """

        return self.scale(1.0 / abs(self))

    @property
    def downvector(self):
        """
        Calculate new vector pointing downwards
        """

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

    def slope_radians(self):

        return atan(self.z / self.length_horiz)

    def as_geolaxis(self):

        if abs(self) < MINIMUM_VECTOR_MAGNITUDE:
            return None

        unit_vect = self.versor

        plunge = - degrees(asin(unit_vect.z))  # upward negative, downward positive

        trend = 90.0 - degrees(atan2(unit_vect.y, unit_vect.x))
        if trend < 0.0:
            trend += 360.0
        elif trend > 360.0:
            trend -= 360.0

        assert 0.0 <= trend < 360.0
        assert -90.0 <= plunge <= 90.0

        return GDirect(trend, plunge)

    def scalar_product(self, another):

        return self.x * another.x + self.y * another.y + self.z * another.z

    def vectors_cos_angle(self, another):

        try:
            val = self.scalar_product(another) / (abs(self) * abs(another))
            if val > 1.0:
                return 1.0
            elif val < -1.0:
                return -1.0
            else:
                return val
        except ZeroDivisionError:
            return np.nan

    def angle(self, another):
        """
        angle between two vectors,
        in 0 - pi range
        """

        return degrees(acos(self.vectors_cos_angle(another)))

    def vector_product(self, another):

        x = self.y * another.z - self.z * another.y
        y = self.z * another.x - self.x * another.z
        z = self.x * another.y - self.y * another.x

        return CVect(x, y, z)

    def by_matrix(self, matrix3x3):

        vx = matrix3x3[0, 0] * self.x + matrix3x3[0, 1] * self.y + matrix3x3[0, 2] * self.z
        vy = matrix3x3[1, 0] * self.x + matrix3x3[1, 1] * self.y + matrix3x3[1, 2] * self.z
        vz = matrix3x3[2, 0] * self.x + matrix3x3[2, 1] * self.y + matrix3x3[2, 2] * self.z

        return CVect(vx, vy, vz)


class GDirect(object):
    """
    Geological direction.
    Defined by trend and plunge (both in degrees)
     - trend: [0.0, 360.0[ clockwise, from 0 (North)
     - plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """

    def __init__(self, srcTrend, srcPlunge):
        """

        :param srcTrend: Trend range: [0.0, 360.0[ clockwise, from 0 (North)
        :param srcPlunge: Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis

        Example:
          >>> a = GDirect(120, -27)
          >>> b = GDirect(54, -320)
          Traceback (most recent call last):
          ...
          AssertionError: plunge must be between -90° and +90° (comprised)
        """

        assert -90.0 <= srcPlunge <= 90.0, "plunge must be between -90° and +90° (comprised)"
        self._trend = srcTrend % 360.0
        self._plunge = float(srcPlunge)

    @property
    def trend(self):
        """
        Returns trend of the geological direction.
        Range is [0, 360[

        Example:
          >>> GDirect(420, -17).trend
          60.0
          >>> GDirect(-20, 49).trend
          340.0
        """

        return self._trend

    @property
    def plunge(self):
        """
        Returns plugne of the geological direction.
        Range is [-90, 90]

        Example:
          >>> GDirect(420, -17).plunge
          -17.0
        """

        return self._plunge

    @property
    def tp(self):
        """
        Returns trend and plunge of the geological direction

        Example:
          >>> GDirect(-90, -45).tp
          (270.0, -45.0)
        """
        
        return self.trend, self.plunge

    @property
    def versor(self):
        """
        Returns the CVect corresponding to the direction

        Examples:
          >>> print GDirect(0, 90).versor
          CVect(0.000, 0.000, -1.000)
          >>> print GDirect(0, -90).versor
          CVect(0.000, 0.000, 1.000)
        """

        north_coord = cos(radians(self.plunge)) * cos(radians(self.trend))
        east_coord = cos(radians(self.plunge)) * sin(radians(self.trend))
        down_coord = sin(radians(self.plunge))

        return CVect(east_coord, north_coord, -down_coord)

    def as_downgeolaxis(self):

        trend, plunge = self.trend, self.plunge
        if plunge < 0.0:
            trend = (trend + 180.0) % 360.0
            plunge = - plunge

        return GDirect(trend, plunge)

    def as_normalgeolplane(self):

        down_axis = self.as_downgeolaxis()
        dipdir = (down_axis.trend + 180.0) % 360.0
        dipangle = 90.0 - down_axis.plunge

        return GPlane(dipdir, dipangle)


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    """

    def __init__(self, a=None, b=None, c=None, d=None):
        self._a = a
        self._b = b
        self._c = c
        self._d = d

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def d(self):
        return self._d

    @classmethod
    def from_points(cls, pt1, pt2, pt3):
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

    @property
    def normal_versor3d(self):
        """
        return the normal versor to the cartesian plane
        """

        return CVect(self.a, self.b, self.c).versor

    def as_geolplane_and_point_3d(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.normal_versor3d.as_geolaxis().as_normalgeolplane()
        point = CPoint(point_solution(np.array([[self.a, self.b, self.c]]),
                                      np.array([-self.d])))
        return geol_plane, point

    def intersection_versor3d(self, another):
        """
        return intersection versor for two intersecting planes
        """

        return self.normal_versor3d.vector_product(another.normal_versor3d).versor

    def intersection_point3dt(self, another):
        """
        return point on intersection line (obviously non-unique solution)
        for two planes
        """

        # find a point lying on the intersection line (this is a non-unique solution)    
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return CPoint(x, y, z)

    def set_point_inside(self, pt):
        return self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d

    def angle_degr(self, another):
        angle_degr = self.normal_versor3d.angle(another.normal_versor3d)

        assert angle_degr > 0.0

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr



class GPlane(object):
    """
    Geological plane.
    Defined by dip direction and dip angle.
    """

    def __init__(self, srcAzimuth, srcDipAngle, isRHRStrike=False):
        """
        Class constructor

        @param  srcAzimuth:  Azimuth of the plane (RHR strike or dip direction).
        @type  srcAzimuth:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90°).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  GeolPlane.

        Example:
          >>> gp = GPlane(0, 90)
    
        """

        if isRHRStrike:
            self._dipdir = rhrstrk2dd(srcAzimuth)
        else:
            self._dipdir = srcAzimuth % 360.0
        self._dipangle = srcDipAngle

    @property
    def dd(self):
        """

        Returns the dip direction of the geological plane.

        Example:
          >>> GPlane(34.2, 89.7).dd
          34.2

        """

        return self._dipdir

    @property
    def da(self):
        """
        Returns the dip angle of the geological plane.

        Example:
          >>> GPlane(183, 77).da
          77

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Returns a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = GPlane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)

        """
        
        return self.dd, self.da

    @property
    def normal(self):
        """
        Returns the normal to the plane, as a geological (downward) axis.

        Example:
            >>> ga = GPlane(90, 55).normal
            >>> ga.trend
            270.0
            >>> ga.plunge
            35.0
        """
        
        trend = (self.dd + 180.0) % 360.0
        plunge = 90.0 - self.da

        return GDirect(trend, plunge)

    def plane_x_coeff(self):
        """
        Calculate the slope of a given plane along the x direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.    
        """
        return - sin(radians(self.dd)) * tan(radians(self.da))

    def plane_y_coeff(self):
        """
        Calculate the slope of a given plane along the y direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.     
        """
        return - cos(radians(self.dd)) * tan(radians(self.da))

    def plane_from_geo(self, or_Pt):
        """
        Closure that embodies the analytical formula for a given, non-vertical plane.
        This closure is used to calculate the z value from given horizontal coordinates (x, y).
    
        @param  or_Pt:  CartesianPoint3DT instance expressing a location point contained by the plane.
        @type  or_Pt:  CartesianPoint3DT.
        
        @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.
        """

        x0 = or_Pt.x
        y0 = or_Pt.y
        z0 = or_Pt.z

        # slope of the line parallel to the x axis and contained by the plane
        a = self.plane_x_coeff()

        # slope of the line parallel to the y axis and contained by the plane
        b = self.plane_y_coeff()

        return lambda x, y: a * (x - x0) + b * (y - y0) + z0

    def as_cartesplane(self, point):

        normal_versor = self.normal.as_downgeolaxis().versor
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return CPlane(a, b, c, d)

    def angle_degr(self, another):
        """
        calculate angle (in degrees) between two planes

        >>> p1 = GPlane(100.0, 50.0)
        >>> p1.angle_degr(p1)
        0.0

        >>> p2 = GPlane(300.0, 10.0)
        >>> p3 = GPlane(300.0, 90.0)
        >>> p2.angle_degr(p3)
        80.0
        """

        vec0 = self.normal.versor
        vec1 = another.normal.versor

        return vec0.angle(vec1)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
