# -*- coding: utf-8 -*-

from __future__ import division

from math import sqrt, sin, cos, tan, radians, asin, acos, atan, atan2, degrees
import numpy as np

from array_utils import point_solution
from geosurf_utils import rhrstrk2dd


MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10
MINIMUM_SCALAR_VALUE = 1e-15
MINIMUM_VECTOR_MAGN_DIFF = MINIMUM_SCALAR_VALUE


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D + time
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan, t=np.nan):

        self._p = np.array([x, y, z, t], dtype=np.float64)

    def __repr__(self):

        return "Point({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z, self.t)

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a point from a numpy 1x4 array.

        Example:
          >>> Point.from_array(np.array([1, 0, 1]))
          Point(1.0000, 0.0000, 1.0000, nan)
        """

        obj = cls()

        assert 3 <= a.size <= 4
        b = a.astype(np.float64)
        if b.size == 3:
            c = np.append(b, [np.nan])
        else:
            c = b
        obj._p = c
        return obj

    @property
    def v(self):
        """
        Return values as array

        Example:
          >>> Point(1, 0, 0).v
          array([  1.,   0.,   0.,  nan])
        """

        return self._p

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.v[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.v[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.v[2]

    @property
    def t(self):
        """
        Return time value

        Example:
          >>> Point(1.5, 3.2, 41., 22.).t
          22.0
        """
        return self.v[3]

    def clone(self):
        """
        Clone the point.

        Example:
          >>> Point(1, 1, 1).clone()
          Point(1.0000, 1.0000, 1.0000, nan)
        """

        return Point.from_array(self.v)

    def __sub__(self, another):
        """Return point difference

        Example:
          >>> Point(1., 1., 1.) - Point(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000, nan)
        """

        return Point.from_array(self.v - another.v)

    def __abs__(self):
        """
        Point distance from frame origin.

        Example:
          >>> abs(Point(3, 4, 0))
          5.0
          >>> abs(Point(0, 12, 5))
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def dist_3d(self, another):
        """
        Calculate Euclidean spatial distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist_3d(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1, 4).dist_3d(Point(4, 5, 1, 14))
          5.0
        """

        return abs(self - another)

    def dist_2d(self, another):
        """
        Calculate horizontal (2D) distance between two points.

        Examples:
          >>> Point(1.,1.,1.).dist_2d(Point(4.,5.,7.))
          5.0
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def coincident(self, another, tolerance=MINIMUM_SEPARATION_THRESHOLD):
        """
        Check spatial coincidence of two points

        Example:
          >>> Point(1., 0., -1.).coincident(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).coincident(Point(1., 0., 0.))
          True
        """

        if self.dist_3d(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0, st=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).translate(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000, nan)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz, self.t + st)

    def vect_offset(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector.

        Example:
          >>> Point(1, 2, 0).vect_offset(Vect(10, 5, 0))
          Point(11.0000, 7.0000, 0.0000, nan)
        """

        return Point(self.x + displ_vect.x,
                     self.y + displ_vect.y,
                     self.z + displ_vect.z,
                     self.t)

    @property
    def vector(self):
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0, 5).vector
          Vect(1.0000, 1.0000, 0.0000)
        """

        return Vect(self.x, self.y, self.z)

    def delta_time(self, another):
        """
        Calculate the time difference between two points

        Example:
          >>> Point(1,1,1,4).delta_time(Point(1,1,2,5))
          1.0
        """

        return another.t - self.t

    def speed(self, another):
        """
        Calculate the speed required to displace self to another.

        Example:
          >>> Point(1, 1, 1, 4).speed(Point(4, 5, 1, 14))
          0.5
        """

        try:
            return self.dist_3d(another) / self.delta_time(another)
        except:
            return np.Infinity


class Vect(object):
    """
    Cartesian vector, 3D
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    z axis -> Up
    """

    def __init__(self, x=np.nan, y=np.nan, z=np.nan):
        """
        Vect constructor
        """

        self._v = np.array([x, y, z], dtype=np.float64)

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a vector from a numpy 1x3 array.

        Example:
          >>> Vect.from_array(np.array([1,0,1]))
          Vect(1.0000, 0.0000, 1.0000)
        """

        obj = cls()

        assert a.size == 3
        b = a.astype(np.float64)
        obj._v = b
        return obj

    @property
    def v(self):
        """
        Return the vector values as array

        Example:
          >>> Vect(1, 1, 0).v
          array([ 1.,  1.,  0.])
        """

        return self._v

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Vect(1, 2, 0).x
          1.0
        """

        return self.v[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Vect(1, 2, 0).y
          2.0
        """

        return self.v[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Vect(1, 2, 0).z
          0.0
        """

        return self.v[2]

    def __sub__(self, another):
        """Return vector difference.

        Example:
          >>> Vect(1., 1., 1.) - Vect(1., 1., 1.)
          Vect(0.0000, 0.0000, 0.0000)
        """

        return Vect.from_array(self.v - another.v)

    def __eq__(self, another):
        """Returns True if vectors are equal.

        Example:
          >>> Vect(1., 1., 1.) == Vect(1., 1., 1.)
          True
        """
        return bool(abs(self - another) < MINIMUM_VECTOR_MAGN_DIFF)

    def __ne__(self, another):
        """Returns False if vectors are equal.

        Example:
          >>> Vect(1., 1., 1.) != Vect(0., 0., 0.)
          True
        """

        return not self == another

    def __repr__(self):

        return "Vect({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __add__(self, another):
        """
        Sum of two vectors

        Example:
          >>> Vect(1, 0, 0) + Vect(0, 1, 1)
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(1, 1, 1) + Vect(-1, -1, -1)
          Vect(0.0000, 0.0000, 0.0000)
        """

        return Vect.from_array(self.v + another.v)

    def clone(self):
        """
        Clone the vector.

        Example:
          >>> Vect(1, 1, 1).clone()
          Vect(1.0000, 1.0000, 1.0000)
        """
        return Vect.from_array(self.v)

    def __abs__(self):
        """
        Vector magnitude.

        Example:
          >>> abs(Vect(3, 4, 0))
          5.0
          >>> abs(Vect(0, 12, 5))
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def hlen(self):
        """
        Vector length projected on the horizontal (xy) plane.

        Example:
          >>> Vect(3, 4, 7).hlen
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    def scale(self, scale_factor):
        """
        Create a scaled vector.

        Example;
          >>> Vect(1,0,1).scale(2.5)
          Vect(2.5000, 0.0000, 2.5000)
        """

        return Vect.from_array(self.v * scale_factor)

    @property
    def versor(self):
        """
        Calculate versor.

        Example:
          >>> Vect(5, 0, 0).versor
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor
          Vect(0.0000, 0.0000, -1.0000)
        """

        return self.scale(1.0 / abs(self))

    @property
    def downvect(self):
        """
        Calculate a new vector downward-pointing

        Example:
          >>> Vect(1, 1, 1).downvect
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(-1, -1, -1).downvect
          Vect(-1.0000, -1.0000, -1.0000)
        """

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

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
        """

        hlen = self.hlen
        if hlen == 0.0:
            if self.z > 0.:
                return -90.
            elif self.z < 0.:
                return 90.
            else:
                raise Exception("Zero-valued vector")
        else:
            slope = - degrees(atan(self.z / self.hlen))
            if abs(slope) > MINIMUM_SCALAR_VALUE:
                return slope
            else:
                return 0.

    @property
    def gvect(self):
        """
        Calculate the geological axis parallel to the original vector.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(0, 1, 1).gvect
          GVect(0.00, -45.00)
          >>> Vect(1, 0, 1).gvect
          GVect(90.00, -45.00)
          >>> Vect(0, 0, 1).gvect
          GVect(0.00, -90.00)
          >>> Vect(0, 0, -1).gvect
          GVect(0.00, 90.00)
          >>> Vect(-1, 0, 0).gvect
          GVect(270.00, 0.00)
          >>> Vect(0, -1, 0).gvect
          GVect(180.00, 0.00)
          >>> Vect(-1, -1, 0).gvect
          GVect(225.00, 0.00)
        """

        if abs(self) < MINIMUM_VECTOR_MAGNITUDE:
            raise Exception("Provided vector has near-zero magnitude")

        plunge = self.slope  # upward pointing -> negative value, downward -> positive

        unit_vect = self.versor
        if unit_vect.y == 0. and unit_vect.x == 0:
            trend = 0.
        else:
            trend = (90. - degrees(atan2(unit_vect.y, unit_vect.x))) % 360.

        return GVect(trend, plunge)

    def sp(self, another):
        """
        Vector scalar product.

        Examples:
          >>> Vect(1, 0, 0).sp(Vect(1, 0, 0))
          1.0
          >>> Vect(1, 0, 0).sp(Vect(0, 1, 0))
          0.0
          >>> Vect(1, 0, 0).sp(Vect(-1, 0, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y + self.z * another.z

    def cos_angle(self, another):
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect(1,0,0).cos_angle(Vect(0,0,1))
          0.0
          >>> Vect(1,0,0).cos_angle(Vect(-1,0,0))
          -1.0
          >>> Vect(1,0,0).cos_angle(Vect(1,0,0))
          1.0
        """

        try:
            val = self.sp(another) / (abs(self) * abs(another))
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
        angle between two vectors, as degrees
        in 0 - 180 range

        Example:
          >>> Vect(1, 0, 0).angle(Vect(0, 0, 1))
          90.0
          >>> Vect(1, 0, 0).angle(Vect(-1, 0, 0))
          180.0
          >>> Vect(0, 0, 1).angle(Vect(0, 0, -1))
          180.0
          >>> Vect(1, 1, 1).angle(Vect(1, 1,1 ))
          0.0
        """

        return degrees(acos(self.cos_angle(another)))

    def vp(self, another):
        """
        Vector product

        Examples:
          >>> Vect(1, 0, 0).vp(Vect(0, 1, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> Vect(1, 0, 0).vp(Vect(1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
          >>> Vect(1, 0, 0).vp(Vect(-1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
        """

        return Vect.from_array(np.cross(self.v, another.v))

    def by_matrix(self, array3x3):
        """
        Matrix multiplication of a vector

        """

        return Vect.from_array(array3x3.dot(self.v))


class GVect(object):
    """
    Geological vector.
    Defined by trend and plunge (both in degrees)
     - trend: [0.0, 360.0[ clockwise, from 0 (North)
     - plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """

    def __init__(self, srcTrend, srcPlunge):
        """
        Geological vector constructor.
        srcTrend: Trend range: [0.0, 360.0[ clockwise, from 0 (North)
        srcPlunge: Plunge: [-90.0, 90.0], negative value: upward pointing axis, positive values: downward axis

        Example:
          >>> a = GVect(120, -27)
          >>> b = GVect(54, -320)
          Traceback (most recent call last):
          ...
          AssertionError: plunge must be between -90° and +90° (comprised)
        """

        assert -90.0 <= srcPlunge <= 90.0, "plunge must be between -90° and +90° (comprised)"
        self._trend = srcTrend % 360.0
        self._plunge = float(srcPlunge)

    @property
    def tr(self):
        """
        Returns trend of the geological direction.
        Range is [0, 360[

        Example:
          >>> GVect(420, -17).tr
          60.0
          >>> GVect(-20, 49).tr
          340.0
        """

        return self._trend

    @property
    def tp(self):
        """
        Returns trend and plunge of the geological direction

        Example:
          >>> GVect(-90, -45).tp
          (270.0, -45.0)
        """

        return self.tr, self.pl

    @property
    def pl(self):
        """
        Returns plugne of the geological direction.
        Range is [-90, 90]

        Example:
          >>> GVect(420, -17).pl
          -17.0
        """

        return self._plunge

    def __repr__(self):

        return "GVect({:.2f}, {:.2f})".format(*self.tp)

    @property
    def versor(self):
        """
        Returns the Vect corresponding to the geological vector

        Examples:
          >>> print GVect(0, 90).versor
          Vect(0.0000, 0.0000, -1.0000)
          >>> print GVect(0, -90).versor
          Vect(0.0000, 0.0000, 1.0000)
        """

        north_coord = cos(radians(self.pl)) * cos(radians(self.tr))
        east_coord = cos(radians(self.pl)) * sin(radians(self.tr))
        down_coord = sin(radians(self.pl))

        return Vect(east_coord, north_coord, -down_coord)

    @property
    def dwngvect(self):
        """
        Return downward-point geological vector

        Examples:
          >>> GVect(90, -45).dwngvect
          GVect(270.00, 45.00)
          >>> GVect(180, 45).dwngvect
          GVect(180.00, 45.00)
          >>> GVect(0, 0).dwngvect
          GVect(0.00, 0.00)
          >>> GVect(0, 90).dwngvect
          GVect(0.00, 90.00)
        """

        trend, plunge = self.tr, self.pl
        if plunge < 0.0:
            trend = (trend + 180.0) % 360.0
            plunge = - plunge

        return GVect(trend, plunge)

    @property
    def ngplane(self):
        """
        Return the geological plane that is normal to the geological vector.

        Examples:
          >>> GVect(0, 45).ngplane
          GPlane(180.00, 45.00)
          >>> GVect(0, -45).ngplane
          GPlane(0.00, 45.00)
          >>> GVect(0, 90).ngplane
          GPlane(180.00, 0.00)
        """

        down_axis = self.dwngvect
        dipdir = (down_axis.tr + 180.0) % 360.0
        dipangle = 90.0 - down_axis.pl

        return GPlane(dipdir, dipangle)


class Plane(object):
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

        return Vect(self.a, self.b, self.c).versor

    def as_geolplane_and_point_3d(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.normal_versor3d.gvect.ngplane()
        point = Point(point_solution(np.array([[self.a, self.b, self.c]]),
                                     np.array([-self.d])))
        return geol_plane, point

    def intersection_versor3d(self, another):
        """
        return intersection versor for two intersecting planes
        """

        return self.normal_versor3d.vp(another.normal_versor3d).versor

    def intersection_point3dt(self, another):
        """
        return point on intersection line (obviously non-unique solution)
        for two planes
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point(x, y, z)

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

    def __repr__(self):

        return "GPlane({:.2f}, {:.2f})".format(*self.dda)

    @property
    def normal(self):
        """
        Returns the normal to the plane, as a geological (downward) axis.

        Example:
            >>> ga = GPlane(90, 55).normal
            >>> ga.tr
            270.0
            >>> ga.pl
            35.0
        """
        
        trend = (self.dd + 180.0) % 360.0
        plunge = 90.0 - self.da

        return GVect(trend, plunge)

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

        normal_versor = self.normal.dwngvect.versor
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return Plane(a, b, c, d)

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
