# -*- coding: utf-8 -*-

from __future__ import division

from math import sqrt, sin, cos, tan, radians, asin, acos, atan, atan2, degrees
import numpy as np

from apsg import Vec3

from .array_utils import point_solution


MINIMUM_SEPARATION_THRESHOLD = 1e-10
MINIMUM_VECTOR_MAGNITUDE = 1e-10


class CartesianPoint3DT(object):

    def __init__(self, x=np.nan, y=np.nan, z=np.nan, t=None):

        self._x = x
        self._y = y
        self._z = z
        self._t = t

    @property
    def p_x(self):

        return self._x

    @property
    def p_y(self):

        return self._y

    @property
    def p_z(self):

        return self._z

    @property
    def p_t(self):

        return self._t

    def clone(self):

        return CartesianPoint3DT(self.p_x, self.p_y, self.p_z, self.p_t)

    def spat_distance(self, another):
        """
        Calculate Euclidean spatial distance between two points.

        @param  another:  the CartesianPoint3DT instance for which the spatial distance should be calculated
        @type  another:  CartesianPoint3DT.

        @return:  spatial distance between the two points - float.
        """

        return sqrt((self.p_x - another.p_x) ** 2 + (self.p_y - another.p_y) ** 2 + (self.p_z - another.p_z) ** 2)

    def distance_2d(self, another):

        return sqrt((self.p_x - another.p_x) ** 2 + (self.p_y - another.p_y) ** 2)

    def spat_coincident_with(self, another, tolerance=MINIMUM_SEPARATION_THRESHOLD):

        if self.spat_distance(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        @param  sx:  the shift to be applied along the x axis.
        @type  sx:  float.
        @param  sy:  the shift to be applied along the y axis.
        @type  sy:  float.
        @param  sz:  the shift to be applied along the z axis.
        @type  sz:  float.

        @return:  a new CartesianPoint3DT instance shifted by the given amounts with respect to the original one.
        """

        return CartesianPoint3DT(self.p_x + sx, self.p_y + sy, self.p_z + sz, self.p_t)

    def translate_with_vector(self, displacement_vector):

        return CartesianPoint3DT(self.p_x + displacement_vector.x, self.p_y + displacement_vector.y,
                                 self.p_z + displacement_vector.z, self.p_t)

    def as_vector3d(self):

        return CartesianVector3D(self.p_x, self.p_y, self.p_z)

    def delta_time(self, another):

        return another.p_t - self.p_t

    def speed(self, another):

        try:
            return self.spat_distance(another) / self.delta_time(another)
        except:
            return np.nan


class CartesianVector3D(object):
    def __init__(self, x=np.nan, y=np.nan, z=np.nan):

        self._x = x
        self._y = y
        self._z = z

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def z(self):

        return self._z

    def clone(self):

        return CartesianVector3D(self.x, self.y, self.z)

    @property
    def length(self):

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def length_horiz(self):

        return sqrt(self.x * self.x + self.y * self.y)

    def scale(self, scale_factor):

        return CartesianVector3D(self.x * scale_factor,
                                 self.y * scale_factor,
                                 self.z * scale_factor)

    def as_versor3d(self):

        return self.scale(1.0 / self.length)

    def as_downvector3d(self):

        if self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.clone()

    def add(self, another):

        return CartesianVector3D(self.x + another.x,
                                 self.y + another.y,
                                 self.z + another.z)

    def slope_radians(self):

        return atan(self.z / self.length_horiz)

    def as_geolaxis(self):

        if self.length < MINIMUM_VECTOR_MAGNITUDE:
            return None

        unit_vect = self.as_versor3d()

        plunge = - degrees(asin(unit_vect.z))  # upward negative, downward positive

        trend = 90.0 - degrees(atan2(unit_vect.y, unit_vect.x))
        if trend < 0.0:
            trend += 360.0
        elif trend > 360.0:
            trend -= 360.0

        assert 0.0 <= trend < 360.0
        assert -90.0 <= plunge <= 90.0

        return GeolAxis(trend, plunge)

    def scalar_product(self, another):

        return self.x * another.x + self.y * another.y + self.z * another.z

    def vectors_cos_angle(self, another):

        try:
            val = self.scalar_product(another) / (self.length * another.length)
            if val > 1.0:
                return 1.0
            elif val < -1.0:
                return -1.0
            else:
                return val
        except ZeroDivisionError:
            return np.nan

    def angle_degr(self, another):
        """
        angle between two vectors,
        in 0 - pi range
        """

        return degrees(acos(self.vectors_cos_angle(another)))

    def vector_product(self, another):

        x = self.y * another.z - self.z * another.y
        y = self.z * another.x - self.x * another.z
        z = self.x * another.y - self.y * another.x

        return CartesianVector3D(x, y, z)

    def by_matrix(self, matrix3x3):

        vx = matrix3x3[0, 0] * self.x + matrix3x3[0, 1] * self.y + matrix3x3[0, 2] * self.z
        vy = matrix3x3[1, 0] * self.x + matrix3x3[1, 1] * self.y + matrix3x3[1, 2] * self.z
        vz = matrix3x3[2, 0] * self.x + matrix3x3[2, 1] * self.y + matrix3x3[2, 2] * self.z

        return CartesianVector3D(vx, vy, vz)


class GeolAxis(object):
    """
    Structural axis,
    defined by trend and plunge (both in degrees)
    Trend range: [0.0, 360.0[ clockwise, from 0 (North) 
    Plunge: [-90.0, 90.0], negative value: upward axis, positive values: downward axis
    """

    def __init__(self, srcTrend, srcPlunge):

        assert 0.0 <= srcTrend < 360.0
        assert -90.0 <= srcPlunge <= 90.0

        self._trend = srcTrend
        self._plunge = srcPlunge

    @property
    def vals(self):
        
        return self._trend, self._plunge
    
    @property
    def trend(self):

        return self._trend

    @property
    def plunge(self):

        return self._plunge

    def versor_3d(self):

        north_coord = cos(radians(self.plunge)) * cos(radians(self.trend))
        east_coord = cos(radians(self.plunge)) * sin(radians(self.trend))
        down_coord = sin(radians(self.plunge))

        return CartesianVector3D(east_coord, north_coord, -down_coord)

    def as_downgeolaxis(self):

        trend, plunge = self.trend, self.plunge
        if plunge < 0.0:
            trend = (trend + 180.0) % 360.0
            plunge = - plunge

        return GeolAxis(trend, plunge)

    def as_normalgeolplane(self):

        down_axis = self.as_downgeolaxis()
        dipdir = (down_axis.trend + 180.0) % 360.0
        dipangle = 90.0 - down_axis.plunge

        return GeolPlane(dipdir, dipangle)


class CartesianPlane(object):
    """
    Cartesian plane, expressed by equation:
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
        matr_a = np.array([[pt1.p_y, pt1.p_z, 1],
                           [pt2.p_y, pt2.p_z, 1],
                           [pt3.p_y, pt3.p_z, 1]])

        matr_b = - np.array([[pt1.p_x, pt1.p_z, 1],
                             [pt2.p_x, pt2.p_z, 1],
                             [pt3.p_x, pt3.p_z, 1]])

        matr_c = np.array([[pt1.p_x, pt1.p_y, 1],
                           [pt2.p_x, pt2.p_y, 1],
                           [pt3.p_x, pt3.p_y, 1]])

        matr_d = - np.array([[pt1.p_x, pt1.p_y, pt1.p_z],
                             [pt2.p_x, pt2.p_y, pt2.p_z],
                             [pt3.p_x, pt3.p_y, pt3.p_z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    def normal_versor3d(self):
        """
        return the normal versor to the cartesian plane
        """

        return CartesianVector3D(self.a, self.b, self.c).as_versor3d()

    def as_geolplane_and_point_3d(self):
        """
        converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution)
        """

        geol_plane = self.normal_versor3d().as_geolaxis().as_normalgeolplane()
        point = CartesianPoint3DT(point_solution(np.array([[self.a, self.b, self.c]]),
                                                 np.array([-self.d])))
        return geol_plane, point

    def intersection_versor3d(self, another):
        """
        return intersection versor for two intersecting planes
        """

        return self.normal_versor3d().vector_product(another.normal_versor3d()).as_versor3d()

    def intersection_point3dt(self, another):
        """
        return point on intersection line (obviously non-unique solution)
        for two planes
        """

        # find a point lying on the intersection line (this is a non-unique solution)    
        a = np.array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = np.array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return CartesianPoint3DT(x, y, z)

    def set_point_inside(self, pt):
        return self.a * pt.p_x + self.b * pt.p_y + self.c * pt.p_z + self.d

    def angle_degr(self, another):
        angle_degr = self.normal_versor3d().angle_degr(another.normal_versor3d())

        assert angle_degr > 0.0

        if angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr


def rhrstrike2dipdir(rhr_strk):

    return (rhr_strk + 90.0) % 360.0


class GeolPlane(object):
    """
    Structural plane, following geological conventions:
    dip direction and dip angle.
    
    """

    def __init__(self, srcAzimuth, srcDipAngle, isRHRStrike=False):
        """
        Class constructor

        @param  srcAzimuth:  Azimuth of the plane (RHR strike or dip direction).
        @type  srcAzimuth:  number or string convertible to float.
        @param  srcDipAngle:  Dip angle of the plane (0-90�).
        @type  srcDipAngle:  number or string convertible to float.
           
        @return:  GeolPlane.
    
        """

        if isRHRStrike:
            self._dipdir = rhrstrike2dipdir(srcAzimuth)
        else:
            self._dipdir = srcAzimuth % 360.0
        self._dipangle = srcDipAngle

    @property
    def vals(self):
        
        return self._dipdir, self._dipangle
    
    @property
    def dipdir(self):
        
        return self._dipdir

    @property
    def dipangle(self):
        
        return self._dipangle

    def as_normalgeolaxis(self):
        
        trend = (self.dipdir + 180.0) % 360.0
        plunge = 90.0 - self.dipangle

        return GeolAxis(trend, plunge)

    def plane_x_coeff(self):
        """
        Calculate the slope of a given plane along the x direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.    
        """
        return - sin(radians(self.dipdir)) * tan(radians(self.dipangle))

    def plane_y_coeff(self):
        """
        Calculate the slope of a given plane along the y direction.
        The plane orientation  is expressed following the geological convention. 
               
        @return:  slope - float.     
        """
        return - cos(radians(self.dipdir)) * tan(radians(self.dipangle))

    def plane_from_geo(self, or_Pt):
        """
        Closure that embodies the analytical formula for a given, non-vertical plane.
        This closure is used to calculate the z value from given horizontal coordinates (x, y).
    
        @param  or_Pt:  CartesianPoint3DT instance expressing a location point contained by the plane.
        @type  or_Pt:  CartesianPoint3DT.
        
        @return:  lambda (closure) expressing an analytical formula for deriving z given x and y values.
        """

        x0 = or_Pt.p_x
        y0 = or_Pt.p_y
        z0 = or_Pt.p_z

        # slope of the line parallel to the x axis and contained by the plane
        a = self.plane_x_coeff()

        # slope of the line parallel to the y axis and contained by the plane
        b = self.plane_y_coeff()

        return lambda x, y: a * (x - x0) + b * (y - y0) + z0

    def as_cartesplane(self, point):
        normal_versor = self.as_normalgeolaxis().as_downgeolaxis().versor_3d()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.p_x + b * point.p_y + c * point.p_z)
        return CartesianPlane(a, b, c, d)

    def angle_degr(self, another):
        """
        calculate angle (in degrees) between two planes

        >>> p1 = GeolPlane(100.0, 50.0)
        >>> p1.angle_degr(p1)
        0.0

        >>> p2 = GeolPlane(0.0, 0.0)
        >>> p3 = GeolPlane(300.0, 90.0)
        >>> p2.angle_degr(p3)
        90.0
        """

        vec0 = self.as_normalgeolaxis().versor_3d()
        vec1 = another.as_normalgeolaxis().versor_3d()

        return vec0.angle_degr(vec1)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
