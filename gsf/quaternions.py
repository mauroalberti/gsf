
from math import sqrt, degrees, acos
import numpy as np

from .geometry import Vect
from .transformations import Rotation
from .errors import QuaternionException


quat_normaliz_tolerance = 1.0e-6
quat_division_tolerance = 1.0e-10
vect_magn_thresh = 1.0e-6


class Quaternion(object):
    """
    Quaternions.
    """

    def __init__(self, w=np.nan, x=np.nan, y=np.nan, z=np.nan):
        """
        Construct a Quaternion instance.

        Example;
          >>> Quaternion(1, 0, 1, 0)
          Quaternion(1.00000, 0.00000, 1.00000, 0.00000)
          >>> Quaternion()
          Quaternion(nan, nan, nan, nan)
        """

        self.q = np.array([w, x, y, z], dtype=np.float64)

    def __repr__(self):

        return "Quaternion({:.5f}, {:.5f}, {:.5f}, {:.5f})".format(self.q[0], self.q[1], self.q[2], self.q[3])

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a quaternion from a numpy 1x4 array.

        Example:
          >>> Quaternion.from_array(np.array([1, 0, 1, 0]))
          Quaternion(1.00000, 0.00000, 1.00000, 0.00000)
        """

        assert a.size == 4

        obj = cls()
        obj.q = a.astype(np.float64)

        return obj

    @classmethod
    def from_cart_matr(cls, matr):
        """
        Class method to construct a quaternion from a 3x3 matrix.
        """

        q0 = sqrt(1 + matr[0, 0] + matr[1, 1] + matr[2, 2]) / 2.0
        q1 = sqrt(1 + matr[0, 0] - matr[1, 1] - matr[2, 2]) / 2.0
        q2 = sqrt(1 - matr[0, 0] + matr[1, 1] - matr[2, 2]) / 2.0
        q3 = sqrt(1 - matr[0, 0] - matr[1, 1] + matr[2, 2]) / 2.0

        q0q1 = (matr[2, 1] - matr[1, 2]) / 4.0
        q0q2 = (matr[0, 2] - matr[2, 0]) / 4.0
        q0q3 = (matr[1, 0] - matr[0, 1]) / 4.0
        q1q2 = (matr[0, 1] + matr[1, 0]) / 4.0
        q1q3 = (matr[0, 2] + matr[2, 0]) / 4.0
        q2q3 = (matr[1, 2] + matr[2, 1]) / 4.0

        if (3 * q0) > (q1 + q2 + q3):

            q1 = q0q1 / q0
            q2 = q0q2 / q0
            q3 = q0q3 / q0

        elif (3 * q1) > (q0 + q2 + q3):

            q0 = q0q1 / q1
            q2 = q1q2 / q1
            q3 = q1q3 / q1

        elif (3 * q2) > (q0 + q1 + q3):

            q0 = q0q2 / q2
            q1 = q1q2 / q2
            q3 = q2q3 / q2

        else:

            q0 = q0q3 / q3
            q1 = q1q3 / q3
            q2 = q2q3 / q3

        w, x, y, z = q0, q1, q2, q3

        obj = cls()
        q = np.array([w, x, y, z], dtype=np.float64)
        obj.q = q

        return obj

    def __eq__(self, another):
        """
        Quaternion equality.

        :param another: a Quaternion instance
        :return: Boolean

        Example:
          >>> Quaternion(1, 1, 3, 0) == Quaternion(0, 7, -2, 4)
          False
          >>> Quaternion(1.0, 1.0, 3.0, 0.0) == Quaternion(1.0, 1.0, 3.0, 0.0)
          True
          >>> Quaternion(1.0, 1.0, 3.0, np.nan) == Quaternion(1.0, 1.0, 3.0, np.nan)
          True
          >>> Quaternion(1.0, 1.0, 3.0, 0.0) == Quaternion(1.0, 1.0, 3.0, -1.0e-20)
          False
        """

        return  ((self.q == another.q) | (np.isnan(self.q) & np.isnan(another.q))).all()

    def __ne__(self, another):
        """
        Quaternion inequality.

        :param another: a Quaternion instance
        :return: Boolean

        Example:
          >>> Quaternion(1, 1, 3, 0) != Quaternion(0, 7, -2, 4)
          True
          >>> Quaternion(1.0, 1.0, 3.0, np.nan) != Quaternion(1.0, 1.0, 3.0, np.nan)
          False
        """

        return not (self == another)

    def __add__(self, another):
        """
        Quaternion sum.

        :param another: Quaternion instance.
        :return: Quaternion instance.

        Example:
          >>> Quaternion(1, 1, 3, 0) + Quaternion(0, 7, -2, 4)
          Quaternion(1.00000, 8.00000, 1.00000, 4.00000)
        """

        return Quaternion.from_array(self.q + another.q)

    def __sub__(self, another):
        """
        Quaternion difference.

        :param another: Quaternion instance.
        :return: Quaternion instance.

        Example:
          >>> Quaternion(1, 1, 3, 0) - Quaternion(0, 7, -2, 4)
          Quaternion(1.00000, -6.00000, 5.00000, -4.00000)
        """

        return Quaternion.from_array(self.q - another.q)

    def scalar_mult(self, val):
        """
        Multiplication of a quaternion by a scalar value.

        :param val: Integer or Float
        :return: Quaternion instance

        Example:
          >>> Quaternion(1, 1, 3, 0).scalar_mult(4)
          Quaternion(4.00000, 4.00000, 12.00000, 0.00000)
        """

        return Quaternion.from_array(self.q * val)

    def quater_mult(self, another):
        """
        Quaternion multiplication.

        :param another: Quaternion instance .
        :return: Quaternion instance.

        Example:
          >>> Quaternion(3, 1, -2, 1).quater_mult(Quaternion(2, -1, 2, 3))
          Quaternion(8.00000, -9.00000, -2.00000, 11.00000)
        """
        
        a = + (self.q[0] * another.q[0]) \
            - (self.q[1] * another.q[1]) \
            - (self.q[2] * another.q[2]) \
            - (self.q[3] * another.q[3])
        
        b = + (self.q[0] * another.q[1]) \
            + (self.q[1] * another.q[0]) \
            + (self.q[2] * another.q[3]) \
            - (self.q[3] * another.q[2])

        c = + (self.q[0] * another.q[2]) \
            - (self.q[1] * another.q[3]) \
            + (self.q[2] * another.q[0]) \
            + (self.q[3] * another.q[1])
        
        d = + (self.q[0] * another.q[3]) \
            + (self.q[1] * another.q[2]) \
            - (self.q[2] * another.q[1]) \
            + (self.q[3] * another.q[0])
                
        return Quaternion(a, b, c, d)

    def __mul__(self, another):
        """
        Wrapper for quaternion multiplication.

        :param another: Quaternion or Float/Integer
        :return: Quaternion instance

        Example:
          >>> Quaternion(1, 1, 3, 0) * 3
          Quaternion(3.00000, 3.00000, 9.00000, 0.00000)
          >>> Quaternion(3, 1, -2, 1) * Quaternion(2, -1, 2, 3)
          Quaternion(8.00000, -9.00000, -2.00000, 11.00000)
          >>> Quaternion(3, 1, -2, 1) * "a"
          Traceback (most recent call last):
          ...
          QuaternionException: Multiplicand is not number or quaternion
        """

        if isinstance(another, (int, long, float)):
            return self.scalar_mult(another)
        elif isinstance(another, Quaternion):
            return self.quater_mult(another)
        else:
            raise QuaternionException("Multiplicand is not number or quaternion")

    def __div__(self, denominator):
        """
        Division of a quaternion by a scalar.

        :param denominator: Float
        :return: Quaternion instance.

        Example:
          >>> Quaternion(1, 1, 3, 0) / 3
          Quaternion(0.33333, 0.33333, 1.00000, 0.00000)
          >>> Quaternion(1, 1, 3, 0) / 1e-11
          Traceback (most recent call last):
          ...
          QuaternionException: Quaternion division by almost zero value
        """

        if abs(denominator) < quat_division_tolerance:
            raise QuaternionException("Quaternion division by almost zero value")
        else:
            return Quaternion.from_array(self.q / denominator)

    def conjugate(self):
        """
        Quaternion conjugate.

        :return: Quaternion instance.

        Example:
          >>> Quaternion(1, 1, 3, 0).conjugate()
          Quaternion(1.00000, -1.00000, -3.00000, -0.00000)
          >>> Quaternion(2.0, 0.0, -3.3, 17.09).conjugate()
          Quaternion(2.00000, -0.00000, 3.30000, -17.09000)
          >>> Quaternion(2.0, 0.0, np.nan, 17.09).conjugate()
          Quaternion(2.00000, -0.00000, nan, -17.09000)
        """
        
        a = + self.q[0]
        b = - self.q[1]
        c = - self.q[2]
        d = - self.q[3]

        return Quaternion(a, b, c, d)

    def sqrd_norm(self):
        """
        Squared norm of a quaternion.

        :return: Float

        Example:
          >>> Quaternion(1, 0, 0, 0).sqrd_norm()
          1.0
          >>> Quaternion(1, 1, 0, 2).sqrd_norm()
          6.0
          >>> Quaternion(2, -1, 2, 3).sqrd_norm()
          18.0
          >>> Quaternion(2, np.nan, 2, 3).sqrd_norm()
          nan
        """

        return self.q[0]**2 + self.q[1]**2 + self.q[2]**2 + self.q[3]**2


    def inverse(self):
        """
        Quaternion inverse.

        :return: Quaternion instance.

        Example:
          >>> Quaternion(0, 1, 0, 0).inverse()
          Quaternion(0.00000, -1.00000, -0.00000, -0.00000)
        """

        return self.conjugate() / self.sqrd_norm()

    def is_normalized(self):
        """
        Check if a quaternion is unitary.

        :return: Boolean

        Example:
          >>> Quaternion(0, 1, 0, 0).is_normalized()
          True
          >>> Quaternion(1, 4, 0, -4).is_normalized()
          False
        """

        return abs(1.0 - sqrt(self.sqrd_norm())) < quat_normaliz_tolerance

    def normalize(self):
        """
        Normalize a quaternion.

        :return: Quaternion instance.

        Example:
          >>> Quaternion(0, 4, 0, 0).normalize()
          Quaternion(0.00000, 1.00000, 0.00000, 0.00000)
          >>> Quaternion(0, 4, 0, 8).normalize()
          Quaternion(0.00000, 0.44721, 0.00000, 0.89443)
        """

        return self / sqrt(self.sqrd_norm())

    def to_rotation(self):
        """
        Calculates the rotation expressed by a quaternion.
        The resulting rotation vector is set to point downward.
        Examples are taken from Kuipers, 2002, chp. 5.

        :return: Rotation instance.

        Examples:
          >>> Quaternion(0.5, 0.5, 0.5, 0.5).to_rotation()
          Rotation(225.0000, 35.2644, 240.0000)
          >>> Quaternion(sqrt(2)/2, 0.0, 0.0, sqrt(2)/2).to_rotation()
          Rotation(180.0000, 90.0000, 270.0000)
          >>> Quaternion(sqrt(2)/2, sqrt(2)/2, 0.0, 0.0).to_rotation()
          Rotation(90.0000, 0.0000, 90.0000)
        """

        norm_quat = self.normalize()

        qvect = Vect.from_array(norm_quat.q[1:])
        qvect_magn = qvect.len_3d

        if qvect_magn < vect_magn_thresh:

            rot_ang = 0.0
            rot_axis_tr = 0.0
            rot_axis_pl = 0.0

        else:

            rot_ang = (2 * degrees(acos(norm_quat.q[0]))) % 360.0
            rot_gvect = qvect.gvect()
            if rot_gvect.is_upward:
                rot_gvect = rot_gvect.downward()
                rot_ang = - rot_ang
            rot_axis_tr = rot_gvect.tr
            rot_axis_pl = rot_gvect.pl

        return Rotation(
            rot_axis_tr,
            rot_axis_pl,
            rot_ang)



if __name__ == "__main__":

    import doctest
    import numtest  # external module, used in doctest float checks
    doctest.testmod()

