# -*- coding: utf-8 -*-


from typing import List

from pygsf.mathematics.quaternions import *
from .orientations import *


class RotationAxis(object):
    """
    Rotation axis, expressed by an Orientation and a rotation angle.
    """

    def __repr__(self):

        return "RotationAxis({:.4f}, {:.4f}, {:.4f})".format(*self.gv.d, self.a)

    def __init__(self, trend=float('nan'), plunge=float('nan'), rot_ang=float('nan')):
        """
        Constructor.

        :param trend: Float/Integer
        :param plunge: Float/Integer
        :param rot_ang: Float/Integer

        Example:
        >> RotationAxis(0, 90, 120)
        RotationAxis(0.0000, 90.0000, 120.0000)
        """

        self.gv = Direct.fromAzPl(trend, plunge)
        self.a = float(rot_ang)

    @classmethod
    def fromQuaternion(cls, quat: Quaternion):
        """
        Calculates the Rotation Axis expressed by a quaternion.
        The resulting rotation asVect is set to point downward.
        Examples are taken from Kuipers, 2002, chp. 5.

        :return: RotationAxis instance.

        Examples:
          >>> RotationAxis.fromQuaternion(Quaternion(0.5, 0.5, 0.5, 0.5))
          RotationAxis(45.0000, -35.2644, 120.0000)
          >>> RotationAxis.fromQuaternion(Quaternion(sqrt(2)/2, 0.0, 0.0, sqrt(2)/2))
          RotationAxis(0.0000, -90.0000, 90.0000)
          >>> RotationAxis.fromQuaternion(Quaternion(sqrt(2)/2, sqrt(2)/2, 0.0, 0.0))
          RotationAxis(90.0000, 0.0000, 90.0000)
        """

        if abs(quat) < quat_magn_thresh:

            rot_ang = 0.0
            rot_gvect = Direct(0.0, 0.0)

        elif are_close(quat.scalar, 1):

            rot_ang = 0.0
            rot_gvect = Direct(0.0, 0.0)

        else:

            unit_quat = quat.normalize()
            rot_ang = unit_quat.rotation_angle()
            rot_gvect = Direct.fromVect(unit_quat.vector)

        obj = cls()
        obj.gv = rot_gvect
        obj.a = rot_ang

        return obj

    @classmethod
    def fromDirect(cls, direct: Direct, angle: float):
        """
        Class constructor from a Direct instance and an angle value.

        :param gvector: a Direct instance
        :param angle: float value
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.fromDirect(Direct(320, 12), 30)
          RotationAxis(320.0000, 12.0000, 30.0000)
          >>> RotationAxis.fromDirect(Direct(315.0, -0.0), 10)
          RotationAxis(315.0000, -0.0000, 10.0000)
        """

        obj = cls()
        obj.gv = direct
        obj.a = float(angle)

        return obj

    @property
    def rot_ang(self) -> float:
        """
        Returns the rotation angle of the rotation axis.

        :return: rotation angle (Float)

        Example:
          >>> RotationAxis(10, 15, 230).rot_ang
          230.0
        """

        return self.a

    @property
    def rot_vect(self) -> Direct:
        """
        Returns the rotation axis, expressed as a Direct.

        :return: Direct instance

        Example:
          >>> RotationAxis(320, 40, 15).rot_vect
          Direct(az: 320.00°, pl: 40.00°)
          >>> RotationAxis(135, 0, -10).rot_vect
          Direct(az: 135.00°, pl: 0.00°)
          >>> RotationAxis(45, 10, 10).rot_vect
          Direct(az: 45.00°, pl: 10.00°)
        """

        return self.gv

    @classmethod
    def fromVect(cls, vector: Vect, angle: float):
        """
        Class constructor from a Vect instance and an angle value.

        :param vector: a Vect instance
        :param angle: float value
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.fromVect(Vect(0, 1, 0), 30)
          RotationAxis(0.0000, 0.0000, 30.0000)
          >>> RotationAxis.fromVect(Vect(1, 0, 0), 30)
          RotationAxis(90.0000, 0.0000, 30.0000)
          >>> RotationAxis.fromVect(Vect(0, 0, -1), 30)
          RotationAxis(0.0000, 90.0000, 30.0000)
        """

        obj = cls()
        obj.gv = vector.asOrien()
        obj.a = angle

        return obj

    @property
    def versor(self) -> Vect:
        """
        Return the versor equivalent to the Rotation geological asVect.

        :return: Vect
        """

        return self.gv.versor()

    def specular(self):
        """
        Derives the rotation axis with opposite asVect direction
        and rotation angle that is the complement to 360°.
        The resultant rotation is equivalent to the original one.

        :return: RotationAxis instance.

        Example
          >>> RotationAxis(90, 45, 320).specular()
          RotationAxis(270.0000, -45.0000, 40.0000)
          >>> RotationAxis(135, 0, -10).specular()
          RotationAxis(315.0000, -0.0000, 10.0000)
          >>> RotationAxis(45, 10, 10).specular()
          RotationAxis(225.0000, -10.0000, 350.0000)
        """

        gvect_opp = self.rot_vect.opposite()
        opposite_angle = (360.0 - self.rot_ang) % 360.0

        return RotationAxis.fromDirect(gvect_opp, opposite_angle)

    def compl180(self):
        """
        Creates a new rotation axis that is the complement to 180 of the original one.

        :return: RotationAxis instance.

        Example:
          >>> RotationAxis(90, 45, 120).compl180()
          RotationAxis(90.0000, 45.0000, 300.0000)
          >>> RotationAxis(117, 34, 18).compl180()
          RotationAxis(117.0000, 34.0000, 198.0000)
          >>> RotationAxis(117, 34, -18).compl180()
          RotationAxis(117.0000, 34.0000, 162.0000)
        """

        rot_ang = - (180.0 - self.rot_ang) % 360.0
        return RotationAxis.fromDirect(self.gv, rot_ang)

    def strictly_equival(self, another, angle_tolerance: float = VECTOR_ANGLE_THRESHOLD) -> bool:
        """
        Checks if two RotationAxis are almost equal, based on a strict checking
        of the Direct component and of the rotation angle.

        :param another: another RotationAxis instance, to be compared with
        :type another: RotationAxis
        :return: the equivalence (true/false) between the two compared RotationAxis
        :rtype: bool

        Examples:
          >>> ra_1 = RotationAxis(180, 10, 10)
          >>> ra_2 = RotationAxis(180, 10, 10.5)
          >>> ra_1.strictly_equival(ra_2)
          True
          >>> ra_3 = RotationAxis(180.2, 10, 10.4)
          >>> ra_1.strictly_equival(ra_3)
          True
          >>> ra_4 = RotationAxis(184.9, 10, 10.4)
          >>> ra_1.strictly_equival(ra_4)
          False
        """

        if not self.gv.isAlmostParallel(another.gv, angle_tolerance):
            return False

        if not are_close(self.a, another.a, atol = 1.0):
            return False

        return True

    def to_rotation_quaternion(self) -> Quaternion:
        """
        Converts the RotationAxis instance to the corresponding rotation quaternion.

        :return: the rotation quaternion.
        :rtype: Quaternion
        """

        rotation_angle_rad = radians(self.a)
        rotation_vector = self.gv.versor()

        w = cos(rotation_angle_rad / 2.0)
        x, y, z = rotation_vector.scale(sin(rotation_angle_rad / 2.0)).toXYZ()

        return Quaternion(w, x, y, z).normalize()

    def to_rotation_matrix(self):
        """
        Derives the rotation matrix from the RotationAxis instance.

        :return: 3x3 numpy array
        """

        rotation_versor = self.versor
        phi = radians(self.a)

        l = rotation_versor.x
        m = rotation_versor.y
        n = rotation_versor.z

        cos_phi = cos(phi)
        sin_phi = sin(phi)

        a11 = cos_phi + ((l * l) * (1 - cos_phi))
        a12 = ((l * m) * (1 - cos_phi)) - (n * sin_phi)
        a13 = ((l * n) * (1 - cos_phi)) + (m * sin_phi)

        a21 = ((l * m) * (1 - cos_phi)) + (n * sin_phi)
        a22 = cos_phi + ((m * m) * (1 - cos_phi))
        a23 = ((m * n) * (1 - cos_phi)) - (l * sin_phi)

        a31 = ((l * n) * (1 - cos_phi)) - (m * sin_phi)
        a32 = ((m * n) * (1 - cos_phi)) + (l * sin_phi)
        a33 = cos_phi + ((n * n) * (1 - cos_phi))

        return np.array([(a11, a12, a13),
                         (a21, a22, a23),
                         (a31, a32, a33)])

    def to_min_rotation_axis(self):
        """
        Calculates the minimum rotation axis from the given quaternion.

        :return: RotationAxis instance.
        """

        return self if abs(self.rot_ang) <= 180.0 else self.specular()


def sort_rotations(rotation_axes: List[RotationAxis]) -> List[RotationAxis]:
    """
    Sorts a list or rotation axes, based on the rotation angle (absolute value),
    in an increasing order.

    :param rotation_axes: o list of RotationAxis objects.
    :return: the sorted list of RotationAxis

    Example:
      >>> rots = [RotationAxis(110, 14, -23), RotationAxis(42, 13, 17), RotationAxis(149, 87, 13)]
      >>> sort_rotations(rots)
      [RotationAxis(149.0000, 87.0000, 13.0000), RotationAxis(42.0000, 13.0000, 17.0000), RotationAxis(110.0000, 14.0000, -23.0000)]
    """

    return sorted(rotation_axes, key=lambda rot_ax: abs(rot_ax.rot_ang))


def quat_rot_vect(quat: Quaternion, vect: Vect) -> Vect:
    """
    Calculates a rotated solution of a asVect given a normalized quaternion.
    Original formula in Ref. [1].
    Eq.6: R(qv) = q qv q(-1)

    :param quat: a Quaternion instance
    :param vect: a Vect instance
    :return: a rotated Vect instance

    Example:
      >>> q = Quaternion.i()  # rotation of 180° around the x axis
      >>> quat_rot_vect(q, Vect(0, 1, 0))
      Vect(0.0000, -1.0000, 0.0000)
      >>> quat_rot_vect(q, Vect(0, 1, 1))
      Vect(0.0000, -1.0000, -1.0000)
      >>> q = Quaternion.k()  # rotation of 180° around the z axis
      >>> quat_rot_vect(q, Vect(0, 1, 1))
      Vect(0.0000, -1.0000, 1.0000)
      >>> q = Quaternion.j()  # rotation of 180° around the y axis
      >>> quat_rot_vect(q, Vect(1, 0, 1))
      Vect(-1.0000, 0.0000, -1.0000)
    """

    q = quat.normalize()
    qv = Quaternion.from_vect(vect)

    rotated_v = q * (qv * q.inverse)

    return rotated_v.vector


def quat_rot_quat(q1: Quaternion, q2: Quaternion) -> Quaternion:
    """
    Calculates the rotated quaternion given by:
    q_rot = q2 * q1
    Formula (13) in Ref. [1].

    :param q1: Quaternion instance.
    :param q2: Quaternion instance.
    :return: Quaternion instance.
    """

    return q2 * q1


if __name__ == "__main__":
    import doctest

    doctest.testmod()
