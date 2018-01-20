# -*- coding: utf-8 -*-

from math import radians, sin, cos

import numpy as np

from .geometry import Vect, GVect, GAxis


class RotationAxis(object):
    """
    Rotation axis, expressed by a geological vector and a rotation angle.
    """

    def __repr__(self):

        return "RotationAxis({:.4f}, {:.4f}, {:.4f})".format(self.gv.tr, self.gv.pl, self.a)

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

        self.gv = GVect(trend, plunge)
        self.a = float(rot_ang)

    @classmethod
    def from_gvect(cls, gvector, angle):
        """
        Class constructor from a GVect instance and an angle value.

        :param gvector: a GVect instance
        :param angle: float value
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.from_gvect(GVect(320, 12), 30)
          RotationAxis(320.0000, 12.0000, 30.0000)
          >>> RotationAxis.from_gvect(GVect(315.0, -0.0), 10)
          RotationAxis(315.0000, -0.0000, 10.0000)
        """

        obj = cls()
        obj.gv = gvector
        obj.a = float(angle)

        return obj

    @property
    def rot_ang(self):
        """
        Returns the rotation angle of the rotation axis.

        :return: rotation angle (Float)

        Example:
          >>> RotationAxis(10, 15, 230).rot_ang
          230.0
        """

        return self.a

    @property
    def rot_vect(self):
        """
        Returns the rotation axis, expressed as a GVect.

        :return: GVect instance

        Example:
          >>> RotationAxis(320, 40, 15).rot_vect
          GVect(320.00, +40.00)
          >>> RotationAxis(135, 0, -10).rot_vect
          GVect(135.00, +00.00)
          >>> RotationAxis(45, 10, 10).rot_vect
          GVect(045.00, +10.00)
        """

        return self.gv

    @classmethod
    def from_vect(cls, vector, angle):
        """
        Class constructor from a Vect instance and an angle value.

        :param vector: a Vect instance
        :param angle: float value
        :return: RotationAxis instance

        Example:
          >>> RotationAxis.from_vect(Vect(0, 1, 0), 30)
          RotationAxis(0.0000, 0.0000, 30.0000)
          >>> RotationAxis.from_vect(Vect(1, 0, 0), 30)
          RotationAxis(90.0000, 0.0000, 30.0000)
          >>> RotationAxis.from_vect(Vect(0, 0, -1), 30)
          RotationAxis(0.0000, 90.0000, 30.0000)
        """

        obj = cls()
        obj.gv = vector.gvect()
        obj.a = angle

        return obj

    @property
    def versor(self):
        """
        Return the versor equivalent to the Rotation geological vector.

        :return: Vect
        """

        return self.gv.versor()

    def specular(self):
        """
        Derives the rotation axis with opposite vector orientation
        and rotation angle that is the complement to 360°.

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

        return RotationAxis.from_gvect(gvect_opp, opposite_angle)

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



if __name__ == "__main__":

    import doctest
    doctest.testmod()
