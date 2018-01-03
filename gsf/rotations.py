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
        self.a = rot_ang % 360.0

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
        """

        obj = cls()
        obj.gv = gvector
        obj.a = angle

        return obj

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
    import numtest  # external module, used in doctest float checks
    doctest.testmod()
