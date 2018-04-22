# -*- coding: utf-8 -*-

import numpy as np

from .default_parameters import *
from .arrays import arrays_are_close
from .geography import Vect, GAxis, GVect, GPlane
from .faults import Slick, GFault
from .quaternions import Quaternion


class PTBAxes(object):
    """
    Represent the triad of P, T and B kinematic axes.
    It can also calculate the M plane.
    """

    def __repr__(self):
        return "PTBAxes(P: {}, T: {})".format(
            self._p_versor.asAxis(),
            self._t_versor.asAxis())

    def __init__(self, p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)):
        """
        Create a new PTBAxes instances, given the two
        P and T axes (provided as GAxis instances).
        T and P axes are recalculated to be strictly orthogonal,
        based on fixed T orientation.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00))
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(80, 0))
          Traceback (most recent call last):
          ...
          PTBAxesInputException: P and T axes must be sub-orthogonal
        """

        if not p_axis.isSubOrthogonal(t_axis):
            raise PTBAxesInputException("P and T axes must be sub-orthogonal")

        b_vect = t_axis.normVersor(p_axis)

        self._t_versor = t_axis.as_versor()
        self._p_versor = b_vect.vCross(self._t_versor).versor()

    @classmethod
    def from_faultslick(cls, fault_slick):
        """
        Class method to calculate P-T axes from a FaultSlick instance.
        Return P and T axes and a third Boolean variable,
        indicating whether the P-T derivation is from a slickenline with a known movement sense (True)
        or with unknown/uncertain movement sense (False).

        Example:
          >>> PTBAxes.from_faultslick(GFault(PPlane(90, 45), Slick(GVect(90, 45))))
          PTBAxes(P: GAxis(000.00, -90.00), T: GAxis(090.00, +00.00))
        """

        s_versor = fault_slick.slick.geom.versor()
        f_versor = fault_slick.gplane._normal_orien_frwrd().versor()

        obj = cls()
        obj._p_versor = (f_versor - s_versor).versor()
        obj._t_versor = (f_versor + s_versor).versor()

        return obj

    @classmethod
    def from_vectors(cls, t_vector, p_vector):
        """
        Class method to create a PTBAxes instance from T and P axis vectors.
        Vectors are not required to be normalized but are required to be
        sub-orthogonal.

        :param t_vector: the asVect representing the T axis (Vect instance).
        :param p_vector: the asVect representing the P axis (Vect instance).
        :return: a PTBAxes instance.

        Example:
          >>> PTBAxes.from_vectors(t_vector=Vect(1,0,0), p_vector=Vect(0,1,0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00))
          >>> PTBAxes.from_vectors(t_vector=Vect(0,0,-1), p_vector=Vect(1,0,0))
          PTBAxes(P: GAxis(090.00, +00.00), T: GAxis(000.00, +90.00))
          >>> PTBAxes.from_vectors(t_vector=Vect(1,1,0), p_vector=Vect(-1,1,0))
          PTBAxes(P: GAxis(315.00, +00.00), T: GAxis(045.00, +00.00))
          >>> PTBAxes.from_vectors(t_vector=Vect(1, 1, 0), p_vector=Vect(0.5, 1, 0))
          Traceback (most recent call last):
          ...
          PTBAxesInputException: P and T vectors must be sub-orthogonal
        """

        if not t_vector.isSubOrthogonal(p_vector):
            raise PTBAxesInputException("P and T vectors must be sub-orthogonal")

        t_versor = t_vector.versor()
        p_versor = p_vector.versor()
        b_versor = t_versor.vCross(p_versor).versor()

        obj = cls()
        obj._p_versor = b_versor.vCross(t_versor).versor()
        obj._t_versor = t_versor

        return obj

    @classmethod
    def from_quaternion(cls, quaternion):
        """
        Creates a PTBAxes instance from a given quaternion.
        Formula extracted from eq. 10 in:
        Kagan, Y.Y, 1991. 3-D rotation of double-couple earthquake sources.

        :param quaternion: a Quaternion instance.
        :return:a PTBAxes instance.
        """

        q0, q1, q2, q3 = quaternion.normalize().toXYZ()

        q0q0 = q0*q0
        q0q1 = q0*q1
        q0q2 = q0*q2
        q0q3 = q0*q3

        q1q1 = q1*q1
        q1q2 = q1*q2
        q1q3 = q1*q3

        q2q2 = q2*q2
        q2q3 = q2*q3

        q3q3 = q3*q3

        t1 = q0q0 + q1q1 - q2q2 - q3q3
        t2 = 2*(q1q2 + q0q3)
        t3 = 2*(q1q3 - q0q2)

        p1 = 2*(q1q2 - q0q3)
        p2 = q0q0 - q1q1 + q2q2 - q3q3
        p3 = 2*(q2q3 + q0q1)

        t_vector = Vect(t1, t2, t3)
        p_vector = Vect(p1, p2, p3)

        return PTBAxes.from_vectors(t_vector=t_vector, p_vector=p_vector)

    @property
    def p_versor(self):
        """
        Return the P versor component of the PTBAxes instance.

        :return: P versor instance

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).p_versor
          Vect(-0.0000, 1.0000, 0.0000)
        """

        return self._p_versor

    @property
    def t_versor(self):
        """
        Return the T versor component of the PTBAxes instance.

        :return: T versor instance

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).t_versor
          Vect(1.0000, 0.0000, -0.0000)
        """

        return self._t_versor

    @property
    def b_versor(self):
        """
        Return the B versor component of the PTBAxes instance.

        :return: B versor instance

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).b_versor
          Vect(0.0000, 0.0000, 1.0000)
        """

        return self.t_versor.vCross(self.p_versor)

    @property
    def p_axis(self):
        """
        Return the P axis.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).p_axis
          GAxis(000.00, +00.00)
          >>> PTBAxes(p_axis=GAxis(0, 90), t_axis=GAxis(90, 0)).p_axis
          GAxis(000.00, +90.00)
        """

        return self.p_versor.asAxis()

    @property
    def t_axis(self):
        """
        Return the T axis.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).t_axis
          GAxis(090.00, +00.00)
          >>> PTBAxes(p_axis=GAxis(0, -90), t_axis=GAxis(90, 0)).t_axis
          GAxis(090.00, +00.00)
        """

        return self.t_versor.asAxis()

    @property
    def b_axis(self):
        """
        Calculate the B axis.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).b_axis
          GAxis(000.00, -90.00)
          >>> PTBAxes(p_axis=GAxis(0, 90), t_axis=GAxis(0, 0)).b_axis
          GAxis(270.00, +00.00)
        """

        return self.b_versor.asAxis()

    @property
    def m_plane(self):
        """
        Calculate M plane.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 90), t_axis=GAxis(90, 0)).m_plane.isAlmostParallel(PPlane(0.0, 90.0))
          True
          >>> (PTBAxes(p_axis=GAxis(45, 45), t_axis=GAxis(225, 45)).m_plane).isAlmostParallel(PPlane(315.00, 90.00))
          True
        """

        return self.p_axis.commonGPlane(self.t_axis)

    def almost_equal(self, another, tolerance_angle=VECTOR_ANGLE_THRESHOLD):
        """
        Checks for equivalence between two PTBAXes instances
        within a given tolerance angle (default is VECTOR_ANGLE_THRESHOLD)

        :param another: a PTBAXes instance.
        :param tolerance_angle: the tolerance angle for the equality check (float)
        :return: Boolean.

        Examples:
          >>> fm1 = PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0))
          >>> fm2 = PTBAxes(p_axis=GAxis(0, 0.5), t_axis=GAxis(90, 0))
          >>> fm1.almost_equal(fm2)
          True
          >>> fm3 = PTBAxes(p_axis=GAxis(180.5, 0), t_axis=GAxis(90.5, 0))
          >>> fm1.almost_equal(fm3)
          True
          >>> fm3.almost_equal(fm2)
          True
          >>> fm4 = PTBAxes(p_axis=GAxis(181.5, 0), t_axis=GAxis(91.5, 0))
          >>> fm1.almost_equal(fm4)
          False
        """

        if not self.p_axis.isAlmostParallel(another.p_axis, tolerance_angle):
            return False

        if not self.t_axis.isAlmostParallel(another.t_axis, tolerance_angle):
            return False

        return True

    def to_matrix(self):
        """
        Creates a rotation matrix from the PTB as_vect xyz.
        Formula as in:
        - eq. 3 in Kagan, Y. Y., 2007. Simplified algorithms for calculating double-couple rotation.
        - eq. 31 in Kagan, Y. Y., 2008. On geometric complexity of earthquake focal zone and fault system.

        :return: a 3x3 numpy arrays fo floats.

        Example:
          >>> arrays_are_close(PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).to_matrix(), np.identity(3))
          True
        """

        t = self.t_versor
        p = self.p_versor
        b = self.b_versor

        return np.array([
            [t.x, p.x, b.x],
            [t.y, p.y, b.y],
            [t.z, p.z, b.z]
        ])

    def to_quaternion(self):
        """
        Transforms the focal mechanism into a quaternion.

        :return: a Quaternion instance.

        Example:
          >>> PTBAxes(p_axis=GAxis(232, 41), t_axis=GAxis(120, 24)).to_quaternion()
          Quaternion(-0.41567, 0.85017, -0.31120, -0.08706)
          >>> PTBAxes(p_axis=GAxis(51, 17), t_axis=GAxis(295, 55)).to_quaternion()
          Quaternion(0.38380, 0.30459, 0.80853, -0.32588)
        """

        return Quaternion.from_rot_matr(self.to_matrix())


class PTBAxesInputException(Exception):
    """
    Exception for PTBAxes input
    """


if __name__ == "__main__":

    import doctest
    doctest.testmod()
