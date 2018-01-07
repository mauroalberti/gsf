# -*- coding: utf-8 -*-

from .arrays import arrays_are_close
from .geometry import *
from .faults import Slickenline, FaultSlick
from .quaternions import Quaternion


class PTBAxes(object):
    """
    Represent the triad of P, T and B kinematic axes.
    It can also calculate the M plane.
    """

    def __repr__(self):
        return "PTBAxes(P: {}, T: {})".format(
            self._p_versor.gaxis(),
            self._t_versor.gaxis())

    def __init__(self, p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)):
        """
        Create a new PTBAxes instances, given the two
        P and T axes (provided as GAxis instances).

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00))
        """

        assert p_axis.is_suborthogonal(t_axis)

        self._p_versor = p_axis.as_versor()
        self._t_versor = t_axis.as_versor()

    @classmethod
    def from_faultslick(cls, fault_slick):
        """
        Class method to calculate P-T axes from a FaultSlick instance.
        Return P and T axes and a third Boolean variable,
        indicating whether the P-T derivation is from a slickenline with a known movement sense (True)
        or with unknown/uncertain movement sense (False).

        Example:
          >>> PTBAxes.from_faultslick(FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))))
          PTBAxes(P: GAxis(000.00, -90.00), T: GAxis(090.00, +00.00))
        """

        s_versor = fault_slick.sl.lin.versor()
        f_versor = fault_slick.fp.normal().versor()

        obj = cls()
        obj._p_versor = (f_versor - s_versor).versor()
        obj._t_versor = (f_versor + s_versor).versor()

        return obj

    @property
    def p_versor(self):
        """
        Return the P versor component of the PTBAxes instance.

        :return: P versor instance

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 0), t_axis=GAxis(90, 0)).p_versor
          Vect(0.0000, 1.0000, -0.0000)
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

        return self.t_versor.vp(self.p_versor)

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

        return self.p_versor.gaxis()

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

        return self.t_versor.gaxis()

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

        return self.b_versor.gaxis()

    @property
    def m_plane(self):
        """
        Calculate M plane.

        Example:
          >>> PTBAxes(p_axis=GAxis(0, 90), t_axis=GAxis(90, 0)).m_plane.almost_parallel(GPlane(0.0, 90.0))
          True
          >>> (PTBAxes(p_axis=GAxis(45, 45), t_axis=GAxis(225, 45)).m_plane).almost_parallel(GPlane(315.00, 90.00))
          True
        """

        return self.p_axis.common_plane(self.t_axis)

    def to_matrix(self):
        """
        Creates a rotation matrix from the PTB vector components.
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

    def calculate_rotations(self, another):
        """
        Calculate the rotations between two focal mechanisms, sensu Kagan.
        See Kagan Y.Y. papers for theoretical basis.
        Practical implementation derive from Alberti, 2010:
        Analysis of kinematic correlations in faults and focal mechanisms with GIS and Fortran programs.
        Examples are from Kagan Y. Y., 1991. 3-D rotation of double-couple earthquake sources.

        :param another: a PTBAxes instance
        :return:

        Example:
          >>> fm1 = PTBAxes(p_axis=GAxis(232, 41), t_axis=GAxis(120, 24))
          >>> fm2 = PTBAxes(p_axis=GAxis(51, 17), t_axis=GAxis(295, 55))
          >>> rots = fm1.calculate_rotations(fm2)
        """

        fm1 = self
        fm2 = another

        # processing of equal focal mechanisms
        t_axes_angle = fm1.t_axis.angle(fm2.t_axis)
        p_axes_angle = fm1.p_axis.angle(fm2.p_axis)

        if t_axes_angle < 0.5 and p_axes_angle < 0.5:
          return []

        # transformation of XYZ axes cartesian components (fm1,2) into quaternions q1,2

        focmec1_matrix = fm1.to_matrix()
        focmec2_matrix = fm2.to_matrix()

        fm1_quaternion = Quaternion.from_rot_matr(focmec1_matrix)
        fm2_quaternion = Quaternion.from_rot_matr(focmec2_matrix)

        # calculation of quaternion inverse q1,2[-1]

        fm1_inversequatern = fm1_quaternion.inverse
        fm2_inversequatern = fm2_quaternion.inverse

        # calculation of rotation quaternion : q' = q2*q1[-1]

        base_rot_quater = fm2_quaternion * fm1_inversequatern

        # calculation of secondary rotation pure quaternions: a(i,j,k) = q2*(i,j,k)*q2[-1]
        axes_quats = [
            Quaternion.i(),
            Quaternion.j(),
            Quaternion.k()]

        suppl_prod2quat = map(lambda ax_quat: fm2_quaternion * (ax_quat * fm2_inversequatern), axes_quats)

        # calculation of the other 3 rotation quaternions: q'(i,j,k) = a(i,j,k)*q'

        rotations_quaternions = list(map(lambda quat: quat * base_rot_quater, suppl_prod2quat))
        rotations_quaternions.append(base_rot_quater)
        rotations_axes = map(lambda quat: quat.to_rotation_axis(), rotations_quaternions)

        return sorted(rotations_axes, key=lambda rot_ax: rot_ax.a)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
