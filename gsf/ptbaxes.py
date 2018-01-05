# -*- coding: utf-8 -*-

from .geometry import *
from .faults import Slickenline, FaultSlick


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
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00))
        """

        p_versor = p_axis.as_versor()
        t_versor = t_axis.as_versor()
        assert p_versor.is_suborthogonal(t_versor)

        self._p_versor = p_versor
        self._t_versor = t_versor

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
    def p_axis(self):
        """
        Return the P axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).p_axis
          GAxis(000.00, +00.00)
        """

        return self._p_versor.gaxis()

    @property
    def t_axis(self):
        """
        Return the T axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).t_axis
          GAxis(090.00, +00.00)
        """

        return self._t_versor.gaxis()

    @property
    def b_axis(self):
        """
        Calculate the B axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).b_axis
          GAxis(000.00, +90.00)
        """

        return self._p_versor.vp(self._t_versor).gaxis()

    @property
    def m_plane(self):
        """
        Calculate M plane.

        Example:
          >>> PTBAxes(GAxis(0, 90), GAxis(90, 0)).m_plane.almost_parallel(GPlane(0.0, 90.0))
          True
          >>> (PTBAxes(GAxis(45, 45), GAxis(225, 45)).m_plane).almost_parallel(GPlane(315.00, 90.00))
          True
        """

        return self.p_axis.common_plane(self.t_axis)


if __name__ == "__main__":

    import doctest
    doctest.testmod()
