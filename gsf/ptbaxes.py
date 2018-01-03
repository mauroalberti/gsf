# -*- coding: utf-8 -*-

from .geometry import *
from .faults import Slickenline, FaultSlick


class PTBAxes(object):
    """
    Represent the triad of P, T and B kinematic axes.
    It can also calculate the M plane.
    """

    def __init__(self, p_axis=None, t_axis=None, known=True):
        """
        Create a new PTBAxes instances, given the two
        P and T axes.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0))
          PTBAxes(P: GAxis(000.00, +00.00), T: GAxis(090.00, +00.00), True)
        """

        self._pax = p_axis
        self._tax = t_axis
        self._known = known

    @classmethod
    def from_faultslick(cls, fault_slick):
        """
        Class method to calculate P-T axes from a FaultSlick instance.
        Return P and T axes and a third Boolean variable,
        indicating whether the P-T derivation is from a slickenline with a known movement sense (True)
        or with unknown/uncertain movement sense (False).

        Example:
          >>> PTBAxes.from_faultslick(FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))))
          PTBAxes(P: GAxis(000.00, -90.00), T: GAxis(090.00, +00.00), True)
        """

        s_versor = fault_slick.sl.lin.versor()
        f_versor = fault_slick.fp.normal().versor()
        known = fault_slick.known_sense

        obj = cls()
        obj._pax = (f_versor - s_versor).gaxis()
        obj._tax = (f_versor + s_versor).gaxis()
        obj._known = known

        return obj

    @property
    def p_axis(self):
        """
        Return the P axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).p_axis
          GAxis(000.00, +00.00)
        """

        return self._pax

    @property
    def t_axis(self):
        """
        Return the T axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).t_axis
          GAxis(090.00, +00.00)
        """

        return self._tax

    @property
    def known(self):
        """
        Indicate if the movement sense is known.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).known
          True
        """

        return self._known

    def __repr__(self):
        return "PTBAxes(P: {}, T: {}, {})".format(self.p_axis, self.t_axis, self.known)

    @property
    def b_axis(self):
        """
        Calculate the B axis.

        Example:
          >>> PTBAxes(GAxis(0, 0), GAxis(90, 0)).b_axis
          GAxis(000.00, +90.00)
        """

        return self.p_axis.vp(self.t_axis).as_axis()

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
    import numtest  # external module, used in doctest float checks
    doctest.testmod()
