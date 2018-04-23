# -*- coding: utf-8 -*-

from pygsf.orientations.orientations import *


class Slick(Orien):
    """
    Slickeline.
    It can be defined through a OrienM instance, in which case it has a movement sense,
    or via a GAxis, when the movement sense is unknown or not sure.
    When the movement sense is known, the OrienM instance indicates the displacement of the block that is:
    - for a horizontal or a dipping, non vertical fault: the upper block
    - for a vertical fault: the block individuated by the (formal) dip direction.
    """

    def __init__(self, trend: float, plunge: float, known=True):
        """"
        Class constructors from trend, plunge and optional known movement sense flag.

        :param trend: the trend of the slickenline
        :type trend: float
        :param plunge: the slickenline plunge
        :type plunge: float
        :param known: the known movement sense flag
        :type known: bool
        :return: the Slick instance

        Example:
          >>> Slick(90, 10)
          Slick(090.00, +10.00, known_dir: True)
          >>> Slick(90, 10, known=False)
          Slick(090.00, +10.00, known_dir: False)
          >>> Slick("90", 10, False)
          Traceback (most recent call last):
          ...
          SlickInputTypeException: Trend must be a number
        """

        if not isinstance(trend, (int, float)):
            raise SlickInputTypeException("Trend must be a number")
        if not isinstance(plunge, (int, float)):
            raise SlickInputTypeException("Plunge must be a number")
        if not isinstance(known, bool):
            raise SlickInputTypeException("Known movement sense must be a boolean")

        super().__init__(trend, plunge, is_axis=known)

    def __repr__(self):

        return "Slick({:06.2f}, {:+06.2f}, known_dir: {})".format(self.tr, self.pl, self.ax)

    def has_known_sense(self):
        """
        Check whether the slickenline has known movement sense.

        Example:
          >>> Slick(90, 45).has_known_sense()
          True
          >>> Slick(90, 45, False).has_known_sense()
          False
        """

        return self.ax

    def has_unknown_sense(self):
        """
        Check whether the slickenline has unknown/uncertain movement sense.

        Example:
          >>> Slick(90, 45, False).has_unknown_sense()
          True
          >>> Slick(90, 45).has_unknown_sense()
          False
        """

        return not self.has_known_sense()

    def set_known_sense(self):
        """
        Return a new slickenline with movement sense set to known.

        Example:
          >>> Slick(180, -30, False).set_known_sense()
          Slick(180.00, -30.00, known_dir: True)
        """

        return Slick(*self.tp)

    def set_unknown_sense(self):
        """
        Set to unknown/uncertain the movement sense for the current Slickline instance.

        Example:
          >>> Slick(180, -30).set_unknown_sense()
          Slick(180.00, -30.00, known_dir: False)
        """

        return Slick(*self.tpoa)

    def invert(self):
        """
        Invert the slickenline sense, when known, otherwise raise SlickSenseException.

        Example:
         >>> Slick(30, 45, False).invert()
         Traceback (most recent call last):
         ...
         SlickSenseException: Slick must have know movement sense
         >>> Slick(30, 45).invert()
         Slick(210.00, -45.00, known_dir: True)
        """

        if not self.has_known_sense():
            raise SlickSenseException("Slick must have know movement sense")

        return Slick(*self.opposite().tp)


class GFault(PPlane):
    """
    Represent a fault plane, represented by a PPlane instance,
    and zero, one or more slickenlines, represented by a list of Slick instances (None when no slickenlines).
    """

    def __init__(self, azim: float, dip_ang: float, is_rhr_strike=False, slickenlines=None):
        """
        Create an instance of a GFault.

        :param  azim:  azimuth of the plane (RHR strike or dip direction).
        :type  azim:  number or string convertible to float.
        :param  dip_ang:  Dip angle of the plane (0-90°).
        :type  dip_ang:  number or string convertible to float.
        :param is_rhr_strike: if the source azimuth is RHR strike (default is False, i.e. it is dip direction)
        :return: the instantiated geological plane.

        Example:
          >>> GFault(90, 45, [Slick(90, 45)])
          GFault(090.00, +45.00) with 1 slickenline(s))
          >>> GFault(90, 45. [Slick(80, 45)])
          Traceback (most recent call last):
          ...
          GFaultInputTypeException: Slick is not within fault plane
        """

        if not isinstance(azim, (int, float)):
            raise GFaultInputTypeException("Azim must be a number")
        if not isinstance(dip_ang, (int, float)):
            raise GFaultInputTypeException("Dip angle must be a number")

        if slickenlines is None:
            slickenlines = []
        temp_gplane = PPlane(azim, dip_ang, is_rhr_strike)
        for slickenline in slickenlines:
            if not temp_gplane.isAlmostParallel(slickenline):
                raise GFaultInputTypeException("All slickenlines must lie on the plane")

        super().__init__(azim, dip_ang, is_rhr_strike)

        self.slickenlines = slickenlines

    def __repr__(self):

        return "GFault({}, {}) with {} slickenine(s)".format(*self.dda, len(self.slick))


class SlickSenseException(Exception):
    """
    Exception for slickenline movement sense.
    """

    pass


class SlickInputTypeException(Exception):
    """
    Exception for slickenline type.
    """

    pass


class GFaultInputTypeException(Exception):
    """
    Exception for GFault input type.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
