# -*- coding: utf-8 -*-

from .geometry import *
from .mathematics import are_close


class Slickenline(object):
    """
    Slickeline.
    It can be defined through a GVect instance, in which case it has a movement sense,
    or via a GAxis, when the movement sense is unknown or not sure.
    When the movement sense is known, the GVect instance indicates the displacement of the block that is:
    - for a horizontal or a dipping, non vertical fault: the upper block
    - for a vertical fault: the block individuated by the (formal) dip direction.
    """

    def __init__(self, mov_lin):
        """"
        Slickenline constructor.
        The 'mov_lin' argument is a GVect or a GAxis instance. 
        Depending on that, the movement sense will be known
        when a GVect provided, and unknown/uncertain when a GAxis provided.

        Example:
          >>> Slickenline(GVect(90, 10))
          Slickenline(090.00, +10.00, True)
          >>> Slickenline(GAxis(90, 10))
          Slickenline(090.00, +10.00, False)
          >>> Slickenline(GPlane(180,20))
          Traceback (most recent call last):
          ...
          SlickelineInputTypeException: Input movement is not of the correct type
        """

        if not isinstance(mov_lin, (GVect, GAxis)):
            raise SlickelineInputTypeException("Input movement is not of the correct type")

        self._mov_lin = mov_lin

    @classmethod
    def from_trpl(cls, trend: float, plunge: float, known=True):
        """
        Class constructors from trend, plunge and optional known movement sense flag.

        :param trend: the trend of the slickenline
        :type trend: float
        :param plunge: the slickenline plunge
        :type plunge: float
        :param known: the known movement sense flag
        :type known: bool
        :return: the Slickenline instance

        Examples:
          >>> Slickenline.from_trpl(90, 10)
          Slickenline(090.00, +10.00, True)
          >>> Slickenline.from_trpl(90, 10, False)
          Slickenline(090.00, +10.00, False)
          >>> Slickenline.from_trpl("90", 10, False)
          Traceback (most recent call last):
          ...
          SlickelineInputTypeException: Trend must be a number
        """

        if not isinstance(trend, (int, float)):
            raise SlickelineInputTypeException("Trend must be a number")
        if not isinstance(plunge, (int, float)):
            raise SlickelineInputTypeException("Plunge must be a number")
        if not isinstance(known, bool):
            raise SlickelineInputTypeException("Known movement sense must be a boolean")

        val = GVect(trend, plunge) if known else GAxis(trend, plunge)

        return cls(val)

    def has_known_sense(self):
        """
        Check whether the slickenline has known movement sense.

        Example:
          >>> Slickenline(GVect(90, 45)).has_known_sense()
          True
          >>> Slickenline(GAxis(90, 45)).has_known_sense()
          False
        """

        if isinstance(self._mov_lin, GAxis):
            return False
        elif isinstance(self._mov_lin, GVect):
            return True
        else:
            raise SlickelineInputTypeException("Error with provided slickeline type")

    def has_unknown_sense(self):
        """
        Check whether the slickenline has unknown/uncertain movement sense.

        Example:
          >>> Slickenline(GAxis(90, 45)).has_unknown_sense()
          True
          >>> Slickenline(GVect(90, 45)).has_unknown_sense()
          False
        """

        return not self.has_known_sense()

    def set_known_sense(self):
        """
        Set (formal) movement sense to Slickline instance without known/certain movement sense.

        Example:
          >>> Slickenline(GAxis(180, -30)).set_known_sense() 
          Slickenline(180.00, -30.00, True)
        """

        return Slickenline(self.geom.as_gvect())

    def set_unknown_sense(self):
        """
        Set to unknown/uncertain the movement sense for the current Slickline instance. 

        Example:
          >>> Slickenline(GVect(180, -30)).set_unknown_sense() 
          Slickenline(180.00, -30.00, False)
        """

        return Slickenline(self.geom.as_axis())

    @property
    def geom(self):
        """
        Return the slickenline orientation value,
        as a GVect (known movement sense)
        or a GAxis instance (unknown movement sense).

        Example:
          >>> Slickenline(GVect(90, 45)).geom
          GVect(090.00, +45.00)
          >>> Slickenline(GAxis(90, 45)).geom
          GAxis(090.00, +45.00)
        """

        return self._mov_lin

    @property
    def vals(self):
        """
        The slickenline parameters.
        """

        known_mov = self.has_known_sense()

        return self._mov_lin.tr, self._mov_lin.pl, known_mov

    def __repr__(self):

        return "Slickenline({:06.2f}, {:+06.2f}, {})".format(*self.vals)

    def invert(self):
        """
        Invert the slickenline sense, when known, otherwise raise SlickelineSenseException.

        Example:
         >>> Slickenline(GAxis(30, 45)).invert()
         Traceback (most recent call last):
         ...
         SlickelineSenseException: Slickenline must have know movement sense
         >>> Slickenline(GVect(30, 45)).invert()
         Slickenline(210.00, -45.00, True)
        """

        if not self.has_known_sense():
            raise SlickelineSenseException("Slickenline must have know movement sense")

        return Slickenline(self.geom.opposite())


class FaultSlick(object):
    """
    Represent a couple of geological observations,
    made up by a fault plane, represented by a GPlane instance,
    and a slickenline observation, represented by a Slickenline instance.
    """

    def __init__(self, fault_plane, slickenline):
        """
        Create an instance of a FaultSlick.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45)))
          FaultSlick(GPlane(090.00, +45.00), Slickenline(090.00, +45.00, False))
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(80, 45)))
          Traceback (most recent call last):
          ...
          FaultSlickInputTypeException: Slickenline is not within fault plane
        """

        if not isinstance(fault_plane, GPlane):
            raise FaultSlickInputTypeException("Provided fault plane must be a GPlane instance")
        elif not isinstance(slickenline, Slickenline):
            raise FaultSlickInputTypeException("Provided slickenline must be a Slickenline instance")
        elif not are_close(fault_plane.normal().angle(slickenline.geom), 90.):
            raise FaultSlickInputTypeException("Slickenline is not within fault plane")

        self._fltpln = fault_plane
        self._slick = slickenline

    @property
    def gplane(self):
        """
        Return fault plane, as a GPlane instance.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).gplane
          GPlane(090.00, +45.00)
        """

        return self._fltpln

    @property
    def slick(self):
        """
        Return the slickenline associated with the fault. 

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).slick
          Slickenline(090.00, +45.00, False)
        """

        return self._slick

    def slick_geom(self):
        """
        Return the geometric object (GVect or GAxis) associated with slickenline.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).slick_geom()
          GAxis(090.00, +45.00)
        """

        return self.slick.geom

    @property
    def known_sense(self):
        """
        Check if the Slickenline instance in the FaultSlick instance has a known movement sense.

        Example: 
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).known_sense
          False
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).known_sense
          True
        """

        return self.slick.has_known_sense()

    def set_known_sense(self):
        """
        Create FaultSlick instance with known movement sense from another instance.

        Example:
          >>> FaultSlick(GPlane(0, 45), Slickenline(GAxis(0, 45))).set_known_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, True))
        """

        return FaultSlick(self.gplane, self.slick.set_known_sense())

    def set_unknown_sense(self):
        """
        Create FaultSlick instance with unknown/uncertain movement sense.

        Example:
          >>> FaultSlick(GPlane(0, 45), Slickenline(GVect(0, 45))).set_unknown_sense()
          FaultSlick(GPlane(000.00, +45.00), Slickenline(000.00, +45.00, False))
        """

        return FaultSlick(self.gplane, self.slick.set_unknown_sense())

    def __repr__(self):

        return "FaultSlick({}, {})".format(self.gplane, self.slick)

    def opposite_mov(self):
        """
        Create FaultSlick instance with opposite movement, when the source instance
        has defined movement sense, otherwise raise SlickelineSenseException.

        Example:
          >>> FaultSlick(GPlane(90, 45), Slickenline(GAxis(90, 45))).opposite_mov()
          Traceback (most recent call last):
          ...
          SlickelineSenseException: Fault slickenline must have known movement sense
          >>> FaultSlick(GPlane(90, 45), Slickenline(GVect(90, 45))).opposite_mov()
          FaultSlick(GPlane(090.00, +45.00), Slickenline(270.00, -45.00, True))
        """

        if not self.known_sense:
            raise SlickelineSenseException("Fault slickenline must have known movement sense")

        return FaultSlick(self.gplane, self.slick.invert())

    def rake(self):
        """
        Calculates the rake (sensu Aki & Richards, 1980) of the slickenline.
        The slickenlines must have known sense movement.

        :return: the rake value
        :rtype: double

        Examples:
          >>> fault_plane = GPlane(180, 45)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 0)).rake()
          0.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, -45)).rake()
          90.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, 0)).rake()
          180.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 45)).rake()
          -90.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 45, False)).rake()
          Traceback (most recent call last):
          ...
          FaultSlickInputTypeException: Slickeline must have known movement sense
          >>> fault_plane = GPlane(90, 90)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).rake()
          0.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 90)).rake()
          -90.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, -90)).rake()
          90.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 1)).rake()
          -179.0000000000001
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, -1)).rake()
          179.0000000000001
          >>> fault_plane = GPlane(90, 0)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).rake()
          0.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 0)).rake()
          -90.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(135, 0)).rake()
          -135.0
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, 0)).rake()
          90.00000000000001
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(225, 0)).rake()
          135.00000000000003
        """

        if not self.known_sense:
            raise FaultSlickInputTypeException("Slickeline must have known movement sense")

        sl_gv = self.slick_geom()
        angle = sl_gv.angle(self.gplane.strk_rhr_gv())

        if self.gplane.dipdir_gv().angle(sl_gv) < 90.0:
            return -angle
        else:
            return angle

    def is_normal(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has normal movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if normal, False if not applicable

        Examples:
          >>> fault_plane = GPlane(90, 90)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_normal()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 90)).is_normal()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_normal()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 45)).is_normal()
          True
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, -45)).is_normal()
          False
        """

        if self.gplane.is_very_high_angle(dip_angle_threshold) or self.gplane.is_very_low_angle(dip_angle_threshold):
            return False

        if - rk_threshold >= self.rake() >= -(180.0 - rk_threshold):
            return True
        else:
            return False

    def is_reverse(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has reverse movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if reverse, False if not applicable

        Examples:
          >>> fault_plane = GPlane(90, 90)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_reverse()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 90)).is_reverse()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_reverse()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, -45)).is_reverse()
          True
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(90, 45)).is_reverse()
          False
        """

        if self.gplane.is_very_high_angle(dip_angle_threshold) or self.gplane.is_very_low_angle(dip_angle_threshold):
            return False

        if rk_threshold <= self.rake() <= (180.0 - rk_threshold):
            return True
        else:
            return False

    def is_right_lateral(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has right-lateral movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if right-lateral, False if not applicable

        Examples:
          >>> fault_plane = GPlane(90, 90)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_right_lateral()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 0)).is_right_lateral()
          True
          >>> fault_plane = GPlane(90, 45)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_right_lateral()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 0)).is_right_lateral()
          True
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, -45)).is_right_lateral()
          False
          >>> fault_plane = GPlane(90, 2)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 0)).is_right_lateral()
          False
        """

        if self.gplane.is_very_low_angle(dip_angle_threshold):
            return False

        rake = self.rake()
        if rake >= (90.0 + rk_threshold) or rake <= (-90.0 - rk_threshold):
            return True
        else:
            return False

    def is_left_lateral(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has left-lateral movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if left-lateral, False if not applicable

        Examples:
          >>> fault_plane = GPlane(90, 90)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_left_lateral()
          True
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 0)).is_left_lateral()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_left_lateral()
          True
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(180, 0)).is_left_lateral()
          False
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(270, -45)).is_left_lateral()
          False
          >>> fault_plane = GPlane(90, 2)
          >>> FaultSlick(fault_plane, Slickenline.from_trpl(0, 0)).is_left_lateral()
          False
        """

        if self.gplane.is_very_low_angle(dip_angle_threshold):
            return False

        if (-90.0 + rk_threshold) <= self.rake() <= (90.0 - rk_threshold):
            return True
        else:
            return False


class SlickelineSenseException(Exception):
    """
    Exception for slickenline movement sense.
    """

    pass


class SlickelineInputTypeException(Exception):
    """
    Exception for slickenline type.
    """

    pass


class FaultSlickInputTypeException(Exception):
    """
    Exception for FaultSlick input type.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
