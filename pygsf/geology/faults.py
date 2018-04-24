# -*- coding: utf-8 -*-


from ..orientations.orientations import *
from ..exceptions.geology import *


class Slick(Direct):
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


class FaultMS(PPlane):
    """
    Represent a fault plane, represented by a PPlane instance,
    and zero, one or more slickenlines, represented by a list of Slick instances (None when no slickenlines).
    """

    def __init__(self, azim: float, dip_ang: float, is_rhr_strike=False, slickenlines=None):
        """
        Create an instance of a FaultMS.

        :param  azim:  azimuth of the plane (RHR strike or dip direction).
        :type  azim:  number or string convertible to float.
        :param  dip_ang:  Dip angle of the plane (0-90°).
        :type  dip_ang:  number or string convertible to float.
        :param is_rhr_strike: if the source azimuth is RHR strike (default is False, i.e. it is dip direction)
        :return: the instantiated geological plane.

        Example:
          >>> FaultMS(90, 45, [Slick(90, 45)])
          FaultMS(090.00, +45.00) with 1 slickenline(s))
          >>> FaultMS(90, 45. [Slick(80, 45)])
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

        return "FaultMS({}, {}) with {} slickenine(s)".format(*self.dda, len(self.slick))


class GFault(object):
    """
    Represent a couple of geological observations,
    made up by a fault plane, represented by a PPlane instance,
    and a slickenline observation, represented by a Slick instance.
    """

    def __init__(self, fault_plane, slickenline):
        """
        Create an instance of a GFault.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45)))
          GFault(PPlane(090.00, +45.00), Slick(090.00, +45.00, False))
          >>> GFault(PPlane(90, 45), Slick(GAxis(80, 45)))
          Traceback (most recent call last):
          ...
          GFaultInputTypeException: Slick is not within fault plane
        """

        if not isinstance(fault_plane, PPlane):
            raise GFaultInputTypeException("Provided fault plane must be a PPlane instance")
        elif not isinstance(slickenline, Slick):
            raise GFaultInputTypeException("Provided slickenline must be a Slick instance")
        elif not are_close(fault_plane._normal_orien_frwrd().angle(slickenline.geom), 90.):
            raise GFaultInputTypeException("Slick is not within fault plane")

        self._fltpln = fault_plane
        self._slick = slickenline

    @property
    def gplane(self):
        """
        Return fault plane, as a PPlane instance.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45))).gplane
          PPlane(090.00, +45.00)
        """

        return self._fltpln

    @property
    def slick(self):
        """
        Return the slickenline associated with the fault.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45))).slick
          Slick(090.00, +45.00, False)
        """

        return self._slick

    def slick_geom(self):
        """
        Return the geometric object (GVect or GAxis) associated with slickenline.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45))).slick_geom()
          GAxis(090.00, +45.00)
        """

        return self.slick.geom

    @property
    def known_sense(self):
        """
        Check if the Slick instance in the GFault instance has a known movement sense.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45))).known_sense
          False
          >>> GFault(PPlane(90, 45), Slick(GVect(90, 45))).known_sense
          True
        """

        return self.slick.has_known_sense()

    def set_known_sense(self):
        """
        Create GFault instance with known movement sense from another instance.

        Example:
          >>> GFault(PPlane(0, 45), Slick(GAxis(0, 45))).set_known_sense()
          GFault(PPlane(000.00, +45.00), Slick(000.00, +45.00, True))
        """

        return GFault(self.gplane, self.slick.set_known_sense())

    def set_unknown_sense(self):
        """
        Create GFault instance with unknown/uncertain movement sense.

        Example:
          >>> GFault(PPlane(0, 45), Slick(GVect(0, 45))).set_unknown_sense()
          GFault(PPlane(000.00, +45.00), Slick(000.00, +45.00, False))
        """

        return GFault(self.gplane, self.slick.set_unknown_sense())

    def __repr__(self):

        return "GFault({}, {})".format(self.gplane, self.slick)

    def opposite_mov(self):
        """
        Create GFault instance with opposite movement, when the source instance
        has defined movement sense, otherwise raise SlickSenseException.

        Example:
          >>> GFault(PPlane(90, 45), Slick(GAxis(90, 45))).opposite_mov()
          Traceback (most recent call last):
          ...
          SlickSenseException: Fault slickenline must have known movement sense
          >>> GFault(PPlane(90, 45), Slick(GVect(90, 45))).opposite_mov()
          GFault(PPlane(090.00, +45.00), Slick(270.00, -45.00, True))
        """

        if not self.known_sense:
            raise SlickSenseException("Fault slickenline must have known movement sense")

        return GFault(self.gplane, self.slick.invert())

    def rake(self):
        """
        Calculates the rake (sensu Aki & Richards, 1980) of the slickenline.
        The slickenlines must have known sense movement.

        :return: the rake value
        :rtype: double

        Examples:
          >>> fault_plane = PPlane(180, 45)
          >>> GFault(fault_plane, Slick.from_trpl(90, 0)).rake()
          0.0
          >>> GFault(fault_plane, Slick.from_trpl(0, -45)).rake()
          90.0
          >>> GFault(fault_plane, Slick.from_trpl(270, 0)).rake()
          180.0
          >>> GFault(fault_plane, Slick.from_trpl(180, 45)).rake()
          -90.0
          >>> GFault(fault_plane, Slick.from_trpl(180, 45, False)).rake()
          Traceback (most recent call last):
          ...
          GFaultInputTypeException: Slickeline must have known movement sense
          >>> fault_plane = PPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).rake()
          0.0
          >>> GFault(fault_plane, Slick.from_trpl(90, 90)).rake()
          -90.0
          >>> GFault(fault_plane, Slick.from_trpl(90, -90)).rake()
          90.0
          >>> GFault(fault_plane, Slick.from_trpl(180, 1)).rake()
          -179.0000000000001
          >>> GFault(fault_plane, Slick.from_trpl(180, -1)).rake()
          179.0000000000001
          >>> fault_plane = PPlane(90, 0)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).rake()
          0.0
          >>> GFault(fault_plane, Slick.from_trpl(90, 0)).rake()
          -90.0
          >>> GFault(fault_plane, Slick.from_trpl(135, 0)).rake()
          -135.0
          >>> GFault(fault_plane, Slick.from_trpl(270, 0)).rake()
          90.00000000000001
          >>> GFault(fault_plane, Slick.from_trpl(225, 0)).rake()
          135.00000000000003
        """

        if not self.known_sense:
            raise GFaultInputTypeException("Slickeline must have known movement sense")

        sl_gv = self.slick_geom()
        angle = sl_gv.angle(self.gplane.rhrStrikeOrien())

        if self.gplane.dipDirOrien().angle(sl_gv) < 90.0:
            return -angle
        else:
            return angle

    def is_normal(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has _normal_orien_frwrd movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if _normal_orien_frwrd, False if not applicable

        Examples:
          >>> fault_plane = PPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_normal()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 90)).is_normal()
          False
          >>> fault_plane = PPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_normal()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 45)).is_normal()
          True
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_normal()
          False
        """

        if self.gplane.isVHighAngle(dip_angle_threshold) or self.gplane.isVLowAngle(dip_angle_threshold):
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
          >>> fault_plane = PPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_reverse()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 90)).is_reverse()
          False
          >>> fault_plane = PPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_reverse()
          False
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_reverse()
          True
          >>> GFault(fault_plane, Slick.from_trpl(90, 45)).is_reverse()
          False
        """

        if self.gplane.isVHighAngle(dip_angle_threshold) or self.gplane.isVLowAngle(dip_angle_threshold):
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
          >>> fault_plane = PPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_right_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          True
          >>> fault_plane = PPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_right_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_right_lateral()
          False
          >>> fault_plane = PPlane(90, 2)
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          False
        """

        if self.gplane.isVLowAngle(dip_angle_threshold):
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
          >>> fault_plane = PPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_left_lateral()
          False
          >>> fault_plane = PPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_left_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_left_lateral()
          False
          >>> fault_plane = PPlane(90, 2)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          False
        """

        if self.gplane.isVLowAngle(dip_angle_threshold):
            return False

        if (-90.0 + rk_threshold) <= self.rake() <= (90.0 - rk_threshold):
            return True
        else:
            return False




if __name__ == "__main__":

    import doctest
    doctest.testmod()
