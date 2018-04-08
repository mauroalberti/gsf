
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
        :return: the Slick instance

        Examples:
          >>> Slick.from_trpl(90, 10)
          Slick(090.00, +10.00, True)
          >>> Slick.from_trpl(90, 10, False)
          Slick(090.00, +10.00, False)
          >>> Slick.from_trpl("90", 10, False)
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

        val = GVect(trend, plunge) if known else GAxis(trend, plunge)

        return cls(val)

    def has_known_sense(self):
        """
        Check whether the slickenline has known movement sense.

        Example:
          >>> Slick(GVect(90, 45)).has_known_sense()
          True
          >>> Slick(GAxis(90, 45)).has_known_sense()
          False
        """

        if isinstance(self._mov_lin, GAxis):
            return False
        elif isinstance(self._mov_lin, GVect):
            return True
        else:
            raise SlickInputTypeException("Error with provided slickeline type")

    def has_unknown_sense(self):
        """
        Check whether the slickenline has unknown/uncertain movement sense.

        Example:
          >>> Slick(GAxis(90, 45)).has_unknown_sense()
          True
          >>> Slick(GVect(90, 45)).has_unknown_sense()
          False
        """

        return not self.has_known_sense()

    def set_known_sense(self):
        """
        Set (formal) movement sense to Slickline instance without known/certain movement sense.

        Example:
          >>> Slick(GAxis(180, -30)).set_known_sense()
          Slick(180.00, -30.00, True)
        """

        return Slick(self.geom.as_gvect())

    def set_unknown_sense(self):
        """
        Set to unknown/uncertain the movement sense for the current Slickline instance.

        Example:
          >>> Slick(GVect(180, -30)).set_unknown_sense()
          Slick(180.00, -30.00, False)
        """

        return Slick(self.geom.as_gaxis())

    @property
    def geom(self):
        """
        Return the slickenline orientation value,
        as a GVect (known movement sense)
        or a GAxis instance (unknown movement sense).

        Example:
          >>> Slick(GVect(90, 45)).geom
          GVect(090.00, +45.00)
          >>> Slick(GAxis(90, 45)).geom
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



    def invert(self):
        """
        Invert the slickenline sense, when known, otherwise raise SlickSenseException.

        Example:
         >>> Slick(GAxis(30, 45)).invert()
         Traceback (most recent call last):
         ...
         SlickSenseException: Slick must have know movement sense
         >>> Slick(GVect(30, 45)).invert()
         Slick(210.00, -45.00, True)
        """

        if not self.has_known_sense():
            raise SlickSenseException("Slick must have know movement sense")

        return Slick(self.geom.opposite())


class GFault(object):
    """
    Represent a couple of geological observations,
    made up by a fault plane, represented by a GPlane instance,
    and a slickenline observation, represented by a Slick instance.
    """

    def __init__(self, fault_plane, slickenline):
        """
        Create an instance of a GFault.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45)))
          GFault(GPlane(090.00, +45.00), Slick(090.00, +45.00, False))
          >>> GFault(GPlane(90, 45), Slick(GAxis(80, 45)))
          Traceback (most recent call last):
          ...
          GFaultInputTypeException: Slick is not within fault plane
        """

        if not isinstance(fault_plane, GPlane):
            raise GFaultInputTypeException("Provided fault plane must be a GPlane instance")
        elif not isinstance(slickenline, Slick):
            raise GFaultInputTypeException("Provided slickenline must be a Slick instance")
        elif not are_close(fault_plane._normal_gv_frwrd().angle(slickenline.geom), 90.):
            raise GFaultInputTypeException("Slick is not within fault plane")

        self._fltpln = fault_plane
        self._slick = slickenline

    @property
    def gplane(self):
        """
        Return fault plane, as a GPlane instance.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45))).gplane
          GPlane(090.00, +45.00)
        """

        return self._fltpln

    @property
    def slick(self):
        """
        Return the slickenline associated with the fault.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45))).slick
          Slick(090.00, +45.00, False)
        """

        return self._slick

    def slick_geom(self):
        """
        Return the geometric object (GVect or GAxis) associated with slickenline.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45))).slick_geom()
          GAxis(090.00, +45.00)
        """

        return self.slick.geom

    @property
    def known_sense(self):
        """
        Check if the Slick instance in the GFault instance has a known movement sense.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45))).known_sense
          False
          >>> GFault(GPlane(90, 45), Slick(GVect(90, 45))).known_sense
          True
        """

        return self.slick.has_known_sense()

    def set_known_sense(self):
        """
        Create GFault instance with known movement sense from another instance.

        Example:
          >>> GFault(GPlane(0, 45), Slick(GAxis(0, 45))).set_known_sense()
          GFault(GPlane(000.00, +45.00), Slick(000.00, +45.00, True))
        """

        return GFault(self.gplane, self.slick.set_known_sense())

    def set_unknown_sense(self):
        """
        Create GFault instance with unknown/uncertain movement sense.

        Example:
          >>> GFault(GPlane(0, 45), Slick(GVect(0, 45))).set_unknown_sense()
          GFault(GPlane(000.00, +45.00), Slick(000.00, +45.00, False))
        """

        return GFault(self.gplane, self.slick.set_unknown_sense())

    def __repr__(self):

        return "GFault({}, {})".format(self.gplane, self.slick)

    def opposite_mov(self):
        """
        Create GFault instance with opposite movement, when the source instance
        has defined movement sense, otherwise raise SlickSenseException.

        Example:
          >>> GFault(GPlane(90, 45), Slick(GAxis(90, 45))).opposite_mov()
          Traceback (most recent call last):
          ...
          SlickSenseException: Fault slickenline must have known movement sense
          >>> GFault(GPlane(90, 45), Slick(GVect(90, 45))).opposite_mov()
          GFault(GPlane(090.00, +45.00), Slick(270.00, -45.00, True))
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
          >>> fault_plane = GPlane(180, 45)
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
          >>> fault_plane = GPlane(90, 90)
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
          >>> fault_plane = GPlane(90, 0)
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
        angle = sl_gv.angle(self.gplane.strk_rhr_gv())

        if self.gplane.dipdir_gv().angle(sl_gv) < 90.0:
            return -angle
        else:
            return angle

    def is_normal(self, rk_threshold=rake_threshold, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a fault has _normal_gv_frwrd movements.

        :param rk_threshold: the threshold, in degrees, for the rake angle
        :param dip_angle_threshold: the threshold, in degrees, for the dip angle of the geological plane
        :return: True if _normal_gv_frwrd, False if not applicable

        Examples:
          >>> fault_plane = GPlane(90, 90)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_normal()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 90)).is_normal()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_normal()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 45)).is_normal()
          True
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_normal()
          False
        """

        if self.gplane.is_vhigh_angle(dip_angle_threshold) or self.gplane.is_vlow_angle(dip_angle_threshold):
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
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_reverse()
          False
          >>> GFault(fault_plane, Slick.from_trpl(90, 90)).is_reverse()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_reverse()
          False
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_reverse()
          True
          >>> GFault(fault_plane, Slick.from_trpl(90, 45)).is_reverse()
          False
        """

        if self.gplane.is_vhigh_angle(dip_angle_threshold) or self.gplane.is_vlow_angle(dip_angle_threshold):
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
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_right_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          True
          >>> fault_plane = GPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_right_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_right_lateral()
          False
          >>> fault_plane = GPlane(90, 2)
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_right_lateral()
          False
        """

        if self.gplane.is_vlow_angle(dip_angle_threshold):
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
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_left_lateral()
          False
          >>> fault_plane = GPlane(90, 45)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          True
          >>> GFault(fault_plane, Slick.from_trpl(180, 0)).is_left_lateral()
          False
          >>> GFault(fault_plane, Slick.from_trpl(270, -45)).is_left_lateral()
          False
          >>> fault_plane = GPlane(90, 2)
          >>> GFault(fault_plane, Slick.from_trpl(0, 0)).is_left_lateral()
          False
        """

        if self.gplane.is_vlow_angle(dip_angle_threshold):
            return False

        if (-90.0 + rk_threshold) <= self.rake() <= (90.0 - rk_threshold):
            return True
        else:
            return False


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