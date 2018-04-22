


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

        if not isinstance(fault_plane, GPlane):
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



