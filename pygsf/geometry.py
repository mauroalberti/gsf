# -*- coding: utf-8 -*-


from __future__ import division
#from __future__ import annotations

from math import sqrt, sin, cos, radians, acos, atan, isnan, pi
import numpy as np

from typing import Dict, Tuple, List

from .default_parameters import *
from .arrays import *
from .mathematics import *
from .geometry_utils import *


isfinite = np.isfinite
array = np.array


class Azim(object):
    """
    Azim class
    """

    def __init__(self, val: [int, float], unit: str='d'):
        """
        Creates an azimuth instance.

        :param val: azimuth value
        :param unit: angle measurement unit, 'd' (default, stands for decimal degrees) or 'r' (stands for radians)

        Examples:
          >>> Azim(10)
          Azimuth(10.00°)
          >>> Azim(370)
          Azimuth(10.00°)
          >>> Azim(pi/2, unit='r')
          Azimuth(90.00°)
          >>> Azim("10")
          Traceback (most recent call last):
          ...
          GeomInputException: Input azimuth value must be int/float
          >>> Azim(np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Input azimuth value must be finite
        """

        # unit check
        if unit not in ("d", "r"):
            raise GeomInputException("Unit input must be 'd' or 'r'")

        if not (isinstance(val, (int, float))):
            raise GeomInputException("Input azimuth value must be int/float")
        elif not isfinite(val):
            raise GeomInputException("Input azimuth value must be finite")

        if unit == 'd':
            val = radians(val)

        self.a = val % (2*pi)

    @property
    def d(self):
        """
        Returns the angle in decimal degrees.

        :return: angle in decimal degrees

        Example:
          >>> Azim(10).d
          10.0
          >>> Azim(pi/2, unit='r').d
          90.0
        """

        return degrees(self.a)

    @property
    def r(self):
        """
        Returns the angle in radians.

        :return: angle in radians

        Example:
          >>> Azim(180).r
          3.141592653589793
        """

        return self.a

    @classmethod
    def fromXY(cls, x: [int, float], y: [int, float]) -> 'Azim':
        """
        Calculates azimuth given cartesian components.

        :param cls: class
        :param x: x component
        :param y: y component
        :return: Azimuth instance

        Examples:
          >>> Azim.fromXY(1, 1)
          Azimuth(45.00°)
          >>> Azim.fromXY(1, -1)
          Azimuth(135.00°)
          >>> Azim.fromXY(-1, -1)
          Azimuth(225.00°)
          >>> Azim.fromXY(-1, 1)
          Azimuth(315.00°)
          >>> Azim.fromXY(0, 0)
          Azimuth(0.00°)
          >>> Azim.fromXY(0, np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Input x and y values must be finite
          >>> Azim.fromXY("10", np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Input x and y values must be integer or float
        """

        # input vals checks
        vals = [x, y]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            raise GeomInputException("Input x and y values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise GeomInputException("Input x and y values must be finite")

        angle = atan2(x, y)
        return cls(angle, unit='r')

    def __repr__(self) -> str:

        return "Azimuth({:.2f}°)".format(self.d)

    def toXY(self) -> Tuple[float, float]:
        """
        Converts an azimuth to x-y components.

        :return: a tuple storing x and y values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Azim(0).toXY())
          (0.0, 1.0)
          >>> apprFTuple(Azim(90).toXY())
          (1.0, 0.0)
          >>> apprFTuple(Azim(180).toXY())
          (0.0, -1.0)
          >>> apprFTuple(Azim(270).toXY())
          (-1.0, 0.0)
          >>> apprFTuple(Azim(360).toXY())
          (0.0, 1.0)
        """

        return sin(self.a), cos(self.a)


class Plunge(object):
    """
    Class representing a plunge
    """

    def __init__(self, val: [int, float], unit: str='d'):
        """
        Creates a Plunge instance.

        :param val: plunge value
        :param unit: angle measurement unit, decimal degrees ('d') or radians ('r')

        Examples:
          >>> Plunge(10)
          Plunge(10.00°)
          >>> Plunge("10")
          Traceback (most recent call last):
          ...
          GeomInputException: Input plunge value must be int/float
          >>> Plunge(np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Input plunge value must be finite
          >>> Plunge(-100)
          Traceback (most recent call last):
          ...
          GeomInputException: Input value in degrees must be between -90° and 90°
         """

        # unit check
        if unit not in ('d', 'r'):
            raise GeomInputException("Unit input must be 'd' (for degrees) or 'r' (for radians)")

        # val check
        if not (isinstance(val, (int, float))):
            raise GeomInputException("Input plunge value must be int/float")
        elif not isfinite(val):
            raise GeomInputException("Input plunge value must be finite")
        if unit == 'd' and not (-90.0 <= val <= 90.0):
            raise GeomInputException("Input value in degrees must be between -90° and 90°")
        elif unit == 'r' and not (-pi/2 <= val <= pi/2):
            raise GeomInputException("Input value in radians must be between -pi/2 and pi/2")

        if unit == 'd':
            val = radians(val)

        self.p = val

    @property
    def d(self):
        """
        Returns the angle in decimal degrees.

        :return: angle in decimal degrees

        Example:
          >>> Plunge(10).d
          10.0
          >>> Plunge(-pi/2, unit='r').d
          -90.0
        """

        return degrees(self.p)

    @property
    def r(self):
        """
        Returns the angle in radians.

        :return: angle in radians

        Example:
          >>> Plunge(90).r
          1.5707963267948966
          >>> Plunge(45).r
          0.7853981633974483
        """

        return self.p

    @classmethod
    def fromHZ(cls, h: [int, float], z: [int, float]) -> 'Plunge':
        """
        Calculates plunge from h and z components.

        :param cls: class
        :param h: horizontal component (always positive)
        :param z: vertical component (positive upward)
        :return: Plunge instance

        Examples:
          >>> Plunge.fromHZ(1, 1)
          Plunge(-45.00°)
          >>> Plunge.fromHZ(1, -1)
          Plunge(45.00°)
          >>> Plunge.fromHZ(0, 1)
          Plunge(-90.00°)
          >>> Plunge.fromHZ(0, -1)
          Plunge(90.00°)
          >>> Plunge.fromHZ(-1, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Horizontal component cannot be negative
          >>> Plunge.fromHZ(0, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Input h and z values cannot be both zero
        """

        # input vals check

        vals = [h, z]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            raise GeomInputException("Input h and z values must be integer or float")
        elif not all(map(isfinite, vals)):
            raise GeomInputException("Input h and z values must be finite")

        if h == 0.0 and z == 0.0:
            raise GeomInputException("Input h and z values cannot be both zero")
        elif h < 0.0:
            raise GeomInputException("Horizontal component cannot be negative")

        angle = atan2(-z, h)

        return cls(angle, unit='r')

    def __repr__(self) -> str:

        return "Plunge({:.2f}°)".format(self.d)

    def toHZ(self):

        """
        Converts an azimuth to h-z components.

        :return: a tuple storing h (horizontal) and z values:
        :type: tuple of two floats

        Examples:
          >>> apprFTuple(Plunge(0).toHZ())
          (1.0, 0.0)
          >>> apprFTuple(Plunge(90).toHZ())
          (0.0, -1.0)
          >>> apprFTuple(Plunge(-90).toHZ())
          (0.0, 1.0)
          >>> apprFTuple(Plunge(-45).toHZ(), ndec=6)
          (0.707107, 0.707107)
          >>> apprFTuple(Plunge(45).toHZ(), ndec=6)
          (0.707107, -0.707107)
        """

        return cos(self.p), -sin(self.p)

    @property
    def isUpward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Plunge(10).isUpward
          False
          >>> Plunge(0.0).isUpward
          False
          >>> Plunge(-45).isUpward
          True
        """

        return self.r < 0.0

    @property
    def isDownward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Plunge(15).isDownward
          True
          >>> Plunge(0.0).isDownward
          False
          >>> Plunge(-45).isDownward
          False
        """

        return self.r > 0.0


class Orien(object):
    """
    Class describing an orientation, expressed as a polar direction.
    """

    def __init__(self, az: 'Azim', pl: 'Plunge'):
        """
        Creates a polar direction instance.

        :param az: the azimuth value
        :param pl: the plunge value
        """

        if not isinstance(az, Azim):
            raise GeomInputException("First input value must be of type Azim")

        if not isinstance(pl, Plunge):
            raise GeomInputException("Second input value must be of type Plunge")

        self._az = az
        self._pl = pl

    @property
    def d(self):
        """
        Returns azimuth and plunge in decimal degrees as a tuple.

        :return: tuple of azimuth and plunge in decimal degrees

        Example:
          >>> Orien.fromAzPl(100, 20).d
          (100.0, 20.0)
          >>> Orien.fromAzPl(-pi/2, -pi/4, unit='r').d
          (270.0, -45.0)
        """

        return self.az.d, self.pl.d

    @property
    def r(self):
        """
        Returns azimuth and plunge in radians as a tuple.

        :return: tuple of azimuth and plunge in radians

        Example:
          >>> Orien.fromAzPl(90, 45).r
          (1.5707963267948966, 0.7853981633974483)
        """

        return self.az.r, self.pl.r

    @property
    def az(self):
        """
        Returns the azimuth instance.

        :return: Azimuth
        """

        return self._az

    @property
    def pl(self):
        """
        Returns the plunge instance.

        :return: Plunge
        """

        return self._pl

    @classmethod
    def fromAzPl(cls, az: [int, float], pl: [int, float], unit='d') -> 'Orien':
        """
        Class constructor from trend and plunge.

        :param az: trend value
        :param pl: plunge value
        :param unit: measurement unit, in degrees ('d') or radians ('r')
        :return: Orientation instance

        Examples:
          >>> Orien.fromAzPl(30, 40)
          Orien(az: 30.00°, pl: 40.00°)
          >>> Orien.fromAzPl(370, 80)
          Orien(az: 10.00°, pl: 80.00°)
          >>> Orien.fromAzPl(pi/2, pi/4, unit='r')
          Orien(az: 90.00°, pl: 45.00°)
          >>> Orien.fromAzPl(280, -100)
          Traceback (most recent call last):
          ...
          GeomInputException: Input value in degrees must be between -90° and 90°
          >>> Orien.fromAzPl("10", 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Input azimuth value must be int/float
          >>> Orien.fromAzPl(100, np.nan)
          Traceback (most recent call last):
          ...
          GeomInputException: Input plunge value must be finite
        """

        azim = Azim(az, unit=unit)
        plng = Plunge(pl, unit=unit)

        return cls(azim, plng)

    @classmethod
    def _from_xyz(cls, x: [int, float], y: [int, float], z: [int, float]) -> 'Orien':
        """
        Private class constructor from three Cartesian values. Note: norm of components is unit.

        :param x: x component
        :param y: y component
        :param z: z component
        :return: Orientation instance
        """

        h = sqrt(x*x + y*y)

        az = Azim.fromXY(x, y)
        pl = Plunge.fromHZ(h, z)

        return cls(az, pl)

    @classmethod
    def fromXYZ(cls, x: [int, float], y: [int, float], z: [int, float]) -> 'Orien':
        """
        Class constructor from three generic Cartesian values.

        :param x: x component
        :param y: y component
        :param z: z component
        :return: Orientation instance

        Examples:
          >>> Orien.fromXYZ(1, 0, 0)
          Orien(az: 90.00°, pl: -0.00°)
          >>> Orien.fromXYZ(0, 1, 0)
          Orien(az: 0.00°, pl: -0.00°)
          >>> Orien.fromXYZ(0, 0, 1)
          Orien(az: 0.00°, pl: -90.00°)
          >>> Orien.fromXYZ(0, 0, -1)
          Orien(az: 0.00°, pl: 90.00°)
          >>> Orien.fromXYZ(1, 1, 0)
          Orien(az: 45.00°, pl: -0.00°)
          >>> Orien.fromXYZ(0.5, -0.5, -0.7071067811865476)
          Orien(az: 135.00°, pl: 45.00°)
          >>> Orien.fromXYZ(-0.5, 0.5, 0.7071067811865476)
          Orien(az: 315.00°, pl: -45.00°)
          >>> Orien.fromXYZ(0, 0, 0)
          Traceback (most recent call last):
          ...
          GeomInputException: Input components have near-zero values
        """

        mag, norm_xyz = normXYZ(x, y, z)

        if norm_xyz is None:
            raise GeomInputException("Input components have near-zero values")

        return Orien._from_xyz(*norm_xyz)

    def __repr__(self) -> str:

        return "Orien(az: {:.2f}°, pl: {:.2f}°)".format(*self.d)

    def toXYZ(self) -> Tuple[float, float, float]:
        """
        Converts an orientation to a tuple of x, y and z cartesian components (with unit norm).

        :return: tuple of x, y and z components.

        Examples:
          >>> az, pl = Azim(90), Plunge(0)
          >>> apprFTuple(Orien(az, pl).toXYZ())
          (1.0, 0.0, 0.0)
          >>> az, pl = Azim(135), Plunge(45)
          >>> apprFTuple(Orien(az, pl).toXYZ(), ndec=6)
          (0.5, -0.5, -0.707107)
          >>> az, pl = Azim(135), Plunge(0)
          >>> apprFTuple(Orien(az, pl).toXYZ(), ndec=6)
          (0.707107, -0.707107, 0.0)
          >>> az, pl = Azim(180), Plunge(45)
          >>> apprFTuple(Orien(az, pl).toXYZ(), ndec=6)
          (0.0, -0.707107, -0.707107)
          >>> az, pl = Azim(225), Plunge(-45)
          >>> apprFTuple(Orien(az, pl).toXYZ(), ndec=6)
          (-0.5, -0.5, 0.707107)
          >>> az, pl = Azim(270), Plunge(90)
          >>> apprFTuple(Orien(az, pl).toXYZ(), ndec=6)
          (0.0, 0.0, -1.0)
        """

        x, y = self.az.toXY()
        h, z = self.pl.toHZ()

        return x*h, y*h, z

    def copy(self):
        """
        Return a copy of the Orientation instance.

        Example:
          >>> Orien.fromAzPl(10, 20).copy()
          Orien(az: 10.00°, pl: 20.00°)
        """

        return Orien(self.az, self.pl)

    def opposite(self):
        """
        Return the opposite orientation.

        Example:
          >>> Orien.fromAzPl(0, 30).opposite()
          Orien(az: 180.00°, pl: -30.00°)
          >>> Orien.fromAzPl(315, 10).opposite()
          Orien(az: 135.00°, pl: -10.00°)
          >>> Orien.fromAzPl(135, 0).opposite()
          Orien(az: 315.00°, pl: -0.00°)
        """

        az, pl = self.r

        az = (az + pi) % (2*pi)
        pl = -pl

        return Orien.fromAzPl(az, pl, unit='r')

    def mirrorHoriz(self):
        """
        Return the mirror Orientation using a horizontal plane.

        Example:
          >>> Orien.fromAzPl(0, 30).mirrorHoriz()
          Orien(az: 0.00°, pl: -30.00°)
          >>> Orien.fromAzPl(315, 10).mirrorHoriz()
          Orien(az: 315.00°, pl: -10.00°)
          >>> Orien.fromAzPl(135, 0).mirrorHoriz()
          Orien(az: 135.00°, pl: -0.00°)
        """

        az = self.az.r
        pl = -self.pl.r

        return Orien.fromAzPl(az, pl, unit='r')

    @property
    def colatNorth(self) -> float:
        """
        Calculates the colatitude from the North (top).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: float

        Examples:
          >>> Orien.fromAzPl(320, 90).colatNorth
          180.0
          >>> Orien.fromAzPl(320, 45).colatNorth
          135.0
          >>> Orien.fromAzPl(320, 0).colatNorth
          90.0
          >>> Orien.fromAzPl(320, -45).colatNorth
          45.0
          >>> Orien.fromAzPl(320, -90).colatNorth
          0.0
        """

        return plng2colatTop(self.pl.d)

    @property
    def colatSouth(self) -> float:
        """
        Calculates the colatitude from the South (bottom).

        :return: an angle between 0 and 180 (in degrees).
        :rtype: float

        Examples:
          >>> Orien.fromAzPl(320, 90).colatSouth
          0.0
          >>> Orien.fromAzPl(320, 45).colatSouth
          45.0
          >>> Orien.fromAzPl(320, 0).colatSouth
          90.0
          >>> Orien.fromAzPl(320, -45).colatSouth
          135.0
          >>> Orien.fromAzPl(320, -90).colatSouth
          180.0
        """

        return plng2colatBottom(self.pl.d)

    def asVersor(self):
        """
        Return the unit vector corresponding to the Orien instance.

        Examples:
          >>> Orien.fromAzPl(0, 90).asVersor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> Orien.fromAzPl(0, -90).asVersor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> Orien.fromAzPl(90, 90).asVersor()
          Vect(0.0000, 0.0000, -1.0000)
        """

        az, pl = self.r
        cos_az, cos_pl = cos(az), cos(pl)
        sin_az, sin_pl = sin(az), sin(pl)
        north_coord = cos_pl * cos_az
        east_coord = cos_pl * sin_az
        down_coord = sin_pl

        return Vect(east_coord, north_coord, -down_coord)

    @property
    def isUpward(self):
        """
        Check whether the instance is pointing upward or horizontal.

        Examples:
          >>> Orien.fromAzPl(10, 15).isUpward
          False
          >>> Orien.fromAzPl(257.4, 0.0).isUpward
          False
          >>> Orien.fromAzPl(90, -45).isUpward
          True
        """

        return self.pl.isUpward

    @property
    def isDownward(self):
        """
        Check whether the instance is pointing downward or horizontal.

        Examples:
          >>> Orien.fromAzPl(10, 15).isDownward
          True
          >>> Orien.fromAzPl(257.4, 0.0).isDownward
          False
          >>> Orien.fromAzPl(90, -45).isDownward
          False
        """

        return self.pl.isDownward

    def upward(self):
        """
        Return upward-point geological vector.

        Examples:
          >>> Orien.fromAzPl(90, -45).upward().isAlmostParallel(Orien.fromAzPl(90.0, -45.0))
          True
          >>> Orien.fromAzPl(180, 45).upward().isAlmostParallel(Orien.fromAzPl(0.0, -45.0))
          True
          >>> Orien.fromAzPl(0, 0).upward().isAlmostParallel(Orien.fromAzPl(0.0, 0.0))
          True
          >>> Orien.fromAzPl(0, 90).upward().isAlmostParallel(Orien.fromAzPl(180.0, -90.0))
          True
          >>> Orien.fromAzPl(90, -45).upward().isAlmostParallel(Orien.fromAzPl(90.0, -35.0))
          False
          >>> Orien.fromAzPl(180, 45).upward().isAlmostParallel(Orien.fromAzPl(10.0, -45.0))
          False
          >>> Orien.fromAzPl(0, 0).upward().isAlmostParallel(Orien.fromAzPl(170.0, 0.0))
          False
          >>> Orien.fromAzPl(0, 90).upward().isAlmostParallel(Orien.fromAzPl(180.0, -80.0))
          False
        """

        if not self.isDownward:
            return self.copy()
        else:
            return self.opposite()

    def downward(self):
        """
        Return downward-pointing geological vector.

        Examples:
          >>> Orien.fromAzPl(90, -45).downward().isAlmostParallel(Orien.fromAzPl(270.0, 45.0))
          True
          >>> Orien.fromAzPl(180, 45).downward().isAlmostParallel(Orien.fromAzPl(180.0, 45.0))
          True
          >>> Orien.fromAzPl(0, 0).downward().isAlmostParallel(Orien.fromAzPl(180.0, 0.0))
          False
          >>> Orien.fromAzPl(0, 90).downward().isAlmostParallel(Orien.fromAzPl(0.0, 90.0))
          True
          >>> Orien.fromAzPl(90, -45).downward().isAlmostParallel(Orien.fromAzPl(270.0, 35.0))
          False
          >>> Orien.fromAzPl(180, 45).downward().isAlmostParallel(Orien.fromAzPl(170.0, 45.0))
          False
          >>> Orien.fromAzPl(0, 0).downward().isAlmostParallel(Orien.fromAzPl(180.0, 10.0))
          False
          >>> Orien.fromAzPl(0, 90).downward().isAlmostParallel(Orien.fromAzPl(0.0, 80.0))
          False
        """

        if not self.isUpward:
            return self.copy()
        else:
            return self.opposite()

    def isAbsDipInRange(self, min_val, max_val, min_val_incl=False, max_value_incl=True):
        """
        Check whether the absolute value of the dip angle of an Orien instance is within a given range
        (default: minimum value is not included, maximum value is included).

        :param min_val: the minimum dip angle, positive, domain: 0-90°.
        :param max_val: the maximum dip angle, positive, domain: 0-90°.
        :param min_val_incl: is minimum value included, boolean.
        :param max_value_incl: is maximum value included, boolean.
        :return: Boolean

        Examples:
          >>> Orien.fromAzPl(90, -45).isAbsDipInRange(30, 60)
          True
          >>> Orien.fromAzPl(120, 0).isAbsDipInRange(0, 60)
          False
          >>> Orien.fromAzPl(120, 0).isAbsDipInRange(0, 60, min_val_incl=True)
          True
          >>> Orien.fromAzPl(120, 60).isAbsDipInRange(0, 60)
          True
        """

        abs_dip = abs(self.pl.d)

        if abs_dip < min_val or abs_dip > max_val:
            return False
        elif abs_dip == min_val:
            if min_val_incl:
                return True
            else:
                return False
        elif abs_dip == max_val:
            if max_value_incl:
                return True
            else:
                return False
        else:
            return True

    def isSubHorizontal(self, max_dip_angle=DIP_ANGLE_THRESHOLD):
        """
        Check whether the instance is almost horizontal.

        Examples:
          >>> Orien.fromAzPl(10, 15).isSubHorizontal()
          False
          >>> Orien.fromAzPl(257, 2).isSubHorizontal()
          True
          >>> Orien.fromAzPl(90, -5).isSubHorizontal()
          False
        """

        return abs(self.pl.d) < max_dip_angle

    def isSubVertical(self, min_dip_angle=90.0 - DIP_ANGLE_THRESHOLD):
        """
        Check whether the instance is almost vertical.

        Examples:
          >>> Orien.fromAzPl(10, 15).isSubVertical()
          False
          >>> Orien.fromAzPl(257, 89).isSubVertical()
          True
        """

        return abs(self.pl.d) > min_dip_angle

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two Vect instances or a Vect and a GAXis instances.
        Range is 0°-180°.

        Examples:
          >>> are_close(Orien.fromAzPl(0, 90).angle(Orien.fromAzPl(90, 0)), 90)
          True
          >>> are_close(Orien.fromAzPl(0, 0).angle(Orien.fromAzPl(270, 0)), 90)
          True
          >>> are_close(Orien.fromAzPl(0, 0).angle(Orien.fromAzPl(0, 0)), 0)
          True
          >>> are_close(Orien.fromAzPl(0, 0).angle(Orien.fromAzPl(180, 0)), 180)
          True
          >>> are_close(Orien.fromAzPl(90, 0).angle(Orien.fromAzPl(270, 0)), 180)
          True
        """

        angle_vers = self.asVersor().angle(another.asVersor())

        return angle_vers

    def isAlmostParallel(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two Orien instances are sub-parallel,

        :param another: an Orien instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Orien.fromAzPl(0, 90).isAlmostParallel(Orien.fromAzPl(90, 0))
          False
          >>> Orien.fromAzPl(0, 0).isAlmostParallel(Orien.fromAzPl(0, 1e-6))
          True
          >>> Orien.fromAzPl(0, 90).isAlmostParallel(Orien.fromAzPl(180, 0))
          False
          >>> Orien.fromAzPl(0, 90).isAlmostParallel(Orien.fromAzPl(0, -90))
          False
        """

        fst_gvect = self

        snd_geoelem = another

        angle = fst_gvect.angle(snd_geoelem)

        if isinstance(another, PPlane):
            return angle > (90.0 - angle_tolerance)
        else:
            return angle <= angle_tolerance

    def isAlmostAntiParallel(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two Vect instances are almost anti-parallel,

        :param another: a Vect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Orien.fromAzPl(0, 90).isAlmostAntiParallel(Orien.fromAzPl(90, -89.5))
          True
          >>> Orien.fromAzPl(0, 0).isAlmostAntiParallel(Orien.fromAzPl(180, 1e-6))
          True
          >>> Orien.fromAzPl(90, 45).isAlmostAntiParallel(Orien.fromAzPl(270, -45.5))
          True
          >>> Orien.fromAzPl(45, 90).isAlmostAntiParallel(Orien.fromAzPl(0, -90))
          True
          >>> Orien.fromAzPl(45, 72).isAlmostAntiParallel(Orien.fromAzPl(140, -38))
          False
        """

        return self.angle(another) > (180.0 - angle_tolerance)

    def isSubOrthogonal(self, another, angle_tolerance=VECTOR_ANGLE_THRESHOLD):
        """
        Check that two Vect instance are sub-orthogonal

        :param another: a Vect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees) from orthogonality
        :return: Boolean

         Examples:
          >>> Orien.fromAzPl(0, 90).isSubOrthogonal(Orien.fromAzPl(90, 0))
          True
          >>> Orien.fromAzPl(0, 0).isSubOrthogonal(Orien.fromAzPl(0, 1.e-6))
          False
          >>> Orien.fromAzPl(0, 0).isSubOrthogonal(Orien.fromAzPl(180, 0))
          False
          >>> Orien.fromAzPl(90, 0).isSubOrthogonal(Orien.fromAzPl(270, 89.5))
          True
          >>> Orien.fromAzPl(0, 90).isSubOrthogonal(Orien.fromAzPl(0, 0.5))
          True
        """

        return 90.0 - angle_tolerance <= self.angle(another) <= 90.0 + angle_tolerance

    def normVersor(self, another):
        """
        Calculate the versor (Vect) defined by the vector product of two Orien instances.

        Examples:
          >>> Orien.fromAzPl(0, 0).normVersor(Orien.fromAzPl(90, 0))
          Vect(0.0000, 0.0000, -1.0000)
          >>> Orien.fromAzPl(45, 0).normVersor(Orien.fromAzPl(310, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> Orien.fromAzPl(0, 0).normVersor(Orien.fromAzPl(90, 90))
          Vect(-1.0000, 0.0000, -0.0000)
          >>> Orien.fromAzPl(315, 45).normVersor(Orien.fromAzPl(315, 44.5)) is None
          True
        """

        if self.isAlmostParallel(another):
            return None
        else:
            return self.asVersor().vCross(another.asVersor()).versor()

    def normPPlane(self):
        """
        Return the geological plane that is normal to the orientation.

        Examples:
          >>> Orien.fromAzPl(0, 45).normPPlane()
          PPlane(180.00, +45.00)
          >>> Orien.fromAzPl(0, -45).normPPlane()
          PPlane(000.00, +45.00)
          >>> Orien.fromAzPl(0, 90).normPPlane()
          PPlane(180.00, +00.00)
        """

        down_orien = self.downward()
        dipdir = (down_orien.az.d + 180.0) % 360.0
        dipangle = 90.0 - down_orien.pl.d

        return PPlane(dipdir, dipangle)

    def commonPPlane(self, another):
        """
        Calculate PPlane instance defined by the two Vect instances.

        Examples:
          >>> Orien.fromAzPl(0, 0).commonPPlane(Orien.fromAzPl(90, 0)).isAlmostParallel(PPlane(180.0, 0.0))
          True
          >>> Orien.fromAzPl(0, 0).commonPPlane(Orien.fromAzPl(90, 90)).isAlmostParallel(PPlane(90.0, 90.0))
          True
          >>> Orien.fromAzPl(45, 0).commonPPlane(Orien.fromAzPl(135, 45)).isAlmostParallel(PPlane(135.0, 45.0))
          True
          >>> Orien.fromAzPl(315, 45).commonPPlane(Orien.fromAzPl(135, 45)).isAlmostParallel(PPlane(225.0, 90.0))
          True
          >>> Orien.fromAzPl(0, 0).commonPPlane(Orien.fromAzPl(90, 0)).isAlmostParallel(PPlane(180.0, 10.0))
          False
          >>> Orien.fromAzPl(0, 0).commonPPlane(Orien.fromAzPl(90, 90)).isAlmostParallel(PPlane(90.0, 80.0))
          False
          >>> Orien.fromAzPl(45, 0).commonPPlane(Orien.fromAzPl(135, 45)).isAlmostParallel(PPlane(125.0, 45.0))
          False
          >>> Orien.fromAzPl(315, 45).commonPPlane(Orien.fromAzPl(135, 45)).isAlmostParallel(PPlane(225.0, 80.0))
          False
          >>> Orien.fromAzPl(315, 45).commonPPlane(Orien.fromAzPl(315, 44.5)) is None
          True
        """

        normal_versor = self.normVersor(another)
        if normal_versor is None:
            return None
        else:
            return normal_versor.asOrien().normPPlane()

    def asAxis(self):
        """
        Create Axis instance with the same attitude as the self instance.

        Example:
          >>> Orien.fromAzPl(220, 32).asAxis()
          Axis(az: 220.00°, pl: 32.00°)
        """

        return Axis(self.az, self.pl)

    def normOrien(self, another):
        """
        Calculate the Orien instance that is normal to the two provided sources.
        Angle between sources must be larger than MIN_ANGLE_DEGR_DISORIENTATION,
        otherwise a SubparallelLineationException will be raised.

        Example:
          >>> Orien.fromAzPl(0, 0).normOrien(Orien.fromAzPl(0.5, 0)) is None
          True
          >>> Orien.fromAzPl(0, 0).normOrien(Orien.fromAzPl(179.5, 0)) is None
          True
          >>> Orien.fromAzPl(0, 0).normOrien(Orien.fromAzPl(5.1, 0))
          Orien(az: 0.00°, pl: 90.00°)
          >>> Orien.fromAzPl(90, 45).normOrien(Orien.fromAzPl(90, 0))
          Orien(az: 180.00°, pl: 0.00°)

        """

        if self.isAlmostAntiParallel(another):
            return None
        elif self.isAlmostParallel(another):
            return None
        else:
            return self.normVersor(another).asOrien()


class Axis(Orien):
    """
    Polar Axis. Inherits from Orientation
    """

    def __init__(self, az: Azim, pl: Plunge):

        super().__init__(az, pl)

    def __repr__(self):

        return "Axis(az: {:.2f}°, pl: {:.2f}°)".format(*self.d)

    def asOrien(self):
        """
        Create OrienM instance with the same attitude as the self instance.

        Example:
          >>> Axis.fromAzPl(220, 32).asOrien()
          Orien(az: 220.00°, pl: 32.00°)
        """

        return Orien(self.az, self.pl)

    def normAxis(self, another):
        """
        Calculate the Axis instance that is perpendicular to the two provided.
        The two source Axis must not be subparallel (threshold is MIN_ANGLE_DEGR_DISORIENTATION),
        otherwise a SubparallelLineationException will be raised.

        Example:
          >>> Axis.fromAzPl(0, 0).normAxis(Axis.fromAzPl(0.5, 0)) is None
          True
          >>> Axis.fromAzPl(0, 0).normAxis(Axis.fromAzPl(180, 0)) is None
          True
          >>> Axis.fromAzPl(90, 0).normAxis(Axis.fromAzPl(180, 0))
          Axis(az: 0.00°, pl: 90.00°)
          >>> Axis.fromAzPl(90, 45).normAxis(Axis.fromAzPl(180, 0))
          Axis(az: 270.00°, pl: 45.00°)
          >>> Axis.fromAzPl(270, 45).normAxis(Axis.fromAzPl(180, 90)).isAlmostParallel(Axis.fromAzPl(180, 0))
          True
        """

        norm_orien = self.normOrien(another)
        if norm_orien is None:
            return None
        else:
            return norm_orien.asAxis()

    def angle(self, another):
        """
        Calculate angle (in degrees) between the two Axis instances.
        Range is 0°-90°.

        Examples:
          >>> are_close(Axis.fromAzPl(0, 90).angle(Axis.fromAzPl(90, 0)), 90)
          True
          >>> are_close(Axis.fromAzPl(0, 0).angle(Axis.fromAzPl(270, 0)), 90)
          True
          >>> are_close(Axis.fromAzPl(0, 0).angle(Axis.fromAzPl(0, 0)), 0)
          True
          >>> are_close(Axis.fromAzPl(0, 0).angle(Axis.fromAzPl(180, 0)), 0)
          True
          >>> are_close(Axis.fromAzPl(0, 0).angle(Axis.fromAzPl(179, 0)), 1)
          True
          >>> are_close(Axis.fromAzPl(0, -90).angle(Axis.fromAzPl(0, 90)), 0)
          True
          >>> are_close(Axis.fromAzPl(90, 0).angle(Axis.fromAzPl(315, 0)), 45)
          True
        """

        angle_vers = self.asVersor().angle(another.asVersor())

        return min(angle_vers, 180.0 - angle_vers)


class OrienM(object):
    """
    Polar vector class.
    """

    def __init__(self, orien: Orien, mag: [int, float], unit_m: str='m'):
        """
        Constructs a polar vector object.
        """

        if not isinstance(orien, Orien):
            raise GeomInputException("First vect argument must be of type Orien")

        if not isinstance(mag, (int, float)):
            raise GeomInputException("Second vect arg must be int/float")
        elif not isfinite(mag):
            raise GeomInputException("Second vect arg must be finite")
        elif mag <= 0.0:
            raise GeomInputException("Second vect argument must be positive")

        if not isinstance(unit_m, str):
            raise GeomInputException("Magnitude unit must be string")

        self.o = orien
        self.mag = float(mag)
        self.um = unit_m

    @property
    def d(self):
        """
        Returns azimuth and plunge of orientation in decimal degrees, as a tuple.

        :return: tuple of azimuth and plunge in decimal degrees

        Example:
          >>> OrienM.fromAzPlMg(100, 20, 1).d
          (100.0, 20.0)
          >>> OrienM.fromAzPlMg(-pi/2, -pi/4, 0.5, unit_a='r').d
          (270.0, -45.0)
        """

        return self.o.d

    @property
    def r(self):
        """
        Returns azimuth and plunge of orientation in radians, as a tuple.

        :return: tuple of azimuth and plunge in radians

        Example:
          >>> OrienM.fromAzPlMg(90, 45, 1).r
          (1.5707963267948966, 0.7853981633974483)
        """

        return self.o.r

    @property
    def m(self):
        """
        Returns magnitude and its measurement unit, as a tuple.

        :return: tuple of magnitude and measurement unit

        Example:
          >>> OrienM.fromAzPlMg(90, 45, 1, unit_m='mm').m
          (1.0, 'mm')
        """

        return self.mag, self.um

    @classmethod
    def fromAzPlMg(cls, az: [int, float], pl: [int, float], mag: [int, float], unit_a='d', unit_m='m') -> 'OrienM':
        """
        Class constructor from azimuth, plunge and magnitude.

        :param az: azimuth value
        :param pl: plunge value
        :param magn: magnitude
        :param unit_d: angle measurement unit, in degrees ('d') or radians ('r')
        :param unit_m: magnitude measurement unit. Default is meter ('m')
        :return: OrienM instance

        Examples:
          >>> OrienM.fromAzPlMg(30, 40, 0.1)
          OrienM(az: 30.00°, pl: 40.00°, mag: 0.1000 m)
        """

        az = Azim(az, unit=unit_a)
        pl = Plunge(pl, unit=unit_a)
        orien = Orien(az, pl)

        return cls(orien, mag, unit_m)
        
    @classmethod
    def fromXYZM(cls, x: [int, float], y: [int, float], z: [int, float], unit_m='m') -> 'OrienM':
        """
        Constructs a Vect object given a x-y-z triplet
        :param x: x component
        :param y: y component
        :param z: z component
        :return: the Vect instance

        Examples:
          >>> OrienM.fromXYZM(1, 1, 0, 'mm')
          OrienM(az: 45.00°, pl: -0.00°, mag: 1.4142 mm)
        """

        mag, norm_xyz = normXYZ(x, y, z)

        if norm_xyz is None:
            raise GeomInputException("Input components have near-zero values")

        orientation = Orien._from_xyz(*norm_xyz)

        return cls(orientation, mag, unit_m)

    def __repr__(self) -> str:

        return "OrienM(az: {:.2f}°, pl: {:.2f}°, mag: {:.4f} {})".format(*self.d, *self.m)


    def toXYZM(self) -> Tuple[float, float, float, str]:
        """
        Converts a polar vector to a tuple of x, y and z cartesian components, with measurement unit.

        :return: tuple of x, y and z components, and measurement unit.

        Examples:
          >>> OrienM.fromAzPlMg(90, 0, 0.2).toXYZM()
          (0.2, 1.2246467991473533e-17, -0.0, 'm')
        """

        vals = self.o.toXYZ()
        mag, unit_m = self.m
        x, y, z = map(lambda val: val*mag, vals)
        return x, y, z, unit_m


class PPlane(object):
    """
    Geological plane.
    Defined by dip direction and dip angle (both in degrees):
     - dip direction: [0.0, 360.0[ clockwise, from 0 (North);
     - dip angle: [0, 90.0]: downward-pointing.
    """

    def __init__(self, azim: float, dip_ang: float, is_rhr_strike=False):
        """
        Geological plane constructor.

        :param  azim:  azimuth of the plane (RHR strike or dip direction).
        :type  azim:  number or string convertible to float.
        :param  dip_ang:  Dip angle of the plane (0-90°).
        :type  dip_ang:  number or string convertible to float.
        :param is_rhr_strike: if the source azimuth is RHR strike (default is False, i.e. it is dip direction)
        :return: the instantiated geological plane.
        :rtype: PPlane.

        Example:
          >>> PPlane(0, 90)
          PPlane(000.00, +90.00)
          >>> PPlane(0, 90, is_rhr_strike=True)
          PPlane(090.00, +90.00)
          >>> PPlane(0, 90, True)
          PPlane(090.00, +90.00)
          >>> PPlane(0, "90", True)
          Traceback (most recent call last):
          ...
          GeomInputException: Source dip angle must be number
          >>> PPlane(0, 900)
          Traceback (most recent call last):
          ...
          GeomInputException: Dip angle must be between 0° and 90°
        """

        def rhrstrk2dd(rhr_strk):
            """Converts RHR strike value to dip direction value.

            Example:
                >>> rhrstrk2dd(285.5)
                15.5
            """

            return (rhr_strk + 90.0) % 360.0

        if not isinstance(azim, (int, float)):
            raise GeomInputException("Source azimuth must be number")
        if not isinstance(dip_ang, (int, float)):
            raise GeomInputException("Source dip angle must be number")
        if not isinstance(is_rhr_strike, bool):
            raise GeomInputException("Source azimuth type must be boolean")

        if not (0.0 <= dip_ang <= 90.0):
            raise GeomInputException("Dip angle must be between 0° and 90°")

        if is_rhr_strike:
            self._dipdir = rhrstrk2dd(azim)
        else:
            self._dipdir = azim % 360.0
        self._dipangle = float(dip_ang)

    @property
    def dd(self):
        """
        Return the dip direction of the geological plane.

        Example:
          >>> PPlane(34.2, 89.7).dd
          34.2
        """

        return self._dipdir

    @property
    def da(self):
        """
        Return the dip angle of the geological plane.

        Example:
          >>> PPlane(183, 77).da
          77.0

        """

        return self._dipangle

    @property
    def dda(self):
        """
        Return a tuple storing the dip direction and dip angle values of a geological plane.

        Example:
          >>> gp = PPlane(89.4, 17.2)
          >>> gp.dda
          (89.4, 17.2)
        """

        return self.dd, self.da

    @property
    def rhrStrike(self):
        """
        Return the strike according to the right-hand-rule.

        Examples:
          >>> PPlane(90, 45).rhrStrike
          0.0
          >>> PPlane(45, 89).rhrStrike
          315.0
          >>> PPlane(275, 38).rhrStrike
          185.0
          >>> PPlane(0, 38).rhrStrike
          270.0
        """

        return (self.dd - 90.0) % 360.0

    @property
    def srda(self):
        """
        Return a tuple storing the right-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> PPlane(100, 17.2).srda
          (10.0, 17.2)
          >>> PPlane(10, 87).srda
          (280.0, 87.0)
        """

        return self.rhrStrike, self.da

    @property
    def lhrStrike(self):
        """
        Return the strike according to the left-hand-rule.

        Examples:
          >>> PPlane(90, 45).lhrStrike
          180.0
          >>> PPlane(45, 89).lhrStrike
          135.0
          >>> PPlane(275, 38).lhrStrike
          5.0
          >>> PPlane(0, 38).lhrStrike
          90.0
        """

        return (self.dd + 90.0) % 360.0

    @property
    def slda(self):
        """
        Return a tuple storing the left-hand-rule strike and dip angle values of a geological plane.

        Example:
          >>> PPlane(100, 17.2).slda
          (190.0, 17.2)
          >>> PPlane(10, 87).slda
          (100.0, 87.0)
        """

        return self.lhrStrike, self.da

    def __repr__(self):

        return "PPlane({:06.2f}, {:+06.2f})".format(*self.dda)

    def rhrStrikeOrien(self):
        """
        Creates a OrienM instance that is parallel to the right-hand rule strike.

        :return: OrienM instance,

        Examples:
          >>> PPlane(90, 45).rhrStrikeOrien()
          Orien(az: 0.00°, pl: 0.00°)
          >>> PPlane(45, 17).rhrStrikeOrien()
          Orien(az: 315.00°, pl: 0.00°)
          >>> PPlane(90, 0).rhrStrikeOrien()
          Orien(az: 0.00°, pl: 0.00°)
        """

        return Orien.fromAzPl(
            az=self.rhrStrike,
            pl=0.0)

    def lhrStrikeOrien(self):
        """
        Creates an Orientation instance that is parallel to the left-hand rule strike.

        :return: OrienM instance.

        Examples:
          >>> PPlane(90, 45).lhrStrikeOrien()
          Orien(az: 180.00°, pl: 0.00°)
          >>> PPlane(45, 17).lhrStrikeOrien()
          Orien(az: 135.00°, pl: 0.00°)
        """

        return Orien.fromAzPl(
            az=self.lhrStrike,
            pl=0.0)

    def dipDirOrien(self):
        """
        Creates a OrienM instance that is parallel to the dip direction.

        :return: OrienM instance.

        Examples:
          >>> PPlane(90, 45).dipDirOrien()
          Orien(az: 90.00°, pl: 45.00°)
          >>> PPlane(45, 17).dipDirOrien()
          Orien(az: 45.00°, pl: 17.00°)
        """

        return Orien.fromAzPl(
            az=self.dd,
            pl=self.da)

    def dipDirOppOrien(self):
        """
        Creates a OrienM instance that is anti-parallel to the dip direction.

        :return: OrienM instance.

        Examples:
          >>> PPlane(90, 45).dipDirOppOrien()
          Orien(az: 270.00°, pl: -45.00°)
          >>> PPlane(45, 17).dipDirOppOrien()
          Orien(az: 225.00°, pl: -17.00°)
        """

        return self.dipDirOrien().opposite()

    def mirrorVertPPlane(self):
        """
        Mirror a geological plane around a vertical plane
        creating a new one that has a dip direction opposite
        to the original one but with downward plunge.

        :return: geological plane
        :rtype: PPlane

        Examples:
          >>> PPlane(0, 45).mirrorVertPPlane()
          PPlane(180.00, +45.00)
          >>> PPlane(225, 80).mirrorVertPPlane()
          PPlane(045.00, +80.00)
          >>> PPlane(90, 90).mirrorVertPPlane()
          PPlane(270.00, +90.00)
          >>> PPlane(270, 0).mirrorVertPPlane()
          PPlane(090.00, +00.00)
        """

        return PPlane(
            azim=opposite_trend(self.dd),
            dip_ang=self.da)

    def _normal_orien_frwrd(self):
        """
        Return the geological vector normal to the geological plane,
        pointing in the same direction as the geological plane.

        Example:
            >>> PPlane(90, 55)._normal_orien_frwrd()
            Orien(az: 90.00°, pl: -35.00°)
            >>> PPlane(90, 90)._normal_orien_frwrd()
            Orien(az: 90.00°, pl: 0.00°)
            >>> PPlane(90, 0)._normal_orien_frwrd()
            Orien(az: 90.00°, pl: -90.00°)
        """

        tr = self.dd % 360.0
        pl = self.da - 90.0

        return Orien.fromAzPl(
            az=tr,
            pl=pl)

    def _normal_orien_anti(self):
        """
        Return the geological vector normal to the geological plane,
        pointing in the opposite direction to the geological plane.

        Example:
            >>> PPlane(90, 55)._normal_orien_anti()
            Orien(az: 270.00°, pl: 35.00°)
            >>> PPlane(90, 90)._normal_orien_anti()
            Orien(az: 270.00°, pl: -0.00°)
            >>> PPlane(90, 0)._normal_orien_anti()
            Orien(az: 270.00°, pl: 90.00°)
        """

        return self._normal_orien_frwrd().opposite()

    def downNormOrien(self):
        """
        Return the geological vector normOrien to the geological plane,
        pointing downward.

        Example:
            >>> PPlane(90, 55).downNormOrien()
            Orien(az: 270.00°, pl: 35.00°)
            >>> PPlane(90, 90).downNormOrien()
            Orien(az: 90.00°, pl: 0.00°)
            >>> PPlane(90, 0).downNormOrien()
            Orien(az: 270.00°, pl: 90.00°)
        """

        return self._normal_orien_frwrd().downward()

    def upNormOrien(self):
        """
        Return the orientation normal to the polar plane,
        pointing upward.

        Example:
            >>> PPlane(90, 55).upNormOrien()
            Orien(az: 90.00°, pl: -35.00°)
            >>> PPlane(90, 90).upNormOrien()
            Orien(az: 90.00°, pl: 0.00°)
            >>> PPlane(90, 0).upNormOrien()
            Orien(az: 90.00°, pl: -90.00°)
        """

        return self._normal_orien_frwrd().upward()

    def normOrien(self):
        """
        Wrapper to down_normal_gv.

        :return: OrienM normOrien to the PPlane self instance
        """

        return self.downNormOrien()

    def normAxis(self):
        """
        Normal Axis.

        :return: Axis normal to the PPlane self instance
        """

        return self.downNormOrien().asAxis()

    def angle(self, another):
        """
        Calculate angle (in degrees) between two geoplanes.
        Range is 0°-90°.

        Examples:
          >>> PPlane(100.0, 50.0).angle(PPlane(100.0, 50.0))
          0.0
          >>> PPlane(300.0, 10.0).angle(PPlane(300.0, 90.0))
          80.0
          >>> PPlane(90.0, 90.0).angle(PPlane(270.0, 90.0))
          0.0
          >>> are_close(PPlane(90.0, 90.0).angle(PPlane(130.0, 90.0)), 40)
          True
          >>> are_close(PPlane(90, 70).angle(PPlane(270, 70)), 40)
          True
          >>> are_close(PPlane(90.0, 10.0).angle(PPlane(270.0, 10.0)), 20.0)
          True
          >>> are_close(PPlane(90.0, 10.0).angle(PPlane(270.0, 30.0)), 40.0)
          True
        """

        gpl_axis = self._normal_orien_frwrd().asAxis()
        if isinstance(another, PPlane):
            an_axis = another._normal_orien_frwrd().asAxis()
        else:
            raise GeomInputException("Provided another instance for angle is of {} type".format(type(another)))

        angle = gpl_axis.angle(an_axis)

        if isinstance(another, PPlane):
            return angle
        else:
            return 90.0 - angle

    def isAlmostParallel(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two GPlanes are sub-parallel

        :param another: a PPlane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> PPlane(0, 90).isAlmostParallel(PPlane(270, 90))
          False
          >>> PPlane(0, 90).isAlmostParallel(PPlane(180, 90))
          True
          >>> PPlane(0, 90).isAlmostParallel(PPlane(0, 0))
          False
          >>> PPlane(0, 0).isAlmostParallel(PPlane(0, 1e-6))
          True
          >>> PPlane(0, 0).isAlmostParallel(PPlane(0, 1.1))
          False
        """

        return self.angle(another) < angle_tolerance

    def isSubOrthogonal(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two GPlanes are sub-orthogonal.

        :param another: a PPlane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> PPlane(0, 90).isSubOrthogonal(PPlane(270, 90))
          True
          >>> PPlane(0, 90).isSubOrthogonal(PPlane(180, 90))
          False
          >>> PPlane(0, 90).isSubOrthogonal(PPlane(0, 0))
          True
          >>> PPlane(0, 0).isSubOrthogonal(PPlane(0, 88))
          False
          >>> PPlane(0, 0).isSubOrthogonal(PPlane(0, 45))
          False
        """

        fst_axis = self.normOrien().asAxis()

        if isinstance(another, PPlane):
            snd_gaxis = another.normOrien().asAxis()
        else:
            raise GeomInputException("Not accepted argument type for isSubOrthogonal method")

        angle = fst_axis.angle(snd_gaxis)

        if isinstance(another, PPlane):
            return angle > 90.0 - angle_tolerance
        else:
            return angle < angle_tolerance

    def rakeToOrien(self, rake):
        """
        Calculate OrienM given a PPlane instance and a rake value.
        The rake is defined according to the Aki and Richards, 1980 conventions:
        rake = 0° -> left-lateral
        rake = 90° -> reverse
        rake = +/- 180° -> right-lateral
        rake = -90° -> normal

        Examples:
          >>> PPlane(180, 45).rakeToOrien(0.0)
          Orien(az: 90.00°, pl: 0.00°)
          >>> PPlane(180, 45).rakeToOrien(90.0)
          Orien(az: 0.00°, pl: -45.00°)
          >>> PPlane(180, 45).rakeToOrien(-90.0)
          Orien(az: 180.00°, pl: 45.00°)
          >>> PPlane(180, 45).rakeToOrien(180.0).isAlmostParallel(Orien.fromAzPl(270.00, 0.00))
          True
          >>> PPlane(180, 45).rakeToOrien(-180.0)
          Orien(az: 270.00°, pl: 0.00°)
        """

        rk = radians(rake)
        strk = radians(self.rhrStrike)
        dip = radians(self.da)

        x = cos(rk) * sin(strk) - sin(rk) * cos(dip) * cos(strk)
        y = cos(rk) * cos(strk) + sin(rk) * cos(dip) * sin(strk)
        z = sin(rk) * sin(dip)

        return Vect(x, y, z).asOrien()

    def isVLowAngle(self, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a geological plane is very low angle.

        :param threshold: the limit for the plane angle, in degrees
        :type threshold: float
        :return: bool flag indicating if it is very low angle

        Examples:
          >>> PPlane(38.9, 1.2).isVLowAngle()
          True
          >>> PPlane(38.9, 7.4).isVLowAngle()
          False
        """

        return self.da < dip_angle_threshold

    def isVHighAngle(self, dip_angle_threshold=angle_gplane_thrshld):
        """
        Checks if a geological plane is very high angle.

        :param threshold: the limit for the plane angle, in degrees
        :type threshold: float
        :return: bool flag indicating if it is very high angle

        Examples:
          >>> PPlane(38.9, 11.2).isVHighAngle()
          False
          >>> PPlane(38.9, 88.4).isVHighAngle()
          True
        """

        return self.da > (90.0 - dip_angle_threshold)


class Point(object):
    """
    Cartesian point.
    Dimensions: 3D (space)
    """

    def __init__(self, x: [int, float], y: [int, float], z: [int, float]):
        """
        Construct a Point instance.
        """

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, (int, float)), vals)):
            raise GeomInputException("Input values must be integer of float")
        elif not all(map(isfinite, vals)):
            raise GeomInputException("Input values must be finite")
        else:
            self._a = array(vals, dtype=np.float64)

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self, another: 'Point') -> Optional[bool]:
        """
        Return True if objects are equal.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, -1)
          False
        """

        if not isinstance(another, Point):
            raise GeomInputException("Variables must be of the same type")
        else:
            return all([
                self.x == another.x,
                self.y == another.y,
                self.z == another.z])

    def __ne__(self, another: 'Point') -> Optional[bool]:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
        """

        if not isinstance(another, Point):
            return None
        else:
            return not (self == another)

    @property
    def a(self) -> 'numpy.array':
        """
        Return a copy of the object inner array.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7).a
          array([4., 3., 7.])
        """

        return np.copy(self._a)

    @property
    def x(self) -> float:
        """
        Return x value

        Example:
          >>> Point(1.5, 1, 1).x
          1.5
        """

        return self.a[0]

    @property
    def y(self) -> float:
        """
        Return y value

        Example:
          >>> Point(1.5, 3.0, 1).y
          3.0
        """
        return self.a[1]

    @property
    def z(self) -> float:
        """
        Return z value

        Example:
          >>> Point(1.5, 3.2, 41.).z
          41.0
        """
        return self.a[2]

    def toXYZ(self) -> Tuple[float, float, float]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> 'numpy.array':
        """
        Return a double Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> Point(1, 2, 3).toArray()
          array([1., 2., 3.])
        """

        return self.a

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return self.__class__(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return self.__class__(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return self.__class__(0.0, self.y, self.z)

    @property
    def len3D(self) -> float:
        """
        Spatial distance of the point from the axis origin.

        :return: distance
        :rtype: float

        Examples:
          >>> Point(4.0, 3.0, 0.0).len3D
          5.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @property
    def len2D(self) -> float:
        """
        2D distance of the point from the axis origin.

        Example:
          >>> Point(3, 4, 0).len2D
          5.0
          >>> Point(12, 5, 3).len2D
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y)

    def deltaX(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.x - self.x

    def deltaY(self, another: 'Point') -> Optional[float]:
        """
        Delta between y components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.y - self.y

    def deltaZ(self, another: 'Point') -> Optional[float]:
        """
        Delta between x components of two Point Instances.

        :return: difference value
        :rtype: float

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return another.z - self.z

    def dist3DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate Euclidean spatial distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist3DWith(Point(4., 5., 1,))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist2DWith(self, another: 'Point') -> Optional[float]:
        """
        Calculate horizontal (2D) distance between two points.

        Examples:
          >>> Point(1., 1., 1.).dist2DWith(Point(4., 5., 7.))
          5.0
        """

        if not isinstance(another, Point):
            return None
        else:
            return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self, scale_factor: [int, float]) -> Optional['Point']:
        """
        Create a scaled object.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(np.nan) is None
          True
          >>> Point(1, 0, 1).scale(np.inf) is None
          True
        """

        if not isinstance(scale_factor, (int, float)):
            return None
        elif not isfinite(scale_factor):
            return None
        else:
            x, y, z = arr2tuple(self.a * scale_factor)
            return self.__class__(x, y, z)

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def isCoinc(self, another: 'Point', tolerance: [int, float] = MIN_SEPARATION_THRESHOLD) -> Optional[bool]:
        """
        Check spatial coincidence of two points

        Example:
          >>> Point(1., 0., -1.).isCoinc(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc(Point(1., 0., 0.))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4))
          True
          >>> Point(1.2, 7.4, 1.4).isCoinc(Point(1.2, 7.4, 1.4), tolerance=np.nan) is None
          True
        """

        if not isinstance(another, Point):
            return None
        elif not isinstance(tolerance, (int, float)):
            return None
        elif not isfinite(tolerance):
            return None
        else:
            distance_2d = self.dist2DWith(another)
            if np.isnan(distance_2d) or distance_2d > tolerance:
                return False
            else:
                distance_3d = self.dist3DWith(another)
                if np.isnan(distance_3d) or distance_3d > tolerance:
                    return False
                else:
                    return True

    def shift(self, sx: float, sy: float, sz: float) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
          >>> Point(1, 2, -1).shift(0.5, np.nan, 1.5) is None
          True
       """

        vals = [sx, sy, sz]
        if not all(map(lambda val: isinstance(val, (int, float)), vals)):
            return None
        elif not all(map(isfinite, vals)):
            return None
        else:
            return self.__class__(self.x + sx, self.y + sy, self.z + sz)

    def asVect(self) -> 'Vect':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0).asVect()
          Vect(1.0000, 1.0000, 0.0000)
          >>> Point(0.2, 1, 6).asVect()
          Vect(0.2000, 1.0000, 6.0000)
        """

        return Vect(self.x, self.y, self.z)


class Vect(Point):
    """
    Cartesian 3D vector.
    Right-handed rectangular Cartesian coordinate system (ENU):
    x axis -> East
    y axis -> North
    z axis -> Up
    """

    def __init__(self, x: [int, float], y: [int, float], z: [int, float]):
        """
        Vect constructor.

        Example;
          >>> Vect(1, 0, 1)
          Vect(1.0000, 0.0000, 1.0000)
          >>> Vect(1, np.nan, 1)
          Traceback (most recent call last):
          ...
          GeomInputException: Input values must be finite
          >>> Vect(1, 0, np.inf)
          Traceback (most recent call last):
          ...
          GeomInputException: Input values must be finite
          >>> Vect(0, 0, 0)
          Vect(0.0000, 0.0000, 0.0000)
        """

        super().__init__(x, y, z)

    def __repr__(self) -> str:

        return "Vect({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __add__(self, another: 'Vect') -> 'Vect':
        """
        Sum of two vectors.

        Example:
          >>> Vect(1, 0, 0) + Vect(0, 1, 1)
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(1, 1, 1) + Vect(-1, -1, -1)
          Vect(0.0000, 0.0000, 0.0000)
        """

        x, y, z = arr2tuple(self.a + another.a)
        return self.__class__(x, y, z)

    def __sub__(self, another: 'Vect') -> 'Vect':
        """Return object difference

        Example:
          >>> Vect(1., 1., 1.) - Vect(1., 1., 1.)
          Vect(0.0000, 0.0000, 0.0000)
          >>> Vect(1., 1., 3.) - Vect(1., 1., 2.2)
          Vect(0.0000, 0.0000, 0.8000)
        """

        x, y, z = arr2tuple(self.a - another.a)
        return self.__class__(x, y, z)

    @property
    def isAlmostZero(self) -> bool:
        """
        Check if the Vect instance length is near zero.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isAlmostZero
          False
          >>> Vect(0.0, 0.0, 0.0).isAlmostZero
          True
        """

        return are_close(self.len3D, 0)

    @property
    def isAlmostUnit(self) -> bool:
        """
        Check if the Vect instance length is near unit.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isAlmostUnit
          False
          >>> Vect(0.0, 1.0, 0.0).isAlmostUnit
          True
        """

        return are_close(self.len3D, 1)

    @property
    def isValid(self) -> bool:
        """
        Check if the Vect instance components are not all valid and the xyz not all zero-valued.

        :return: Boolean

        Example:
          >>> Vect(1, 2, 0).isValid
          True
          >>> Vect(0.0, 0.0, 0.0).isValid
          False
        """

        return not self.isAlmostZero

    def versor(self) -> Optional['Vect']:
        """
        Calculate versor in xyz space.

        Example:
          >>> Vect(5, 0, 0).versor()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, -1).versor()
          Vect(0.0000, 0.0000, -1.0000)
          >>> Vect(0, 0, 0).versor() is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.scale(1.0 / self.len3D)

    def versor2D(self) -> Optional['Vect']:
        """
        Create 2D versor version of the current vector

        :return: unit vector

        Example:
          >>> Vect(7, 0, 10).versor2D()
          Vect(1.0000, 0.0000, 0.0000)
          >>> Vect(0, 0, 10).versor2D() is None
          True
        """

        vXY = self.pXY()
        if vXY.isValid:
            return self.pXY().versor()
        else:
            return None

    def trend(self) -> Optional[float]:
        """
        Trend of a vector
        (degrees, clockwise from North, range 0°-360°)

        :return: an optional float value, representing trend, in decimal degrees.

        Examples:
          >>> Vect(1, 0, 0).trend()
          90.0
          >>> Vect(0, 1, 0).trend()
          0.0
          >>> Vect(1, 1, 0).trend()
          45.0
          >>> Vect(1, -1, 0).trend()
          135.0
          >>> Vect(0, -1, 0).trend()
          180.0
          >>> Vect(-1, -1, 0).trend()
          225.0
          >>> Vect(-1, 0, 0).trend()
          270.0
          >>> Vect(-1, 1, 0).trend()
          315.0
          >>> Vect(1, 1, 10).trend()
          45.0
          >>> Vect(0, 0, 0).trend() is None
          True
         """

        return angle_north_clock(self.x, self.y)

    def slope(self) -> Optional[float]:
        """
        Slope of a vector.
        Degrees, positive: downward-directed, negative: upward-dir., range -90°/90°

        :return: an optional float, representing the vector slope, in decimal degrees.

        Examples:
          >>> Vect(1, 0, -1).slope()
          45.0
          >>> Vect(1, 0, 1).slope()
          -45.0
          >>> Vect(0, 1, 0).slope()
          0.0
          >>> Vect(0, 0, 1).slope()
          -90.0
          >>> Vect(0, 0, -1).slope()
          90.0
          >>> Vect(0, 0, 0).slope() is None
          True
         """
        h = self.len2D
        zv = self.z
        sl = slope(h, abs(zv))

        if sl is None:
            return None
        else:
            if zv <= 0.0:
                return sl
            else:
                return -sl

    def absSlope(self) -> Optional[float]:
        """
        Return the absolute value of the slope.

        :return: optional float, slope in decimal degrees.

          >>> Vect(1, 0, -1).absSlope()
          45.0
          >>> Vect(1, 0, 1).absSlope()
          45.0
          >>> Vect(0, 1, 0).absSlope()
          0.0
          >>> Vect(0, 0, 1).absSlope()
          90.0
          >>> Vect(0, 0, -1).absSlope()
          90.0
          >>> Vect(0, 0, 0).absSlope() is None
          True
        """

        sl = self.slope()
        if sl is None:
            return None
        else:
            return abs(sl)

    @property
    def isUpward(self) -> Optional[bool]:
        """
        Check that a vector is upward-directed.

        :return: boolean

        Example:
          >>> Vect(0,0,1).isUpward
          True
          >>> Vect(0,0,-0.5).isUpward
          False
          >>> Vect(1, 3, 0).isUpward
          False
          >>> Vect(0, 0, 0).isUpward is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.z > 0.0

    @property
    def isDownward(self) -> Optional[bool]:
        """
        Check that a vector is downward-directed.

        :return: boolean

        Example:
          >>> Vect(0,0,1).isDownward
          False
          >>> Vect(0,0,-0.5).isDownward
          True
          >>> Vect(1, 3, 0).isDownward
          False
          >>> Vect(0, 0, 0).isDownward is None
          True
        """

        if not self.isValid:
            return None
        else:
            return self.z < 0.0

    def upward(self) -> Optional['Vect']:
        """
        Calculate a new upward-pointing vector.

        Example:
          >>> Vect(1, 1, 1).upward()
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(-1, -1, -1).upward()
          Vect(1.0000, 1.0000, 1.0000)
          >>> Vect(0, 0, 0).upward() is None
          True
        """

        if not self.isValid:
            return None
        elif self.z < 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def downward(self) -> Optional['Vect']:
        """
        Calculate a new vector downward-pointing.

        Example:
          >>> Vect(1, 1, 1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(-1, -1, -1).downward()
          Vect(-1.0000, -1.0000, -1.0000)
          >>> Vect(0, 0, 0).downward() is None
          True
        """

        if not self.isValid:
            return None
        elif self.z > 0.0:
            return self.scale(-1.0)
        else:
            return self.scale(1.0)

    def asOrien(self) -> Optional[Orien]:
        """
        Calculate the polar orientation parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(1, 1, 1).asOrien()
          Orien(az: 45.00°, pl: -35.26°)
          >>> Vect(0, 1, 1).asOrien()
          Orien(az: 0.00°, pl: -45.00°)
          >>> Vect(1, 0, 1).asOrien()
          Orien(az: 90.00°, pl: -45.00°)
          >>> Vect(0, 0, 1).asOrien()
          Orien(az: 0.00°, pl: -90.00°)
          >>> Vect(0, 0, -1).asOrien()
          Orien(az: 0.00°, pl: 90.00°)
          >>> Vect(-1, 0, 0).asOrien()
          Orien(az: 270.00°, pl: 0.00°)
          >>> Vect(0, -1, 0).asOrien()
          Orien(az: 180.00°, pl: 0.00°)
          >>> Vect(-1, -1, 0).asOrien()
          Orien(az: 225.00°, pl: 0.00°)
          >>> Vect(0, 0, 0).asOrien() is None
          True
        """

        if self.isValid:

            pl = self.slope()  # upward pointing -> negative value, downward -> positive

            unit_vect = self.versor()
            if unit_vect.y == 0. and unit_vect.x == 0:
                az = 0.
            else:
                az = (90. - degrees(atan2(unit_vect.y, unit_vect.x))) % 360.

            return Orien.fromAzPl(az, pl)

        else:

            return None

    def asAxis(self) -> Optional['Axis']:
        """
        Calculate the axis parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Vect(0, 1, 1).asAxis()
          Axis(az: 0.00°, pl: -45.00°)
          >>> Vect(1, 0, 1).asAxis()
          Axis(az: 90.00°, pl: -45.00°)
          >>> Vect(0, 0, 1).asAxis()
          Axis(az: 0.00°, pl: -90.00°)
          >>> Vect(0, 0, -1).asAxis()
          Axis(az: 0.00°, pl: 90.00°)
          >>> Vect(-1, 0, 0).asAxis()
          Axis(az: 270.00°, pl: 0.00°)
          >>> Vect(0, -1, 0).asAxis()
          Axis(az: 180.00°, pl: 0.00°)
          >>> Vect(-1, -1, 0).asAxis()
          Axis(az: 225.00°, pl: 0.00°)
          >>> Vect(0, 0, 0).asAxis() is None
          True
        """

        if self.isValid:
            return self.asOrien().asAxis()
        else:
            return None

    def vDot(self, another: 'Vect') -> float:
        """
        Vector scalar multiplication.

        Examples:
          >>> Vect(1, 0, 0).vDot(Vect(1, 0, 0))
          1.0
          >>> Vect(1, 0, 0).vDot(Vect(0, 1, 0))
          0.0
          >>> Vect(1, 0, 0).vDot(Vect(-1, 0, 0))
          -1.0
        """

        return self.x * another.x + self.y * another.y + self.z * another.z

    def angleCos(self, another: 'Vect') -> Optional[float]:
        """
        Return the cosine of the angle between two vectors.

        Examples:
          >>> Vect(1,0,0).angleCos(Vect(0,0,1))
          0.0
          >>> Vect(1,0,0).angleCos(Vect(-1,0,0))
          -1.0
          >>> Vect(1,0,0).angleCos(Vect(1,0,0))
          1.0
          >>> Vect(0, 0, 0).angleCos(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).angleCos(Vect(0,0,0)) is None
          True
        """

        if not (self.isValid and another.isValid):
            return None
        else:
            val = self.vDot(another) / (self.len3D * another.len3D)
            if val > 1.0:
                return 1.0
            elif val < -1.0:
                return -1.0
            else:
                return val

    def angle(self, another: 'Vect') -> Optional[float]:
        """
        Calculate angle between two vectors, as degrees
        in 0° - 180° range.

        Example:
          >>> Vect(1, 0, 0).angle(Vect(0, 0, 1))
          90.0
          >>> Vect(1, 0, 0).angle(Vect(-1, 0, 0))
          180.0
          >>> Vect(0, 0, 1).angle(Vect(0, 0, -1))
          180.0
          >>> Vect(1, 1, 1).angle(Vect(1, 1,1 ))
          0.0
          >>> Vect(0, 0, 0).angle(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).angle(Vect(0,0,0)) is None
          True
        """

        if not (self.isValid and another.isValid):
            return None
        else:
            return degrees(acos(self.angleCos(another)))

    def isAlmostParallel(self, another: 'Vect', angle_tolerance: [int, float] = VECTOR_ANGLE_THRESHOLD) -> Optional[
        bool]:
        """
        Check that two Vect are sub-parallel,

        :param another: a Vect instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

        Examples:
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(1, 0, 0))
          True
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(0, 0, 1))
          False
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(-1, 0, 0))
          False
          >>> Vect(0, 0, 0).isAlmostParallel(Vect(1,0,0)) is None
          True
          >>> Vect(1, 0, 0).isAlmostParallel(Vect(0,0,0)) is None
          True
        """

        if not isinstance(another, Vect):
            return None
        elif not (self.isValid and another.isValid):
            return None
        elif not isinstance(angle_tolerance, (int, float)):
            return None
        elif not isfinite(angle_tolerance):
            return None
        else:
            return self.angle(another) <= angle_tolerance

    def isSubOrthogonal(self, another: 'Vect', angle_tolerance: [int, float] = VECTOR_ANGLE_THRESHOLD) -> Optional[
        bool]:
        """
        Check whether two vectors are sub-orhogonal.

        :param another: a second Vect instance
        :param angle_tolerance: the tolerance angle, in decimal degrees
        :return: Boolean

        Example:
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 0))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 1, 1))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal(Vect(0, 0.9999999999999, 0))
          True
          >>> Vect(1, 0, 0).isSubOrthogonal((Vect(0, 0, 0))) is None
          True
        """

        if not isinstance(another, Vect):
            return None
        elif not (self.isValid and another.isValid):
            return None
        elif not isinstance(angle_tolerance, (int, float)):
            return None
        elif not isfinite(angle_tolerance):
            return None
        else:
            return are_close(0, self.angleCos(another), atol=cos(angle_tolerance))

    def vCross(self, another: 'Vect') -> 'Vect':
        """
        Vector product (cross product).

        Examples:
          >>> Vect(1, 0, 0).vCross(Vect(0, 1, 0))
          Vect(0.0000, 0.0000, 1.0000)
          >>> Vect(1, 0, 0).vCross(Vect(1, 0, 0))
          Vect(0.0000, 0.0000, 0.0000)
          >>> (Vect(1, 0, 0).vCross(Vect(-1, 0, 0))).isAlmostZero
          True
        """

        x, y, z = arr2tuple(np.cross(self.a[:3], another.a[:3]))
        return Vect(x, y, z)

    def byMatrix(self, array3x3: 'np.array') -> 'Vect':
        """
        Matrix multiplication of a vector.

        """

        x, y, z = arr2tuple(array3x3.dot(self.a))
        return Vect(x, y, z)


class CPlane(object):
    """
    Cartesian plane.
    Expressed by equation:
    ax + by + cz + d = 0

    Note: Plane is locational - its position in space is defined.
    This contrast with PPlane, defined just by its attitude, but with undefined position

    """

    def __init__(self, a, b, c, d):

        self._a = float(a)
        self._b = float(b)
        self._c = float(c)
        self._d = float(d)

    @property
    def a(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).a
          1.0
        """

        return self._a

    @property
    def b(self):
        """
        Return b coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 4, 0, 2).b
          4.0
        """

        return self._b

    @property
    def c(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 5.4, 2).c
          5.4
        """

        return self._c

    @property
    def d(self):
        """
        Return a coefficient of a Plane instance.

        Example:
          >>> CPlane(1, 0, 0, 2).d
          2.0
        """

        return self._d

    @property
    def v(self):
        """
        Return coefficients of a Plane instance.

        Example:
          >>> CPlane(1, 1, 7, -4).v
          (1.0, 1.0, 7.0, -4.0)
        """
        return self.a, self.b, self.c, self.d

    @classmethod
    def fromPoints(cls, pt1, pt2, pt3):
        """
        Create a Plane from three given Point instances.

        Example:
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))
          CPlane(0.0000, 0.0000, 1.0000, 0.0000)
          >>> CPlane.fromPoints(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))
          CPlane(1.0000, 0.0000, 0.0000, 0.0000)
        """

        matr_a = array([[pt1.y, pt1.z, 1],
                           [pt2.y, pt2.z, 1],
                           [pt3.y, pt3.z, 1]])

        matr_b = - array([[pt1.x, pt1.z, 1],
                             [pt2.x, pt2.z, 1],
                             [pt3.x, pt3.z, 1]])

        matr_c = array([[pt1.x, pt1.y, 1],
                           [pt2.x, pt2.y, 1],
                           [pt3.x, pt3.y, 1]])

        matr_d = - array([[pt1.x, pt1.y, pt1.z],
                             [pt2.x, pt2.y, pt2.z],
                             [pt3.x, pt3.y, pt3.z]])

        return cls(np.linalg.det(matr_a),
                   np.linalg.det(matr_b),
                   np.linalg.det(matr_c),
                   np.linalg.det(matr_d))

    @classmethod
    def fromGPlanePt(cls, gplane, point):
        """
        Given a PPlane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> CPlane.fromGPlanePt(PPlane(0, 0), Point(0, 0, 0))
          CPlane(0.0000, 0.0000, 1.0000, -0.0000)
          >>> CPlane.fromGPlanePt(PPlane(90, 45), Point(0, 0, 0))
          CPlane(0.7071, 0.0000, 0.7071, -0.0000)
          >>> CPlane.fromGPlanePt(PPlane(0, 90), Point(0, 0, 0))
          CPlane(0.0000, 1.0000, -0.0000, -0.0000)
        """

        normal_versor = gplane._normal_orien_frwrd().asVersor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return cls(a, b, c, d)

    def __repr__(self):

        return "CPlane({:.4f}, {:.4f}, {:.4f}, {:.4f})".format(*self.v)

    def normVersor(self):
        """
        Return the versor normal to the cartesian plane.

        Examples:
          >>> CPlane(0, 0, 5, -2).normVersor()
          Vect(0.0000, 0.0000, 1.0000)
          >>> CPlane(0, 7, 0, 5).normVersor()
          Vect(0.0000, 1.0000, 0.0000)
        """

        return Vect(self.a, self.b, self.c).versor()

    def toGPlanePoint(self):
        """
        Converts a cartesian plane into a geological plane
        and a point lying in the plane (non-unique solution).

        Examples:
          >>> gpl, pt = CPlane(0, 0, 1, -1).toGPlanePoint()
          >>> gpl
          PPlane(000.00, +00.00)
          >>> pt
          Point(0.0000, 0.0000, 1.0000)
        """

        geol_plane = self.normVersor().asOrien().normPPlane()
        point = Point(*point_solution(array([[self.a, self.b, self.c]]),
                                      array([-self.d])))
        return geol_plane, point

    def intersVersor(self, another):
        """
        Return intersection versor for two intersecting planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersVersor(b)
        Vect(0.0000, -1.0000, 0.0000)
        """

        return self.normVersor().vCross(another.normVersor()).versor()

    def intersPoint(self, another):
        """
        Return point on intersection line (obviously non-unique solution)
        for two planes.

        >>> a = CPlane(1, 0, 0, 0)
        >>> b = CPlane(0, 0, 1, 0)
        >>> a.intersPoint(b)
        Point(0.0000, 0.0000, 0.0000)
        """

        # find a point lying on the intersection line (this is a non-unique solution)
        a = array([[self.a, self.b, self.c], [another.a, another.b, another.c]])
        b = array([-self.d, -another.d])
        x, y, z = point_solution(a, b)

        return Point(x, y, z)

    def isPointInPlane(self, pt):
        """
          Check whether a point lie in a plane.

          >>> pl = CPlane(0, 0, 1, 0)
          >>> pt = Point(0, 1, 0)
          >>> pl.isPointInPlane(pt)
          True
        """

        if abs(self.a * pt.x + self.b * pt.y + self.c * pt.z + self.d) < MIN_SCALAR_VALUE:
            return True
        else:
            return False

    def angle(self, another):
        """
        Calculate angle (in degrees) between two planes.

        Examples:
          >>> CPlane(1,0,0,0).angle(CPlane(0,1,0,0))
          90.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,1,0))
          45.0
          >>> CPlane(1,0,0,0).angle(CPlane(1,0,0,0))
          0.0
        """

        angle_degr = self.normVersor().angle(another.normVersor())
        if abs(angle_degr) < MIN_ANGLE_DEGR_VALUE:
            angle_degr = 0.0
        elif angle_degr > 90.0:
            angle_degr = 180.0 - angle_degr

        return angle_degr

    def isAlmostParallel(self, another, angle_tolerance=PLANE_ANGLE_THRESHOLD):
        """
        Check that two Plane are sub-parallel

        :param another: a Plane instance
        :param angle_tolerance: the maximum allowed divergence angle (in degrees)
        :return: Boolean

         Examples:
          >>> CPlane(1,0,0,0).isAlmostParallel(CPlane(1,0,0,0))
          True
          >>> CPlane(1,0,0,0).isAlmostParallel(CPlane(1,0,1,0))
          False
        """

        return self.angle(another) < angle_tolerance


class GeomInputException(Exception):
    """
    Exception for geometric input.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
