# -*- coding: utf-8 -*-


from __future__ import division

from math import sin, cos, radians, pi

from typing import Tuple

from pygsf.mathematics.arrays import *
from pygsf.spatial.vector.vector import *
from pygsf.orientations.utils import *


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

    @classmethod
    def fromVect(cls, vect: Vect) -> Optional['Orien']:
        """
        Calculate the polar orientation parallel to the Vect instance.
        Trend range: [0°, 360°[
        Plunge range: [-90°, 90°], with negative values for upward-pointing
        geological axes and positive values for downward-pointing axes.

        Examples:
          >>> Orien.fromVect(Vect(1, 1, 1))
          Orien(az: 45.00°, pl: -35.26°)
          >>> Orien.fromVect(Vect(0, 1, 1))
          Orien(az: 0.00°, pl: -45.00°)
          >>> Orien.fromVect(Vect(1, 0, 1))
          Orien(az: 90.00°, pl: -45.00°)
          >>> Orien.fromVect(Vect(0, 0, 1))
          Orien(az: 0.00°, pl: -90.00°)
          >>> Orien.fromVect(Vect(0, 0, -1))
          Orien(az: 0.00°, pl: 90.00°)
          >>> Orien.fromVect(Vect(-1, 0, 0))
          Orien(az: 270.00°, pl: 0.00°)
          >>> Orien.fromVect(Vect(0, -1, 0))
          Orien(az: 180.00°, pl: 0.00°)
          >>> Orien.fromVect(Vect(-1, -1, 0))
          Orien(az: 225.00°, pl: 0.00°)
          >>> Orien.fromVect(Vect(0, 0, 0)) is None
          True
        """

        x, y, z = vect.toXYZ()
        return Orien.fromXYZ(x, y, z)

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

    def toCPlane(self, point):
        """
        Given a PPlane instance and a provided Point instance,
        calculate the corresponding Plane instance.

        Example:
          >>> PPlane(0, 0).toCPlane(Point(0, 0, 0))
          CPlane(0.0000, 0.0000, 1.0000, -0.0000)
          >>> PPlane(90, 45).toCPlane(Point(0, 0, 0))
          CPlane(0.7071, 0.0000, 0.7071, -0.0000)
          >>> PPlane(0, 90).toCPlane(Point(0, 0, 0))
          CPlane(0.0000, 1.0000, -0.0000, -0.0000)
        """

        normal_versor = self._normal_orien_frwrd().asVersor()
        a, b, c = normal_versor.x, normal_versor.y, normal_versor.z
        d = - (a * point.x + b * point.y + c * point.z)
        return CPlane(a, b, c, d)


class GeomInputException(Exception):
    """
    Exception for geometric input.
    """

    pass


if __name__ == "__main__":

    import doctest
    doctest.testmod()
