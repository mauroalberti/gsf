import datetime
from math import sqrt
from typing import Optional

import numpy
from numpy.core._multiarray_umath import absolute

from qygsf import Vect
from qygsf.geometries.shapes.geometry import MIN_SEPARATION_THRESHOLD


class Point4D(object):
    """
    Cartesian point.
    Dimensions: 3D + time
    """

    def __init__(self,
                 x: Optional[np.ndarray] = np.nan,
                 y: Optional[np.ndarray] = np.nan,
                 z: Optional[np.ndarray] = np.nan,
                 t: Optional[datetime.datetime] = None
                 ):
        """
        Construct a Point instance given 3 float values and an optional time value (seconds from ....).
        """

        self._p = np.array([x, y, z], dtype=np.float64)
        self._t = t

    def __repr__(self):

        return f"Point({self.x:.4f}, {self.y:.4f}, {self.z:.4f}, {self.t:s})"

    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a point from a numpy 1x4 array.

        Example:
          >>> Point4D.from_array(np.array([1, 0, 1]))
          Point(1.0000, 0.0000, 1.0000, None)
        """

        obj = cls()

        assert 2 <= a.size <= 3
        b = a.astype(np.float64)
        if b.size == 2:
            c = np.append(b, [np.nan])
        else:
            c = b
        obj._p = c
        obj._t = None
        return obj

    '''
    @property
    def v(self):
        """
        Return values as array

        Example:
          >>> Point(1, 0, 0).v
          array([  1.,   0.,   0.,  nan])
        """

        return self._p
    '''

    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Point4D(1.5, 1, 1).x
          1.5
        """

        return self._p[0]

    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Point4D(1.5, 3.0, 1).y
          3.0
        """
        return self._p[1]

    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Point4D(1.5, 3.2, 41.).z
          41.0
        """
        return self._p[2]

    @property
    def t(self):
        """
        Return time value

        Example:
          >>> Point4D(1.5, 3.2, 41., 22.).t
          22.0
        """
        return self._t

    def clone(self):
        """
        Clone the point.

        Example:
          >>> Point4D(1, 1, 1).clone()
          Point(1.0000, 1.0000, 1.0000, nan)
        """

        return Point4D(
            x=self.x,
            y=self.y,
            z=self.z,
            t=self.t
        )

    def __sub__(self, another):
        """Return point difference

        Example:
          >>> Point4D(1., 1., 1.) - Point4D(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000, nan)
        """

        return Point4D(
            x=self.x - another.x,
            y=self.y - another.y,
            z=self.z - another.z if np.isfinite(self.z) and np.isfinite(another.z) else np.nan,
            t=self.t - another.t if self.t is not None and another.t is not None else None
        )

    def __abs__(self):
        """
        Point distance from frame origin.
        todo: make sure it works as intended with nan values

        Example:
          >>> abs(Point4D(3, 4, 0))
          5.0
          >>> abs(Point4D(0, 12, 5))
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def dist_3d(self, another):
        """
        Calculate Euclidean spatial distance between two points.
        todo: make sure it works as intended with nan values

        Examples:
          >>> Point4D(1., 1., 1.).dist_3d(Point4D(4., 5., 1,))
          5.0
          >>> Point4D(1, 1, 1, 4).dist_3d(Point4D(4, 5, 1, 14))
          5.0
        """

        return abs(self - another)

    def dist_2d(self, another):
        """
        Calculate horizontal (2D) distance between two points.
        todo: make sure it works as intended with nan values

        Examples:
          >>> Point4D(1.,1.,1.).dist_2d(Point4D(4.,5.,7.))
          5.0
        """

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def coincident(self, another, tolerance=MIN_SEPARATION_THRESHOLD):
        """
        Check spatial coincidence of two points
        todo: make sure it works as intended with Points with nan values

        Example:
          >>> Point4D(1., 0., -1.).coincident(Point4D(1., 1.5, -1.))
          False
          >>> Point4D(1., 0., 0.).coincident(Point4D(1., 0., 0.))
          True
        """

        if self.dist_2d(another) > tolerance:
            return False
        elif self.dist_3d(another) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0, st=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        Example:
          >>> Point4D(1, 1, 1).translate(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000, nan)
       """

        return Point4D(self.x + sx, self.y + sy, self.z + sz, self.t + st)

    def vect_offset(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector.

        Example:
          >>> Point4D(1, 2, 0).vect_offset(Vect(10, 5, 0))
          Point(11.0000, 7.0000, 0.0000, nan)
        """

        return Point4D(self.x + displ_vect.x,
                       self.y + displ_vect.y,
                       self.z + displ_vect.z,
                       self.t)

    @property
    def vector(self):
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point4D(1, 1, 0, 5).vector
          Vect(1.0000, 1.0000, 0.0000)
        """

        return Vect(self.x, self.y, self.z)

    def delta_time(self, another):
        """
        Calculate the time difference between two points

        Example:
          >>> Point4D(1,1,1,4).delta_time(Point4D(1,1,2,5))
          1.0
        """

        return another.t - self.t

    def speed(self, another):
        """
        Calculate the speed required to displace self to another.

        Example:
          >>> Point4D(1, 1, 1, 4).speed(Point4D(4, 5, 1, 14))
          0.5
        """

        try:
            return self.dist_3d(another) / self.delta_time(another)
        except:
            return np.Infinity