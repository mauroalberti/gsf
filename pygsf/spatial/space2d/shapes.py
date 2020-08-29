import abc

import numbers
from math import pi


class Shape2D(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    @property
    def center(self):
        """Calculate shape center"""

    @abc.abstractmethod
    @property
    def area(self):
        """Calculate shape area"""
        
    @abc.abstractmethod
    @property
    def length(self):
        """Calculate shape area"""

    @abc.abstractmethod
    @property
    def clone(self):
        """Create a clone of the shape"""


class Point2D(Shape2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real
                 ):

        self._x = x
        self._y = y

    def clone(self):
        return Point2D(self._x, self._y)

    @property
    def center(self):
        return self.clone()

    @property
    def area(self):
        return 0.0

    @property
    def length(self):
        return 0.0


class Circle2D(Shape2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 r: numbers.Real
                 ):

        self._x = x
        self._y = y
        self._r = r

    @property
    def center(self):
        return Point2D(self._x, self._y)

    @property
    def radius(self):
        return self._r

    @property
    def area(self):
        return pi * self._r * self._r

    @property
    def length(self):
        return 2.0 * pi * self._r

    def clone(self):
        return Circle2D(self._x, self._y, self._r)


class Square2D(Shape2D):

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 side: numbers.Real,
                 rotation: numbers.Real
                 ):
        """

        :param x: x coordinate of square center
        :param y: y coordinate of square center
        :param side: square side
        :param rotation: square rotation, counter-clockwise, decimal degrees
        """

        self._x = x
        self._y = y
        self._side = side
        self._cc_rotat = rotation % 360.0

    @property
    def center(self):
        return Point2D(self._x, self._y)

    @property
    def area(self):
        return self._side * self._side

    @property
    def length(self):
        return 4.0 * self._side

    def clone(self):
        return Square2D(self._x, self._y, self._side, self._cc_rotat)
