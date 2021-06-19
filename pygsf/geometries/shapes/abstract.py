
from typing import List, Optional
import numbers
import abc

import math

from ...mathematics.defaults import *
from ...utils.types import *


class Shape(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def area(self):
        """Calculate shape area"""

    @abc.abstractmethod
    def length(self):
        """Calculate shape length"""

    @abc.abstractmethod
    def clone(self) -> 'Shape':
        """
        Clone a shape."""


class Point(Shape, metaclass=abc.ABCMeta):

    def __init__(self,
        x: numbers.Real,
        y: numbers.Real
    ):

        check_type(x, "Input x value", numbers.Real)
        check_type(y, "Input y value", numbers.Real)

        if not all(map(math.isfinite, [x, y])):
            raise Exception("Both x and y input values must be finite")

        self._x = float(x)
        self._y = float(y)

    @property
    def x(self):
        """X coordinate"""

        return self._x

    @property
    def y(self):
        """Y coordinate"""

        return self._y

    @abc.abstractmethod
    def as_point2d(self):
        """Converts a point to a point 2D"""

    def area(self):
        """Calculate shape area"""

        return 0.0

    def length(self):
        """Calculate shape length"""

        return 0.0

    @abc.abstractmethod
    def distance(self, other: 'Point') -> numbers.Real:
        """Calculate distance between two points"""

    @abc.abstractmethod
    def is_coincident(self,
                      other: 'Point',
                      tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                      ) -> bool:
        """Check whether two points are coincident"""


class Segment(Shape, metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def start_pt(self) -> Point:
        """Return start point"""

    @property
    @abc.abstractmethod
    def end_pt(self) -> Point:
        """Return end point"""

    @property
    @abc.abstractmethod
    def length(self) -> numbers.Real:
        """Returns horizontal length of the segment."""

    @abc.abstractmethod
    def as_segment2d(self):
        """Return a Segment2D instance"""

    def area(self):
        """Calculate shape area"""

        return 0.0


class Line(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0

    @abc.abstractmethod
    def pt(self,
           ndx: numbers.Integral) -> Optional[Point]:
        """Extract point at given index"""

    @abc.abstractmethod
    def num_pts(self) -> numbers.Integral:
        """
        Returns the number of line points.
        """

    @abc.abstractmethod
    def as_segment(self) -> Segment:
        """Return the segment defined by line start and end points"""

    @abc.abstractmethod
    def segments(self) -> List[Segment]:
        """Convert to a list of segments"""

    @abc.abstractmethod
    def add_pt(self, pt: Point):
        """
        In-place transformation of the original Line instance
        by adding a new point at the end."""

    @abc.abstractmethod
    def add_pts(self, pt_list: List[Point]):
        """
        In-place transformation of the original Line instance
        by adding a new set of points at the end."""

    @abc.abstractmethod
    def remove_coincident_points(self) -> Optional['Line']:
        """
        Remove coincident successive points.
        """

    @abc.abstractmethod
    def reversed(self) -> 'Line':
        """Reverse line"""

    '''
    @abc.abstractmethod
    def length(self) -> numbers.Real:
        """ The line length"""
    '''

    """
    def x_list(self) -> List[numbers.Real]:
        return list(map(lambda pt: pt.x, self))

    def y_list(self) -> List[numbers.Real]:
        return list(map(lambda pt: pt.y, self))

    def x_array(self):
        return np.asarray(self.x_list())

    def y_array(self):
        return np.asarray(self.y_list())

    def xy_arrays(self):
        return self.x_array, self.y_array

    def x_min(self):
        return np.nanmin(self.x_array())

    def x_max(self):
        return np.nanmax(self.x_array())

    def y_min(self):
        return np.nanmin(self.y_array())

    def y_max(self):
        return np.nanmax(self.y_array())
    """


class Polygon(Shape, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def num_side(self):
        """Return number of sides"""

    def area(self):

        return NotImplemented

    def length(self):

        return NotImplemented

    def clone(self):

        return NotImplemented

