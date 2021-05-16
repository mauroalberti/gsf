
from typing import List, Optional, Tuple
import numbers
import abc

import numpy as np


class Shape(object, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def area(self):
        """Calculate shape area"""

    @abc.abstractmethod
    def length(self):
        """Calculate shape length"""


class Point(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0

    def length(self):
        """Calculate shape length"""

        return 0.0

    @abc.abstractmethod
    def is_coincident(self, other: 'Point') -> bool:
        """Check whether two points are coincident"""

    @abc.abstractmethod
    def distance(self, other: 'Point') -> bool:
        """Calculate distaance between two points"""


class Segment(Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0


class Line(list, Shape, metaclass=abc.ABCMeta):

    def area(self):
        """Calculate shape area"""

        return 0.0

    def pt(self,
           ndx: numbers.Integral) -> Point:

        return self[ndx]

    def num_pts(self) -> numbers.Integral:
        return len(self)

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

    def start_pt(self) -> Optional[Point]:
        """
        Return the first point of a Line or None when no points.

        :return: the first point or None.
        """

        return self.pt(0) if len(self) > 0 else None

    def end_pt(self) -> Optional[Point]:
        """
        Return the last point of a Line or None when no points.

        :return: the last point or None.
        """

        return self.pt(-1) if len(self) > 0 else None

    @abc.abstractmethod
    def as_segment(self) -> Segment:
        """Return the segment defined by line start and end points"""

    @abc.abstractmethod
    def as_segments(self) -> List[Segment]:
        """Convert to a list of segments"""

    def xy_zipped(self) -> List[Tuple[numbers.Real, numbers.Real]]:

        return [(x, y) for x, y in zip(self.x_list(), self.y_list())]

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

    @abc.abstractmethod
    def clone(self) -> 'Line':
        """
        Clone a line."""

    @abc.abstractmethod
    def length(self) -> numbers.Real:
        """ The line length"""


class Polygon(Shape, metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def num_side(self):
        """Return number of sides"""


