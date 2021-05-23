
from typing import List, Optional, Tuple
import numbers
import abc

import numpy as np

from ...geometries.shapes.abstract import Shape, Line


class GeoShape(object, metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def shape(self) -> Shape:
        """Return shape"""

    @shape.setter
    @abc.abstractmethod
    def shape(self,
                  shp: Shape):
        """Set shape"""

    @property
    @abc.abstractmethod
    def epsg_code(self):
        """Return EPSG code"""

    @epsg_code.setter
    @abc.abstractmethod
    def epsg_code(self,
                  epsg_cd: numbers.Integral):
        """Set EPSG code"""


class GeoLine(GeoShape, metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def shape(self) -> Line:
        """Return shape"""

    @shape.setter
    @abc.abstractmethod
    def shape(self,
                  shp: Line):
        """Set shape"""
