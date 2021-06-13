
import numbers
import abc

from ...geometries.shapes.abstract import Shape
from ...utils.types import *


class GeoShape(object, metaclass=abc.ABCMeta):

    def __init__(self,
                 shape: Shape,
                 epsg_cd: numbers.Integral):

        check_type(shape, "Shape", Shape)
        check_type(epsg_cd, "EPSG code", numbers.Integral)

        self._shape = shape
        self._epsg_code = epsg_cd

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

'''
class GeoLine(GeoShape, metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def geometry(self) -> Line:
        """Return shape"""

    @geometry.setter
    @abc.abstractmethod
    def shape(self,
              shp: Line):
        """Set shape"""
'''