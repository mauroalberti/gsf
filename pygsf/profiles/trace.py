import numbers

from ..geometries.shapes.abstract import Line
from ..geometries.shapes.space2d import Line2D
from ..geometries.shapes.space3d import Line3D
from ..geometries.shapes.space4d import Line4D
from ..georeferenced.shapes.abstract import GeoLine


class Trace(GeoLine):
    """
    The profile trace.
    """

    def __init__(self,
                 line: Line,
                 epsg_code: numbers.Integral = -1
                 ):

        super(Trace, self).__init__(line, epsg_code)



