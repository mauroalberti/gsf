
import numbers

from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon

from .geometries import *


class MPolygon:
    """
    A shapely (multi)polygon with EPSG code.

    """

    def __init__(self,
                 shapely_geom: Union[Polygon, MultiPolygon],
                 epsg_code: numbers.Integral
                 ):
        """
        :param shapely_geom: the (multi)polygon
        :type shapely_geom: Union[Polygon, MultiPolygon]
        :param epsg_code: the EPSG code of the two geometries
        :type epsg_code: numbers.Integral
        """

        self._geom = shapely_geom
        self._epsg_code = epsg_code

    @property
    def geom(self):
        return self._geom

    @property
    def epsg_code(self):
        return self._epsg_code

    def intersect_line(self,
                       line: LineString,
                       ) -> Lines:
        """
        Determine the intersections between a mpolygon and a line.

        :param line: the line
        :type line: shapely.geometry.LineString
        :return: the intersecting lines
        :rtype: Lines
        """

        lines = Lines()

        intersections = line.intersection(self.geom)

        if intersections:

            if intersections.geom_type == "LineString":

                inters_ln = line_from_shapely(
                    shapely_geom=intersections,
                    epsg_code=self.epsg_code
                )

                lines.append(inters_ln)

            elif intersections.geom_type == "MultiLineString":

                for intersection_line in intersections:

                    inters_ln = line_from_shapely(
                        shapely_geom=intersection_line,
                        epsg_code=self.epsg_code
                    )

                    lines.append(inters_ln)

            else:

                pass

        return lines

