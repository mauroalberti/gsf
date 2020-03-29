
import numbers

import shapely

from .geometries import *


def polygon_line_intersection(
        mpolygon: Union[shapely.geometry.Polygon, shapely.geometry.MultiPolygon],
        line: shapely.geometry.LineString,
        epsg_code: numbers.Integral
) -> Lines:
    """
    Determine the intersections between a polygon and a line.

    :param mpolygon: the (multi)polygon to intersect
    :type mpolygon: Union[shapely.geometry.Polygon, shapely.geometry.MultiPolygon]
    :param line: the line
    :type line: shapely.geometry.LineString
    :param epsg_code: the EPSG code of the two geometries
    :type epsg_code: numbers.Integral
    :return: the intersecting lines
    :rtype: Lines
    """

    lines = Lines()

    intersections = line.intersection(mpolygon)

    if intersections:

        if intersections.geom_type == "LineString":

            inters_ln = line_from_shapely(
                shapely_linestring=intersections,
                epsg_code=epsg_code
            )

            lines.append(inters_ln)

        elif intersections.geom_type == "MultiLineString":

            for intersection_line in intersections:
                inters_ln = line_from_shapely(
                    shapely_linestring=intersection_line,
                    epsg_code=epsg_code
                )

                lines.append(inters_ln)

        else:

            pass

    return lines
