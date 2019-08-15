# -*- coding: utf-8 -*-


from typing import Optional

from .rasters.geoarray import GeoArray
from .vectorial.geometries import Line
from pygsf.spatial.vectorial.geometries import Point


def line_on_grid(ga: GeoArray, profile_line: Line) -> Optional[Line]:
    """
    Calculates a line draped on a grid.

    :param ga: geoarray
    :type ga: GeoArray.
    :param profile_line: the profile line.
    :type profile_line: Line
    :return: the profile.
    :rtype: Optional[Line].
    """

    if ga.crs() != profile_line.crs():
        return None

    epsg_line = profile_line.epsg()

    lnProfile = Line(epsg_cd=epsg_line)

    for point in profile_line.pts():

        z = ga.interpolate_bilinear(point.x, point.y)
        if z:
            lnProfile.add_pt(
                Point(
                    x=point.x,
                    y=point.y,
                    z=z,
                    epsg_cd=epsg_line)
            )

    return lnProfile
