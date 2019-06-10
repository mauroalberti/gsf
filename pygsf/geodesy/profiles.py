
from typing import List, Tuple

from math import asin, cos, pi

import numpy as np

from ..spatial.vectorial.geometries import Line

from .geodetic import epsg_4326_str


def profile_parameters(profile: Line) -> Tuple[List[float], List[float], List[float]]:
    """
    Calculates profile parameters for polar projections source datasets.

    :param profile: the profile line.
    :type profile: Line.
    :return: three profile parameters: horizontal distances, 3D distances, directional slopes
    :rtype: Tuple of three floats lists.
    """

    # calculate 3D distances between consecutive points

    if profile.crs == epsg_4326_str:

        # convert original values into ECEF values (x, y, z, time in ECEF global coordinate system)
        ecef_ln = profile.wgs842ecef()

        dist_3d_values = ecef_ln.step_lengths_3d()

    else:

        dist_3d_values = profile.step_lengths_3d()

    # calculate delta elevations between consecutive points

    delta_elev_values = profile.step_delta_z()

    # calculate slope along section

    dir_slopes_rads = []
    for delta_elev, dist_3D in zip(delta_elev_values, dist_3d_values):
        if dist_3D == 0.0:
            if delta_elev == 0.0:
                slope_rads = 0.0
            elif delta_elev < 0.0:
                slope_rads = - 0.5 * pi
            else:
                slope_rads = 0.5 * pi
        else:
            slope_rads = asin(delta_elev / dist_3D)

        dir_slopes_rads.append(slope_rads)

    # calculate horizontal distance along section

    horiz_dist_values = []
    for slope_rads, dist_3D in zip(dir_slopes_rads, dist_3d_values):

        horiz_dist_values.append(dist_3D * cos(slope_rads))

    return horiz_dist_values, dist_3d_values, dir_slopes_rads


