# -*- coding: utf-8 -*-

from typing import Optional

from ...mathematics.defaults import *

from ...utils.time import *


WGS84 = {'semi-major axis': 6378137.0,
         'first eccentricity squared': 6.69437999014e-3}


def n_phi(phi_rad):
    """

    :param phi_rad:
    :return:
    """

    a = WGS84['semi-major axis']
    e_squared = WGS84['first eccentricity squared']
    return a / sqrt(1.0 - e_squared * sin(phi_rad) ** 2)


def geodetic2ecef(lat, lon, height):
    """

    :param lat:
    :param lon:
    :param height:
    :return:
    """

    e_squared = WGS84['first eccentricity squared']

    lat_rad, lon_rad = radians(lat), radians(lon)

    nphi = n_phi(lat_rad)

    x = (nphi + height) * cos(lat_rad) * cos(lon_rad)
    y = (nphi + height) * cos(lat_rad) * sin(lon_rad)
    z = (nphi * (1 - e_squared) + height) * sin(lat_rad)

    return x, y, z


class PolarSTPoint(object):

    def __init__(self, lat: float, lon: float, elev: float, time: Optional[float], crs_id: str):
        """
        Creates a space-time point expressed with polar (latitude-longitude) coordinates.

        :param lat: latitude.
        :type lat: float.
        :param lon: longitude.
        :type lon: float.
        :param elev: elevation.
        :type elev: float.
        :param time: time.
        :type time: optional float.
        :param crs_id: CRS code id.
        :type crs_id: basestring.
        """

        self.lat = lat
        self.lon = lon
        self.elev = elev
        self.time = time
        self.crs_id = crs_id

    def ecef_stpt(self):
        """
        Converts to ECEF (Earth-centered, Earth-fixed) Cartesian reference frame.

        :return:
        """

        x, y, _ = geodetic2ecef(self.lat, self.lon, self.elev)
        t = standard_gpstime_to_seconds(self.time)

        return x, y, self.elev, t
