
from typing import List, Tuple, Optional, Union

from geopandas import GeoDataFrame

from pygsf.spatial.vectorial.geometries import Plane
from pygsf.spatial.vectorial.geodataframes import *

from .base import GeorefAttitude


def try_extract_georeferenced_attitudes(
        geodataframe: GeoDataFrame,
        azim_fldnm: str,
        dip_ang_fldnm: str,
        id_fldnm: Optional[str] = None,
        is_rhrstrike: bool = False
) -> Tuple[bool, Union[str, List[GeorefAttitude]]]:
    """
    Try extracting the georeferenced _attitudes from a geopandas GeoDataFrame instance representing point records.

    :param geodataframe: the source geodataframe.
    :type geodataframe: GeoDataFrame.
    :param azim_fldnm: the name of the azimuth field in the geodataframe.
    :type azim_fldnm: basestring.
    :param dip_ang_fldnm: the name of the dip angle field in the geodataframe.
    :type dip_ang_fldnm: basestring.
    :param is_rhrstrike: whether the dip azimuth is strike RHR
    :type is_rhrstrike: bool.
    :return: the success status and an error message or a collection of georeferenced attitudes, one for each source record.
    :rtype: Tuple[bool, Union[str, List[GeorefAttitude]]]
    """

    try:

        epsg = get_epsg(geodataframe)

        attitudes = []

        for ndx, row in geodataframe.iterrows():

            pt = row['geometry']
            x, y = pt.x, pt.y

            if id_fldnm:
                azimuth, dip_ang, id = row[azim_fldnm], row[dip_ang_fldnm], row[id_fldnm]
            else:
                azimuth, dip_ang, id = row[azim_fldnm], row[dip_ang_fldnm], ndx + 1

            if is_rhrstrike:
                azimuth = (azimuth + 90.0) % 360.0

            attitudes.append(GeorefAttitude(id, Point(x, y, epsg_code=epsg), Plane(azimuth, dip_ang)))

        return True, attitudes

    except Exception as e:

        return False, str(e)


