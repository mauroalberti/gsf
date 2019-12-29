
from typing import Set

import geopandas as gpd


def geodataframe_geom_types(
    geodataframe: gpd.GeoDataFrame
) -> Set[str]:
    # Side effects: none
    """
    Return a set storing the geometric types in a GeoDataFrame instance.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: the set of geometric types values
    :rtype: Set[str]
    """

    return set(geodataframe.geom_type)


def containsLines(
        geodataframe: gpd.GeoDataFrame
) -> bool:
    # Side effects: none
    """
    Check if a GeoDataFrame instance contains lines.

    :param geodataframe: the input geodataframe
    :type geodataframe: gpd.GeoDataFrame
    :return: if a GeoDataFrame instance contains lines
    :rtype: bool
    """

    return 'LineString' in geodataframe_geom_types(geodataframe)

