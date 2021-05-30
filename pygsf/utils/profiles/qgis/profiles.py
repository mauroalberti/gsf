from math import ceil
from typing import Union, List, Tuple

from pygsf.geometries.shapes.collections import NamedLines
from pygsf.geometries.shapes.space2d import Line2D
from pygsf.geometries.shapes.space3d import Line3D
from pygsf.io.vectors.gps import TrackSource, GPXElevationUsage
from pygsf.utils.profiles.qgis.traces import lines3dnamed_from_dems
from pygsf.utils.rasters.qgis import get_min_dem_resolution


def try_prepare_grids_profile(
    profile_line: Union[Line2D, Line3D],
    track_source: TrackSource,
    gpx_elevation_usage: GPXElevationUsage,
    selected_grids: List,
    selected_grids_parameters: List,
    gpx_track_name: str
    ) -> Tuple[bool, Union[str, NamedLines]]:

    pt_num_threshold = 1e4

    try:

        # get DEMs resolutions in project CRS and choose the min value

        if track_source != TrackSource.GPX_FILE or \
           gpx_elevation_usage != GPXElevationUsage.ONLY_TO_USE:

            sample_distance = get_min_dem_resolution(
                selected_grids,
                selected_grids_parameters
            )

            # check total number of points in line(s) to create

            estimated_total_num_pts = 0

            profile_length = profile_line.length_2d()
            profile_num_pts = profile_length / sample_distance
            estimated_total_num_pts += profile_num_pts

            estimated_total_num_pts = int(ceil(estimated_total_num_pts))

            if estimated_total_num_pts > pt_num_threshold:
                return False, f"There are {estimated_total_num_pts} estimated points (limit is {pt_num_threshold}) in profile(s) to create.\nTry increasing sample distance value"

            dem_named_3dlines = lines3dnamed_from_dems(
                source_profile_line=profile_line,
                sample_distance=sample_distance,
                selected_dems=selected_grids,
                selected_dem_parameters=selected_grids_parameters
            )

            if dem_named_3dlines is None:
                return False, "Debug: profile not created"

        else:

            dem_named_3dlines = None

        # Manage GPX data

        if track_source == TrackSource.GPX_FILE and \
           gpx_elevation_usage != GPXElevationUsage.NOT_USED:

            gpx_named_3dlines = [(gpx_track_name, profile_line)]

        else:

            gpx_named_3dlines = None

        if dem_named_3dlines is None and gpx_named_3dlines is None:
            named_3dlines = None
        elif dem_named_3dlines is None:
            named_3dlines = gpx_named_3dlines
        elif gpx_named_3dlines is None:
            named_3dlines = dem_named_3dlines
        else:
            named_3dlines = dem_named_3dlines + gpx_named_3dlines

        if named_3dlines is None:
            return False, "Unable to create profiles"

        grids_profile = named_3dlines

        return True, grids_profile

    except Exception as e:

        return False, str(e)