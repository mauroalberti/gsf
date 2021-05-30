from enum import Enum, auto
from typing import Tuple, Union
from xml.dom.minidom import parse

from pygsf.geometries.shapes.space4d import Line4D
from pygsf.utils.vectors.qgis.points import TrackPointGPX
from pygsf.utils.qgis_utils.project import projectCrs
from pygsf.utils.time import standard_gpstime_to_datetime


def try_extract_track_from_gpxfile(
        source_gpx_path: str,
        invert_profile: bool
) -> Tuple[bool, Union[str, Tuple[str, Line4D]]]:

    try:

        doc = parse(source_gpx_path)

        # define track name
        try:
            profile_name = doc.getElementsByTagName('trk')[0].getElementsByTagName('name')[0].firstChild.data
        except:
            profile_name = ''

        # get raw track point values (lat, lon, elev, time)
        track_raw_data = []
        for trk_node in doc.getElementsByTagName('trk'):
            for trksegment in trk_node.getElementsByTagName('trkseg'):
                for tkr_pt in trksegment.getElementsByTagName('trkpt'):
                    track_raw_data.append((tkr_pt.getAttribute("lat"),
                                           tkr_pt.getAttribute("lon"),
                                           tkr_pt.getElementsByTagName("ele")[0].childNodes[0].data,
                                           tkr_pt.getElementsByTagName("time")[0].childNodes[0].data))

        # create list of TrackPointGPX elements
        track_points = []
        for val in track_raw_data:
            lat, lon, ele, time = val
            time = standard_gpstime_to_datetime(time)
            gpx_trackpoint = TrackPointGPX(
                lat,
                lon,
                ele,
                time
            )
            track_points.append(gpx_trackpoint)

        # check for the presence of track points
        if len(track_points) == 0:
            return False, "No track point found in this file"

        # project track points to QGIS project CRS

        projected_pts = []
        for track_pt in track_points:

            projected_pt = track_pt.project(
                    dest_crs=projectCrs()
            )
            projected_pts.append(projected_pt)

        profile_line = Line4D(pts=projected_pts)

        if invert_profile:

            profile_line = profile_line.invert_direction()

        return True, (profile_name, profile_line)

    except Exception as e:

        return False, str(e)


class TrackSource(Enum):
    """
    The profile source type.
    """

    UNDEFINED  = auto()
    LINE_LAYER = auto()
    DIGITATION = auto()
    POINT_LIST = auto()
    GPX_FILE   = auto()


class GPXElevationUsage(Enum):
    """
    The profile source type.
    """

    NOT_USED = auto()
    USE_WITH_DEMS = auto()
    ONLY_TO_USE = auto()