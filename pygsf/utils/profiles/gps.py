from math import degrees, asin, cos, radians
from xml.dom.minidom import parse

from pygsf.profiles.sets import ProfileElevations
from pygsf.utils.vectors.qgis.points import TrackPointGPX


def topoprofiles_from_gpxfile(
        source_gpx_path: str,
        invert_profile: bool
) -> ProfileElevations:

    doc = parse(source_gpx_path)

    # define track name
    try:
        trkname = doc.getElementsByTagName('trk')[0].getElementsByTagName('name')[0].firstChild.data
    except:
        trkname = ''

    # get raw track point values (lat, lon, elev, time)
    track_raw_data = []
    for trk_node in doc.getElementsByTagName('trk'):
        for trksegment in trk_node.getElementsByTagName('trkseg'):
            for tkr_pt in trksegment.getElementsByTagName('trkpt'):
                track_raw_data.append((tkr_pt.getAttribute("lat"),
                                       tkr_pt.getAttribute("lon"),
                                       tkr_pt.getElementsByTagName("ele")[0].childNodes[0].data,
                                       tkr_pt.getElementsByTagName("time")[0].childNodes[0].data))

    # reverse profile orientation if requested
    if invert_profile:
        track_data = track_raw_data[::-1]
    else:
        track_data = track_raw_data

    # create list of TrackPointGPX elements
    track_points = []
    for val in track_data:
        gpx_trackpoint = TrackPointGPX(*val)
        track_points.append(gpx_trackpoint)

    # check for the presence of track points
    if len(track_points) == 0:
        raise Exception("No track point found in this file")

    # calculate delta elevations between consecutive points
    delta_elev_values = [np.nan]
    for ndx in range(1, len(track_points)):
        delta_elev_values.append(track_points[ndx].elev - track_points[ndx - 1].elev)

    # convert original values into ECEF values (x, y, z in ECEF global coordinate system)
    trk_ECEFpoints = [trck.as_pt4d() for trck in track_points]

    # calculate 3D distances between consecutive points
    dist_3D_values = [np.nan]
    for ndx in range(1, len(trk_ECEFpoints)):
        dist_3D_values.append(trk_ECEFpoints[ndx].distance(trk_ECEFpoints[ndx - 1]))

    # calculate slope along track
    dir_slopes = []
    for delta_elev, dist_3D in zip(delta_elev_values, dist_3D_values):
        try:
            slope = degrees(asin(delta_elev / dist_3D))
        except:
            slope = 0.0
        dir_slopes.append(slope)

    # calculate horizontal distance along track
    horiz_dist_values = []
    for slope, dist_3D in zip(dir_slopes, dist_3D_values):
        try:
            horiz_dist_values.append(dist_3D * cos(radians(slope)))
        except:
            horiz_dist_values.append(np.nan)

    # defines the cumulative 2D distance values
    cum_distances_2D = [0.0]
    for ndx in range(1, len(horiz_dist_values)):
        cum_distances_2D.append(cum_distances_2D[-1] + horiz_dist_values[ndx])

    # defines the cumulative 3D distance values
    cum_distances_3D = [0.0]
    for ndx in range(1, len(dist_3D_values)):
        cum_distances_3D.append(cum_distances_3D[-1] + dist_3D_values[ndx])

    lat_values = [track.lat for track in track_points]
    lon_values = [track.lon for track in track_points]
    time_values = [track.time for track in track_points]
    elevations = [track.elev for track in track_points]

    topo_profiles = ProfileElevations()

    #topo_profiles.line_source = gpx_source
    topo_profiles.inverted = invert_profile

    topo_profiles.lons = np.asarray(lon_values)
    topo_profiles.lats = np.asarray(lat_values)
    topo_profiles.times = time_values
    topo_profiles.surface_names = [trkname]  # [] required for compatibility with DEM case
    topo_profiles.incr_len_2d = np.asarray(cum_distances_2D)
    topo_profiles.incr_len_3d = [np.asarray(cum_distances_3D)]  # [] required for compatibility with DEM case
    topo_profiles.z_array = [np.asarray(elevations)]  # [] required for compatibility with DEM case
    topo_profiles.dir_slopes = [np.asarray(dir_slopes)]  # [] required for compatibility with DEM case

    return topo_profiles