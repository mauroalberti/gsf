# -*- coding: utf-8 -*-


from math import asin
import copy
#import xml.dom.minidom

from .exceptions import GPXIOException

from ..mathematics.statistics import get_statistics
from ..geography.earth import projectionType
from ..spatial.vectorial.vectorial import *
from ..spatial.vectorial.geodetic import *

#from ..libs_utils.qgis.projections import multiline_project
#from ..libs_utils.qgis.qgs_tools import get_project_crs, line_geoms_attrs, get_zs_from_dem


class GeoProfilesSet(object):
    """
    Represents a set of Geoprofile elements,
    stored as a list
    """

    def __init__(self, name=""):
        """

        :param name:
        """

        self._name = name
        self._geoprofiles = []  # a list of ProfileElements instances
        self.profiles_created = False
        self.plot_params = None

    @property
    def name(self):
        """

        :return:
        """
        return self._name

    @name.setter
    def name(self, new_name):
        """

        :param new_name:
        :return:
        """

        self._name = new_name

    @property
    def geoprofiles(self):
        """

        :return:
        """

        return self._geoprofiles

    @property
    def geoprofiles_num(self):
        """

        :return:
        """

        return len(self._geoprofiles)

    def geoprofile(self, ndx):
        """

        :param ndx:
        :return:
        """

        return self._geoprofiles[ndx]

    def append(self, geoprofile):
        """

        :param geoprofile:
        :return:
        """

        self._geoprofiles.append(geoprofile)

    def insert(self, ndx, geoprofile):
        """

        :param ndx:
        :param geoprofile:
        :return:
        """

        self._geoprofiles.insert(ndx, geoprofile)

    def move(self, ndx_init, ndx_final):
        """

        :param ndx_init:
        :param ndx_final:
        :return:
        """

        geoprofile = self._geoprofiles.pop(ndx_init)
        self.insert(ndx_final, geoprofile)

    def move_up(self, ndx):
        """

        :param ndx:
        :return:
        """

        self.move(ndx, ndx -1)

    def move_down(self, ndx):
        """

        :param ndx:
        :return:
        """

        self.move(ndx, ndx + 1)

    def remove(self, ndx):
        """

        :param ndx:
        :return:
        """

        _ = self._geoprofiles.pop(ndx)


class GeoProfile(object):
    """
    Class representing the topographic and geological elements
    embodying a single geological profile.
    """

    def __init__(self):
        """

        """

        self.source_data_type = None
        self.original_line = None
        self.sample_distance = None  # max spacing along profile; float
        self.resampled_line = None

        self.topoprofiles = None  # instance of TopoProfiles
        self.geoplane_attitudes = []
        self.geosurfaces = []
        self.geosurfaces_ids = []
        self.lineaments = []
        self.outcrops = []

    def set_topo_profiles(self, topo_profiles):
        """

        :param topo_profiles:
        :return:
        """

        self.topoprofiles = topo_profiles

    def add_intersections_pts(self, intersection_list):
        """

        :param intersection_list:
        :return:
        """

        self.lineaments += intersection_list

    def add_intersections_lines(self, formation_list, intersection_line3d_list, intersection_polygon_s_list2):
        """

        :param formation_list:
        :param intersection_line3d_list:
        :param intersection_polygon_s_list2:
        :return:
        """

        self.outcrops = zip(formation_list, intersection_line3d_list, intersection_polygon_s_list2)

    def get_current_dem_names(self):
        """

        :return:
        """

        return self.topoprofiles.surface_names

    def max_s(self):
        """

        :return:
        """

        return self.topoprofiles.max_s()

    def min_z_topo(self):
        """

        :return:
        """

        return self.topoprofiles.min_z()

    def max_z_topo(self):
        """

        :return:
        """

        return self.topoprofiles.max_z()

    def min_z_plane_attitudes(self):
        """

        :return:
        """

        # TODO:  manage case for possible nan p_z values
        return min([plane_attitude.pt_3d.p_z for plane_attitude_set in self.geoplane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.max_s()])

    def max_z_plane_attitudes(self):
        """

        :return:
        """

        # TODO:  manage case for possible nan p_z values
        return max([plane_attitude.pt_3d.p_z for plane_attitude_set in self.geoplane_attitudes for plane_attitude in
                    plane_attitude_set if 0.0 <= plane_attitude.sign_hor_dist <= self.max_s()])

    def min_z_curves(self):
        """

        :return:
        """

        return min([pt_2d.p_y for multiline_2d_list in self.geosurfaces for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.max_s()])

    def max_z_curves(self):
        """

        :return:
        """

        return max([pt_2d.p_y for multiline_2d_list in self.geosurfaces for multiline_2d in multiline_2d_list for line_2d in
                    multiline_2d.lines for pt_2d in line_2d.pts if 0.0 <= pt_2d.p_x <= self.max_s()])

    def min_z(self):
        """

        :return:
        """

        min_z = self.min_z_topo()

        if len(self.geoplane_attitudes) > 0:
            min_z = min([min_z, self.min_z_plane_attitudes()])

        if len(self.geosurfaces) > 0:
            min_z = min([min_z, self.min_z_curves()])

        return min_z

    def max_z(self):
        """

        :return:
        """

        max_z = self.max_z_topo()

        if len(self.geoplane_attitudes) > 0:
            max_z = max([max_z, self.max_z_plane_attitudes()])

        if len(self.geosurfaces) > 0:
            max_z = max([max_z, self.max_z_curves()])

        return max_z

    def add_plane_attitudes(self, plane_attitudes):
        """

        :param plane_attitudes:
        :return:
        """

        self.geoplane_attitudes.append(plane_attitudes)

    def add_curves(self, lMultilines, lIds):
        """

        :param lMultilines:
        :param lIds:
        :return:
        """

        self.geosurfaces.append(lMultilines)
        self.geosurfaces_ids.append(lIds)


class TopoProfiles(object):
    """
    A set of topographic profiles created from a
    single planar profile trace and a set of source elevations
    data (for instance, a pair of dems).
    """

    def __init__(self,
                 crs_authid: str,
                 profile_source: str,
                 source_names: List[str],
                 xs: List[float],
                 ys: List[float],
                 zs: List[List[float]],
                 times: List[float],
                 inverted: bool):
        """
        Creates the topoprofile instance.

        :param crs_authid: the CRS authid. E.g.: "EPSG:4326".
        :type crs_authid: basestring.
        :param profile_source: the source type of the elevations, i.e., DEMs or GPS.
        :type profile_source: basestring.
        :param source_names: list of data sources names.
        :type source_names: list of basestrings.
        :param xs:
        :param ys:
        :param zs:
        :param times:
        :param inverted:
        """


        self.crs_authid = crs_authid
        self.profile_source = profile_source
        self.source_names = source_names
        self.xs = xs
        self.ys = ys
        self.zs = zs
        self.times = times
        self.inverted = inverted

        ###

        if projectionType(self.crs_authid) == "polar":
            self.profile_s, self.profiles_s3ds, self.profiles_dirslopes = self.params_polar()
        else:
            self.profile_s = np.asarray(source_trace.incremental_length_2d())

            self.profiles_dirslopes = [np.asarray(line.slopes()) for line in topo_lines]

            self.profiles_s3ds = [np.asarray(line.incremental_length_3d()) for line in topo_lines]

            """
            self.dem_params = [DEMParams(dem, params) for (dem, params) in
                               zip(source_dems, dem_parameters)]
            self.gpx_params = None                  
            """

    @property
    def num_pnts(self) -> int:
        """
        Returns the number of points in the profile.
        """

        return len(self.xs)

    def params_polar(self):
        """
        Calculates parameters for polar projections source datasets.
        """

        xs = self.xs.tolist()
        ys = self.ys.tolist()

        for profile_ndx in range(len(self.zs)):

            zs = self.zs[profile_ndx].tolist()
            times = self.times.tolist()

            # convert original values into ECEF values (x, y, z, time in ECEF global coordinate system)

            ecef_pt = [PolarSTPoint(lat, lon, elev, time).ecef_stpt() for lat, lon, elev, time in zip(xs, ys, zs, times)]

            # calculate 3D distances between consecutive points

            dist_3D_values = [np.nan]
            for ndx in range(1, self.num_pnts):
                dist_3D_values.append(ecef_pt[ndx].dist3DWith(ecef_pt[ndx - 1]))

            # calculate delta elevations between consecutive points

            delta_elev_values = [np.nan]
            for ndx in range(1, self.num_pnts):
                delta_elev_values.append(zs[ndx] - zs[ndx - 1])

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

            """
            # defines the cumulative 2D distance values
            
            cum_distances_2D = [0.0]
            for ndx in range(1, len(horiz_dist_values)):
                cum_distances_2D.append(cum_distances_2D[-1] + horiz_dist_values[ndx])
    
            # defines the cumulative 3D distance values
            
            cum_distances_3D = [0.0]
            for ndx in range(1, len(dist_3D_values)):
                cum_distances_3D.append(cum_distances_3D[-1] + dist_3D_values[ndx])
            """

        return horiz_dist_values, dist_3D_values, dir_slopes


    def max_s(self) -> float:
        """
        Return the horizontal length of the topographic profile trace.

        :return: the length of the profile.
        :rtype: float
        """

        return self.profile_s[-1]

    def min_z(self) -> float:
        """
        Returns the minimum elevation value in the
        set of stored topographic profiles.

        :return:
        """

        return float(min([np.nanmin(prof_elev) for prof_elev in self.profiles_elevs]))

    def max_z(self) -> float:
        """
        Return the maximum elevation in the set of stored topographic profiles.

        :return: the maximum elevation in the profiles.
        :rtype: float.
        """

        return float(max([np.nanmax(prof_elev) for prof_elev in self.profiles_elevs]))

    @property
    def absolute_slopes(self) -> List[np.ndarray]:
        """
        Returns a list of the absolute slopes of the stored topographic profiles.

        :return: a list of the absolute slopes of the stored topographic profiles.
        :rtype: list of nd.arrays.
        """

        return [np.fabs(prof_dirslopes) for prof_dirslopes in self.profiles_dirslopes]

    @property
    def statistics_elev(self) -> List(Dict):
        """
        Return statistics (min, max, mean, var, std) for elevations in stored topographic profiles.

        :return: statistics values.
        :rtype: List of dictionaries, one for each topographic profile.
        """

        return [get_statistics(prof_elev) for prof_elev in self.profiles_elevs]

    @property
    def statistics_dirslopes(self) -> List(Dict):
        """
        Return statistics (min, max, mean, var, std) for directional slopes of stored topographic profiles.

        :return: statistics values.
        :rtype: List of dictionaries, one for each topographic profile.
        """

        return [get_statistics(prof_dirslopes) for prof_dirslopes in self.profiles_dirslopes]

    @property
    def statistics_slopes(self) -> List(Dict):
        """
        Return statistics (min, max, mean, var, std) for absolute slopes of stored topographic profiles.

        :return: statistics values.
        :rtype: List of dictionaries, one for each topographic profile.
        """

        return [get_statistics(prof_abs_slopes) for prof_abs_slopes in self.absolute_slopes]


class PlaneAttitude(object):

    def __init__(self, rec_id, source_point_3d, source_geol_plane, point_3d, slope_rad, dwnwrd_sense, sign_hor_dist):
        """

        :param rec_id:
        :param source_point_3d:
        :param source_geol_plane:
        :param point_3d:
        :param slope_rad:
        :param dwnwrd_sense:
        :param sign_hor_dist:
        """

        self.id = rec_id
        self.src_pt_3d = source_point_3d
        self.src_geol_plane = source_geol_plane
        self.pt_3d = point_3d
        self.slope_rad = slope_rad
        self.dwnwrd_sense = dwnwrd_sense
        self.sign_hor_dist = sign_hor_dist


def topoprofiles_from_gpxsource(source_gpx_path, invert_profile, gpx_source):
    """

    :param source_gpx_path:
    :param invert_profile:
    :param gpx_source:
    :return:
    """

    doc = xml.dom.minidom.parse(source_gpx_path)

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

    # create list of PolarSTPoint elements
    track_points = []
    for val in track_data:
        gpx_trackpoint = PolarSTPoint(*val)
        track_points.append(gpx_trackpoint)

    # check for the presence of track points
    if len(track_points) == 0:
        raise GPXIOException("No track point found in this file")

    lats =       np.asarray([track.lat  for track in track_points])
    lons =       np.asarray([track.lon  for track in track_points])
    times =      np.asarray([track.time for track in track_points])
    elevations = np.asarray([track.elev for track in track_points])

    return TopoProfiles(
        crs_authid=4236,
        profile_source=gpx_source,
        xs=lons,
        ys=lats,
        times=times,
        source_names=[trkname],
        topo_lines=None,
        zs=[elevations],
        inverted=invert_profile)


    profile_s=np.asarray(cum_distances_2D),
    profiles_s3ds=[np.asarray(cum_distances_3D)]

    profiles_dirslopes=[np.asarray(dir_slopes)]



def calculate_profile_lines_intersection(multilines2d_list, id_list, profile_line2d):
    """

    :param multilines2d_list:
    :param id_list:
    :param profile_line2d:
    :return:
    """

    profile_segment2d_list = profile_line2d.as_segments()

    profile_segment2d = profile_segment2d_list[0]

    intersection_list = []
    for ndx, multiline2d in enumerate(multilines2d_list):
        if id_list is None:
            multiline_id = ''
        else:
            multiline_id = id_list[ndx]
        for line2d in multiline2d.lines:
            for line_segment2d in line2d.as_segments():
                try:
                    intersection_point2d = profile_segment2d.intersection_2d_pt(line_segment2d)
                except ZeroDivisionError:
                    continue
                if intersection_point2d is None:
                    continue
                if line_segment2d.contains_2d_pt(intersection_point2d) and \
                   profile_segment2d.contains_2d_pt(intersection_point2d):
                    intersection_list.append([intersection_point2d, multiline_id])

    return intersection_list


def intersection_distances_by_profile_start_list(profile_line, intersections):
    """

    :param profile_line:
    :param intersections:
    :return:
    """

    # convert the profile line from a CartesianLine2DT to a CartesianSegment2DT
    profile_segment2d_list = profile_line.as_segments()
    # debug
    assert len(profile_segment2d_list) == 1
    profile_segment2d = profile_segment2d_list[0]

    # determine distances for each point in intersection list
    # creating a list of float values
    distance_from_profile_start_list = []
    for intersection in intersections:
        distance_from_profile_start_list.append(profile_segment2d.start_pt.dist2DWith(intersection[0]))

    return distance_from_profile_start_list


def profile_polygon_intersection(profile_qgsgeometry, polygon_layer, inters_polygon_classifaction_field_ndx):
    """

    :param profile_qgsgeometry:
    :param polygon_layer:
    :param inters_polygon_classifaction_field_ndx:
    :return:
    """

    intersection_polyline_polygon_crs_list = []

    if polygon_layer.selectedFeatureCount() > 0:
        features = polygon_layer.selectedFeatures()
    else:
        features = polygon_layer.getFeatures()

    for polygon_feature in features:
        # retrieve every (selected) feature with its geometry and attributes

        # fetch geometry
        poly_geom = polygon_feature.geometry()

        intersection_qgsgeometry = poly_geom.intersection(profile_qgsgeometry)

        try:
            if intersection_qgsgeometry.isEmpty():
                continue
        except:
            try:
                if intersection_qgsgeometry.isGeosEmpty():
                    continue
            except:
                return False, "Missing function for checking empty geometries.\nPlease upgrade QGIS"

        if inters_polygon_classifaction_field_ndx >= 0:
            attrs = polygon_feature.attributes()
            polygon_classification = attrs[inters_polygon_classifaction_field_ndx]
        else:
            polygon_classification = None

        if intersection_qgsgeometry.isMultipart():
            lines = intersection_qgsgeometry.asMultiPolyline()
        else:
            lines = [intersection_qgsgeometry.asPolyline()]

        for line in lines:
            intersection_polyline_polygon_crs_list.append([polygon_classification, line])

    return True, intersection_polyline_polygon_crs_list


def extract_multiline2d_list(structural_line_layer, project_crs):
    """

    :param structural_line_layer:
    :param project_crs:
    :return:
    """

    line_orig_crs_geoms_attrs = line_geoms_attrs(structural_line_layer)

    line_orig_geom_list3 = [geom_data[0] for geom_data in line_orig_crs_geoms_attrs]
    line_orig_crs_MultiLine2D_list = [xytuple_l2_to_MultiLine(xy_list2) for xy_list2 in line_orig_geom_list3]
    line_orig_crs_clean_MultiLine2D_list = [multiline_2d.remove_coincident_points() for multiline_2d in
                                            line_orig_crs_MultiLine2D_list]

    # get CRS information
    structural_line_layer_crs = structural_line_layer.crs()

    # project input line layer to project CRS
    line_proj_crs_MultiLine2D_list = [multiline_project(multiline2d, structural_line_layer_crs, project_crs) for
                                          multiline2d in line_orig_crs_clean_MultiLine2D_list]

    return line_proj_crs_MultiLine2D_list


def define_plot_structural_segment(structural_attitude, profile_length, vertical_exaggeration, segment_scale_factor=70.0):
    """

    :param structural_attitude:
    :param profile_length:
    :param vertical_exaggeration:
    :param segment_scale_factor:
    :return:
    """

    ve = float(vertical_exaggeration)
    intersection_point = structural_attitude.pt_3d
    z0 = intersection_point.z

    h_dist = structural_attitude.sign_hor_dist
    slope_rad = structural_attitude.slope_rad
    intersection_downward_sense = structural_attitude.dwnwrd_sense
    length = profile_length / segment_scale_factor

    s_slope = sin(float(slope_rad))
    c_slope = cos(float(slope_rad))

    if c_slope == 0.0:
        height_corr = length / ve
        structural_segment_s = [h_dist, h_dist]
        structural_segment_z = [z0 + height_corr, z0 - height_corr]
    else:
        t_slope = s_slope / c_slope
        width = length * c_slope

        length_exag = width * sqrt(1 + ve*ve * t_slope*t_slope)

        corr_width = width * length / length_exag
        corr_height = corr_width * t_slope

        structural_segment_s = [h_dist - corr_width, h_dist + corr_width]
        structural_segment_z = [z0 + corr_height, z0 - corr_height]

        if intersection_downward_sense == "left":
            structural_segment_z = [z0 - corr_height, z0 + corr_height]

    return structural_segment_s, structural_segment_z


def calculate_projected_3d_pts(canvas, struct_pts, structural_pts_crs, demObj):
    """

    :param canvas:
    :param struct_pts:
    :param structural_pts_crs:
    :param demObj:
    :return:
    """

    demCrs = demObj.params.crs

    # check if on-the-fly-projection is set on
    project_crs = get_project_crs(canvas)

    # set points in the project crs
    if structural_pts_crs != project_crs:
        struct_pts_in_prj_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, project_crs)
    else:
        struct_pts_in_prj_crs = copy.deepcopy(struct_pts)

        # project the source points from point layer crs to DEM crs
    # if the two crs are different
    if structural_pts_crs != demCrs:
        struct_pts_in_dem_crs = calculate_pts_in_projection(struct_pts, structural_pts_crs, demCrs)
    else:
        struct_pts_in_dem_crs = copy.deepcopy(struct_pts)

    # - 3D structural points, with x, y, and z extracted from the current DEM
    struct_pts_z = get_zs_from_dem(struct_pts_in_dem_crs, demObj)

    assert len(struct_pts_in_prj_crs) == len(struct_pts_z)

    return [Point(pt.x, pt.y, z) for (pt, z) in zip(struct_pts_in_prj_crs, struct_pts_z)]

