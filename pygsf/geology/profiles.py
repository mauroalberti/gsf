# -*- coding: utf-8 -*-


from math import asin

from ..spatial.vectorial.vectorial import *
from ..spatial.geodetic import *
from ..orientations.orientations import Axis


def profile_params(profile: Line) -> Tuple[List[float], List[float], List[float]]:
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
        try:
            slope_rads = asin(delta_elev / dist_3D)
        except:
            slope_rads = 0.0
        dir_slopes_rads.append(slope_rads)

    # calculate horizontal distance along section

    horiz_dist_values = []
    for slope_rads, dist_3D in zip(dir_slopes_rads, dist_3d_values):
        try:
            horiz_dist_values.append(dist_3D * cos(slope_rads))
        except:
            horiz_dist_values.append(np.nan)

    return horiz_dist_values, dist_3d_values, dir_slopes_rads


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



def calculate_distance_with_sign(projected_point, section_init_pt, section_vector):
    """

    :param projected_point:
    :param section_init_pt:
    :param section_vector:
    :return:
    """

    assert projected_point.z != np.nan
    assert projected_point.z is not None

    projected_vector = Segment(section_init_pt, projected_point).vector()
    cos_alpha = section_vector.cos_angle(projected_vector)

    return projected_vector.len3D * cos_alpha


def get_intersection_slope(intersection_versor_3d, section_vector):
    """

    :param intersection_versor_3d:
    :param section_vector:
    :return:
    """

    slope_radians = abs(radians(intersection_versor_3d.slope))
    scalar_product_for_downward_sense = section_vector.sp(intersection_versor_3d.downward)
    if scalar_product_for_downward_sense > 0.0:
        intersection_downward_sense = "right"
    elif scalar_product_for_downward_sense == 0.0:
        intersection_downward_sense = "vertical"
    else:
        intersection_downward_sense = "left"

    return slope_radians, intersection_downward_sense


def calculate_intersection_versor(section_cartes_plane, structural_cartes_plane):
    """

    :param section_cartes_plane:
    :param structural_cartes_plane:
    :return:
    """

    return section_cartes_plane.inters_versor(structural_cartes_plane)


def calculate_nearest_intersection(intersection_versor_3d, section_cartes_plane, structural_cartes_plane,
                                   structural_pt):
    """

    :param intersection_versor_3d:
    :param section_cartes_plane:
    :param structural_cartes_plane:
    :param structural_pt:
    :return:
    """

    dummy_inters_point = section_cartes_plane.inters_point(structural_cartes_plane)
    dummy_structural_vector = Segment(dummy_inters_point, structural_pt).vector()
    dummy_distance = dummy_structural_vector.vDot(intersection_versor_3d)
    offset_vector = intersection_versor_3d.scale(dummy_distance)

    return Point(dummy_inters_point.x + offset_vector.x,
                 dummy_inters_point.y + offset_vector.y,
                 dummy_inters_point.z + offset_vector.z)


def calculate_axis_intersection(map_axis, section_cartes_plane, structural_pt):
    """

    :param map_axis:
    :param section_cartes_plane:
    :param structural_pt:
    :return:
    """

    axis_versor = map_axis.as_vect().versor
    l, m, n = axis_versor.x, axis_versor.y, axis_versor.z
    axis_param_line = ParamLine3D(structural_pt, l, m, n)
    return axis_param_line.intersect_cartes_plane(section_cartes_plane)


def map_measure_to_section(structural_rec, section_data, map_axis=None):
    """

    :param structural_rec:
    :param section_data:
    :param map_axis:
    :return:
    """

    # extract source data
    structural_pt, structural_plane, structural_pt_id = structural_rec
    section_init_pt, section_cartes_plane, section_vector = section_data['init_pt'], section_data['cartes_plane'], \
                                                            section_data['vector']

    # transform geological plane attitude into Cartesian plane
    structural_cartes_plane = structural_plane.plane(structural_pt)

    ## intersection versor
    intersection_versor_3d = calculate_intersection_versor(section_cartes_plane, structural_cartes_plane)

    # calculate slope of geological plane onto section plane
    slope_radians, intersection_downward_sense = get_intersection_slope(intersection_versor_3d, section_vector)

    # intersection point
    if map_axis is None:
        intersection_point_3d = calculate_nearest_intersection(intersection_versor_3d, section_cartes_plane,
                                                               structural_cartes_plane, structural_pt)
    else:
        intersection_point_3d = calculate_axis_intersection(map_axis, section_cartes_plane, structural_pt)

    # horizontal spat_distance between projected structural point and profile start
    signed_distance_from_section_start = calculate_distance_with_sign(intersection_point_3d, section_init_pt,
                                                                      section_vector)

    # solution for current structural point
    return PlaneAttitude(structural_pt_id,
                         structural_pt,
                         structural_plane,
                         intersection_point_3d,
                         slope_radians,
                         intersection_downward_sense,
                         signed_distance_from_section_start)


def map_struct_pts_on_section(structural_data, section_data, mapping_method):
    """
    defines:
        - 2D x-y location in section
        - plane-plane segment intersection

    :param structural_data:
    :param section_data:
    :param mapping_method:
    :return:
    """


    if mapping_method['method'] == 'nearest':
        return [map_measure_to_section(structural_rec, section_data) for structural_rec in structural_data]

    if mapping_method['method'] == 'common axis':
        map_axis = Axis(mapping_method['trend'], mapping_method['plunge'])
        return [map_measure_to_section(structural_rec, section_data, map_axis) for structural_rec in structural_data]

    if mapping_method['method'] == 'individual axes':
        assert len(mapping_method['individual_axes_values']) == len(structural_data)
        result = []
        for structural_rec, (trend, plunge) in zip(structural_data, mapping_method['individual_axes_values']):
            try:
                map_axis = Axis(trend, plunge)
                result.append(map_measure_to_section(structural_rec, section_data, map_axis))
            except:
                continue
        return result

