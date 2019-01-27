
from ...geology.profile import TopoProfiles
from ...config.profiles import demline_source, gpxfile_source
from .qgs_tools import topoline_from_dem, get_project_crs


def profile_elevations(canvas, source_profile_line, sample_distance, selected_dems, selected_dem_parameters,
                       invert_profile) -> TopoProfiles:
    """

    :param canvas:
    :param source_profile_line:
    :param sample_distance:
    :param selected_dems:
    :param selected_dem_parameters:
    :param invert_profile:
    :return: TopoProfiles.
    """

    # get project CRS information

    project_crs = get_project_crs(canvas)

    if invert_profile:
        line = source_profile_line.reverse_direction()
    else:
        line = source_profile_line

    resampled_line = line.densify_2d_line(sample_distance)  # line resampled by sample distance

    # calculate 3D profiles from DEMs

    llnsTopoLines = []
    for dem, dem_params in zip(selected_dems, selected_dem_parameters):
        dem_topoline3d = topoline_from_dem(
            resampled_line,
            project_crs,
            dem,
            dem_params)
        llnsTopoLines.append(dem_topoline3d)

    # setup topoprofiles properties

    return TopoProfiles(
        crs_authid=crs_authid,
        profile_source=demline_source,
        source_names=[dem.name() for dem in selected_dems],
        xs=np.asarray(resampled_line.x_list),
        ys=np.asarray(resampled_line.y_list),
        zs=[np.asarray([line.z_array() for line in topo_lines])],
        times=None,
        inverted=invert_profile)

