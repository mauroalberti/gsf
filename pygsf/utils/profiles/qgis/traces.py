import numbers
from typing import List, Tuple

from qgis._core import QgsRasterLayer

from pygsf.geometries.shapes.space2d import Line2D
from pygsf.geometries.shapes.space3d import Line3D, Point3D
from pygsf.geometries.shapes.space4d import Line4D
from pygsf.utils.vectors.qgis.lines import project_line2d
from pygsf.utils.qgis_utils.project import projectCrs
from pygsf.utils.rasters.qgis import QGisRasterParameters, interpolate_z


def line3d_from_dem(
    src_line2d: Line2D,
    qgs_raster_layer: QgsRasterLayer,
    qgis_raster_parameters: QGisRasterParameters
) -> Line3D:

    project_crs = projectCrs()
    if qgs_raster_layer.crs() != projectCrs():
        line2d_in_dem_crs = project_line2d(
            src_line2d=src_line2d,
            src_crs=project_crs,
            dest_crs=qgs_raster_layer.crs()
        )
    else:
        line2d_in_dem_crs = src_line2d

    line_3d = Line3D()

    for point2d_dem_crs, point2d_project_crs in zip(line2d_in_dem_crs.pts(), src_line2d.pts()):

        interpolated_z = interpolate_z(
            qgs_raster_layer=qgs_raster_layer,
            qgis_raster_parameters=qgis_raster_parameters,
            point2d=point2d_dem_crs
        )

        pt3d = Point3D(
            x=point2d_project_crs.x,
            y=point2d_project_crs.y,
            z=interpolated_z
        )

        line_3d.add_pt(pt3d)

    return line_3d


def lines3dnamed_from_dems(
        source_profile_line: Line2D,
        sample_distance: numbers.Real,
        selected_dems: List,
        selected_dem_parameters
) -> List[Tuple[str, Line3D]]:

    if isinstance(source_profile_line, Line4D):
        template_line2d = source_profile_line.as_line2d()
    elif isinstance(source_profile_line, Line3D):
        raise Exception("Source profile line is a Line3D instance")
    else:
        template_line2d = source_profile_line

    resampled_line = template_line2d.densify_2d_line(sample_distance)  # line resampled by sample distance

    # calculate 3D profiles from DEMs

    lines3d = []

    for dem, dem_params in zip(selected_dems, selected_dem_parameters):

        line3d = line3d_from_dem(
            resampled_line,
            dem,
            dem_params
        )

        lines3d.append((dem.name(), line3d))

    return lines3d