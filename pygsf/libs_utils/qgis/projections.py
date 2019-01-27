
from qgis.core import QgsCoordinateReferenceSystem

from .pygsf.spatial.vectorial.vectorial import Point, Segment, Line, MultiLine
from .pygsf.libs_utils.qgis.qgs_tools import qgs_project_xy


def line_project(line: Line, srcCrs: QgsCoordinateReferenceSystem, destCrs: QgsCoordinateReferenceSystem) -> Line:
    """
    Projects a line from a source to a destination CRS.

    :param line: the original line, to be projected.
    :type line: Line.
    :param srcCrs: the CRS of the original line.
    :type srcCrs: QgsCoordinateReferenceSystem.
    :param destCrs: the final CRS of the line.
    :type destCrs: QgsCoordinateReferenceSystem.
    :return: the projected line.
    :rtype: Line.
    """

    points = []
    for point in line.pts:
        x, y, z = point.toXYZ()
        x, y = qgs_project_xy(
            x=x,
            y=y,
            srcCrs=srcCrs,
            destCrs=destCrs)
        points.append(Point(x, y, z))

    return Line(points)


def multiline_project(multiline: MultiLine, srcCrs: QgsCoordinateReferenceSystem, destCrs: QgsCoordinateReferenceSystem) -> MultiLine:
    """
    Projects a multiline from a source to a destination CRS.

    :param multiline: the original multiline, to be projected.
    :type multiline: MultiLine.
    :param srcCrs: the CRS of the original multiline.
    :type srcCrs: QgsCoordinateReferenceSystem.
    :param destCrs: the final CRS of the multiline.
    :type destCrs: QgsCoordinateReferenceSystem.
    :return: the projected multiline.
    :rtype: MultiLine.
    """

    lines = []
    for line in multiline.lines:
        lines.append(multiline_project(line, srcCrs, destCrs))

    return MultiLine(lines)


def calculate_azimuth_correction(src_pt: Point, crs: QgsCoordinateReferenceSystem) -> float:
    """
    Calculates the empirical azimuth correction (angle between y-axis direction and geographic North)
    for a given point.

    :param src_pt: the point for which to calculate the correction.
    :type src_pt: Point.
    :param crs: the considered coordinate reference system.
    :type crs: QgsCoordinateReferenceSystem.
    :return: the azimuth angle.
    :rtype: float.
    """

    # Calculates dip direction correction with respect to project CRS y-axis orientation

    srcpt_prjcrs_x = src_pt.x
    srcpt_prjcrs_y = src_pt.y

    srcpt_epsg4326_lon, srcpt_epsg4326_lat = qgs_project_xy(
        x=srcpt_prjcrs_x,
        y=srcpt_prjcrs_y,
        srcCrs=crs)

    north_dummpy_pt_lon = srcpt_epsg4326_lon  # no change
    north_dummpy_pt_lat = srcpt_epsg4326_lat + (1.0 / 1200.0)  # add 3 minute-seconds (approximately 90 meters)

    dummypt_prjcrs_x, dummypt_prjcrs_y = qgs_project_xy(
        x=north_dummpy_pt_lon,
        y=north_dummpy_pt_lat,
        destCrs=crs)

    start_pt = Point(
        srcpt_prjcrs_x,
        srcpt_prjcrs_y)

    end_pt = Point(
        dummypt_prjcrs_x,
        dummypt_prjcrs_y)

    north_vector = Segment(
        start_pt=start_pt,
        end_pt=end_pt).vector()

    azimuth_correction = north_vector.azimuth

    return azimuth_correction
