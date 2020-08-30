
import numbers
from typing import Optional, List, Union, Tuple

import numpy
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon

from pygsf.spatial.space3d.vectorial.geometries import *

from pygsf.crs.crs import Crs, check_epsg, check_crs
from pygsf.geometries.geom3d.shapes import Point, Segment, Points, analizeJoins
from pygsf.utils.types import check_type


class MGeoPolygon:
    """
    A shapely (multi)polygon with EPSG code.

    """

    def __init__(self,
                 shapely_geom: Union[Polygon, MultiPolygon],
                 epsg_code: numbers.Integral
                 ):
        """
        :param shapely_geom: the (multi)polygon
        :type shapely_geom: Union[Polygon, MultiPolygon]
        :param epsg_code: the EPSG code of the two geometries
        :type epsg_code: numbers.Integral
        """

        check_type(
            shapely_geom,
            "Polygon",
            (Polygon, MultiPolygon)
        )

        self._geom = shapely_geom
        self._epsg_code = epsg_code

    @property
    def geom(self):
        return self._geom

    @property
    def epsg_code(self):
        return self._epsg_code

    def intersect_line(self,
                       line: LineString,
                       ) -> Lines:
        """
        Determine the intersections between a mpolygon and a line.

        :param line: the line
        :type line: shapely.geometry.LineString
        :return: the intersecting lines
        :rtype: Lines
        """

        lines = Lines()

        intersections = line.intersection(self.geom)

        if intersections:

            if intersections.geom_type == "LineString":

                inters_ln = line_from_shapely(
                    shapely_geom=intersections,
                    epsg_code=self.epsg_code
                )

                lines.append(inters_ln)

            elif intersections.geom_type == "MultiLineString":

                for intersection_line in intersections:

                    inters_ln = line_from_shapely(
                        shapely_geom=intersection_line,
                        epsg_code=self.epsg_code
                    )

                    lines.append(inters_ln)

            else:

                pass

        return lines


class Points:
    """
    Collection of points.
    """

    def __init__(self,
                 epsg_code: numbers.Integral,
                 x_array: np.ndarray,
                 y_array: np.ndarray,
                 z_array: Optional[np.ndarray] = None
                 #t_array: Optional[np.ndarray] = None
                 ):
        """
        Construct a point list from a set of array values and an EPSG code.

        :param epsg_code: the EPSG code of the points
        :type epsg_code: numbers.Integral
        :param x_array: the array storing the x values
        :type x_array: np.ndarray
        :param y_array: the array storing the y values
        :type y_array: np.ndarray
        :param z_array: the optional array storing the z values
        :type z_array: np.ndarray
        :param t_array: the optional array storing the t values
        :type t_array: np.ndarray
        """

        check_type(
            var=epsg_code,
            name="EPSG code",
            expected_types=numbers.Integral
        )

        check_type(
            var=x_array,
            name="X array",
            expected_types=np.ndarray
        )

        check_type(
            var=y_array,
            name="Y array",
            expected_types=np.ndarray
        )

        array_length = len(x_array)

        if len(y_array) != array_length:
            raise Exception(f"Y array has length {len(y_array)} while X array has length {len(x_array)}")

        if z_array is not None:

            check_type(
                var=z_array,
                name="Z array",
                expected_types=np.ndarray
            )

            if len(z_array) != array_length:
                raise Exception(f"Z array has length {len(z_array)} while X array has length {len(x_array)}")

        else:

            z_array = np.zeros_like(x_array)

        if t_array is not None:

            check_type(
                var=t_array,
                name="T array",
                expected_types=np.ndarray
            )

            if len(t_array) != array_length:
                raise Exception(f"T array has length {len(t_array)} while X array has length {len(x_array)}")

        else:

            t_array = np.zeros_like(x_array)

        self._epsg_code = epsg_code
        self._x_array = x_array
        self._y_array = y_array
        self._z_array = z_array
        self._t_array = t_array

    @classmethod
    def fromPoints(cls,
                   points: List[Point],
                   epsg_code: numbers.Integral = None,
                   crs_check: bool = True
                   ):
        """

        :param points: list of points
        :type points: List[Point]
        :param epsg_code: optional EPSG code
        :type epsg_code: numbers.Integral
        :param crs_check: whether to check points crs
        :type crs_check: bool
        """

        for ndx, point in enumerate(points):

            check_type(point, "Input point {}".format(ndx), Point)

        if not epsg_code:
            epsg_code = points[0].epsg_code()

        if crs_check:

            for ndx, point in enumerate(points):

                if point.epsg_code() != epsg_code:

                    raise Exception("Point {} has EPSG code {} but {} required".format(ndx, point.epsg_code(), epsg_code))

        return Points(
            epsg_code=epsg_code,
            x_array=np.array([p.x for p in points]),
            y_array=np.array([p.y for p in points]),
            z_array=np.array([p.z for p in points]),
            t_array=np.array([p.t for p in points])
        )

    @property
    def xs(self):
        """
        The points x values.

        :return: points x values
        :rtype: float
        """

        return self._x_array

    @property
    def ys(self):
        """
        The points y values.

        :return: points y values
        :rtype: float
        """

        return self._y_array

    @property
    def zs(self):
        """
        The points z values.

        :return: points z values
        :rtype: float
        """

        return self._z_array


    @property
    def ts(self):
        """
        The points t values.

        :return: points t values
        :rtype: float
        """

        return self._t_array

    def epsg_code(self) -> numbers.Integral:
        """
        The points EPSG code.

        :return: the points EPSG code
        :rtype: numbers.Integral
        """

        return self._epsg_code

    @property
    def crs(self) -> Crs:
        """
        The points CRS.

        :return: the points CRS
        :rtype: Crs
        """

        return Crs(self.epsg_code())

    def asXyzArray(self):
        """
        Convert to a Numpy x-y-z array
        """

        return np.vstack(
            (
                self.xs,
                self.ys,
                self.zs
            )
        ).transpose()

    def x_min(self):
        """
        The minimum x value.
        """

        return np.nanmin(self.xs)

    def x_max(self):
        """
        The maximum x value.
        """

        return np.nanmax(self.xs)

    def x_mean(self):
        """
        The mean x value.
        """

        return np.nanmean(self.xs)

    def y_min(self):
        """
        The minimum y value.
        """

        return np.nanmin(self.ys)

    def y_max(self):
        """
        The maximum y value.
        """

        return np.nanmax(self.ys)

    def y_mean(self):
        """
        The mean y value.
        """

        return np.nanmean(self.ys)

    def z_min(self):
        """
        The minimum z value.
        """

        return np.nanmin(self.zs)

    def z_max(self):
        """
        The maximum z value.
        """

        return np.nanmax(self.zs)

    def z_mean(self):
        """
        The mean z value.
        """

        return np.nanmean(self.zs)

    def t_min(self):
        """
        The minimum t value.
        """

        return np.nanmin(self.ts)

    def t_max(self):
        """
        The maximum t value.
        """

        return np.nanmax(self.ts)

    def t_mean(self):
        """
        The mean t value.
        """

        return np.nanmean(self.ts)

    def nanmean_point(self) -> Point:
        """
        Returns the nan- excluded mean point of the collection.
        It is the mean point for a collection of point in a x-y-z frame (i.e., not lat-lon).

        :return: the nan- excluded mean point of the collection.
        :rtype: Point
        """

        return Point(
            x=np.nanmean(self.xs),
            y=np.nanmean(self.ys),
            z=np.nanmean(self.zs),
            t=np.nanmean(self.ts),
            epsg_code=self.epsg_code()
        )


class PointSegmentCollection(list):
    """
    Collection of point or segment elements.

    """

    def __init__(
            self,
            geoms: Optional[List[Union[Point, Segment]]] = None,
            epsg_code: Optional[numbers.Integral] = None
    ):

        if geoms is not None:

            for geom in geoms:
                check_type(geom, "Spatial element", (Point, Segment))

        if epsg_code is not None:
            check_type(
                var=epsg_code,
                name="EPSG code",
                expected_types=numbers.Integral
            )

        if geoms is not None and epsg_code is not None:

            for geom in geoms:
                check_epsg(
                    spatial_element=geom,
                    epsg_code=epsg_code
                )

        elif geoms is not None and len(geoms) > 0:

            epsg_code = geoms[0].epsg_code()

        if geoms is not None and len(geoms) > 0:

            super(PointSegmentCollection, self).__init__(geoms)

        else:

            super(PointSegmentCollection, self).__init__()

        self.epsg_code = epsg_code

    def append(self,
               spatial_element: Union[Point, Segment]
               ) -> None:

        check_type(
            var=spatial_element,
            name="Spatial element",
            expected_types=(Point, Segment)
        )

        if self.epsg_code is not None:

            check_epsg(
                spatial_element=spatial_element,
                epsg_code=self.epsg_code
            )

        else:

            self.epsg_code = spatial_element.epsg_code()

        self.append(spatial_element)


class PointSegmentCollections(list):

    def __init__(self, atts: List[Tuple[Union[str, numbers.Integral], PointSegmentCollection]]):

        check_type(atts, "Point-segment collections", List)
        for label, spat_element in atts:
            check_type(label, "Label", (str, numbers.Integral))
            check_type(spat_element, "Point-segment collection", PointSegmentCollection)

        super(PointSegmentCollections, self).__init__(atts)


class MultiLine(object):
    """
    MultiLine is a list of Line objects, each one with the same CRS.
    """

    def __init__(self, lines: Optional[List[Points]] = None, epsg_cd: numbers.Integral = -1):

        if lines is None:
            lines = []

        if lines and epsg_cd == -1:
            epsg_cd = lines[0].epsg_code()

        for ndx in range(len(lines)):
            if lines[ndx].epsg_code() != epsg_cd:
                raise Exception("Input line with index {} should have EPSG code {} but has {}".format(
                    ndx,
                    epsg_cd,
                    lines[ndx].epsg_code()
                ))

        self._lines = lines
        self._crs = Crs(epsg_cd)

    def lines(self):

        return self._lines

    @property
    def crs(self) -> Crs:

        return self._crs

    def epsg(self) -> numbers.Integral:

        return self._crs.epsg_code()

    def num_lines(self):

        return len(self.lines())

    def num_tot_pts(self) -> numbers.Integral:

        num_points = 0
        for line in self._lines:
            num_points += line.num_pts()

        return num_points

    def line(self, ln_ndx: numbers.Integral = 0) -> Optional[Points]:
        """
        Extracts a line from the multiline instance, based on the provided index.

        :return: Line instance or None when ln_ndx is out-of-range.
        :rtype: Optional[Line].
        """

        num_lines = self.num_lines()
        if num_lines == 0:
            return None

        if ln_ndx not in range(num_lines):
            return None

        return self.lines()[ln_ndx]

    def __iter__(self):
        """
        Return the elements of a MultiLine, i.e., its lines.
        """

        return (self.line(i) for i in range(0, self.num_lines()-1))

    def __repr__(self) -> str:
        """
        Represents a MultiLine instance as a shortened text.

        :return: a textual shortened representation of a MultiLine instance.
        :rtype: basestring.
        """

        num_lines = self.num_lines()
        num_tot_pts = self.num_tot_pts()
        epsg = self.epsg()

        txt = "MultiLine with {} line(s) and {} total point(s) - EPSG: {}".format(num_lines, num_tot_pts, epsg)

        return txt

    def __len__(self):
        """
        Return number of lines.

        :return: number of lines
        :rtype: numbers.Integral
        """

        return self.num_lines()

    def add_line(self, line) -> bool:
        """
        In-place addition of a Line instance (that is not cloned).

        :param line: the line to add.
        :type line: Line.
        :return: status of addition. True when added, False otherwise.
        :rtype: bool.
        """

        if self.num_lines() == 0 and not self.crs.valid():
            self._crs = line.crs

        if self.num_lines() > 0 and line.crs != self.crs:
            return False

        self._lines += [line]
        return True

    def clone(self) -> 'MultiLine':

        return MultiLine(
            lines=[line.clone() for line in self._lines],
            epsg_cd=self.epsg()
        )

    def x_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.x_min() for line in self.lines()]))

    def x_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.x_max() for line in self.lines()]))

    def y_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.y_min() for line in self.lines()]))

    def y_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.y_max() for line in self.lines()]))

    def z_min(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmin([line.z_min() for line in self.lines()]))

    def z_max(self) -> Optional[numbers.Real]:

        if self.num_tot_pts() == 0:
            return None
        else:
            return float(np.nanmax([line.z_max() for line in self.lines()]))

    def is_continuous(self) -> bool:
        """
        Checks whether all lines in a multiline are connected.

        :return: whether all lines are connected.
        :rtype: bool.
        """

        if len(self._lines) <= 1:
            return False

        for line_ndx in range(len(self._lines) - 1):
            first = self._lines[line_ndx]
            second = self._lines[line_ndx + 1]
            if not analizeJoins(first, second):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines()) - 1):
            if not self.lines()[line_ndx].pt(-1).isCoinc3D(self.lines()[line_ndx + 1].pt(0)):
                return False

        return True

    def to_line(self):

        return Points([point for line in self._lines for point in line.pts()], epsg_cd=self.epsg())

    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines():
            lDensifiedLines.append(line.densify_2d_line(sample_distance))

        return MultiLine(lDensifiedLines, self.epsg())

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines():
            cleaned_lines.append(line.remove_coincident_points())

        return MultiLine(cleaned_lines, self.epsg())

    def intersectSegment(self,
        segment: Segment
    ) -> List[Optional[Union[Point, 'Segment']]]:
        """
        Calculates the possible intersection between the multiline and a provided segment.

        :param segment: the input segment
        :type segment: Segment
        :return: the possible intersections, points or segments
        :rtype: List[List[Optional[Union[Point, 'Segment']]]]
        """

        check_type(segment, "Input segment", Segment)
        check_crs(self, segment)

        intersections = []
        for line in self:
            intersections.extend(line.intersectSegment(segment))

        return intersections


class MultiLines(list):
    """
    Collection of multilines, inheriting from list.

    """

    def __init__(self,
                 multilines: List[MultiLine] = None
                 ):

        if multilines:

            check_type(multilines, "MultiLines", List)
            for el in multilines:
                check_type(el, "MultiLine", MultiLine)

            super(MultiLines, self).__init__(multilines)

        else:

            super(MultiLines, self).__init__()

    def append(self,
               item: MultiLine
               ) -> None:

        check_type(item, "MultiLine", MultiLine)
        super(MultiLines, self).append(item)


class SimpleGeometryCollections(list):

    def __init__(self,
                 epsg_code: numbers.Integral = 4326):

        super(SimpleGeometryCollections, self).__init__()

        self.epsg_code = epsg_code

    def append(self, geom):

        if not isinstance(geom, (Point, Segment, Points)):
            raise Exception(f"Expected Point, Segment or Line but got {type(geom)}")

        if geom.epsg_code() != self.epsg_code:
            raise Exception(f"Expected {self.epsg_code} EPSG code but got {geom.epsg_code()}")

        self.append(geom)