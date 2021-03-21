
from copy import copy

import numbers
from typing import Optional, List, Union, Tuple

from shapely.geometry import LineString
from shapely.geometry.polygon import Polygon
from shapely.geometry.multipolygon import MultiPolygon

from pygsf.georeferenced.crs import *
from pygsf.geometries.shapes.space2d import *
from pygsf.geometries.shapes.space3d import *
from pygsf.utils.types import *


class GeoPoints2D:
    """
    Collection of points.
    """

    def __init__(self,
                 epsg_code: numbers.Integral,
                 x_array: array,
                 y_array: array
                 ):
        """
        Construct a point list from a set of array values and an EPSG code.

        :param epsg_code: the EPSG code of the points
        :param x_array: the array storing the x values
        :param y_array: the array storing the y values
        """

        check_type(
            var=epsg_code,
            name="EPSG code",
            expected_types=numbers.Integral
        )

        check_type(
            var=x_array,
            name="X array",
            expected_types=array
        )

        check_type(
            var=y_array,
            name="Y array",
            expected_types=array
        )

        array_length = len(x_array)

        if len(y_array) != array_length:
            raise Exception(f"Y array has length {len(y_array)} while X array has length {len(x_array)}")

        self._epsg_code = epsg_code
        self._x_array = x_array
        self._y_array = y_array

    def num_pts(self
                ) -> int:
        """
        Numbers of points.
        """

        return len(self._x_array)

    @classmethod
    def fromPoints(cls,
                   points: List[Point2D],
                   epsg_code: numbers.Integral
                   ):
        """

        :param points: list of points
        :param epsg_code: optional EPSG code
        """

        for ndx, point in enumerate(points):

            check_type(point, "Input point {}".format(ndx), Point2D)

        return GeoPoints2D(
            epsg_code=epsg_code,
            x_array=array('d', [p.x for p in points]),
            y_array=array('d', [p.y for p in points])
        )

    @property
    def xs(self
           ) -> array:
        """
        Returns a copy of the points x values.

        :return: points x values
        """

        return copy(self._x_array)

    @property
    def ys(self
           ) -> array:
        """
        Returns a copy of the points y values.

        :return: points y values
        """

        return copy(self._y_array)

    def pt(self, pt_ndx: numbers.Integral) -> Point2D:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :return: the extracted Point instance.

        Examples:
        """

        return Point2D(
            x=self._x_array[pt_ndx],
            y=self._y_array[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
                  ) -> Tuple[float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x and y values
        """

        return (
            self._x_array[ndx],
            self._y_array[ndx]
        )

    def pts(self):

        return [Point2D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def __repr__(self) -> str:
        """
        Represents a GeoPoints2D instance as a shortened text.

        :return: a textual shortened representation of a Points instance.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty GeoPoints2D"
        else:
            x1, y1 = self.values_at(0)
            if num_points == 1:
                txt = "GeoPoints2D with unique point: {.4f}.{.4f}".format(x1, y1)
            else:
                x2, y2 = self.values_at(self.num_pts()-1)
                txt = "GeoPoints2D with {} points: ({:.4f}, {:.4f}) ... ({:.4f}, {:.4f})".format(
                    num_points, x1, y1, x2, y2)

        return txt

    def __iter__(self):
        """
        Return each point.
        """

        return (self.pt(ndx) for ndx in range(self.num_pts()))

    @property
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

        return Crs(self.epsg_code)

    def asXyArray(self):
        """
        Convert to a Numpy x-y array
        """

        return np.vstack(
            (
                self.xs,
                self.ys
            )
        ).transpose()

    def add_pt(self, pt) -> None:
        """
        In-place transformation of the original Line instance
        by adding a new point at the end.

        :param pt: the point to add
        :return: nothing
        """

        self._x_array.append(pt.x)
        self._y_array.append(pt.y)

    def add_pts(self,
                pts: List[Point2D]):
        """
        In-place transformation of the original Points instance
        by adding a new set of points at the end.

        :param pts: list of Points.
        """

        check_type(pts, "GeoPoints2D", GeoPoints2D)
        check_crs(self, pts)

        self._x_array.extend(pts.xs)
        self._y_array.extend(pts.ys)

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = GeoPoints2D([[0, 0], [1, 0], [0, 1]])
          >>> l.x_min()
          0.0
          >>> m = GeoPoints2D([])
          >>> m.x_min()
          None
        """

        return np.nanmin(self._x_array) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum x value.
        """

        return np.nanmax(self._x_array) if self.num_pts() > 0 else None

    def x_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean x value.
        """

        return np.nanmean(self._x_array) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum y value.
        """

        return np.nanmin(self._y_array) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum y value.
        """

        return np.nanmax(self._y_array) if self.num_pts() > 0 else None

    def y_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean y value.
        """

        return np.nanmean(self._y_array) if self.num_pts() > 0 else None

    def nanmean_point(self) -> Point2D:
        """
        Returns the nan- excluded mean point of the collection.
        It is the mean point for a collection of point in a x-y frame (i.e., not lat-lon).

        :return: the nan- excluded mean point of the collection.
        :rtype: Point2D
        """

        return Point2D(
            x=np.nanmean(self._x_array),
            y=np.nanmean(self._y_array)
        )

    def segment(self,
        ndx: int
    ) -> Optional[Segment2D]:
        """
        Returns the optional segment starting at index ndx.

        :param ndx: the segment index.
        :return: the optional segment
        """

        if ndx < 0 or ndx >= self.num_pts() - 1:
            return None

        return Segment2D(
            start_pt=self.pt(ndx),
            end_pt=self.pt(ndx + 1)
        )

    def reversed(self) -> 'GeoPoints2D':
        """
        Return a Points instance with reversed point list.

        :return: a new Points instance.
        """

        xs = self._x_array.reverse()
        ys = self._y_array.reverse()

        return GeoPoints2D(
            epsg_code=self.epsg_code,
            x_array=xs,
            y_array=ys
        )


class GeoPoints3D:
    """
    Collection of points.
    """

    def __init__(self,
                 epsg_code: numbers.Integral,
                 x_array: array,
                 y_array: array,
                 z_array: Optional[array] = None
                 ):
        """
        Construct a point list from a set of array values and an EPSG code.

        :param epsg_code: the EPSG code of the points
        :param x_array: the array storing the x values
        :param y_array: the array storing the y values
        :param z_array: the optional array storing the z values
        """

        check_type(
            var=epsg_code,
            name="EPSG code",
            expected_types=numbers.Integral
        )

        check_type(
            var=x_array,
            name="X array",
            expected_types=array
        )

        check_type(
            var=y_array,
            name="Y array",
            expected_types=array
        )

        array_length = len(x_array)

        if len(y_array) != array_length:
            raise Exception(f"Y array has length {len(y_array)} while X array has length {len(x_array)}")

        if z_array is not None:

            check_type(
                var=z_array,
                name="Z array",
                expected_types=array
            )

            if len(z_array) != array_length:
                raise Exception(f"Z array has length {len(z_array)} while X array has length {len(x_array)}")

        else:

            z_array = np.zeros_like(x_array)

        self._epsg_code = epsg_code
        self._x_array = x_array
        self._y_array = y_array
        self._z_array = z_array

    def num_pts(self
                ) -> int:
        """
        Numbers of points.
        """

        return len(self._x_array)

    @classmethod
    def fromPoints(cls,
                   points: List[Point3D],
                   epsg_code: numbers.Integral
                   ):
        """

        :param points: list of points
        :param epsg_code: optional EPSG code
        """

        for ndx, point in enumerate(points):

            check_type(point, "Input point {}".format(ndx), Point3D)

        return GeoPoints3D(
            epsg_code=epsg_code,
            x_array=array('d', [p.x for p in points]),
            y_array=array('d', [p.y for p in points]),
            z_array=array('d', [p.z for p in points])
        )

    @property
    def xs(self
           ) -> array:
        """
        Returns a copy of the points x values.

        :return: points x values
        """

        return copy(self._x_array)

    @property
    def ys(self
           ) -> array:
        """
        Returns a copy of the points y values.

        :return: points y values
        """

        return copy(self._y_array)

    @property
    def zs(self
           ) -> array:
        """
        Returns a copy of the points z values.

        :return: points z values
        """

        return copy(self._z_array)

    def pt(self, pt_ndx: numbers.Integral) -> Point3D:
        """
        Extract the point at index pt_ndx.

        :param pt_ndx: point index.
        :type pt_ndx: numbers.Integral.
        :return: the extracted Point instance.
        :rtype: Point.

        Examples:
        """

        return Point3D(
            x=self._x_array[pt_ndx],
            y=self._y_array[pt_ndx],
            z=self._z_array[pt_ndx]
        )

    def values_at(self,
        ndx: numbers.Integral
                  ) -> Tuple[float, float, float]:
        """
        Return the values at given index.

        :param ndx: the index of the point values to extract
        :type ndx: numbers.Integral
        :return: the x, y and z values
        """

        return (
            self._x_array[ndx],
            self._y_array[ndx],
            self._z_array[ndx]
        )

    def pts(self):

        return [Point3D(*self.values_at(ndx)) for ndx in range(self.num_pts())]

    def __repr__(self) -> str:
        """
        Represents a GeoPoints3D instance as a shortened text.

        :return: a textual shortened representation of a Points instance.
        """

        num_points = self.num_pts()

        if num_points == 0:
            txt = "Empty GeoPoints3D"
        else:
            x1, y1, z1 = self.values_at(0)
            if num_points == 1:
                txt = "GeoPoints3D with unique point: {.4f}.{.4f},{.4f}".format(x1, y1, z1)
            else:
                x2, y2, z2 = self.values_at(self.num_pts()-1)
                txt = "GeoPoints3D with {} points: ({:.4f}, {:.4f}, {:.4f}) ... ({:.4f}, {:.4f}, {:.4f})".format(
                    num_points, x1, y1, z1, x2, y2, z2)

        return txt

    def __iter__(self):
        """
        Return each point.
        """

        return (self.pt(ndx) for ndx in range(self.num_pts()))

    @property
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

        return Crs(self.epsg_code)

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

    def add_pt(self, pt) -> None:
        """
        In-place transformation of the original Line instance
        by adding a new point at the end.

        :param pt: the point to add
        :return: nothing
        """

        self._x_array.append(pt.x)
        self._y_array.append(pt.y)
        self._z_array.append(pt.z)

    def add_pts(self,
                pts: List[Point3D]):
        """
        In-place transformation of the original Points instance
        by adding a new set of points at the end.

        :param pts: list of Points.
        """

        check_type(pts, "GeoPoints3D", GeoPoints3D)
        check_crs(self, pts)

        self._x_array.extend(pts.xs)
        self._y_array.extend(pts.ys)
        self._z_array.extend(pts.zs)

    def x_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum of x values.

        :return: the optional minimum of x values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = GeoPoints3D([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
          >>> l.x_min()
          0.0
          >>> m = GeoPoints3D([])
          >>> m.x_min()
          None
        """

        return np.nanmin(self._x_array) if self.num_pts() > 0 else None

    def x_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum x value.
        """

        return np.nanmax(self._x_array) if self.num_pts() > 0 else None

    def x_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean x value.
        """

        return np.nanmean(self._x_array) if self.num_pts() > 0 else None

    def y_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum y value.
        """

        return np.nanmin(self._y_array) if self.num_pts() > 0 else None

    def y_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum y value.
        """

        return np.nanmax(self._y_array) if self.num_pts() > 0 else None

    def y_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean y value.
        """

        return np.nanmean(self._y_array) if self.num_pts() > 0 else None

    def z_min(self) -> Optional[numbers.Real]:
        """
        Optional minimum z value.
        """

        return np.nanmin(self._z_array) if self.num_pts() > 0 else None

    def z_max(self) -> Optional[numbers.Real]:
        """
        Optional maximum z value.
        """

        return np.nanmax(self._z_array) if self.num_pts() > 0 else None

    def z_mean(self) -> Optional[numbers.Real]:
        """
        Optional mean z value.
        """

        return np.nanmean(self._z_array) if self.num_pts() > 0 else None

    def z_var(self) -> Optional[numbers.Real]:
        """
        Optional variance of z values.

        :return: the optional variance of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = GeoPoints3D.fromPoints([[0, 0, 2], [1, 0, 2], [0, 1, 2]], epsg_code=32633)
          >>> l.z_var()
          0.0
        """

        return np.nanvar(self._z_array) if self.num_pts() > 0 else None

    def z_std(self) -> Optional[numbers.Real]:
        """
        Optional standard deviation of z values.

        :return: the optional standard deviation of z values.
        :rtype: Optional[numbers.Real]

        Examples:
          >>> l = GeoPoints3D.fromPoints([[0, 0, 2], [1, 0, 2], [0, 1, 2]], epsg_code=32633)
          >>> l.z_std()
          0.0
        """

        return np.nanstd(self._z_array) if self.num_pts() > 0 else None

    def nanmean_point(self) -> Point3D:
        """
        Returns the nan- excluded mean point of the collection.
        It is the mean point for a collection of point in a x-y-z frame (i.e., not lat-lon-z).

        :return: the nan- excluded mean point of the collection.
        :rtype: Point3D
        """

        return Point3D(
            x=np.nanmean(self._x_array),
            y=np.nanmean(self._y_array),
            z=np.nanmean(self._z_array)
        )

    def segment(self,
        ndx: int
    ) -> Optional[Segment3D]:
        """
        Returns the optional segment starting at index ndx.

        :param ndx: the segment index.
        :return: the optional segment
        """

        if ndx < 0 or ndx >= self.num_pts() - 1:
            return None

        return Segment3D(
            start_pt=self.pt(ndx),
            end_pt=self.pt(ndx + 1)
        )

    def reversed(self) -> 'GeoPoints3D':
        """
        Return a Points instance with reversed point list.

        :return: a new Points instance.
        """

        xs = self._x_array.reverse()
        ys = self._y_array.reverse()
        zs = self._z_array.reverse()

        return GeoPoints3D(
            epsg_code=self.epsg_code,
            x_array=xs,
            y_array=ys,
            z_array=zs
        )


class GeoPointSegmentCollection(list):
    """
    Collection of point or segment elements.

    """

    def __init__(
            self,
            geoms: Optional[List[Union[Point3D, Segment3D]]] = None,
            epsg_code: Optional[numbers.Integral] = None
    ):

        if geoms is not None:
            for geom in geoms:
                check_type(geom, "Spatial element", (Point3D, Segment3D))

        if epsg_code is not None:
            check_type(
                var=epsg_code,
                name="EPSG code",
                expected_types=numbers.Integral
            )

        if geoms is not None and len(geoms) > 0:

            super(GeoPointSegmentCollection, self).__init__(geoms)

        else:

            super(GeoPointSegmentCollection, self).__init__()

        self.epsg_code = epsg_code

    def append(self,
               spatial_element: Union[Point3D, Segment3D]
               ) -> None:

        check_type(
            var=spatial_element,
            name="Spatial element",
            expected_types=(Point3D, Segment3D)
        )

        self.append(spatial_element)


class GeoPointSegmentCollections(list):

    def __init__(self, atts: List[Tuple[Union[str, numbers.Integral], GeoPointSegmentCollection]]):

        check_type(atts, "Point-segment collections", List)
        for label, spat_element in atts:
            check_type(label, "Label", (str, numbers.Integral))
            check_type(spat_element, "Point-segment collection", GeoPointSegmentCollection)

        super(GeoPointSegmentCollections, self).__init__(atts)


class GeoLines3D(list):
    """
    Collection of lines.

    """

    def __init__(self,
                 lines: Optional[List[Line3D]] = None
                 ):

        if lines:

            check_type(lines, "Lines", List)
            for line in lines:
                check_type(line, "Line", Line3D)
            first_line = lines[0]
            for line in lines[1:]:
                check_crs(first_line, line)

            super(GeoLines3D, self).__init__(lines)

        else:

            super(GeoLines3D, self).__init__()

    def append(self,
               item: Line3D
               ) -> None:

        check_type(item, "Line", Line3D)
        if len(self) > 0:
            check_crs(self[0], item)

        super(GeoLines3D, self).append(item)

    def intersectSegment(self,
        segment: Segment3D
    ) -> List[Optional[Union[Point3D, 'Segment3D']]]:
        """
        Calculates the possible intersection between the multiline and a provided segment.

        :param segment: the input segment
        :type segment: Segment3D
        :return: the possible intersections, points or segments
        :rtype: List[List[Optional[Union[Point, 'Segment']]]]
        """

        check_type(segment, "Input segment", Segment3D)
        #check_crs(self, segment)

        intersections = []
        for line in self:
            intersections.extend(line.intersectSegment(segment))

        return intersections


class GeoMultiLine(object):
    """
    MultiLine is a list of Line objects, each one with the same CRS.
    """

    def __init__(self, lines: Optional[List[Line3D]] = None, epsg_cd: numbers.Integral = -1):

        if lines is None:
            lines = []

        self._lines = lines
        self._crs = Crs(epsg_cd)

    def lines(self):

        return self._lines

    @property
    def crs(self) -> Crs:

        return self._crs

    @property
    def epsg_code(self) -> numbers.Integral:

        return self._crs.epsg_code

    def num_lines(self):

        return len(self.lines())

    def num_tot_pts(self) -> numbers.Integral:

        num_points = 0
        for line in self._lines:
            num_points += line.num_pts()

        return num_points

    def line(self, ln_ndx: numbers.Integral = 0) -> Optional[GeoPoints3D]:
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
        epsg = self.epsg_code

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

    def clone(self) -> 'GeoMultiLine':

        return GeoMultiLine(
            lines=[line.clone() for line in self._lines],
            epsg_cd=self.epsg_code
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

        return Line3D([point for line in self._lines for point in line.pts()])

    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines():
            lDensifiedLines.append(line.densify_2d_line(sample_distance))

        return GeoMultiLine(lDensifiedLines, self.epsg_code)

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines():
            cleaned_lines.append(line.remove_coincident_points())

        return GeoMultiLine(cleaned_lines, self.epsg_code)

    def intersectSegment(self,
        segment: Segment3D
    ) -> List[Optional[Union[Point3D, 'Segment3D']]]:
        """
        Calculates the possible intersection between the multiline and a provided segment.

        :param segment: the input segment
        :type segment: Segment3D
        :return: the possible intersections, points or segments
        :rtype: List[List[Optional[Union[Point, 'Segment']]]]
        """

        check_type(segment, "Input segment", Segment3D)
        check_crs(self, segment)

        intersections = []
        for line in self:
            intersections.extend(line.intersectSegment(segment))

        return intersections


class GeoMultiLines(list):
    """
    Collection of multilines, inheriting from list.

    """

    def __init__(self,
                 multilines: List[GeoMultiLine] = None
                 ):

        if multilines:

            check_type(multilines, "MultiLines", List)
            for el in multilines:
                check_type(el, "MultiLine", GeoMultiLine)

            super(GeoMultiLines, self).__init__(multilines)

        else:

            super(GeoMultiLines, self).__init__()

    def append(self,
               item: GeoMultiLine
               ) -> None:

        check_type(item, "MultiLine", GeoMultiLine)
        super(GeoMultiLines, self).append(item)


class GeoMPolygon:
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
                       ) -> GeoLines3D:
        """
        Determine the intersections between a mpolygon and a line.

        :param line: the line
        :type line: shapely.geometry.LineString
        :return: the intersecting lines
        :rtype: GeoLines3D
        """

        lines = GeoLines3D()

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


class GeoSimpleGeometryCollections(list):

    def __init__(self,
                 epsg_code: numbers.Integral = 4326):

        super(GeoSimpleGeometryCollections, self).__init__()

        self.epsg_code = epsg_code

    def append(self, geom):

        if not isinstance(geom, (Point3D, Segment3D, Line3D)):
            raise Exception(f"Expected Point, Segment or Line but got {type(geom)}")

        self.append(geom)


def line_from_shapely(
        shapely_geom: LineString,
        epsg_code: numbers.Integral
) -> GeoPoints3D:
    # Side effects: none
    """
    Create a Line instance from a shapely Linestring instance.

    :param shapely_geom: the shapely input LineString instance
    :type shapely_geom: shapely.geometry.linestring.LineString
    :param epsg_code: the EPSG code of the LineString instance
    :type epsg_code: numbers.Integral
    :return: the converted Line instance
    :rtype: Line3D
    """

    x_array, y_array = shapely_geom.xy

    return Line3D.fromArrays(
        x_array,
        y_array
    )


def line_to_shapely(
        src_line: GeoPoints3D
) -> LineString:
    """
    Create a shapely.LineString instance from a Line one.

    :param src_line: the source line to convert to the shapely format
    :type src_line: Line3D
    :return: the shapely LineString instance and the EPSG code
    :rtype: Tuple[LineString, numbers.Integral]
    """

    return LineString(src_line.xy_zipped()), src_line.epsg_code


class GeoSegments(list):
    """
    Collection of segments, inheriting from list.

    """

    def __init__(self, segments: List[Segment3D]):

        check_type(segments, "Segments", List)
        for el in segments:
            check_type(el, "Segment", Segment3D)

        super(GeoSegments, self).__init__(segments)

