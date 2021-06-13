
from typing import Union
from math import sqrt
import datetime
from copy import deepcopy

from ...mathematics.vectors3d import *
from .space2d import Point2D, Line2D
from .space3d import Point3D, Line3D, Segment3D


class Point4D(Point3D):
    """
    Cartesian point.
    Dimensions: 3D + time
    """

    def __init__(self,
                 x: numbers.Real,
                 y: numbers.Real,
                 z: numbers.Real,
                 t: Optional[datetime.datetime] = None
                 ):
        """
        Construct a Point instance given 3 float values and an optional time value (seconds from ....).
        """

        super(Point4D, self).__init__(x, y, z)

        if t is not None:
            check_type(t, "Input t value", datetime.datetime)

        self._t = t

    def __repr__(self):

        return f"Point4D({self.x:.4f}, {self.y:.4f}, {self.z:.4f}, {self.t})"

    '''
    @classmethod
    def from_array(cls, a):
        """
        Class method to construct a point from a numpy 1x4 array.

        Example:
          >>> Point4D.from_array(np.array([1, 0, 1]))
          Point4D(1.0000, 0.0000, 1.0000, None)
        """

        obj = cls(
            x=0.0,
            y=0.0,
            z=0.0,
            t=None
        )

        assert 2 <= a.size <= 3
        b = a.astype(np.float64)
        if b.size == 2:
            c = np.append(b, [np.nan])
        else:
            c = b
        obj._p = c
        obj._t = None
        return obj
    '''

    '''
    @property
    def x(self):
        """
        Return x value

        Example:
          >>> Point4D(1.5, 1, 1).x
          1.5
        """

        return self._p[0]
    '''

    '''
    @property
    def y(self):
        """
        Return y value

        Example:
          >>> Point4D(1.5, 3.0, 1).y
          3.0
        """
        return self._p[1]
    '''

    '''
    @property
    def z(self):
        """
        Return z value

        Example:
          >>> Point4D(1.5, 3.2, 41.).z
          41.0
        """
        return self._p[2]
    '''

    @property
    def t(self):
        """
        Return time value

        Example:
          >>> Point4D(1.5, 3.2, 41.).t
          None
        """
        return self._t

    def a(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real, datetime.datetime]:
        """
        Return the individual values of the point.

        :return: x, y, z and t values

        Examples:
          >>> Point4D(4, 3, 7, datetime.datetime(2020, 12, 4)).a()
          (4.0, 3.0, 7.0, datetime.datetime(2020, 12, 4, 0, 0))
        """

        return self.x, self.y, self.z, self.t

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y, z, t = Point4D(4, 3, 7, datetime.datetime(2020, 12, 4))
          >>> x == 4
          True
          >>> y == 3
          True
          >>> z == 7
          True
          >>> t == datetime.datetime(2020, 12, 4)
          True

        """

        return (i for i in self.a())

    def clone(self):
        """
        Clone the point.

        Example:
          >>> Point4D(1, 1, 1).clone()
          Point4D(1.0000, 1.0000, 1.0000, None)
        """

        return Point4D(
            x=self.x,
            y=self.y,
            z=self.z,
            t=self.t
        )

    def as_point2d(self) -> Point2D:
        """
        Convert a 4D point to a 2D point.
        """

        return Point2D(
            x=self.x,
            y=self.y
        )

    def as_point3d(self) -> Point3D:
        """
        Convert a 4D point to a 3D point.
        """

        return Point3D(
            x=self.x,
            y=self.y,
            z=self.z
        )

    def __sub__(self, another):
        """Return point difference

        Example:
          >>> Point4D(1., 1., 1.) - Point4D(1., 1., 1.)
          Point4D(0.0000, 0.0000, 0.0000, None)
        """

        return Point4D(
            x=self.x - another.x,
            y=self.y - another.y,
            z=self.z - another.z if np.isfinite(self.z) and np.isfinite(another.z) else np.nan,
            t=self.t - another.t if self.t is not None and another.t is not None else None
        )

    def __abs__(self):
        """
        Point distance from frame origin.
        todo: make sure it works as intended with nan values

        Example:
          >>> abs(Point4D(3, 4, 0))
          5.0
          >>> abs(Point4D(0, 12, 5))
          13.0
        """

        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def distance(self, another):
        """
        Calculate Euclidean spatial distance between two points.
        todo: make sure it works as intended with nan values

        Examples:
          >>> Point4D(1., 1., 1.).distance(Point4D(4., 5., 1,))
          5.0
          >>> Point4D(1, 1, 1).distance(Point4D(4, 5, 1))
          5.0
        """

        return self.distance_3d(another)

    def is_coincident(self,
                      other,
                      tolerance=MIN_SEPARATION_THRESHOLD
                      ) -> bool:
        """
        Check spatial coincidence of two points
        todo: make sure it works as intended with Points with nan values

        Example:
          >>> Point4D(1., 0., -1.).is_coincident(Point4D(1., 1.5, -1.))
          False
          >>> Point4D(1., 0., 0.).is_coincident(Point4D(1., 0., 0.))
          True
        """

        if self.distance_2d(other) > tolerance:
            return False
        elif self.distance(other) > tolerance:
            return False
        else:
            return True

    def translate(self, sx=0.0, sy=0.0, sz=0.0):
        """
        Create a new point shifted by given amount from the self instance.

        Example:
          >>> Point4D(1, 1, 1).translate(0.5, 1., 1.5)
          Point4D(1.5000, 2.0000, 2.5000, None)
       """

        return Point4D(
            self.x + sx,
            self.y + sy,
            self.z + sz,
            self.t
        )

    def vect_offset(self, displ_vect):
        """
        Create a new point from the self, with offsets defined by a vector.
        """

        return Point4D(self.x + displ_vect.x,
                       self.y + displ_vect.y,
                       self.z + displ_vect.z,
                       self.t)

    def vector(self):
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point4D(1, 1, 0).vector()
          Vect3D(1.0000, 1.0000, 0.0000)
        """

        return Vect3D(self.x, self.y, self.z)

    def delta_time(self, another):
        """
        Calculate the time difference between two points
        """

        return another.t - self.t

    def speed(self, another):
        """
        Calculate the speed required to displace self to another.
        """

        try:

            return self.distance(another) / self.delta_time(another)

        except:

            return np.nan


class Segment4D(Segment3D):
    """
    Segment is a geometric object defined by a straight line between
    two points.
    """

    def __init__(self,
                 start_pt: Point4D,
                 end_pt: Point4D):

        check_type(start_pt, "Start point", Point4D)
        check_type(end_pt, "End point", Point4D)

        super(Segment4D, self).__init__(start_pt, end_pt)

    def clone(self):

        return Segment4D(self.start_pt, self.end_pt)

    def increasing_x(self):

        if self.end_pt.x < self.start_pt.x:
            return Segment4D(self.end_pt, self.start_pt)
        else:
            return self.clone()

    @property
    def x_range(self):

        if self.start_pt.x < self.end_pt.x:
            return self.start_pt.x, self.end_pt.x
        else:
            return self.end_pt.x, self.start_pt.x

    @property
    def y_range(self):

        if self.start_pt.y < self.end_pt.y:
            return self.start_pt.y, self.end_pt.y
        else:
            return self.end_pt.y, self.start_pt.y

    @property
    def z_range(self):

        if self.start_pt.z < self.end_pt.z:
            return self.start_pt.z, self.end_pt.z
        else:
            return self.end_pt.z, self.start_pt.z

    @property
    def delta_x(self):

        return self.end_pt.x - self.start_pt.x

    @property
    def delta_y(self):

        return self.end_pt.y - self.start_pt.y

    @property
    def delta_z(self):

        return self.end_pt.z - self.start_pt.z

    @property
    def length(self):

        return self.length_3d()

    def vector(self):

        return Vect3D(self.delta_x,
                      self.delta_y,
                      self.delta_z
                      )

    """
    def segment_2d_m(self):

        return (self.end_pt.y - self.start_pt.y) / (self.end_pt.x - self.start_pt.x)

    def segment_2d_p(self):

        return self.start_pt.y - self.segment_2d_m() * self.start_pt.x

    def intersection_2d_pt(self, another):

        assert self.length_2d > 0.0
        assert another.length_2d > 0.0

        if self.start_pt.x == self.end_pt.x:  # self segment parallel to y axis
            x0 = self.start_pt.x
            try:
                m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            except:
                return None
            y0 = m1 * x0 + p1
        elif another.start_pt.x == another.end_pt.x:  # another segment parallel to y axis
            x0 = another.start_pt.x
            try:
                m1, p1 = self.segment_2d_m(), self.segment_2d_p()
            except:
                return None
            y0 = m1 * x0 + p1
        else:  # no segment parallel to y axis
            m0, p0 = self.segment_2d_m(), self.segment_2d_p()
            m1, p1 = another.segment_2d_m(), another.segment_2d_p()
            x0 = (p1 - p0) / (m0 - m1)
            y0 = m0 * x0 + p0

        return Point4D(x0, y0)

    def contains_2d_pt(self, pt2d):

        segment_length2d = self.length_2d
        segmentstart_pt2d_distance = self.start_pt.dist_2d(pt2d)
        segmentend_pt2d_distance = self.end_pt.dist_2d(pt2d)

        if segmentstart_pt2d_distance > segment_length2d or \
           segmentend_pt2d_distance > segment_length2d:
            return False
        else:
            return True

    def fast_2d_contains_pt(self, pt2d):
        '''
        to work properly, this function requires that the pt lies on the line defined by the segment
        '''

        range_x = self.x_range
        range_y = self.y_range

        if range_x[0] <= pt2d.x <= range_x[1] or \
           range_y[0] <= pt2d.y <= range_y[1]:
            return True
        else:
            return False
    """

    def scale(self, scale_factor):
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: float
        :return: Segment instance
        """

        delta_x = self.delta_x * scale_factor
        delta_y = self.delta_y * scale_factor
        delta_z = self.delta_z * scale_factor

        end_pt = Point4D(self.start_pt.x + delta_x,
                         self.start_pt.y + delta_y,
                         self.start_pt.z + delta_z)

        return Segment4D(self.start_pt,
                         end_pt)

    '''
    def densify_2d_segment(self, densify_distance):
        """
        Densify a segment by adding additional points
        separated a distance equal to densify_distance.
        The result is no longer a Segment instance, instead it is a Line instance.

        :param densify_distance: float
        :return: Line
        """

        assert densify_distance > 0.0

        length2d = self.length_2d

        assert length2d > 0.0

        vect = self.vector()
        vers_2d = vect.versor_2d()
        generator_vector = vers_2d.scale(densify_distance)

        assert generator_vector.len_2d > 0.0

        interpolated_line = Line4D([self.start_pt])
        n = 0
        while True:
            n += 1
            new_pt = self.start_pt.vect_offset(generator_vector.scale(n))
            distance = self.start_pt.dist_2d(new_pt)
            if distance >= length2d:
                break
            interpolated_line.add_pt(new_pt)
        interpolated_line.add_pt(self.end_pt)

        return interpolated_line
    '''


class Line4D(Line3D):
    """
    A list of Point objects.
    """

    def __init__(self,
                 x_seq: Optional[Union[Sequence[float], np.ndarray]] = None,
                 y_seq: Optional[Union[Sequence[float], np.ndarray]] = None,
                 z_seq: Optional[Union[Sequence[float], np.ndarray]] = None,
                 t_seq: Optional[Sequence[datetime.datetime]] = None
                 ):
        """
        Creates the Line4D instance.
        """

        super(Line4D, self).__init__(x_seq, y_seq, z_seq)

        if t_seq is None and self.num_pts() != 0:
            raise Exception(f"T input is None while points are {self.num_pts()}")

        if t_seq is not None and self.num_pts() == 0:
            raise Exception(f"T input is not None while points are {self.num_pts()}")

        self._times = t_seq

    def clone(self):

        return Line4D(
            x_seq=self.x_array().copy(),
            y_seq=self.y_array().copy(),
            z_seq=self.z_array().copy(),
            t_seq=deepcopy(self._times)
        )

    def add_pt(self, pt: Point4D):
        """
        In-place transformation of the original Line4D instance
        by adding a new point at the end.

        :param pt: the point to add
        """

        check_type(pt, "Point", Point4D)

        self._x_array = np.append(self._x_array, pt.x)
        self._y_array = np.append(self._y_array, pt.y)
        self._z_array = np.append(self._z_array, pt.z)
        self._times.append(pt.t)

    def add_pts(self, pt_list: List[Point4D]):
        """
        In-place transformation of the original Line instance
        by adding a new set of points at the end.

        :param pt_list: list of Points.
        """

        for pt in pt_list:
            check_type(pt, "Point", Point4D)

        x_vals = [pt.x for pt in pt_list]
        y_vals = [pt.y for pt in pt_list]
        z_vals = [pt.z for pt in pt_list]
        t_vals = [pt.t for pt in pt_list]

        self._x_array = np.append(self._x_array, x_vals)
        self._y_array = np.append(self._y_array, y_vals)
        self._z_array = np.append(self._z_array, z_vals)
        self._times.extend(t_vals)

    def as_segment(self) -> Segment4D:
        """Return the segment defined by line start and end points"""

        return Segment4D(self.start_pt(), self.end_pt())

    def segments(self):
        """
        Convert to a list of segments.

        :return: list of Segment objects
        """

        pts_pairs = list(zip(self[:-1], self[1:]))

        segments = [Segment4D(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    def as_line2d(self) -> Line2D:

        pts2d = []
        for pt4d in self:
            x, y, _, _ = pt4d
            pts2d.append(
                Point2D(
                    x=x,
                    y=y
                )
            )

        return Line2D(pts=pts2d)

    def join(self, another):
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line4D(self + another)

    def length(self) -> numbers.Real:

        return self.length_3d()

    def invert_direction(self):

        new_line = self.clone()
        new_line.reverse()  # in-place operation on new_line

        return new_line

    def absolute_slopes(self) -> np.ndarray:

        return np.asarray(list(map(abs, self.dir_slopes())))

    def segment(self,
        ndx: numbers.Integral
    ) -> Optional[Segment4D]:
        """
        Returns the optional segment at index ndx.

        :param ndx: the segment index.
        :return: the optional segment
        """

        start_pt = self.pt(ndx)
        end_pt = self.pt(ndx + 1)

        if start_pt.is_coincident(end_pt):
            return None
        else:
            return Segment4D(
                start_pt=self.pt(ndx),
                end_pt=self.pt(ndx + 1)
            )

    def reversed(self) -> 'Line4D':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        """

        pts = [pt.clone() for pt in self]
        pts.reverse()

        return Line4D(
            pts=pts
        )


class MultiLine4D(object):
    """
    MultiLine4D is a list of Line4D objects
    """

    def __init__(self, lines_list=None):

        if lines_list is None:
            lines_list = []
        self._lines = lines_list

    @property
    def lines(self):

        return self._lines

    def add(self, line):

        return MultiLine4D(self.lines + [line])

    def clone(self):

        return MultiLine4D(self.lines)

    @property
    def num_parts(self):

        return len(self.lines)

    @property
    def num_points(self):

        num_points = 0
        for line in self.lines:
            num_points += line.num_pts()

        return num_points

    @property
    def x_min(self):

        return np.nanmin([line.x_min for line in self.lines])

    @property
    def x_max(self):

        return np.nanmax([line.x_max for line in self.lines])

    @property
    def y_min(self):

        return np.nanmin([line.y_min for line in self.lines])

    @property
    def y_max(self):

        return np.nanmax([line.y_max for line in self.lines])

    @property
    def z_min(self):

        return np.nanmin([line.z_min for line in self.lines])

    @property
    def z_max(self):

        return np.nanmax([line.z_max for line in self.lines])

    def is_continuous(self):

        for line_ndx in range(len(self._lines) - 1):
            if not self.lines[line_ndx].pts()[-1].is_coincident(self.lines[line_ndx + 1].pts()[0]) or \
               not self.lines[line_ndx].pts()[-1].is_coincident(self.lines[line_ndx + 1].pts()[-1]):
                return False

        return True

    def is_unidirectional(self):

        for line_ndx in range(len(self.lines) - 1):
            if not self.lines[line_ndx].pts()[-1].is_coincident(self.lines[line_ndx + 1].pts()[0]):
                return False

        return True

    def to_line(self):

        return Line4D([point for line in self.lines for point in line.pts()])

    '''
    def crs_project(self, srcCrs, destCrs):

        lines = []
        for line in self.lines:
            lines.append(line.crs_project(srcCrs, destCrs))

        return MultiLine4D(lines)

    def densify_2d_multiline(self, sample_distance):

        lDensifiedLines = []
        for line in self.lines:
            lDensifiedLines.append(line.densify_2d_line(sample_distance))

        return MultiLine4D(lDensifiedLines)
    '''

    def remove_coincident_points(self):

        cleaned_lines = []
        for line in self.lines:
            cleaned_lines.append(line.remove_coincident_points())

        return MultiLine4D(cleaned_lines)


