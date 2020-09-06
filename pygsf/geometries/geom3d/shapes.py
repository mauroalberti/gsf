
from typing import List, Optional, Union

import numbers
from array import array

from pygsf.orientations.orientations import *
from pygsf.mathematics.statistics import *
from pygsf.mathematics.quaternions import *
from pygsf.geometries.geom2d.shapes import *
from pygsf.utils.types import *


class Point:
    """
    Cartesian point.
    Dimensions: 3D
    """

    def __init__(
        self,
        x: numbers.Real,
        y: numbers.Real,
        z: numbers.Real = 0.0
    ):
        """
        Construct a Point instance.

        :param x: point x coordinate.
        :type x: numbers.Real.
        :param y: point y coordinate.
        :type y: numbers.Real.
        :param z: point z coordinate.
        :type z: numbers.Real.
        """

        vals = [x, y, z]
        if any(map(lambda val: not isinstance(val, numbers.Real), vals)):
            raise Exception("Input values must be integer or float type")
        if not all(map(math.isfinite, vals)):
            raise Exception("Input values must be finite")

        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    @classmethod
    def fromVect(cls,
        vect: Vect) -> 'Point':
        """

        :param vect:
        :return:
        """

        return cls(
            x=vect.x,
            y=vect.y,
            z=vect.z
        )

    @property
    def x(self) -> numbers.Real:
        """
        Return the x coordinate of the current point.

        :return: x coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7).x
          4.0
          >>> Point(-0.39, 3, 7).x
          -0.39
        """

        return self._x

    @property
    def y(self) -> numbers.Real:
        """
        Return the y coordinate of the current point.

        :return: y coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7).y
          3.0
          >>> Point(-0.39, 17.42, 7).y
          17.42
        """

        return self._y

    @property
    def z(self) -> numbers.Real:
        """
        Return the z coordinate of the current point.

        :return: z coordinate.
        :rtype: numbers.Real

        Examples:
          >>> Point(4, 3, 7).z
          7.0
          >>> Point(-0.39, 17.42, 8.9).z
          8.9
        """

        return self._z

    def __iter__(self):
        """
        Return the elements of a Point.

        :return:

        Examples;
          >>> x, y, z = Point(1,1)
          >>> x == 1
          True
          >>> y == 1
          True

        """

        return (i for i in self.a())

    def __repr__(self) -> str:

        return "Point({:.4f}, {:.4f}, {:.4f})".format(self.x, self.y, self.z)

    def __eq__(self,
        another: 'Point'
    ) -> bool:
        """
        Return True if objects are equal.

        :param another: another point.
        :type another: Point.
        :raise: Exception.

        Example:
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, 1)
          True
          >>> Point(1., 1., 1.) == Point(1, 1, -1)
          False
        """

        if not isinstance(another, Point):
            raise Exception("Another instance must be a Point")

        return all([
            self.x == another.x,
            self.y == another.y,
            self.z == another.z
            ]
        )

    def __ne__(self,
        another: 'Point'
    ) -> bool:
        """
        Return False if objects are equal.

        Example:
          >>> Point(1., 1., 1.) != Point(0., 0., 0.)
          True
          >>> Point(1., 1., 1.) != Point(1, 1, 1)
          True
        """

        return not (self == another)

    def a(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Return the individual values of the point.

        :return: double array of x, y, z values

        Examples:
          >>> Point(4, 3, 7).a()
          (4.0, 3.0, 7.0)
        """

        return self.x, self.y, self.z

    def __add__(self, another: 'Point') -> 'Point':
        """
        Sum of two points.

        :param another: the point to add
        :type another: Point
        :return: the sum of the two points
        :rtype: Point
        :raise: Exception

        Example:
          >>> Point(1, 0, 0) + Point(0, 1, 1)
          Point(1.0000, 1.0000, 1.0000)
          >>> Point(1, 1, 1) + Point(-1, -1, -1)
          Point(0.0000, 0.0000, 0.0000)
        """

        check_type(another, "Second point", Point)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point(
            x=x0+x1,
            y=y0+y1,
            z=z0+z1
        )

    def __sub__(self,
        another: 'Point'
    ) -> 'Point':
        """Subtract two points.

        :param another: the point to subtract
        :type another: Point
        :return: the difference between the two points
        :rtype: Point
        :raise: Exception

        Example:
          >>> Point(1., 1., 1.) - Point(1., 1., 1.)
          Point(0.0000, 0.0000, 0.0000)
          >>> Point(1., 1., 3.) - Point(1., 1., 2.2)
          Point(0.0000, 0.0000, 0.8000)
        """

        check_type(another, "Second point", Point)

        x0, y0, z0 = self
        x1, y1, z1 = another

        return Point(
            x=x0 - x1,
            y=y0 - y1,
            z=z0 - z1
        )

    def clone(self) -> 'Point':
        """
        Clone a point.

        :return: a new point.
        :rtype: Point.
        """

        return Point(*self.a())

    def toXYZ(self) -> Tuple[numbers.Real, numbers.Real, numbers.Real]:
        """
        Returns the spatial components as a tuple of three values.

        :return: the spatial components (x, y, z).
        :rtype: a tuple of three floats.

        Examples:
          >>> Point(1, 0, 3).toXYZ()
          (1.0, 0.0, 3.0)
        """

        return self.x, self.y, self.z

    def toArray(self) -> np.ndarray:
        """
        Return a Numpy array representing the point values.

        :return: Numpy array

        Examples:
          >>> np.allclose(Point(1, 2, 3).toArray(), np.array([ 1., 2., 3., 0.]))
          True
        """

        return np.asarray(self.toXYZ())

    def pXY(self) -> 'Point':
        """
        Projection on the x-y plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXY()
          Point(2.0000, 3.0000, 0.0000)
        """

        return Point(self.x, self.y, 0.0)

    def pXZ(self) -> 'Point':
        """
        Projection on the x-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pXZ()
          Point(2.0000, 0.0000, 4.0000)
        """

        return Point(self.x, 0.0, self.z)

    def pYZ(self) -> 'Point':
        """
        Projection on the y-z plane

        :return: projected object instance

        Examples:
          >>> Point(2, 3, 4).pYZ()
          Point(0.0000, 3.0000, 4.0000)
        """

        return Point(0.0, self.y, self.z)

    def deltaX(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between x components of two Point Instances.

        :return: x coordinates difference value.
        :rtype: optional numbers.Real.
        :raise: Exception

        Examples:
          >>> Point(1, 2, 3).deltaX(Point(4, 7, 1))
          3.0
        """

        check_crs(self, another)

        return another.x - self.x

    def deltaY(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between y components of two Point Instances.

        :return: y coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3).deltaY(Point(4, 7, 1))
          5.0
        """

        check_crs(self, another)

        return another.y - self.y

    def deltaZ(self,
        another: 'Point'
    ) -> Optional[numbers.Real]:
        """
        Delta between z components of two Point Instances.

        :return: z coordinates difference value.
        :rtype: optional numbers.Real.

        Examples:
          >>> Point(1, 2, 3).deltaZ(Point(4, 7, 1))
          -2.0
        """

        check_crs(self, another)

        return another.z - self.z

    def dist3DWith(self,
        another: 'Point'
    ) -> numbers.Real:
        """
        Calculate Euclidean spatial distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point(1., 1., 1.).dist3DWith(Point(4., 5., 1))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
          >>> Point(1, 1, 1).dist3DWith(Point(4, 5, 1))
          5.0
        """

        check_type(another, "Point", Point)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2 + (self.z - another.z) ** 2)

    def dist2DWith(self,
        another: 'Point'
    ) -> numbers.Real:
        """
        Calculate horizontal (2D) distance between two points.
        TODO: consider case of polar CRS

        :param another: another Point instance.
        :type another: Point.
        :return: the 2D distance (when the two points have the same CRS).
        :rtype: numbers.Real.
        :raise: Exception.

        Examples:
          >>> Point(1., 1., 1.).dist2DWith(Point(4., 5., 7.))
          5.0
        """

        check_type(another, "Second point", Point)

        return sqrt((self.x - another.x) ** 2 + (self.y - another.y) ** 2)

    def scale(self,
        scale_factor: numbers.Real
    ) -> 'Point':
        """
        Create a scaled object.
        Note: it does not make sense for polar coordinates.
        TODO: manage polar coordinates cases OR deprecate and remove - after dependency check.

        Example;
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
          >>> Point(1, 0, 1).scale(2.5)
          Point(2.5000, 0.0000, 2.5000)
        """

        x, y, z = self.x * scale_factor, self.y * scale_factor, self.z * scale_factor
        return Point(x, y, z)

    def invert(self) -> 'Point':
        """
        Create a new object with inverted direction.
        Note: it depends on scale method, that could be deprecated/removed.

        Examples:
          >>> Point(1, 1, 1).invert()
          Point(-1.0000, -1.0000, -1.0000)
          >>> Point(2, -1, 4).invert()
          Point(-2.0000, 1.0000, -4.0000)
        """

        return self.scale(-1)

    def reflect_vertical(self) -> 'Point':
        """
        Reflect a point along a vertical axis.

        :return: reflected point.
        :rtype: Point

        Examples:
          >>> Point(1,1,1).reflect_vertical()
          Point(-1.0000, -1.0000, 1.0000)
        """

        x, y, z = self

        return Point(
            x=-x,
            y=-y,
            z=z
        )

    def isCoinc3D(self,
                  another: 'Point',
                  tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
                  ) -> bool:
        """
        Check spatial coincidence of two points

        :param another: the point to compare.
        :type another: Point.
        :param tolerance: the maximum allowed distance between the two points.
        :type tolerance: numbers.Real.
        :return: whether the two points are coincident.
        :rtype: bool.
        :raise: Exception.

        Example:
          >>> Point(1., 0., -1.).isCoinc3D(Point(1., 1.5, -1.))
          False
          >>> Point(1., 0., 0.).isCoinc3D(Point(1., 0., 0.))
          True
        """

        check_type(another, "Second point", Point)

        return self.dist3DWith(another) <= tolerance

    def already_present(self,
        pt_list: List['Point'],
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> Optional[bool]:
        """
        Determines if a point is already in a given point list, using an optional distance separation,

        :param pt_list: list of points. May be empty.
        :type pt_list: List of Points.
        :param tolerance: optional maximum distance between near-coincident point pair.
        :type tolerance: numbers.Real.
        :return: True if already present, False otherwise.
        :rtype: optional boolean.
        """

        for pt in pt_list:
            if self.isCoinc3D(pt, tolerance=tolerance):
                return True
        return False

    def shift(self,
        sx: numbers.Real,
        sy: numbers.Real,
        sz: numbers.Real
    ) -> Optional['Point']:
        """
        Create a new object shifted by given amount from the self instance.

        Example:
          >>> Point(1, 1, 1).shift(0.5, 1., 1.5)
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shift(0.5, 1., 1.5)
          Point(1.5000, 3.0000, 0.5000)
       """

        return Point(self.x + sx, self.y + sy, self.z + sz)

    def shiftByVect(self,
        v: Vect
    ) -> 'Point':
        """
        Create a new point shifted from the self instance by given vector.

        :param v: the shift vector.
        :type v: Vect.
        :return: the shifted point.
        :rtype: Point.
        :raise: Exception

        Example:
          >>> Point(1, 1, 1).shiftByVect(Vect(0.5, 1., 1.5))
          Point(1.5000, 2.0000, 2.5000)
          >>> Point(1, 2, -1).shiftByVect(Vect(0.5, 1., 1.5))
          Point(1.5000, 3.0000, 0.5000)
       """

        x, y, z = self

        sx, sy, sz = v.toXYZ()

        return Point(x + sx, y + sy, z + sz)

    def asVect(self) -> 'Vect':
        """
        Create a vector based on the point coordinates

        Example:
          >>> Point(1, 1, 0).asVect()
          Vect(1.0000, 1.0000, 0.0000)
          >>> Point(0.2, 1, 6).asVect()
          Vect(0.2000, 1.0000, 6.0000)
        """

        return Vect(self.x, self.y, self.z)

    def rotate(self,
        rotation_axis: 'RotationAxis',
        center_point: 'Point' = None
        ) -> 'Point':
        """
        Rotates a point.
        :param rotation_axis:
        :param center_point:
        :return: the rotated point
        :rtype: Point

        Examples:
          >>> pt = Point(0,0,1)
          >>> rot_axis = RotationAxis(0,0,90)
          >>> center_pt = Point(0,0,0.5)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(0.5000, 0.0000, 0.5000)
          >>> center_pt = Point(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(0.0000, 0.0000, 1.0000)
          >>> center_pt = Point(0, 0, 2)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-1.0000, 0.0000, 2.0000)
          >>> rot_axis = RotationAxis(0,0,180)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-0.0000, 0.0000, 3.0000)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(0.0000, 0.0000, -1.0000)
          >>> pt = Point(1,1,1,5)
          >>> rot_axis = RotationAxis(0,90,90)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(1.0000, -1.0000, 1.0000)
          >>> rot_axis = RotationAxis(0,90,180)
          >>> pt.rotate(rotation_axis=rot_axis)
          Point(-1.0000, -1.0000, 1.0000)
          >>> center_pt = Point(1,1,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(1.0000, 1.0000, 1.0000)
          >>> center_pt = Point(2,2,10)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(3.0000, 3.0000, 1.0000)
          >>> pt = Point(1, 1, 2)
          >>> rot_axis = RotationAxis(135, 0, 180)
          >>> center_pt = Point(0,0,1)
          >>> pt.rotate(rotation_axis=rot_axis, center_point=center_pt)
          Point(-1.0000, -1.0000, 0.0000)
        """

        if not center_point:

            center_point = Point(
                x=0.0,
                y=0.0,
                z=0.0
            )

        check_type(center_point, "Center point", Point)

        p_diff = self - center_point

        p_vect = p_diff.asVect()

        rot_vect = rotVectByAxis(
            v=p_vect,
            rot_axis=rotation_axis
        )

        x, y, z, epsg_cd = rot_vect

        rot_pt = Point(
            x=x,
            y=y,
            z=z
        )

        transl_pt = center_point + rot_pt

        return transl_pt

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE
    ):
        """
        Creates a random point.

        :return: random point
        :rtype: Point
        """

        vals = [random.uniform(lower_boundary, upper_boundary) for _ in range(3)]
        return cls(*vals)


def pack_to_points(
    xs: array,
    ys: array,
    zs: Optional[array] = None,
) -> List[Point]:
    # Side effects: None
    """
    Create a list of points given a set
    of input arrays.

    :param xs: array of x values
    :type xs: array
    :param ys: array of y values
    :type ys: array
    :param zs: optional array of z values
    :type zs: Optional[array]
    :return: a list of Point instances
    :rtype: List[Point]
    """

    if zs is None:
        zs = [0.0] * len(xs)
    pts = []
    for x, y, z, t in zip(xs, ys, zs):
        pts.append(
            Point(
                x,
                y,
                z
            )
        )

    return pts


class Segment:
    """
    Segment is a geometric object defined by the straight line between
    two vertices.
    """

    def __init__(self,
                 start_pt: Point,
                 end_pt: Point):
        """
        Creates a segment instance provided the two points have the same CRS code.

        :param start_pt: the start point.
        :type: Point.
        :param end_pt: the end point.
        :type end_pt: Point.
        :return: the new segment instance if both points have the same crs.
        :raises: CRSCodeException.
        """

        check_type(start_pt, "Start point", Point)

        check_type(end_pt, "End point", Point)

        if start_pt.dist3DWith(end_pt) == 0.0:
            raise Exception("Source points cannot be coincident")

        self._start_pt = start_pt.clone()
        self._end_pt = end_pt.clone()

    @classmethod
    def fromVector(cls,
                   point: Point,
                   dir_vector: Vect):

        check_type(point, "Input point", Point)
        check_type(dir_vector, "Directional vector", Vect)

        start_pt = point
        end_pt = start_pt.shiftByVect(dir_vector)

        return cls(
            start_pt=start_pt,
            end_pt=end_pt
        )

    def __repr__(self) -> str:
        """
        Represents a Segment instance.

        :return: the Segment representation.
        :rtype: str.
        """

        return "Segment(start_pt={}, end_pt={})".format(
            self.start_pt,
            self.end_pt
        )

    @property
    def start_pt(self) -> Point:

        return self._start_pt

    @property
    def end_pt(self) -> Point:

        return self._end_pt

    def __iter__(self):
        """
        Return the elements of a Segment, i.e., start and end point.
        """

        return (i for i in [self.start_pt, self.end_pt])

    def clone(self) -> 'Segment':

        return Segment(self._start_pt, self._end_pt)

    def increasing_x(self) -> 'Segment':

        if self.end_pt.x < self.start_pt.x:
            return Segment(self.end_pt, self.start_pt)
        else:
            return self.clone()

    def x_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.x < self.end_pt.x:
            return self.start_pt.x, self.end_pt.x
        else:
            return self.end_pt.x, self.start_pt.x

    def y_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.y < self.end_pt.y:
            return self.start_pt.y, self.end_pt.y
        else:
            return self.end_pt.y, self.start_pt.y

    def z_range(self) -> Tuple[numbers.Real, numbers.Real]:

        if self.start_pt.z < self.end_pt.z:
            return self.start_pt.z, self.end_pt.z
        else:
            return self.end_pt.z, self.start_pt.z

    def delta_x(self) -> numbers.Real:

        return self.end_pt.x - self.start_pt.x

    def delta_y(self) -> numbers.Real:

        return self.end_pt.y - self.start_pt.y

    def delta_z(self) -> numbers.Real:
        """
        Z delta between segment end point and start point.

        :return: numbers.Real.
        """

        return self.end_pt.z - self.start_pt.z

    def length2D(self) -> numbers.Real:

        return self.start_pt.dist2DWith(self.end_pt)

    def length3D(self) -> numbers.Real:

        return self.start_pt.dist3DWith(self.end_pt)

    def deltaZS(self) -> Optional[numbers.Real]:
        """
        Calculates the delta z - delta s ratio of a segment.

        :return: optional numbers.Real.
        """

        len2d = self.length2D()

        if len2d == 0.0:
            return None

        return self.delta_z() / len2d

    def slope_rad(self) -> Optional[numbers.Real]:
        """
        Calculates the slope in radians of the segment.
        Positive is downward point, negative upward pointing.

        :return: optional numbers.Real.
        """

        delta_zs = self.deltaZS()

        if delta_zs is None:
            return None
        else:
            return - math.atan(delta_zs)

    def vector(self) -> Vect:

        return Vect(self.delta_x(),
                    self.delta_y(),
                    self.delta_z()
        )

    def antivector(self) -> Vect:
        """
        Returns the vector pointing from the segment end to the segment start.

        :return: the vector pointing from the segment end to the segment start.
        :rtype: Vect.
        """

        return self.vector().invert()

    def contains_pt(self,
        pt: Point
    ) -> bool:
        """
        Checks whether a point is contained in a segment.

        :param pt: the point for which to check containement.
        :return: bool.
        :raise: Exception.

        Examples:
          >>> segment = Segment(Point(0, 0, 0), Point(1, 0, 0))
          >>> segment.contains_pt(Point(0, 0, 0))
          True
          >>> segment.contains_pt(Point(1, 0, 0))
          True
          >>> segment.contains_pt(Point(0.5, 0, 0))
          True
          >>> segment.contains_pt(Point(0.5, 0.00001, 0))
          False
          >>> segment.contains_pt(Point(0.5, 0, 0.00001))
          False
          >>> segment.contains_pt(Point(1.00001, 0, 0))
          False
          >>> segment.contains_pt(Point(0.000001, 0, 0))
          True
          >>> segment.contains_pt(Point(-0.000001, 0, 0))
          False
          >>> segment.contains_pt(Point(0.5, 1000, 1000))
          False
          >>> segment = Segment(Point(0, 0, 0), Point(0, 1, 0))
          >>> segment.contains_pt(Point(0, 0, 0))
          True
          >>> segment.contains_pt(Point(0, 0.5, 0))
          True
          >>> segment.contains_pt(Point(0, 1, 0))
          True
          >>> segment.contains_pt(Point(0, 1.5, 0))
          False
          >>> segment = Segment(Point(0, 0, 0), Point(1, 1, 1))
          >>> segment.contains_pt(Point(0.5, 0.5, 0.5))
          True
          >>> segment.contains_pt(Point(1, 1, 1))
          True
          >>> segment = Segment(Point(1,2,3), Point(9,8,2))
          >>> segment.contains_pt(segment.pointAt(0.745))
          True
          >>> segment.contains_pt(segment.pointAt(1.745))
          False
          >>> segment.contains_pt(segment.pointAt(-0.745))
          False
          >>> segment.contains_pt(segment.pointAt(0))
          True
        """

        check_type(pt, "Point", Point)

        segment_length = self.length3D()
        length_startpt_pt = self.start_pt.dist3DWith(pt)
        length_endpt_pt = self.end_pt.dist3DWith(pt)

        return areClose(
            a=segment_length,
            b=length_startpt_pt + length_endpt_pt
        )

    def pointAt(self,
        scale_factor: numbers.Real
    ) -> Point:
        """
        Returns a point aligned with the segment
        and lying at given scale factor, where 1 is segment length
        ans 0 is segment start.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point

        Examples:
          >>> s = Segment(Point(0,0,0), Point(1,0,0))
          >>> s.pointAt(0)
          Point(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point(0.5000, 0.0000, 0.0000)
          >>> s.pointAt(1)
          Point(1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-1)
          Point(-1.0000, 0.0000, 0.0000)
          >>> s.pointAt(-2)
          Point(-2.0000, 0.0000, 0.0000)
          >>> s.pointAt(2)
          Point(2.0000, 0.0000, 0.0000)
          >>> s = Segment(Point(0,0,0), Point(0,0,1))
          >>> s.pointAt(0)
          Point(0.0000, 0.0000, 0.0000)
          >>> s.pointAt(0.5)
          Point(0.0000, 0.0000, 0.5000)
          >>> s.pointAt(1)
          Point(0.0000, 0.0000, 1.0000)
          >>> s.pointAt(-1)
          Point(0.0000.0000, 0.0000, -1)
          >>> s.pointAt(-2)
          Point(0.0000, 0.0000, -2.0000)
          >>> s.pointAt(2)
          Point(0.0000, 0.0000, 2.0000)
          >>> s = Segment(Point(0,0,0), Point(1,1,1))
          >>> s.pointAt(0.5)
          Point(0.5000, 0.5000, 0.5000)
          >>> s = Segment(Point(0,0,0), Point(4,0,0))
          >>> s.pointAt(7.5)
          Point(30.0000, 0.0000, 0.0000)
        """

        dx = self.delta_x() * scale_factor
        dy = self.delta_y() * scale_factor
        dz = self.delta_z() * scale_factor

        return Point(
            x=self.start_pt.x + dx,
            y=self.start_pt.y + dy,
            z=self.start_pt.z + dz
        )

    def pointProjection(self,
        point: Point
    ) -> Point:
        """
        Return the point projection on the segment.

        Examples:
          >>> s = Segment(start_pt=Point(0,0,0), end_pt=Point(1,0,0))
          >>> p = Point(0.5, 1, 4)
          >>> s.pointProjection(p)
          Point(0.5000, 0.0000, 0.0000)
          >>> s = Segment(start_pt=Point(0,0,0), end_pt=Point(4,0,0))
          >>> p = Point(7.5, 19.2, -14.72)
          >>> s.pointProjection(p)
          Point(7.5000, 0.0000, 0.0000)
        """

        check_type(point, "Input point", Point)

        check_crs(self, point)

        other_segment = Segment(
            self.start_pt,
            point
        )
        
        scale_factor = self.vector().scalarProjection(other_segment.vector()) / self.length3D()
        return self.pointAt(scale_factor)

    def pointDistance(self,
        point: Point
    ) -> numbers.Real:
        """
        Returns the point distance to the segment.

        :param point: the point to calculate the distance with
        :type point: Point
        :return: the distance of the point to the segment
        :rtype: numbers.Real

        Examples:
          >>> s = Segment(Point(0,0,0), Point(0,0,4))
          >>> s.pointDistance(Point(-17.2, 0.0, -49,3))
          17.2
          >>> s.pointDistance(Point(-17.2, 1.22, -49,3))
          17.24321315764553
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        point_projection = self.pointProjection(point)

        return point.dist3DWith(point_projection)

    def pointS(self,
        point: Point
    ) -> Optional[numbers.Real]:
        """
        Calculates the optional distance of the point along the segment.
        A zero value is for a point coinciding with the start point.
        Returns None if the point is not contained in the segment.

        :param point: the point to calculate the optional distance in the segment.
        :type point: Point
        :return: the the optional distance of the point along the segment.
        """

        check_type(point, "Input point", Point)

        #check_crs(self, point)

        if not self.contains_pt(point):
            return None

        return self.start_pt.dist3DWith(point)

    def scale(self,
        scale_factor
    ) -> 'Segment':
        """
        Scale a segment by the given scale_factor.
        Start point does not change.

        :param scale_factor: the scale factor, where 1 is the segment length.
        :type scale_factor: numbers.Real
        :return: Point at scale factor
        :rtype: Point
        """

        end_pt = self.pointAt(scale_factor)

        return Segment(
            self.start_pt,
            end_pt)

    def vertical_plane(self) -> Optional[CPlane]:
        """
        Returns the vertical Cartesian plane containing the segment.

        :return: the vertical Cartesian plane containing the segment.
        :rtype: Optional[CPlane].
        """

        if self.length2D() == 0.0:
            return None

        # arbitrary point on the same vertical as end point

        section_final_pt_up = self.end_pt.shift(
            sx=0.0,
            sy=0.0,
            sz=1000.0)

        return CPlane.fromPoints(
            pt1=self.start_pt,
            pt2=self.end_pt,
            pt3=section_final_pt_up)

    def same_start(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the two segments have the same start point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(0,0,0), Point(0,1,0))
          >>> s1.same_start(s2)
          True
        """

        return self.start_pt.isCoinc3D(
            another=another.start_pt,
            tolerance=tol
        )

    def same_end(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the two segments have the same end point.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the two segments have the same end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(2,0,0), Point(1,0,0))
          >>> s1.same_end(s2)
          True
        """

        return self.end_pt.isCoinc3D(
            another=another.end_pt,
            tolerance=tol)

    def conn_to_other(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the first segment is sequentially connected to the second one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the first segment is sequentially connected to the second one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(1,0,0), Point(2,0,0))
          >>> s1.conn_to_other(s2)
          True
        """

        return self.end_pt.isCoinc3D(
            another=another.start_pt,
            tolerance=tol)

    def other_connected(self,
        another: 'Segment',
        tol: numbers.Real = 1e-12
    ) -> bool:
        """
        Check whether the second segment is sequentially connected to the first one.

        :param another: a segment to check for.
        :type another: Segment.
        :param tol: tolerance for distance between points.
        :type tol: numbers.Real.
        :return: whether the second segment is sequentially connected to the first one.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-1,0,0), Point(0,0,0))
          >>> s1.other_connected(s2)
          True
        """

        return another.end_pt.isCoinc3D(
            another=self.start_pt,
            tolerance=tol)

    def segment_start_in(self,
        another: 'Segment'
    ) -> bool:
        """
        Check whether the second segment contains the first segment start point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment start point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-0.5,0,0), Point(0.5,0,0))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment(Point(0,0,0), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          True
          >>> s1 = Segment(Point(0,1,0), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          False
          >>> s1 = Segment(Point(-1,-1,-1), Point(1,1,1))
          >>> s1.segment_start_in(s2)
          False
        """

        return another.contains_pt(self.start_pt)

    def segment_end_in(self,
        another: 'Segment'
    ) -> bool:
        """
        Check whether the second segment contains the first segment end point.

        :param another: a segment to check for.
        :type another: Segment.
        :return: whether the second segment contains the first segment end point.
        :rtype: bool.

        Examples:
          >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
          >>> s2 = Segment(Point(-0.5,0,0), Point(0.5,0,0))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment(Point(0,0,0), Point(1,1,1))
          >>> s1.segment_end_in(s2)
          False
          >>> s1 = Segment(Point(0,1,0), Point(1,1,1))
          >>> s2 = Segment(Point(1,1,1), Point(0.5,0,0))
          >>> s1.segment_end_in(s2)
          True
          >>> s1 = Segment(Point(-1,-1,3), Point(1,1,3))
          >>> s2 = Segment(Point(0,2,3), Point(2,0,3))
          >>> s1.segment_end_in(s2)
          True
        """

        return another.contains_pt(self.end_pt)

    def rotate(self,
        rotation_axis: 'RotationAxis',
        center_point: 'Point' = None
        ) -> 'Segment':
        """
        Rotates a segment.
        :param rotation_axis:
        :param center_point:
        :return: the rotated segment
        :rtype: Segment

        Examples:
        >>> seg = Segment(Point(0,0,0), Point(0,0,1))
        >>> rot_ax = RotationAxis(0, 0, 90)
        >>> seg.rotate(rot_ax)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
        >>> rot_ax = RotationAxis(0, 0, 180)
        >>> seg.rotate(rot_ax)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.0000.0000))
        >>> centr_pt = Point(0,0,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(-0.0000, 0.0000, 1.0000), end_pt=Point(0.0000, 0.0000, 0.0000))
        >>> seg = Segment(Point(0,0,0), Point(1,1,0))
        >>> centr_pt = Point(1,0,0)
        >>> rot_ax = RotationAxis(0, 90, 90)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(1.0000, 1.0000, 0.0000), end_pt=Point(2.0000, 0.0000, -0.0000))
        >>> seg = Segment(Point(1,1,1), Point(0,0,0))
        >>> rot_ax = RotationAxis(135, 0, 180)
        >>> centr_pt = Point(0.5,0.5,0.5)
        >>> seg.rotate(rotation_axis=rot_ax, center_point=centr_pt)
        Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 1.0000, 1.0000))
        """

        start_pt, end_pt = self

        rotated_start_pt = start_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        rotated_end_pt = end_pt.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        return Segment(
            start_pt=rotated_start_pt,
            end_pt=rotated_end_pt
        )

    @classmethod
    def random(cls,
        lower_boundary: float = -MAX_SCALAR_VALUE,
        upper_boundary: float = MAX_SCALAR_VALUE):
        """
        Creates a random segment.

        :return: random segment
        :rtype: Segment
        """

        return cls(
            start_pt=Point.random(lower_boundary, upper_boundary),
            end_pt=Point.random(lower_boundary, upper_boundary)
        )


def point_or_segment(
        point1: Point,
        point2: Point,
        tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Union[Point, Segment]:
    """
    Creates a point or segment based on the points distance.

    :param point1: first input point.
    :type point1: Point.
    :param point2: second input point.
    :type point2: Point.
    :param tol: distance tolerance between the two points.
    :type tol: numbers.Real.
    :return: point or segment based on their distance.
    :rtype: PointOrSegment.
    :raise: Exception.
    """

    check_type(point1, "First point", Point)
    check_type(point2, "Second point", Point)

    #check_crs(point1, point2)

    if point1.dist3DWith(point2) <= tol:
        return Points.fromPoints([point1, point2]).nanmean_point()
    else:
        return Segment(
            start_pt=point1,
            end_pt=point2
        )


def intersect_segments(
    segment1: Segment,
    segment2: Segment,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Point, Segment]]:
    """
    Determines the optional point or segment intersection between the segment pair.

    :param segment1: the first segment
    :type segment1: Segment
    :param segment2: the second segment
    :type segment2: Segment
    :param tol: the distance tolerance for collapsing a intersection segment into a point
    :type tol: numbers.Real
    :return: the optional point or segment intersection between the segment pair.
    :rtype: Optional[Union[Point, Segment]]

    Examples:
      >>> s2 = Segment(Point(0,0,0), Point(1,0,0))
      >>> s1 = Segment(Point(0,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(-1,0,0))
      >>> intersect_segments(s1, s2) is None
      True
      >>> s1 = Segment(Point(-2,0,0), Point(0,0,0))
      >>> intersect_segments(s1, s2)
      Point(0.0000, 0.0000, 0.0000)
      >>> s1 = Segment(Point(-2,0,0), Point(0.5,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(-2,0,0), Point(2,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0,0,0), Point(0.5,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(0.5000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(0.75,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(0.7500, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(1,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0.25,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(0,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.0000, 0.0000, 0.0000), end_pt=Point(1.0000, 0.0000, 0.0000))
      >>> s1 = Segment(Point(1,0,0), Point(1.25,0,0))
      >>> intersect_segments(s1, s2)
      Point(1.0000, 0.0000, 0.0000)
      >>> s2 = Segment(Point(0,0,0), Point(1,1,1))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(0.75,0.75,0.75))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.2500, 0.2500), end_pt=Point(0.7500, 0.7500, 0.7500))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,1.75,1.75))
      >>> intersect_segments(s1, s2)
      Segment(start_pt=Point(0.2500, 0.2500, 0.2500), end_pt=Point(1.0000, 1.0000, 1.0000))
      >>> s1 = Segment(Point(0.25,0.25,0.25), Point(1.75,0,1.75))
      >>> intersect_segments(s1, s2)
      Point(0.2500, 0.2500, 0.2500)
      >>> s1 = Segment(Point(0.25,1,0.25), Point(0.75,0.75,0.75))
      >>> intersect_segments(s1, s2)
      Point(0.7500, 0.7500, 0.7500)
      >>> s2 = Segment(Point(-1,-1,-1), Point(1,1,1))
      >>> s1 = Segment(Point(-1,1,1), Point(1,-1,-1))
      >>> intersect_segments(s1, s2)
      Point(-0.0000, 0.0000, 0.0000)
    """

    check_type(segment1, "First segment", Segment)
    check_type(segment2, "Second segment", Segment)

    #check_crs(segment1, segment2)

    s1_startpt_inside = segment1.segment_start_in(segment2)
    s2_startpt_inside = segment2.segment_start_in(segment1)

    s1_endpt_inside = segment1.segment_end_in(segment2)
    s2_endpt_inside = segment2.segment_end_in(segment1)

    elements = [s1_startpt_inside, s2_startpt_inside, s1_endpt_inside, s2_endpt_inside]

    if all(elements):
        return segment1.clone()

    if s1_startpt_inside and s1_endpt_inside:
        return segment1.clone()

    if s2_startpt_inside and s2_endpt_inside:
        return segment2.clone()

    if s1_startpt_inside and s2_startpt_inside:
        return point_or_segment(
            segment1.start_pt,
            segment2.start_pt,
            tol=tol
        )

    if s1_startpt_inside and s2_endpt_inside:
        return point_or_segment(
            segment1.start_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_startpt_inside:
        return point_or_segment(
            segment2.start_pt,
            segment1.end_pt,
            tol=tol
        )

    if s1_endpt_inside and s2_endpt_inside:
        return point_or_segment(
            segment1.end_pt,
            segment2.end_pt,
            tol=tol
        )

    if s1_startpt_inside:
        return segment1.start_pt.clone()

    if s1_endpt_inside:
        return segment1.end_pt.clone()

    if s2_startpt_inside:
        return segment2.start_pt.clone()

    if s2_endpt_inside:
        return segment2.end_pt.clone()

    cline1 = Segment.fromSegment(segment1)
    cline2 = Segment.fromSegment(segment2)

    shortest_segm_or_pt = cline1.shortest_segment_or_point(
        cline2,
        tol=tol
    )

    if not shortest_segm_or_pt:
        return None

    if not isinstance(shortest_segm_or_pt, Point):
        return None

    inters_pt = shortest_segm_or_pt

    if not segment1.contains_pt(inters_pt):
        return None

    if not segment2.contains_pt(inters_pt):
        return None

    return inters_pt


class Line:
    """
    A line.
    """

    def __init__(self,
                 pts: List[Point]):
        """

        """

        check_type(pts, "List", list)
        for el in pts:
            check_type(el, "Point", Point)

        self._pts = pts

    def pts(self):
        return self._pts

    def pt(self,
           ndx: numbers.Integral):
        """

        """

        return self._pts[ndx]

    def add_pt(self,
               pt: Point):

        self._pts.append(pt)

    def num_pts(self):
        return len(self._pts)

    def x_min(self):
        return min(map(lambda pt: pt.x, self._pts))

    def x_max(self):
        return max(map(lambda pt: pt.x, self._pts))

    def y_min(self):
        return min(map(lambda pt: pt.y, self._pts))

    def y_max(self):
        return max(map(lambda pt: pt.y, self._pts))

    def z_min(self):
        return min(map(lambda pt: pt.z, self._pts))

    def z_max(self):
        return max(map(lambda pt: pt.z, self._pts))

    def as_segments(self):
        """
        Convert to a list of segments.

        :return: list of Segment objects
        """

        pts_pairs = zip(self.pts()[:-1], self.pts()[1:])

        segments = [Segment(pt_a, pt_b) for (pt_a, pt_b) in pts_pairs]

        return segments

    '''
    def densify_2d_line(self, sample_distance) -> 'Points':
        """
        Densify a line into a new line instance,
        using the provided sample distance.
        Returned Line instance has coincident successive points removed.

        :param sample_distance: numbers.Real
        :return: Line instance
        """

        if sample_distance <= 0.0:
            raise Exception(f"Sample distance must be positive. {sample_distance} received")

        segments = self.as_segments()

        densified_line_list = [segment.densify2d_asLine(sample_distance) for segment in segments]

        densifyied_multiline = MultiLine(densified_line_list)

        densifyied_line = densifyied_multiline.to_line()

        densifyied_line_wo_coinc_pts = densifyied_line.remove_coincident_points()

        return densifyied_line_wo_coinc_pts
    '''

    def join(self, another) -> 'Points':
        """
        Joins together two lines and returns the join as a new line without point changes,
        with possible overlapping points
        and orientation mismatches between the two original lines
        """

        return Line(self.pts() + another.pts())

    def length_3d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).dist3DWith(self.pt(ndx + 1))
        return length

    '''
    def length_2d(self) -> numbers.Real:

        length = 0.0
        for ndx in range(self.num_pts() - 1):
            length += self.pt(ndx).dist2DWith(self.pt(ndx + 1))
        return length
    '''

    def step_delta_z(self) -> List[numbers.Real]:
        """
        Return the difference in elevation between consecutive points:
        z[ndx+1] - z[ndx]

        :return: a list of height differences.
        :rtype: list of floats.
        """

        delta_z = [0.0]

        for ndx in range(1, self.num_pts()):
            delta_z.append(self.pt(ndx).z - self.pt(ndx - 1).z)

        return delta_z

    def step_lengths_3d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 3D distances.
        It is the distance between a point and its previous one.
        The list has the same lenght as the source point list.

        :return: the individual 3D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self.pt(ndx).dist3DWith(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list

    '''
    def step_lengths_2d(self) -> List[numbers.Real]:
        """
        Returns the point-to-point 2D distances.
        It is the distance between a point and its previous one.
        The list has the same length as the source point list.

        :return: the individual 2D segment lengths.
        :rtype: list of floats.

        Examples:
        """

        step_length_list = [0.0]
        for ndx in range(1, self.num_pts()):
            length = self.pt(ndx).dist2DWith(self.pt(ndx - 1))
            step_length_list.append(length)

        return step_length_list
    '''

    def incremental_length_3d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 3D segment lengths.

        :return: accumulated 3D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_3d()))

    '''
    def incremental_length_2d(self) -> List[numbers.Real]:
        """
        Returns the accumulated 2D segment lengths.

        :return: accumulated 2D segment lenghts
        :rtype: list of floats.
        """

        return list(itertools.accumulate(self.step_lengths_2d()))
    '''

    def reversed(self) -> 'Points':
        """
        Return a Line instance with reversed point list.

        :return: a new Line instance.
        :rtype: Line.
        """

        pts = [pt.clone() for pt in self.pts()]
        pts.reverse()

        return Line(
            pts=pts
        )

    def slopes_degr(self) -> List[Optional[numbers.Real]]:
        """
        Calculates the slopes (in degrees) of each Line segment.
        The first value is the slope of the first segment.
        The last value, always None, is the slope of the segment starting at the last point.
        The number of elements is equal to the number of points in the Line.

        :return: list of slopes (degrees).
        :rtype: List[Optional[numbers.Real]].
        """

        lSlopes = []

        segments = self.as_segments()
        for segment in segments:
            vector = segment.vector()
            lSlopes.append(-vector.slope_degr())  # minus because vector convention is positive downward

        lSlopes.append(None)  # None refers to the slope of the Segment starting with the last point

        return lSlopes

    def slopes_stats(self) -> Dict:
        """
        Returns the line directional slope statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        """

        return get_statistics(self.slopes_degr())

    def abs_slopes_degr(self) -> List[Optional[numbers.Real]]:

        return [abs(val) for val in self.slopes_degr()]

    def abs_slopes_stats(self) -> Dict:
        """
        Returns the line absolute slopes statistics.

        :return: the statistics parameters: min, max, mean, var, std.
        :rtype: Dictionary.
        """

        return get_statistics(self.abs_slopes_degr())

    def extremes_distance_3d(self) -> numbers.Real:
        """
        Calculate the 3D distance between start and end points.

        :return: the 3D distance between start and end points
        :rtype: numbers.Real
        """

        return self.end_pt().dist3DWith(self.start_pt())

    '''
    def extremes_distance_2d(self) -> numbers.Real:
        """
        Calculate the 2D distance between start and end points.

        :return: the 2D distance between start and end points
        """

        return self.end_pt().dist2DWith(self.start_pt())
    '''

    def isClosed_3d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 3D-closed.

        :param tolerance: the tolerance for considering the line closed
        :type tolerance: numbers.Real
        :return: whether the line is to be considered 3D-closed
        :rtype: bool
        """

        return self.end_pt().isCoinc3D(self.start_pt(), tolerance=tolerance)

    '''
    def isClosed_2d(self,
        tolerance: numbers.Real = MIN_SEPARATION_THRESHOLD
    ) -> bool:
        """
        Determine if the line is 2D-closed.

        :param tolerance: the tolerance for considering the line closed
        :return: whether the line is to be considered 2D-closed
        """

        return self.end_pt().isCoinc2D(self.start_pt(), tolerance=tolerance)
    '''

    def walk_backward(self) -> 'Points':
        """
        Create a new line by walking the line backward from the last point up to the first and thus closing it.

        :return: a closed line with zero area
        :rtype: 'Line'
        """

        return Line(self.pts() + self.reversed().pts()[1:])

    def clone(self) -> 'Points':
        """
        Clone a line.

        :return: the cloned line
        :rtype: Points
        """

        return Line(self.pts())

    '''
    def close_2d(self) -> 'Points':
        """
        Return a line that is 2D-closed.

        :return: a 2D-closed line
        :rtype: Points
        """

        line = self.clone()

        if not line.isClosed_2d():

            line.add_pt(line.start_pt())

        return line
    '''

    def close_3d(self) -> 'Points':
        """
        Return a line that is 3D-closed.

        :return: a 3D-closed line
        :rtype: Points
        """

        line = self.clone()

        if not line.isClosed_3d():

            line.add_pt(line.start_pt())

        return line

    def remove_coincident_points(self) -> Optional['Points']:
        """
        Remove coincident successive points

        :return: Line instance
        :rtype: Optional[Points]
        """

        if self.num_pts() == 0:
            return

        new_line = Line(
            pts=[self.pt(0)]
        )

        for ndx in range(1, self.num_pts()):
            if not self.pt(ndx).isCoinc3D(new_line.pt(-1)):
                new_line.add_pt(self.pt(ndx))

        return new_line


class JoinTypes(Enum):
    """
    Enumeration for Line and Segment type.
    """

    START_START = 1  # start point coincident with start point
    START_END   = 2  # start point coincident with end point
    END_START   = 3  # end point coincident with start point
    END_END     = 4  # end point coincident with end point


def analizeJoins(first: Union[Line, Segment], second: Union[Line, Segment]) -> List[Optional[JoinTypes]]:
    """
    Analyze join types between two lines/segments.

    :param first: a line or segment.
    :type first: Line or Segment.
    :param second: a line or segment.
    :param second: Line or Segment.
    :return: a list of join types.
    :rtype: List[Optional[JoinTypes]].

    Examples:
      >>> first = Segment(Point(x=0,y=0), Point(x=1,y=0))
      >>> second = Segment(Point(x=1,y=0), Point(x=0,y=0))
      >>> analizeJoins(first, second)
      [<JoinTypes.START_END: 2>, <JoinTypes.END_START: 3>]
      >>> first = Segment(Point(x=0,y=0), Point(x=1,y=0))
      >>> second = Segment(Point(x=2,y=0), Point(x=3,y=0))
      >>> analizeJoins(first, second)
      []
    """

    join_types = []

    if first.start_pt.isCoinc3D(second.start_pt):
        join_types.append(JoinTypes.START_START)

    if first.start_pt.isCoinc3D(second.end_pt):
        join_types.append(JoinTypes.START_END)

    if first.end_pt.isCoinc3D(second.start_pt):
        join_types.append(JoinTypes.END_START)

    if first.end_pt.isCoinc3D(second.end_pt):
        join_types.append(JoinTypes.END_END)

    return join_types


def shortest_segment_or_point(
    first_segment: Segment,
    second_segment: Segment,
    tol: numbers.Real = PRACTICAL_MIN_DIST
) -> Optional[Union[Segment, Point]]:

    """
    Calculates the optional shortest segment - or the intersection point - between two lines represented by two segments.

    Adapted from:
        http://paulbourke.net/geometry/pointlineplane/

    C code from:
        http://paulbourke.net/geometry/pointlineplane/lineline.c
[
    typedef struct {
    double x,y,z;
    } XYZ;

    /*
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. Calculate also the values of mua and mub where
      Pa = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - P3)
    Return FALSE if no solution exists.
    */
    int LineLineIntersect(
    XYZ p1,XYZ p2,XYZ p3,XYZ p4,XYZ *pa,XYZ *pb,
    double *mua, double *mub)
    {
    XYZ p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;

    p13.x = p1.x - p3.x;
    p13.y = p1.y - p3.y;
    p13.z = p1.z - p3.z;
    p43.x = p4.x - p3.x;
    p43.y = p4.y - p3.y;
    p43.z = p4.z - p3.z;
    if (ABS(p43.x) < EPS && ABS(p43.y) < EPS && ABS(p43.z) < EPS)
      return(FALSE);
    p21.x = p2.x - p1.x;
    p21.y = p2.y - p1.y;
    p21.z = p2.z - p1.z;
    if (ABS(p21.x) < EPS && ABS(p21.y) < EPS && ABS(p21.z) < EPS)
      return(FALSE);

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

    denom = d2121 * d4343 - d4321 * d4321;
    if (ABS(denom) < EPS)
      return(FALSE);
    numer = d1343 * d4321 - d1321 * d4343;

    *mua = numer / denom;
    *mub = (d1343 + d4321 * (*mua)) / d4343;

    pa->x = p1.x + *mua * p21.x;
    pa->y = p1.y + *mua * p21.y;
    pa->z = p1.z + *mua * p21.z;
    pb->x = p3.x + *mub * p43.x;
    pb->y = p3.y + *mub * p43.y;
    pb->z = p3.z + *mub * p43.z;

    return(TRUE);
    }

    :param first_segment: the first segment
    :param second_segment: the second segment
    :param tol: tolerance value for collapsing a segment into the midpoint.
    :return: the optional shortest segment or an intersection point.
    """

    check_type(second_segment, "Second Cartesian line", Segment)

    p1 = first_segment.start_pt
    p2 = first_segment.end_pt

    p3 = second_segment.start_pt
    p4 = second_segment.end_pt

    p13 = Point(
        x=p1.x - p3.x,
        y=p1.y - p3.y,
        z=p1.z - p3.z
    )

    p43 = Point(
        x=p4.x - p3.x,
        y=p4.y - p3.y,
        z=p4.z - p3.z
    )

    if p43.asVect().isAlmostZero:
        return None

    p21 = Point(
        x=p2.x - p1.x,
        y=p2.y - p1.y,
        z=p2.z - p1.z,
    )

    if p21.asVect().isAlmostZero:
        return None

    d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z
    d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z
    d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z
    d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z
    d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z

    denom = d2121 * d4343 - d4321 * d4321

    if fabs(denom) < MIN_SCALAR_VALUE:
        return None

    numer = d1343 * d4321 - d1321 * d4343

    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    pa = Point(
        x=p1.x + mua * p21.x,
        y=p1.y + mua * p21.y,
        z=p1.z + mua * p21.z
    )

    pb = Point(
        x=p3.x + mub * p43.x,
        y=p3.y + mub * p43.y,
        z=p3.z + mub * p43.z
    )

    intersection = point_or_segment(
        point1=pa,
        point2=pb,
        tol=tol
    )

    return intersection




