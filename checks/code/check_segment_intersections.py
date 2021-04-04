# coding: utf-8

# # Check segment intersections
from pygsf.orientations.orientations import RotationAxis
from pygsf.geometries.shapes.space3d import *


def check_segment_intersections(n=100):

    for i in range(n):

        random_segment = Segment3D.random()

        random_interval_value = random.uniform(0, 1)

        rotation_axis = RotationAxis.randomNaive()

        center_point = random_segment.pointAt(random_interval_value)

        rotated_segment = random_segment.rotate(
            rotation_axis=rotation_axis,
            center_point=center_point
        )

        pair = SegmentPair(
            random_segment,
            rotated_segment
        )

        intersection = pair.intersect(tol=1)

        if not intersection:
            print("\nCase None: {}".format(i))
            print(random_segment)
            print(random_segment.length())
            print(random_interval_value)
            print(rotation_axis)
            print(center_point)
            print(rotated_segment)
            print(rotated_segment.length())
            print(random_segment.contains_pt(center_point))
            print(rotated_segment.contains_pt(center_point))
            connection_segment = shortest_segment_or_point(random_segment, rotated_segment)
            print(connection_segment)
            print(connection_segment.length())

        elif isinstance(intersection, Segment3D):
            print("\nCase Segment: {}".format(i))
            print(random_segment)
            print(random_segment.length())
            print(random_interval_value)
            print(rotation_axis)
            print(center_point)
            print(rotated_segment)
            print(rotated_segment.length())
            print(random_segment.contains_pt(center_point))
            print(rotated_segment.contains_pt(center_point))
            print(intersection)
            print(intersection.length())

        elif not intersection.is_coincident(center_point, tolerance=1):
            print("\nCase not coincident: {}".format(i))
            print(random_segment)
            print(random_segment.length())
            print(random_interval_value)
            print(rotation_axis)
            print(center_point)
            print(rotated_segment)
            print(rotated_segment.length())
            print(random_segment.contains_pt(center_point))
            print(rotated_segment.contains_pt(center_point))
            print(intersection)
            print(random_segment.contains_pt(intersection))
            print(rotated_segment.contains_pt(intersection))
            print(center_point.distance(intersection))

        else:
            pass


if __name__ == "__main__":

    check_segment_intersections(10000)
