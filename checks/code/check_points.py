# coding: utf-8

# # Check point collections

from pygsf.geometries.shapes.space3d import *


def check_points():

    pts = [
        Point3D(0, 0, 0),
        Point3D(1, 1, 1),
        Point3D(2, 2, 2)
    ]

    point_coll = Points3D(pts)
    print(point_coll.nanmean_point())


if __name__ == '__main__':

    check_points()