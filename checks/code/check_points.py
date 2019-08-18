# coding: utf-8

# # Check point collections

import math
import numpy as np

from pygsf.spatial.vectorial.geometries import *


def check_points():

    pts = [
        Point(0,0,0),
        Point(1,1,1),
        Point(2,2,2)
    ]

    point_coll = Points(pts)
    print(point_coll.nanmean_point())


if __name__ == '__main__':

    check_points()