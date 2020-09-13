# -*- coding: utf-8 -*-


import unittest

import pygsf.orientations.orientations
import pygsf.geometries.space3d.shapes
import pygsf.orientations.direct_utils
from pygsf.geometries.space3d.shapes import *


class TestOrientations(unittest.TestCase):

    def setUp(self):

        pass

    def test_direct_general(self):
        """
        Check expected OrienM results for downward dip.
        """

        assert pygsf.orientations.orientations.Direct.fromAzPl(90, 90).isDownward
        assert pygsf.orientations.orientations.Direct.fromAzPl(90, -45).isUpward
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(90, 90).asVersor().z, -1.0)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(90, -90).asVersor().z, 1.0)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(0, 90).upward().asVersor().z, 1.0)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(0, -90).downward().asVersor().z, -1.0)

    def test_direct_angle(self):

        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(90, 45).angle(
            pygsf.orientations.orientations.Direct.fromAzPl(90, 55)), 10.)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(90, 45).angle(
            pygsf.orientations.orientations.Direct.fromAzPl(270, 10)), 125.)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(90, 90).angle(
            pygsf.orientations.orientations.Direct.fromAzPl(135, 90)), 0.)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(0, 0).angle(
            pygsf.orientations.orientations.Direct.fromAzPl(135, 0)), 135.)
        assert areClose(pygsf.orientations.orientations.Direct.fromAzPl(0, 80).angle(
            pygsf.orientations.orientations.Direct.fromAzPl(180, 80)), 20.)

    def test_axis_angle(self):

        assert areClose(pygsf.orientations.orientations.Axis.fromAzPl(90, 0).angle(
            pygsf.orientations.orientations.Axis.fromAzPl(270, 0)), 0.)

    def test_plane_normal(self):

        assert areClose(
            pygsf.orientations.orientations.Plane(90, 45).normDirectFrwrd().angle(
                pygsf.orientations.orientations.Direct.fromAzPl(90, -45)), 0.)

    def test_plane2cplane(self):

        pl = pygsf.orientations.orientations.Plane(90, 45).toCPlane(Point(0, 0, 0))
        assert areClose(pl.angle(CPlane(1, 0, 1, 0), 0.))

    def test_plane_angle(self):

        assert areClose(
            pygsf.orientations.orientations.Plane(90, 45).angle(
                pygsf.orientations.orientations.Plane(90, 45)), 0.)
        assert areClose(
            pygsf.orientations.orientations.Plane(90, 45).angle(
                pygsf.orientations.orientations.Plane(90, 55)), 10.)
        assert areClose(
            pygsf.orientations.orientations.Plane(90, 5).angle(
                pygsf.orientations.orientations.Plane(270, 5)), 10.)
        assert areClose(
            pygsf.orientations.orientations.Plane(90, 85).angle(
                pygsf.orientations.orientations.Plane(270, 85)), 10.)
        assert areClose(
            pygsf.orientations.orientations.Plane(0, 0).angle(
                pygsf.orientations.orientations.Plane(0, 10)), 10.)
        assert areClose(
            pygsf.orientations.orientations.Plane(0, 0).angle(
                pygsf.orientations.orientations.Plane(180, 0)), 0.)

    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()
