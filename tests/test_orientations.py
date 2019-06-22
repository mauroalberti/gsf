# -*- coding: utf-8 -*-


import unittest

from pygsf.geology import orientations as orient
from pygsf.mathematics.scalars import areClose
from pygsf.spatial.vectorial.geometries import Point, CPlane


class TestOrientations(unittest.TestCase):

    def setUp(self):

        pass

    def test_direct_general(self):
        """
        Check expected OrienM results for downward dip.
        """

        assert orient.Direct.fromAzPl(90, 90).isDownward
        assert orient.Direct.fromAzPl(90, -45).isUpward
        assert areClose(orient.Direct.fromAzPl(90, 90).asVersor().z, -1.0)
        assert areClose(orient.Direct.fromAzPl(90, -90).asVersor().z, 1.0)
        assert areClose(orient.Direct.fromAzPl(0, 90).upward().asVersor().z, 1.0)
        assert areClose(orient.Direct.fromAzPl(0, -90).downward().asVersor().z, -1.0)

    def test_direct_angle(self):

        assert areClose(orient.Direct.fromAzPl(90, 45).angle(orient.Direct.fromAzPl(90, 55)), 10.)
        assert areClose(orient.Direct.fromAzPl(90, 45).angle(orient.Direct.fromAzPl(270, 10)), 125.)
        assert areClose(orient.Direct.fromAzPl(90, 90).angle(orient.Direct.fromAzPl(135, 90)), 0.)
        assert areClose(orient.Direct.fromAzPl(0, 0).angle(orient.Direct.fromAzPl(135, 0)), 135.)
        assert areClose(orient.Direct.fromAzPl(0, 80).angle(orient.Direct.fromAzPl(180, 80)), 20.)

    def test_axis_angle(self):

        assert areClose(orient.Axis.fromAzPl(90, 0).angle(orient.Axis.fromAzPl(270, 0)), 0.)

    def test_plane_normal(self):

        assert areClose(orient.Plane(90, 45).normDirectFrwrd().angle(orient.Direct.fromAzPl(90, -45)), 0.)

    def test_plane2cplane(self):

        pl = orient.Plane(90, 45).toCPlane(Point(0, 0, 0, epsg_cd=2000))
        assert areClose(pl.angle(CPlane(1, 0, 1, 0, epsg_cd=2000)), 0.)

    def test_plane_angle(self):

        assert areClose(orient.Plane(90, 45).angle(orient.Plane(90, 45)), 0.)
        assert areClose(orient.Plane(90, 45).angle(orient.Plane(90, 55)), 10.)
        assert areClose(orient.Plane(90, 5).angle(orient.Plane(270, 5)), 10.)
        assert areClose(orient.Plane(90, 85).angle(orient.Plane(270, 85)), 10.)
        assert areClose(orient.Plane(0, 0).angle(orient.Plane(0, 10)), 10.)
        assert areClose(orient.Plane(0, 0).angle(orient.Plane(180, 0)), 0.)

    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()
