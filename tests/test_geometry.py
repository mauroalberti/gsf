

import unittest

from ..geography import *


class TestGeometry(unittest.TestCase):

    def setUp(self):

        pass

    def test_gvect_general(self):
        """
        Check expected OrienM results for downward dip.
        """

        assert GVect(90, 90).isDownward
        assert GVect(90, -45).isUpward
        assert are_close(GVect(90, 90).versor().z, -1.0)
        assert are_close(GVect(90, -90).versor().z, 1.0)
        assert are_close(GVect(0, 90).upward().versor().z, 1.0)
        assert are_close(GVect(0, -90).downward().versor().z, -1.0)

    def test_gvect_angle(self):

        assert are_close(GVect(90, 45).angle(GVect(90, 55)), 10.)
        assert are_close(GVect(90, 45).angle(GVect(270, 10)), 125.)
        assert are_close(GVect(90, 90).angle(GVect(135, 90)), 0.)
        assert are_close(GVect(0, 0).angle(GVect(135, 0)), 135.)
        assert are_close(GVect(0, 80).angle(GVect(180, 80)), 20.)

    def test_gaxis_angle(self):

        assert are_close(GAxis(90, 0).angle(GAxis(270, 0)), 0.)

    def test_gplane_normal(self):

        assert are_close(GPlane(90, 45).normalOrienFrwrd().angle(GVect(90, -45)), 0.)

    def test_gplane_plane(self):

        pl = GPlane(90, 45).plane(Point(0, 0, 0))
        assert are_close(pl.angle(CPlane(1, 0, 1, 0)), 0.)

    def test_gplane_angle(self):

        assert are_close(GPlane(90, 45).angle(GPlane(90, 45)), 0.)
        assert are_close(GPlane(90, 45).angle(GPlane(90, 55)), 10.)
        assert are_close(GPlane(90, 5).angle(GPlane(270, 5)), 10.)
        assert are_close(GPlane(90, 85).angle(GPlane(270, 85)), 10.)
        assert are_close(GPlane(0, 0).angle(GPlane(0, 10)), 10.)
        assert are_close(GPlane(0, 0).angle(GPlane(180, 0)), 0.)


    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()