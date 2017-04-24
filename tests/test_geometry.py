

import unittest

from gsf.math_utils import *
from gsf.geometry import *


class TestGsf(unittest.TestCase):

    def setUp(self):

        pass

    def test_gvect_general(self):
        """
        Check expected GVect results for downward dip.
        """

        assert GVect(90, 90).is_downward
        assert GVect(90, -45).is_upward
        assert isclose(GVect(90, 90).versor.z, -1.0)
        assert isclose(GVect(90, -90).versor.z, 1.0)
        assert isclose(GVect(0, 90).upward.versor.z, 1.0)
        assert isclose(GVect(0, -90).downward.versor.z, -1.0)

    def test_gvect_angle(self):

        assert isclose(GVect(90, 45).angle(GVect(90, 55)), 10.)
        assert isclose(GVect(90, 45).angle(GVect(270, 10)), 125.)
        assert isclose(GVect(90, 90).angle(GVect(135, 90)), 0.)
        assert isclose(GVect(0, 0).angle(GVect(135, 0)), 135.)
        assert isclose(GVect(0, 80).angle(GVect(180, 80)), 20.)

    def test_gaxis_angle(self):

        assert isclose(GAxis(90, 0).angle(GAxis(270, 0)), 0.)

    def test_gplane_normal(self):

        assert GPlane(90, 45).normal.versor.is_downward

    def test_gplane_plane(self):

        pl = GPlane(90, 45).plane(Point(0, 0, 0))
        assert isclose(pl.angle(Plane(1, 0, 1, 0)), 0.)

    def test_gplane_angle(self):

        assert isclose(GPlane(90, 45).angle(GPlane(90, 45)), 0.)
        assert isclose(GPlane(90, 45).angle(GPlane(90, 55)), 10.)
        assert isclose(GPlane(90, 5).angle(GPlane(270, 5)), 10.)
        assert isclose(GPlane(90, 85).angle(GPlane(270, 85)), 10.)
        assert isclose(GPlane(0, 0).angle(GPlane(0, 10)), 10.)
        assert isclose(GPlane(0, 0).angle(GPlane(180, 0)), 0.)


    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()