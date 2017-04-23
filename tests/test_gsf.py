

import unittest

from gsf.math_utils import *
from gsf.geosurfaces import *


class TestGsf(unittest.TestCase):

    def setUp(self):

        pass

    def test_dip_convention(self):
        """
        Check that predefined parameters for dip convention
        are as expected.
        """

        assert len(dip_conventions) == 2
        assert "positive_up" in dip_conventions
        assert "positive_down" in dip_conventions
        assert dip_angle_convention == default_dip_angle_convention

    def test_gvect_up(self):
        """
        Check expected GVect results for upward dip convention.
        """

        set_dip_angle_positive_up()
        assert GVect(90, 90).is_upward
        assert GVect(90, -45).is_downward
        assert isclose(GVect(90, 90).versor.z, 1.0)
        assert isclose(GVect(90, -90).versor.z, -1.0)
        assert isclose(GVect(0, -90).upward.versor.z, 1.0)
        assert isclose(GVect(0, 90).downward.versor.z, -1.0)

    def test_gvect_down(self):
        """
        Check expected GVect results for downward dip convention.
        """

        set_dip_angle_positive_down()
        assert GVect(90, 90).is_downward
        assert GVect(90, -45).is_upward
        assert isclose(GVect(90, 90).versor.z, -1.0)
        assert isclose(GVect(90, -90).versor.z, 1.0)
        assert isclose(GVect(0, 90).upward.versor.z, 1.0)
        assert isclose(GVect(0, -90).downward.versor.z, -1.0)

    def test_gvect_angle(self):

        set_dip_angle_positive_up()
        assert isclose(GVect(90, 45).angle(GVect(90, 55)), 10.)
        assert isclose(GVect(90, 45).angle(GVect(270, 10)), 125.)
        assert isclose(GVect(90, 90).angle(GVect(135, 90)), 0.)
        assert isclose(GVect(0, 0).angle(GVect(135, 0)), 135.)
        assert isclose(GVect(0, 80).angle(GVect(180, 80)), 20.)
        set_dip_angle_positive_down()
        assert isclose(GVect(90, 45).angle(GVect(90, 55)), 10.)
        assert isclose(GVect(90, 45).angle(GVect(270, 10)), 125.)
        assert isclose(GVect(90, 90).angle(GVect(135, 90)), 0.)
        assert isclose(GVect(0, 0).angle(GVect(135, 0)), 135.)
        assert isclose(GVect(0, 80).angle(GVect(180, 80)), 20.)

    def test_gplane_normal(self):

        set_dip_angle_positive_up()
        assert GPlane(90, 45).normal.versor.is_upward
        set_dip_angle_positive_down()
        assert GPlane(90, 45).normal.versor.is_downward

    def test_gplane_plane(self):

        set_dip_angle_positive_up()
        pl = GPlane(90, 45).plane(Point(0, 0, 0))
        assert isclose(pl.angle(Plane(-1, 0, 1, 0)), 0.)
        set_dip_angle_positive_down()
        pl = GPlane(90, 45).plane(Point(0, 0, 0))
        assert isclose(pl.angle(Plane(1, 0, 1, 0)), 0.)

    def test_gplane_angle(self):

        set_dip_angle_positive_up()
        assert isclose(GPlane(90, 45).angle(GPlane(90, 45)), 0.)
        assert isclose(GPlane(90, 45).angle(GPlane(90, 55)), 10.)
        assert isclose(GPlane(90, 5).angle(GPlane(270, 5)), 10.)
        assert isclose(GPlane(90, 85).angle(GPlane(270, 85)), 10.)


    """
    def test_axes_angle(self):

        self.assertTrue(GAxis(90, 0).angle.(GAxis(270, 0)) == 0.0)
    """


    def tearDown(self):

        pass


if __name__ == '__main__':

    unittest.main()