

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
        assert isclose(GVect(90, 90).versor.z, 1.0)
        assert isclose(GVect(90, -90).versor.z, -1.0)

    def test_gvect_down(self):
        """
        Check expected GVect results for downward dip convention.
        """

        set_dip_angle_positive_down()
        assert isclose(GVect(90, 90).versor.z, -1.0)
        assert isclose(GVect(90, -90).versor.z, 1.0)

    """
    def test_axes_angle(self):

        self.assertTrue(GAxis(90, 0).angle.(GAxis(270, 0)) == 0.0)
    """


    def tearDown(self):

        pass


if __name__ == '__main__':
    unittest.main()