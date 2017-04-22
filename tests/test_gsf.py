
import unittest

from gsf.config import dip_conventions
from gsf import *


class TestGsf(unittest.TestCase):

    def setUp(self):

        pass

    def test_dip_convention(self):

        assert len(dip_conventions) == 2
        assert "positive_up" in dip_conventions
        assert "positive_down" in dip_conventions

    """
    def test_axes_angle(self):

        self.assertTrue(GAxis(90, 0).angle.(GAxis(270, 0)) == 0.0)
    """


    def tearDown(self):

        pass


if __name__ == '__main__':
    unittest.main()