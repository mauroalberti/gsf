# -*- coding: utf-8 -*-


import unittest

from ..geometry import GAxis
from ..ptbaxes import PTBAxes
from ..quaternions import *
from ..rotations import RotationAxis


class TestFocalMechamismRotations(unittest.TestCase):

    def test_kagan_examples(self):
        """
        Test of focal mechanims rotations examples
        as described in Kagan Y.Y. 3-D rotation of double-couple earthquake sources

        :return:
        """

        fm1 = PTBAxes(p_axis=GAxis(232, 41), t_axis=GAxis(120, 24))
        fm2 = PTBAxes(p_axis=GAxis(51, 17), t_axis=GAxis(295, 55))

        rots = fm1.calculate_rotations(fm2)

        for rot in rots:
            print(rot)

