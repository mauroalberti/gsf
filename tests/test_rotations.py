# -*- coding: utf-8 -*-


import unittest

from .test_datasets import *

from pygsf.orientations.rotations import *


class TestFocalMechamismRotations(unittest.TestCase):

    def test_rotate_focal_mechanism(self):
        """
        Test whether a focal mechanism is correctly rotated.

        :return:
        """

        src_fm = k91_fs_PTBaxes
        rot_axes = map(sols2rotaxis, k91_rot_sols)
        rot_fm = k91_ss_PTBaxes

        print("\nTest forward focal mechanism rotation")
        for rot_axis in rot_axes:
            calc_rot_fm = focmech_rotate(src_fm, rot_axis)
            print("calculated rotated fm: {}".format(calc_rot_fm))
            print("expected rotated fm: {}".format(rot_fm))
            assert calc_rot_fm.almostEqual(rot_fm)

    def test_inversion_kagan_examples_1(self):
        """
        Test of focal mechanims rotations examples
        as described in Kagan Y.Y. 3-D rotation of double-couple earthquake sources

        :return:
        """

        rots = focmechs_invert_rotations(k91_fs_PTBaxes, k91_ss_PTBaxes)

        print("\nTest inverse focal mechanism rotation")
        for ndx, rot in enumerate(rots):
            print("calculated solution: {}".format(rot))
            k91_sol = sols2rotaxis(k91_rot_sols[ndx])
            print("Kagan 1991 solution: {}".format(k91_sol))
            assert rot.strictly_equival(k91_sol)


if __name__ == '__main__':

    unittest.main()

