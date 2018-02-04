# -*- coding: utf-8 -*-


import unittest

from .test_datasets import *

from ..rotations import *


class TestFocalMechamismRotations(unittest.TestCase):

    def test_kagan_examples_1(self):
        """
        Test of focal mechanims rotations examples
        as described in Kagan Y.Y. 3-D rotation of double-couple earthquake sources

        :return:
        """

        def sols2rotaxis(rot_sol: dict) -> RotationAxis:
            az = rot_sol["az"]
            colat = rot_sol["colat_b"]
            rot_ang = rot_sol["rot_ang"]

            pl_from_N = plng2colatBottom(colat)

            gv = GVect(az, pl_from_N)
            return RotationAxis.from_gvect(gv, rot_ang)

        rots = focmechs_disorientations(k91_fs_PTBaxes, k91_ss_PTBaxes)

        for ndx, rot in enumerate(rots):
            print("calculated solution: {}".format(rot))
            k91_sol = sols2rotaxis(k91_rot_sols[ndx])
            print("Kagan 1991 solution: {}".format(k91_sol))
            assert rot.strictly_equival(k91_sol)

    def test_kagan_examples_2(self):
        """
        Test of focal mechanims rotations calculations
        as described in Kagan Y.Y. 3-D rotation of double-couple earthquake sources

        :return:
        """

        qr_o = k91_ss_quater * k91_fs_quater.inverse

        a_i = k91_ss_quater * (Quaternion.i() * k91_ss_quater.inverse)
        a_j = k91_ss_quater * (Quaternion.j() * k91_ss_quater.inverse)
        a_k = k91_ss_quater * (Quaternion.k() * k91_ss_quater.inverse)

        qr_i = a_i * qr_o
        qr_j = a_j * qr_o
        qr_k = a_k * qr_o

        """
        print(qr_o, RotationAxis.from_quaternion(qr_o).to_min_rotation_axis())
        print(qr_i, RotationAxis.from_quaternion(qr_i).to_min_rotation_axis())
        print(qr_j, RotationAxis.from_quaternion(qr_j).to_min_rotation_axis())
        print(qr_k, RotationAxis.from_quaternion(qr_k).to_min_rotation_axis())
        """

if __name__ == '__main__':
    unittest.main()

