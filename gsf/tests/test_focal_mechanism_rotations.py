
from pprint import pprint

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

        fs_T_axis_vals, fs_P_axis_vals = (120, 24), (232, 41)
        ss_T_axis_vals, ss_P_axis_vals = (295, 55), ( 51, 17)

        fs_T_gaxis, fs_P_gaxis = GAxis(*fs_T_axis_vals), GAxis(*fs_P_axis_vals)
        ss_T_gaxis, ss_P_gaxis = GAxis(*ss_T_axis_vals), GAxis(*ss_P_axis_vals)

        fs_PTBaxes = PTBAxes(
            p_axis=fs_P_gaxis,
            t_axis=fs_T_gaxis)

        ss_PTBaxes = PTBAxes(
            p_axis=ss_P_gaxis,
            t_axis=ss_T_gaxis)

        fs_quater = fs_PTBaxes.to_quaternion()
        ss_quater = ss_PTBaxes.to_quaternion()

        qr_o = ss_quater * fs_quater.inverse

        a_i = ss_quater * (Quaternion.i() * ss_quater.inverse)
        a_j = ss_quater * (Quaternion.j() * ss_quater.inverse)
        a_k = ss_quater * (Quaternion.k() * ss_quater.inverse)

        qr_i = a_i * qr_o
        qr_j = a_j * qr_o
        qr_k = a_k * qr_o

        print(qr_o, qr_o.to_rotation_axis())
        print(qr_i, qr_i.to_rotation_axis())
        print(qr_j, qr_j.to_rotation_axis())
        print(qr_k, qr_k.to_rotation_axis())


if __name__ == '__main__':
    unittest.main()

