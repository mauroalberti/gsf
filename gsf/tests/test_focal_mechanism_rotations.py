
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

        first_solution_T_axis_vals = (120, 24)
        first_solution_P_axis_vals = (232, 41)

        second_solution_T_axis_vals = (295, 55)
        second_solution_P_axis_vals = (51, 17)

        fs_taxis = GAxis(*first_solution_T_axis_vals)
        fs_paxis = GAxis(*first_solution_P_axis_vals)

        fst_sol_ptbaxes = PTBAxes(
            p_axis=fs_paxis,
            t_axis=fs_taxis)

        ss_taxis = GAxis(*second_solution_T_axis_vals)
        ss_paxis = GAxis(*second_solution_P_axis_vals)

        snd_sol_ptbaxes = PTBAxes(
            p_axis=ss_paxis,
            t_axis=ss_taxis)

        fst_sol_quater = fst_sol_ptbaxes.to_quaternion()
        snd_sol_quater = snd_sol_ptbaxes.to_quaternion()

        qr_o = snd_sol_quater * (fst_sol_quater.inverse)

        a_i = snd_sol_quater * (Quaternion.i() * snd_sol_quater.inverse)
        a_j = snd_sol_quater * (Quaternion.j() * snd_sol_quater.inverse)
        a_k = snd_sol_quater * (Quaternion.k() * snd_sol_quater.inverse)

        qr_i = a_i * qr_o
        qr_j = a_j * qr_o
        qr_k = a_k * qr_o


if __name__ == '__main__':
    unittest.main()

