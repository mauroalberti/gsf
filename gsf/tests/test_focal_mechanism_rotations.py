
from pprint import pprint

import unittest

from ..geometry import GAxis
from ..ptbaxes import PTBAxes
from ..rotations import RotationAxis



class TestFocalMechamismRotations(unittest.TestCase):

    def test_kagan_examples(self):
        """
        Test of focal mechanims rotations examples
        as described in Kaga Y.Y. 3-D rotation of double-couple earthquake sources

        :return:
        """

        first_solution_T_axis_vals = (120, 24)
        first_solution_P_axis_vals = (232, 41)

        second_solution_T_axis_vals = (295, 55)
        second_solution_P_axis_vals = (51, 17)

        fs_taxis = GAxis(*first_solution_T_axis_vals)
        fs_paxis = GAxis(*first_solution_P_axis_vals)

        first_solution = PTBAxes(
            p_axis=fs_paxis,
            t_axis=fs_taxis)

        print(first_solution)

        ss_taxis = GAxis(*second_solution_T_axis_vals)
        ss_paxis = GAxis(*second_solution_P_axis_vals)

        second_solution = PTBAxes(
            p_axis=ss_paxis,
            t_axis=ss_taxis)

        print(second_solution)

        fs_quaternion = first_solution.to_quaternion()

        print(fs_quaternion)
        print(fs_quaternion.to_rotation_axis())
        print(fs_quaternion.inverse)
        print(fs_quaternion.inverse.to_rotation_axis())
        print(fs_quaternion * (-1))
        print((fs_quaternion * (-1)).to_rotation_axis())

        """
        rotations = first_solution.calculate_rotations(second_solution)

        for rotation in rotations:
            print(rotation)
        """


if __name__ == '__main__':
    unittest.main()

