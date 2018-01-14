
import unittest

from ..quaternions import Quaternion 


q_case_1 = Quaternion(3.2, 17.4, 9.25, -8.47)


class TestQuaternions(unittest.TestCase):

    def test_sqrd_norm(self):

        self.assertAlmostEqual(Quaternion.zero().sqrd_norm(), 0.0)
        self.assertAlmostEqual(Quaternion.identity().sqrd_norm(), 1.0)
        self.assertAlmostEqual(Quaternion.i().sqrd_norm(), 1.0)
        self.assertAlmostEqual(Quaternion.j().sqrd_norm(), 1.0)
        self.assertAlmostEqual(Quaternion.k().sqrd_norm(), 1.0)

    def test_normalized(self):

        norm_quat = q_case_1.normalize()

        self.assertAlmostEqual(norm_quat.sqrd_norm(), 1.0)

        cnj_norm = norm_quat.conjugate
        inv_norm = norm_quat.inverse
        assert cnj_norm.is_close_to(inv_norm)

    def test_rotation_axis(self):

        quat = Quaternion.identity()
        print(quat)
        print(quat.to_rotation_axis())


if __name__ == '__main__':

    unittest.main()


