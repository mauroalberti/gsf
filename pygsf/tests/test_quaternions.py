
import unittest

from pygsf.mathematics.scalars import are_close
from pygsf.mathematics.quaternions import Quaternion


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

        quat_1 = Quaternion(0.696, 0.322, -0.152, 0.624)
        assert are_close(quat_1.normalize().norm, 1.0)


if __name__ == '__main__':

    unittest.main()


