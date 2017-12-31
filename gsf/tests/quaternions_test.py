
import unittest

from ..quaternions import Quaternion 

q_null = Quaternion(0, 0, 0, 0)
q_one = Quaternion(1, 0, 0, 0)
q_x = Quaternion(0, 1, 0, 0)
q_y = Quaternion(0, 0, 1, 0)
q_z = Quaternion(0, 0, 0, 1)

class TestQuaternions(unittest.TestCase):

    def test_sqrd_norm(self):

        self.assertAlmostEqual(q_null.sqrd_norm(), 0.0)
        self.assertAlmostEqual(q_one.sqrd_norm(), 1.0)
        self.assertAlmostEqual(q_y.sqrd_norm(), 1.0)