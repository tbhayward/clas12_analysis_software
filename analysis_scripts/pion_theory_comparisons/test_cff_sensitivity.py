import unittest
import numpy as np
from analyze_cff_sensitivity import integral_abs
class IntegrationTests(unittest.TestCase):
 def test_sign_crossings(self):
  self.assertAlmostEqual(integral_abs(np.array([0.,1.]),np.array([-1.,1.])),.5)
  self.assertAlmostEqual(integral_abs(np.array([0.,2.]),np.array([-2.,1.])),5/3)
 def test_constant_and_zero(self):
  self.assertAlmostEqual(integral_abs(np.array([.05,.65,1.2]),np.ones(3)*.2),.23)
  self.assertEqual(integral_abs(np.array([0.,1.]),np.zeros(2)),0.)
if __name__=='__main__':unittest.main()
