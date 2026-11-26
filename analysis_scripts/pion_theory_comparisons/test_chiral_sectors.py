import unittest
import numpy as np
from chiral_sectors import classify,grouped
class SectorTests(unittest.TestCase):
 def test_physical_token_categories_and_mixed_interference(self):
  for t in ['F11xG14','HTxcalE','EtL_polexcalE','HTxcalEt','F11xHtL']:self.assertEqual(classify(t),'added')
  for t in ['HTxEtL_pole','EbarTxHtL','HTxEtL_np']:self.assertEqual(classify(t),'chiral_odd')
  for t in ['HtLxEtL_pole','EtL_npxEtL_pole']:self.assertEqual(classify(t),'longitudinal_pole')
  self.assertEqual(classify('other'),'unknown')
  with self.assertRaises(ValueError):classify('HTxNotDefined')
 def test_grouping_is_additive_and_keeps_other_separate(self):
  b={'terms':['F11xG14','HTxEtL_pole','HTxcalE','other'],'arrays':{'total':np.array([.5,-.3]),'F11xG14':np.array([.2,-.4]),'HTxEtL_pole':np.array([.4,.2]),'HTxcalE':np.array([-.2,-.2]),'other':np.array([.1,.1])}}
  g=grouped(b);np.testing.assert_allclose(sum(g.values()),b['arrays']['total']);np.testing.assert_allclose(g['added'],[0,-.6]);np.testing.assert_allclose(g['unknown'],[.1,.1])
if __name__=='__main__':unittest.main()
