"""Checks for experimental-bin routing and supplied CFF normalization."""
from pathlib import Path
import unittest
import numpy as np
from theory_io import load_theory,load_measurements,OBS
from plot_comparisons import components,limits,state_data,focused_limits
from cff_style import configure_styles,BLUE,ORANGE,GREEN,RED,PURPLE
HERE=Path(__file__).parent
INPUTS=HERE/'inputs' if (HERE/'inputs').exists() else HERE.parent/'data/theory'
MEASUREMENTS=HERE/'measurements/measurements.csv' if (HERE/'measurements/measurements.csv').exists() else HERE.parent/'data/measurements/measurements.csv'
class ComparisonTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.theory,cls.audit=load_theory(INPUTS);cls.data=load_measurements(MEASUREMENTS)
 def test_all_three_species_have_seven_ratios(self):
  self.assertEqual(len(self.theory),21)
  self.assertTrue(all(r['max_abs_closure']<.00004 for r in self.audit if 'max_abs_closure' in r))
 def test_clarity_palette_is_unambiguous_per_observable(self):
  registry=configure_styles(self.theory)
  for obs,styles in registry.items():
   products=[v for t,v in styles.items() if t not in ['other','remainder']]
   self.assertEqual(len(set(products)),len(products))
  self.assertEqual(registry['ul1']['HTxEtL_pole'],(RED,'solid'))
  self.assertEqual(registry['ul1']['EtL_polexcalE'],(ORANGE,'solid'))
  self.assertEqual(registry['ul1']['EbarTxHtL'],(GREEN,'solid'))
  self.assertEqual(registry['ul1']['F11xHtL'],(BLUE,'solid'))
  self.assertEqual(registry['ul1']['HTxEtL_np'],(PURPLE,'solid'))
  self.assertEqual(registry['ul2']['F11xG14'],(BLUE,'solid'))
 def test_full_and_leading_decompositions_close(self):
  for blocks in self.theory.values():
   leading=sorted([t for t in blocks[0]['terms'] if t!='other'],key=lambda t:-max(np.max(np.abs(b['arrays'][t])) for b in blocks))
   for b in blocks:
    b['leading_terms']=leading
    for n in [None,5]:np.testing.assert_allclose(sum(v for _,v in components(b,n)),b['arrays']['total'],rtol=0,atol=.00004)
 def test_pi_minus_uses_six_lu_and_two_target_bins(self):
  for xi in range(4):
   self.assertEqual(len(state_data(self.data,'pim','lu1',xi)),6)
   for o in ['ul1','ul2']:self.assertIn(len(state_data(self.data,'pim',o,xi)),[2,3])
   for o in ['ll0','ll1']:self.assertIn(len(state_data(self.data,'pim',o,xi)),[0,2])
   for o in ['u1','u2']:self.assertEqual(len(state_data(self.data,'pim',o,xi)),0)
 def test_fixed_axes_respect_requested_ranges(self):
  for species in ['pip','pim']:
   for obs in OBS:
    for decomp in [False,True]:
     low,high=limits(self.theory,self.data,species,obs,decomp)
     self.assertTrue(-1<=low<high<=1)
  expected={'u1':(-.4,.8),'u2':(-.8,.6),'lu1':(-.2,.3),'ul1':(-.3,.3),'ul2':(-.4,.1),'ll0':(-.2,1.),'ll1':(-.4,.8)}
  for species in ['pip','pim']:
   for obs,bounds in expected.items():
    actual=limits(self.theory,self.data,species,obs,True)
    self.assertEqual(actual,(-.2,.3) if species=='pim' and obs=='lu1' else (-1.,1.) if species=='pim' and obs=='ul2' else bounds)
    self.assertEqual(focused_limits(species,obs,True,actual),actual)
if __name__=='__main__':unittest.main()
