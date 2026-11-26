"""The conference snapshot must preserve fine LU and the three-bin free UL baseline."""
import unittest
from pathlib import Path
from theory_io import load_measurements
HERE=Path(__file__).parent
MEASUREMENTS=HERE/'measurements/measurements.csv' if (HERE/'measurements').exists() else HERE.parent/'data/measurements/measurements.csv'


class ReleaseInputTests(unittest.TestCase):
 def test_conference_binning(self):
  data=load_measurements(MEASUREMENTS)
  pim=[r for r in data if r['channel']=='pim']
  # Legacy data may be loaded explicitly for historical comparisons.
  if len(pim)!=48:self.skipTest('Historical measurement snapshot')
  self.assertEqual({r['observable'] for r in pim},{'lu1','ul1','ul2'})
  for obs,n,edges in [('lu1',6,[.05,.25,.45,.65,.85,1.05,1.25]),('ul1',3,[.05,.45,.85,1.25]),('ul2',3,[.05,.45,.85,1.25])]:
   for xi in range(4):
    rows=sorted([r for r in pim if r['observable']==obs and r['x_index']==xi],key=lambda r:r['bin'])
    self.assertEqual(len(rows),n)
    self.assertEqual([r['x_low'] for r in rows],edges[:-1])
    self.assertEqual([r['x_high'] for r in rows],edges[1:])
    self.assertTrue(all(r['status']=='measured' and r['stat']>0 for r in rows))
  self.assertEqual(len([r for r in data if r['channel']=='pip']),168)


if __name__=='__main__':unittest.main()
