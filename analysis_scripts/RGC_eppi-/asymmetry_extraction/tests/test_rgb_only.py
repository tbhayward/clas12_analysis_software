import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from rgb_only_check import RGBObjective

class RGBSupportTest(unittest.TestCase):
 def test_legacy_check_can_allow_negative_other_state(self):
  # Ten run/helicity states: an observed-run probability below one does not
  # guarantee every other run/helicity probability is nonnegative.
  states=np.array([[.1,.8 if i%2 else -.8,0.,0.] for i in range(10)])
  events=dict(phi=np.array([np.pi/2]),rb=np.array([.5]),rc=np.array([.8]),rv=np.array([1.]),rw=np.array([1.]),h=np.array([1.]),pb=np.array([.8]),pt=np.array([0.]),q_fraction=np.array([.1]))
  groups=[dict(events=events,states=states)]
  x=np.array([0.,1.5,.5])
  strict=RGBObjective(groups,legacy=False);legacy=RGBObjective(groups,legacy=True)
  self.assertLess(np.min(1+strict.core.allowed@x),0)
  self.assertEqual(strict(*x),1e30)
  self.assertLess(legacy(*x),1e29)

if __name__=='__main__':unittest.main()
