import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_production_results import axes_style,PANELS

class ProductionStyle(unittest.TestCase):
 def test_pi_plus_layout_and_limits(self):
  self.assertEqual(PANELS,[0,1,3,4,5,6,7])
  for j in range(7):
   fig,ax=plt.subplots();axes_style(ax,j)
   self.assertEqual(ax.get_ylim(),(-.4,.4) if j in [2,3,4] else (-1.,1.))
   self.assertEqual(ax.get_xlim(),(0.,1.3))
   self.assertEqual(len(ax.patches),0)
   plt.close(fig)

if __name__=='__main__':unittest.main()
