"""Regression checks for updated CFF totals in the existing pi+ entry point."""
from pathlib import Path
import shutil,tempfile,unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_kroll_comparisons as plot
class CFFDirectoryTests(unittest.TestCase):
 def setUp(self):
  self.path=plot.DEFAULT_CFF_DIRECTORY
  if not self.path.is_dir():
   candidates=[Path(__file__).parents[1]/'inputs/pip',Path(__file__).parents[2]/'data/theory/pip']
   self.path=next(p for p in candidates if p.is_dir())
 def test_new_totals_are_copied_exactly(self):
  result=plot.load_new_theory(self.path)
  self.assertFalse(result.attrs['uncertainty_supplied'])
  for obs,stem in plot.NEW_THEORY_COLUMNS.items():
   raw=pd.read_csv(self.path/(stem+'.csv')).sort_values(['xB','minus_tprime'])
   np.testing.assert_array_equal(result[obs],raw.total)
 def test_no_band_is_invented(self):
  result=plot.load_new_theory(self.path);fig,ax=plt.subplots()
  try:
   plot.draw_new_theory_curve(ax,result[result.xB==.2],'ul2','black')
   self.assertEqual(len(ax.collections),0)
  finally:plt.close(fig)
 def test_nonclosing_decomposition_is_rejected(self):
  with tempfile.TemporaryDirectory() as tmp:
   dest=Path(tmp)/'pip';shutil.copytree(self.path,dest)
   p=dest/'A_UL_sinphi.csv';a=pd.read_csv(p);a.loc[0,'total']+=.01;a.to_csv(p,index=False)
   with self.assertRaisesRegex(ValueError,'sum'):plot.load_new_theory(dest)
if __name__=='__main__':unittest.main()
