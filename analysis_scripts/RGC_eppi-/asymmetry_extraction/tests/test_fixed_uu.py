import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from likelihood import fit
from extract_structure_function_ratios import read_uu_inputs
# Installed tests provide a reproducible physical toy.
sys.path.insert(0,str(Path(__file__).resolve().parent))
from test_likelihood import toy

class FixedUU(unittest.TestCase):
 def test_nonzero_inputs_are_held_fixed(self):
  g,truth=toy(n=20000)
  result=fit([g],(.35,.02),fixed_uu=(.2,-.12))
  self.assertEqual(result['status'],'valid')
  np.testing.assert_array_equal(result['values'][:2],[.2,-.12])
  np.testing.assert_array_equal(result['errors'][:2],[0,0])
  self.assertEqual(len(result['free_parameters']),5)
  self.assertTrue(np.all(np.isfinite(result['correlation'])))
 def test_rgb_only_lu_is_free(self):
  g,truth=toy(n=10000,polarized=False)
  result=fit([g],(.35,.02),fixed_uu=(.2,-.12))
  self.assertEqual(result['status'],'valid')
  self.assertEqual(result['free_parameters'],['LU_sin_phi'])
 def test_import_requires_complete_unique_bins(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'uu.csv'
   p.write_text('bin,UU_cos_phi,UU_cos_2phi\n'+''.join(f'{b},.2,-.12\n' for b in range(1,25)))
   self.assertEqual(read_uu_inputs(p)[24],(.2,-.12))
   p.write_text(p.read_text()+'1,0,0\n')
   with self.assertRaises(ValueError):read_uu_inputs(p)
   p.write_text('bin,UU_cos_phi,UU_cos_2phi\n1,0,0\n')
   with self.assertRaises(ValueError):read_uu_inputs(p)

if __name__=='__main__':unittest.main()
