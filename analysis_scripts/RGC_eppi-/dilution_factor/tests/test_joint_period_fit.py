import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from joint_period_fit import fit_bin, replica_bin
from model import COEFFICIENTS as c
class JointTest(unittest.TestCase):
 def test_known_common_fraction(self):
  inputs=[];f=.35
  for p in ['summer2022','fall2022','spring2023']:
   r=dict(C=.002,CH2=.003,He=.0002,ET=.00005,CD2=.004)
   s=r['CD2']+c['A']*r['C']+c['B']*r['He']+c['C']*r['ET']
   d=r['CD2']+c['D']*r['C']+c['E']*r['He']+c['F']*r['ET']
   r['ND3']=r['He']*s/(s-f*d)
   for target,rate in r.items():
    if target=='CD2' and p!='spring2023':continue
    inputs.append(dict(dataset='RGC_'+p+'_'+target,charge=1e8,selected_counts=[rate*1e8]*24))
  result=fit_bin(inputs,0)
  self.assertAlmostEqual(result.x[16],f,places=5)
  self.assertLess(abs(result.fun),1e-5)
  at=fit_bin(inputs,0,fixed=f,start=result.x)
  shifted=fit_bin(inputs,0,fixed=.45,start=result.x)
  self.assertGreater(shifted.fun,at.fun+1)
  row,values=replica_bin(({'inputs':inputs},0,16,42))
  row2,values2=replica_bin(({'inputs':inputs},0,16,42))
  np.testing.assert_array_equal(values,values2)
  self.assertAlmostEqual(row['stat_std'],np.std(values,ddof=1))
  self.assertEqual(row['replicas_successful']+row['replicas_failed'],16)
if __name__=='__main__':unittest.main()
