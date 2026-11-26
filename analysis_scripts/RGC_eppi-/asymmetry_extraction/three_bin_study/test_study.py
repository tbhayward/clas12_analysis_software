"""Check rebin conservation, fixed-ratio covariance, and spawn-worker execution."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key]='1'
from pathlib import Path
import sys
import unittest
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[0]))
from constrained_fit import fit
from uu_model_study.likelihood import fit as original
from run_three_bin_study import MEMBERS,initialize,worker


def synthetic():
    rng=np.random.default_rng(191)
    groups=[]
    for period in range(2):
        n=2000
        h=rng.choice([-1,1],n)
        ev=dict(phi=rng.uniform(0,2*np.pi,n),rb=np.full(n,.7),rc=np.full(n,.6),rv=np.full(n,.9),rw=np.full(n,.8),h=h,pb=np.full(n,.84),pt=rng.choice([-.22,.22],n))
        states=np.array([[1,h*.84,t,h*.84*t] for h in [-1,1] for t in [-.22,.22]])
        groups.append(dict(events=ev,states=states))
    return groups


class StudyTests(unittest.TestCase):
    def test_membership(self):
        self.assertEqual([i for group in MEMBERS for i in group],list(range(1,25)))
        fine=np.arange(24)
        coarse=3*(fine//6)+(fine%6)//2+1
        for b,members in enumerate(MEMBERS,1):
            self.assertEqual((fine[coarse==b]+1).tolist(),members)

    def test_baseline_and_dilution_response(self):
        groups=synthetic()
        base=fit(groups,(.35,.03))
        old=original(groups,(.35,.03),fixed_uu=(0,0))
        self.assertEqual(base['status'],'valid')
        np.testing.assert_allclose(base['values'],old['values'],atol=3e-4,rtol=0)
        np.testing.assert_allclose(base['errors'],old['errors'],atol=1e-4,rtol=0)
        fixed=fit(groups,(.35,.03),.2)
        plus=fit(groups,(.3501,0),.2)
        minus=fit(groups,(.3499,0),.2)
        self.assertAlmostEqual(fixed['values'][6],.2)
        self.assertEqual(fixed['errors'][6],0)
        numeric=(np.array(plus['values'])-minus['values'])/.0002
        np.testing.assert_allclose(fixed['dilution_derivative'],numeric,atol=.003,rtol=0)
        self.assertGreater(np.min(np.linalg.eigvalsh(np.array(fixed['covariance'])[2:6,2:6])),0)

    def test_spawn_worker(self):
        groups={str(i):dict(events=dict(g['events'],bin=np.ones(len(g['events']['phi']),dtype=int)),states=g['states']) for i,g in enumerate(synthetic())}
        settings={'step1_2':{'1':dict(fixed_LLcos=.2)}}
        with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn'),initializer=initialize,initargs=(groups,{1:(.35,.03)},settings)) as pool:
            result=pool.submit(worker,('step1_2',1)).result(timeout=60)
        self.assertEqual(result['status'],'valid')
        self.assertAlmostEqual(result['values'][6],.2)


if __name__=='__main__':
    unittest.main()
