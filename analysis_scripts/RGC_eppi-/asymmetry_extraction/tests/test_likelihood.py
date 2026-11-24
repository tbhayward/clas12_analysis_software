import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from likelihood import ConditionalLikelihood,basis,fit
from inputs import states_for


def toy(seed=51,n=3000,polarized=True):
 rng=np.random.default_rng(seed)
 # Deliberately unequal charges and signed target polarizations.
 states=np.array([[2,.84,.27,.84*.27],[1,-.84,.27,-.84*.27],[1.3,.84,-.22,-.84*.22],[1.8,-.84,-.22,.84*.22]])
 if not polarized:states[:,2:]=0
 phi=rng.uniform(0,2*np.pi,n)
 # A nonuniform common acceptance, independent of the spin state.
 phi=phi[rng.random(n)<np.exp(.5*np.cos(phi)-.5)]
 n=len(phi)
 rb=np.full(n,.6);rc=np.full(n,.8);rv=np.full(n,1.1);rw=np.full(n,.5)
 truth=np.array([.2,-.12,.25,.14,-.07,.21,.08])
 probabilities=np.column_stack([s[0]*(1+basis(phi,rb,rc,rv,rw,*s[1:])@truth) for s in states])
 probabilities/=probabilities.sum(axis=1,keepdims=True)
 choice=(rng.random(n)[:,None]>np.cumsum(probabilities,axis=1)).sum(axis=1)
 chosen=states[choice]
 h=np.sign(chosen[:,1]);pb=np.abs(chosen[:,1]);pt=chosen[:,2]
 return dict(events=dict(phi=phi,rb=rb,rc=rc,rv=rv,rw=rw,h=h,pb=pb,pt=pt),states=states),truth


class TestLikelihood(unittest.TestCase):
 def test_charge_normalization_and_gradient(self):
  g,truth=toy();cost=ConditionalLikelihood([g])
  eps=1e-6
  numeric=[]
  for j in range(7):
   step=np.eye(7)[j]*eps
   numeric.append((cost(*(truth+step))-cost(*(truth-step)))/(2*eps))
  np.testing.assert_allclose(cost.gradient(*truth),numeric,atol=1e-6)
  # Explicit state-by-state normalization must equal charge-moment reduction.
  a=g['events'];s=g['states'];q=s[:,0]/s[:,0].sum()
  factors=np.column_stack([1+basis(*[a[k] for k in ['phi','rb','rc','rv','rw']],*row[1:])@truth for row in s])
  np.testing.assert_allclose(factors@q,1+cost.average@truth)

 def test_independent_period_normalization(self):
  rgc,truth=toy(seed=1)
  rgb,_=toy(seed=2,polarized=False)
  joint=ConditionalLikelihood([rgc,rgb])
  self.assertAlmostEqual(joint(*truth),ConditionalLikelihood([rgc])(*truth)+ConditionalLikelihood([rgb])(*truth))
  rgb_scaled=dict(events=rgb['events'],states=rgb['states'].copy())
  rgb_scaled['states'][:,0]*=10
  self.assertAlmostEqual(joint(*truth),ConditionalLikelihood([rgc,rgb_scaled])(*truth))

 def test_rgb_has_no_target_sensitivity(self):
  g,truth=toy(polarized=False)
  c=ConditionalLikelihood([g],target=False)
  self.assertEqual(c.npar,3)
  full=ConditionalLikelihood([g],target=True)
  np.testing.assert_array_equal(full.observed[:,3:],0)
  self.assertAlmostEqual(full(*truth),c(*truth[:3]))
  # With all spin terms zero, UU factors cancel exactly.
  self.assertAlmostEqual(c(.2,-.1,0),0,places=10)
  self.assertAlmostEqual(c(-.2,.1,0),0,places=10)

 def test_closure_and_shared_dilution_covariance(self):
  g,truth=toy(n=40000)
  result=fit([g],(.35,.02))
  self.assertEqual(result['status'],'valid')
  expected=truth.copy();expected[3:]/=.35
  values=np.array(result['values']);errors=np.array(result['errors'])
  self.assertTrue(np.all(np.abs(values-expected)<4*errors))
  covariance=np.array(result['covariance']);cb=np.array(result['diluted_covariance'])
  jac=np.r_[np.ones(3),np.ones(4)/.35]
  d=np.r_[np.zeros(3),-values[3:]/.35]
  np.testing.assert_allclose(covariance,cb*jac[:,None]*jac[None,:]+.02**2*np.outer(d,d))

 def test_rgb_total_charge_assumption(self):
  states,records=states_for('RGB_spring2019_D2',[6265],{6265:[100.,0,0,0,0,0]})
  np.testing.assert_array_equal(states[:,0],[50.,50.])
  np.testing.assert_array_equal(states[:,2:],0)
  self.assertEqual(records[6265]['pt'],0)

if __name__=='__main__':unittest.main()
