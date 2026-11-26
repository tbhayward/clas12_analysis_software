"""Period-conditional spin likelihood with event-level depolarization ratios.

Fit diluted target amplitudes first. A shared Gaussian dilution nuisance is
then propagated analytically through beta_target = Df * physics_target.
This is the local-Hessian equivalent of profiling that nuisance, without an
unnecessary flat direction in the numerical minimization.
"""
import numpy as np
from iminuit import Minuit

NAMES=['UU_cos_phi','UU_cos_2phi','LU_sin_phi','UL_lab_sin_phi','UL_lab_sin_2phi','LL_lab','LL_lab_cos_phi']


def basis(phi, rb, rc, rv, rw, hb, pt, hbpt):
    return np.column_stack([rv*np.cos(phi),rb*np.cos(2*phi),hb*rw*np.sin(phi),
      pt*rv*np.sin(phi),pt*rb*np.sin(2*phi),hbpt*rc,hbpt*rw*np.cos(phi)])


def extreme_states(states):
    """Keep only convex-hull vertices. This preserves the exact positivity domain."""
    from scipy.spatial import ConvexHull
    x=np.unique(states,axis=0)
    if len(x)<3:return x
    centered=x-x[0]
    _,singular,vectors=np.linalg.svd(centered,full_matrices=False)
    rank=int(np.sum(singular>max(singular[0]*1e-12,1e-14)))
    if rank==0:return x[:1]
    projected=centered@vectors[:rank].T
    if rank==1:return x[[np.argmin(projected[:,0]),np.argmax(projected[:,0])]]
    return x[ConvexHull(projected).vertices]


class ConditionalLikelihood:
    def __init__(self, groups, target=True):
        self.npar=7 if target else 3
        observed=[];average=[];allowed=[]
        for g in groups:
            a=g['events'];states=g['states'];q=states[:,0];w=q/q.sum()
            # State columns are Q, h*Pb, Pt, h*Pb*Pt.
            args=[a[k] for k in ['phi','rb','rc','rv','rw']]
            observed.append(basis(*args,a['h']*a['pb'],a['pt'],a['h']*a['pb']*a['pt'])[:,:self.npar])
            moments=w@states[:,1:]
            average.append(basis(*args,*moments)[:,:self.npar])
            for state in extreme_states(states[:,1:]):
                allowed.append(basis(*args,*state)[:,:self.npar])
        self.observed=np.concatenate(observed);self.average=np.concatenate(average)
        self.allowed=np.concatenate(allowed)
        self.events=len(self.observed)

    def __call__(self,*beta):
        beta=np.asarray(beta)
        if np.min(1+self.allowed@beta)<=1e-10:return 1e30
        n=1+self.observed@beta;d=1+self.average@beta
        if np.any(d<=0):return 1e30
        return float(np.log(d).sum()-np.log(n).sum())

    def gradient(self,*beta):
        beta=np.asarray(beta)
        if np.min(1+self.allowed@beta)<=1e-10:return np.zeros(self.npar)
        n=1+self.observed@beta;d=1+self.average@beta
        return np.sum(self.average/d[:,None]-self.observed/n[:,None],axis=0)


def fit(groups,dilution,starts=3,fixed_uu=None):
    target=any(np.any(g['events']['pt']!=0) for g in groups)
    cost=ConditionalLikelihood(groups,target)
    seeds=[np.array([0,0,.05,.01,.01,.01,.01]),np.array([.3,-.15,.08,-.01,.01,.01,-.01]),np.array([-.3,.15,-.05,.01,-.01,-.01,.01])]
    if fixed_uu is not None:
        fixed_uu=np.asarray(fixed_uu,dtype=float)
        if fixed_uu.shape!=(2,) or not np.all(np.isfinite(fixed_uu)):raise ValueError('Need two finite fixed UU ratios')
        baseline=np.zeros(cost.npar);baseline[:2]=fixed_uu
        if np.min(1+cost.allowed@baseline)<=1e-10:raise ValueError('Imported UU values give a nonpositive unpolarized factor')
        seeds[0]=np.zeros(7)
    fits=[]
    for seed in seeds[:starts]:
        seed=seed.copy()
        if fixed_uu is not None:seed[:2]=fixed_uu
        m=Minuit(cost,*seed[:cost.npar],name=NAMES[:cost.npar],grad=cost.gradient)
        m.errordef=.5;m.strategy=1
        if fixed_uu is not None:m.fixed[0]=True;m.fixed[1]=True
        m.migrad(ncall=4000);m.hesse()
        if not (m.valid and m.accurate):
            m.simplex(ncall=2000)
            m.strategy=2;m.migrad(ncall=6000);m.hesse()
        fits.append(m)
    valid=[m for m in fits if m.valid and m.accurate and m.covariance is not None]
    if not valid:
        return dict(status='invalid',events=cost.events,candidates=[dict(nll=float(m.fval),valid=bool(m.valid),accurate=bool(m.accurate),edm=float(m.fmin.edm),values=list(m.values),minimum_state_factor=float(np.min(1+cost.allowed@np.array(m.values)))) for m in fits])
    best=min(valid,key=lambda m:m.fval)
    beta=np.array(best.values);cb=np.asarray(best.covariance)
    values=beta.copy();jac=np.ones(len(beta));df,sigma=dilution
    if target:values[3:]/=df;jac[3:]/=df
    cov=cb*jac[:,None]*jac[None,:]
    if target:
        derivative=np.zeros(7);derivative[3:]=-values[3:]/df
        cov+=sigma*sigma*np.outer(derivative,derivative)
    free=np.array([not best.fixed[i] for i in range(cost.npar)])
    eigen=np.linalg.eigvalsh(cb[np.ix_(free,free)])
    normalization=np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
    correlation=np.divide(cov,normalization,out=np.zeros_like(cov),where=normalization>0)
    warnings=[]
    if eigen[-1]/eigen[0]>1e8:warnings.append('weak_parameter_identification')
    if np.min(1+cost.allowed@beta)<1e-5:warnings.append('near_positivity_boundary')
    if len(valid)>1 and max(m.fval for m in valid)-best.fval>.1:warnings.append('multistart_minima_disagree')
    derivative=np.zeros(len(beta))
    if target:derivative[3:]=-values[3:]/df
    return dict(status='valid',fixed_parameters={NAMES[i]:float(beta[i]) for i in range(cost.npar) if not free[i]},free_parameters=[NAMES[i] for i in range(cost.npar) if free[i]],events=cost.events,names=NAMES[:cost.npar],values=values.tolist(),
      errors=np.sqrt(np.diag(cov)).tolist(),covariance=cov.tolist(),correlation=correlation.tolist(),
      diluted_values=beta.tolist(),diluted_covariance=cb.tolist(),
      dilution=df if target else None,dilution_stat=sigma if target else None,dilution_derivative=derivative.tolist(),
      nll=float(best.fval),edm=float(best.fmin.edm),
      minimum_state_factor=float(np.min(1+cost.allowed@beta)),
      covariance_condition=float(eigen[-1]/eigen[0]),
      multistart_nll=[float(m.fval) for m in fits],
      warnings=warnings)
