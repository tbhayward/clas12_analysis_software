"""Nominal preliminary models: Gaussian mass peak and asymmetric momentum peak."""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import gammainc
from scipy.stats import gamma
from mass_model import MASS_EDGES,ANGLE_EDGES,fit_with_checks as mass_checks,fit_angle,components as mass_components
from mass_model import deviance_residual,merge_pairs

MOM_EDGES=np.linspace(0,1.0,201)
MOM_RANGE=(0,.6)
PROB_LOW=.02275013194817921
PROB_HIGH=1-PROB_LOW

def momentum_expected(edges,p,background_shape=3):
    signal_yield,mode,shape,bg_yield,gap=p
    scale=mode/(shape-1)
    bg_scale=(mode+gap)/(background_shape-1)
    return signal_yield*np.diff(gammainc(shape,edges/scale))+bg_yield*np.diff(gammainc(background_shape,edges/bg_scale))

def momentum_components(x,p,width,background_shape=3):
    ns,mode,shape,nb,gap=p
    return width*ns*gamma.pdf(x,shape,scale=mode/(shape-1)),width*nb*gamma.pdf(x,background_shape,scale=(mode+gap)/(background_shape-1))

def fit_momentum(hist,edges,background_shape=3,reference=None):
    centers=(edges[:-1]+edges[1:])/2
    select=(centers>=MOM_RANGE[0])&(centers<MOM_RANGE[1])
    ids=np.flatnonzero(select);edges_local=edges[ids[0]:ids[-1]+2]
    y=np.asarray(hist,float)[select]
    result=dict(kind='momentum',model='gamma_signal_plus_gamma_tail',background_shape=background_shape,usable=False,warnings=[],fit_range=list(MOM_RANGE),fit_entries=int(y.sum()),cut_method='central 95.45% of fitted signal gamma distribution')
    if y.sum()<100:return dict(result,reason='fewer than 100 entries in fit range')
    total=float(y.sum())
    low=[0,.02,2.0,0,.04]
    high=[total*5,.20,15.,total*5,.50]
    def residual(p):return deviance_residual(y,momentum_expected(edges_local,p,background_shape))
    fits=[]
    if reference is not None:
        q=reference['parameters'];start=[total*.7,float(np.clip(q[1],.0201,.1999)),float(np.clip(q[2],2.01,14.99)),total*.3,float(np.clip(q[4],.0401,.4999))]
        fits.append(least_squares(residual,start,bounds=(low,high),max_nfev=1500))
    for mode in [.06,.10,.14]:
        for shape in [3.,6.]:
            p=[total*.7,mode,shape,total*.3,.15]
            fits.append(least_squares(residual,p,bounds=(low,high),max_nfev=1500))
    opt=min(fits,key=lambda f:np.dot(f.fun,f.fun));p=opt.x
    expected=momentum_expected(edges_local,p,background_shape)
    scale=p[1]/(p[2]-1)
    lower,upper=gamma.ppf([PROB_LOW,PROB_HIGH],p[2],scale=scale)
    chi=float(np.sum((y-expected)**2/np.maximum(expected,1e-12)))
    result.update(parameters=p.tolist(),mode=float(p[1]),shape=float(p[2]),mu=float(p[2]*scale),sigma=float(np.sqrt(p[2])*scale),cut_low=float(lower),cut_high=float(upper),deviance=float(np.dot(opt.fun,opt.fun)),chi2=chi,ndf=len(y)-5,chi2_ndf=chi/(len(y)-5))
    if not opt.success:return dict(result,reason='fit did not converge')
    if any(min(p[i]-low[i],high[i]-p[i])/(high[i]-low[i])<.005 for i in [1,2]):
        return dict(result,reason='signal mode or shape reached a bound')
    # A broad single tail must not substitute for evidence of the narrow component.
    def bg_residual(q):
        nb,mode=q
        pred=nb*np.diff(gammainc(background_shape,edges_local/(mode/(background_shape-1))))
        return deviance_residual(y,pred)
    bg=least_squares(bg_residual,[total,.20],bounds=([0,.06],[total*5,.70]),max_nfev=1500)
    result['deviance_improvement']=float(np.dot(bg.fun,bg.fun)-result['deviance'])
    if result['deviance_improvement']<9:return dict(result,reason='narrow component not sufficiently resolved above broad tail')
    result['usable']=True
    if result['chi2_ndf']>5:result['warnings'].append('poor signal-plus-tail description')
    if upper>MOM_RANGE[1]:result['warnings'].append('cut extends outside fitted range')
    result['reason']='nominal preliminary cut'
    return result

def fit_with_checks(hist,edges,kind,reference=None):
    if kind=='mass':return mass_checks(hist,edges,kind,initial_peak=(reference['mu'],reference['sigma']) if reference else None)
    nominal=fit_momentum(hist,edges,reference=reference)
    coarse,ce=merge_pairs(hist,edges)
    nominal['checks']=dict(tail_shape_2=fit_momentum(hist,edges,2),coarser_histogram=fit_momentum(coarse,ce))
    shifts=[]
    if nominal['usable']:
        for name,f in nominal['checks'].items():
            if not f['usable']:nominal['warnings'].append(name+' check unresolved');continue
            shifts.append(max(abs(f['cut_low']-nominal['cut_low']),abs(f['cut_high']-nominal['cut_high'])))
        threshold=max(2*(edges[1]-edges[0]),.1*(nominal['cut_high']-nominal['cut_low']))
        nominal['max_boundary_shift']=float(max(shifts,default=0));nominal['stability_threshold']=float(threshold)
        if any(s>threshold for s in shifts):nominal['warnings'].append('cut unstable under model or histogram change')
    return nominal

def components(x,fit,width):
    if fit['kind']=='mass':return mass_components(x,fit['parameters'],fit['fit_range'])
    return momentum_components(x,fit['parameters'],width,fit['background_shape'])
