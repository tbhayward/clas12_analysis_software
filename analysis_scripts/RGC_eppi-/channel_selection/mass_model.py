"""Positive linear background baseline with concave-down quadratic checks."""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import xlogy, ndtr

MASS_EDGES = np.linspace(.30, 1.60, 131)  # 0.010 GeV^2, twice the original width
MOM_EDGES = np.linspace(0, .50, 101)  # 0.005 GeV
ANGLE_EDGES = np.linspace(-6, 6, 51)  # degrees, measured minus predicted
FIT_RANGES = {'mass': (.70, 1.30), 'momentum': (0, .25)}
RANGE_LIMITS = {'mass': (.40, 1.50), 'momentum': (0, .45)}

def components(x, p, region):
    amp, mu, sigma, b0, b1, b2 = p
    u = (x-region[0])/(region[1]-region[0])
    return amp*np.exp(-.5*((x-mu)/sigma)**2), b0*(1-u)**2+2*b1*u*(1-u)+b2*u*u

def expected_bins(edges, p, region):
    """Bin-averaged model, avoiding a coarse-bin center approximation."""
    amp,mu,sigma,b0,b1,b2=p
    lo,hi=edges[:-1],edges[1:]
    gaussian=amp*sigma*np.sqrt(2*np.pi)*(ndtr((hi-mu)/sigma)-ndtr((lo-mu)/sigma))/(hi-lo)
    u0,u1=(lo-region[0])/(region[1]-region[0]),(hi-region[0])/(region[1]-region[0])
    bg=b0+2*(b1-b0)*(u0+u1)/2+(b0-2*b1+b2)*(u0*u0+u0*u1+u1*u1)/3
    return gaussian+bg

def canonical(p, model):
    amp,mu,sigma,b0,b2 = p[:5]
    return np.array([amp,mu,sigma,b0,(b0+b2)/2+(p[5] if model=='concave_quadratic' else 0),b2])

def deviance_residual(obs, expected):
    expected = np.maximum(expected, 1e-12)
    return np.sign(expected-obs)*np.sqrt(np.maximum(0, 2*(expected-obs+xlogy(obs, obs/expected))))

def fit_histogram(hist, edges, kind, model='linear', region=None, initial_peak=None):
    region=tuple(region or FIT_RANGES[kind])
    xall=(edges[:-1]+edges[1:])/2
    sel=(xall>=region[0])&(xall<region[1])
    x,y=xall[sel],np.asarray(hist,float)[sel]
    indices=np.flatnonzero(sel)
    local_edges=edges[indices[0]:indices[-1]+2]
    result=dict(kind=kind,model=model,usable=False,warnings=[],fit_range=list(region),fit_entries=int(y.sum()))
    if y.sum()<100:
        return dict(result,reason='fewer than 100 entries in fit range')
    upper_y=max(10.,y.max()*10)
    mb,sb,seeds=((.78,1.05),(.01,.20),[.86,.92,.98]) if kind=='mass' else ((.005,.18),(.005,.15),[.05,.08,.12])
    if initial_peak is not None:seeds=[float(np.clip(initial_peak[0],mb[0]+1e-5,mb[1]-1e-5)),*seeds]
    low=[0,mb[0],sb[0],0,0]+([0] if model=='concave_quadratic' else [])
    high=[upper_y,mb[1],sb[1],upper_y,upper_y]+([upper_y] if model=='concave_quadratic' else [])
    def residual(p):
        return deviance_residual(y,expected_bins(local_edges,canonical(p,model),region))
    fits=[]
    for seed in seeds:
        initial=[max(1,y.max()*.7),seed,float(np.clip(initial_peak[1],sb[0]+1e-5,sb[1]-1e-5)) if initial_peak is not None else (.07 if kind=='mass' else .035),max(.1,y[0]),max(.1,y[-1])]+([max(.1,y.max()*.05)] if model=='concave_quadratic' else [])
        fits.append(least_squares(residual,initial,bounds=(low,high),max_nfev=1500))
    opt=min(fits,key=lambda f:np.dot(f.fun,f.fun))
    pars=canonical(opt.x,model)
    expected=expected_bins(local_edges,pars,region)
    ndf=len(x)-len(opt.x)
    chi=float(np.sum((y-expected)**2/np.maximum(expected,1e-12)))
    errors=np.sqrt(np.maximum(0,np.diag(np.linalg.pinv(opt.jac.T@opt.jac))))
    result.update(parameters=pars.tolist(),mu=float(pars[1]),sigma=float(pars[2]),
                  mu_error=float(errors[1]),sigma_error=float(errors[2]),
                  cut_low=float(max(0,pars[1]-2*pars[2]) if kind=='momentum' else pars[1]-2*pars[2]),cut_high=float(pars[1]+2*pars[2]),
                  deviance=float(np.dot(opt.fun,opt.fun)),chi2=chi,ndf=ndf,chi2_ndf=chi/ndf)
    if not opt.success:return dict(result,reason='fit did not converge')
    if any(min(opt.x[i]-low[i],high[i]-opt.x[i])/(high[i]-low[i])<.005 for i in [1,2]):
        return dict(result,reason='peak mean or width reached a bound')
    def bg_residual(b):
        p=np.r_[0,pars[1],pars[2],b]
        return deviance_residual(y,expected_bins(local_edges,canonical(p,model),region))
    bgfit=least_squares(bg_residual,opt.x[3:],bounds=(0,upper_y),max_nfev=1500)
    result['deviance_improvement']=float(np.dot(bgfit.fun,bgfit.fun)-result['deviance'])
    if result['deviance_improvement']<9:return dict(result,reason='peak not sufficiently resolved above background')
    result['usable']=True
    if chi/ndf>5:result['warnings'].append('poor Gaussian-plus-background description')
    if result['cut_low']<region[0]-1e-10 or result['cut_high']>region[1]+1e-10:
        result['warnings'].append('cut extends outside fitted range')
    result['reason']='provisional cut'
    return result

def adaptive_fit(hist,edges,kind,initial_peak=None):
    region=FIT_RANGES[kind]
    history=[]
    step=float(edges[1]-edges[0])
    for _ in range(4):
        fit=fit_histogram(hist,edges,kind,region=region,initial_peak=initial_peak)
        history.append(list(region))
        if not fit['usable'] or 'cut extends outside fitted range' not in fit['warnings']:break
        limits=RANGE_LIMITS[kind]
        # Expand only as needed, with bin-aligned edges and two histogram bins of margin.
        lo=max(limits[0],min(region[0],np.floor((fit['cut_low']-2*step)/step)*step))
        hi=min(limits[1],max(region[1],np.ceil((fit['cut_high']+2*step)/step)*step))
        new=(round(lo,8),round(hi,8))
        if np.allclose(new,region):break
        region=new
    fit['range_history']=history
    return fit

def merge_pairs(hist,edges):
    # Preserve both edges even for an odd number of bins.
    starts=np.arange(0,len(hist),2)
    return np.add.reduceat(hist,starts),np.r_[edges[starts],edges[-1]]

def fit_with_checks(hist,edges,kind,initial_peak=None):
    nominal=adaptive_fit(hist,edges,kind,initial_peak)
    region=nominal['fit_range']
    alternate=fit_histogram(hist,edges,kind,model='concave_quadratic',region=region)
    rebinned,redges=merge_pairs(hist,edges)
    rebin_fit=fit_histogram(rebinned,redges,kind,region=region)
    nominal['checks']=dict(concave_quadratic=alternate,coarser_histogram=rebin_fit)
    shifts=[]
    if nominal['usable']:
        for label,fit in nominal['checks'].items():
            if not fit['usable']:
                nominal['warnings'].append(label+' check unresolved')
                continue
            shift=max(abs(nominal['cut_low']-fit['cut_low']),abs(nominal['cut_high']-fit['cut_high']))
            nominal[label+'_max_boundary_shift']=float(shift)
            shifts.append(shift)
        threshold=max(2*(edges[1]-edges[0]),.1*(nominal['cut_high']-nominal['cut_low']))
        nominal['stability_threshold']=float(threshold)
        nominal['max_boundary_shift']=float(max(shifts,default=0))
        if any(s>threshold for s in shifts):nominal['warnings'].append('cut unstable under background or histogram change')
    return nominal

def fit_angle(hist):
    """Gaussian visual guide in degrees. No angular cut is supplied."""
    from scipy.optimize import curve_fit
    x=(ANGLE_EDGES[:-1]+ANGLE_EDGES[1:])/2
    sel=np.abs(x)<3.5
    x,y=x[sel],hist[sel].astype(float)
    if y.sum()<100:return dict(usable=False,reason='insufficient entries',unit='degree')
    def gaussian(x,a,mu,sigma):return a*np.exp(-.5*((x-mu)/sigma)**2)
    try:
        p,_=curve_fit(gaussian,x,y,p0=[max(1,y.max()),0,.86],sigma=np.sqrt(np.maximum(1,y)),absolute_sigma=True,bounds=([0,-2.3,.05],[max(10,y.max()*10),2.3,3.5]),maxfev=5000)
        return dict(usable=True,parameters=p.tolist(),mu=float(p[1]),sigma=float(p[2]),unit='degree',role='visual guide only, no angular cut')
    except (RuntimeError,ValueError):return dict(usable=False,reason='Gaussian guide did not converge',unit='degree')
