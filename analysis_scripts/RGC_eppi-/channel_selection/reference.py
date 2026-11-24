"""Explicit cross-period shape assistance for weak Su22 calibration bins."""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import gamma
from mass_model import expected_bins,deviance_residual
from models import MOM_RANGE,momentum_expected,PROB_LOW,PROB_HIGH,fit_with_checks

REFERENCE_GROUP='RGC_fall2022_spring2023'

def adjustment(own,reference,kind):
    donors=[i for i,(a,b) in enumerate(zip(own,reference)) if a['usable'] and b['usable'] and not a['warnings'] and not b['warnings']]
    if len(donors)<6:donors=[i for i,(a,b) in enumerate(zip(own,reference)) if a['usable'] and b['usable']]
    if len(donors)<6:return dict(donor_bins=[i+1 for i in donors],mean_shift=0.,width_ratio=1.,mode_ratio=1.,correction_status='insufficient donors, unshifted reference',shift_mad=None,ratio_mad=None)
    if kind=='mass':
        shifts=np.array([own[i]['mu']-reference[i]['mu'] for i in donors]);ratios=np.array([own[i]['sigma']/reference[i]['sigma'] for i in donors])
        shift=float(np.median(shifts));ratio=float(np.median(ratios))
        return dict(donor_bins=[i+1 for i in donors],mean_shift=shift,width_ratio=ratio,mode_ratio=1.,correction_status='median Su22 minus reference mean and width ratio',shift_mad=float(np.median(np.abs(shifts-shift))),ratio_mad=float(np.median(np.abs(ratios-ratio))))
    ratios=np.array([own[i]['mode']/reference[i]['mode'] for i in donors]);ratio=float(np.median(ratios))
    return dict(donor_bins=[i+1 for i in donors],mean_shift=0.,width_ratio=1.,mode_ratio=ratio,correction_status='median Su22/reference signal-mode ratio',shift_mad=None,ratio_mad=float(np.median(np.abs(ratios-ratio))))

def transfer(hist,edges,kind,ref,correction,original,bin_number):
    x=(edges[:-1]+edges[1:])/2
    if kind=='mass':
        mu=ref['mu']+correction['mean_shift'];sigma=ref['sigma']*correction['width_ratio']
        region=(max(float(edges[0]),min(.70,mu-2*sigma-.02)),min(float(edges[-1]),max(1.30,mu+2*sigma+.02)))
    else:region=MOM_RANGE
    ids=np.flatnonzero((x>=region[0])&(x<region[1]));local_edges=edges[ids[0]:ids[-1]+2];y=hist[ids].astype(float);total=max(float(y.sum()),1.)
    if kind=='mass':
        def unpack(q):return [q[0],mu,sigma,q[1],(q[1]+q[2])/2,q[2]]
        def model(q):return expected_bins(local_edges,unpack(q),region)
        bounds=(np.zeros(3),np.full(3,max(10.,y.max()*10)))
        initial=np.array([max(1.,y.max()*.7),max(.1,y[0]),max(.1,y[-1])])
        cut_low,cut_high=mu-2*sigma,mu+2*sigma
    else:
        mode=ref['mode']*correction['mode_ratio'];shape=ref['shape'];scale=mode/(shape-1)
        mu=shape*scale;sigma=np.sqrt(shape)*scale
        def unpack(q):return [q[0],mode,shape,q[1],q[2]]
        def model(q):return momentum_expected(local_edges,unpack(q),3)
        bounds=([0,0,.04],[total*5,total*5,.50]);initial=[total*.7,total*.3,.15]
        cut_low,cut_high=gamma.ppf([PROB_LOW,PROB_HIGH],shape,scale=scale)
    opt=least_squares(lambda q:deviance_residual(y,model(q)),initial,bounds=bounds,max_nfev=1500)
    pred=model(opt.x);ndf=len(y)-len(opt.x);chi=float(np.sum((y-pred)**2/np.maximum(pred,1e-12)))
    result=dict(kind=kind,model='reference_fixed_shape_mass' if kind=='mass' else 'gamma_signal_plus_gamma_tail',background_shape=3,usable=True,origin='reference_assisted',parameters=np.asarray(unpack(opt.x)).tolist(),mu=float(mu),sigma=float(sigma),cut_low=float(cut_low),cut_high=float(cut_high),fit_range=list(region),fit_entries=int(y.sum()),chi2=chi,ndf=ndf,chi2_ndf=chi/ndf,deviance=float(np.dot(opt.fun,opt.fun)),warnings=['reference-assisted shape, not an independent Su22 peak determination'],reason='nominal preliminary reference-assisted cut',reference=dict(group=REFERENCE_GROUP,bin=bin_number,correction=correction,unadjusted_window=[ref['cut_low'],ref['cut_high']]),independent_fit=original,checks=dict(unadjusted_reference=ref))
    if kind=='momentum':result.update(mode=float(mode),shape=float(shape),cut_method='central 95.45% of transferred signal gamma distribution')
    result['reference_window_shift']=float(max(abs(cut_low-ref['cut_low']),abs(cut_high-ref['cut_high'])))
    if not opt.success:result['warnings'].append('Su22 amplitude/background fit did not converge; window still comes from reference')
    return result

def assist(histograms,edges,kind,own,reference):
    # Retry only unresolved bins with the same-bin reference as an optimizer seed.
    for i,f in enumerate(own):
        f.setdefault('origin','independent')
        if not f['usable'] and reference[i]['usable']:
            retry=fit_with_checks(histograms[i],edges,kind,reference=reference[i]);retry['origin']='reference_seeded_independent'
            retry['reference_seed_group']=REFERENCE_GROUP
            if retry['usable']:own[i]=retry
    correction=adjustment(own,reference,kind)
    for i,f in enumerate(own):
        if not f['usable'] and reference[i]['usable']:own[i]=transfer(histograms[i],edges,kind,reference[i],correction,f,i+1)
    return own
