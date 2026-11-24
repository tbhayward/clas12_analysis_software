"""Fit five calibration groups and freeze reusable, target-independent windows."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from candidate_cache import iterate,save_json
from selection import GROUPS,VERSION,fingerprint,seal
from models import MASS_EDGES,MOM_EDGES,ANGLE_EDGES,fit_with_checks,fit_angle,components
from kinematics import XB,TP

EDGES=dict(mass=MASS_EDGES,momentum=MOM_EDGES,angle=ANGLE_EDGES)

def hist_add(hist,values,bins,edges):
    index=np.searchsorted(edges,values,side='right')-1
    good=np.isfinite(values)&(index>=0)&(index<len(edges)-1)
    hist+=np.bincount(bins[good]*(len(edges)-1)+index[good],minlength=hist.size).reshape(hist.shape)

def window(values,bins,fits):
    lo=np.array([f['cut_low'] if f['usable'] else np.inf for f in fits])
    hi=np.array([f['cut_high'] if f['usable'] else -np.inf for f in fits])
    return (values>=lo[bins])&(values<=hi[bins])

def group_histograms(members,caches,stage,mass=None,momentum=None):
    hist=np.zeros((24,len(EDGES[stage])-1),dtype=np.int64)
    counts={name:np.zeros(24,dtype=np.int64) for name in members}
    for name in members:
        path,meta=caches[name]
        for a in iterate(Path(path),meta):
            bins=a['bin'].astype(np.int64)
            good=np.ones(len(bins),dtype=bool)
            if mass is not None:good &= window(a['mass'],bins,mass)
            if momentum is not None:good &= window(a['momentum'],bins,momentum)
            hist_add(hist,a[stage][good],bins[good],EDGES[stage])
            counts[name]+=np.bincount(bins[good],minlength=24)
    return hist,{k:v.tolist() for k,v in counts.items()}

def fit_group(job):
    group,members,caches,output=job[:4]
    reference=job[4] if len(job)>4 else None
    out=Path(output)/'fits'/group;out.mkdir(parents=True,exist_ok=True)
    print('FIT '+group,flush=True)
    save_json(out/'progress.json',dict(stage='mass fits',state='fitting'))
    mh,initial=group_histograms(members,caches,'mass')
    mass=[fit_with_checks(h,MASS_EDGES,'mass') for h in mh]
    if reference is not None:
        from reference import assist
        mass=assist(mh,MASS_EDGES,'mass',mass,reference['mass'])
    ph,after_mass=group_histograms(members,caches,'momentum',mass)
    save_json(out/'progress.json',dict(stage='momentum fits',state='fitting'))
    momentum=[fit_with_checks(h,MOM_EDGES,'momentum') for h in ph]
    if reference is not None:momentum=assist(ph,MOM_EDGES,'momentum',momentum,reference['momentum'])
    ah,after_both=group_histograms(members,caches,'angle',mass,momentum)
    guides=[fit_angle(h) for h in ah]
    result=dict(group=group,members=members,mass=mass,momentum=momentum,angle_guides=guides,counts=dict(initial=initial,after_mass=after_mass,after_both=after_both))
    save_json(out/'fit_results.json',result)
    np.savez_compressed(out/'histograms.npz',mass=mh,momentum=ph,angle=ah)
    write_tables(out,result)
    plot_group(out,result)
    save_json(out/'progress.json',dict(stage='complete',state='complete'))
    print('FITTED '+group,flush=True)
    return result

def write_tables(out,result):
    with (out/'cuts.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['bin','stage','origin','usable','model','mean','width','low','high','chi2_ndf','max_boundary_shift','warnings','reason'])
        for b in range(24):
            for stage in ['mass','momentum']:
                f=result[stage][b]
                w.writerow([b+1,stage,f.get('origin','independent'),f['usable'],f['model'],*[f.get(k,'') for k in ['mu','sigma','cut_low','cut_high','chi2_ndf','max_boundary_shift']],' | '.join(f['warnings']),f['reason']])
    with (out/'stability.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['bin','stage','variant','usable','low','high','chi2_ndf','reason'])
        for b in range(24):
            for stage in ['mass','momentum']:
                f=result[stage][b]
                for name,v in [('nominal',f),*f['checks'].items()]:w.writerow([b+1,stage,name,v['usable'],v.get('cut_low',''),v.get('cut_high',''),v.get('chi2_ndf',''),v['reason']])

def freeze(results,caches,output):
    if {r['group'] for r in results} != set(GROUPS):raise ValueError('All five calibration groups are required')
    groups={}
    for r in results:
        bins=[]
        for i,(m,p) in enumerate(zip(r['mass'],r['momentum'])):
            active=m['usable'] and p['usable']
            v=dict(bin=i+1,active=active,mass_origin=m.get('origin','independent'),momentum_origin=p.get('origin','independent'),reference_assisted=any(f.get('origin')=='reference_assisted' for f in [m,p]),warnings=m['warnings']+p['warnings'],reason='usable nominal preliminary cut' if active else 'mass: '+m['reason']+'; momentum: '+p['reason'])
            if m.get('origin')=='reference_assisted':v['mass_reference']=m['reference']
            if p.get('origin')=='reference_assisted':v['momentum_reference']=p['reference']
            if active:v.update(mass_low=m['cut_low'],mass_high=m['cut_high'],momentum_low=p['cut_low'],momentum_high=p['cut_high'])
            bins.append(v)
        groups[r['group']]=dict(calibration_members=r['members'],bins=bins)
    payload=dict(fit_implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name('models.py').read_bytes()+Path(__file__).with_name('mass_model.py').read_bytes()+Path(__file__).with_name('reference.py').read_bytes()).hexdigest(),version=VERSION,group_ids={name:i+1 for i,name in enumerate(GROUPS)},partial=any(meta['partial'] for _,meta in caches.values()),base_selection_sha256=fingerprint(),xB_edges=XB.tolist(),minus_tprime_edges=TP.tolist(),angle_unit='degree',angle_convention='measured minus predicted',groups=groups,calibration_sources={n:m for n,(_,m) in caches.items()},policy='Inactive bins have no nominal cut and must be reported as unavailable, never as physical zero yields. Identical windows apply to ND3 and all control targets of each period.')
    unresolved={g:[b['bin'] for b in v['bins'] if not b['active']] for g,v in groups.items() if any(not b['active'] for b in v['bins'])}
    if unresolved and not payload['partial']:raise ValueError('Cannot publish nominal cut set with unavailable bins: '+str(unresolved))
    payload['cutset_id']=seal(payload)
    path=Path(output)/'nominal_cuts.json';save_json(path,payload)
    return path,payload

def plot_group(out,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for stage,label,limits,subtitle in [('mass',r'$M_{X,d}^{2}$ (GeV$^2$)',(.65,1.35),'Before exclusivity selection'),('momentum',r'$|\vec p_{\rm miss}|$ (GeV)',(0,.60),'After missing-mass selection'),('angle',r'$\theta_{\rm measured}-\theta_{\rm predicted}$ (deg)',(-6,6),'After both selections')]:
        with np.load(out/'histograms.npz') as data:all_hist=data[stage]
        edges=EDGES[stage];width=edges[1]-edges[0]
        for xi in range(4):
            fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
            for ti,ax in enumerate(axes.flat):
                b=xi*6+ti;hist=all_hist[b]
                ax.stairs(hist,edges,color='black',lw=1)
                notes=[];shown=limits
                if stage!='angle':
                    f=result[stage][b]
                    if 'parameters' in f:
                        x=np.linspace(*f['fit_range'],500);g,bg=components(x,f,width)
                        ax.plot(x,g+bg,color='tab:red',lw=1.4);ax.plot(x,bg,color='gray',ls=':',lw=1)
                    if f['usable']:
                        lo,hi=f['cut_low'],f['cut_high'];ax.axvspan(lo,hi,color='tab:green',alpha=.09)
                        for v in [lo,hi]:ax.axvline(v,color='tab:green',ls='--',lw=1)
                        notes=[f'Cut: [{lo:.3f}, {hi:.3f}]']
                        if f.get('origin')=='reference_assisted':notes.append('Reference-assisted shape')
                        elif f['warnings']:notes.append('Model sensitivity: see table')
                        shown=(max(edges[0],min(limits[0],lo-.01,f['fit_range'][0]-.01)),min(edges[-1],max(limits[1],hi+.01,f['fit_range'][1]+.01)))
                    else:notes=['NO NOMINAL CUT','Unavailable bin']
                else:
                    f=result['angle_guides'][b];ax.axvline(0,color='gray',ls='--',lw=1)
                    notes=[f'N = {sum(v[b] for v in result["counts"]["after_both"].values()):,}']
                    if f['usable']:
                        amp,mu,sigma=f['parameters'];x=np.linspace(-3.5,3.5,400)
                        ax.plot(x,amp*np.exp(-.5*((x-mu)/sigma)**2),color='tab:red',lw=1.3)
                        notes.append(rf'$\mu$={mu:.2f}, $\sigma$={sigma:.2f} deg')
                ax.text(.97,.96,'\n'.join(notes),ha='right',va='top',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                ax.set(title=f'Bin {b+1}: {TP[ti]:.2f} ≤ −t′ < {TP[ti+1]:.2f}',xlabel=label,ylabel='Combinations / bin',xlim=shown)
                center=(edges[:-1]+edges[1:])/2;h=hist[(center>=shown[0])&(center<=shown[1])]
                ax.set_ylim(0,max(1,float(h.max())*1.28));ax.spines[['top','right']].set_visible(False)
            fig.suptitle(f'{result["group"]} | {XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f}\n{subtitle}',fontsize=13)
            fig.savefig(out/f'{stage}_xb{xi+1}.png',dpi=150);plt.close(fig)
