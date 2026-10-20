#!/usr/bin/env python3
"""Production-bin |t| template-fit diagnostics. No raw-data efficiency claims.
Run in dvcs_cross_section/external_scripts; requires the original fitter module.
"""
import os
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,csv,json,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
import plot_exclusivity_data_dvcs_pi0_mc as fitlib
MP=0.9382720813
BEAMS={'fa18_inb':10.604,'fa18_out':10.604,'sp18_inb':10.594,'sp18_out':10.594,'sp19_inb':10.200}
VARS=fitlib.VARIABLES

def edges_from_csv(path):
    with open(path,newline='') as f:
        rows=list(csv.DictReader(f))
    return np.array(sorted({float(r[k]) for r in rows for k in ('t_abs_min','t_abs_max') if r.get(k)}))

def vector_t(a,beam):
    # q=k-k', q'=measured photon; t_gamma=(q-q')^2.
    # p' is measured proton, t_proton=(p'-p)^2=-2 M T_p.
    def xyz(prefix):
        p=np.asarray(a[prefix+'_p'],float); th=np.asarray(a[prefix+'_theta'],float); ph=np.asarray(a[prefix+'_phi'],float)
        return np.stack((p*np.sin(th)*np.cos(ph),p*np.sin(th)*np.sin(ph),p*np.cos(th)),axis=0)
    e=xyz('e'); g=xyz('p2'); p=xyz('p1')
    ee=np.sqrt(np.sum(e*e,axis=0)+0.00051099895**2)
    pg=np.sqrt(np.sum(g*g,axis=0)); pp=np.sqrt(np.sum(p*p,axis=0))
    beam_p=np.sqrt(beam**2-0.00051099895**2)
    d=np.stack((-e[0]-g[0],-e[1]-g[1],beam_p-e[2]-g[2]),axis=0)
    tg=(beam-ee-pg)**2-np.sum(d*d,axis=0)
    tp=-2*MP*(np.sqrt(pp*pp+MP*MP)-MP)
    return tp-tg

def read_sample(job):
    pk,tk,kind,path,edges,step,global_json=job
    fitlib.ACTIVE_GLOBAL_CUTS=fitlib.load_global_cuts_config(global_json)
    topo=next(t for t in fitlib.TOPOLOGIES if t.key==tk)
    with uproot.open(f'{path}:{fitlib.TREE_NAME}') as tree:
        resolved=fitlib.resolve_variable_branches(path,VARS)
        required=set(fitlib.global_cut_required_branches())|set(resolved.values())
        required.update(f'{p}_{v}' for p in ('e','p1','p2') for v in ('p','theta','phi'))
        missing=required-set(tree.keys())
        if missing: raise RuntimeError(f'{path}: missing {sorted(missing)}')
    nb=len(edges)-1
    h={v.branch:np.zeros((nb,v.bins),dtype=np.int64) for v in VARS}
    dt=[[] for _ in range(nb)]
    n=np.zeros(nb,dtype=np.int64)
    for a in uproot.iterate(f'{path}:{fitlib.TREE_NAME}',expressions=sorted(required),step_size=step,library='np'):
        d1=np.asarray(a['detector1']);d2=np.asarray(a['detector2'])
        mask=fitlib.global_cut_mask_uproot(a,pk,d1,d2)&(d1==topo.detector1)&(d2==topo.detector2)
        t=-np.asarray(a['t1'],float)
        ix=np.searchsorted(edges,t,side='right')-1
        mask &= np.isfinite(t)&(ix>=0)&(ix<nb)
        if not np.any(mask):continue
        delta=vector_t(a,BEAMS[pk])
        for i in np.unique(ix[mask]):
            sel=mask&(ix==i)
            n[i]+=np.count_nonzero(sel)
            for v in VARS:
                values=np.asarray(a[resolved[v.branch]][sel],float)
                h[v.branch][i]+=np.histogram(values[np.isfinite(values)],bins=v.bins,range=(v.xmin,v.xmax))[0]
            ds=delta[sel];dt[i].append(ds[np.isfinite(ds)].astype(np.float32))
    dt=[np.concatenate(chunks) if chunks else np.array([],dtype=np.float32) for chunks in dt]
    return pk,tk,kind,h,dt,n

def cut_bounds(spec):
    if 'cut_low' in spec and 'cut_high' in spec:return float(spec['cut_low']),float(spec['cut_high'])
    return float(spec['mean'])-3*float(spec['std']),float(spec['mean'])+3*float(spec['std'])

def shape_width(shape,v):
    shape=np.asarray(shape,float)
    if not np.all(np.isfinite(shape)) or np.sum(shape)<=0:return np.nan
    centers=np.linspace(v.xmin+(v.xmax-v.xmin)/(2*v.bins),v.xmax-(v.xmax-v.xmin)/(2*v.bins),v.bins)
    c=np.cumsum(shape)/np.sum(shape)
    return float((np.interp(.84,c,centers)-np.interp(.16,c,centers))/2)

def contained(shape,v,low,high):
    if shape is None:return np.nan
    centers=v.xmin+(np.arange(v.bins)+.5)*(v.xmax-v.xmin)/v.bins
    a=np.asarray(shape,float)
    return float(a[(centers>=low)&(centers<=high)].sum()/a.sum()) if a.sum()>0 else np.nan

def fit_one(job):
    pk,tk,i,edges,inputs,cutmap,opts=job
    topo=next(t for t in fitlib.TOPOLOGIES if t.key==tk)
    dh={v.branch:inputs['data'][0][v.branch][i] for v in VARS}
    mh={v.branch:inputs['mc'][0][v.branch][i] for v in VARS}
    ph={v.branch:inputs['pi0'][0][v.branch][i] for v in VARS}
    result=fitlib.fit_shared_two_templates(dh,mh,ph,topo,
        opts['shift'],opts['smear'],opts['min_counts'],('Delta_phi','theta_gamma_gamma','pTmiss'),
        4.,8.,True,.90,.95,.90,.95)
    rows=[]
    for v in VARS:
        f=(result.variable_results or {}).get(v.branch)
        if f is None or not f.success:continue
        lo,hi=cut_bounds(cutmap[v.branch]) if v.branch in cutmap else (np.nan,np.nan)
        mc=mh[v.branch]/max(1,np.sum(mh[v.branch]))
        data_shape=f.transformed_dvcs_shape
        rows.append(dict(period=pk,topology=tk,t_min=edges[i],t_max=edges[i+1],variable=v.branch,
                         fit_success=result.success,f_pi0=result.f_pi0,f_pi0_error=result.f_pi0_err,
                         deviance=result.deviance,ndf=result.ndf,
                         n_data=inputs['data'][2][i],n_mc=inputs['mc'][2][i],n_pi0=inputs['pi0'][2][i],
                         fitted_data_dvcs_width=shape_width(data_shape,v),mc_dvcs_width=shape_width(mc,v),
                         fitted_data_single_containment=contained(data_shape,v,lo,hi),
                         mc_single_containment=contained(mc,v,lo,hi),
                         cut_low=lo,cut_high=hi))
    return pk,tk,i,result,rows,dh,mh,ph

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--periods',nargs='+',default=['fa18_inb','sp18_out'])
    ap.add_argument('--topologies',nargs='+',default=['CD_FD','CD_FT','FD_FD'])
    ap.add_argument('--step-size',default='100 MB')
    ap.add_argument('--binning',default=str(ROOT/'imports/all_bin_v3.csv'))
    ap.add_argument('--cuts',default=str(ROOT/'output/jsons/combined_cuts.json'))
    ap.add_argument('--global-cuts-json',default=str(ROOT/'output/jsons/global_cuts_config.json'))
    ap.add_argument('--output',default=str(ROOT/'output/t_dependent_exclusivity_diagnostics'))
    ap.add_argument('--fit-min-counts',type=int,default=100)
    ap.add_argument('--max-shift-bins',type=float,default=12.)
    ap.add_argument('--max-smear-bins',type=float,default=20.)
    args=ap.parse_args();edges=edges_from_csv(args.binning)
    with open(args.cuts) as f:cuts=json.load(f)
    periods={p.key:p for p in fitlib.PERIODS}
    jobs=[]
    for pk in args.periods:
        p=periods[pk]
        for tk in args.topologies:
            for kind,path in [('data',p.data_file),('mc',p.dvcs_mc_file),('pi0',p.pi0_as_dvcs_mc_file)]:
                if not Path(path).is_file():raise FileNotFoundError(path)
                jobs.append((pk,tk,kind,path,edges,args.step_size,args.global_cuts_json))
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    inputs={}
    with ProcessPoolExecutor(max_workers=min(8,max(1,args.workers))) as pool:
        for future in as_completed([pool.submit(read_sample,j) for j in jobs]):
            pk,tk,kind,h,dt,n=future.result();inputs[(pk,tk,kind)]=(h,dt,n)
            print(f'[read] {pk} {tk} {kind}: {n.sum():,}',flush=True)
    fitjobs=[]
    for pk in args.periods:
        for tk in args.topologies:
            key=fitlib.cut_json_key('dvcs',periods[pk],next(t for t in fitlib.TOPOLOGIES if t.key==tk))
            if key not in cuts:raise KeyError(key)
            cmap=cuts[key]
            sample={kind:inputs[(pk,tk,kind)] for kind in ('data','mc','pi0')}
            for i in range(len(edges)-1):
                fitjobs.append((pk,tk,i,edges,sample,cmap['data'],dict(shift=args.max_shift_bins,smear=args.max_smear_bins,min_counts=args.fit_min_counts)))
    allrows=[];fits={}
    with ProcessPoolExecutor(max_workers=min(8,max(1,args.workers))) as pool:
        for future in as_completed([pool.submit(fit_one,j) for j in fitjobs]):
            pk,tk,i,result,rows,dh,mh,ph=future.result()
            fits[(pk,tk,i)]=(result,dh,mh,ph)
            allrows.extend(rows)
            print(f'[fit] {pk} {tk} bin {i+1}: success={result.success} f_pi0={result.f_pi0:.3g}',flush=True)
    if allrows:
        with open(out/'template_fit_by_t.csv','w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    for pk in args.periods:
        for tk in args.topologies:
            fig,axs=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
            for ax,v in zip(axs.flat,VARS):
                for kind,label,color in [('fitted_data_dvcs_width','Fitted data DVCS','tab:orange'),('mc_dvcs_width','DVCS MC','tab:blue')]:
                    rr=sorted((r for r in allrows if r['period']==pk and r['topology']==tk and r['variable']==v.branch),key=lambda r:r['t_min'])
                    ax.plot([.5*(r['t_min']+r['t_max']) for r in rr],[r[kind] for r in rr],'o-',label=label,color=color)
                ax.set_title(v.branch);ax.set_xlabel(r'$|t|$ (GeV$^2$)');ax.grid(alpha=.2)
            for ax in axs.flat[len(VARS):]:ax.axis('off')
            axs.flat[0].legend(fontsize=8)
            fig.suptitle(f'{pk} / {tk}: fitted signal widths (not raw data widths)')
            fig.savefig(out/f'fitted_widths_{pk}_{tk}.png',dpi=150);plt.close(fig)
            fig,axs=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
            for ax,v in zip(axs.flat,VARS):
                rr=sorted((r for r in allrows if r['period']==pk and r['topology']==tk and r['variable']==v.branch),key=lambda r:r['t_min'])
                x=[.5*(r['t_min']+r['t_max']) for r in rr]
                a=np.array([r['fitted_data_dvcs_width'] for r in rr]);b=np.array([r['mc_dvcs_width'] for r in rr])
                ax.plot(x,np.divide(a,b,out=np.full_like(a,np.nan),where=b>0),'o-');ax.axhline(1,color='gray',ls='--')
                ax.set_title(v.branch);ax.set_xlabel(r'$|t|$ (GeV$^2$)');ax.grid(alpha=.2)
            for ax in axs.flat[len(VARS):]:ax.axis('off')
            fig.suptitle(f'{pk} / {tk}: fitted DVCS data / MC width')
            fig.savefig(out/f'width_ratios_{pk}_{tk}.png',dpi=150);plt.close(fig)
            fig,axs=plt.subplots(3,3,figsize=(15,11),constrained_layout=True)
            for ax,v in zip(axs.flat,VARS):
                rr=sorted((r for r in allrows if r['period']==pk and r['topology']==tk and r['variable']==v.branch),key=lambda r:r['t_min'])
                x=[.5*(r['t_min']+r['t_max']) for r in rr]
                ax.plot(x,[r['fitted_data_single_containment'] for r in rr],'o-',label='Fitted data DVCS')
                ax.plot(x,[r['mc_single_containment'] for r in rr],'o-',label='DVCS MC')
                ax.set_title(v.branch);ax.set_xlabel(r'$|t|$ (GeV$^2$)');ax.set_ylim(0,1.05);ax.grid(alpha=.2)
            for ax in axs.flat[len(VARS):]:ax.axis('off')
            axs.flat[0].legend(fontsize=8)
            fig.suptitle(f'{pk} / {tk}: one-dimensional nominal-cut containment (NOT joint efficiency)')
            fig.savefig(out/f'single_cut_containment_{pk}_{tk}.png',dpi=150);plt.close(fig)
            fig,axs=plt.subplots(2,4,figsize=(17,8),constrained_layout=True)
            for i,ax in enumerate(axs.flat):
                if i>=len(edges)-1:ax.axis('off');continue
                for kind,color in [('data','black'),('mc','tab:blue')]:
                    dt=inputs[(pk,tk,kind)][1][i]
                    if len(dt)>5:
                        q=np.quantile(dt,[.01,.99]);ax.hist(dt,bins=65,range=(max(-2,q[0]),min(2,q[1])),density=True,histtype='step',label=kind,color=color)
                ax.set_title(f'{edges[i]:.2f}–{edges[i+1]:.2f} GeV²');ax.set_xlabel(r'$t_p-t_\gamma$ (GeV$^2$)');ax.grid(alpha=.2)
            axs.flat[0].legend(fontsize=8);fig.suptitle(f'{pk} / {tk}: four-vector momentum-transfer closure')
            fig.savefig(out/f'delta_t_{pk}_{tk}.png',dpi=150);plt.close(fig)
    (out/'README.txt').write_text('Reuses fit_shared_two_templates from the production exclusivity fitter, with the same default nuisance priors and shared pi0 fraction drivers.\n'
        'Each fit uses the production |t| bin edges pooled over xB/Q2/phi and the production global-cut configuration.\n'
        'Widths are central 68% quantile half-widths of the MORPHED fitted DVCS component, not fitted Gaussian sigma parameters.\n'
        'Containment is ONE-DIMENSIONAL only, evaluated on the nominal data cut window. This is NOT joint signal efficiency and NOT a cross-section correction.\n'
        'No iterative multi-variable conditional fitting or data event-level signal weighting is performed.\n'
        'The pi0 MC is included in every fit; fit failures are not silently replaced by raw widths.\n'
        'Delta t uses measured e, p and gamma four-vectors with period-dependent beam energies.\n'
        'No caches, no PDF, no production modifications.\n')
    print('[output]',out)
if __name__=='__main__':main()
