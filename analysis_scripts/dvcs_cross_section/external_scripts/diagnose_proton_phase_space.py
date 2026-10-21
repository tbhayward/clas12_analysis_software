#!/usr/bin/env python3
"""Pass-2 proton phase-space diagnostics. Run from dvcs_cross_section root.
Requires numpy, matplotlib, uproot. Five period jobs, max eight workers.
Generated MC is plotted at generator level; no event-level truth matching is assumed.
Reconstructed preselection approximates global preselection; optional fitted exclusivity JSON applies topology-specific data/MC windows.
All yield ratios are explicitly shape-normalized, not absolute efficiencies.
"""
import argparse
import concurrent.futures as cf
import csv
import json
import time
import os
from pathlib import Path
import re
import numpy as np
os.environ.setdefault('MPLBACKEND','Agg')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot

PERIODS=['Sp18_inb','Sp18_out','Fa18_inb','Fa18_out','Sp19_inb']
BASE=Path('/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2')
ENERGY={'Sp18_inb':'10594','Sp18_out':'10594','Fa18_inb':'10604','Fa18_out':'10604','Sp19_inb':'10200'}
CURRENT={'Sp18_inb':'50','Sp18_out':'45','Fa18_inb':'50','Fa18_out':'50','Sp19_inb':'50'}
PBINS=np.linspace(0.15,2.5,95); TBINS=np.linspace(30,70,81)
CUTS=[(.30,64.23),(.35,62),(.40,60),(.45,60),(.50,58),(.50,50)]
T_EDGES=[0,.15,.2,.3,.4,.6,1.0]
P_REGIONS=[(.30,.40),(.40,.50),(.50,.70),(.70,1.5)]
TH_REGIONS=[(37,50),(50,58),(58,64.23)]

def path_for(period,kind):
    if kind=='data': return BASE/'data/dvcs'/f'rga_{period.lower()}_epgamma.root'
    stem='gen' if kind=='gen' else 'rec'
    return BASE/'mc/dvcsgen'/f'{stem}_dvcsgen_rga_{period.lower()}_{CURRENT[period]}nA_{ENERGY[period]}MeV.root'

def angle(a):
    a=np.asarray(a,dtype=float)
    return np.degrees(a) if np.nanpercentile(np.abs(a),95)<3.5 else a

def status(period,message):
    print(f"[{time.strftime('%H:%M:%S')}] [{period}] {message}",flush=True)

def load_exclusivity(path):
    if not path:
        return None
    with open(path) as f:
        payload=json.load(f)
    if not isinstance(payload,dict):
        raise ValueError("Exclusivity JSON must be a dictionary")
    return payload

def windows_for(payload,period,detector2,kind):
    if payload is None:
        return {}
    key=f"DVCS_{period.split('_')[0].capitalize()}_{period.split('_')[1].capitalize()}_CD_{'FD' if detector2==1 else 'FT'}"
    # Canonical period labels: Fa18_Inb, Sp19_Inb, etc.
    if key not in payload:
        raise KeyError(f"Missing exclusivity key {key}; example keys: {list(payload)[:5]}")
    block=payload[key]
    if kind not in block:
        raise KeyError(f"Missing '{kind}' cuts in {key}")
    return block[kind]

def collect(period,kind,step,exclusivity_path,strict):
    file=path_for(period,kind)
    if not file.exists(): raise FileNotFoundError(str(file))
    status(period,f"Opening {kind}: {file}")
    start=time.monotonic()
    payload=load_exclusivity(exclusivity_path) if kind!='gen' else None
    with uproot.open(file) as f:
        tree=f['PhysicsEvents']; keys=set(tree.keys())
        required=['p1_p','p1_theta']
        if any(x not in keys for x in required): raise RuntimeError(f'{file}: missing {set(required)-keys}')
        optional=['detector1','detector2','t1','open_angle_ep2','pTmiss','e_theta','p1_phi','e_phi','p2_phi']
        cuts_by_top={}
        if payload is not None:
            for det2 in (0,1):
                cuts_by_top[det2]=windows_for(payload,period,det2,'data' if kind=='data' else 'mc')
            required_cut_branches=set().union(*(v.keys() for v in cuts_by_top.values()))
            missing=required_cut_branches-keys
            if missing: raise RuntimeError(f'{file}: exclusivity branches absent: {sorted(missing)}')
            if 'detector2' not in keys: raise RuntimeError(f'{file}: missing detector2 for topology cuts')
        branches=list(dict.fromkeys(required+[x for x in optional if x in keys]+[x for x in set().union(*(v.keys() for v in cuts_by_top.values())) if x in keys]))
        h2=np.zeros((len(PBINS)-1,len(TBINS)-1)); hlow=np.zeros_like(h2); hmid=np.zeros_like(h2); hhigh=np.zeros_like(h2)
        ht=np.zeros((len(T_EDGES)-1,len(PBINS)-1,len(TBINS)-1))
        counts=np.zeros(len(CUTS)); total=0; lowtotal=0
        status(period,f"{kind}: {tree.num_entries:,} ROOT entries; {len(branches)} branches; exclusivity={'ON' if payload else 'OFF'}")
        scanned=0
        for chunk_index,a in enumerate(tree.iterate(branches,library='np',step_size=step),1):
            scanned+=len(a['p1_p'])
            if chunk_index==1 or chunk_index%10==0:
                status(period,f"{kind}: scanned {scanned:,}/{tree.num_entries:,} ({100*scanned/max(1,tree.num_entries):.1f}%), {time.monotonic()-start:.0f}s")
            p=np.asarray(a['p1_p'],float); th=angle(a['p1_theta']); mask=np.isfinite(p)&np.isfinite(th)
            if kind!='gen':
                if 'detector1' in a: mask &= a['detector1']==2
                if 't1' in a: mask &= (-a['t1']<1.)
                if 'open_angle_ep2' in a: mask &= angle(a['open_angle_ep2'])>5
                if payload is not None:
                    exclusive=np.zeros(len(p),dtype=bool)
                    for det2,cuts in cuts_by_top.items():
                        m=np.asarray(a['detector2'])==det2
                        for branch,definition in cuts.items():
                            low=float(definition['cut_low']); high=float(definition['cut_high'])
                            v=np.asarray(a[branch],float)
                            m &= np.isfinite(v)&(v>=low)&(v<=high)
                        exclusive |= m
                    mask &= exclusive
            else: mask &= (th>=30)&(th<70)
            p=p[mask]; th=th[mask]
            t=np.asarray(a['t1'],float)[mask] if 't1' in a else None
            h2+=np.histogram2d(p,th,bins=(PBINS,TBINS))[0]
            if t is not None:
                abs_t=np.abs(t)
                regions=[abs_t<.2,(abs_t>=.2)&(abs_t<.4),abs_t>=.4]
                for h,sel in zip((hlow,hmid,hhigh),regions): h+=np.histogram2d(p[sel],th[sel],bins=(PBINS,TBINS))[0]
                ht+=np.histogramdd((abs_t,p,th),bins=(T_EDGES,PBINS,TBINS))[0]
                lowtotal+=int(np.sum(regions[0]))
            total+=len(p)
            for i,(pc,tc) in enumerate(CUTS): counts[i]+=np.count_nonzero((p>pc)&(th<tc))
        status(period,f"{kind}: DONE, selected {total:,}, elapsed {time.monotonic()-start:.1f}s")
        return dict(kind=kind,period=period,h=h2,low=hlow,mid=hmid,high=hhigh,ht=ht,counts=counts,total=total,lowtotal=lowtotal)

def plot1d(results,period,out):
    fig,axs=plt.subplots(2,2,figsize=(13,9)); colors={'data':'black','rec':'tab:blue','gen':'tab:orange'}
    for kind,r in results.items():
        h=r['h']; color=colors[kind]
        for ax,v,bins in ((axs[0,0],h.sum(axis=1),PBINS),(axs[0,1],h.sum(axis=0),TBINS)):
            # Unit-area shape comparison; no claim of absolute data/MC normalization.
            y=v/max(v.sum(),1); ax.stairs(y,bins,label=kind,color=color,linewidth=1.5)
        for ax,key in ((axs[1,0],'low'),(axs[1,1],'high')):
            v=r[key].sum(axis=0); ax.stairs(v/max(v.sum(),1),TBINS,label=kind,color=color,linewidth=1.5)
    for ax in axs.flat: ax.grid(alpha=.2); ax.legend(fontsize=8)
    axs[0,0].set(xlabel='Proton momentum (GeV)',ylabel='Unit-area fraction',title='Momentum: all selected events')
    axs[0,1].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='Angle: all selected events')
    axs[1,0].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='Low |t| (<0.2 GeV²)')
    axs[1,1].set(xlabel='Proton angle (degrees)',ylabel='Unit-area fraction',title='High |t| (≥0.4 GeV²)')
    for ax in axs[:,1]:
        for x in (50,58,64.23): ax.axvline(x,color='gray',linestyle='--',alpha=.5)
    axs[0,0].axvline(.5,color='gray',linestyle='--')
    fig.suptitle(f'{period}: generated MC / reconstructed MC / reconstructed data (shape only)')
    fig.tight_layout();fig.savefig(out/'distributions.png',dpi=155);plt.close(fig)

def maps(results,period,out):
    fig,axs=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    for j,kind in enumerate(('gen','rec','data')):
        h=results[kind]['h']; im=axs[0,j].pcolormesh(PBINS,TBINS,np.log10(h.T+1),shading='auto',cmap='viridis');fig.colorbar(im,ax=axs[0,j]);axs[0,j].set_title(f'{kind}: log10(count+1)')
    rec=results['rec']['h'];dat=results['data']['h'];gen=results['gen']['h']
    # Only relative shape ratios. Not detector efficiency or absolute data/MC ratio.
    for j,(num,den,title) in enumerate(((rec,gen,'rec/gen shape ratio (NOT efficiency)'),(dat,rec,'data/recMC shape ratio'),(results['data']['low'],results['rec']['low'],'low-|t| data/recMC shape ratio'))):
        norm=(den.sum()/num.sum()) if num.sum()>0 else np.nan
        val=np.divide(num*norm,den,out=np.full_like(num,np.nan),where=den>=15)
        im=axs[1,j].pcolormesh(PBINS,TBINS,val.T,shading='auto',vmin=0,vmax=2,cmap='coolwarm');fig.colorbar(im,ax=axs[1,j]);axs[1,j].set_title(title+' (unit-area normalized)')
    for ax in axs.flat:
        ax.axvline(.5,color='white',ls='--',lw=1);ax.axhline(58,color='white',ls='--',lw=1)
        ax.set(xlabel='Proton momentum (GeV)',ylabel='Proton angle (degrees)')
    fig.suptitle(period);fig.savefig(out/'maps.png',dpi=145);plt.close(fig)

def ratio_map(data,mc,refdata,refmc,min_mc=15):
    # Reference-region normalized D/MC: unity is the shape baseline in the high-|t| reference.
    scale=refmc.sum()/refdata.sum() if refdata.sum()>0 else np.nan
    return np.divide(data*scale,mc,out=np.full(data.shape,np.nan),where=mc>=min_mc)

def plot_t_slices(results,period,out):
    d=results['data']['ht'];m=results['rec']['ht']
    reference_d=d[-2:].sum(axis=0);reference_m=m[-2:].sum(axis=0)
    fig,axs=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    for i,ax in enumerate(axs.flat):
        ratio=ratio_map(d[i],m[i],reference_d,reference_m)
        im=ax.pcolormesh(PBINS,TBINS,ratio.T,shading='auto',cmap='coolwarm',vmin=0,vmax=2)
        ax.axvline(.5,color='white',ls='--');ax.axhline(58,color='white',ls='--')
        ax.set(title=f'{T_EDGES[i]:.2f} < |t| < {T_EDGES[i+1]:.2f} (GeV²)',xlabel='Proton momentum (GeV)',ylabel='Proton angle (degrees)')
        fig.colorbar(im,ax=ax)
    fig.suptitle(f'{period}: exclusive data/recMC ratios; common high-|t| normalization' if d.sum() else period)
    fig.savefig(out/'data_mc_ratio_by_t.png',dpi=150);plt.close(fig)
    # Momentum-sliced theta and angle-sliced momentum; ratios use SAME high-|t| reference scale.
    for mode,regions,axis in [('theta_by_momentum',P_REGIONS,0),('momentum_by_angle',TH_REGIONS,1)]:
        fig,axs=plt.subplots(len(regions),2,figsize=(13,3.2*len(regions)),layout='constrained')
        for i,(lo,hi) in enumerate(regions):
            edges=PBINS if axis==0 else TBINS
            select=(edges[:-1]>=lo)&(edges[1:]<=hi)
            for j,idx in enumerate((0,len(T_EDGES)-2)):
                hd=d[idx];hm=m[idx]
                dd=hd[select,:].sum(axis=0) if axis==0 else hd[:,select].sum(axis=1)
                mm=hm[select,:].sum(axis=0) if axis==0 else hm[:,select].sum(axis=1)
                rd=reference_d[select,:].sum() if axis==0 else reference_d[:,select].sum()
                rm=reference_m[select,:].sum() if axis==0 else reference_m[:,select].sum()
                ratio=np.divide(dd*(rm/rd if rd>0 else np.nan),mm,out=np.full(dd.shape,np.nan),where=mm>=10)
                bins=TBINS if axis==0 else PBINS
                axs[i,j].stairs(ratio,bins,color='black');axs[i,j].axhline(1,color='gray',ls='--')
                axs[i,j].set(xlabel='Proton angle (degrees)' if axis==0 else 'Proton momentum (GeV)',ylabel='Data/recMC (reference normalized)',title=f'{lo:g}–{hi:g} '+('GeV' if axis==0 else 'degrees')+('  low |t|' if j==0 else '  high |t|'))
                axs[i,j].grid(alpha=.2)
        fig.suptitle(period+' — '+mode)
        fig.savefig(out/(mode+'.png'),dpi=145);plt.close(fig)
    # CSV grid for inspection, including counts and normalized ratio.
    with (out/'data_mc_ratio_by_t_grid.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['period','t_low_GeV2','t_high_GeV2','p_low_GeV','p_high_GeV','theta_low_deg','theta_high_deg','data','rec_mc','ratio_reference_normalized'])
        for i in range(len(T_EDGES)-1):
            r=ratio_map(d[i],m[i],reference_d,reference_m)
            for ip in range(len(PBINS)-1):
                for it in range(len(TBINS)-1):
                    if d[i,ip,it]+m[i,ip,it]==0:continue
                    w.writerow([period,T_EDGES[i],T_EDGES[i+1],PBINS[ip],PBINS[ip+1],TBINS[it],TBINS[it+1],int(d[i,ip,it]),int(m[i,ip,it]),r[ip,it]])

def worker(args):
    period,outdir,step,exclusivity_path,strict=args; out=Path(outdir)/period;out.mkdir(parents=True,exist_ok=True)
    status(period,'START')
    results={k:collect(period,k,step,exclusivity_path,strict) for k in ('data','rec','gen')}
    status(period,'Making distribution plots, |t|-sliced maps, 1D comparisons and CSV')
    plot1d(results,period,out);maps(results,period,out);plot_t_slices(results,period,out)
    rows=[]
    for kind,r in results.items():
        for (pc,tc),n in zip(CUTS,r['counts']):
            rows.append(dict(period=period,sample=kind,p_min_GeV=pc,theta_max_deg=tc,n_survive=int(n),n_baseline=int(r['total']),fraction=float(n/max(r['total'],1))))
    with (out/'survival.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    status(period,'DONE: plots and CSV saved')
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',default='output/proton_phase_space_root_diagnostics')
    ap.add_argument('--workers',type=int,default=5)
    ap.add_argument('--step-size',default='100 MB')
    ap.add_argument('--periods',nargs='+',choices=PERIODS,default=PERIODS)
    ap.add_argument('--exclusivity-json',default=None,help='Path to iterative_cuts/jsons/combined_cuts_95.json; data and MC receive their respective fitted cuts')
    ap.add_argument('--strict',action='store_true',help='Reserved for future exact-global-cut validation; currently rejects execution rather than silently claiming exact equivalence')
    a=ap.parse_args()
    if a.strict:ap.error('Exact full C++ global selection has not been ported. --strict is not supported.')
    if a.exclusivity_json:
        payload=load_exclusivity(a.exclusivity_json)
        for period in a.periods:
            for det2 in (0,1):
                for kind in ('data','mc'):windows_for(payload,period,det2,kind)
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    jobs=[(p,str(out),a.step_size,a.exclusivity_json,a.strict) for p in a.periods];rows=[]
    status('main',f'Launching {len(jobs)} periods with {min(8,max(1,a.workers),len(jobs))} workers, step={a.step_size}, exclusivity={a.exclusivity_json or "OFF"}')
    with cf.ProcessPoolExecutor(max_workers=min(8,max(1,a.workers),len(jobs))) as pool:
        futures={pool.submit(worker,job):job[0] for job in jobs}
        for future in cf.as_completed(futures):
            period=futures[future]
            try:rows.extend(future.result())
            except Exception as exc:
                status('main',f'ERROR in {period}: {exc}')
                raise
    with (out/'survival_all_periods.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    status('main',f'Completed: {out}')
if __name__=='__main__':main()
