#!/usr/bin/env python3
"""Published RGA Pass-1 versus Pass-2 BSA, in Pass-1 (Q2,xB,-t') cells.

Run from dvcs_cross_section/external_scripts:
  python3 compare_pass1_published_bsa.py

Pass-1 and Pass-2 phi samples are NEVER rebinned. Pass-1/Pass-2 reuse events:
independent-error pulls are descriptive, not calibrated significance tests.
"""
from __future__ import annotations
import argparse, ast, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
QEDGES=[0,1.4,1.8,2.4,3.25,5,np.inf]
XEDGES=[[0,.13,.21,np.inf],[0,.13,.21,np.inf],[0,.16,.26,np.inf],
        [0,.21,.33,np.inf],[0,.33,np.inf],[0,.55,np.inf]]
TEDGES=[0,.2,.4,.8,np.inf]
OFFSETS=[0,12,24,36,48,56]
BLUE='#0072B2'; ORANGE='#E69F00'
M=0.9382720813

def tprime(q,x,t):
    eps=4*M*M*x*x/q
    minus_tmin=q*(2*(1-x)*(1-np.sqrt(1+eps))+eps)/(4*x*(1-x)+eps)
    return np.abs(t)-minus_tmin
#enddef

def cell_id(q,x,tp):
    if not all(np.isfinite(v) for v in (q,x,tp)) or tp<0:return -1
    iq=np.searchsorted(QEDGES,q,side='right')-1
    if iq<0 or iq>=6:return -1
    ix=np.searchsorted(XEDGES[iq],x,side='right')-1
    it=np.searchsorted(TEDGES,tp,side='right')-1
    if ix<0 or ix>=len(XEDGES[iq])-1 or it<0 or it>=4:return -1
    return OFFSETS[iq]+(len(XEDGES[iq])-1)*it+ix+1
#enddef

def cell_info(cid):
    z=cid-1
    iq=max(i for i,o in enumerate(OFFSETS) if o<=z)
    nx=len(XEDGES[iq])-1
    it=(z-OFFSETS[iq])//nx
    ix=(z-OFFSETS[iq])%nx
    return iq,ix,it
#enddef

def parse_pair(v):
    try:
        z=ast.literal_eval(str(v))
        # The analysis exports (value, statistical_error, auxiliary_flag).
        # Only the first two entries are used for this comparison.
        if not isinstance(z,(tuple,list)) or len(z)<2:return np.nan,np.nan
        return float(z[0]),float(z[1])
    except (ValueError,TypeError,SyntaxError,IndexError):return np.nan,np.nan
    #endtry
#enddef

def read_pass1(path):
    d=pd.read_csv(path,sep=r'\s+',comment='#',header=None,
                  names=['phi_rad','Q2','xB','t','energy','A','err'])
    d['phi']=np.rad2deg(d.phi_rad)%360
    d['tp']=tprime(d.Q2,d.xB,d.t)
    d['cell']=[cell_id(q,x,t) for q,x,t in zip(d.Q2,d.xB,d.tp)]
    d['source']='Pass-1'
    return d[np.isfinite(d.A)&np.isfinite(d.err)&(d.err>0)&(d.cell>0)].copy()
#enddef

def read_pass2(path):
    src=pd.read_csv(path,low_memory=False)
    out=[]
    for energy,label in [(10.6,'10.6 GeV'),(10.2,'Sp19 Inb')]:
        bsa_col=f'BSA, counts, {label}'
        required=['Bin Name','phimin','phimax','Q2min','Q2max','xBmin','xBmax',
                  't_abs_min','t_abs_max',bsa_col]
        missing=[c for c in required if c not in src.columns]
        if missing:raise ValueError(f'Missing Pass-2 CSV columns: {missing}')
        #endif
        for _,r in src.iterrows():
            a,e=parse_pair(r[bsa_col])
            if not np.isfinite(a) or not np.isfinite(e) or e<=0:continue
            #endif
            cid=pd.to_numeric(r['Bin Name'],errors='coerce')
            low=pd.to_numeric(r['phimin'],errors='coerce')
            high=pd.to_numeric(r['phimax'],errors='coerce')
            if not all(np.isfinite(z) for z in (cid,low,high)):continue
            #endif
            cid=int(cid)
            if cid<1 or cid>64:continue
            #endif
            # In this dedicated CSV t_abs_avg stores event-weighted -t',
            # matching its t_abs_min/max boundaries. Use real mean phi.
            def mean_value(prefix):
                col=f'{prefix}, {label}'
                return pd.to_numeric(r[col],errors='coerce') if col in src.columns else np.nan
            phiavg=mean_value('phiavg')
            qavg=mean_value('Q2avg')
            xavg=mean_value('xBavg')
            tpavg=mean_value('t_abs_avg')
            phi=float(phiavg) if np.isfinite(phiavg) else (float(low)+float(high))/2.0
            complete=all(np.isfinite(z) for z in (phiavg,qavg,xavg,tpavg))
            out.append(dict(phi=phi%360,Q2=qavg,xB=xavg,t=np.nan,tp=tpavg,
                            energy=energy,A=float(a),err=float(e),cell=cid,
                            source='Pass-2',mean_kinematics_cell=cid if complete else -1,
                            phi_min=float(low),phi_max=float(high),mean_available=complete))
        #endfor
    #endfor
    return pd.DataFrame(out,columns=['phi','Q2','xB','t','tp','energy','A','err',
                                     'cell','source','mean_kinematics_cell',
                                     'phi_min','phi_max','mean_available'])
#enddef

def coverage(g,maxgap,minpts):
    ph=np.sort(np.unique(g.phi.to_numpy(float)%360))
    if len(ph)<2:return False,360.,0,0
    #endif
    gap=float(np.max(np.diff(np.r_[ph,ph[0]+360])))
    quadrants=len(set(np.floor(ph/90).astype(int)))
    halves=min(int(np.sum(ph<180)),int(np.sum(ph>=180)))
    return bool(len(ph)>=minpts and gap<=maxgap and quadrants==4 and halves>=3),gap,quadrants,halves
#enddef

def sinonly(ph,a):return a*np.sin(np.deg2rad(ph))
#enddef

def twopar(ph,a,b):
    p=np.deg2rad(ph)
    return a*np.sin(p)/(1+b*np.cos(p))
#enddef

def fit(g,model):
    ph=g.phi.to_numpy(float); y=g.A.to_numpy(float);err=g.err.to_numpy(float)
    npar=1 if model=='sin' else 2
    if len(ph)<=npar+1:return None
    #endif
    try:
        if model=='sin':
            popt,pcov=curve_fit(sinonly,ph,y,sigma=err,absolute_sigma=True,p0=[0.1],maxfev=15000)
        else:
            popt,pcov=curve_fit(twopar,ph,y,sigma=err,absolute_sigma=True,
                                p0=[0.1,0],bounds=([-1.5,-.95],[1.5,.95]),maxfev=15000)
        #endif
        se=np.sqrt(np.diag(pcov))
        if not np.all(np.isfinite(popt)) or not np.all(np.isfinite(se)) or np.any(se<=0):return None
        #endif
        if model=='two' and (abs(popt[1])>.90 or abs(popt[0])>1.45):return None
        #endif
        pred=sinonly(ph,*popt) if model=='sin' else twopar(ph,*popt)
        chi=float(np.sum(((y-pred)/err)**2))
        corr=float(pcov[0,1]/(se[0]*se[1])) if model=='two' else np.nan
        return dict(a=float(popt[0]),sa=float(se[0]),b=float(popt[1]) if model=='two' else np.nan,
                    sb=float(se[1]) if model=='two' else np.nan,chi2=chi,ndf=len(ph)-npar,corr=corr)
    except (RuntimeError,ValueError,TypeError,FloatingPointError):return None
    #endtry
#enddef

def create_plots(p1,p2,fit_df,pairs,out,energy,no_fits):
    base=out/('10p6_GeV' if energy==10.6 else '10p2_GeV')
    canvas=base/'xb_canvases';canvas.mkdir(parents=True,exist_ok=True)
    # Remove legacy xB-grouped plots so old canvases are not mistaken for new output.
    for stale in canvas.glob('bsa_xB_*.png'):stale.unlink()
    aa=p1[np.isclose(p1.energy,energy)];bb=p2[np.isclose(p2.energy,energy)]
    common=sorted(set(aa.cell)&set(bb.cell))
    # Six canvases per energy, one per published Q2 region. The xB
    # subdivisions vary with Q2 and must not be grouped across Q2 regions.
    for iq in range(6):
        cids=[c for c in common if cell_info(c)[0]==iq]
        if not cids:continue
        nx=len(XEDGES[iq])-1
        fig,axes=plt.subplots(4,nx,figsize=(4.5*nx+1,12.0),
                              squeeze=False,sharex=True,sharey=True)
        for ax in axes.flat:
            ax.set_xlim(0,360);ax.set_ylim(-.8,.8)
            ax.axhline(0,color='.65',lw=.75)
            ax.set_xticks([0,90,180,270,360]);ax.grid(alpha=.14)
        for cid in cids:
            _,ix,it=cell_info(cid);ax=axes[it,ix]
            for data,color,marker,name in [(aa,BLUE,'s','Published Pass-1'),
                                            (bb,ORANGE,'o','Pass-2')]:
                g=data[data.cell==cid].sort_values('phi')
                if g.empty:continue
                ax.errorbar(g.phi,g.A,yerr=g.err,fmt=marker,ms=3.5,
                            capsize=1.5,lw=.9,color=color,alpha=.9)
                if not no_fits:
                    row=fit_df[(fit_df.energy==energy)&(fit_df.cell==cid)&
                               (fit_df.source==name.replace('Published ',''))]
                    if len(row):
                        r=row.iloc[0]
                        if r['good'] and np.isfinite(r['a_two']):
                            x=np.linspace(0,360,361)
                            ax.plot(x,twopar(x,r['a_two'],r['b_two']),color=color,lw=1.35)
                #endif
            ax.text(.035,.95,f'Cell {cid}',transform=ax.transAxes,
                    va='top',fontsize=9)
        for ix in range(nx):
            lo,hi=XEDGES[iq][ix:ix+2]
            txt=rf'$x_B\in[{lo:g},{hi:g})$' if np.isfinite(hi) else rf'$x_B\geq {lo:g}$'
            axes[0,ix].set_title(txt,fontsize=11,pad=9)
            axes[-1,ix].set_xlabel(r'$\phi$ (deg)')
        for it in range(4):
            lo,hi=TEDGES[it:it+2]
            lab=f'{lo:g}–{hi:g}' if np.isfinite(hi) else f'>{lo:g}'
            axes[it,0].set_ylabel(r'$A_{LU}$'+'\n'+rf'$-t^\prime$: {lab} GeV$^2$',fontsize=10)
        qlo,qhi=QEDGES[iq:iq+2]
        qtxt=f'{qlo:g}–{qhi:g}' if np.isfinite(qhi) else f'>{qlo:g}'
        fig.suptitle(rf'CLAS12 DVCS BSA: Pass-1 vs Pass-2 | {energy:g} GeV | '
                     +rf'$Q^2$: {qtxt} GeV$^2$',fontsize=15,y=.986)
        from matplotlib.lines import Line2D
        handles=[Line2D([0],[0],color=BLUE,marker='s',ls='-',lw=1.3,ms=5),
                 Line2D([0],[0],color=ORANGE,marker='o',ls='-',lw=1.3,ms=5)]
        fig.legend(handles,['Published Pass-1','Pass-2'],loc='upper center',
                   bbox_to_anchor=(.5,.956),ncol=2,frameon=False,fontsize=11)
        fig.subplots_adjust(left=.105,right=.985,bottom=.075,top=.885,
                            hspace=.22,wspace=.12)
        fig.savefig(canvas/f'bsa_Q2_{iq+1:02d}.png',dpi=170)
        plt.close(fig)
    #endfor
    sub=pairs[np.isclose(pairs.energy,energy)] if len(pairs) else pd.DataFrame()
    if len(sub):
        fig,ax=plt.subplots(figsize=(6.5,5.7))
        ax.errorbar(sub.a_pass1,sub.a_pass2,xerr=sub.sa_pass1,yerr=sub.sa_pass2,
                    fmt='o',color=BLUE,ms=4,alpha=.8,capsize=1)
        lo=min(sub.a_pass1.min(),sub.a_pass2.min())-.05
        hi=max(sub.a_pass1.max(),sub.a_pass2.max())+.05
        ax.plot([lo,hi],[lo,hi],'--',color='.4');ax.set(xlim=(lo,hi),ylim=(lo,hi),
            xlabel='Published Pass-1 sine amplitude',ylabel='Pass-2 sine amplitude',
            title=f'{energy:g} GeV | qualified modulation fits')
        ax.grid(alpha=.2);fig.tight_layout();fig.savefig(base/'amplitude_comparison.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(max(8,len(sub)*.12),4.5))
        ax.errorbar(np.arange(len(sub)),sub.pull,fmt='o',color=BLUE,ms=3)
        for y in [-2,0,2]:ax.axhline(y,color='.5',ls='--' if y else '-',lw=.8)
        #endfor
        ax.set(xlabel='Qualified matched kinematic cell (ordered)',ylabel='Naive statistical pull',
               title=f'{energy:g} GeV | Pass-2 minus Pass-1')
        fig.tight_layout();fig.savefig(base/'amplitude_pulls.png',dpi=180);plt.close(fig)
        fig,ax=plt.subplots(figsize=(6,4.5));ax.hist(sub.pull,bins=np.linspace(-4,4,25),edgecolor='white',color=BLUE)
        ax.set(xlabel='Naive statistical pull',ylabel='Number of cells',title=f'{energy:g} GeV | amplitude pulls')
        fig.tight_layout();fig.savefig(base/'pull_distribution.png',dpi=180);plt.close(fig)
    #endif
    return len(common),len({cell_info(c)[0] for c in common})
#enddef

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pass1',type=Path,default=HERE/'import/rga_pass1_bsa.txt')
    ap.add_argument('--pass2',type=Path,default=HERE.parent/'output/csvs/dvcs_pass2_pass1_bsa_binning.csv')
    ap.add_argument('--output-dir',type=Path,default=HERE.parent/'output/pass1_bsa_comparison')
    ap.add_argument('--min-phi-points',type=int,default=8)
    ap.add_argument('--max-phi-gap',type=float,default=90.)
    ap.add_argument('--fit-model',choices=['both','sin','two'],default='both')
    ap.add_argument('--no-fits',action='store_true')
    ap.add_argument('--verbose',action='store_true')
    args=ap.parse_args();out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    p1=read_pass1(args.pass1);p2=read_pass2(args.pass2)
    print(f'[inputs] Pass-1 {len(p1)} valid points; Pass-2 {len(p2)} valid points')
    if len(p2)==0:raise RuntimeError('No valid Pass-2 BSA values; check the alternate CSV was populated.')
    #endif
    mismatches=p2[(p2.mean_kinematics_cell>0)&(p2.cell!=p2.mean_kinematics_cell)]
    if nmeans==0:
        print('[warning] No event-weighted means found. Rebuild the C++ suite and rerun '
              'the --pass1-bsa-binning workflow before publication plots.')
    nmeans=int(p2.mean_available.sum())
    print(f'[bin audit] Pass-2 mean kinematics available for {nmeans}/{len(p2)} valid points. '
          'Using mean phi when present, nominal phi center only when absent.')
    allpoints=pd.concat([p1,p2],ignore_index=True)
    allpoints.to_csv(out/'csv_matched_source_points.csv',index=False)
    records=[]
    for (energy,cid,source),g in allpoints.groupby(['energy','cell','source']):
        ok,gap,quadrants,halves=coverage(g,args.max_phi_gap,args.min_phi_points)
        s=fit(g,'sin') if not args.no_fits else None
        t=fit(g,'two') if not args.no_fits else None
        good=bool(ok and t is not None and s is not None and
                  (not np.isfinite(t['corr']) or abs(t['corr'])<.98) and
                  abs(t['a']-s['a'])<=2*math.hypot(t['sa'],s['sa']))
        records.append(dict(energy=energy,cell=cid,source=source,n=len(g),max_phi_gap=gap,
                            quadrants=quadrants,min_half_points=halves,coverage_ok=ok,
                            sin_fit_ok=s is not None,two_fit_ok=t is not None,good=good,
                            a_sin=s['a'] if s else np.nan,sa_sin=s['sa'] if s else np.nan,
                            a_two=t['a'] if t else np.nan,sa_two=t['sa'] if t else np.nan,
                            b_two=t['b'] if t else np.nan,sb_two=t['sb'] if t else np.nan,
                            corr=t['corr'] if t else np.nan,chi2_two=t['chi2'] if t else np.nan,
                            ndf_two=t['ndf'] if t else np.nan))
    #endfor
    fits=pd.DataFrame(records);fits.to_csv(out/'modulation_fits.csv',index=False)
    pairs=[]
    for (energy,cid),g in fits.groupby(['energy','cell']):
        if set(g.source)!={'Pass-1','Pass-2'}:continue
        #endif
        r=g.set_index('source');a=r.loc['Pass-1'];b=r.loc['Pass-2']
        if not (a.good and b.good):continue
        #endif
        use_two=args.fit_model!='sin'
        av='a_two' if use_two else 'a_sin';sv='sa_two' if use_two else 'sa_sin'
        delta=b[av]-a[av];sig=math.hypot(a[sv],b[sv])
        pairs.append(dict(energy=energy,cell=cid,model='two' if use_two else 'sin',
                          a_pass1=a[av],sa_pass1=a[sv],a_pass2=b[av],sa_pass2=b[sv],
                          delta=delta,naive_sigma=sig,pull=delta/sig))
    #endfor
    pairs=pd.DataFrame(pairs,columns=['energy','cell','model','a_pass1','sa_pass1','a_pass2','sa_pass2','delta','naive_sigma','pull'])
    pairs.to_csv(out/'modulation_comparison.csv',index=False)
    # Report both functional forms, independently of which is primary.
    model_rows=[]
    for (energy,cid),g in fits.groupby(['energy','cell']):
        if set(g.source)!={'Pass-1','Pass-2'}:continue
        z=g.set_index('source'); a=z.loc['Pass-1'];b=z.loc['Pass-2']
        if not (a.good and b.good):continue
        for model in ('sin','two'):
            key='a_'+model;errkey='sa_'+model
            delta=b[key]-a[key];sig=math.hypot(a[errkey],b[errkey])
            model_rows.append(dict(energy=energy,cell=cid,model=model,
                                   a_pass1=a[key],a_pass2=b[key],
                                   sa_pass1=a[errkey],sa_pass2=b[errkey],
                                   delta=delta,pull=delta/sig,
                                   ratio=b[key]/a[key] if abs(a[key])>1e-8 else np.nan))
    pd.DataFrame(model_rows).to_csv(out/'modulation_comparison_both_models.csv',index=False)
    # Compare actual event-weighted Pass-2 kinematics with the published means.
    kinematics=[]
    for (energy,cid),g in p2.groupby(['energy','cell']):
        if not g.mean_available.any():continue
        v=g[g.mean_available];old=p1[(p1.energy==energy)&(p1.cell==cid)]
        if old.empty:continue
        kinematics.append(dict(energy=energy,cell=cid,
            Q2_pass1=float(old.Q2.mean()),Q2_pass2=float(v.Q2.mean()),
            xB_pass1=float(old.xB.mean()),xB_pass2=float(v.xB.mean()),
            tprime_pass1=float(old.tp.mean()),tprime_pass2=float(v.tp.mean()),
            delta_Q2=float(v.Q2.mean()-old.Q2.mean()),
            delta_xB=float(v.xB.mean()-old.xB.mean()),
            delta_tprime=float(v.tp.mean()-old.tp.mean())))
    pd.DataFrame(kinematics).to_csv(out/'mean_kinematics_comparison.csv',index=False)

    summary=[f'Pass-2 populated mean-kinematics rows: {nmeans}/{len(p2)}',
             'Pass-1 vs Pass-2 BSA; same underlying events imply correlated statistical errors.',
             'Naive independent-error pulls are diagnostic only.',
             f'Pass-2 mean-kinematics cell mismatches: {len(mismatches)}']
    for energy in [10.6,10.2]:
        ncell,ncanvas=create_plots(p1,p2,fits,pairs,out,energy,args.no_fits)
        sub=pairs[np.isclose(pairs.energy,energy)]
        chisq=float(np.sum(sub.pull**2)) if len(sub) else np.nan
        line=(f'{energy:g} GeV: {ncell} matched 3D cells, {ncanvas} Q2 canvases, '
              f'{len(sub)} qualified amplitude pairs, naive chi2/N={chisq/len(sub) if len(sub) else np.nan:.3f}')
        print('[summary]',line);summary.append(line)
    #endfor
    (out/'comparison_summary.txt').write_text('\n'.join(summary)+'\n')
    print('[output]',out)
#enddef

if __name__=='__main__':main()
