#!/usr/bin/env python3
"""Direct experimental BSA comparison: published Pass-1, Sebastian Pass-2, our Pass-2.

No KM15/CFF model. Sebastian Appendix-H points were digitized from a PDF; their
error bars are graphical estimates. Sebastian kinematics are Pass-1 proxies.
Same-run Pass-2 datasets are statistically correlated; independent-error pull
numbers are diagnostics, NOT calibrated compatibility probabilities.
"""
from __future__ import annotations
import argparse, ast, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
M=0.9382720813

def tmin(Q,x):
    e=4*M*M*x*x/Q
    return -Q*(2*(1-x)*(1-np.sqrt(1+e))+e)/(4*x*(1-x)+e)
#enddef

def cell_assign(d):
    regions=[(0,1.4,[0,.13,.21,np.inf]),(1.4,1.8,[0,.13,.21,np.inf]),
             (1.8,2.4,[0,.16,.26,np.inf]),(2.4,3.25,[0,.21,.33,np.inf]),
             (3.25,5,[0,.33,np.inf]),(5,np.inf,[0,.55,np.inf])]
    edges=[0,.2,.4,.8,np.inf]
    out=d.copy(); out['cell']=-1
    tp=tmin(out.Q2.to_numpy(float),out.xB.to_numpy(float))-out.t.to_numpy(float)
    out['tprime']=tp
    k=0
    for q0,q1,xedges in regions:
        for x0,x1 in zip(xedges[:-1],xedges[1:]):
            for t0,t1 in zip(edges[:-1],edges[1:]):
                k+=1
                sel=(out.Q2>=q0)&(out.Q2<q1)&(out.xB>=x0)&(out.xB<x1)&(tp>=t0)&(tp<t1)
                out.loc[sel,'cell']=k
            #endfor
        #endfor
    #endfor
    return out
#enddef

def parse_pair(x):
    try:
        a=ast.literal_eval(str(x));return float(a[0]),float(a[1])
    except (ValueError,TypeError,SyntaxError,IndexError):return np.nan,np.nan
    #endtry
#enddef

def read_our(path,labels):
    src=pd.read_csv(path,low_memory=False); result=[]
    for label in labels:
        b=f'BSA, counts, {label}'
        cols=[b,f'xBavg, {label}',f'Q2avg, {label}',f't_abs_avg, {label}',f'phiavg, {label}']
        missing=[c for c in cols if c not in src]
        if missing:raise ValueError(f'Missing {label}: {missing}')
        #endif
        for _,row in src.iterrows():
            a,stat=parse_pair(row[b]);x=row[cols[1]];q=row[cols[2]];t=row[cols[3]];phi=row[cols[4]]
            if not all(np.isfinite(v) for v in (a,stat,x,q,t,phi)) or stat<=0:continue
            #endif
            result.append(dict(source='Our Pass-2',subset=label,Q2=q,xB=x,t=-abs(t),phi_deg=phi,A=a,err=stat,error_kind='statistical'))
        #endfor
    #endfor
    return cell_assign(pd.DataFrame(result))
#enddef

def read_pass1(path):
    d=pd.read_csv(path,sep=r'\s+',comment='#',header=None,names=['phi','Q2','xB','t','E','A','err'])
    d=d[np.isclose(d.E,10.6)].copy();d['phi_deg']=np.rad2deg(d.phi)
    d['source']='Published Pass-1';d['subset']='Fa18 published combined';d['error_kind']='published uncertainty'
    return cell_assign(d)
#enddef

def read_seb(path,proxy):
    d=pd.read_csv(path)
    d=d[(d.series=='sebastian_pass2')&(d.topology=='egamma_p')].copy()
    k=pd.read_csv(proxy)
    d=d.merge(k[['published_rga_bin','Q2','xB','t']],on='published_rga_bin',validate='many_to_one')
    d['sebastian_published_bin']=d.published_rga_bin.astype(int)
    d=cell_assign(d)  # Map published Table-10 numbering to analysis cells via physical kinematics.
    d['source']='Sebastian Pass-2';d['subset']=d.torus.map({'inbending':'Fa18 Inb','outbending':'Fa18 Out'})
    d['A']=d.bsa;d['err']=d.bsa_error_graphical;d['error_kind']='digitized graphical uncertainty'
    return d
#enddef

def amplitude(d):
    """Weighted one-harmonic sin(phi) fit, no imposed intercept; diagnostic only."""
    d=d[np.isfinite(d.A)&np.isfinite(d.err)&(d.err>0)&np.isfinite(d.phi_deg)]
    if len(d)<4:return None
    #endif
    x=np.sin(np.deg2rad(d.phi_deg.to_numpy(float)));y=d.A.to_numpy(float);e=d.err.to_numpy(float)
    w=1/e**2;den=np.sum(w*x*x)
    if den<=0:return None
    #endif
    a=np.sum(w*x*y)/den;se=np.sqrt(1/den)
    chi=float(np.sum(((y-a*x)/e)**2))
    return dict(amplitude=float(a),error=float(se),n=int(len(d)),chi2=chi,ndf=len(d)-1,
                Q2=float(d.Q2.mean()),xB=float(d.xB.mean()),tabs=float((-d.t).mean()))
#enddef

def run(args):
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    p1=read_pass1(args.pass1);our=read_our(args.our,['Fa18 Inb','Fa18 Out','Fa18']);seb=read_seb(args.sebastian,args.proxy)
    print(f'[input] Pass-1 {len(p1)}; our Pass-2 {len(our)}; Sebastian Pass-2 {len(seb)}')
    if not set(seb.cell.unique()).issubset(set(p1.cell.unique())):
        print('[warning] Some Sebastian cells have no published Pass-1 points')
    #endif
    # Published Table-10 IDs and the fit script's sequential cell IDs differ.
    # Match by Q2/xB/t-prime region, never by raw integer ID.
    groups=[]
    for name,df in [('Published Pass-1',p1),('Sebastian Pass-2',seb),('Our Pass-2',our)]:
        for (sub,cell),g in df.groupby(['subset','cell']):
            if cell<1:continue
            #endif
            r=amplitude(g)
            if r:groups.append(dict(source=name,subset=sub,cell=int(cell),**r))
            #endif
        #endfor
    #endfor
    fits=pd.DataFrame(groups);fits.to_csv(out/'sinphi_amplitudes.csv',index=False)
    pairs=[]
    for subset in ['Fa18 Inb','Fa18 Out']:
        s=fits[(fits.source=='Sebastian Pass-2')&(fits.subset==subset)].set_index('cell')
        o=fits[(fits.source=='Our Pass-2')&(fits.subset==subset)].set_index('cell')
        p=fits[fits.source=='Published Pass-1'].set_index('cell')
        for cell in sorted(set(s.index)&set(o.index)&set(p.index)):
            sr=s.loc[cell];orr=o.loc[cell];pr=p.loc[cell]
            if isinstance(sr,pd.DataFrame) or isinstance(orr,pd.DataFrame) or isinstance(pr,pd.DataFrame):continue
            #endif
            pairs.append(dict(subset=subset,cell=cell,Q2_proxy=sr.Q2,xB_proxy=sr.xB,t_abs_proxy=sr.tabs,
                Q2_our=orr.Q2,xB_our=orr.xB,t_abs_our=orr.tabs,
                sebastian_A=sr.amplitude,sebastian_err=sr.error,our_A=orr.amplitude,our_err=orr.error,
                pass1_A=pr.amplitude,pass1_err=pr.error,
                our_minus_sebastian=orr.amplitude-sr.amplitude,
                our_minus_pass1=orr.amplitude-pr.amplitude,
                sebastian_minus_pass1=sr.amplitude-pr.amplitude,
                independent_error_pull_our_sebastian=(orr.amplitude-sr.amplitude)/math.hypot(orr.error,sr.error),
                note='Same Pass-2 events may overlap: independent-error pull is not a significance'))
        #endfor
    #endfor
    pair=pd.DataFrame(pairs);pair.to_csv(out/'three_way_amplitude_comparison.csv',index=False)
    for subset in ['Fa18 Inb','Fa18 Out']:
        sel=pair[pair.subset==subset]
        print(f'[three-way {subset}] {len(sel)} common cells')
        if len(sel):
            print(f'  median(our - Sebastian) sinphi amplitude = {sel.our_minus_sebastian.median():+.4f}')
            print(f'  median(our - Pass1) sinphi amplitude = {sel.our_minus_pass1.median():+.4f}')
            print(f'  median(Sebastian - Pass1) sinphi amplitude = {sel.sebastian_minus_pass1.median():+.4f}')
        #endif
        directory=out/subset.replace(' ','_');directory.mkdir(exist_ok=True)
        for cell in sel.cell:
            fig,ax=plt.subplots(figsize=(8.2,5.4))
            sources=[('Published Pass-1',p1,'Fa18 published combined','o'),('Sebastian Pass-2',seb,subset,'s'),('Our Pass-2',our,subset,'^')]
            for name,df,sub,marker in sources:
                g=df[(df.cell==cell)&(df.subset==sub)].sort_values('phi_deg')
                if g.empty:continue
                #endif
                ax.errorbar(g.phi_deg,g.A,yerr=g.err,fmt=marker,ms=3.3,capsize=1.5,alpha=.75,label=f'{name} ({len(g)} pts)')
                row=fits[(fits.source==name)&(fits.subset==sub)&(fits.cell==cell)]
                if len(row):
                    ph=np.linspace(0,360,361);ax.plot(ph,float(row.iloc[0].amplitude)*np.sin(np.deg2rad(ph)),lw=1.2,alpha=.8)
                #endif
            #endfor
            rr=sel[sel.cell==cell].iloc[0]
            ax.set_title(f'Fa18 {subset[-3:]} | $\\langle Q^2\\rangle_{{proxy}}$={rr.Q2_proxy:.2f} GeV$^2$, '
                         f'$\\langle x_B\\rangle_{{proxy}}$={rr.xB_proxy:.3f}, '
                         f'$\\langle -t\\rangle_{{proxy}}$={rr.t_abs_proxy:.3f} GeV$^2$')
            ax.set_xlabel(r'$\phi$ (degrees)');ax.set_ylabel(r'$A_{LU}$');ax.set_xlim(0,360)
            ax.axhline(0,color='gray',lw=.7);ax.legend(fontsize=8);ax.grid(alpha=.15)
            fig.text(.5,.01,'Sebastian: digitized PDF points; kinematics are published Pass-1 proxies. Overlapping Pass-2 samples are correlated.',ha='center',fontsize=7)
            fig.tight_layout(rect=(0,.035,1,1));fig.savefig(directory/f'cell_{int(cell):02d}.png',dpi=145);plt.close(fig)
        #endfor
        if len(sel):
            fig,ax=plt.subplots(figsize=(8.8,5.2))
            ss=sel.sort_values(['Q2_proxy','xB_proxy','t_abs_proxy'])
            xx=np.arange(len(ss));ax.errorbar(xx,ss.our_minus_sebastian,yerr=np.hypot(ss.our_err,ss.sebastian_err),fmt='o',capsize=2)
            ax.axhline(0,color='gray',lw=.8);ax.set_xticks(xx,[f'{r.Q2_proxy:.1f}/{r.xB_proxy:.2f}/{r.t_abs_proxy:.2f}' for r in ss.itertuples()],rotation=90,fontsize=7)
            ax.set_ylabel(r'Our $\sin\phi$ amplitude $-$ Sebastian amplitude')
            ax.set_xlabel(r'Proxy $\langle Q^2\rangle$ / $\langle x_B\rangle$ / $\langle -t\rangle$ (GeV$^2$ where applicable)')
            ax.set_title(f'{subset}: direct fitted BSA amplitude differences (error bars ignore event covariance)')
            fig.tight_layout();fig.savefig(out/f'amplitude_difference_{subset.replace(" ","_")}.png',dpi=170);plt.close(fig)
        #endif
    #endfor
    (out/'README.txt').write_text('Direct BSA comparison, no KM15. Fits A_LU(phi)=a*sin(phi) to each experimental dataset independently.\n'
       'Published Pass-1 data are Fa18 combined, NOT separated by torus polarity.\n'
       'Sebastian Pass-2 Appendix H: egamma_p only; digitized graphical errors, not source numerical tables.\n'
       'Sebastian mean Q2/xB/t are published Pass-1 proxies, NOT independently measured Sebastian means.\n'
       'Our Pass-2 is native-bin data assigned to published 64-bin regions, not truly rebinned event yields.\n'
       'Sin(phi)-only fit is descriptive and can differ from the original extraction harmonic model.\n'
       'Our and Sebastian Pass-2 likely share events: no covariance, so difference pulls are not significance.\n'
       'Per-bin kinematic variations and phi coverage can also create apparent differences.\n')
    print('[write]',out)
#enddef

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--our',type=Path,default=HERE.parent/'output/csvs/dvcs_pass2_analysis.csv')
    ap.add_argument('--pass1',type=Path,default=HERE/'import/rga_pass1_bsa.txt')
    ap.add_argument('--sebastian',type=Path,default=HERE/'import/sebastian_pass2_bsa_from_note.csv')
    ap.add_argument('--proxy',type=Path,default=HERE/'import/sebastian_appendixH_kinematics_proxy.csv')
    ap.add_argument('--output',type=Path,default=HERE.parent/'output/bsa/three_way_pass1_sebastian_our')
    run(ap.parse_args())
#endif
