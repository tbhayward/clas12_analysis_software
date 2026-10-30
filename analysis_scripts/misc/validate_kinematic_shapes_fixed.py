#!/usr/bin/env python3
"""Compare legacy and rebuilt selected-data kinematics with named or positional branches."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from schema import SCHEMAS
PERIODS=('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb')
VARIABLES=('Q2','x','W','t1','t2','phi1','phi2','e_p','e_theta','e_phi','p1_p','p1_theta','p1_phi','p2_p','p2_theta','p2_phi','Mx2','Mx2_1','Mx2_2','Emiss2','pTmiss','theta_gamma_gamma','Mh_gammagamma')
TOPOLOGIES={'all':None,'FD-FD':(1,1),'CD-FD':(2,1),'CD-FT':(2,0)}

def legacy_path(period,mode,base):
    if mode=='dvcs_data':return base/f'dvcs/data/pass2/data/dvcs/rga_{period}_epgamma.root'
    return base/f'eppi0/data/pass2/data/rga_{period}_eppi0.root'

def resolve(tree,mode,variable):
    names=set(tree.keys())
    if variable in names:return variable
    schema=SCHEMAS[mode]
    if variable in schema:
        c=f'col{schema.index(variable)+1:02d}'
        if c in names:return c
    return None

def collect(paths,tree_name,mode,variables,topology,limit):
    output={v:[] for v in variables};count=0
    for path in paths:
        if count>=limit:break
        with uproot.open(path) as f:
            tree=f[tree_name] if tree_name in f else f["DataEvents" if tree_name=="PhysicsEvents" else "PhysicsEvents"]
            resolved={v:resolve(tree,mode,v) for v in variables}
            if topology is not None:
                resolved['detector1']=resolve(tree,mode,'detector1')
                resolved['detector2']=resolve(tree,mode,'detector2')
            missing=[k for k,v in resolved.items() if v is None]
            if missing:raise RuntimeError(f'{path}: missing branches {missing}')
            unique=list(set(resolved.values()))
            for batch in tree.iterate(unique,step_size='100 MB',library='np'):
                if count>=limit:break
                n=len(batch[unique[0]])
                mask=np.ones(n,dtype=bool)
                if topology is not None:
                    mask=(batch[resolved['detector1']]==topology[0])&(batch[resolved['detector2']]==topology[1])
                for variable in variables:
                    a=np.asarray(batch[resolved[variable]])[mask]
                    a=a[np.isfinite(a)]
                    if len(a):output[variable].append(a[:max(0,limit-count)])
                #endfor
                count+=int(mask.sum())
            #endfor
        #endif
    #endfor
    return {v:np.concatenate(a) if a else np.array([]) for v,a in output.items()}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    ap.add_argument('--legacy-base',default='/work/clas12/thayward/CLAS12_exclusive')
    ap.add_argument('--out',default='/volatile/clas12/thayward/dvcs_pass2_validation_v4')
    ap.add_argument('--sample',type=int,default=200000)
    ap.add_argument('--periods',nargs='+',default=list(PERIODS))
    ap.add_argument('--modes',nargs='+',default=['dvcs_data','pi0_data'])
    a=ap.parse_args();root=Path(a.root);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for period in a.periods:
        for mode in a.modes:
            paths=sorted(p for p in (root/'data'/period/mode/'parts').glob('*.root') if not p.name.endswith('.partial.root') and Path(str(p)+'.ok').exists())
            old=legacy_path(period,mode,Path(a.legacy_base))
            if not paths or not old.exists():print('[SKIP]',period,mode,'missing input');continue
            variables=[v for v in VARIABLES if v in SCHEMAS[mode]]
            for topology,det in TOPOLOGIES.items():
                if mode=='pi0_data' and topology!='all':continue
                try:
                    new=collect(paths,'DataEvents',mode,variables,det,a.sample)
                    olddata=collect([old],'PhysicsEvents',mode,variables,det,a.sample)
                except Exception as exc:
                    print('[ERROR]',period,mode,topology,str(exc),flush=True);continue
                #endtry
                for variable in variables:
                    x=new[variable];y=olddata[variable]
                    if len(x)<5 or len(y)<5:continue
                    bounds=np.quantile(np.concatenate([x,y]),[0.005,0.995]);lo,hi=map(float,bounds)
                    if hi<=lo:continue
                    hnew,edges=np.histogram(x,bins=60,range=(lo,hi));hold,_=np.histogram(y,bins=edges)
                    if not hnew.sum() or not hold.sum():continue
                    pn=hnew/hnew.sum();po=hold/hold.sum();tv=float(np.abs(pn-po).sum()/2)
                    centers=(edges[1:]+edges[:-1])/2
                    fig,(ax,ratio)=plt.subplots(2,1,figsize=(8,6),sharex=True,gridspec_kw={'height_ratios':[3,1]})
                    ax.step(centers,po,where='mid',label=f'Legacy ({len(y):,} sampled)',color='#D55E00')
                    ax.step(centers,pn,where='mid',label=f'Rebuild ({len(x):,} sampled)',color='#0072B2')
                    ax.legend(fontsize=9);ax.set_ylabel('Fraction per bin');ax.set_title(f'{period} / {mode} / {topology} — {variable}')
                    r=np.divide(pn,po,out=np.full_like(pn,np.nan),where=po>0)
                    ratio.step(centers,r,where='mid',color='black');ratio.axhline(1,color='gray',linestyle='--')
                    ratio.set_ylabel('New / old');ratio.set_xlabel(variable)
                    fig.tight_layout()
                    filename=f'{period}_{mode}_{topology}_{variable}.png';fig.savefig(out/filename,dpi=140);plt.close(fig)
                    rows.append(dict(period=period,mode=mode,topology=topology,variable=variable,new_sample=len(x),legacy_sample=len(y),total_variation=tv,plot=str(out/filename)))
                #endfor
                print('[COMPARE]',period,mode,topology,'variables',len(variables),flush=True)
            #endfor
        #endfor
    #endfor
    with (out/'shape_comparisons.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['period','mode','topology','variable','new_sample','legacy_sample','total_variation','plot']);w.writeheader();w.writerows(rows)
    (out/'summary.json').write_text(json.dumps({'plots':len(rows),'comparison':'shape normalized to unit integral, selected data only','warning':'Sampled first N selected events per source; does not establish run-by-run HIPO equivalence.'},indent=2))
    print('[OUTPUT]',len(rows),'plots',out,flush=True)
if __name__=='__main__':main()
