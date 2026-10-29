#!/usr/bin/env python3
"""Read-only CLAS12 DVCS Pass-2 rebuild inventory and legacy comparison.
Run from any directory; requires uproot, numpy, matplotlib.
"""
import argparse, csv, json, re
from pathlib import Path
import numpy as np

PERIODS = ('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb')
LEGACY_BASE = Path('/work/clas12/thayward/CLAS12_exclusive')
KIN = {'Q2':('Q2','Q2','q2'), 'xB':('x','xB','xb'), 't':('t1','t','minus_t'), 'phi':('phi','phi_trento'), 'W':('W','w')}

def legacy_path(period, mode):
    p=period
    energy='10200' if p=='sp19_inb' else '10594' if p.startswith('sp18') else '10604'
    current='45' if p=='sp18_out' else '50'
    if mode=='dvcs_data':return LEGACY_BASE/f'dvcs/data/pass2/data/dvcs/rga_{p}_epgamma.root'
    if mode=='pi0_data':return LEGACY_BASE/f'eppi0/data/pass2/data/rga_{p}_eppi0.root'
    if mode=='dvcs_mc':return LEGACY_BASE/f'dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_{p}_{current}nA_{energy}MeV.root'
    if mode=='pi0_mc':return LEGACY_BASE/f'eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_{p}_{current}nA_{energy}MeV.root'
    if mode=='pi0_as_dvcs':return LEGACY_BASE/f'eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_{p}_{current}nA_{energy}MeV_epgamma.root'
    raise ValueError(mode)

def tree_info(path,preferred=None):
    import uproot
    with uproot.open(path) as f:
        if preferred and preferred in f:t=f[preferred]
        else:
            trees=[f[k] for k in f.keys(cycle=False) if getattr(f[k],'classname','')=='TTree']
            if not trees:raise ValueError('no TTree found')
            t=trees[0]
        return int(t.num_entries),list(t.keys()),t.name

def write_csv(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def compare_hist(files,legacy,branch,legacy_branch,tree_name,legacy_tree,out,limit):
    import uproot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def read(paths,tname,bname):
        arrays=[];remaining=limit
        for p in paths:
            if remaining<=0:break
            with uproot.open(p) as f:
                t=f[tname]; a=t[bname].array(entry_stop=remaining,library='np')
            a=np.asarray(a).ravel();a=a[np.isfinite(a)];arrays.append(a);remaining-=len(a)
        return np.concatenate(arrays) if arrays else np.array([])
    a=read(files,tree_name,branch);b=read([legacy],legacy_tree,legacy_branch)
    if len(a)<2 or len(b)<2:return None
    lo,hi=np.quantile(np.concatenate([a,b]),[0.005,0.995]);
    if not hi>lo:return None
    ha,edges=np.histogram(a,bins=60,range=(lo,hi));hb,_=np.histogram(b,bins=edges)
    if not ha.sum() or not hb.sum():return None
    pa=ha/ha.sum();pb=hb/hb.sum();tv=0.5*np.abs(pa-pb).sum()
    fig,ax=plt.subplots(figsize=(7,4.4));centers=(edges[:-1]+edges[1:])/2
    ax.step(centers,pa,where='mid',color='#0072B2',label=f'Rebuild (sample {len(a):,})')
    ax.step(centers,pb,where='mid',color='#D55E00',label=f'Legacy (sample {len(b):,})')
    ax.set_xlabel(branch);ax.set_ylabel('Fraction per bin');ax.legend(fontsize=9)
    fig.tight_layout();fig.savefig(out,dpi=150);plt.close(fig)
    return float(tv)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    ap.add_argument('--out',default='/volatile/clas12/thayward/dvcs_pass2_validation_v3')
    ap.add_argument('--sample',type=int,default=250000,help='Maximum entries per side per distribution')
    ap.add_argument('--no-plots',action='store_true')
    args=ap.parse_args();root=Path(args.root);out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    modes={'dvcs_data':('data','DataEvents'),'pi0_data':('data','DataEvents'),'dvcs_mc':('mc','MCEvents'),'pi0_mc':('mc','MCEvents'),'pi0_as_dvcs':('leakage','DataEvents')}
    inventory=[];shapes=[]
    for period in PERIODS:
      for mode,(section,preferred) in modes.items():
        base=root/section/period/mode
        # MC has current / nominal subdirectories; others have parts directly.
        paths=sorted(p for p in base.rglob('*.root') if p.parent.name=='parts' and not p.name.endswith('.partial.root')) if base.exists() else []
        validated=[p for p in paths if Path(str(p)+'.ok').exists()]
        unfinished=len(paths)-len(validated)
        n_new=0;bad=[];new_branches=set();tree_name=preferred
        for p in validated:
            try:
                n,br,tn=tree_info(p,preferred);n_new+=n;new_branches.update(br);tree_name=tn
            except Exception as exc:bad.append(f'{p}: {exc}')
        lp=legacy_path(period,mode);n_old=None;old_branches=set();legacy_tree=None;legacy_error=''
        if lp.exists():
            try:n_old,br,legacy_tree=tree_info(lp,'PhysicsEvents');old_branches=set(br)
            except Exception as exc:legacy_error=str(exc)
        else:legacy_error='legacy file missing'
        row=dict(period=period,mode=mode,rebuilt_files=len(paths),validated_files=len(validated),unvalidated_files=unfinished,unreadable_files=len(bad),rebuilt_entries=n_new,legacy_entries=n_old if n_old is not None else '',entry_ratio=(n_new/n_old if n_old else ''),legacy_path=str(lp),legacy_status=legacy_error or 'ok',rebuild_status='; '.join(bad[:2]) or 'ok')
        inventory.append(row)
        print(f'{period:9} {mode:12} valid={len(validated):3} new={n_new:12,} legacy={str(n_old):>12} ratio={str(round(n_new/n_old,4)) if n_old else "n/a"}',flush=True)
        if args.no_plots or not validated or n_old is None or bad:continue
        # Only compare genuinely matching branch names. Do not guess positional data columns.
        for variable,candidates in KIN.items():
            shared=next((b for b in candidates if b in new_branches and b in old_branches),None)
            if shared is None:continue
            # In MCEvents generated and reconstructed branches can have different meanings;
            # use only exact shared names, and mark this as a diagnostic, not event-matched closure.
            target=out/f'{period}_{mode}_{variable}.png'
            try:
                tv=compare_hist(validated,lp,shared,shared,tree_name,legacy_tree,target,args.sample)
                if tv is not None:shapes.append(dict(period=period,mode=mode,variable=variable,branch=shared,total_variation=tv,plot=str(target)))
            except Exception as exc:print(f'[PLOT_SKIPPED] {period} {mode} {shared}: {exc}',flush=True)
    fields=['period','mode','rebuilt_files','validated_files','unvalidated_files','unreadable_files','rebuilt_entries','legacy_entries','entry_ratio','legacy_path','legacy_status','rebuild_status']
    write_csv(out/'inventory.csv',inventory,fields)
    write_csv(out/'shape_comparisons.csv',shapes,['period','mode','variable','branch','total_variation','plot'])
    (out/'summary.json').write_text(json.dumps({'inventory':inventory,'shape_comparisons':shapes,'warnings':['MCEvents contains all generated HIPO events, while legacy reconstructed PhysicsEvents is selected-event only: raw entry ratios are NOT reconstruction efficiency.','New DataEvents colXX columns are not automatically mapped to legacy branches.','Partial two-node campaigns cannot be compared as complete periods.','MC current variants must be matched before interpreting yield ratios.']},indent=2))
    print('[OUTPUT]',out/'inventory.csv',out/'shape_comparisons.csv',out/'summary.json',flush=True)
if __name__=='__main__':main()
