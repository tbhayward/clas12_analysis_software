#!/usr/bin/env python3
"""Fixed add-pT response: FD–FD leakage budgets and kinematic cut stability."""
from __future__ import annotations
import os
for _n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_n]='1'
import argparse,csv,html,json
from pathlib import Path
import numpy as np
import study_particle_common_cut_stability as common
base=common.base;cuts=common.cuts;closure=common.closure
AXES={'xB':(12,1.),'Q2':(13,1.),'electron_p':(14,1.),'electron_theta_deg':(15,180/np.pi),'proton_p':(16,1.),'proton_theta_deg':(17,180/np.pi),'photon_p':(18,1.),'photon_theta_deg':(19,180/np.pi)}


def partitions(data,axis):
    col,factor=AXES[axis];values=data[:,col]*factor;finite=values[np.isfinite(values)]
    if len(finite)<3:return []
    edges=np.unique(np.quantile(finite,[0,1/3,2/3,1]))
    if len(edges)<2:return []
    return [(float(lo),float(hi),i==len(edges)-2) for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:]))]


def subset(arrays,axis,lo,hi,last):
    col,factor=AXES[axis]
    return {s:a[np.isfinite(a[:,col])&(a[:,col]*factor>=lo)&((a[:,col]*factor<=hi) if last else (a[:,col]*factor<hi))] for s,a in arrays.items()}


def shell_rows(arrays,block,thresholds,scale,prob,variable='pTmiss',center=0.):
    """Disjoint annuli, not differences treated as independent cumulative counts."""
    counts=[];col=1+base.OBS.index(variable)
    for sample in ('data','background','mc'):
        a=arrays[sample];sel=closure.mask(a,0,block,'data',omit=variable)
        x=np.sort(np.abs(a[sel,col]-center));counts.append(np.diff(np.searchsorted(x,thresholds,side='right')))
    alpha=scale['background_scale'];weight=alpha/prob['background']
    ctrl_d=scale['control_data_anchor'];ctrl_m=scale['control_mc_anchor']*prob['pi0_mc']
    va=alpha**2*(1/ctrl_d+1/ctrl_m) if ctrl_d>0 and ctrl_m>0 else float('nan');result=[]
    for i,(lo,hi) in enumerate(zip(thresholds[:-1],thresholds[1:])):
        d,b,s=[int(c[i]) for c in counts];bg=weight*b;res=d-bg;err=np.sqrt(d+weight**2*b+(b/prob['background'])**2*va)
        result.append(dict(shell_low=float(lo),shell_high=float(hi),data_count=d,background_mc_count=b,signal_mc_count=s,predicted_background=float(bg),subtracted_data=float(res),conditional_stat_error=float(err),residual_over_stat_error=float(res/err) if err>0 else float('nan'),data_over_predicted_background=float(d/bg) if bg>0 else float('nan'),status='negative residual' if res<0 else 'nonnegative residual'))
    return result


def render_group(period,bounds,axis,groups,block,scale,prob,out):
    import matplotlib.pyplot as plt
    folder=out/'plots'/period/'FD_FD'/f't{bounds[0]:.2f}_{bounds[1]:.2f}';folder.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');scanrows=[];shells=[];summaries=[]
    colors=['#205cab','#22946b','#d45c18']
    for gi,(label,arrays,low,high) in enumerate(groups):
        color=colors[gi%len(colors)]
        for vi,var in enumerate(('pTmiss','Mx2')):
            st=block['data'].get(var);center=0. if var=='pTmiss' else (float(st['cut_low'])+float(st['cut_high']))/2 if st else 0.
            ref=float(st['cut_high'])-center if st else (.2 if var=='pTmiss' else .02)
            th=np.unique(np.r_[np.linspace(.04,.42,25) if var=='pTmiss' else np.linspace(.002,.045,25),ref])
            scan=common.calculate_scan(arrays,0,block,var,bounds,th,ref,center,scale,prob)
            meta=dict(period=period,t_low=bounds[0],t_high=bounds[1],axis=axis,subdivision=label,kinematic_low=low,kinematic_high=high,variable=var)
            for r in scan:r.update(meta)
            scanrows+=scan
            ax=axes[1,vi];ax.errorbar(th,[r['relative_signal_scale'] for r in scan],yerr=[r['conditional_stat_error'] for r in scan],label=label,color=color,lw=1,elinewidth=.5)
            ax.axhline(1,color='black',lw=.7);ax.axvline(ref,color='black',ls=':',lw=.7)
            ax.set(xlabel=var+(' upper cut (GeV)' if vi==0 else ' half-width (GeV²)'),ylabel='Signal scale / reference');ax.grid(alpha=.2)
            target=.15 if vi==0 else .01;near=min(scan,key=lambda r:abs(r['threshold']-target));atref=min(scan,key=lambda r:abs(r['threshold']-ref))
            summaries.append(dict(**meta,test_threshold=near['threshold'],relative_signal_scale=near['relative_signal_scale'],conditional_stat_error=near['conditional_stat_error'],reference_data=atref['data_count'],reference_signal_mc=atref['signal_mc_count'],reference_background_fraction=atref['background_estimate']/atref['data_count'] if atref['data_count'] else float('nan'),support='low statistics' if min(near['data_count'],atref['data_count'])<100 or min(near['signal_mc_count'],atref['signal_mc_count'])<50 else 'adequate counts; model uncertainties excluded'))
            edges=np.linspace(0,.45,121) if vi==0 else np.linspace(-.045,.045,151)
            h=[]
            for s in ('data','mc','background'):
                a=arrays[s];sel=closure.mask(a,0,block,'data',omit=var);h.append(closure.histogram(a[:,1+base.OBS.index(var)],sel,edges))
            d,s,b=h;ss,bb=common.scaled_shape(d,s,b,scale['signal_scale'],scale['background_scale'],prob,True)
            ax=axes[0,vi]
            if d.sum():ax.stairs(d[1:-1]/d.sum()/np.diff(edges),edges,color=color,ls='--',label=label+' DATA')
            ax.stairs((ss+bb)[1:-1]/np.diff(edges),edges,color=color,label=label+' summed MC');ax.set(xlabel=var,ylabel='Equal-area density');ax.grid(alpha=.2)
            if vi==0:
                thresholds=np.unique(np.r_[0,.08,.12,.16,.20,.24,.30,.36,.42,ref])
                for r in shell_rows(arrays,block,thresholds,scale,prob):r.update(meta);shells.append(r)
    for ax in axes.flat:ax.legend(fontsize=7)
    fig.suptitle(f'{period} · FD–FD · {bounds[0]:.2f}<|t|<{bounds[1]:.2f} · {axis}\nFixed add-pT response; common DATA windows; other cuts retained')
    path=folder/(axis+'.png');fig.savefig(path,dpi=110);plt.close(fig)
    return scanrows,shells,summaries,path.relative_to(out).as_posix()


def plot_shell_budget(period,bounds,rows,out):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(11,7),layout='constrained',sharex=True)
    x=np.arange(len(rows));data=np.array([r['data_count'] for r in rows]);background=np.array([r['predicted_background'] for r in rows]);error=np.array([r['conditional_stat_error'] for r in rows])
    axes[0].bar(x-.18,data,width=.36,color='#205cab',label='DATA');axes[0].bar(x+.18,background,width=.36,color='#22946b',label='Predicted π⁰')
    axes[0].errorbar(x-.18,data,yerr=np.sqrt(data),fmt='none',color='black',lw=.6)
    axes[0].errorbar(x+.18,background,yerr=np.sqrt(np.maximum(0,error**2-data)),fmt='none',color='black',lw=.6)
    axes[0].set_ylabel('Events in added interval');axes[0].legend()
    residual=data-background;axes[1].bar(x,residual,color=['#bd4242' if y<0 else '#205cab' for y in residual]);axes[1].errorbar(x,residual,yerr=error,fmt='none',color='black',lw=.8);axes[1].axhline(0,color='black',lw=.8)
    axes[1].set(ylabel='DATA − predicted π⁰',xlabel='Disjoint missing-pT interval (GeV)');axes[1].set_xticks(x,[f"({r['shell_low']:.3g}, {r['shell_high']:.3g}]" for r in rows],rotation=40,ha='right')
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle(f'{period} · FD–FD · {bounds[0]:.2f}<|t|<{bounds[1]:.2f} · background budget\nFixed control-derived π⁰ scale; other common cuts retained')
    path=out/'plots'/period/'FD_FD'/f't{bounds[0]:.2f}_{bounds[1]:.2f}'/'shell_budget.png';fig.savefig(path,dpi=110);plt.close(fig);return path.relative_to(out).as_posix()


def analyze(exports,spec,normalizations,periods,family,closed,out):
    scans=[];shells=[];summaries=[];links=[];coverage=[]
    for period in periods:
        cfg=next(p for p in cuts.production.PERIODS if p.key==period);key=cuts.production.cut_json_key('dvcs',cfg,cuts.production.TOPOLOGIES[0])
        blocks=json.loads((closed/'global_cut_refits'/period/'add_pt'/'combined_cuts_95.json').read_text()) if family=='refitted_DATA' else spec['exclusivity']
        block=common.common_blocks(blocks)[key];scale,prob=normalizations[period]
        for bounds in common.BINS:
            a={s:v[(v[:,9]>=bounds[0])&(v[:,9]<bounds[1])] for s,v in exports[period].items()}
            coverage.append(dict(period=period,t_low=bounds[0],t_high=bounds[1],global_data=len(a['data']),global_signal_mc=len(a['mc']),global_background_mc=len(a['background'])))
            if not len(a['data']):continue
            groups=[('inclusive',a,float('nan'),float('nan'))]
            for axis in ('inclusive',*AXES):
                if axis!='inclusive':groups=[(f'{lo:.3g}–{hi:.3g}',subset(a,axis,lo,hi,last),lo,hi) for lo,hi,last in partitions(a['data'],axis)]
                if not groups:continue
                rr,ss,tt,link=render_group(period,bounds,axis,groups,block,scale,prob,out);scans+=rr;shells+=ss;summaries+=tt;links.append(link)
                if axis=='inclusive':links.append(plot_shell_budget(period,bounds,ss,out))
        print('[done] '+period+' FD–FD shell budgets and kinematic subdivisions',flush=True)
    for name,rows in [('scans',scans),('shell_budgets',shells),('summary',summaries),('coverage',coverage)]:base.table(out/(name+'.csv'),rows)
    gallery=''
    for p in periods:
        gallery+='<details open><summary>'+html.escape(p)+'</summary>'
        for b in common.BINS:
            stem=f'plots/{p}/FD_FD/t{b[0]:.2f}_{b[1]:.2f}/';ll=[l for l in links if l.startswith(stem)]
            gallery+=f'<p><b>{b[0]:.2f}–{b[1]:.2f} GeV²</b>: '+(' · '.join('<a href="'+l+'">'+html.escape(Path(l).stem)+'</a>' for l in ll) if ll else 'No DATA in this cell')+'</p>'
        gallery+='</details>'
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>FD–FD background transfer</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}summary{cursor:pointer;margin:12px}</style><h1>FD–FD background-transfer and kinematic checks</h1><p>Fixed add-pT response, no fit or extra histogram morph. All samples receive identical '+family+' windows and saved full globals. The old common-cut study supplies fixed two-photon/background and topology signal anchors; no new per-subdivision background normalization is fitted.</p><p>Read shell_budget and inclusive first, then each particle momentum/angle or xB/Q² subdivision. Top: DATA dashed, summed DVCS+π⁰ MC solid, same color per subdivision, independently equal total areas. Bottom: (DATA−predicted π⁰)/DVCS MC relative to the dotted reference threshold, separately in each subdivision. Flatness is the test; crossing one at the reference is guaranteed.</p><p>Subdivisions are DATA terciles before exclusivity in each |t| bin, with the same numerical edges on MC. Angles are degrees. These are one-dimensional subdivisions, not a joint kinematic matching or proof of which particle is responsible. Momentum transfer and all particle coordinates of MC are evaluated after the fixed response. Other exclusivity cuts remain applied when plotting/scanning a variable.</p><p>shell_budget.png shows DATA and predicted π⁰ counts in disjoint intervals above, and DATA−π⁰ below. A negative residual is highlighted in red; no DVCS signal shape or signal scale enters that budget. shell_budgets.csv contains disjoint pTmiss intervals (low,high], DATA, predicted π⁰, and their difference. Negative signal residuals identify intervals where this model predicts more background than total DATA; do not clip them. data_over_predicted_background is a point diagnostic, not a confidence limit or a fitted transfer correction. Errors include conditional counting/control statistics, not response/cut/transfer-model uncertainties. summary.csv flags low counts. No cross sections or absolute acceptances are inferred.</p><p>'+ ' · '.join('<a href="'+n+'.csv">'+n+'</a>' for n in ('summary','shell_budgets','scans','coverage'))+'</p>'+gallery)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--common-study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--periods',nargs='+');p.add_argument('--workers',type=int,default=10);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--reuse-exports',type=Path);p.add_argument('--window-family',choices=['saved_DATA','refitted_DATA'],default='refitted_DATA')
    args=p.parse_args(argv)
    if not 1<=args.workers<=10:p.error('Workers must be 1–10 total')
    origin=args.common_study.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a new directory')
    cm=json.loads((origin/'manifest.json').read_text());source=Path(cm['arguments']['followup_study'])/'cut_redetermination';closed=Path(cm['arguments']['closure_study']);old=json.loads((source/'manifest.json').read_text());spec=json.loads((source/'spec.json').read_text());periods=list(dict.fromkeys(args.periods or spec['periods']))
    if not set(periods)<=set(spec['periods']):p.error('Requested period absent from source')
    if base.digest(source/'manifest.json')!=cm['source_sha256'] or base.digest(closed/'manifest.json')!=cm['closure_sha256'] or old['candidate_hashes']!=cm['candidate_hashes']:p.error('Common-study provenance differs from source/cuts')
    ref=Path(old['arguments']['reference_study']);run=Path(old['arguments']['response_run']);training=Path(old['arguments']['training_study'])
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    cache=json.loads((run/'cache_manifest.json').read_text())
    if base.digest(run/'cache_manifest.json')!=old['response_cache_sha256']:p.error('Response cache metadata changed')
    paths=cuts.candidate_paths({'add_pt':run},periods,training);vectors={}
    for (period,c),path in paths.items():
        if base.digest(path)!=cm['candidate_hashes'][c][period]:p.error('Candidate parameters changed')
        vectors[period]=cuts.load_candidate_parameters(path)
    norms={}
    for r in csv.DictReader((origin/'normalizations.csv').open()):
        if r['scenario']=='add_pt' and r['topology']=='FD_FD':norms[r['period']]={k:float(v) for k,v in r.items() if k not in ('period','scenario','topology','status')}
    if not set(periods)<=norms.keys():p.error('Common-study FD–FD normalization missing')
    fingerprint=dict(common_sha256=base.digest(origin/'manifest.json'),normalizations_sha256=base.digest(origin/'normalizations.csv'),candidate_hashes={period:cm['candidate_hashes']['add_pt'][period] for period in periods},spec_sha256=base.digest(source/'spec.json'),worker_sources={n:base.digest(base.SOURCE/n) for n in ('particle_response_cut_arrays.h','particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h','exclusivity_resolution_core.h','exclusivity_resolution_worker.cpp')},saved_globals={n:base.digest(ref/'backend'/n) for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h')})
    exports={period:{} for period in periods};records=[]
    if args.reuse_exports:
        saved=json.loads((args.reuse_exports/'export_manifest.json').read_text())
        if saved['fingerprint']!=fingerprint:p.error('Reused exports have different inputs/worker')
        records=saved['records']
        for r in records:
            if common.export_digest(r['path'])!=r['sha256']:p.error('Reused export changed')
    else:
        # Validate every source before creating report/scratch or starting a ROOT scan.
        jobs=[];meta=json.loads((ref/'manifest.json').read_text())
        if base.digest(ref/'manifest.json')!=old['reference_sha256']:p.error('Reference inventory changed')
        for period in periods:
            entries=[r for r in meta['inventory'] if r['period']==period and r['sample']=='data']
            if not entries:p.error('Missing DATA inventory')
            for r in entries:
                st=Path(r['path']).stat()
                if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):p.error('DATA input changed')
            common_job=dict(period=period,export_kinematics=True,export_topology=0,export_global_only=True)
            jobs.append(dict(**common_job,sample='data',files=[r['path'] for r in entries]))
            for sample in ('mc','background'):
                r=cache['samples'][period+'__'+sample];path=(run/r['file']).resolve()
                if common.export_digest(path)!=r['sha256']:p.error('MC cache changed')
                jobs.append(dict(**common_job,sample=sample,cache=str(path),parameters=vectors[period]))
        out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));base.dump(out/'spec.json',spec)
        print('[export] FD–FD/full globals only; DATA ROOT once per period, existing signal/leakage MC caches; no control reread or response fit',flush=True)
        worker=base.compile_worker(out,ref);pool=base.Pool(worker,out/'spec.json',work,args.workers)
        for i,job,path in pool.run('cut-arrays',jobs,'fd_kinematics'):
            shape=json.loads(Path(str(path)+'.json').read_text())['shape']
            records.append(dict(period=job['period'],sample=job['sample'],path=str(path),shape=shape,sha256=common.export_digest(path)))
        for r in meta['inventory']:
            if r['period'] in periods and r['sample']=='data':
                st=Path(r['path']).stat()
                if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('DATA input changed during export')
    for r in records:
        shape=r['shape'];path=Path(r['path'])
        if len(shape)!=2 or shape[1]!=20 or path.stat().st_size!=int(np.prod(shape))*8:p.error('Invalid extended export layout')
        exports[r['period']][r['sample']]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,20))
    if any(set(samples)!=set(('data','mc','background')) for samples in exports.values()):p.error('Incomplete export samples')
    out.mkdir(parents=True,exist_ok=True);base.dump(out/'export_manifest.json',dict(fingerprint=fingerprint,records=records));base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},fingerprint=fingerprint,window_family=args.window_family))
    normprob={period:(norms[period],{s:cache['samples'][period+'__'+s]['probability'] for s in ('mc','background','pi0_mc')}) for period in periods}
    analyze(exports,spec,normprob,periods,args.window_family,closed,out)
    print('[done] '+str(out/'REPORT.html')+'; transfer report directory only; extended arrays remain on scratch',flush=True)

if __name__=='__main__':main()
