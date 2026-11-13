#!/usr/bin/env python3
"""Order-independent two-photon control shapes from fresh DATA and existing MC caches."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,csv,html,json,shutil
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts
from study_particle_common_cut_stability import export_digest
from study_cd_proton_angle_response import SOURCES,replace_once

GROUPS=[p+'_'+g for p in ('FDp','CDp') for g in ('FD_FD','FD_FT','FT_FT')]

def isolated_sources(source,target):
    target.mkdir()
    for name in SOURCES:shutil.copy2(source/name,target/name)
    p=target/'particle_response_cut_arrays.h';s=p.read_text()
    s=replace_once(s,'const int ncols=photons?28:extended?20:12;','const int ncols=29;')
    s=replace_once(s,'const std::array<double,8>* photon_values=nullptr){','const std::array<double,8>* photon_values=nullptr,double pair_mass=NAN){')
    s=replace_once(s,'int d1=int(v[3]),d2=int(v[4]),top=d1==1&&d2==1?0:d1==2&&d2==1?1:d1==2&&d2==0?2:-1;if(top<0)return;', '''int d1=int(v[3]),d2=int(v[4]);
  if(!control||!photons||!photon_values)throw std::runtime_error("Unordered export requires control photon tags");
  if(!std::isfinite((*photon_values)[3])||!std::isfinite((*photon_values)[7]))throw std::runtime_error("Nonfinite photon detector tag");
  int g1=int((*photon_values)[3]),g2=int((*photon_values)[7]);
  if((d1!=1&&d1!=2)||(g1!=0&&g1!=1)||(g2!=0&&g2!=1))throw std::runtime_error("Invalid unordered-control detector tag");
  int top=(d1==2?3:0)+2-g1-g2;''')
    old='bool global=!is_excluded_run(int(v[23]),cfg)&&passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],cfg)&&passes_proton_momentum_cut(d1,period,v[8],cfg)&&passes_global_cuts(v[2],v[14],v[15],d1,d2,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],cfg);'
    new='''auto assigned_global=[&](int tag){return passes_global_cuts(v[2],v[14],v[15],d1,tag,period,v[5],v[6],v[7],v[9],v[10],v[11],v[12],v[13],cfg);};
  // Conservative common domain: both legacy photon-detector assignments pass.
  // Existing control globals still use the composite pi0 coordinates, explicitly.
  bool global=!is_excluded_run(int(v[23]),cfg)&&passes_analysis_phase_space_cuts(period,v[0],v[1],v[2],v[5],cfg)&&passes_proton_momentum_cut(d1,period,v[8],cfg)&&assigned_global(g1)&&assigned_global(g2);'''
    s=replace_once(s,old,new)
    s=replace_once(s,'std::array<double,28> row{};','std::array<double,29> row{};row[28]=pair_mass;')
    s=replace_once(s,'emit(vary_control(c,pars,beam,base,noise,mass),c.base.key,photons?&pp:nullptr);','auto shifted=vary_control(c,pars,beam,base,noise,mass);emit(shifted,c.base.key,photons?&pp:nullptr,mass);')
    s=replace_once(s,'emit(r.v,mix(pathkey(path)^uint64_t(i)),&pp);','emit(r.v,mix(pathkey(path)^uint64_t(i)),&pp,c.photons[8]);')
    p.write_text(s)

def compile_worker(out,reference,work):
    source=work/'unordered_sources';isolated_sources(base.SOURCE,source);old=base.SOURCE
    try:
        base.SOURCE=source
        return base.compile_worker(out,reference)
    finally:base.SOURCE=old

def blocks(spec,period,group):
    cfg=next(c for c in cuts.production.PERIODS if c.key==period)
    # One common intersection per proton region, independent of photon pair/order.
    tops=[0] if group<3 else [1,2]
    return [spec['exclusivity'][cuts.production.cut_json_key('eppi0',cfg,cuts.production.TOPOLOGIES[t])]['data'] for t in tops]

def selection(a,group,windows,omit=None):
    keep=(a[:,0]==group)&(a[:,8]>.5)
    for block in windows:
        for name,w in block.items():
            if name==omit:continue
            col=1+base.OBS.index(name);v=a[:,col];keep&=np.isfinite(v)&(v<=float(w['cut_high']))
            if w.get('mode')!='upper_quantile':keep&=v>=float(w['cut_low'])
    return keep

def histogram(a,keep,col,edges,half=None):
    if half is not None:keep=keep&(a[:,11]==half)
    values=a[keep,col];counts=np.zeros(len(edges)+1,dtype=np.int64);finite=np.isfinite(values);np.add.at(counts,np.searchsorted(edges,values[finite],side='right'),1);return counts

def render(period,arrays,spec,out,min_entries,dpi):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    variables=[*spec['variables'],dict(name='m_gamma_gamma',label='Two-photon mass (GeV)',edges=np.linspace(.10,.17,141).tolist(),bins=140)]
    rows=[];flows=[];links=[];metrics=[]
    for scope in ('global_only','common_windows_Nminus1'):
        for group,label in enumerate(GROUPS):
            windows=[] if scope=='global_only' else blocks(spec,period,group)
            fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained')
            for scenario,a in arrays.items():
                keep=selection(a,group,windows);flows.append(dict(period=period,scope=scope,group=label,sample=scenario,selected=int(keep.sum()),unresolved_photons=int((keep&~np.isfinite(a[:,20])).sum())))
            for vi,v in enumerate(variables):
                name=v['name'];edges=base.diag.variable_edges(v);col=28 if name=='m_gamma_gamma' else 1+base.OBS.index(name);ax=axes.flat[vi]
                d=histogram(arrays['DATA'],selection(arrays['DATA'],group,windows,name),col,edges);curves=[('DATA',d,'#205cab')];models={s:[] for s in ('nominal','frozen')}
                for scenario in models:
                    a=arrays[scenario];keep=selection(a,group,windows,name)
                    for half in (0,1):models[scenario].append(histogram(a,keep,col,edges,half))
                for scenario,color in [('nominal','#888888'),('frozen','#22946b')]:
                    for half in (0,1):
                        m=models[scenario][half].astype(float);ref=models['nominal'][half].astype(float);score,bins=base.shape_score(d,m,m,int(ref.sum()),min_entries,ref,ref);dm=base.diag.shape_metrics(d,edges,0,0);mm=base.diag.shape_metrics(m,edges,0,0);width=mm['core_spread']
                        rows.append(dict(period=period,scope=scope,group=label,variable=name,scenario=scenario,mc_half=half,data_count=int(d.sum()),raw_mc=int(m.sum()),score=score,usable_bins=bins,data_over_mc_width=dm['core_spread']/width if width>0 else float('nan'),median_shift_mc_widths=(dm['median']-mm['median'])/width if width>0 else float('nan'),status='usable' if min(d.sum(),m.sum(),ref.sum())>=min_entries and np.isfinite(score) else 'insufficient support'))
                    curves.append((scenario+' MC',models[scenario][0]+models[scenario][1],color))
                for text,h,color in curves:
                    if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,label=text,color=color)
                if name=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_ in curves]))
                elif name!='m_gamma_gamma':
                    qs=[base.diag.quantile(h[1:-1],edges,q) for _,h,_ in curves if h.sum()>0 for q in (.01,.99)];qs=[q for q in qs if np.isfinite(q)]
                    if qs:
                        lo,hi=min(qs),max(qs);pad=max(.001,(hi-lo)*.08);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                ax.set(xlabel=v['label'],ylabel='Unit-area density',title=name);ax.legend(fontsize=7);ax.grid(alpha=.2)
            for ax in axes.flat[len(variables):]:ax.axis('off')
            fig.suptitle(period+' · '+label+' · '+scope+'\nUnordered photon detectors; fixed response, no refitting');folder=out/'plots'/period/scope;folder.mkdir(parents=True,exist_ok=True);path=folder/(label+'.png');fig.savefig(path,dpi=dpi);plt.close(fig);links.append(path.relative_to(out).as_posix())
    # Physical tag/angle audit before any exclusivity windows: distributions, not inferred tags.
    for scenario,a in arrays.items():
        keep=a[:,8]>.5
        for role,col in [(1,20),(2,24)]:
            for tag in (0,1):
                m=keep&(a[:,col+3]==tag)&np.isfinite(a[:,col+1]);theta=a[m,col+1]*180/np.pi
                metrics.append(dict(period=period,sample=scenario,photon_role=role,detector='FT' if tag==0 else 'FD',events=len(theta),theta_p01=float(np.quantile(theta,.01)) if len(theta) else float('nan'),theta_median=float(np.median(theta)) if len(theta) else float('nan'),theta_p99=float(np.quantile(theta,.99)) if len(theta) else float('nan')))
    return rows,flows,links,metrics

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joint-run',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--workers',type=int,default=10);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--min-entries',type=int,default=100);p.add_argument('--dpi',type=int,default=100);args=p.parse_args(argv)
    if not 1<=args.workers<=10 or args.min_entries<10 or args.dpi<50:p.error('Invalid workers, entries or dpi')
    out=args.output.resolve();run=args.joint_run.resolve()
    if out.exists():p.error('Output exists; choose a fresh directory')
    meta=json.loads((run/'manifest.json').read_text());spec=json.loads((run/'spec.json').read_text());cm=json.loads((run/'cache_manifest.json').read_text());reference=Path(meta['arguments']['reference_study']);ref=json.loads((reference/'manifest.json').read_text())
    if not spec.get('joint_control_fit') or cm['reference_sha256']!=base.digest(reference/'manifest.json') or meta['reference_manifest_sha256']!=base.digest(reference/'manifest.json'):p.error('Joint cache/reference mismatch')
    for name in ('particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h'):
        if base.digest(base.SOURCE/name)!=meta['source_hashes'][name]:p.error('Response evaluator changed: '+name)
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/name)!=spec['input_hashes'][key]:p.error('Saved global source changed')
    jobs=[]
    for period in spec['periods']:
        entries=[r for r in ref['inventory'] if r['period']==period and r['sample']=='pi0_data']
        if not entries:p.error('Missing DATA inventory')
        for r in entries:
            st=Path(r['path']).stat()
            if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):p.error('DATA input changed')
        jobs.append(dict(period=period,sample='pi0_data',files=[r['path'] for r in entries],scenario='DATA',export_individual_photons=True,export_global_only=True))
        cache=cm['samples'][period+'__pi0_mc'];path=(run/cache['file']).resolve()
        if export_digest(path)!=cache['sha256']:p.error('MC cache changed')
        pars=cuts.load_candidate_parameters(run/period/'candidate_parameters.json')
        for scenario,parameters in [('nominal',[0.]*20),('frozen',pars)]:jobs.append(dict(period=period,sample='pi0_mc',cache=str(path),parameters=parameters,scenario=scenario,export_individual_photons=True,export_global_only=True))
    out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));base.dump(out/'spec.json',spec);worker=compile_worker(out,reference,work);pool=base.Pool(worker,out/'spec.json',work,args.workers)
    print('[study] DATA control reread once; both MC hypotheses reuse caches; order-independent global intersection and six proton/pair groups; production unchanged',flush=True)
    records=[];arrays={period:{} for period in spec['periods']}
    for _,job,path in pool.run('cut-arrays',jobs,'unordered_control'):
        shape=json.loads(Path(str(path)+'.json').read_text())['shape']
        if len(shape)!=2 or shape[1]!=29 or path.stat().st_size!=int(np.prod(shape))*8:raise ValueError('Invalid unordered control export')
        arrays[job['period']][job['scenario']]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,29));records.append(dict(period=job['period'],scenario=job['scenario'],path=str(path),shape=shape,sha256=export_digest(path)))
    # Catch changed DATA inputs during streaming.
    for r in ref['inventory']:
        if r['sample']=='pi0_data':
            st=Path(r['path']).stat()
            if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('DATA changed during export')
    base.dump(out/'export_manifest.json',dict(records=records));base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},joint_manifest_sha256=base.digest(run/'manifest.json'),candidates={period:base.digest(run/period/'candidate_parameters.json') for period in spec['periods']},sources={name:base.digest(out/'backend'/name) for name in SOURCES},script_sha256=base.digest(Path(__file__))))
    scores=[];flows=[];links=[];metrics=[]
    for period,a in arrays.items():
        rr,ff,ll,mm=render(period,a,spec,out,args.min_entries,args.dpi);scores+=rr;flows+=ff;links+=ll;metrics+=mm;print('[done] '+period,flush=True)
    for name,rows in [('shape_scores',scores),('selection_counts',flows),('photon_tag_angles',metrics)]:base.table(out/(name+'.csv'),rows)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Unordered control response</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Control shapes by unordered photon detectors</h1><p>FD or CD proton, crossed with photon pairs FD–FD, FD–FT and FT–FT. Photon order never changes the group. DATA, nominal AAOGEN MC and frozen morphed AAOGEN MC; no response refits or absolute-rate normalization. The pi0 mass is now shown explicitly.</p><p>Fresh DATA export and cached MC exports reevaluate global cuts using BOTH possible legacy photon-detector assignments, retaining their intersection. Existing globals still use the stored composite pi0 coordinates; this conservative diagnostic domain is not a new physical two-photon efficiency model. It is permutation invariant and includes all six proton/pair groups. Upstream processing losses remain unrecoverable from these trees.</p><p>Start with global_only plots (no exclusivity windows). common_windows_Nminus1 applies an identical intersection of saved DATA windows within each proton region, omitting the displayed variable. CD uses CD–FD and CD–FT window intersection; FD uses its available FD–FD windows, including the newly displayed FD-proton mixed/FT–FT groups. These are conservative comparison windows, not newly optimized 95% cuts. Mass plots apply the other existing windows; no additional pi0 mass selection is made. Unresolved photon coordinates remain in composite-variable plots and are counted; individual angle summaries exclude them.</p><p>Score CSV uses common nominal-MC support/variance, both MC halves and low-statistics flags; half 1 is a check on MC sampling, with shared DATA. Shape curves use both halves for precision. None of these counts or shape scores determines a production contamination or cross-section correction. Main DVCS acceptance, pi0 subtraction, production cuts and current calibration are unchanged.</p><p><a href="shape_scores.csv">Shape scores</a> · <a href="selection_counts.csv">Counts and unresolved pairs</a> · <a href="photon_tag_angles.csv">Recovered angles by detector/role</a></p>'+''.join('<p><a href="'+x+'">'+html.escape(x)+'</a></p>' for x in links));print('[done] '+str(out/'REPORT.html')+'; arrays stay on scratch',flush=True)

if __name__=='__main__':main()
