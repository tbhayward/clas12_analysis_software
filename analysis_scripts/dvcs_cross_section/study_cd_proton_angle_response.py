#!/usr/bin/env python3
"""One extra CD-proton momentum-width dependence, using completed joint caches."""
from __future__ import annotations
import argparse,copy,csv,json,shutil,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import study_particle_response as base
import study_frozen_response_kinematics as frozen
import study_response_cut_redetermination as cuts
from study_particle_common_cut_stability import export_digest

PLOT_LOCK=threading.Lock()
EXTRA='proton_CD_low_angle_p_sigma_fraction'
SOURCES=('exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h','particle_response_core.h','particle_response_worker.cpp','particle_control_core.h','particle_control_worker.h','particle_response_cut_arrays.h')

def replace_once(s,old,new):
    if s.count(old)!=1:raise ValueError('Evaluator changed; cannot safely apply isolated angle extension: '+old)
    return s.replace(old,new)

def angle_sources(source,target):
    """Compile-time isolated extension: production's 20-parameter law is untouched."""
    target.mkdir()
    for name in SOURCES:shutil.copy2(source/name,target/name)
    p=target/'particle_response_core.h';s=p.read_text()
    s=replace_once(s,'NC=24, NP=20','NC=24, NP=21')
    anchor='inline std::array<double,NC> vary('
    helper='''// Additional fractional momentum width: full at theta<=35 deg, zero at >=55 deg.
inline double cd_angle_sigma(double theta,const std::array<double,NP>& pars){
 return pars[20]*std::clamp((55.-theta*180/PI)/20.,0.,1.);
}
'''
    s=replace_once(s,anchor,helper+anchor)
    s=replace_once(s,'if(pars[off]==0&&pars[off+1]==0&&pars[off+2]==0&&pars[off+3]==0)continue;',
      'double extra=regions[i]==2?cd_angle_sigma(v[col+1],pars):0.;\n  if(pars[off]==0&&pars[off+1]==0&&pars[off+2]==0&&pars[off+3]==0&&extra==0)continue;\n  active=active||extra!=0;')
    s=replace_once(s,'pars[off]*z(3*i)','(extra==0?pars[off]:std::hypot(pars[off],extra))*z(3*i)')
    p.write_text(s)
    p=target/'particle_control_core.h';s=p.read_text()
    s=replace_once(s,'int j=4*region;v[col]*=1+pars[j+3]+pars[j]*z[first];',
      'int j=4*region;double extra=region==2?cd_angle_sigma(v[col+1],pars):0.;v[col]*=1+pars[j+3]+(extra==0?pars[j]:std::hypot(pars[j],extra))*z[first];')
    p.write_text(s)
    p=target/'particle_response_worker.cpp';s=p.read_text()
    s=replace_once(s,'if(!affected)continue;','affected=affected||(int(r.v[3])==2&&parameters[s][20]!=0);if(!affected)continue;');p.write_text(s)

def compile_angle_worker(out,reference,work):
    source=work/'angle_sources';angle_sources(base.SOURCE,source)
    original=base.SOURCE
    try:
        base.SOURCE=source
        return base.compile_worker(out,reference)
    finally:base.SOURCE=original

def pick(parameters,scores,start):
    eligible=[i for i,p in enumerate(parameters) if list(p[:20])==list(start)]
    finite=[i for i in eligible if np.isfinite(scores[i]['train'])]
    if not finite:raise ValueError('No finite joint training objective')
    return min(finite,key=lambda i:(scores[i]['train'],parameters[i][20]))

def check_source_provenance(meta,fingerprint):
    # Fitting does not call cut-arrays. Its helper gained the optional 28-column
    # photon export after the joint fit; compare it with the exporting study.
    evaluators=('particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h')
    for name in evaluators:
        if base.digest(base.SOURCE/name)!=meta['source_hashes'][name]:raise ValueError('Frozen evaluator differs: '+name)
    for name in ('particle_response_cut_arrays.h',*evaluators):
        if base.digest(base.SOURCE/name)!=fingerprint['export_sources'][name]:raise ValueError('Frozen diagnostic exporter differs: '+name)


def load_inputs(run,diagnostic):
    meta=json.loads((run/'manifest.json').read_text());spec=json.loads((run/'spec.json').read_text());cm=json.loads((run/'cache_manifest.json').read_text())
    if not spec.get('joint_control_fit') or spec.get('joint_control_weight',.5)!=.5:raise ValueError('Requires default-weight joint fit')
    reference=Path(meta['arguments']['reference_study']);ref=json.loads((reference/'manifest.json').read_text())
    if base.digest(reference/'manifest.json')!=meta['reference_manifest_sha256'] or cm['reference_sha256']!=base.digest(reference/'manifest.json'):raise ValueError('Reference provenance changed')
    em=json.loads((diagnostic/'export_manifest.json').read_text());fp=em['fingerprint']
    check_source_provenance(meta,fp)
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/name)!=spec['input_hashes'][key]:raise ValueError('Saved global selection changed')
    cache={k:dict(r,path=str((run/r['file']).resolve())) for k,r in cm['samples'].items()}
    for r in cache.values():
        if export_digest(r['path'])!=r['sha256']:raise ValueError('Cache checksum changed: '+r['path'])
    starts={p:cuts.load_candidate_parameters(run/p/'candidate_parameters.json') for p in spec['periods']}
    if fp['joint_manifest']!=base.digest(run/'manifest.json') or fp['cache_manifest']!=base.digest(run/'cache_manifest.json') or fp['candidates']!={p:base.digest(run/p/'candidate_parameters.json') for p in spec['periods']}:raise ValueError('Frozen diagnostic uses different fit/cache')
    data_arrays={p:{} for p in spec['periods']}
    for r in em['records']:
        if r['scenario']!='shared':continue
        if export_digest(r['path'])!=r['sha256']:raise ValueError('DATA export checksum changed')
        data_arrays[r['period']][r['sample']]=np.memmap(r['path'],dtype='<f8',mode='r',shape=tuple(r['shape']))
    if any(set(v)!={'data','pi0_data'} for v in data_arrays.values()):raise ValueError('Missing frozen DATA arrays')
    for r in ref['inventory']:
        if r['sample'] in ('data','pi0_data'):
            st=Path(r['path']).stat()
            if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('Reference DATA changed')
    return meta,spec,cm,reference,ref,cache,starts,data_arrays

def low_high(period,engine,worker,pool,parameters,indices,multipliers,data_arrays,out):
    import matplotlib.pyplot as plt
    cfg=next(p for p in cuts.production.PERIODS if p.key==period)
    split=next(r for r in frozen.choose_splits(data_arrays['data'],data_arrays['pi0_data']) if r['category']=='proton_CD' and r['axis']=='angle')
    arrays={};jobs=[]
    for name,pars in parameters.items():
        for sample in ('mc','background','pi0_mc'):jobs.append(dict(period=period,sample=sample,cache=engine.cache[period+'__'+sample]['path'],parameters=pars,scenario=name,export_kinematics=True,export_individual_photons=sample=='pi0_mc',export_global_only=True))
    for _,job,path in pool.run('cut-arrays',jobs,'low_high_arrays'):
        shape=json.loads(Path(str(path)+'.json').read_text())['shape'];arrays[job['scenario'],job['sample']]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape))
    audit=[]
    for label,a in [('DATA control',data_arrays['pi0_data']),('Frozen MC control',arrays['frozen','pi0_mc'])]:
        keep=frozen.category_mask(a,'photon1_FD')&frozen.category_mask(a,'photon2_FD')
        first=int((a[keep,20]>a[keep,24]).sum());second=int((a[keep,24]>a[keep,20]).sum())
        audit.append(dict(period=period,sample=label,both_fd=int(keep.sum()),photon1_harder=first,photon2_harder=second,photon1_harder_fraction=first/int(keep.sum()) if keep.any() else float('nan')))
    base.table(out/period/'photon_order_audit.csv',audit)
    rows=[];links=[]
    for channel in ('epgamma','epgammagamma'):
        for top in (1,2):
            block=engine.spec['exclusivity'][cuts.production.cut_json_key('dvcs' if channel=='epgamma' else 'eppi0',cfg,cuts.production.TOPOLOGIES[top])]
            fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
            for si,side in enumerate(('low','high')):
                for vi,name in enumerate(('pTmiss','Emiss2','Mx2_2')):
                    variable=next(v for v in engine.spec['variables'] if v['name']==name);edges=base.diag.variable_edges(variable)
                    d=frozen.histogram(data_arrays['data' if channel=='epgamma' else 'pi0_data'],top,block,name,split,side,edges);ax=axes[si,vi];curves=[('DATA',d,'#205cab')]
                    for scenario in ('frozen','angle'):
                        for half in (0,1):
                            if channel=='epgammagamma':
                                m=frozen.histogram(arrays[scenario,'pi0_mc'],top,block,name,split,side,edges,half);variance=m.copy();raw=m.sum()
                            else:
                                beta,alpha,_=engine.normalization(period,indices[scenario],half,top,multipliers[scenario]);m=np.zeros_like(d,dtype=float);variance=m.copy();raw=0
                                for sample,scale in [('mc',beta),('background',alpha)]:
                                    h=frozen.histogram(arrays[scenario,sample],top,block,name,split,side,edges,half);weight=2*scale/engine.cache[period+'__'+sample]['probability'];m+=weight*h;variance+=weight*weight*h;raw+=h.sum()
                            # Fixed frozen reference for score denominator; same bins/support across models.
                            if scenario=='frozen':
                                if half==0:refs={}
                                refs[half]=(m.copy(),variance.copy(),raw)
                            refm,refv,refraw=refs[half];score,bins=base.shape_score(d,m,variance,int(refraw),engine.min_mc,refm,refv)
                            dm=base.diag.shape_metrics(d,edges,0,0);mm=base.diag.shape_metrics(m,edges,0,0);w=mm['core_spread'];rows.append(dict(period=period,channel=channel,topology=base.TOPS[top],side=side,angle_split_deg=split['threshold'],variable=name,scenario=scenario,mc_half=half,data_count=d.sum(),raw_mc=raw,shape_score=score,usable_bins=bins,data_over_mc_width=dm['core_spread']/w if w>0 else float('nan'),median_shift_mc_widths=(dm['median']-mm['median'])/w if w>0 else float('nan'),status='usable' if d.sum()>=100 and raw>=100 and refraw>=100 and np.isfinite(score) else 'insufficient support'))
                            if half==1:curves.append((scenario+' · other MC half',m,'#d45c18' if scenario=='frozen' else '#22946b'))
                    for label,h,color in curves:
                        if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,label=label,color=color)
                    qs=[base.diag.quantile(h[1:-1],edges,q) for _,h,_ in curves if h.sum()>0 for q in (.01,.99)];qs=[q for q in qs if np.isfinite(q)]
                    if qs:
                        lo,hi=min(qs),max(qs);pad=max(.001,(hi-lo)*.08);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                    ax.set(title=side+' CD proton angle · '+name,xlabel=variable['label'],ylabel='Unit-area density');ax.legend(fontsize=7);ax.grid(alpha=.2)
            fig.suptitle(period+' · '+base.TOPS[top]+' · '+channel+f' · proton-angle split {split["threshold"]:.2f}°\nSame observable above/below; fixed common DATA windows, no subdivision refits')
            folder=out/'plots'/period/'angle_subdivisions';folder.mkdir(parents=True,exist_ok=True);path=folder/(channel+'_'+base.TOPS[top]+'.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    return rows,links

def fit_period(period,worker,budget,work,out,spec,cache,data,norm,start,data_arrays,args):
    import matplotlib.pyplot as plt
    one=copy.deepcopy(spec);one['periods']=[period];folder=out/period;folder.mkdir();base.dump(folder/'spec.json',one);pw=work/period;pw.mkdir();pool=base.Pool(worker,folder/'spec.json',pw,args.workers,budget)
    engine=base.Engine(pool,one,cache,data,norm,base.TOPS,args.min_mc,args.background_transfer_uncertainty)
    trials=[start+[float(x)] for x in np.linspace(0,args.max_sigma,17)]
    pars,scores,_,hist,_=engine.evaluate(trials,'coarse');chosen=pick(pars,scores,start);best=pars[chosen][20];scan=[dict(extra_sigma=p[20],**s) for p,s in zip(pars,scores) if list(p[:20])==start];base.discard_histograms(hist)
    step=args.max_sigma/16;trials=[start+[float(x)] for x in sorted(set([0.,best,*np.linspace(max(0,best-step),min(args.max_sigma,best+step),9)]))]
    pars,scores,_,hist,_=engine.evaluate(trials,'refine');chosen=pick(pars,scores,start);best=pars[chosen][20];scan += [dict(extra_sigma=p[20],**s) for p,s in zip(pars,scores) if list(p[:20])==start];base.discard_histograms(hist)
    final={'frozen':start+[0.],'angle':start+[best]};pars,scores,records,hist,selected=engine.evaluate(list(final.values()),'final');idx={s:pars.index(tuple(p)) for s,p in final.items()};mult={s:scores[i]['background_multiplier'] for s,i in idx.items()}
    PLOT_LOCK.acquire()
    try:
        base.table(folder/'scan.csv',scan);base.table(folder/'projection_scores.csv',records)
        base.dump(folder/'angle_candidate.json',dict(period=period,frozen_parameters=dict(zip(base.PARAMS,start)),extra_parameter=EXTRA,value=best,boundary=bool(best==args.max_sigma),law='sigma_p_extra/p = value * clamp((55-theta_deg)/20,0,1); independent widths added in quadrature; original reconstructed angle',production_adopted=False))
        replay=base.replay_check(engine,hist,selected,pars);base.table(folder/'nominal_replay.csv',replay)
        if any(abs(r['pull'])>7 for r in replay):raise ValueError('Nominal cache replay failed')
        summary=[];links=[];fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
        for ax,ch in zip(axes,('epgamma','epgammagamma')):
            for j,s in enumerate(('frozen','angle')):
                values=[scores[idx[s]][ch+'_'+metric] for metric in ('train','validation_mc','heldout_observables')];ax.bar(np.arange(3)+(j-.5)*.35,values,.35,label=s,color=['#d45c18','#22946b'][j])
                for metric,value in zip(('train','validation_mc','heldout_observables'),values):summary.append(dict(period=period,scenario=s,channel=ch,metric=metric,score=value,extra_sigma=final[s][20],at_upper_bound=best==args.max_sigma))
            ax.set(xticks=range(3),xticklabels=['Training','Other MC half','Reserved'],ylabel='Shape score (lower is better)',title=ch);ax.legend();ax.grid(axis='y',alpha=.2)
        path=out/'plots'/period/'scores.png';path.parent.mkdir(parents=True,exist_ok=True);fig.suptitle(period+' · frozen joint response versus one CD-angle parameter');fig.savefig(path,dpi=120);plt.close(fig);links.append(path.relative_to(out).as_posix())
        # Inclusive paired shape plots include FD–FD as an unchanged-track control.
        for ch in ('epgamma','epgammagamma'):
            for top in range(3):
                fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained')
                for vi,v in enumerate(one['variables']):
                    ax=axes.flat[vi];n=v['bins']+2;edges=base.diag.variable_edges(v);d=data[period+('__data' if ch=='epgamma' else '__pi0_data')][top,0,vi,:n];curves=[('DATA',d,'#205cab')]
                    for s in ('frozen','angle'):
                        m=hist[period+'__pi0_mc'][idx[s],1,top,0,vi,:n] if ch=='epgammagamma' else engine.weighted(hist,period,idx[s],1,top,0,vi,mult[s])[0]
                        curves.append((s+' · other MC half',m,'#d45c18' if s=='frozen' else '#22946b'))
                    for label,h,color in curves:
                        if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,label=label,color=color)
                    qs=[base.diag.quantile(h[1:-1],edges,q) for _,h,_ in curves if h.sum()>0 for q in (.01,.99)];qs=[q for q in qs if np.isfinite(q)]
                    if qs:
                        lo,hi=min(qs),max(qs);pad=max(.001,(hi-lo)*.08);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                    if v['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_ in curves]))
                    ax.set(xlabel=v['label'],ylabel='Unit-area density');ax.legend(fontsize=7)
                for ax in axes.flat[7:]:ax.axis('off')
                fig.suptitle(period+' · '+base.TOPS[top]+' · '+ch);path=out/'plots'/period/'inclusive'/f'{ch}_{base.TOPS[top]}.png';path.parent.mkdir(exist_ok=True);fig.savefig(path,dpi=100);plt.close(fig);links.append(path.relative_to(out).as_posix())
        rows,ll=low_high(period,engine,worker,pool,final,idx,mult,data_arrays,out);links+=ll;base.table(folder/'angle_residuals.csv',rows)
        return summary,links
    finally:
        base.discard_histograms(hist)
        PLOT_LOCK.release()

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joint-run',required=True,type=Path);p.add_argument('--frozen-study',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--workers',default=10,type=int);p.add_argument('--scratch-root',default=Path('/scratch/thayward/dvcs_particle_response'),type=Path);p.add_argument('--max-sigma',default=.08,type=float);p.add_argument('--min-mc',default=100,type=int);p.add_argument('--background-transfer-uncertainty',default=.1,type=float);args=p.parse_args(argv)
    if not 1<=args.workers<=10 or not 0<args.max_sigma<=.24 or args.min_mc<1 or not np.isfinite(args.background_transfer_uncertainty) or args.background_transfer_uncertainty<0:p.error('Invalid workers, width bound, MC support or background uncertainty')
    out=args.output.resolve();run=args.joint_run.resolve()
    if out.exists():p.error('Output exists; choose a fresh directory')
    meta,spec,cm,reference,ref,cache,starts,data_arrays=load_inputs(run,args.frozen_study.resolve())
    if args.min_mc!=meta['arguments']['min_mc']:p.error('Keep original MC support for a controlled comparison')
    if meta['arguments']['fit_topologies']!=base.TOPS:p.error('Frozen joint response must cover all three topologies')
    if args.background_transfer_uncertainty!=meta['arguments']['background_transfer_uncertainty']:p.error('Keep the original background nuisance uncertainty for a controlled comparison')
    base.TRAIN_OBS=[0,3,4,5];base.VALID_OBS=[2,6]
    spec['reference_flows']=ref['flows'];spec['joint_control_weight']=.5
    with np.load(reference/'histograms.npz',allow_pickle=False) as ar:data={p+'__data':ar[p+'__data'].copy() for p in spec['periods']}
    data.update(base.load_joint_control_data(Path(meta['arguments']['joint_control_study']),spec,reference))
    norm={(r['period'],base.TOPS.index(r['topology'])):r for r in csv.DictReader((reference/'normalizations.csv').open())}
    out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));worker=compile_angle_worker(out,reference,work)
    base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},joint_manifest_sha256=base.digest(run/'manifest.json'),frozen_manifest_sha256=base.digest(args.frozen_study/'manifest.json'),sources={n:base.digest(base.SOURCE/n) for n in SOURCES},isolated_sources={n:base.digest(out/'backend'/n) for n in SOURCES},study_source_sha256=base.digest(Path(__file__)),extra_parameter=EXTRA,angle_edges_deg=[35,55],frozen_response=True))
    budget=base.Budget(args.workers);summary=[];links=[]
    print('[study] Existing MC caches and DATA arrays only; 17-point batched scan plus local refinement; one extra width per period; no production changes',flush=True)
    # Native evaluations share the global worker budget; serialize Matplotlib only.
    with ThreadPoolExecutor(max_workers=min(args.workers,len(spec['periods']))) as executor:
        futures={executor.submit(fit_period,period,worker,budget,work,out,spec,cache,data,norm,starts[period],data_arrays[period],args):period for period in spec['periods']}
        try:
            for future in as_completed(futures):
                ss,ll=future.result();summary+=ss;links+=ll;print('[done] '+futures[future],flush=True)
        except BaseException:
            budget.cancel_all()
            for future in futures:future.cancel()
            raise
    base.table(out/'score_comparison.csv',summary)
    audit=[]
    for period in spec['periods']:
        audit.extend(csv.DictReader((out/period/'photon_order_audit.csv').open()))
    base.table(out/'photon_order_audit.csv',audit)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>CD-proton angle test</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>One CD-proton angle-dependent width</h1><p>The saved joint response is frozen. Only one extra momentum-width amplitude per running period is fitted, simultaneously using epgamma signal+leakage and epgammagamma AAOGEN shapes, across all three topologies. The extra fractional width is amplitude × clamp((55−theta_deg)/20,0,1), evaluated at the original reconstructed CD-proton angle and added in quadrature to its existing width. This smooth ramp is a diagnostic choice motivated by the preceding angle split, not an established detector calibration. All biases and photon responses remain fixed; photon roles share their detector response.</p><p>Only the training objective (Delta_phi, pTmiss, Emiss2, Mx2, MC half 0) selects the amplitude. Other-half and reserved theta_gamma_gamma/Mx2_2 scores are validation. DATA are shared, so this is not independent DATA validation. The existing constrained background nuisance is profiled identically for both alternatives. Zero amplitude is included, preventing a worse training optimum; improvement alone does not establish validation success. Inclusive and low/high overlays use MC half 1. Compare CD–FT improvements against possible CD–FD regressions and boundary flags. Fixed saved cuts are retained: no production acceptance or cross-section correction is adopted.</p><p>Photon-role energy ordering remains a diagnostic caveat: photon_order_audit.csv tabulates both-FD DATA/control-MC energy rankings, after the same full-global selection, before exclusivity. A ranking reversal needs source-level investigation; it is not evidence for separate photon-role responses. This CD-proton test does not use role-specific energy subdivisions.</p><p><a href="score_comparison.csv">Scores</a> · <a href="photon_order_audit.csv">Photon ordering check</a></p>'+''.join('<p><a href="'+x+'">'+x+'</a></p>' for x in links))
    print('[done] '+str(out/'REPORT.html')+'; transfer this report directory; large arrays remain on scratch',flush=True)

if __name__=='__main__':main()
