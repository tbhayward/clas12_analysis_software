#!/usr/bin/env python3
"""Fit independent joint particle responses from this production run's inputs."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,concurrent.futures,copy,fcntl,hashlib,json,shutil,tempfile,uuid
from pathlib import Path
import numpy as np
import uproot
import study_selected_particle_response as selected
from production_response_stream import stream_fit_cache
base=selected.base
SOURCE=Path(__file__).resolve().parent
PERIODS=('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb')

def digest_json(x):return hashlib.sha256(json.dumps(x,sort_keys=True).encode()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def stamp(p):
    s=Path(p).stat();return dict(path=os.path.abspath(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
def atomic(p,x):
    p=Path(p);temp=p.with_name(p.name+'.'+uuid.uuid4().hex+'.tmp');base.dump(temp,x);temp.replace(p)
def period_for(tag):
    return next(p for p in PERIODS if p in tag.lower())
def beam(period):return 10.1998 if period=='sp19_inb' else 10.594 if period.startswith('sp18') else 10.604

def intersect(blocks):
    result={}
    for kind in ('data','mc'):
        result[kind]={}
        for name in sorted(set().union(*(b[kind] for b in blocks))):
            ws=[b[kind][name] for b in blocks if name in b[kind]];lo=max(float(w['cut_low']) for w in ws);hi=min(float(w['cut_high']) for w in ws)
            if hi<=lo:raise ValueError('Initial control windows do not overlap: '+name)
            result[kind][name]=dict(cut_low=lo,cut_high=hi,mode='upper_quantile' if all(w.get('mode')=='upper_quantile' for w in ws) else 'quantile_window')
    return result

def windows_from_production(cuts,period):
    p=next(p for p in selected.cuts.production.PERIODS if p.key==period);tops=selected.cuts.production.TOPOLOGIES
    dvcs=[copy.deepcopy(cuts[selected.cuts.production.cut_json_key('dvcs',p,t)]) for t in tops];ctrl=[copy.deepcopy(cuts[selected.cuts.production.cut_json_key('eppi0',p,t)]) for t in tops]
    # Pure detector pairs inherit their actual production window. The mixed
    # CD pair uses their intersection so neither photon ordering chooses a cut.
    return dict(dvcs=dvcs,control=[ctrl[0],ctrl[1],intersect(ctrl[1:]),ctrl[2]])

def fit_coordinates(engine,folder,coarse=2,refine=1):
    best=[0.]*20;grid=base.parameter_grid(base.TOPS,'fit','wide');steps=[]
    for phase,rounds in [('coarse',coarse),('refine',refine)]:
      for iteration in range(rounds):
        for k,levels in grid.items():
            values=levels
            if phase=='refine':
                lo=max((v for v in levels if v<best[k]),default=best[k]);hi=min((v for v in levels if v>best[k]),default=best[k]);values=sorted(set([lo,(lo+best[k])/2,best[k],(hi+best[k])/2,hi]))
            trials=[]
            for value in values:p=list(best);p[k]=value;trials.append(p)
            pars,scores,_,hist,_=engine.evaluate(trials,f'{phase}_{iteration}_{k}');chosen=selected.select_training(pars,scores,trials);best=list(pars[chosen]);steps.append(dict(phase=phase,round=iteration,parameter=base.PARAMS[k],value=best[k],train=scores[chosen]['train']));base.discard_histograms(hist)
            print(f"[response fit] {engine.spec['periods'][0]} {phase} {iteration+1} {base.PARAMS[k]}={best[k]:.6g}, score={scores[chosen]['train']:.4g}",flush=True)
    pars,scores,records,hist,counts=engine.evaluate([[0.]*20,best],'final');index=pars.index(tuple(best));base.table(folder/'fit_steps.csv',steps);base.table(folder/'projection_scores.csv',records)
    atomic(folder/'candidate_parameters.json',dict(period=engine.spec['periods'][0],parameters=dict(zip(base.PARAMS,best)),bias_mode='fit',fit_topologies=base.TOPS,joint_control_fit=True,selection_first=True,unordered_control=True,control_mass_trained=True,training_variables=[base.OBS[k] for k in selected.TRAIN],reserved_variables=[base.OBS[k] for k in selected.RESERVED],background_multiplier=scores[index]['background_multiplier'],coarse_rounds=coarse,refinement_rounds=refine,initial_parameters='zero; no imported fit'))
    atomic(folder/'final_scores.json',dict(parameters=pars,scores=scores));base.discard_histograms(hist);return best

def build_caches(inventory,scratch,worker,globals_,seed,events,workers):
    sources={}
    for entry in inventory:
        if entry['sample'] not in ('data','pi0_data','mc','background','pi0_mc'):continue
        key=(period_for(entry['tag']),entry['sample']);sources.setdefault(key,[]).extend(os.path.abspath(p) for p in entry['files'])
    if set(sources)!={(p,s) for p in PERIODS for s in ('data','pi0_data','mc','background','pi0_mc')}:raise ValueError('Fresh joint fit needs DATA, signal, leakage and both control samples for all five periods')
    def make(item):
        (period,sample),files=item;files=sorted(set(files));identities=[stamp(p) for p in files];control=sample.startswith('pi0');data=sample in ('data','pi0_data');recipe=dict(inputs=identities,period=period,control=control,seed=seed,events=None if data else events,worker_sha256=base.digest(worker));folder=scratch/'fit_caches'/digest_json(recipe);folder.mkdir(parents=True,exist_ok=True)
        with (folder/'lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX);target=folder/'events.bin';marker=folder/'complete.json'
            if marker.is_file():
                saved=read(marker)
                if saved['recipe']!=recipe or stamp(target)!=saved['file']:raise ValueError('Original-event fitting cache changed')
                return (period,sample),dict(path=str(target),probability=saved['probability'])
            entries=0
            for p in files:
                with uproot.open(p) as f:entries+=f['PhysicsEvents'].num_entries
            probability=1. if data else min(1.,events/max(1,entries));temp=folder/('events.'+uuid.uuid4().hex+'.partial');records=0
            with temp.open('wb') as merged,(folder/'cache.log').open('w') as log:
                for i,path in enumerate(files):
                    part=folder/f'part_{i}.partial';job=dict(source=path,target=str(part),period=period,control=control,beam_energy=beam(period),seed=seed,global_cuts=globals_,cache_probability=probability)
                    meta=stream_fit_cache(job,worker,log);records+=meta['records']
                    with part.open('rb') as src:shutil.copyfileobj(src,merged,8*1024*1024)
                    part.unlink()
            if [stamp(p) for p in files]!=identities:raise ValueError('Input changed during fit-cache extraction')
            temp.replace(target);atomic(marker,dict(recipe=recipe,file=stamp(target),records=records,probability=probability));print(f'[response fit cache] {period}/{sample}: {records:,} original events; probability={probability:.5g}',flush=True)
            return (period,sample),dict(path=str(target),probability=probability)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:return dict(executor.map(make,sources.items()))

def run_fit(inventory_path,initial_cuts,globals_path,scratch,out,workers=10,events=800000,seed=1729):
    inventory=read(inventory_path);globals_=read(globals_path);initial=read(initial_cuts);out=Path(out);scratch=Path(scratch);out.mkdir(parents=True,exist_ok=True);scratch.mkdir(parents=True,exist_ok=True)
    files=sorted({os.path.abspath(p) for r in inventory if r['sample'] in ('data','pi0_data','mc','background','pi0_mc') for p in r['files']});signature=dict(inputs=[stamp(p) for p in files],initial_cuts_sha256=base.digest(initial_cuts),globals=globals_,seed=seed,events=events,coarse_rounds=2,refinement_rounds=1,sources={n:base.digest(SOURCE/n) for n in ('fit_production_particle_response.py','study_selected_particle_response.py','selected_response_worker.cpp','particle_response_core.h','particle_control_core.h','production_response_writer.cpp','production_response_stream.py','study_particle_response.py','study_unordered_control_response.py','particle_response_cut_arrays.h','exclusivity_resolution_core.h','exclusivity_resolution_worker.cpp','global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h')},worker_sha256=base.digest(SOURCE/'production_response_writer'))
    marker=out/'fit_signature.json'
    if marker.is_file() and read(marker)!=signature:raise ValueError('Fresh-fit inputs/settings changed; use a new analysis run directory')
    atomic(marker,signature);atomic(out/'source_spec.json',dict(seed=seed,periods=PERIODS,control_to_leakage_exposure_ratio=1.));atomic(out/'signature.json',dict(sources=signature['sources'],pipeline=base.digest(Path(__file__))))
    if all((out/'response'/p/'fit_completed.json').is_file() for p in PERIODS):
        for p in PERIODS:
            folder=out/'response'/p;completed=read(folder/'fit_completed.json')
            if any(base.digest(folder/n)!=h for n,h in completed.items()):raise ValueError('Completed fitted response changed: '+p)
        print('[response fit] all five completed fits verified and reused',flush=True);render_scores(out);return out
    cache=build_caches(inventory,scratch,SOURCE/'production_response_writer',globals_,seed,events,workers)
    reference=out/'reference';(reference/'backend').mkdir(parents=True,exist_ok=True)
    for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(SOURCE/n,reference/'backend'/n)
    work=Path(tempfile.mkdtemp(prefix='production_fit-',dir=scratch));build=out/'worker';build.mkdir(exist_ok=True)
    # Compile once per invocation; scoring/export uses only flat caches, no ROOT I/O.
    if (build/'backend').is_dir():shutil.rmtree(build/'backend')
    worker=selected.compile_worker(build,reference,work);budget=base.Budget(workers);steps=[]
    variables=[dict(name=n,label=n,bins=600 if j==5 else 240,low=lo,high=hi) for j,(n,lo,hi) in enumerate(zip(base.OBS,[np.pi-.6,0,0,0,-1.5,-.04,-1.5],[np.pi+.6,np.pi,5,.8,2.5,.04,4.5]))];mass=dict(name='m_gamma_gamma',label='Two-photon mass (GeV)',bins=140,low=.10,high=.17)
    def period_fit(period):
        folder=out/'response'/period;folder.mkdir(parents=True,exist_ok=True);done=folder/'fit_completed.json'
        if done.is_file():
            if any(base.digest(folder/n)!=h for n,h in read(done).items()):raise ValueError('Completed fit changed: '+period)
            return
        windows=windows_from_production(initial,period);spec=dict(periods=[period],global_cuts=globals_,variables=variables,mass_variable=mass,dvcs_windows={period:windows['dvcs']},control_windows={period:windows['control']});atomic(folder/'spec.json',spec);atomic(folder/'initial_windows.json',windows);pw=work/period;pw.mkdir();pool=base.Pool(worker,folder/'spec.json',pw,3,budget=budget);jobs=[]
        for sample in ('data','pi0_data'):
            jobs.append(dict(period=period,sample='pi0_mc' if sample=='pi0_data' else 'mc',cache=cache[period,sample]['path'],parameters=[0.]*20,export_individual_photons=sample=='pi0_data',export_global_only=False,original_sample=sample))
        arrays={job['original_sample']:selected.arrays_from_export(path) for _,job,path in pool.run('export-selected',jobs,'data_reference')};d,counts=selected.data_histograms(period,arrays,windows,spec);mc={s:cache[period,s] for s in ('mc','background','pi0_mc')};engine=selected.Engine(pool,spec,mc,d,counts,100,.1);fit_coordinates(engine,folder);atomic(done,{n:base.digest(folder/n) for n in ('candidate_parameters.json','final_scores.json','projection_scores.csv')});print('[response fit] '+period+' joint fit complete',flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(5,workers)) as executor:
        futures=[executor.submit(period_fit,p) for p in PERIODS]
        try:
            for f in concurrent.futures.as_completed(futures):f.result()
        except BaseException:
            budget.cancel_all()
            for f in futures:f.cancel()
            raise
    render_scores(out);return out

def render_scores(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    links=[]
    for period in PERIODS:
        summary=read(out/'response'/period/'final_scores.json');scores=summary['scores'];index=0 if len(scores)==1 else 1;fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        for ax,ch in zip(axes,('epgamma','epgammagamma')):
            for i,(label,score,color) in enumerate((('Nominal',scores[0],'#888888'),('Fitted',scores[index],'#22946b'))):ax.bar(np.arange(3)+(i-.5)*.35,[float(score[ch+'_'+k]) if score[ch+'_'+k] is not None else np.nan for k in ('train','validation_mc','reserved')],.35,label=label,color=color)
            ax.set(xticks=range(3),xticklabels=['Training','Other MC half','Reserved'],ylabel='Shape score (lower is better)',title=ch);ax.legend();ax.grid(axis='y',alpha=.2)
        fig.suptitle(period+' — response fitted inside this analysis run');folder=out/'plots';folder.mkdir(exist_ok=True);path=folder/(period+'_fit_scores');fig.savefig(path.with_suffix('.pdf'));fig.savefig(path.with_suffix('.png'),dpi=220);plt.close(fig);links.append('plots/'+path.name+'.png')
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:17px system-ui;max-width:1100px;margin:30px auto}img{max-width:100%}</style><h1>Joint response fitted in the production chain</h1><p>Each period/polarity starts from zero response and has an independent shared-particle fit. Two broad coordinate sweeps and one refinement use training observables/MC half 0 only. Both channels constrain the same response, including the two-photon mass. Other-MC-half and reserved scores do not choose parameters. The DATA are shared. Initial production DATA windows are frozen during fitting; mixed control detector pairs use their order-independent intersection. No prior study or saved fitted parameters are required. Background normalization nuisance is constrained and is not a production event correction.</p>'+''.join('<p><img src="'+p+'"></p>' for p in links)+'<p>Parameter files, selection definitions, projection scores and checkpoints are under response/PERIOD.</p>')
