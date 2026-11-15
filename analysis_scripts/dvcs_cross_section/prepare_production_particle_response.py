#!/usr/bin/env python3
"""Materialize a frozen fitted response once for all production reconstructed MC."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,concurrent.futures,fcntl,hashlib,json,time,uuid
from pathlib import Path
import numpy as np
import uproot
import study_particle_response as base
SOURCE=Path(__file__).resolve().parent
VERSION='production-response-v1'

def digest_json(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def atomic_json(p,value):
    p=Path(p);tmp=p.with_name(p.name+'.'+uuid.uuid4().hex+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def identity(path):
    p=Path(path);s=p.stat();return dict(path=os.path.abspath(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
def period_for(tag):
    for p in ('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb'):
        if p in tag.lower():return p
    raise ValueError('Unknown response period: '+tag)
def models(run):
    run=Path(run);spec=read(run/('source_spec.json' if (run/'source_spec.json').is_file() else 'spec.json'))
    periods=('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb');result={}
    for p in periods:
        f=run/'response'/p/'candidate_parameters.json'
        if not f.is_file():f=run/p/'candidate_parameters.json'
        x=read(f);pars=x['parameters']
        marker=f.parent/'fit_completed.json'
        if marker.is_file() and read(marker).get('candidate_parameters.json')!=base.digest(f):raise ValueError('Fitted response changed after completion: '+str(f))
        if set(pars)!=set(base.PARAMS):raise ValueError('Expected exactly 20 shared particle response parameters: '+str(f))
        a=[float(pars[k]) for k in base.PARAMS]
        if not np.all(np.isfinite(a)) or any(v<0 for i,v in enumerate(a) if i%4!=3):raise ValueError('Invalid particle response: '+str(f))
        result[p]=dict(parameters=a,parameter_names=base.PARAMS,source=identity(f),sha256=base.digest(f),background_multiplier=float(x.get('background_multiplier',1.)))
    # The imported fit must have used the same response laws as the full-statistics writer.
    manifest=read(run/('signature.json' if (run/'signature.json').is_file() else 'manifest.json'));hashes=manifest.get('sources',manifest.get('source_hashes',{}))
    for name in ('particle_response_core.h','particle_control_core.h'):
        if name not in hashes or hashes[name]!=base.digest(SOURCE/name):raise ValueError('Fitted response core differs from current source: '+name)
    return dict(periods=result,seed=int(spec['seed']),source_run=str(run.resolve()),source_spec_sha256=base.digest(run/('source_spec.json' if (run/'source_spec.json').is_file() else 'spec.json')),version=VERSION,control_to_leakage_exposure_ratio=float(spec.get('control_to_leakage_exposure_ratio',1.)))

def prepare(args):
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=True);model=models(args.response_run.resolve());globals_=read(args.global_cuts)
    exe=SOURCE/'production_response_writer'
    if not exe.is_file():raise FileNotFoundError('Build production_response_writer with make -j10 dvcs_analysis')
    sources=('production_response_writer.cpp','production_response_kinematics.h','production_response_config.h','production_response_stream.py','particle_response_core.h','particle_control_core.h','global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h')
    evaluator={n:base.digest(SOURCE/n) for n in sources};evaluator['executable']=base.digest(exe)
    global_snapshot=out/'global_cuts_config.json';atomic_json(global_snapshot,globals_)
    inv=read(args.inventory);entries={};preserved=set();all_sources={};diagnostic=set()
    for item in inv:
        sample=item['sample'];transform=sample not in ('data','pi0_data') and 'gen' not in sample
        for path in item['files']:
            path=os.path.abspath(path);all_sources[path]=identity(path)
            if not transform:preserved.add(path);continue
            period=period_for(item['tag']);control=sample=='pi0_mc';existing=entries.get(path)
            if existing and (existing['period'],existing['control'])!=(period,control):raise ValueError('Conflicting reconstructed source roles: '+path)
            if not existing:
                with uproot.open(path) as f:
                    if 'PhysicsEvents' not in f:raise ValueError('Reconstructed source has no PhysicsEvents: '+path)
                    if 'response_applied' in f['PhysicsEvents']:raise ValueError('Input is already modified: '+path)
                    count=f['PhysicsEvents'].num_entries
                entries[path]=dict(source=path,period=period,control=control,entries=count,roles=[])
            entries[path]['roles'].append(dict(sample=sample,tag=item['tag']))
            if sample in ('mc','background','pi0_mc'):diagnostic.add(path)
    if not entries:raise ValueError('No reconstructed MC in production inventory')
    if set(entries)&preserved:raise ValueError('An input cannot be both generated/DATA and reconstructed')
    totals={}
    for path in diagnostic:
        r=entries[path];role=next(x['sample'] for x in r['roles'] if x['sample'] in ('mc','background','pi0_mc'));key=(r['period'],role);totals[key]=totals.get(key,0)+r['entries']
    scratch=args.scratch_root.resolve();scratch.mkdir(parents=True,exist_ok=True)
    model_id=digest_json(dict(parameters={p:r['parameters'] for p,r in model['periods'].items()},seed=model['seed'],evaluator=evaluator,notes=args.note_events))
    period_ids={p:digest_json(dict(period=p,parameters=r['parameters'],seed=model['seed'],evaluator=evaluator,notes=args.note_events)) for p,r in model['periods'].items()}
    def run_one(r):
        path=r['source'];params=model['periods'][r['period']]['parameters'];source_id=all_sources[path];prob=0.
        if path in diagnostic:
            role=next(x['sample'] for x in r['roles'] if x['sample'] in ('mc','background','pi0_mc'));prob=min(1.,args.note_events/max(1,totals[r['period'],role]))
        cache_id=period_ids[r['period']];parent=scratch/cache_id;parent.mkdir(exist_ok=True);key=digest_json(dict(source=source_id,model=cache_id,note_probability=prob));folder=parent/key;folder.mkdir(exist_ok=True)
        with (folder/'lock').open('w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            root=folder/'PhysicsEvents.root';notes=folder/'note_sample.bin';summary=folder/'summary.json';marker=folder/'complete.json'
            reuse=False
            if marker.is_file():
                m=read(marker);reuse=m['source']==source_id and m['model_id']==cache_id and all(identity(folder/n)==v for n,v in m['outputs'].items())
                if not reuse:raise ValueError('Completed response cache changed: '+str(folder))
            if not reuse:
                token='.'+uuid.uuid4().hex+'.partial';rr=folder/('PhysicsEvents'+token+'.root');nn=folder/('notes'+token);ss=folder/('summary'+token+'.json');jp=folder/('job'+token+'.json')
                job=dict(source=path,target=str(rr),notes=str(nn),summary=str(ss),period=r['period'],control=r['control'],parameters=params,seed=model['seed'],beam_energy=10.1998 if r['period']=='sp19_inb' else 10.594 if r['period'].startswith('sp18') else 10.604,global_cuts=globals_,note_probability=prob)
                atomic_json(jp,job)
                from production_response_stream import stream_file
                print(f"[response] processing {r['period']} {Path(path).name}: {r['entries']:,} events",flush=True)
                last=[time.monotonic()]
                def progress(done,total):
                    now=time.monotonic()
                    if now-last[0]>=30:
                        print(f"[response] {Path(path).name}: {done:,}/{total:,} ({100*done/max(1,total):.1f}%)",flush=True);last[0]=now
                try:
                    with (folder/'writer.log').open('w') as log:stream_file(job,exe,log,progress=progress)
                except Exception as error:
                    for partial in (rr,nn,ss,jp):partial.unlink(missing_ok=True)
                    raise RuntimeError('Response writer failed: '+str(folder/'writer.log')+'\n'+(folder/'writer.log').read_text()[-3000:]) from error
                sm=read(ss)
                if sm['written']!=r['entries'] or sm['entries']!=r['entries']:raise ValueError('Response changed reconstructed entry count')
                with uproot.open(rr) as f:
                    t=f['PhysicsEvents']
                    if t.num_entries!=r['entries'] or 'response_applied' not in t:raise ValueError('Invalid response ROOT output')
                if nn.stat().st_size!=int(np.prod(sm['note_shape']))*8:raise ValueError('Truncated response note sample')
                if identity(path)!=source_id:raise ValueError('Source changed during response production: '+path)
                rr.replace(root);nn.replace(notes);ss.replace(summary);jp.unlink()
                atomic_json(marker,dict(source=source_id,model_id=cache_id,outputs={n:identity(folder/n) for n in ('PhysicsEvents.root','note_sample.bin','summary.json')}))
            sm=read(summary);return path,dict(target=str(root),source_identity=source_id,target_identity=identity(root),notes=str(notes),summary=sm,roles=r['roles'],period=r['period'],control=r['control'],cache_model_id=cache_id,reused=reuse)
    routes={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        ordered=sorted(entries.values(),key=lambda r:(-r['entries']*(1.5 if r['control'] else 1.),r['source']))
        futures=[pool.submit(run_one,r) for r in ordered]
        for future in concurrent.futures.as_completed(futures):
            path,value=future.result();routes[path]=value;print(f'[response] {len(routes)}/{len(entries)} '+('reused ' if value['reused'] else 'built ')+Path(path).name,flush=True)
    manifest=dict(version=VERSION,data_array_cache=os.environ.get('DVCS_PARTICLE_RESPONSE_DATA_CACHE'),model_id=model_id,model=model,evaluator=evaluator,global_cuts=globals_,global_cuts_identity=identity(global_snapshot),routes=routes,preserved_files={p:all_sources[p] for p in sorted(preserved)},inventory_sha256=base.digest(args.inventory))
    atomic_json(out/'response_model.json',model);atomic_json(out/'input_manifest.json',manifest);print('[response] Frozen full-statistics manifest: '+str(out/'input_manifest.json'),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inventory',required=True,type=Path);p.add_argument('--response-run',required=True,type=Path);p.add_argument('--global-cuts',required=True,type=Path);p.add_argument('--scratch-root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--workers',type=int,default=10);p.add_argument('--note-events',type=int,default=50000);a=p.parse_args()
    if not 1<=a.workers<=10 or a.note_events<1000:p.error('workers must be 1–10; note-events at least 1000')
    prepare(a)
if __name__=='__main__':main()
