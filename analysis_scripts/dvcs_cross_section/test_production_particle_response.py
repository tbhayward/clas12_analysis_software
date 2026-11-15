#!/usr/bin/env python3
"""Production response: native math, real ROOT streaming, routing and cache reuse."""
from __future__ import annotations
import argparse,copy,json,os,shlex,shutil,subprocess,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import uproot
import prepare_production_particle_response as prep
from production_response_stream import stream_file
import production_response_io as routing
SOURCE=Path(__file__).resolve().parent
REPO=prep.base.SOURCE

FIXTURE=r'''
#include "production_response_kinematics.h"
#include "global_cuts.h"
#include <nlohmann/json.hpp>
#include <iostream>
using namespace production_response;
int main(int argc,char** argv){
 if(argc<2)return 1;GlobalCutConfig cfg;cfg.enable_auxiliary_fiducial_cuts=false;cfg.enable_analysis_phase_space_cuts=false;cfg.enable_pTmiss_cut=false;cfg.t1_abs_max=100;cfg.open_angle_min_deg=0;cfg.enable_reconstructed_electron_p_cut=false;write_global_cuts_config_json(argv[1],cfg);
 nlohmann::json result=nlohmann::json::object();
 for(bool control:{false,true})for(bool swap:{false,true}){
  std::array<double,NC> v{};v[3]=2;v[4]=1;v[5]=5.3;v[6]=.22;v[7]=.1;v[8]=.65;v[9]=1.;v[10]=3.;v[11]=3.;v[12]=.25;v[13]=2.;F one=particle(2.8,.18,2.1,0),two=particle(.8,.055,2.5,0);F sum=one+two;
  if(control){v[11]=sum.p.norm();v[12]=angle(sum.p,{0,0,1});v[13]=wrap(std::atan2(sum.p.y,sum.p.x));}
  auto values=calculate(v,10.604,false,control);nlohmann::json row;
  for(size_t i=0;i<values.size();++i)row[names()[i]]=values[i];
  const std::array<std::string,9> primary={"e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi"};for(int i=0;i<9;++i)row[primary[i]]=v[5+i];row["detector1"]=2.;row["detector2"]=control?-1.:1.;
  if(control){if(swap)std::swap(one,two);F beam={{0,0,std::sqrt(10.604*10.604-ME*ME)},10.604},e=particle(v[5],v[6],v[7],ME);for(int i=0;i<2;++i){F g=i?two:one;row["gamma_phi"+std::to_string(i+1)]=trento((beam-e).p,e.p,g.p);row["open_angle_egamma"+std::to_string(i+1)]=angle(e.p,g.p)*180/PI;}row["Mh_gammagamma"]=std::sqrt(sum.m2());row["detector_gamma1"]=swap?0.:1.;row["detector_gamma2"]=swap?1.:0.;}
  row["Mx2"]=row["Mx2"].get<double>()+.000123;row["phi2"]=wrap(row["phi2"].get<double>()+.00005);result[(control?"control":"signal")+std::string(swap?"_swapped":"")]=row;
 }
 std::cout<<result.dump();
}
'''

def native_build(tmp):
    rc=os.environ.get('ROOT_CONFIG','root-config');flags=shlex.split(subprocess.check_output([rc,'--cflags'],text=True));libs=shlex.split(subprocess.check_output([rc,'--libs'],text=True));cxx=os.environ.get('CXX','c++');includes=['-I'+str(SOURCE),'-I'+str(REPO)]
    fixture=tmp/'fixture.cpp';fixture.write_text(FIXTURE);exe=tmp/'fixture';worker=tmp/'production_response_writer'
    subprocess.run([cxx,*flags,'-O2','-std=c++17',*includes,str(fixture),str(REPO/'global_cuts.cpp'),'-o',str(exe),*libs],check=True)
    subprocess.run([cxx,*flags,'-O2','-std=c++17',*includes,str(SOURCE/'production_response_writer.cpp'),str(REPO/'global_cuts.cpp'),'-o',str(worker),*libs],check=True)
    text=subprocess.check_output([str(exe),str(tmp)],text=True);rows=json.loads(text[text.index('{'):]);return worker,rows,json.loads((tmp/'global_cuts_config.json').read_text())

def write_tree(path,row,n=12):
    path.parent.mkdir(parents=True,exist_ok=True);a={k:np.full(n,v,dtype='<f8') for k,v in row.items()};a.update(runnum=np.full(n,5000.,dtype='<f8'),weight=np.linspace(.2,2,n),event_sequence=np.arange(n,dtype=np.uint64)+np.uint64(2**63+11),evnum=np.arange(n,dtype=np.int64)+100,helicity=np.arange(n,dtype=np.int32)%3-1,gen_reference=np.linspace(.01,.99,n))
    with uproot.recreate(path) as f:
        tree=f.mktree('PhysicsEvents',{k:v.dtype for k,v in a.items()});tree.extend(a)
    return a

def arrays(path):
    with uproot.open(path) as f:return f['PhysicsEvents'].arrays(library='np')

def execute(tmp,worker,source,globals_,pars,control=False,label='modified',batch=4,seed=1729,period='fa18_inb',beam=10.604):
    job=dict(source=str(source),target=str(tmp/(label+'.root')),notes=str(tmp/(label+'.bin')),summary=str(tmp/(label+'.json')),period=period,control=control,beam_energy=beam,parameters=pars,seed=seed,global_cuts=globals_,note_probability=1.)
    with (tmp/(label+'.log')).open('w') as log:meta=stream_file(job,worker,log,batch)
    return arrays(job['target']),meta,job

def minkowski(a):return a[0]*a[0]-np.dot(a[1:],a[1:])
def vector(p,t,f,m):return np.r_[np.sqrt(p*p+m*m),p*np.array([np.sin(t)*np.cos(f),np.sin(t)*np.sin(f),np.cos(t)])]
def independent(a,beam=10.604):
    p=lambda prefix,m:vector(a[prefix+'_p'][0],a[prefix+'_theta'][0],a[prefix+'_phi'][0],m)
    k=vector(np.sqrt(beam**2-.00051099895**2),0,0,.00051099895);t=np.array([.9382720813,0,0,0]);e=p('e',.00051099895);pr=p('p1',.9382720813);g=p('p2',0);q=k-e;miss=k+t-e-pr-g
    x=-minkowski(q)/(2*t[0]*q[0]);normal1=np.cross(q[1:],e[1:]);normal2=np.cross(q[1:],g[1:]);phi=np.arctan2(np.dot(np.cross(e[1:],g[1:]),q[1:])*np.linalg.norm(q[1:]),np.dot(normal1,normal2))%(2*np.pi)
    return dict(Q2=-minkowski(q),x=x,W=np.sqrt(minkowski(q+t)),y=q[0]/k[0],t1=minkowski(t-pr),Emiss2=miss[0],Mx2=minkowski(miss),Mx2_1=minkowski(k+t-e-pr),Mx2_2=minkowski(k+t-e-g),pTmiss=np.linalg.norm(miss[1:3]),phi2=phi)

def run_native(tmp):
    worker,rows,globals_=native_build(tmp);source=tmp/'mc/fa18_inb/dvcs_mc/50nA/parts/signal.rec_only.root';original=write_tree(source,rows['signal']);zero=[0.]*20
    out,meta,_=execute(tmp,worker,source,globals_,zero,label='zero')
    for k in original:assert np.array_equal(original[k],out[k],equal_nan=True),k
    assert not np.any(out['response_applied']) and meta['written']==len(original['x'])
    pars=zero.copy();pars[3]=.01;pars[11]=-.02;pars[15]=-.025;pars[13]=.002;pars[14]=.003
    shifted,meta,job=execute(tmp,worker,source,globals_,pars);second,_,other=execute(tmp,worker,source,globals_,pars,label='other_batch',batch=7)
    for k in shifted:assert np.array_equal(shifted[k],second[k],equal_nan=True),k
    assert Path(job['notes']).read_bytes()==Path(other['notes']).read_bytes()
    for k in ('event_sequence','evnum','helicity','weight','runnum','detector1','detector2','gen_reference'):assert np.array_equal(original[k],shifted[k]),k
    before,after=independent(original),independent(shifted)
    for k in before:
        expected=original[k][0]+after[k]-before[k]
        if k=='phi2':expected%=2*np.pi
        assert np.isclose(shifted[k][0],expected,atol=2e-11), (k,shifted[k][0],expected)
    assert not np.isclose(shifted['phi2'][0],original['phi2'][0],atol=1e-7)
    for period,beam in [('sp18_inb',10.594),('sp18_out',10.594),('fa18_out',10.604),('sp19_inb',10.1998)]:
        changed,_,_=execute(tmp,worker,source,globals_,pars,label=period,period=period,beam=beam)
        old_values,new_values=independent(original,beam),independent(changed,beam)
        for name in ('x','Q2','W','y','phi2'):
            expected=original[name][0]+new_values[name]-old_values[name]
            if name=='phi2':expected%=2*np.pi
            assert np.isclose(changed[name][0],expected,atol=2e-11),(period,name)

    for name in ('e_phi','p1_phi','p2_phi','phi1','phi2','phih','phiR','Delta_phi'):assert np.all((shifted[name]>=0)&(shifted[name]<2*np.pi)),name
    notes=np.fromfile(job['notes'],dtype='<f8').reshape(-1,22,31);assert notes.shape[0]==len(original['x'])
    for j,k in enumerate(prep.base.OBS):assert np.allclose(notes[:,1,j+1],shifted[k]),k
    # Shared detector response is independent of photon role for deterministic biases.
    pars=[0.]*20;pars[15]=-.025;pars[19]=.01;pars[3]=.005;pars[11]=-.02
    controls=[]
    for label in ('control','control_swapped'):
        path=tmp/(label+'.rec_only.root');write_tree(path,rows[label]);a,m,j=execute(tmp,worker,path,globals_,pars,True,label+'_changed');assert m['unresolved']==0;controls.append(a)
    for k in ('p2_p','p2_theta','p2_phi','Mh_gammagamma','Mx2','Emiss2'):assert np.allclose(controls[0][k],controls[1][k],atol=3e-12),k
    assert np.all(controls[0]['Mh_gammagamma']<rows['control']['Mh_gammagamma'])
    broken=copy.deepcopy(rows['control']);broken['gamma_phi1']=float('nan');path=tmp/'unresolved.rec_only.root';old=write_tree(path,broken);a,m,j=execute(tmp,worker,path,globals_,pars,True,'unresolved');assert m['unresolved']==len(a['x']) and not a['response_applied'].any()
    for k in old:assert np.array_equal(old[k],a[k],equal_nan=True),k
    # A future physics branch cannot silently remain stale.
    bad=dict(rows['signal'],unimplemented_physics=2.);path=tmp/'unknown.rec_only.root';write_tree(path,bad)
    try:execute(tmp,worker,path,globals_,pars,label='unknown')
    except (RuntimeError,BrokenPipeError):pass
    else:raise AssertionError('Unknown physics branch accepted')
    print('PASS: real ROOT streaming, exact zero replay, independent Q2/x/W/y/t/phi/missing checks, all 22 note scenarios, batch invariance, uint64 identity/weights, independent photon detectors/order, unresolved preservation and unknown-branch rejection')
    return worker,rows,globals_,source

def cache_and_routing(tmp,worker,rows,globals_,source):
    root=tmp/'driver';root.mkdir();model=root/'fitted';model.mkdir();globals_path=root/'globals.json';globals_path.write_text(json.dumps(globals_))
    for p in ('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb'):
        folder=model/'response'/p;folder.mkdir(parents=True);pars=dict.fromkeys(prep.base.PARAMS,0.);pars['proton_CD_p_bias_fraction']=-.02;(folder/'candidate_parameters.json').write_text(json.dumps(dict(parameters=pars)))
    (model/'source_spec.json').write_text(json.dumps(dict(seed=1729)));hashes={n:prep.base.digest(REPO/n) for n in ('particle_response_core.h','particle_control_core.h')};(model/'signature.json').write_text(json.dumps(dict(sources=hashes)))
    gen=root/'generated.root';write_tree(gen,rows['signal']);data=root/'data.root';write_tree(data,rows['signal']);inv=[dict(sample=s,tag=t,tree='PhysicsEvents',files=[str(path)]) for s,t,path in [('mc','DVCS_Fa18_inb_rec',source),('current_mc','DVCS_Fa18_inb_50nA_rec_current',source),('gen','DVCS_Fa18_inb_gen',gen),('data','DVCS_Fa18_inb_data',data)]];inventory=root/'inventory.json';inventory.write_text(json.dumps(inv));args=argparse.Namespace(response_run=model,output=root/'output',scratch_root=root/'scratch',global_cuts=globals_path,inventory=inventory,workers=2,note_events=1000,fit_response=False)
    # SOURCE points at the installed/staged code; only the test worker is copied.
    testsource=root/'source';testsource.mkdir()
    for path in SOURCE.iterdir():
        if path.is_file() and path.suffix in ('.py','.h','.cpp'):shutil.copy2(path,testsource/path.name)
    for n in ('particle_response_core.h','particle_control_core.h','global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(REPO/n,testsource/n)
    shutil.copy2(worker,testsource/'production_response_writer')
    with patch.object(prep,'SOURCE',testsource):
        prep.prepare(args);first=json.loads((args.output/'input_manifest.json').read_text());assert len(first['routes'])==1 and len(next(iter(first['routes'].values()))['roles'])==2
        target=Path(next(iter(first['routes'].values()))['target']);stamp=target.stat().st_mtime_ns
        with patch('production_response_stream.stream_file',side_effect=AssertionError('Cache unexpectedly rebuilt')):prep.prepare(args)
        assert target.stat().st_mtime_ns==stamp
        moved=root/'moved_fitted';shutil.copytree(model,moved);args.response_run=moved
        with patch('production_response_stream.stream_file',side_effect=AssertionError('Moved model snapshot rebuilt cache')):prep.prepare(args)
        assert target.stat().st_mtime_ns==stamp
        change=moved/'response/sp19_inb/candidate_parameters.json';updated=json.loads(change.read_text());updated['parameters']['proton_CD_p_bias_fraction']=-.03;change.write_text(json.dumps(updated))
        with patch('production_response_stream.stream_file',side_effect=AssertionError('Other period rebuilt cache')):prep.prepare(args)
        last=json.loads((args.output/'input_manifest.json').read_text());assert last['model_id']!=first['model_id'] and target.stat().st_mtime_ns==stamp;first=last
        changed_globals=dict(globals_);changed_globals['open_angle_min_deg']=1.;globals_path.write_text(json.dumps(changed_globals))
        with patch('production_response_stream.stream_file',side_effect=AssertionError('Global cuts rebuilt four-vectors')):prep.prepare(args)
        assert target.stat().st_mtime_ns==stamp
    with patch.dict(os.environ,DVCS_PARTICLE_RESPONSE_MANIFEST=str(args.output/'input_manifest.json')):
        routing._state=None;assert routing.routed_files([str(source),str(data),str(gen)])==[str(target),str(data),str(gen)]
        try:routing.routed_files([str(root/'unregistered.root')])
        except ValueError:pass
        else:raise AssertionError('Unregistered input routed silently')
    with patch.dict(os.environ,DVCS_PARTICLE_RESPONSE_DATA_CACHE=str(root/'data_cache')):
        a=routing.cached_data_arrays(dict(input=1),lambda:{'CD_FD':{'Mx2':np.arange(10,dtype=np.float32)}})
        b=routing.cached_data_arrays(dict(input=1),lambda:(_ for _ in ()).throw(AssertionError('DATA reread')));assert np.array_equal(a['CD_FD']['Mx2'],b['CD_FD']['Mx2']) and not b['CD_FD']['Mx2'].flags.writeable
    import finalize_production_particle_response as finalizer
    csv_path=root/'result.csv';csv_path.write_text('bin,value\n1,2\n');cuts_path=root/'cuts.json';cuts_path.write_text('{}')
    result=finalizer.finalize(args.output,csv_path,cuts_path);assert result['model_id']==first['model_id'] and result['csv_sha256']==finalizer.digest(csv_path) and csv_path.with_suffix('.response.json').is_file()
    print('PASS: full driver, frozen model/provenance, deduplicated reference/current source, checked cache reuse, shared routing, untouched DATA/GEN and unchanged DATA-array reuse')

def data_reader_and_note(tmp,rows,globals_):
    import plot_production_particle_response as note
    production=note.production
    path=tmp/'data/fa18_inb/pi0_data/parts/data.root';row=dict(rows['control']);row['detector2']=row['detector_gamma1'];write_tree(path,row)
    production.ACTIVE_GLOBAL_CUTS=production.load_global_cuts_config(tmp/'global_cuts_config.json')
    with patch.dict(os.environ,DVCS_PARTICLE_RESPONSE_DATA_CACHE=str(tmp/'reader_cache')):
        a=production.load_selected_arrays_for_file(str(path),'fa18_inb',production.TOPOLOGIES,'1 MB',production.PI0_SHAPE_VARIABLES)
        with patch.object(production,'load_selected_arrays_uproot',side_effect=AssertionError('DATA ROOT reread')):
            b=production.load_selected_arrays_for_file(str(path),'fa18_inb',production.TOPOLOGIES,'1 MB',production.PI0_SHAPE_VARIABLES)
        assert 'Mh_gammagamma' in b['CD_FD'] and np.array_equal(a['CD_FD']['Mh_gammagamma'],b['CD_FD']['Mh_gammagamma'])
    rng=np.random.default_rng(41);period=production.PERIODS[0];origin=tmp/'note';origin.mkdir();routes={};data={}
    means=[np.pi,2.6,.3,.08,0.,0.,.88];sigmas=[.045,.2,.15,.03,.15,.004,.2]
    for channel in ('epgamma','epgammagamma'):
        data[channel]=[]
        for group in range(3):
            values={name:rng.normal(means[j],sigmas[j],5000 if channel=='epgamma' else 1000) for j,name in enumerate(prep.base.OBS)}
            if channel=='epgammagamma':values['m_gamma_gamma']=rng.normal(.134,.009,len(values['Mx2']))
            data[channel].append(values)
    for sample in ('mc','background','pi0_mc'):
        n=6000;z=np.zeros((n,22,31));z[:,:,0]=np.repeat(np.arange(3),n//3)[:,None];z[:,:,8]=1;z[:,:,10]=(np.arange(n)%2)[:,None];z[:,:,13]=1
        for j,name in enumerate(('x','Q2','t1','detector1','detector2','e_p','e_theta','e_phi','p1_p','p1_theta','p1_phi','p2_p','p2_theta','p2_phi','open_angle_ep2','pTmiss')):
            z[:,:,14+j]= {'x':.2,'Q2':2.,'t1':-.2,'detector1':2.,'detector2':1.,'e_p':5.,'e_theta':.2,'e_phi':.1,'p1_p':.6,'p1_theta':1.,'p1_phi':3.,'p2_p':3.,'p2_theta':.2,'p2_phi':2.,'open_angle_ep2':10.,'pTmiss':.1}[name]
        z[:,:,30]=5000
        for j,name in enumerate(prep.base.OBS):
            latent=rng.normal(means[j],sigmas[j],n);z[:,0,j+1]=means[j]+(latent-means[j])*.65+sigmas[j]*.35
            for scenario in range(1,22):z[:,scenario,j+1]=latent+sigmas[j]*(.02 if scenario==1 else .03+.002*scenario)
        z[:2000,:,17]=1;z[4000:,:,18]=0
        z[:,:,12]=rng.normal(.134,.009,n)[:,None];file=origin/(sample+'.bin');z.astype('<f8').tofile(file);routes[sample]=dict(period=period.key,roles=[dict(sample=sample)],notes=str(file),summary=dict(note_shape=list(z.shape),note_probability=1.,unresolved=0))
    physical_names=('x','Q2','t1','detector1','detector2','e_p','e_theta','e_phi','p1_p','p1_theta','p1_phi','p2_p','p2_theta','p2_phi','open_angle_ep2','pTmiss')
    full_routes={}
    for sample,r in routes.items():
        z=np.fromfile(r['notes'],dtype='<f8').reshape(-1,22,31)
        for scenario in (0,1):
            frame={name:z[:,scenario,14+j] for j,name in enumerate(physical_names)};frame['runnum']=z[:,scenario,30]
            frame.update({('theta_pi0_pi0' if sample=='pi0_mc' and name=='theta_gamma_gamma' else name):z[:,scenario,j+1] for j,name in enumerate(prep.base.OBS)})
            if sample=='pi0_mc':frame.update(detector_gamma1=frame['detector2'],Mh_gammagamma=z[:,scenario,12])
            path=origin/(sample+str(scenario)+'.root')
            with uproot.recreate(path) as f:tree=f.mktree('PhysicsEvents',{k:v.dtype for k,v in frame.items()});tree.extend(frame)
        full_routes[str(origin/(sample+'0.root'))]=dict(r,target=str(origin/(sample+'1.root')))
    routes=full_routes
    model=dict(seed=1729,periods={period.key:dict(parameters=[.01,0,.001,-.01]*5)})
    manifest=dict(model_id='test',model=model,routes=routes,global_cuts=globals_,global_cuts_identity=dict(path=str(tmp/'global_cuts_config.json')));(origin/'input_manifest.json').write_text(json.dumps(manifest));cut={}
    for channel in ('dvcs','eppi0'):
        for top in production.TOPOLOGIES:cut[production.cut_json_key(channel,period,top)]=dict(data={},mc={})
    (origin/'cuts.json').write_text(json.dumps(cut))
    with patch.object(production,'PERIODS',(period,)),patch.object(note,'data_arrays',return_value=data):note.main(['--response-output',str(origin),'--cuts',str(origin/'cuts.json'),'--output',str(origin/'analysis_note')])
    assert len(list((origin/'analysis_note').rglob('*.png')))==10 and not list((origin/'analysis_note').rglob('*.pdf'))
    with np.load(origin/'analysis_note/plot_histograms.npz') as h:assert len(h.files)>60
    assert (origin/'analysis_note/full_statistics_manifest.json').is_file() and (origin/'analysis_note/influence_scores.csv').is_file()
    with patch.object(note,'data_arrays',side_effect=AssertionError('Replot reread DATA')),patch.object(note,'note_histograms',side_effect=AssertionError('Replot reread MC')):
        note.main(['--replot','--output',str(origin/'analysis_note')])
    print('PASS: actual production DATA reader/cache with saved control mass, final-window note scores, all parameter ablations, component overlays, PNG-only output and histogram-only replot')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',action='store_true',help='Compile the native evaluator and test actual ROOT files');a=p.parse_args()
    # Source-level integration checks complement actual numerical/I/O tests.
    text=(SOURCE/'main.cpp').read_text();assert text.index('prepare_production_particle_response(particle_response')<text.index('// --------- Global bin-averaged kinematics')<text.index('// --------- Current-response calibration')
    assert 'current_mc' in text and 'pi0_mc' in text and 'rad_mc' in text
    for path in SOURCE.glob('*production*response*.py'):compile(path.read_text(),str(path),'exec')
    print('PASS: production stage precedes means/current/counting; all reconstructed maps included; Python modules compile')
    if a.native:
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory);worker,rows,globals_,source=run_native(tmp);cache_and_routing(tmp,worker,rows,globals_,source);data_reader_and_note(tmp,rows,globals_)
if __name__=='__main__':main()
