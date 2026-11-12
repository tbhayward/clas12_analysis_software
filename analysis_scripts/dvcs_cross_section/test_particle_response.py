#!/usr/bin/env python3
"""Physics-kernel, independent-period optimization and optional native ROOT tests."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import numpy as np
import study_particle_response as s

CORE_TEST=r'''
#include "particle_response_core.h"
#include <cassert>
#include <iostream>
using namespace particle_response;
int main(){
 assert(std::abs(trento({0,0,4},{1,0,1},{1,1,2})-PI/4)<1e-12);
 assert(std::abs(trento({0,0,8},{1,0,1},{1,1,2})-PI/4)<1e-12);
 const double E=10.604;F k={{0,0,std::sqrt(E*E-ME*ME)},E},target={{0,0,0},MP};
 F e=particle(5,.2,0,ME),h=k+target-e;double W=std::sqrt(h.m2()),mom=(W*W-MP*MP)/(2*W);
 F g=particle(mom,1.2,1.1,0),p={g.p*(-1),std::sqrt(mom*mom+MP*MP)};
 g=g.boost(h.p*(1/h.e));p=p.boost(h.p*(1/h.e));
 F miss=k+target-e-p-g;assert(std::abs(miss.e)<1e-12&&miss.p.norm()<1e-12);
 Record r{};r.key=123456;
 r.v[3]=2;r.v[4]=0;r.v[5]=5;r.v[6]=.2;r.v[7]=0;
 r.v[8]=p.p.norm();r.v[9]=angle(p.p,{0,0,1});r.v[10]=wrap(std::atan2(p.p.y,p.p.x));
 r.v[11]=g.p.norm();r.v[12]=angle(g.p,{0,0,1});r.v[13]=wrap(std::atan2(g.p.y,g.p.x));
 auto original=derived(r.v,E);int dest[12]={0,1,2,14,15,16,17,18,19,20,21,22};
 for(int i=0;i<12;++i)r.v[dest[i]]=original[i];
 assert(std::abs(original[9])<1e-12&&std::abs(original[10])<1e-20&&original[4]<1e-12);
 assert(std::abs(original[7])<1e-5);assert(std::abs(original[5]-PI)<1e-12);
 std::array<double,NP> pars{};auto zero=vary(r,pars,E,false,original);assert(zero==r.v);
 pars[19]=.01;auto biased=vary(r,pars,E,false,original);
 assert(std::abs(biased[11]-r.v[11]*1.01)<1e-12);
 assert(std::abs(biased[20]-r.v[20]+.01*r.v[11])<1e-12);
 auto direct=derived(biased,E);assert(std::abs(biased[21]-direct[10])<1e-11);
 // Existing stored offsets remain intact while physical shifts are recomputed.
 r.v[21]+=.004;auto offset=vary(r,pars,E,false,original);assert(std::abs(offset[21]-direct[10]-.004)<1e-11);
 pars.fill(0);pars[0]=.01;double sum=0,sum2=0;
 for(int i=0;i<20000;++i){r.key=mix(i);double z=gaussian(r.key,0);sum+=z;sum2+=z*z;}
 assert(std::abs(sum/20000)<.03);assert(std::abs(sum2/20000-1)<.04);
 r.key=42;auto a=vary(r,pars,E,false,original),b=vary(r,pars,E,false,original);assert(a==b);
 assert(a[3]==r.v[3]&&a[4]==r.v[4]&&a[23]==r.v[23]);
 // Shared detector response across topology labels at identical kinematics.
 pars.fill(0);pars[8]=.03;Record cdft=r,cdfd=r;cdft.v[3]=cdfd.v[3]=2;cdft.v[4]=0;cdfd.v[4]=1;
 auto cd1=vary(cdft,pars,E,false,original),cd2=vary(cdfd,pars,E,false,original);assert(cd1[8]==cd2[8]);
 pars.fill(0);pars[12]=.02;Record fdfd=cdfd;fdfd.v[3]=1;
 auto ph1=vary(cdfd,pars,E,false,original),ph2=vary(fdfd,pars,E,false,original);assert(ph1[11]==ph2[11]);
 assert(vary(cdft,pars,E,false,original)[11]==cdft.v[11]);
 pars.fill(0);pars[0]=.01;assert(vary(cdft,pars,E,false,original)[5]==vary(fdfd,pars,E,false,original)[5]);
 assert(sizeof(Record)==200);
 std::cout<<"PASS: physical DVCS conservation, Trento orientation, exact nominal, recomputed shifts, offset preservation, deterministic Gaussian width and record layout\n";
}
'''

CONTROL_TEST='#include "particle_control_core.h"\n#include <iostream>\n#include <cassert>\nusing namespace particle_response;\nint main(){\n Record r{};r.v[5]=5;r.v[6]=.2;r.v[7]=.3;r.v[11]=5;r.v[12]=.12;r.v[13]=1.3;\n double beam=10.604;F pi=particle(5,.12,1.3,MPI0),e=particle(5,.2,.3,ME),k={{0,0,std::sqrt(beam*beam-ME*ME)},beam};\n int failed=0,wrong=0;double maxrel=0;\n auto round5=[](double x){return std::round(x*1e5)/1e5;};\n for(int i=0;i<1000;++i){double ct=2*uniform(i*3)-1,phi=2*PI*uniform(i*3+1);F a=particle(MPI0/2,std::acos(ct),phi,0),b={a.p*(-1),a.e};a=a.boost(pi.p*(1/pi.e));b=b.boost(pi.p*(1/pi.e));\n  double f1=round5(trento((k-e).p,e.p,a.p)),f2=round5(trento((k-e).p,e.p,b.p));double o1=round5(angle(e.p,a.p)*180/PI),o2=round5(angle(e.p,b.p)*180/PI);\n  std::array<double,10> recovered{};double error;bool ok=recover_photons(r,beam,f1,f2,o1,o2,round5(MPI0),recovered,error);\n  if(!ok){++failed;continue;}\n  ControlRecord c{};c.base=r;c.base.v[8]=1;c.base.v[9]=.5;c.base.v[3]=2;c.base.v[4]=1;c.photons=recovered;c.photons[6]=c.photons[7]=1;\n  auto original=derived(c.base.v,beam,false,MPI0);std::array<double,NP> pars{};std::array<double,12> z{};double mass;\n  assert(sizeof(ControlRecord)==280);assert(vary_control(c,pars,beam,original,z,mass)==c.base.v);\n  pars[15]=.01;auto changed=vary_control(c,pars,beam,original,z,mass);\n  assert(std::abs(changed[11]/r.v[11]-1.01)<.001);assert(std::abs(mass/std::sqrt((particle(recovered[0],recovered[1],recovered[2],0)+particle(recovered[3],recovered[4],recovered[5],0)).m2())-1.01)<1e-10);\n  pars.fill(0);pars[12]=.01;z[6]=1;z[9]=-1;auto independent=vary_control(c,pars,beam,original,z,mass);\n  assert(std::abs(mass/std::sqrt((particle(recovered[0],recovered[1],recovered[2],0)+particle(recovered[3],recovered[4],recovered[5],0)).m2())-std::sqrt(1-.0001))<1e-10);\ndouble rel=std::abs(recovered[0]/a.e-1);maxrel=std::max(maxrel,rel);if(rel>.02)++wrong;\n }\n std::cout<<"failed "<<failed<<" wrong "<<wrong<<" maxrelative "<<maxrel<<\'\\n\';assert(failed<10&&wrong<10);\n}\n'

def core(tmp):
    source=tmp/'core.cpp';source.write_text(CORE_TEST);exe=tmp/'core'
    subprocess.run(['c++','-O2','-std=c++17','-I'+str(s.SOURCE),str(source),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)

def control_tests(tmp):
    source=tmp/'control.cpp';source.write_text(CONTROL_TEST);exe=tmp/'control'
    subprocess.run(['c++','-O2','-std=c++17','-I'+str(s.SOURCE),str(source),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
    p='fa18_inb';cache={p+'__mc':dict(probability=1),p+'__background':dict(probability=1)}
    norm={(p,t):dict(signal_scale='2',background_scale='3') for t in range(3)}
    e=s.Engine(None,dict(periods=[p],reference_flows={p+'__data':dict(selected=[1000]*3)}),cache,{},norm,s.TOPS,10)
    e.control=np.array([[[100]*3]*2,[[50]*3]*2])
    e.current_selected={p+'__mc':np.array([[[100]*3]*2,[[50]*3]*2]),p+'__background':np.array([[[20]*3]*2,[[10]*3]*2])}
    b0,a0,r0=e.normalization(p,0,0,0);b1,a1,r1=e.normalization(p,1,0,0)
    assert a1==2*a0 and r1==2 and np.isclose(b1,2*b0)
    assert np.isclose(b0*200+a0*40,1000) and np.isclose(b1*100+a1*20,1000)
    assert e.normalization(p,1,0,0,reference=True)==(2.,3.,1.)
    assert not np.isfinite(e.normalization(p,1,0,0,multiplier=100)[0])
    print('PASS: morphed control transfer, selected-sample signal profiling, raw efficiency separation and invalid-background rejection')


def numerical():
    e=np.linspace(-4,4,401);x=(e[:-1]+e[1:])/2
    d=1000*np.exp(-x*x/2);m=1000*np.exp(-.5*((x-.2)/.8)**2);var=np.ones_like(m)
    good,_=s.shape_score(d,d,var,10000,100,d,var)
    bad,_=s.shape_score(d,m,var,10000,100,d,var)
    assert good<1e-10 and bad>1
    assert np.isnan(s.shape_score(d,m,var,10,100)[0])
    assert s.best_training([[0]*20,[.1]+[0]*19],[dict(train=2,validation_mc=0),dict(train=1,validation_mc=100)])==1
    grid=s.parameter_grid(['CD_FT']);assert set(k//4 for k in grid)=={0,2,4}
    assert len(grid)==12
    p='fa18_inb';n=len(d);spec=dict(periods=[p],variables=[dict(bins=n-2) for _ in s.OBS])
    target=np.r_[0,d[1:-1],0];wrong=np.r_[0,m[1:-1],0]
    data={p+'__data':np.tile(target*2,(3,1,7,1))};hist={};cache={}
    for sample in ('mc','background'):
        hist[p+'__'+sample]=np.tile(wrong,(2,2,3,1,7,1))
        hist[p+'__'+sample][1,0]=np.tile(target,(3,1,7,1))
        cache[p+'__'+sample]=dict(probability=1)
    norm={(p,t):dict(signal_scale=.5,background_scale=.5) for t in range(3)}
    engine=s.Engine(None,spec,cache,data,norm,['CD_FT'],10)
    scores,_=engine.score(hist,[[0]*20,[.01]+[0]*19])
    assert scores[1]['train']<scores[0]['train']
    assert scores[1]['validation_mc']==scores[0]['validation_mc']
    assert scores[1]['training_projections']==3

class ToyEngine:
    def __init__(self,period,goal):
        self.spec=dict(periods=[period],reference_flows={period+'__'+sample:dict(selected=[200]*3) for sample in ('mc','background')})
        self.goal=goal;self.history=[];self.cache={period+'__'+sample:dict(probability=1) for sample in ('mc','background')}
    def evaluate(self,candidates,label,detailed=False):
        p=[tuple([0.]*20)]
        for c in candidates:
            if tuple(c) not in p:p.append(tuple(c))
        scores=[dict(train=(c[0]-self.goal)**2+sum(x*x for x in c[1:]),validation_mc=1,heldout_observables=2,training_projections=3) for c in p]
        selected={key:np.full((len(p),2,3),100) for key in self.cache}
        self.history.extend(dict(call=label,scenario=i,parameters=list(c),**scores[i]) for i,c in enumerate(p))
        return p,scores,[],{},selected

def independent_fit(tmp):
    args=SimpleNamespace(fit_topologies=['CD_FT'],coordinate_rounds=1,refinement_rounds=0)
    def run(item):
        period,goal=item;out=tmp/period;out.mkdir();s.fit_candidate(ToyEngine(period,goal),out,args)
        return json.loads((out/'candidate_parameters.json').read_text())
    with ThreadPoolExecutor(max_workers=2) as pool:a,b=list(pool.map(run,[('sp18_inb',.01),('fa18_inb',.02)]))
    assert a['parameters']['electron_FD_p_sigma_fraction']==.01
    assert b['parameters']['electron_FD_p_sigma_fraction']==.02
    print('PASS: separate period fits recover different injected parameters; validation scores do not select the fit')

def budget_test(tmp):
    worker=tmp/'fake_worker';counter=tmp/'counter.json';s.dump(counter,dict(active=0,peak=0))
    worker.write_text("""#!/usr/bin/env python3
import fcntl,json,sys,time
job=json.load(open(sys.argv[3]))
if job.get('fail'):print('intentional worker failure');sys.exit(3)
def update(delta):
 with open(job['counter'],'r+') as f:
  fcntl.flock(f,fcntl.LOCK_EX);j=json.load(f);j['active']+=delta;j['peak']=max(j['peak'],j['active']);f.seek(0);json.dump(j,f);f.truncate()
update(1);time.sleep(.05);update(-1)
open(sys.argv[4],'wb').write(b'OK')
""");worker.chmod(0o755)
    budget=s.Budget(3);pools=[]
    for i in range(2):
        out=tmp/f'pool{i}';out.mkdir();s.dump(out/'spec.json',{})
        pools.append(s.Pool(worker,out/'spec.json',out,3,budget))
    with ThreadPoolExecutor(max_workers=2) as executor:
        results=list(executor.map(lambda pool:pool.run('evaluate',[dict(counter=str(counter)) for _ in range(4)],'test'),pools))
    j=json.loads(counter.read_text());assert j['peak']<=3 and j['active']==0 and len(results[0])==4
    out=tmp/'failure';out.mkdir();s.dump(out/'spec.json',{})
    failed=s.Budget(2);pool=s.Pool(worker,out/'spec.json',out,2,failed)
    try:pool.run('evaluate',[dict(counter=str(counter),fail=True),dict(counter=str(counter))],'fail')
    except RuntimeError:pass
    else:raise AssertionError('Worker failure was ignored')
    assert not failed.active
    print('PASS: shared worker budget across independent periods and cancellation')

def report_test(tmp):
    spec=dict(periods=['fa18_inb'],variables=[dict(name=n,label=n,bins=40,low=-4,high=4) for n in s.OBS],strata=[dict(axis='inclusive',label='Inclusive')],reference_flows={'fa18_inb__data':dict(selected=[1000]*3)})
    n=42;hist={};cache={};norm={};data={}
    x=np.linspace(-4,4,40);curve=np.r_[0,np.round(50*np.exp(-x*x/2)),0].astype(np.int64)
    for sample in ('mc','background'):
        key='fa18_inb__'+sample;hist[key]=np.tile(curve,(2,2,3,1,7,1));cache[key]=dict(probability=1)
    for top in range(3):norm['fa18_inb',top]=dict(signal_scale='0.8',background_scale='0.2')
    data['fa18_inb__data']=np.tile(curve*2,(3,1,7,1))
    engine=s.Engine(None,spec,cache,data,norm,['CD_FT'],10)
    out=tmp/'report';out.mkdir()
    s.table(out/'single_parameter_ranking.csv',[dict(parameter='electron_FD_p_sigma_fraction',train_score=2,nominal_train=4,heldout_observables=3,nominal_heldout_observables=5)])
    selected={k:np.full((2,2,3),500) for k in hist}
    s.detailed_report(engine,out,[[0]*20,[.01]+[0]*19],[dict(background_multiplier=1),dict(background_multiplier=1.1)],hist,selected)
    assert (out/'REPORT.html').exists() and (out/'plots/fa18_inb/CD_FT.png').exists()
    assert (out/'single_parameter_effects.png').exists()
    with (out/'selection_yields.csv').open() as f:yield_rows=list(s.csv.DictReader(f))
    pair=[r for r in yield_rows if r['topology']=='CD_FT']
    assert np.isclose(float(pair[1]['background_selected'])/float(pair[0]['background_selected']),1.1)
    assert float(pair[1]['signal_selected'])==float(pair[0]['signal_selected'])
    print('PASS: static before/after report, attribution plot and selection/joint summaries (CSV string scales)')
    return engine,out,hist,selected


def profiled_report_test(tmp):
    work=tmp/'profiled';work.mkdir();engine,out,hist,selected=report_test(work)
    engine.control=np.array([[[100]*3]*2,[[50]*3]*2]);engine.control_metadata=dict(unresolved_selected=[0,0],half_records=[1000,1000],mass_outside_input_window=[0,0,4,5])
    selected['fa18_inb__mc'][1]//=2;selected['fa18_inb__background'][1]//=2;engine.current_selected=selected
    s.detailed_report(engine,out,[[0]*20,[.01]+[0]*19],[{},{}],hist,selected)
    with (out/'selection_yields.csv').open() as f:rows=list(s.csv.DictReader(f))
    candidate=[r for r in rows if r['scenario']=='candidate']
    assert all(np.isclose(float(r['predicted_selected']),1000) and float(r['raw_signal_selection_ratio'])==.5 and float(r['control_response_ratio'])==2 for r in candidate)
    assert (out/'control_transfer.csv').exists() and (out/'plots/fa18_inb/CD_FT_shapes.png').exists()
    print('PASS: profiled rate/shape reports retain raw efficiency losses and control mass-window diagnostics')


def scratch_test(tmp):
    a=s.scratch_directory(tmp/'scratch',tmp/'run');b=s.scratch_directory(tmp/'scratch',tmp/'run')
    assert a!=b and a.parent==b.parent==(tmp/'scratch').resolve()
    binary=a/'trial.bin';np.arange(12,dtype=np.int64).tofile(binary)
    mapped={'mc':np.memmap(binary,dtype=np.int64,mode='r',shape=(12,))}
    cache=a/'cache.bin';cache.write_bytes(b'keep')
    s.discard_histograms(mapped)
    assert not binary.exists() and not mapped and cache.read_bytes()==b'keep'
    reference=tmp/'cache_reference';reference.mkdir();(reference/'manifest.json').write_text('{}')
    prior=tmp/'prior';prior.mkdir()
    s.dump(prior/'cache_manifest.json',dict(version=s.VERSION,reference_sha256=s.digest(reference/'manifest.json'),seed=17,
        samples={'sample':dict(file=str(cache),sha256=s.digest(cache),probability=1)}))
    loaded=s.build_cache(None,dict(seed=17),{},reference,prior)
    assert loaded['sample']['path']==str(cache)
    print('PASS: isolated scratch directories, removal of expendable trial binaries, cache preservation and absolute-path reuse')


def recovery_test(tmp):
    from unittest.mock import patch
    case=tmp/'recovery_case';case.mkdir();engine,report,hist,selected=report_test(case)
    run=case/'run';run.mkdir();folder=run/'fa18_inb'
    import shutil
    shutil.move(str(report),folder)
    reference=case/'reference';reference.mkdir();s.dump(reference/'manifest.json',{})
    s.dump(run/'spec.json',engine.spec);s.dump(folder/'spec.json',engine.spec)
    candidate=dict(parameters={s.PARAMS[0]:.01},fit_topologies=['CD_FT'],background_multiplier=1.1);s.dump(folder/'candidate_parameters.json',candidate)
    s.dump(run/'cache_manifest.json',dict(version=s.VERSION,reference_sha256=s.digest(reference/'manifest.json'),samples=engine.cache))
    work=case/'scratch';final=work/'fa18_inb/eval_025_final';final.mkdir(parents=True)
    s.dump(run/'scratch_location.json',dict(directory=str(work)))
    pars=[[0]*20,[.01]+[0]*19]
    for i,sample in enumerate(('mc','background')):
        key='fa18_inb__'+sample;path=final/f'{i}.bin';hist[key].tofile(path)
        s.dump(final/f'{i}.json',dict(period='fa18_inb',sample=sample,detailed=True,parameters=pars))
        s.dump(Path(str(path)+'.json'),dict(shape=hist[key].shape,selected=selected[key].ravel().tolist()))
    with patch.object(s.Engine,'evaluate',side_effect=AssertionError('Recovery must not fit')):
        spec,recovered=s.recover_reports(run,reference,engine.data,engine.norm)
        assert len(recovered)==1;s.render_completed(recovered)
        assert (folder/'report_complete.json').exists()
        assert not list(final.glob('*.bin'))
        _,again=s.recover_reports(run,reference,engine.data,engine.norm);assert again==[]
    (run/'REPORT.html').write_text('index');s.dump(run/'manifest.json',{});(run/'reconstruction_audit.csv').write_text('audit')
    s.compact_results(run,['fa18_inb']);s.compact_results(run,['fa18_inb'])
    assert (run/'results/fa18_inb/REPORT.html').exists()
    print('PASS: saved-trial report recovery without refitting, completion markers, cleanup and repeated publication')


def joint_model_tests(tmp):
    from types import MethodType
    from unittest.mock import patch
    assert len(s.parameter_grid(s.TOPS,'zero'))==15
    assert len(s.parameter_grid(s.TOPS,'fit'))==20
    assert all(k%4!=3 for k in s.parameter_grid(s.TOPS,'zero'))
    p='fa18_inb';norm={(p,t):dict(signal_scale=1,background_scale=1,background_scale_stat=0) for t in range(3)}
    engine=s.Engine(None,dict(periods=[p]),{}, {},norm,s.TOPS,10,.1)
    def fake_score(self,hist,i,m,reference=False):return dict(train=10000*(m-1.1)**2,validation_mc=999 if m>1 else 0,heldout_observables=999,training_projections=9),[]
    engine.score_one=MethodType(fake_score,engine)
    scores,_=engine.score({},[[0]*20]);assert abs(scores[0]['background_multiplier']-1.1)<1e-12
    assert abs(scores[0]['background_prior_penalty']-1)<1e-12
    fixed,_=engine.score({},[[0]*20],profile=False);assert fixed[0]['background_multiplier']==1
    try:s.Engine(None,dict(periods=['sp18_inb','fa18_inb']),{},{},norm,s.TOPS,10)
    except ValueError:pass
    else:raise AssertionError('Cross-period pooling allowed')
    reference=tmp/'compare_reference';reference.mkdir();s.dump(reference/'manifest.json',dict(spec=dict(periods=[p])))
    args=SimpleNamespace(finish_reports=False,scratch_root=tmp/'scratch',workers=3,expected_cache_events=1000,seed=17,
        coordinate_rounds=0,refinement_rounds=0,min_mc=10,background_transfer_uncertainty=.1,fit_topologies=s.TOPS,reuse_cache=tmp/'old_cache')
    output=tmp/'comparison';commands=[]
    def fake_run(cmd,check):
        commands.append(cmd);dest=Path(cmd[cmd.index('--output')+1]);dest.mkdir(parents=True);(dest/'results').mkdir()
        (dest/'results/REPORT.html').write_text('report');s.dump(dest/'manifest.json',{})
        folder=dest/p;folder.mkdir()
        s.table(folder/'final_scores.csv',[dict(scenario='nominal',train=10,validation_mc=10,heldout_observables=10),dict(scenario='candidate',train=5,validation_mc=6,heldout_observables=7)])
        s.dump(folder/'candidate_parameters.json',dict(boundary_parameters=[]))
        s.table(folder/'selection_yields.csv',[dict(topology=t,scenario=a,predicted_selected=100,signal_selected=70) for t in s.TOPS for a in ('nominal','candidate')])
    with patch.object(s.subprocess,'run',fake_run):s.compare_models(args,reference,output)
    assert len(commands)==2
    assert commands[0][commands[0].index('--bias-mode')+1]=='zero'
    assert commands[1][commands[1].index('--bias-mode')+1]=='fit'
    assert commands[1][commands[1].index('--reuse-cache')+1]==str(output/'broadening_only')
    assert '--initial-parameters-from' in commands[1]
    assert (output/'results/model_comparison.csv').exists()
    print('PASS: shared detector laws, no period pooling, broadening-only controls, constrained background profiling and paired cached-model orchestration')


def root_smoke(tmp,reference):
    # Reuse the production-cut test's tiny reconstructed ROOT fixture, not farm data.
    import test_exclusivity_resolution as diagnostic_tests
    diagnostic_tests.root_smoke(tmp,reference)
    saved=json.loads((tmp/'native/spec.json').read_text());saved.update(expected_cache_events=1000,seed=17)
    out=tmp/'response';out.mkdir();s.dump(out/'spec.json',saved)
    # compile_worker needs a reference with saved backend global-cut sources.
    worker=s.compile_worker(out,reference)
    job=dict(period='fa18_inb',sample='mc',files=[str(tmp/'fixture.root')]);s.dump(out/'job.json',job)
    cache=out/'events.bin';subprocess.run([str(worker),'cache',str(out/'spec.json'),str(out/'job.json'),str(cache)],check=True)
    meta=json.loads(Path(str(cache)+'.json').read_text());assert meta['records']==21 and meta['probability']==1
    params=[[0.]*20,[.01]+[0.]*19]
    job=dict(period='fa18_inb',sample='mc',cache=str(cache),parameters=params,detailed=True);s.dump(out/'job.json',job)
    result=out/'hist.bin';subprocess.run([str(worker),'evaluate',str(out/'spec.json'),str(out/'job.json'),str(result)],check=True)
    meta=json.loads(Path(str(result)+'.json').read_text());h=np.fromfile(result,dtype='<i8').reshape(meta['shape'])
    selected=np.array(meta['selected']).reshape(2,2,3);assert np.all(selected[0].sum(axis=0)==2)
    assert not np.array_equal(h[0],h[1]);assert np.all(h>=0)
    for grid in saved['joint_grids']:
        assert np.array_equal(h[0,:,:,grid['offset']:grid['offset']+24].sum(axis=2),h[0,:,:,0])
    print('PASS: actual ROOT cache includes pre-cut events, exact nominal selection, nonzero response and joint reassignment')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root-smoke',action='store_true');parser.add_argument('--reference-study',type=Path)
    args=parser.parse_args()
    if args.root_smoke and not args.reference_study:parser.error('--root-smoke requires --reference-study')
    with tempfile.TemporaryDirectory() as name:
        tmp=Path(name);core(tmp);control_tests(tmp);numerical();independent_fit(tmp);budget_test(tmp);report_test(tmp);profiled_report_test(tmp);scratch_test(tmp);recovery_test(tmp);joint_model_tests(tmp)
        if args.root_smoke:
            # Existing fixture accepts a nominal production source directory;
            # reference backend has those same sources and output cuts below.
            ref=args.reference_study.resolve();fixture=tmp/'fixture_reference';fixture.mkdir()
            spec=json.loads((ref/'spec.json').read_text())
            import shutil
            for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(ref/'backend'/n,fixture/n)
            (fixture/'output/jsons').mkdir(parents=True);s.dump(fixture/'output/jsons/global_cuts_config.json',spec['global_cuts']);s.dump(fixture/'output/jsons/combined_cuts.json',spec['exclusivity'])
            (fixture/'backend').mkdir()
            for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(fixture/n,fixture/'backend'/n)
            root_smoke(tmp,fixture)
    print('PASS: particle-response core and independent-fit checks')
if __name__=='__main__':main()
