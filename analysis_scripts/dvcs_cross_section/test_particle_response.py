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
 assert(sizeof(Record)==200);
 std::cout<<"PASS: physical DVCS conservation, Trento orientation, exact nominal, recomputed shifts, offset preservation, deterministic Gaussian width and record layout\n";
}
'''

def core(tmp):
    source=tmp/'core.cpp';source.write_text(CORE_TEST);exe=tmp/'core'
    subprocess.run(['c++','-O2','-std=c++17','-I'+str(s.SOURCE),str(source),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)

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
    with ThreadPoolExecutor(max_workers=2) as pool:a,b=list(pool.map(run,[('sp18_inb',.006),('fa18_inb',.012)]))
    assert a['parameters']['electron_FD_p_sigma_fraction']==.006
    assert b['parameters']['electron_FD_p_sigma_fraction']==.012
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
    for top in range(3):norm['fa18_inb',top]=dict(signal_scale=.8,background_scale=.2)
    data['fa18_inb__data']=np.tile(curve*2,(3,1,7,1))
    engine=s.Engine(None,spec,cache,data,norm,['CD_FT'],10)
    out=tmp/'report';out.mkdir()
    s.table(out/'single_parameter_ranking.csv',[dict(parameter='electron_FD_p_sigma_fraction',train_score=2,nominal_train=4,heldout_observables=3,nominal_heldout_observables=5)])
    selected={k:np.full((2,2,3),500) for k in hist}
    s.detailed_report(engine,out,[[0]*20,[.01]+[0]*19],[],hist,selected)
    assert (out/'REPORT.html').exists() and (out/'plots/fa18_inb/CD_FT.png').exists()
    assert (out/'single_parameter_effects.png').exists()
    print('PASS: static before/after report, attribution plot and selection/joint summaries')


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
        tmp=Path(name);core(tmp);numerical();independent_fit(tmp);budget_test(tmp);report_test(tmp)
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
