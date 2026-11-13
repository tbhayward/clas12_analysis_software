#!/usr/bin/env python3
"""Native response-law checks, training-only selection and isolated worker smoke."""
import argparse,importlib.util,json,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_cd_proton_angle_response as study


def provenance():
    names=('particle_response_cut_arrays.h','particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h')
    current={n:base.digest(base.SOURCE/n) for n in names}
    meta=dict(source_hashes=dict(current,**{'particle_response_cut_arrays.h':'older helper without individual-photon export'}))
    fp=dict(export_sources=current.copy())
    study.check_source_provenance(meta,fp)
    for source,name in [(meta['source_hashes'],'particle_response_core.h'),(meta['source_hashes'],'particle_response_worker.cpp'),(fp['export_sources'],'particle_response_cut_arrays.h'),(fp['export_sources'],'particle_control_core.h')]:
        old=source[name];source[name]='changed'
        try:study.check_source_provenance(meta,fp)
        except ValueError:pass
        else:raise AssertionError('Changed source was incorrectly accepted: '+name)
        finally:source[name]=old
    print('PASS: historical joint-fit export helper accepted only with matching frozen-export provenance; changed evaluators/exporters rejected')


def native_core(tmp):
    isolated=tmp/'isolated';study.angle_sources(base.SOURCE,isolated)
    cpp=tmp/'core.cpp';cpp.write_text(r'''
#include "particle_control_core.h"
#include <cassert>
#include <iomanip>
#include <iostream>
using namespace particle_response;
int main(){Record r{};r.v[3]=2;r.v[4]=1;r.v[5]=6;r.v[6]=.2;r.v[7]=.1;r.v[8]=.7;r.v[9]=40*PI/180;r.v[10]=2;r.v[11]=3;r.v[12]=.3;r.v[13]=4;r.key=23;
std::array<double,NP> p{};p[8]=.02;p[11]=-.01;p[15]=.005;std::array<double,9> noise{};noise[3]=1.;auto original=derived(r.v,10.604,false);auto v=vary(r,p,10.604,false,original,&noise);
ControlRecord c{};c.base=r;c.photons={2,.2,.5,1,.3,1,1,0,.135,1};std::array<double,12> z{};z[3]=1.;double mass;auto co=derived(r.v,10.604,false,MPI0);auto cv=vary_control(c,p,10.604,co,z,mass);
std::cout<<std::setprecision(17);for(double x:v)std::cout<<x<<' ';for(double x:cv)std::cout<<x<<' ';std::cout<<mass<<'\n';
#ifdef ANGLE_TEST
assert(cd_angle_sigma(35*PI/180,p)==0);p[20]=.04;
assert(std::abs(cd_angle_sigma(45*PI/180,p)-.02)<1e-14);assert(cd_angle_sigma(60*PI/180,p)==0);
auto changed=vary(r,p,10.604,false,original,&noise);auto cc=vary_control(c,p,10.604,co,z,mass);
double extra=.04*.75,expected=.7*(1-.01+std::hypot(.02,extra));assert(std::abs(changed[8]-expected)<1e-14);assert(changed[8]==cc[8]);
r.v[9]=60*PI/180;original=derived(r.v,10.604,false);auto high=vary(r,p,10.604,false,original,&noise);assert(high[8]==v[8]);
r.v[3]=1;r.v[9]=.5;original=derived(r.v,10.604,false);auto fd=vary(r,p,10.604,false,original,&noise);p[20]=0;auto fd0=vary(r,p,10.604,false,original,&noise);assert(fd==fd0);
// Photon-role exchange leaves their summed four-vector unchanged when their noises exchange too.
c.photons={2,.2,.5,1,.3,1,1,0,.135,1};p[12]=.01;p[16]=.03;z[6]=.4;z[9]=-.2;auto a=vary_control(c,p,10.604,co,z,mass);
for(int j=0;j<3;++j)std::swap(c.photons[j],c.photons[j+3]);std::swap(c.photons[6],c.photons[7]);for(int j=0;j<3;++j)std::swap(z[6+j],z[9+j]);auto b=vary_control(c,p,10.604,co,z,mass);for(int j=11;j<=13;++j)assert(std::abs(a[j]-b[j])<1e-12);
#endif
}
''')
    outputs=[]
    for label,source,flags in [('old',base.SOURCE,[]),('angle',isolated,['-DANGLE_TEST'])]:
        exe=tmp/label;subprocess.run(['c++','-std=c++17','-O2',*flags,'-I'+str(source),str(cpp),'-o',str(exe)],check=True);outputs.append(subprocess.check_output([str(exe)]))
    assert outputs[0]==outputs[1],'Zero extra parameter must preserve every old output bit for bit'
    print('PASS: exact zero-extension replay, angle taper, quadrature widths, shared one/two-photon proton response, untouched FD response and photon-role symmetry')


def selector():
    start=[.01]*20;parameters=[tuple([0.]*21),tuple(start+[0]),tuple(start+[.02]),tuple(start+[.04])]
    scores=[{'train':0,'validation_mc':0},{'train':3,'validation_mc':1},{'train':2,'validation_mc':100},{'train':2,'validation_mc':0}]
    assert study.pick(parameters,scores,start)==2
    scores[2]['validation_mc']=-1e9;assert study.pick(parameters,scores,start)==2
    scores[2]['train']=4;scores[3]['train']=4;assert study.pick(parameters,scores,start)==1
    print('PASS: training-only selection, zero option, lower-width tie break; validation cannot choose the amplitude')


def plot_smoke(tmp):
    from types import SimpleNamespace
    rng=np.random.default_rng(17);a=np.zeros((1600,28));a[:,0]=np.repeat([1,2],800);a[:,8]=1;a[:,9]=.2;a[:,10]=1;a[:,11]=np.arange(1600)%2;a[:,14]=6;a[:,15]=.2;a[:,16]=.7;a[:,17]=rng.uniform(.6,1.2,1600);a[:,20]=3;a[:,21]=.2;a[:,23]=1;a[:,24]=1;a[:,25]=.3;a[:,27]=1
    for j in range(1,8):a[:,j]=rng.normal(0,.15,len(a))
    variables=[dict(name=n,label=n,bins=80,low=-1,high=1) for n in base.OBS]
    spec=dict(variables=variables,exclusivity={prefix+'_Fa18_Inb_'+t:dict(data={}) for prefix in ('DVCS','eppi0') for t in base.TOPS})
    class Pool:
        def run(self,mode,jobs,label):
            result=[]
            for i,job in enumerate(jobs):
                path=tmp/(str(i)+'.bin');a.tofile(path);base.dump(Path(str(path)+'.json'),dict(shape=list(a.shape)));result.append((i,job,path))
            return result
    cache={'fa18_inb__'+s:dict(path='unused',probability=1.) for s in ('mc','background','pi0_mc')}
    engine=SimpleNamespace(spec=spec,cache=cache,min_mc=10,normalization=lambda *args:(1.,.1,1.))
    (tmp/'fa18_inb').mkdir()
    rows,links=study.low_high('fa18_inb',engine,None,Pool(),{'frozen':[0.]*21,'angle':[0.]*20+[.02]},{'frozen':0,'angle':1},{'frozen':1.,'angle':1.},{'data':a,'pi0_data':a},tmp)
    assert len(rows)==96 and len(links)==4 and all((tmp/x).is_file() for x in links)
    assert (tmp/'fa18_inb/photon_order_audit.csv').is_file()
    paired={}
    for r in rows:paired.setdefault((r['channel'],r['topology'],r['side'],r['variable'],r['mc_half']),[]).append(r['shape_score'])
    assert all(np.allclose(v[0],v[1],equal_nan=True) for v in paired.values())
    print('PASS: complete low/high plotting, paired frozen-reference scores, both channels/topologies and matched-selection photon-order audit')


def root_smoke(tmp):
    # Reuse existing native I/O fixture, then compare 20- and 21-parameter workers on it.
    import test_frozen_response_exports as fixture
    fixture.native(tmp)
    source=base.SOURCE;out=tmp/'angle_build';out.mkdir();work=tmp/'angle_work';work.mkdir()
    worker=study.compile_angle_worker(out,tmp/'reference',work)
    spec=tmp/'spec.json';job=json.loads((tmp/'job.json').read_text());job['parameters']=[0.]*21
    job.update(export_individual_photons=True,export_global_only=True);base.dump(tmp/'angle_job.json',job)
    oldjob=dict(job,parameters=[0.]*20);base.dump(tmp/'old_job.json',oldjob)
    for exe,j,target in [(tmp/'worker/backend/particle_response_worker',tmp/'old_job.json',tmp/'old_export.bin'),(worker,tmp/'angle_job.json',tmp/'angle_export.bin')]:
        subprocess.run([str(exe),'cut-arrays',str(spec),str(j),str(target)],check=True)
    assert (tmp/'old_export.bin').read_bytes()==(tmp/'angle_export.bin').read_bytes()
    job['parameters']=[[0.]*21,[0.]*20+[.04]];base.dump(tmp/'angle_job.json',job)
    subprocess.run([str(worker),'evaluate',str(spec),str(tmp/'angle_job.json'),str(tmp/'angle_hist.bin')],check=True)
    m=json.loads((tmp/'angle_hist.bin.json').read_text());assert m['hist_shape'][0]==2
    job['sample']='mc';base.dump(tmp/'angle_job.json',job)
    # The signal worker uses a different cache layout; build deterministic native records.
    dt=np.dtype([('v','<f8',(24,)),('key','<u8')]);a=np.zeros(20,dtype=dt);a['v'][:,0]=.2;a['v'][:,1]=2.2;a['v'][:,2]=-.2;a['v'][:,3]=2;a['v'][:,4]=1;a['v'][:,5]=6;a['v'][:,6]=.2;a['v'][:,8]=.7;a['v'][:,9]=.7;a['v'][:,11]=3;a['v'][:,12]=.2;a['v'][:,23]=5032;a['key']=np.arange(20)*4;a.tofile(tmp/'signal.bin');job['cache']=str(tmp/'signal.bin');base.dump(tmp/'angle_job.json',job)
    ss=json.loads(spec.read_text());ss['exclusivity'].update({k.replace('eppi0_','DVCS_'):v for k,v in list(ss['exclusivity'].items())});ss['strata']=[dict(name='inclusive',kind='inclusive')];base.dump(spec,ss)
    subprocess.run([str(worker),'evaluate',str(spec),str(tmp/'angle_job.json'),str(tmp/'signal_hist.bin')],check=True)
    print('PASS: actual ROOT worker compilation, zero-parameter export replay and both control/signal extended evaluator paths')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root-smoke',action='store_true');args=p.parse_args();provenance();selector()
    with tempfile.TemporaryDirectory() as folder:
        native_core(Path(folder))
    with tempfile.TemporaryDirectory() as folder:plot_smoke(Path(folder))
    if args.root_smoke:
        with tempfile.TemporaryDirectory() as folder:root_smoke(Path(folder))
