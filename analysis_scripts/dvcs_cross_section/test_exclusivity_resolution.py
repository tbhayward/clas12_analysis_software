#!/usr/bin/env python3
"""Numerical/native-core checks; --root-smoke additionally exercises actual ROOT I/O."""
import argparse
import contextlib
import copy
import io
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
import study_exclusivity_resolution as s


def core_test(tmp):
    source=tmp/'core_test.cpp';exe=tmp/'core_test'
    source.write_text(r'''
#include "exclusivity_resolution_core.h"
#include <cassert>
using namespace resolution_study;
int main(){
 std::array<Window,NV> c{};c[3]={true,true,0,.1};c[5]={true,false,-.02,.02};
 std::array<double,NV> v{};v[3]=.2;v[5]=.03;
 unsigned mask=failures(v,c);assert(!leave_one_out(mask,3));assert(!leave_one_out(mask,5));
 v[5]=0;mask=failures(v,c);assert(leave_one_out(mask,3));assert(!leave_one_out(mask,5));
 v[3]=.1;assert(failures(v,c)==0);v[5]=-.02;assert(failures(v,c)==0);
 v[5]=std::numeric_limits<double>::quiet_NaN();mask=failures(v,c);assert(leave_one_out(mask,5));assert(!leave_one_out(mask,3));
 Histogram h(10,0,1);h.add(-.1,false);h.add(0,true);h.add(.1,true);h.add(1,false);h.add(std::numeric_limits<double>::infinity(),false);
 assert(h.counts[0]==1&&h.counts[1]==1&&h.counts[2]==1&&h.counts[11]==1&&h.nonfinite==1&&h.pass_own_cut==2);
 Stratum group{"p1_theta",45,50};assert(group.contains(45));assert(!group.contains(50));assert(!group.contains(std::numeric_limits<double>::quiet_NaN()));
}
''')
    subprocess.run(['c++','-std=c++17','-O2','-I'+str(s.SOURCE),str(source),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)


def numerical_tests():
    assert s.saved_window(dict(mode='upper_quantile',cut_low=0,cut_high=.2))==[None,.2]
    assert s.saved_window(dict(mean=1,std=2))==[-5,7]
    for bad in [dict(mode='upper_quantile',cut_high=float('nan')),dict(mean=0,std=0)]:
        try:s.saved_window(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid cuts accepted')
    edges=np.linspace(-5,5,1001);centers=(edges[:-1]+edges[1:])/2
    def gaussian(mu,sd):return np.r_[0,np.round(1e5*np.exp(-.5*((centers-mu)/sd)**2)),0].astype(np.uint64)
    mc=gaussian(0,.4);data=gaussian(.12,.6)
    m=s.shape_metrics(mc,edges,0,100);d=s.shape_metrics(data,edges,0,100)
    assert abs(d['core_spread']/m['core_spread']-1.5)<.002
    assert abs(d['median']-.12)<.001
    assert abs(s.outside_fraction(mc,edges,m['q025'],m['q975'])-.05)<1e-9
    assert s.outside_fraction(data,edges,m['q025'],m['q975'])>.19
    clipped=mc.copy();clipped[0]=mc.sum()
    assert abs(s.shape_metrics(clipped,edges,0,0)['in_range_fraction']-.5)<1e-9
    assert math.isnan(s.shape_metrics(np.zeros(1002),edges,0,0)['median'])
    variables,axes,groups,index=s.configuration()
    for g in groups[1:]:
        if g['axis']=='p1_theta':assert g['high']>g['low']
    assert variables[2]['label'].endswith('(deg)') and variables[4]['label'].endswith('(GeV)')


def synthetic_report(tmp):
    variables,axes,strata,index=s.configuration()
    spec=dict(version=s.VERSION,periods=['fa18_inb'],period_codes={'fa18_inb':'Fa18_Inb'},
        variables=variables,axes=axes,strata=strata,strata_axis_index=index,
        global_cuts=dict(enable_pTmiss_cut=False),exclusivity={})
    for t in s.TOPOLOGIES:spec['exclusivity']['DVCS_Fa18_Inb_'+t]={sample:{'Mx2':dict(cut_low=-.02,cut_high=.02)} for sample in ('data','mc')}
    arrays={};maxbins=max(v['bins'] for v in variables)+2
    for sample in ('data','mc'):
        key='fa18_inb__'+sample
        arrays[key]=np.zeros((3,len(strata),7,maxbins),dtype=np.uint64)
        arrays[key+'__stats']=np.zeros((3,len(strata),7,2),dtype=np.uint64)
        for vi,v in enumerate(variables):
            e=np.linspace(v['low'],v['high'],v['bins']+1);c=(e[:-1]+e[1:])/2
            mu=(v['low']+v['high'])/2;sd=(v['high']-v['low'])/(12 if sample=='data' else 18)
            h=np.r_[0,np.round(1000*np.exp(-.5*((c-mu)/sd)**2)),0].astype(np.uint64)
            arrays[key][:,:,vi,:len(h)]=h
    output=tmp/'report';output.mkdir();np.savez_compressed(output/'histograms.npz',**arrays)
    manifest=dict(spec=spec,flows={k:dict(entries=100,selected=[10,20,30]) for k in ['fa18_inb__data','fa18_inb__mc']},histogram_sha256=s.digest(output/'histograms.npz'))
    s.write_json(output/'manifest.json',manifest)
    s.render(output,arrays,manifest,300,60)
    assert (output/'REPORT.html').is_file() and (output/'summary.csv').is_file()
    page=(output/'REPORT.html').read_text();assert '__PAYLOAD__' not in page
    payload=json.loads(page.split('<script id="payload" type="application/json">')[1].split('</script>')[0])
    assert payload['detail']['fa18_inb|FD_FD|0|5']['summary']['width_ratio']>1.49
    replay=tmp/'replay'
    s.main(['--reuse-histograms',str(output),'--output',str(replay),'--dpi','60'])
    assert (replay/'summary.csv').read_bytes()==(output/'summary.csv').read_bytes()
    for args in [ ['--reuse-histograms',str(output),'--output',str(replay)],
                  ['--reuse-histograms',str(output),'--output',str(tmp/'bad'),'--config','unused.json'] ]:
        with contextlib.redirect_stderr(io.StringIO()):
            try:s.main(args)
            except SystemExit as e:assert e.code!=0
            else:raise AssertionError('Overwrite/reuse guard failed')
    (output/'histograms.npz').write_bytes(b'broken')
    with contextlib.redirect_stderr(io.StringIO()):
        try:s.main(['--reuse-histograms',str(output),'--output',str(tmp/'bad_hash')])
        except SystemExit as e:assert e.code!=0
        else:raise AssertionError('Corrupted cache accepted')


def driver_test(tmp):
    from unittest.mock import patch
    variables,axes,strata,index=s.configuration()
    paths={name:tmp/file for name,file in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]}
    for path in paths.values():path.write_text('// test source snapshot\n')
    spec=dict(version=s.VERSION,periods=['fa18_inb','fa18_out'],variables=variables,strata=strata,
              input_hashes={k:s.digest(v) for k,v in paths.items()})
    files={}
    for period in spec['periods']:
        for sample in ('data','mc'):
            files[period,sample]=[]
            for i in range(3):
                path=tmp/f'{period}_{sample}_{i}.root';path.write_bytes(b'fixture')
                files[period,sample].append(path)
    def fake_compile(output,paths):
        build=output/'backend';build.mkdir()
        for path in paths.values():shutil.copy2(path,build/path.name)
        worker=build/'worker';worker.write_text("""#!/usr/bin/env python3
import json,sys
spec=json.load(open(sys.argv[1]));job=json.load(open(sys.argv[2]));n=len(job['files'])
if any('bad.root' in f for f in job['files']):
 print('intentional worker failure');sys.exit(7)
h=[]
for t in range(3):
 for g in spec['strata']:
  for v in spec['variables']:
   h.append(dict(counts=[0,n]+[0]*v['bins'],nonfinite=0,pass_own_cut=n))
json.dump(dict(histograms=h,entries=n*10,global_pass=[n]*3,selected=[n]*3,files=[dict(path=f,entries=10) for f in job['files']]),open(sys.argv[3],'w'))
""");worker.chmod(0o755);return worker
    out=tmp/'driver';out.mkdir()
    with patch.object(s,'compile_worker',fake_compile):arrays,manifest=s.scan(out,spec,paths,files,2)
    for key in manifest['flows']:
        assert manifest['flows'][key]['entries']==30
        assert np.all(arrays[key][:,:,:,1]==3)
    bad=tmp/'bad.root';bad.write_bytes(b'bad');failed=tmp/'failed';failed.mkdir()
    with patch.object(s,'compile_worker',fake_compile):
        try:s.scan(failed,spec,paths,{('fa18_inb','data'):[bad],('fa18_out','mc'):files['fa18_out','mc']},1)
        except RuntimeError as e:assert 'intentional worker failure' in str(e)
        else:raise AssertionError('Worker failure was ignored')
    assert not (failed/'manifest.json').exists()


def root_smoke(tmp,production):
    """Actual worker and production global functions; deterministic tiny ROOT input."""
    import shlex
    spec,paths=s.build_spec(production,['fa18_inb','fa18_out'])
    cfg=spec['global_cuts'];cfg.update(enable_auxiliary_fiducial_cuts=False,
        enable_analysis_phase_space_cuts=False,enable_pTmiss_cut=True,pTmiss_max=.15,
        excluded_runs=[999],enable_topology_filter=False,enable_dvcsgen_ycol_cut=False)
    for k in list(cfg):
        if k.endswith('_sector_filter'):cfg[k]=False
    for p in spec['periods']:
        for t in s.TOPOLOGIES:
            spec['exclusivity'][f'DVCS_{spec["period_codes"][p]}_{t}']={sample:{'Mx2':dict(cut_low=-.02,cut_high=.02),'pTmiss':dict(mode='upper_quantile',cut_high=.1 if sample=='mc' else .2)} for sample in ('data','mc')}
    out=tmp/'native';out.mkdir();worker=s.compile_worker(out,paths)
    generator=tmp/'generate.cpp';exe=tmp/'generate'
    generator.write_text(r'''
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <string>
int main(int argc,char**argv){TFile f(argv[1],"RECREATE");TTree t("PhysicsEvents","");std::map<std::string,double>v;
for(auto n:{"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pTmiss","theta","theta_gamma_gamma","Emiss2","Mx2","Mx2_2","runnum"}){v[n]=0;t.Branch(n,&v[n],(std::string(n)+"/D").c_str());}
v["x"]=.2;v["Q2"]=2;v["t1"]=-.13;v["e_p"]=5;v["e_theta"]=.25;v["p1_p"]=.6;v["p1_theta"]=.8;v["p1_phi"]=3.141592653589793;v["p2_p"]=4;v["p2_theta"]=.2;v["open_angle_ep2"]=10;v["theta"]=3;v["theta_gamma_gamma"]=1;v["Mx2_2"]=.88;
for(int top=0;top<3;++top){v["detector1"]=top==0?1:2;v["detector2"]=top==2?0:1;
for(int i=0;i<7;++i){v["runnum"]=i==6?999:100;v["pTmiss"]=(i==1||i==3)?.18:(i==5?.12:.05);v["Mx2"]=(i==2||i==3)?.04:(i==4?-.13:0);t.Fill();}}
t.Write();}
''')
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
    subprocess.run(['c++',*flags,'-std=c++17',str(generator),'-o',str(exe),*libs],check=True)
    root=tmp/'fixture.root';subprocess.run([str(exe),str(root)],check=True)
    s.write_json(out/'spec.json',spec)
    # Six nonblacklisted events/topology: pT = .05,.18,.05,.18,.05,.12;
    # Mx = 0,0,.04,.04,-.13,0. N-1 Mx: MC 3, DATA 4;
    # N-1 pT: both 3; non-cut theta after all cuts: MC 1, DATA 2.
    for period in spec['periods']:
        for sample in ('data','mc'):
            jp=out/'job.json';result=out/'result.json'
            s.write_json(jp,dict(period=period,sample=sample,files=[str(root)]))
            subprocess.run([str(worker),str(out/'spec.json'),str(jp),str(result)],check=True)
            payload=json.loads(result.read_text());ng=len(spec['strata'])
            assert payload['entries']==21
            assert payload['global_pass']==[4]*3
            assert payload['selected']==([1]*3 if sample=='mc' else [2]*3)
            for t in range(3):
                def hist(v):return payload['histograms'][(t*ng)*7+v]
                assert sum(hist(5)['counts'])==(3 if sample=='mc' else 4)
                assert hist(5)['counts'][0]==1
                assert sum(hist(3)['counts'])==3 # optional GLOBAL pT pre-cut also omitted
                assert sum(hist(1)['counts'])==(1 if sample=='mc' else 2)
                assert sum(hist(0)['counts'])==(1 if sample=='mc' else 2) # Delta_phi fallback
            assert all(f['delta_phi_fallback'] for f in payload['files'])
    print('PASS: native ROOT I/O, both periods/samples, all topologies, global/own cut omission, blacklist, overflow and phi fallback')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root-smoke',action='store_true')
    parser.add_argument('--production-run',type=Path)
    args=parser.parse_args()
    if args.root_smoke and not args.production_run:parser.error('--root-smoke requires --production-run')
    with tempfile.TemporaryDirectory() as name:
        tmp=Path(name);core_test(tmp);numerical_tests();driver_test(tmp);synthetic_report(tmp)
        print('PASS: native leave-one-out/cut/histogram logic, widths/shifts/tails, empty/range flags, reports, worker orchestration/failure handling, replay and overwrite/checksum guards')
        if args.root_smoke:root_smoke(tmp,args.production_run)
    return 0
if __name__=='__main__':raise SystemExit(main())
