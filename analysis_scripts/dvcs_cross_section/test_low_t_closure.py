#!/usr/bin/env python3
"""Algebra/grid regressions; --root-smoke also exercises real paired ROOT input."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
import time
import numpy as np

SOURCE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('low_t',SOURCE/'study_low_t_closure.py')
study=importlib.util.module_from_spec(sp);sp.loader.exec_module(study)


def fixture_bins():
    return [dict(id=str(i),valid=1,bounds=[.204,.268,1.912,2.51,*b]) for i,b in enumerate([
        [.11,.15,0,30],[.15,.25,0,15],[.15,.25,15,30]])]


def unit():
    rows=fixture_bins();merged=study.merged_bins(rows);control=study.coverage_control_bins(rows)
    assert len(merged)==1 and merged[0]['bounds'][4:]==[.11,.25,0,30]
    assert len(control)==2 and [b['bounds'][4:6] for b in control]==[[.11,.15],[.15,.25]]
    for name,items in [('nominal',rows),('control',control),('merged',merged)]:
        grid=study.make_grid(items,name);assert max(grid['lookup'])==len(items)-1
    try:study.make_grid(rows+[rows[0]],'overlap')
    except ValueError:pass
    else:raise AssertionError('Overlapping bins accepted')
    spec=dict(schemes=[study.make_grid(rows,'nominal')],scenarios=study.SCENARIOS)
    a=study.empty(spec)
    for key in ['train_gen']:a['0_'+key][:]=1000
    for key in ['train_reco','train_reco2']:a['0_'+key][:]=800
    a['0_train_diag'][:]=600
    a['0_test_gen'][:]=1000;a['0_test_gen2'][:]=1000
    a['0_test_reco'][:]=800;a['0_test_reco2'][:]=800;a['0_test_diag'][:]=600
    acceptance,corrected,ratio,error=study.calculate(a,0)
    assert np.allclose(acceptance,.8) and np.allclose(ratio,1)
    assert np.all(np.isfinite(error)) and np.all(error>0)
    # The observed 980 includes 100 outside-bin events as well as diagonal and
    # other-bin contributions; the denominator remains the known 1300 truth.
    a['0_test_gen'][1]=1300;a['0_test_reco'][1]=980
    assert np.allclose(study.calculate(a,0)[2][1],980/(.8*1300))
    print('PASS: non-overlapping phi coarsening, separate control, closure algebra, outside-bin contribution, covariance propagation')


def definition(text,name,struct=False):
    pattern=(r'^struct '+name+r'\s*\{') if struct else (r'^(?:static )?[\w:<>, ]+\b'+name+r'\(')
    match=re.search(pattern,text,re.M)
    assert match,name
    brace=text.index('{',match.start());depth=1;i=brace+1
    while depth:
        depth+=(text[i]=='{')-(text[i]=='}');i+=1
    return text[match.start():i]+(';' if struct else '')+'\n'


def native_and_replay(repo):
    """Exercise actual C++ helpers without loading the local ROOT interpreter."""
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);worker=(SOURCE/'low_t_closure_worker.cpp').read_text();global_code=(repo/'global_cuts.cpp').read_text()
        code='''#include "global_cuts.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <limits>
#include <cassert>
#include <fstream>
#include <stdexcept>
#include <cstdint>
using json=nlohmann::json;
constexpr double DEG=180./3.14159265358979323846;
static void validate_cfg_or_fatal(const GlobalCutConfig&){}
static void fatal_global(const std::string& s){throw std::runtime_error(s);}
'''
        for name in ['config','mix','path_hash']:code+=definition(worker,name)
        for name in ['Grid','Scheme','Window']:code+=definition(worker,name,True)
        for name in ['passes_reconstructed_electron_momentum_cut','dvcs_minus_t_min','passes_analysis_phase_space_cuts','passes_generated_channel_cuts','passes_proton_momentum_cut']:
            code+=definition(global_code,name)
        fields=[n for n in re.findall(r'LOAD\((\w+)\)',worker) if n!='name']
        code+='int main(int argc,char**argv){GlobalCutConfig c;c.enable_auxiliary_fiducial_cuts=true;json j;\n'
        code+=''.join(f'j["{n}"]=c.{n};\n' for n in fields)
        code+='''std::ofstream out(argv[1]);out<<j;out.close();auto cfg=config(j);
 assert(!passes_reconstructed_electron_momentum_cut(2.,cfg));
 assert(passes_analysis_phase_space_cuts(10.604,.23,2.2,-.2,5.,cfg));
 assert(!passes_analysis_phase_space_cuts(10.604,.23,2.2,-.01,5.,cfg));
 assert(!passes_analysis_phase_space_cuts(10.604,.23,2.2,-.2,1.,cfg));
 assert(passes_generated_channel_cuts(.23,2.2,cfg));
 assert(passes_proton_momentum_cut(1,"fa18_inb",.45,cfg));
 assert(!passes_proton_momentum_cut(1,"fa18_out",.45,cfg));
 auto upper=json{{"mode","upper_quantile"},{"cut_high",.2},{"cut_low",0.}};
 assert(Window(upper).pass(-.01));assert(!Window(upper).pass(.21));
 auto symmetric=json{{"mode","quantile_window"},{"cut_low",-.02},{"cut_high",.015599999999999996}};
 assert(Window(symmetric).pass(.015599999999999996));
 assert(!Window(symmetric).pass(std::nextafter(.015599999999999996,1.)));
 assert(!Window(symmetric).pass(std::numeric_limits<double>::quiet_NaN()));
 json grid={{"edges",{{.2,.3},{2.,3.},{.11,.15,.25},{0.,15.,30.}}},
            {"lookup",{0,1,2,3}},{"bins",{0,1,2,3}}};
 Grid g(grid);assert(g.find(.23,2.2,.12,7./DEG)==0);
 assert(g.find(.23,2.2,.20,22./DEG)==3);assert(g.find(.23,2.2,.09,7./DEG)==-1);
 Scheme s(grid,2,2);std::vector<double>w={1.,3.};
 s.add(0,0,true,1,2.,w);s.add(-1,0,true,1,2.,w);
 assert(s.train_gen[0]==1);assert(s.train_reco[0]==4);
 assert(s.train_diag[0]==2);assert(s.train_reco2[0]==8);
 s.add(0,0,false,1,2.,w);
 assert(s.test_gen[8]==3);assert(s.test_reco[8]==6);
 assert(s.test_diag[8]==18);
 int train=0;for(uint64_t i=0;i<1000;++i)train+=(mix(i^path_hash("fixture"))&1)==0;
 assert(train>400&&train<600);return 0;}
'''
        src=root/'native.cpp';src.write_text(code);exe=root/'native';cfg=root/'global.json'
        flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
        subprocess.run(['c++',*flags,'-std=c++17','-I'+str(repo),str(src),'-o',str(exe)],check=True)
        subprocess.run([str(exe),str(cfg)],check=True)
        # Re-render a known saved aggregate for both periods without a ROOT read.
        run=root/'production';(run/'output/jsons').mkdir(parents=True)
        shutil.copy2(cfg,run/'output/jsons/global_cuts_config.json')
        for name in ['global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h']:shutil.copy2(repo/name,run/name)
        cuts={};model={'mc':{'ep->epg':{}}}
        for p,(tag,code) in study.PERIODS.items():
            for topo in ['FD_FD','CD_FD','CD_FT']:cuts[f'DVCS_{code}_{topo}']={'mc':{'Mx2':dict(mode='quantile_window',cut_low=-.02,cut_high=.02)}}
            model['mc']['ep->epg'][tag]={region:dict(reference_factor=.85) for region in ['FT',*['S'+str(i) for i in range(1,7)]]}
        (run/'output/jsons/combined_cuts.json').write_text(json.dumps(cuts))
        d=run/'output/dvcs_current_dependence/calibration';d.mkdir(parents=True);(d/'current_response_model.json').write_text(json.dumps(model))
        d=run/'output/csvs';d.mkdir()
        rows=[{'bin index':r['id'],'valid bin':1,**dict(zip(study.FIELDS,r['bounds']))} for r in fixture_bins()]
        with (d/'dvcs_pass2_analysis.csv').open('w',newline='') as handle:
            w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        spec,_=study.build_spec(run,20261111);source=root/'saved';source.mkdir()
        hashes={n:study.digest(SOURCE/n) for n in ['study_low_t_closure.py','low_t_closure_worker.cpp']}
        hashes.update({n:study.digest(repo/n) for n in ['global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h']})
        fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        manifest=dict(version=study.VERSION,input_hashes=spec['input_hashes'],seed=20261111,source_hashes=hashes,
            event_source_hashes={k:v for k,v in hashes.items() if k!='study_low_t_closure.py'},spec_fingerprint=fingerprint,
            files={p:[] for p in study.PERIODS},file_inventory={p:[] for p in study.PERIODS},max_files_per_period=1)
        (source/'manifest.json').write_text(json.dumps(manifest));(source/'cutflow.json').write_text(json.dumps({p:{} for p in study.PERIODS}))
        for period in study.PERIODS:
            a=study.empty(spec)
            for si in range(3):
                a[f'{si}_train_gen'][:]=1000;a[f'{si}_train_reco'][:]=800;a[f'{si}_train_reco2'][:]=800;a[f'{si}_train_diag'][:]=600
                a[f'{si}_test_gen'][:]=1000;a[f'{si}_test_gen2'][:]=1000;a[f'{si}_test_reco'][:]=800;a[f'{si}_test_reco2'][:]=800;a[f'{si}_test_diag'][:]=600
            control={tuple(b['bounds']):i for i,b in enumerate(spec['schemes'][1]['bins'])}
            for i,b in enumerate(spec['schemes'][2]['bins']):
                if abs(b['bounds'][4]-.11)>1e-9:continue
                bounds=list(b['bounds']);bounds[4:6]=[.11,.15];first=control[tuple(bounds)]
                bounds[4:6]=[.15,.25];second=control[tuple(bounds)]
                for key in list(a):
                    if key.startswith('2_'):
                        a[key][...,i]=a['1_'+key[2:]][...,first]+a['1_'+key[2:]][...,second]
            np.savez_compressed(source/f'aggregates_{period}.npz',**a)
        out=root/'replay';common=['--production-run',str(run),'--reuse-aggregates',str(source),'--dpi','40']
        assert study.main(common+['--output',str(out)])==0
        assert (out/'REPORT.md').is_file() and (out/'closure_nominal_vs_merged.png').is_file()
        # Exercise the real batching/subprocess driver with a lightweight worker.
        # The actual C++ numerical helpers were checked above; no ROOT runtime is
        # required to verify orchestration and cancellation.
        dummy_input=root/'input'
        for period in study.PERIODS:
            directory=dummy_input/'mc'/period/'dvcs_mc/50nA/parts';directory.mkdir(parents=True)
            (directory/'paired.root').write_bytes(b'fixture')
        payload=dict(flow={},schemes=[])
        for key in ['migration','topology','angles']:payload[key]=a[key].ravel().tolist()
        for si in range(3):payload['schemes'].append({k[2:]:v.ravel().tolist() for k,v in a.items() if k.startswith(str(si)+'_')})
        payload_path=root/'payload.json';payload_path.write_text(json.dumps(payload))
        original_compile=study.compile_worker
        def fake_compile(destination,production_source,fail=False):
            script=destination/'fake_worker.py'
            body='#!/usr/bin/env python3\nimport json,sys,shutil,time\n'
            if fail:
                body+='job=json.load(open(sys.argv[2]))\nif job["period"]=="fa18_inb":sys.exit(1)\ntime.sleep(15)\n'
            body+=f'shutil.copyfile({str(payload_path)!r},sys.argv[3])\n'
            script.write_text(body);script.chmod(0o755);return script
        try:
            study.compile_worker=fake_compile
            streamed=root/'streamed'
            assert study.main(['--production-run',str(run),'--input-root',str(dummy_input),'--workers','2','--dpi','40','--output',str(streamed)])==0
            with np.load(streamed/'aggregates_fa18_inb.npz') as z:
                assert all(np.array_equal(z[k],v) for k,v in a.items())
            study.compile_worker=lambda destination,production_source:fake_compile(destination,production_source,True)
            start=time.monotonic()
            try:study.main(['--production-run',str(run),'--input-root',str(dummy_input),'--workers','2','--output',str(root/'failure')])
            except RuntimeError:pass
            else:raise AssertionError('Failed worker accepted')
            assert time.monotonic()-start<8,'Failed scan left other workers running'
        finally:study.compile_worker=original_compile
        try:study.main(common+['--output',str(out)])
        except SystemExit:pass
        else:raise AssertionError('Existing study overwritten')
        (run/'output/jsons/combined_cuts.json').write_text(json.dumps(cuts)+'\n')
        try:study.main(common+['--output',str(root/'stale')])
        except ValueError:pass
        else:raise AssertionError('Stale cuts accepted for aggregate reuse')
    print('PASS: actual C++ cut/accumulator helpers, deterministic split, feed-in, both-period reports, matched coverage, driver/cancellation, reuse and overwrite/provenance guards (no ROOT runtime)')


def root_smoke(repo):
    with tempfile.TemporaryDirectory() as temporary:
        root=Path(temporary);run=root/'production';(run/'output/jsons').mkdir(parents=True)
        for name in ['global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h']:shutil.copy2(repo/name,run/name)
        # Serialize actual production defaults; all reconstructed branches are doubles,
        # as in the paired-MC production branch binder.
        fields=[n for n in re.findall(r'LOAD\((\w+)\)',(SOURCE/'low_t_closure_worker.cpp').read_text()) if n!='name']
        branches=['gen_kinematics_valid','rec_kinematics_valid','gen_x','gen_Q2','gen_t1','gen_phi2','gen_p1_theta',
                  'rec_x','rec_Q2','rec_t1','rec_phi2','rec_detector1','rec_detector2','rec_e_p','rec_e_theta','rec_e_phi',
                  'rec_p1_p','rec_p1_theta','rec_p1_phi','rec_p2_p','rec_p2_theta','rec_p2_phi','rec_open_angle_ep2','rec_Mx2','rec_runnum']
        code='''#include "global_cuts.h"
#include <TFile.h>
#include <TTree.h>
#include <nlohmann/json.hpp>
#include <map>
#include <fstream>
#include <random>
#include <cmath>
int main(int argc,char**argv){
 GlobalCutConfig c;c.enable_auxiliary_fiducial_cuts=true;
 c.auxiliary_cd_proton_theta_min_deg=0;
 nlohmann::json j;
'''+''.join(f'j["{n}"]=c.{n};\n' for n in fields)+'''
 std::ofstream config(argv[2]);config<<j;config.close();
 TFile file(argv[1],"RECREATE");TTree tree("MCEvents","paired fixture");
 std::map<std::string,double> v;int detector1=2,detector2=1,runnum=5000;
'''+''.join((f' tree.Branch("{n[4:]}",&{n[4:]});\n' if n in ('rec_runnum','rec_detector1','rec_detector2') else f' tree.Branch("{n}",&v["{n}"]);\n') for n in branches)+'''
 std::mt19937_64 rng(77);std::uniform_real_distribution<double> u(0,1);std::normal_distribution<double> smear(0,.02);
 const double rad=std::acos(-1.)/180.;
 for(int i=0;i<40000;++i){
  v["gen_kinematics_valid"]=1;v["rec_kinematics_valid"]=i%7!=0;
  v["gen_x"]=v["rec_x"]=.23;v["gen_Q2"]=v["rec_Q2"]=2.2;
  v["gen_t1"]=-.08-.22*u(rng);v["rec_t1"]=v["gen_t1"]+smear(rng);
  v["gen_phi2"]=(-10+50*u(rng))*rad;v["rec_phi2"]=v["gen_phi2"]+smear(rng);
  v["gen_p1_theta"]=(45+18*u(rng))*rad;v["rec_p1_theta"]=v["gen_p1_theta"];
  v["rec_e_p"]=i%5==0?1.5:5.;v["rec_e_theta"]=12*rad;v["rec_e_phi"]=0;
  v["rec_detector1"]=2;v["rec_detector2"]=1;v["rec_p1_p"]=.6;
  v["rec_p1_phi"]=120*rad;v["rec_p2_phi"]=240*rad;
  v["rec_p2_p"]=5;v["rec_p2_theta"]=12*rad;
  v["rec_open_angle_ep2"]=10;v["rec_Mx2"]=i%3==0?.01:.001;v["rec_runnum"]=5000;
  tree.Fill();
 }
 tree.Write();file.Close();return 0;
}
'''
        src=root/'fixture.cpp';src.write_text(code);exe=root/'fixture'
        flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
        libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
        subprocess.run(['c++',*flags,'-std=c++17','-I'+str(repo),str(src),'-o',str(exe),*libs],check=True)
        fixture=root/'paired.root';subprocess.run([str(exe),str(fixture),str(run/'output/jsons/global_cuts_config.json')],check=True)
        exclusive={};model={'mc':{'ep->epg':{}}}
        for p,(tag,code) in study.PERIODS.items():
            for topo in ['FD_FD','CD_FD','CD_FT']:
                exclusive[f'DVCS_{code}_{topo}']=dict(data={'Mx2':dict(cut_low=-1e-10,cut_high=1e-10)},
                    mc={'Mx2':dict(mode='quantile_window',cut_low=-.02 if p=='fa18_inb' else -.005,cut_high=.02 if p=='fa18_inb' else .005)})
            model['mc']['ep->epg'][tag]={region:dict(reference_factor=.85) for region in ['FT',*['S'+str(i) for i in range(1,7)]]}
            path=root/'input/mc'/p/'dvcs_mc/50nA/parts';path.mkdir(parents=True);shutil.copy2(fixture,path/'paired.root')
        (run/'output/jsons/combined_cuts.json').write_text(json.dumps(exclusive))
        p=run/'output/dvcs_current_dependence/calibration';p.mkdir(parents=True);(p/'current_response_model.json').write_text(json.dumps(model))
        p=run/'output/csvs';p.mkdir();rows=[{'bin index':r['id'],'valid bin':1,**dict(zip(study.FIELDS,r['bounds']))} for r in fixture_bins()]
        with (p/'dvcs_pass2_analysis.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        before=study.digest(p/'dvcs_pass2_analysis.csv')
        out=root/'study';common=['--production-run',str(run),'--workers','2','--dpi','40']
        assert study.main(common+['--input-root',str(root/'input'),'--output',str(out)])==0
        flow=json.loads((out/'cutflow.json').read_text())
        a,b=flow['fa18_inb'],flow['fa18_out']
        assert a['global_pass']==b['global_pass']
        assert .6*a['selected']<b['selected']<.75*a['selected']
        assert a['global_pass']<a['valid_reco'] # production electron/physical cuts active
        with np.load(out/'aggregates_fa18_inb.npz') as z:
            g=z['0_train_gen'][0];r=z['0_train_reco'][0]
            assert np.all(g>0) and np.all(r>0)
            assert z['migration'][:11,:4].sum()>0 # subthreshold truth feeding selected reco
        assert (out/'closure_by_cell.csv').is_file() and (out/'closure_nominal_vs_merged.png').is_file()
        replay=root/'replay';assert study.main(common+['--reuse-aggregates',str(out),'--output',str(replay)])==0
        assert (replay/'closure_by_cell.csv').read_bytes()==(out/'closure_by_cell.csv').read_bytes()
        assert study.digest(p/'dvcs_pass2_analysis.csv')==before
        try:study.main(common+['--reuse-aggregates',str(out),'--output',str(replay)])
        except SystemExit:pass
        else:raise AssertionError('Existing output overwritten')
        cuts=run/'output/jsons/combined_cuts.json';cuts.write_text(cuts.read_text()+'\n')
        try:study.main(common+['--reuse-aggregates',str(out),'--output',str(root/'stale')])
        except ValueError:pass
        else:raise AssertionError('Stale cut configuration accepted')
        print('PASS: native ROOT end-to-end for both polarities, production phase space, MC rather than DATA windows, feed-in, all schemes, aggregate reuse, provenance/overwrite protection')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root-smoke',action='store_true');parser.add_argument('--source-dir',type=Path,default=SOURCE);args=parser.parse_args()
    unit()
    native_and_replay(args.source_dir)
    if args.root_smoke:root_smoke(args.source_dir)
