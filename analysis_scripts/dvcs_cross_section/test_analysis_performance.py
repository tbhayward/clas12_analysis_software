#!/usr/bin/env python3
"""Equivalence checks for performance changes; no farm event files needed."""
import ast
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

SOURCE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True


def function_source(text, name):
    tree = ast.parse(text)
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    return '\n'.join(text.splitlines()[node.lineno-1:node.end_lineno])


def model_checks(temp):
    # Load the actual batch client without plotting dependencies.
    text = (SOURCE/'external_scripts/plot_pass2_vs_pass1_model_comparison.py').read_text()
    ns = dict(os=os, sys=sys, json=json, math=__import__('math'), subprocess=subprocess,
              to_float=lambda x, default=0.: float(x) if x else default)
    # Scalar is the existing implementation, including its .12g input formatting.
    ns.update(Sequence=list, Optional=__import__('typing').Optional, Dict=__import__('typing').Dict,
              ModelConfig=object, shlex=shlex)
    for name in ('run_command_capture','km15_xs','km15_xs_many'):
        exec(function_source(text,name),ns)
    fixture = temp/'models';fixture.mkdir()
    shutil.copy2(SOURCE/'km15_cli.py',fixture/'km15_cli.py')
    package=fixture/'gepard';package.mkdir()
    (package/'__init__.py').write_text('class DataPoint:\n def __init__(self, **kw): self.__dict__.update(kw)\n def prepare(self): pass\n')
    (package/'fits.py').write_text('''import math
class Theory:
 def predict(self, p):
  if p.Q2 > 50: raise ValueError('fixture model failure')
  if p.Q2 > 40: return float('nan')
  return p.xB*p.Q2 + math.sin(p.phi) + p.in1energy + p.in1polarization*.1234567891234
th_KM15=Theory()
''')
    from types import SimpleNamespace
    cfg=SimpleNamespace(py_km15=sys.executable, km15_cli=str(fixture/'km15_cli.py'),e_beam=10.604)
    points=[(.23456789123456,2.345678912345,.17,angle,h)
            for angle in (0.,17.123456789123,180.,360.) for h in (0,+1,-1)]
    points += [(-.1,2.,.2,0.,0),(.2,45.,.2,0.,0),(.2,55.,.2,0.,0)]
    scalar=lambda: [ns['km15_xs'](*p[:4],cfg,helicity=p[4]) for p in points]
    start=time.perf_counter();reference=scalar();scalar_time=time.perf_counter()-start
    start=time.perf_counter();actual=ns['km15_xs_many'](points,cfg);batch_time=time.perf_counter()-start
    assert actual==reference
    os.environ['DVCS_VERIFY_KM15_BATCH']='1'
    try: assert ns['km15_xs_many'](points,cfg)==reference
    finally: os.environ.pop('DVCS_VERIFY_KM15_BATCH')
    # Check the exact scalar text for XS and XSD, not just float agreement.
    requests=[['.234567891235','2.34567891235','.17',str(a),'10.604',str(h),mode]
              for a in (0,180,360) for h in (0,1,-1) for mode in ('XS','XSD')]
    results=json.loads(subprocess.check_output([sys.executable,cfg.km15_cli,'--batch'],
                       input=json.dumps(requests),text=True))['values']
    assert results==[subprocess.check_output([sys.executable,cfg.km15_cli,*r],text=True).strip() for r in requests]
    # Old custom CLI fallback.
    old=fixture/'old_cli.py';old.write_text('print("2.25")\n');cfg.km15_cli=str(old)
    assert ns['km15_xs_many'](points,cfg)==[2.25]*len(points)
    # Exercise the actual final-results curve loop against its scalar path.
    cfg.km15_cli=str(fixture/'km15_cli.py');cfg.phi_dense=4
    ns['km15_xs_many']=ns['km15_xs_many']
    ns['bh_xs']=lambda xb,q,t,phi,c: 1.+phi/360.
    base=SimpleNamespace(**ns)
    scope=dict(base=base,finite=__import__('math').isfinite)
    exec(function_source((SOURCE/'external_scripts/plot_pass2_final_results_by_t.py').read_text(),
                         'compute_cell_models'),scope)
    key=SimpleNamespace(xb_min=.2,xb_max=.3,q2_min=2.,q2_max=3.,t_min=.11,t_max=.25)
    for kinds in [('xs',),('bsa',),('xs','bsa')]:
        a=scope['compute_cell_models']((key,cfg,kinds))[1]
        os.environ['DVCS_KM15_SCALAR']='1'
        try: b=scope['compute_cell_models']((key,cfg,kinds))[1]
        finally: os.environ.pop('DVCS_KM15_SCALAR')
        assert a==b
    if '--real-km15' in sys.argv:
        cfg.km15_cli = str(SOURCE/'km15_cli.py')
        actual_points = [(.25,2.5,.2,angle,h)
                         for angle in (0.,137.123456789123,360.) for h in (0,1,-1)]
        start = time.perf_counter()
        real_scalar = [ns['km15_xs'](*p[:4],cfg,helicity=p[4]) for p in actual_points]
        scalar_seconds = time.perf_counter()-start
        start = time.perf_counter()
        real_batch = ns['km15_xs_many'](actual_points,cfg)
        batch_seconds = time.perf_counter()-start
        assert any(v > 0 for v in real_scalar), 'Real KM15 unavailable in this Python environment'
        assert real_scalar == real_batch
        print(f'PASS: real KM15 exact equality; {scalar_seconds:.3f}s -> {batch_seconds:.3f}s')
    print(f'PASS: KM15 scalar/batch text, errors, custom CLI fallback and complete XS/BSA curves; fixture timing {scalar_time:.3f}s -> {batch_time:.3f}s')


def cut_checks(temp):
    text=(SOURCE/'current_dependence.cpp').read_text()
    types=text[text.index('struct SigmaStats {'):text.index('static TopoCutMap load_sigma_cuts')]
    helpers=text[text.index('static bool within_cut_window'):text.index('// Branch binding',text.index('static bool within_cut_window'))]
    helpers=helpers[:helpers.rfind('// -----------------------------------------------------------------------------')]
    fields=text[text.index('struct Branches {'):text.index('    void bind(',text.index('struct Branches {'))]
    method=text[text.index('    double delta_phi_value('):text.index('static bool passes_cone_cut')]
    dispatch=text[text.index('static bool passes_sigma_dispatch_legacy'):text.index('// Data and MC aggregations')]
    dispatch=dispatch[:dispatch.rfind('// -----------------------------------------------------------------------------')]
    delta=text[text.index('static inline double delta_phi_rad_from_two_phi'):text.index('static const std::vector<std::string> PERIOD_ORDER')]
    cpp='''#include <array>
#include <cmath>
#include <string>
#include <unordered_map>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <random>
#include <cassert>
#include <iostream>
#include <chrono>
#include <cstdlib>
constexpr double PI=3.14159265358979323846;
[[noreturn]] static void fatal(const std::string& s) {throw std::runtime_error(s);}
enum class Channel {DVCS,EPPI0};
struct ChannelConfig {Channel channel; std::string cut_prefix;};
struct PeriodTags {std::string period_code;};
'''+types+helpers+delta+fields+method+dispatch+r'''
int main() {
 std::mt19937 rng(271828); std::uniform_real_distribution<double> dist(-8.,8.);
 const double nan=std::numeric_limits<double>::quiet_NaN();
 // Compare all fallback/window modes including exact edges and nonfinite values.
 for(int i=0;i<10000;++i){
  SigmaStats s; s.mean=dist(rng);s.std=dist(rng);s.cut_low=dist(rng);s.cut_high=dist(rng);
  if(i%3==0)s.mode="upper_quantile";
  if(i%5==0)s.cut_low=nan;
  if(i%7==0)s.cut_high=nan;
  if(i%11==0)s.mean=nan;
  if(i%13==0)s.std=nan;
  const auto w=prepare_sigma_window(s);
  for(double v:{dist(rng),s.cut_low,s.cut_high,s.mean-3*s.std,s.mean+3*s.std,nan,
      std::numeric_limits<double>::infinity()})
   assert(w.accepts(v)==within_cut_window(v,s));
 }
 TopoCutMap cuts;
 const std::array<std::string,7> vars{{"Delta_phi","theta","theta_gamma_gamma","pTmiss","Emiss2","Mx2","Mx2_2"}};
 for(auto prefix:{"dvcs","eppi0"})for(auto period:{"Fa18Inb","Sp19Inb"})for(auto topo:{"FD_FD","CD_FD","CD_FT"}) {
  CutVarMap vm;
  for(size_t i=0;i<vars.size();++i){SigmaStats s;s.mean=0;s.std=1;s.cut_low=-3;s.cut_high=3;if(i==3){s.mode="upper_quantile";s.cut_high=2;}vm[vars[i]]=s;}
  vm.erase("theta"); cuts[std::string(prefix)+"_"+period+"_"+topo]=vm;
 }
 Branches b;
 b.has_detector1=b.has_detector2=b.has_Delta_phi=b.has_theta=true;
 b.has_theta_gamma_gamma=b.has_theta_pi0_pi0=b.has_pTmiss=b.has_Emiss2=b.has_Mx2=b.has_Mx2_2=true;
 ChannelConfig cfg{Channel::DVCS,"dvcs"};PeriodTags tags{"Fa18Inb"};
 for(int i=0;i<100000;++i){
  cfg.channel=i%2?Channel::DVCS:Channel::EPPI0;cfg.cut_prefix=i%2?"dvcs":"eppi0";
  tags.period_code=i%3?"Fa18Inb":"Sp19Inb";
  b.detector1=i%3==0?1:2;b.detector2=i%3==2?0:1;
  b.Delta_phi=dist(rng);b.theta=dist(rng);b.theta_gamma_gamma=dist(rng);b.theta_pi0_pi0=dist(rng);
  b.pTmiss=dist(rng);b.Emiss2=dist(rng);b.Mx2=dist(rng);b.Mx2_2=dist(rng);
  assert(passes_sigma_dispatch(cfg,tags,cuts,b)==passes_sigma_dispatch_legacy(cfg,tags,cuts,b));
 }
 // Missing branch failures retain exactly the same diagnostics/short-circuiting.
 cfg={Channel::DVCS,"dvcs"};tags={"Fa18Inb"};b.detector1=b.detector2=1;
 b.Delta_phi=b.theta_gamma_gamma=b.pTmiss=b.Emiss2=b.Mx2=b.Mx2_2=0;
 b.has_Mx2=false;
 std::string a,c;
 try{passes_sigma_dispatch(cfg,tags,cuts,b);}catch(const std::exception& e){a=e.what();}
 try{passes_sigma_dispatch_legacy(cfg,tags,cuts,b);}catch(const std::exception& e){c=e.what();}
 assert(!a.empty()&&a==c);b.has_Mx2=true;
 // Independent cuts with equal topology keys must not reuse windows.
 TopoCutMap alternate=cuts;alternate["dvcs_Fa18Inb_FD_FD"]["Mx2"].cut_low=1.;
 assert(passes_sigma_dispatch(cfg,tags,alternate,b)==passes_sigma_dispatch_legacy(cfg,tags,alternate,b));
 volatile int total=0;
 for(bool fast:{false,true}){
  const auto start=std::chrono::steady_clock::now();
  for(int i=0;i<1000000;++i) total+=fast?passes_sigma_dispatch(cfg,tags,cuts,b):passes_sigma_dispatch_legacy(cfg,tags,cuts,b);
  std::cout<<(fast?"prepared":"legacy")<<" cut fixture seconds: "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 }
 std::cout<<"PASS: prepared/reference cuts, boundaries, channels, topologies, missing branches and cache identities\n";
}
'''
    path=temp/'cuts.cpp';path.write_text(cpp);exe=temp/'cuts'
    subprocess.run(['c++','-std=c++17','-O2',str(path),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)


def root_checks(temp):
    if not shutil.which('root-config'):
        print('SKIP: ROOT cache equivalence (root-config unavailable)');return
    text=(SOURCE/'bin_means.cpp').read_text()
    setup=text[text.index('        const bool uncached ='):text.index('        auto bindD =',text.index('        const bool uncached ='))]
    cpp='''#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <iostream>
int main(int argc,char** argv){
 std::string dir=argv[1];
 for(int part=0;part<3;++part){
  TFile f((dir+"/part"+std::to_string(part)+".root").c_str(),"RECREATE");
  TTree tree("PhysicsEvents","");double x=0,y=0,unused[128]{};
  tree.Branch("x",&x);tree.Branch("y",&y);tree.Branch("unused",unused,"unused[128]/D");
  for(int i=0;i<10000;++i){x=part*10000+i;y=x*.123456789;unused[i%128]=x;tree.Fill();}tree.Write();
 }
 std::vector<std::vector<double>> outputs;
 for(bool cache:{false,true}){
  if(cache)unsetenv("DVCS_BIN_MEANS_NO_CACHE");else setenv("DVCS_BIN_MEANS_NO_CACHE","1",1);
  TChain chain("PhysicsEvents");for(int part=0;part<3;++part)chain.Add((dir+"/part"+std::to_string(part)+".root").c_str());
  TTree* t=&chain;t->SetBranchStatus("*",0);t->SetBranchStatus("x",1);t->SetBranchStatus("y",1);
'''+setup+'''
  double x=0,y=0;t->SetBranchAddress("x",&x);t->SetBranchAddress("y",&y);
  std::vector<double> result;double sx=0,sy=0;
  for(Long64_t i=0;i<t->GetEntries();++i){assert(t->GetEntry(i)>0);result.push_back(x);result.push_back(y);sx+=x;sy+=y;}
  result.push_back(sx);result.push_back(sy);outputs.push_back(result);
 }
 assert(outputs[0]==outputs[1]);
 std::cout<<"PASS: every ROOT value and accumulated sum identical across three cached/uncached chain files\\n";
}
'''
    path=temp/'cache.cpp';path.write_text(cpp);exe=temp/'cache'
    flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
    subprocess.run(['c++',str(path),'-std=c++17',*flags,'-o',str(exe)],check=True)
    subprocess.run([str(exe),str(temp)],check=True)


def main():
    with tempfile.TemporaryDirectory() as d:
        temp=Path(d)
        model_checks(temp)
        cut_checks(temp)
        if '--skip-root' in sys.argv:
            print('SKIP: ROOT runtime check explicitly disabled')
        else:
            root_checks(temp)

if __name__=='__main__':main()
