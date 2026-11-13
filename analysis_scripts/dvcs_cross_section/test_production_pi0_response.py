#!/usr/bin/env python3
import argparse,csv,json,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
import study_production_pi0_response as study


def unit(tmp):
    cols=study.audit.EDGES+['valid bin','reconstructed yield, ep->eppi0, mc, Fa18 Inb','untouched']
    row={k:'0' for k in cols};row.update({'valid bin':'1','reconstructed yield, ep->eppi0, mc, Fa18 Inb':'(20,4,0)','untouched':'DATA'})
    baseline=tmp/'base.csv'
    with baseline.open('w') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerow(row)
    target=tmp/'counts.csv';study.prepare_csv(baseline,target)
    with target.open() as f:r=next(csv.DictReader(f))
    assert r['untouched']=='DATA' and study.audit.triple(r[cols[-2]])[0]==0
    try:study.gate(baseline,target,tmp)
    except ValueError:pass
    else:raise AssertionError('Gate accepted changed counts')
    shutil.copy2(baseline,target);study.gate(baseline,target,tmp)
    text=(study.SOURCE/'total_counts.cpp').read_text();modified=study.instrument(text)
    assert modified.count('apply_pi0_test_response(b,work_cfg,tags,i);')==1
    assert (study.SOURCE/'total_counts.cpp').read_text()==text
    counts=tmp/'nominal.csv';counts.write_text('counts')
    assert not study.completed_count(tmp,'nominal')
    (tmp/'nominal.log').write_text('[total_counts] Updated data and MC count columns in: '+str(counts))
    assert study.completed_count(tmp,'nominal')
    study.response.dump(tmp/'nominal_count_complete.json',dict(csv_sha256=study.response.digest(counts)))
    counts.write_text('changed')
    try:study.completed_count(tmp,'nominal')
    except ValueError:pass
    else:raise AssertionError('Accepted changed completed CSV')
    print('PASS: isolated-source instrumentation, zero-before-recount, untouched DATA, nominal agreement and mismatch gate')


def native(tmp):
    backend=tmp/'backend';exe=study.build(backend)
    fixture=backend/'fixture.cpp';fixture.write_text(r'''
#define main comparison_driver_main
#include "counter.cpp"
#undef main
#include <cassert>
int main(){
 BranchBinder b;b.x=.2;b.Q2=2.2;b.t1=-.2;b.detector1=2;b.detector2=1;b.e_p=6;b.e_theta=.2;b.e_phi=.1;b.p1_p=.6;b.p1_theta=.8;b.p1_phi=2;b.p2_p=3;b.p2_theta=.25;b.p2_phi=4;b.open_angle_ep2=10;b.Delta_phi=3;b.Mx2=.001;b.runnum=5032;b.phi2=6.28;
 b.has_e_p=b.has_e_theta=b.has_e_phi=b.has_p1_p=b.has_p1_theta=b.has_p1_phi=b.has_p2_p=b.has_p2_theta=b.has_p2_phi=b.has_Delta_phi=true;
 WorkConfig work;work.sample_kind=SampleKind::MC_REC;work.channel_cfg=eppi0_bkg_as_dvcs_config();PeriodTags tags;tags.period_label="fa18_inb";tags.tree_key="fixture";
 std::array<double,20> p{};pi0_test_parameters["fa18_inb"]=p;
 auto before=b;apply_pi0_test_response(b,work,tags,0);assert(b.p2_p==before.p2_p&&b.phi2==before.phi2&&b.t1==before.t1);
 p[15]=.01;pi0_test_parameters["fa18_inb"]=p;apply_pi0_test_response(b,work,tags,0);assert(std::abs(b.p2_p-3.03)<1e-12);assert(std::abs(b.Emiss2-before.Emiss2+.03)<1e-10);assert(b.phi2>=0&&b.phi2<2*particle_response::PI);
 b=before;work.channel_cfg=eppi0_config();bool rejected=false;try{apply_pi0_test_response(b,work,tags,0);}catch(...){rejected=true;}assert(rejected);
 b=before;b.response_extra_ok.fill(true);b.response_extra.fill(0);apply_pi0_test_response(b,work,tags,1);
 assert(b.response_unresolved&&b.p2_p==before.p2_p&&b.phi2==before.phi2&&b.t1==before.t1);
 audit_pi0_recovery_row(b,work,tags,7);flush_pi0_recovery_audit();
 assert(pi0_recovery_counts[tags.tree_key][1]==1&&pi0_recovery_bins[tags.tree_key+"|7"][1]==1);

}
''')
    flags=subprocess.check_output(['root-config','--cflags'],text=True).split();libs=subprocess.check_output(['root-config','--libs'],text=True).split();binary=backend/'fixture'
    subprocess.run(['c++',*flags,'-O1','-std=c++17','-pthread','-I'+str(backend),str(fixture),str(study.SOURCE/'global_cuts.cpp'),str(study.SOURCE/'load_trees.cpp'),'-o',str(binary),*libs],check=True);subprocess.run([str(binary)],check=True)
    print('PASS: real production-counter backend compilation, exact zero-response identity, particle bias/derived missing energy, phi wrapping and missing-branch rejection and unchanged unresolved-event retention with selected-bin audit')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');args=p.parse_args()
    with tempfile.TemporaryDirectory() as d:
        unit(Path(d))
        if args.native:native(Path(d))
