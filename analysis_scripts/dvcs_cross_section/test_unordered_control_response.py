#!/usr/bin/env python3
import argparse,json,subprocess,tempfile
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_unordered_control_response as study


def selections(tmp):
    a=np.zeros((1200,29));a[:,0]=np.arange(len(a))%6;a[:,8]=1;a[:,11]=(np.arange(len(a))//6)%2;a[:,20]=3;a[:,21]=.1;a[:,23]=1;a[:,24]=1;a[:,25]=.05;a[:,27]=0;a[:,28]=.135
    rng=np.random.default_rng(5)
    for col in range(1,8):a[:,col]=rng.normal(0,.1,len(a))
    spec=dict(periods=['fa18_inb'],variables=[dict(name=n,label=n,bins=80,low=-1,high=1) for n in base.OBS],exclusivity={'eppi0_Fa18_Inb_'+top:dict(data={'pTmiss':dict(cut_low=-.2,cut_high=.2)}) for top in base.TOPS})
    spec['exclusivity']['eppi0_Fa18_Inb_CD_FT']['data']['pTmiss']['cut_high']=.1
    masks=[study.selection(a,g,study.blocks(spec,'fa18_inb',g)) for g in range(6)];assert np.all(np.sum(masks,axis=0)<=1)
    for g in range(3,6):assert np.all(a[masks[g],4]<=.1)
    assert study.selection(a,3,study.blocks(spec,'fa18_inb',3),omit='pTmiss').sum()==200
    h=study.histogram(a,study.selection(a,3,[]),4,np.linspace(-1,1,81));assert h.sum()==200
    rows,flows,links,metrics=study.render('fa18_inb',{'DATA':a,'nominal':a,'frozen':a},spec,tmp,10,50)
    assert len(links)==12 and all((tmp/x).is_file() for x in links) and len(rows)==384 and len(flows)==36 and len(metrics)==12
    print('PASS: six disjoint groups, common window intersection, N-minus-one selection, both MC halves, pi0 mass and complete paired report plots')


def native(tmp,root_io=False):
    import test_frozen_response_exports as fixture
    fixture.native(tmp)
    out=tmp/'unordered';out.mkdir();work=tmp/'scratch';work.mkdir();worker=study.compile_worker(out,tmp/'reference',work)
    spec=json.loads((tmp/'spec.json').read_text());spec['global_cuts']['enable_auxiliary_fiducial_cuts']=False;base.dump(tmp/'spec.json',spec)
    # Every detector pair in both proton regions, including old unsupported FD-FT.
    dt=np.dtype([('v','<f8',(24,)),('key','<u8'),('photons','<f8',(10,))]);records=np.fromfile(tmp/'cache.bin',dtype=dt);base_record=records[1].copy();a=np.repeat(np.array([base_record]),12)
    for i in range(12):
        g=i%6;d1=1 if g<3 else 2;pair=[(1,1),(1,0),(0,0)][g%3]
        if i>=6:pair=pair[::-1]
        a['v'][i,3]=d1;a['v'][i,4]=pair[0];a['photons'][i,6:8]=pair;a['key'][i]=i*4
    cache=tmp/'all_pairs.bin';a.tofile(cache);job=dict(period='fa18_inb',sample='pi0_mc',cache=str(cache),parameters=[0.]*20,export_individual_photons=True,export_global_only=True);base.dump(tmp/'unordered_job.json',job);target=tmp/'pairs.bin'
    subprocess.run([str(worker),'cut-arrays',str(tmp/'spec.json'),str(tmp/'unordered_job.json'),str(target)],check=True);shape=json.loads(Path(str(target)+'.json').read_text())['shape'];rows=np.fromfile(target,dtype='<f8').reshape(shape)
    assert shape==[12,29] and np.array_equal(rows[:,0],np.arange(12)%6) and np.allclose(rows[:,28],.135)
    # A photon-topology filter must treat either mixed ordering identically.
    spec['global_cuts'].update(enable_topology_filter=True,required_detector1=2,required_detector2=1);base.dump(tmp/'spec.json',spec)
    subprocess.run([str(worker),'cut-arrays',str(tmp/'spec.json'),str(tmp/'unordered_job.json'),str(tmp/'filtered.bin')],check=True);filtered=np.fromfile(tmp/'filtered.bin',dtype='<f8').reshape(-1,29);assert len(filtered)==2 and np.all(filtered[:,0]==3)
    print('PASS: cached-MC native export, all six groups, order-independent globals, mass and topology filter')
    if not root_io:return
    # Actual DATA ROOT branch binding/recovery, including both pair orders.
    cpp=tmp/'data.cpp';cpp.write_text(r'''#include "particle_control_core.h"
#include <TFile.h>
#include <TTree.h>
using namespace particle_response;
int main(int argc,char**argv){TFile f(argv[1],"RECREATE");TTree t("PhysicsEvents","");double v[24]={.2,2.2,-.2,2,1,5,.2,.3,.7,.8,2,5,.12,1.3,10,.05,.01,.01,.01,.01,.01,.01,.01,5032};
const char* names[]={"x","Q2","t1","detector1","detector2","e_p","e_theta","e_phi","p1_p","p1_theta","p1_phi","p2_p","p2_theta","p2_phi","open_angle_ep2","pT","Delta_phi","theta","theta_pi0_pi0","pTmiss","Emiss2","Mx2","Mx2_2","runnum"};
for(int i=0;i<24;++i)t.Branch(names[i],&v[i],(std::string(names[i])+"/D").c_str());double mass=MPI0,d1,d2,az1,az2,op1,op2;
for(auto p:std::vector<std::pair<const char*,double*>>{{"Mh_gammagamma",&mass},{"detector_gamma1",&d1},{"detector_gamma2",&d2},{"gamma_phi1",&az1},{"gamma_phi2",&az2},{"open_angle_egamma1",&op1},{"open_angle_egamma2",&op2}})t.Branch(p.first,p.second,(std::string(p.first)+"/D").c_str());
F pi=particle(5,.12,1.3,MPI0),e=particle(5,.2,.3,ME),k={{0,0,std::sqrt(10.604*10.604-ME*ME)},10.604};F a=particle(MPI0/2,std::acos(.5),1,0),b={a.p*(-1),a.e};a=a.boost(pi.p*(1/pi.e));b=b.boost(pi.p*(1/pi.e));
for(int i=0;i<12;++i){int g=i%6;v[3]=g<3?1:2;d1=g%3==2?0:1;d2=g%3==0?1:0;F one=a,two=b;if(i>=6){std::swap(one,two);std::swap(d1,d2);}v[4]=d1;az1=trento((k-e).p,e.p,one.p);az2=trento((k-e).p,e.p,two.p);op1=angle(e.p,one.p)*180/PI;op2=angle(e.p,two.p)*180/PI;t.Fill();}t.Write();}
''')
    flags=subprocess.check_output(['root-config','--cflags'],text=True).split();libs=subprocess.check_output(['root-config','--libs'],text=True).split();exe=tmp/'data_writer';subprocess.run(['c++',*flags,'-std=c++17','-I'+str(base.SOURCE),str(cpp),'-o',str(exe),*libs],check=True);subprocess.run([str(exe),str(tmp/'data.root')],check=True)
    spec['global_cuts']['enable_topology_filter']=False;base.dump(tmp/'spec.json',spec);job=dict(period='fa18_inb',sample='pi0_data',files=[str(tmp/'data.root')],export_individual_photons=True,export_global_only=True);base.dump(tmp/'unordered_job.json',job)
    subprocess.run([str(worker),'cut-arrays',str(tmp/'spec.json'),str(tmp/'unordered_job.json'),str(tmp/'data_rows.bin')],check=True);rows=np.fromfile(tmp/'data_rows.bin',dtype='<f8').reshape(-1,29);assert len(rows)==12 and np.array_equal(rows[:,0],np.arange(12)%6) and np.all(np.isfinite(rows[:,20])) and np.allclose(rows[:,28],.1349766)
    print('PASS: actual DATA and cached-MC ROOT paths, all six groups, permutation-invariant globals, mass column and topology-filter consistency')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root-smoke',action='store_true');p.add_argument('--root-io',action='store_true');args=p.parse_args()
    with tempfile.TemporaryDirectory() as folder:selections(Path(folder))
    if args.root_smoke or args.root_io:
        with tempfile.TemporaryDirectory() as folder:native(Path(folder),args.root_io)
