#!/usr/bin/env python3
import argparse,copy,json,shutil,subprocess,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_particle_response as base


def scoring():
    edges=np.linspace(-1,1,80);x=(edges[:-1]+edges[1:])/2;good=np.r_[0,100*np.exp(-x*x/.08),0];bad=np.roll(good,8);n=len(good)
    spec=dict(periods=['fa18_inb'],variables=[dict(bins=n-2) for _ in base.OBS],joint_control_fit=True)
    data={'fa18_inb__pi0_data':np.tile(good,(3,1,7,1))}
    engine=base.Engine(None,spec,{},data,{},['CD_FT'],10)
    h=np.tile(good,(2,2,3,1,7,1));h[1,0,2,0,0]=bad
    results=dict(train=0.,validation_mc=0.,heldout_observables=0.,training_projections=1)
    with patch.object(base,'TRAIN_OBS',[0]),patch.object(base,'VALID_OBS',[6]):
        nominal,_=engine.joint_control_score({'fa18_inb__pi0_mc':h},0,results.copy(),[])
        candidate,records=engine.joint_control_score({'fa18_inb__pi0_mc':h},1,results.copy(),[])
        assert nominal['train']==0 and candidate['train']>0
        h[1,1,2,0,0]=bad*2
        changed,_=engine.joint_control_score({'fa18_inb__pi0_mc':h},1,results.copy(),[])
        assert changed['train']==candidate['train'] and changed['validation_mc']>candidate['validation_mc']
        assert all(r['channel']=='epgammagamma' for r in records)
    print('PASS: simultaneous control shape penalty, equal channel weighting, control channel records, validation MC cannot change training objective')


def native(tmp):
    source=base.SOURCE;ref=tmp/'reference';(ref/'backend').mkdir(parents=True);buildsource=tmp/'sources';buildsource.mkdir();out=tmp/'worker';out.mkdir()
    for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(source/n,ref/'backend'/n)
    for n in ('exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h','particle_response_core.h','particle_response_worker.cpp','particle_control_core.h','particle_control_worker.h','particle_response_cut_arrays.h'):shutil.copy2(source/n,buildsource/n)
    shutil.copy2(Path(__file__).parent/'particle_control_worker.h',buildsource/'particle_control_worker.h')
    with patch.object(base,'SOURCE',buildsource):worker=base.compile_worker(out,ref)
    # Generate a complete config using the struct defaults, not a user/farm path.
    import re
    fields=[n for n in re.findall(r'LOAD\((\w+)\)',(source/'exclusivity_resolution_worker.cpp').read_text()) if n!='name']
    cpp=tmp/'cfg.cpp';cpp.write_text('#include "global_cuts.h"\n#include <nlohmann/json.hpp>\n#include <fstream>\nint main(int argc,char**argv){GlobalCutConfig c;nlohmann::json j;'+''.join('j["'+n+'"]=c.'+n+';' for n in fields)+'std::ofstream(argv[1])<<j;}')
    flags=subprocess.check_output(['root-config','--cflags'],text=True).split();exe=tmp/'cfg';subprocess.run(['c++',*flags,'-std=c++17','-I'+str(source),str(cpp),'-o',str(exe)],check=True);subprocess.run([str(exe),str(tmp/'cfg.json')],check=True)
    cfg=json.loads((tmp/'cfg.json').read_text());cfg.update(enable_auxiliary_fiducial_cuts=False,enable_analysis_phase_space_cuts=False,enable_reconstructed_electron_p_cut=False,enable_pTmiss_cut=False,enable_topology_filter=False,enable_electron_fd_sector_filter=False,enable_proton_fd_sector_filter=False,enable_proton_cd_sector_filter=False,enable_photon_fd_sector_filter=False,excluded_runs=[])
    spec=dict(period_codes={'fa18_inb':'Fa18_Inb'},global_cuts=cfg,joint_control_fit=True,variables=[dict(name=n,bins=40,low=-4,high=4) for n in base.OBS],exclusivity={'eppi0_Fa18_Inb_'+t:dict(data={},mc={}) for t in base.TOPS})
    dtype=np.dtype([('v','<f8',(24,)),('key','<u8'),('photons','<f8',(10,))]);assert dtype.itemsize==280
    a=np.zeros(20,dtype=dtype);v=a['v'];v[:,0]=.2;v[:,1]=2.2;v[:,2]=-.2;v[:,3]=2;v[:,4]=1;v[:,5]=6;v[:,6]=.2;v[:,7]=0;v[:,8]=.6;v[:,9]=.8;v[:,10]=2;v[:,11]=3;v[:,12]=.25;v[:,13]=4;v[:,14]=10;v[:,15]=.05;v[:,16:23]=.01;v[:,23]=5032;a['key']=np.arange(20)*4
    a['photons'][:]=[2,.2,.5,1,.3,1,1,0,.135,1]
    cache=tmp/'cache.bin';a.tofile(cache)
    for joint in (False,True):
        spec['joint_control_fit']=joint;base.dump(tmp/'spec.json',spec);base.dump(tmp/'job.json',dict(period='fa18_inb',sample='pi0_mc',cache=str(cache),parameters=[[0.]*20]))
        target=tmp/('joint.bin' if joint else 'legacy.bin');subprocess.run([str(worker),'evaluate',str(tmp/'spec.json'),str(tmp/'job.json'),str(target)],check=True)
        meta=json.loads(Path(str(target)+'.json').read_text())
        if joint:
            hist=np.fromfile(target,dtype='<i8').reshape(meta['hist_shape']);assert hist.sum()==20*7 and hist.shape==(1,2,3,1,7,42)
        else:assert target.stat().st_size==6*8
    # Each photon uses its own detector category, not the first photon's tag.
    fixture=tmp/'tags.cpp';fixture.write_text(r'''
#include "particle_control_core.h"
#include <cassert>
using namespace particle_response;
int main(){ControlRecord c{};c.base.v[3]=2;c.base.v[4]=1;c.base.v[5]=6;c.base.v[6]=.2;c.base.v[8]=.6;c.base.v[9]=.8;c.base.v[11]=3;c.base.v[12]=.25;c.base.v[13]=1;c.photons={2,.2,.5,1,.3,1,1,0,.135,1};std::array<double,20> p{};std::array<double,12> z{};double mass;p[19]=.1;auto old=derived(c.base.v,10.604,false,MPI0);auto changed=vary_control(c,p,10.604,old,z,mass);auto expected=particle(2,.2,.5,0)+particle(1.1,.3,1,0);assert(std::abs(changed[11]-expected.p.norm())<1e-12);}
''');binary=tmp/'tags';subprocess.run(['c++','-std=c++17','-I'+str(source),str(fixture),'-o',str(binary)],check=True);subprocess.run([str(binary)],check=True)
    print('PASS: native control histogram accumulation/shape, legacy selected-count output preserved, FD–FT photons receive separate detector responses')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--native',action='store_true');args=parser.parse_args();scoring()
    if args.native:
        with tempfile.TemporaryDirectory() as d:native(Path(d))
