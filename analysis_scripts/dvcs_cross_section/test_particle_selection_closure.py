#!/usr/bin/env python3
"""Check global phase-space selection, physical mixture scales and cut inputs."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import study_particle_selection_closure as study
base=study.base;cuts=study.cuts;production=cuts.production

def check_loader():
    cfg=SimpleNamespace(t1_abs_max=1.,open_angle_min_deg=5.,enable_reconstructed_electron_p_cut=False,enable_analysis_phase_space_cuts=True,enable_analysis_y_min_cut=True,analysis_y_min=.19,analysis_W_min_GeV=2.,enable_pTmiss_cut=False,pTmiss_max=.2,excluded_runs=(999,),enable_topology_filter=False,enable_electron_fd_sector_filter=False,enable_proton_fd_sector_filter=False,enable_proton_cd_sector_filter=False,enable_photon_fd_sector_filter=False,enable_auxiliary_fiducial_cuts=False)
    arrays={k:np.zeros(4) for k in ['t1','open_angle_ep2','pTmiss','e_phi','p1_phi','p2_phi','e_p','x','Q2','runnum','detector1','detector2']}
    arrays.update(t1=np.full(4,-.2),open_angle_ep2=np.full(4,10.),e_p=np.full(4,6.),x=np.array([.2,.8,.2,.2]),Q2=np.full(4,2.),runnum=np.array([5032,5032,999,5032]),detector1=np.ones(4),detector2=np.ones(4))
    for i,v in enumerate(production.VARIABLES):arrays[v.branch]=np.arange(4,dtype=float)+i*.01
    arrays['pTmiss']=np.full(4,.05)
    variables=(production.VARIABLES[0],)
    with patch.object(production,'ACTIVE_GLOBAL_CUTS',cfg),patch.object(production,'resolve_variable_branches',return_value={'Delta_phi':'Delta_phi'}),patch.object(production,'global_cut_required_branches',return_value=list(arrays)),patch.object(production,'iterate_rebuilt',return_value=iter([arrays])):
        result=production.load_selected_arrays_uproot('/fake','fa18_inb',production.TOPOLOGIES,'1 MB',variables)
    assert result['FD_FD']['Delta_phi'].tolist()==[0.,3.]  # W and blacklist cuts applied; no unwanted sector cut.
    rows=np.zeros((4,12));rows[:,0]=0;rows[:,8]=[1,0,0,1];rows[:,10]=1;rows[:,1]=np.arange(4)
    assert cuts.fit_arrays(rows,0)['Delta_phi'].tolist()==[0,3]
    assert cuts.fit_arrays(rows,0,'minimal')['Delta_phi'].tolist()==[0,1,2,3]
    print('PASS: production uproot loader applies W/blacklist/configured globals; no unconditional legacy sector cut; global/minimal cut-fit domains differ explicitly')

def check_mixture(tmp):
    cfg=next(p for p in production.PERIODS if p.key=='fa18_inb');blocks={}
    for top in production.TOPOLOGIES:
        for channel in ['dvcs','eppi0']:blocks[production.cut_json_key(channel,cfg,top)]=dict(data={},mc={})
    def make(n):
        a=np.zeros((n,12));a[:,0]=2;a[:,8]=1;a[:,10]=1;a[:,9]=.13;a[:,1:8]=.001;return a
    arrays=dict(data=make(100),mc=make(40),background=make(10),pi0_data=make(20),pi0_mc=make(10));prob=dict(mc=.5,background=.5,pi0_mc=.5)
    result=study.scales('fa18_inb',arrays,blocks,prob,1.,1.)
    scale=result[2];assert scale['control_scale']==1 and scale['background_scale']==1 and scale['signal_scale']==1
    key=production.cut_json_key('dvcs',cfg,production.TOPOLOGIES[2]);edges=np.linspace(-.02,.02,41)
    d,s,b=study.physical_prediction(arrays,2,5,edges,scale,prob,blocks[key],'global_pre')
    assert d.sum()==100 and s.sum()==80 and b.sum()==20 and np.allclose(s+b,d)
    rows=cuts.comparison_rows('fa18_inb','same',arrays,arrays,blocks,blocks,prob)
    study.add_background_proxy(rows,'fa18_inb',arrays,arrays,blocks,blocks,result,result,prob)
    assert all(np.isclose(r['subtracted_yield_acceptance_proxy_ratio'],1) for r in rows if r['baseline_mc_selected'])
    arrays['data'][:,8]=0;invalid=study.scales('fa18_inb',arrays,blocks,prob,1.,1.)[2]
    assert invalid['status']!='ok' and not np.isfinite(invalid['signal_scale'])
    arrays['data'][:,8]=1
    # Exact exported response association, including candidates with identical parameters.
    source=tmp/'source';source.mkdir();work=tmp/'scratch';(work/'precut_arrays').mkdir(parents=True);base.dump(source/'scratch_location.json',dict(directory=str(work)))
    vectors={('fa18_inb',c):[.01]+[0.]*19 for c in study.MODELS}
    jobs=[dict(period='fa18_inb',sample=s) for s in ['data','pi0_data']]
    jobs += [dict(period='fa18_inb',sample=s,parameters=v) for v in [[0.]*20,vectors['fa18_inb','current']] for s in ['mc','background','pi0_mc']]
    for i,job in enumerate(jobs):
        path=work/'precut_arrays'/str(i);base.dump(path.with_suffix('.json'),job);np.zeros((1,12)).tofile(path.with_suffix('.bin'));base.dump(Path(str(path)+'.bin.json'),dict(shape=[1,12]))
    loaded,files=study.load_exports(source,['fa18_inb'],vectors)
    assert len(loaded)==11 and ('fa18_inb','add_pt','pi0_mc') in loaded
    print('PASS: control-constrained summed MC, component weights/sampling, invalid background residual flags and reuse of exported candidate events')

def plot_check(tmp):
    cfg=next(p for p in production.PERIODS if p.key=='fa18_inb');blocks={}
    for top in production.TOPOLOGIES:
        for channel in ['dvcs','eppi0']:blocks[production.cut_json_key(channel,cfg,top)]=dict(data={},mc={})
    def make(n):
        rows=np.zeros((n,12));rows[:,0]=2;rows[:,8]=1;rows[:,10]=1;rows[:,9]=.13
        for i in range(7):rows[:,i+1]=np.linspace(-.01,.01,n)
        return rows
    a=dict(data=make(100),mc=make(40),background=make(10),pi0_data=make(20),pi0_mc=make(10));prob=dict(mc=.5,background=.5,pi0_mc=.5)
    anchor={c:study.scales('fa18_inb',a,blocks,prob,1.,1.) for c in study.MODELS}
    spec=dict(variables=[dict(name=n,label=n,bins=40,low=-.02,high=.02) for n in base.OBS])
    links,counts=study.render_physical('fa18_inb',spec,{c:a for c in study.MODELS},blocks,{c:blocks for c in study.MODELS},anchor,anchor,prob,tmp)
    assert len(links)==18 and len(counts)==168
    good=[r for r in counts if r['topology']=='CD_FT' and r['channel']=='DVCS'];assert all(np.isclose(r['signal_entries']+r['background_entries'],r['predicted_entries']) for r in good)
    assert not (tmp/'results').exists()
    print('PASS: physical summed/component plots, full-global/N−1/final views and separate two-photon controls')


def main():
    check_loader()
    with tempfile.TemporaryDirectory() as name:
        check_mixture(Path(name));plot_check(Path(name))
if __name__=='__main__':main()
