#!/usr/bin/env python3
"""Shell algebra, partition coverage, plots, and optional real C++ export smoke."""
import argparse,contextlib,copy,io,json,shutil,subprocess,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_fd_background_transfer as study


def arrays():
    a=np.zeros((300,20));a[:,8]=1;a[:,9]=.2;a[:,10]=1
    a[:,1:8]=.01;a[:,4]=np.linspace(.001,.399,300);a[:,6]=np.linspace(-.04,.04,300)
    for col,factor in study.AXES.values():a[:,col]=np.linspace(.1,1,300)/factor
    return dict(data=np.concatenate([a,a]),mc=a,background=a)


def unit(tmp):
    a=arrays();block=dict(data={},mc={});prob=dict(mc=1,background=1,pi0_mc=1);sc=dict(signal_scale=1.,background_scale=1.,control_data_anchor=1000,control_mc_anchor=1000)
    parts=study.partitions(a['data'],'proton_p');sub=[study.subset(a,'proton_p',*part) for part in parts]
    assert len(parts)==3 and sum(len(s['data']) for s in sub)==len(a['data']) and sum(len(s['mc']) for s in sub)==len(a['mc'])
    rows=study.shell_rows(a,block,[0,.1,.2,.4],sc,prob)
    assert sum(r['data_count'] for r in rows)==600 and sum(r['background_mc_count'] for r in rows)==300
    assert all(r['subtracted_data']==r['signal_mc_count'] for r in rows)
    bad=study.shell_rows(a,block,[0,.1,.2,.4],dict(sc,background_scale=3),prob)
    assert all(r['status']=='negative residual' for r in bad) and all(r['subtracted_data']<0 for r in bad)
    cfg=next(p for p in study.cuts.production.PERIODS if p.key=='fa18_inb');key=study.cuts.production.cut_json_key('dvcs',cfg,study.cuts.production.TOPOLOGIES[0])
    blocks={key:block};spec=dict(exclusivity=blocks);closed=tmp/'closed';folder=closed/'global_cut_refits/fa18_inb/add_pt';folder.mkdir(parents=True);study.base.dump(folder/'combined_cuts_95.json',blocks)
    out=tmp/'report';out.mkdir()
    with patch.object(study,'AXES',{'proton_p':study.AXES['proton_p']}),patch.object(study.common,'BINS',((.11,.15),(.15,.25))):
        study.analyze({'fa18_inb':a},spec,{'fa18_inb':(sc,prob)},['fa18_inb'],'refitted_DATA',closed,out)
    assert (out/'REPORT.html').is_file() and (out/'shell_budgets.csv').is_file()
    assert len(list((out/'plots').rglob('*.png')))==3
    assert 'No DATA in this cell' in (out/'REPORT.html').read_text()
    print('PASS: non-overlapping matched DATA-tercile edges, sampled background shell accounting, negative residuals retained, fixed-response kinematic plots/CSV/HTML and empty-cell reporting')


def driver(tmp):
    root=tmp/'driver';root.mkdir();origin=root/'common';origin.mkdir();follow=root/'followup'/'cut_redetermination';follow.mkdir(parents=True);closed=root/'closed';closed.mkdir();ref=root/'reference';(ref/'backend').mkdir(parents=True);run=root/'response';run.mkdir();training=root/'training';training.mkdir()
    period='fa18_inb';data=root/'data.root';data.write_bytes(b'fixture');st=data.stat()
    for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(study.base.SOURCE/n,ref/'backend'/n)
    study.base.dump(ref/'manifest.json',dict(inventory=[dict(period=period,sample='data',path=str(data),size=st.st_size,mtime_ns=st.st_mtime_ns)]))
    cachefile=root/'cache.bin';cachefile.write_bytes(b'cache')
    study.base.dump(run/'cache_manifest.json',dict(samples={period+'__'+s:dict(probability=1,file=str(cachefile),sha256=study.base.digest(cachefile)) for s in ('mc','background','pi0_mc')}))
    candidate=training/'candidate.json';study.base.dump(candidate,dict(background_multiplier=1));hashes={'add_pt':{period:study.base.digest(candidate)}}
    cfg=next(p for p in study.cuts.production.PERIODS if p.key==period);key=study.cuts.production.cut_json_key('dvcs',cfg,study.cuts.production.TOPOLOGIES[0]);blocks={key:dict(data={},mc={})}
    study.base.dump(follow/'spec.json',dict(periods=[period],exclusivity=blocks))
    study.base.dump(follow/'manifest.json',dict(candidate_hashes=hashes,reference_sha256=study.base.digest(ref/'manifest.json'),response_cache_sha256=study.base.digest(run/'cache_manifest.json'),arguments=dict(reference_study=str(ref),response_run=str(run),training_study=str(training))))
    study.base.dump(closed/'manifest.json',{})
    study.base.dump(origin/'manifest.json',dict(arguments=dict(followup_study=str(follow.parent),closure_study=str(closed)),source_sha256=study.base.digest(follow/'manifest.json'),closure_sha256=study.base.digest(closed/'manifest.json'),candidate_hashes=hashes))
    norm=dict(period=period,scenario='add_pt',topology='FD_FD',signal_scale=1.,background_scale=1.,control_data_anchor=1000,control_mc_anchor=1000,status='ok');study.base.table(origin/'normalizations.csv',[norm])
    class FakePool:
        def __init__(self,worker,spec,work,workers):self.work=work
        def run(self,mode,jobs,label):
            assert mode=='cut-arrays' and len(jobs)==3
            result=[]
            for i,job in enumerate(jobs):
                assert job['export_kinematics'] and job['export_topology']==0 and job['export_global_only']
                a=arrays()[job['sample']];path=self.work/(str(i)+'.bin');a.tofile(path);study.base.dump(Path(str(path)+'.json'),dict(shape=list(a.shape)));result.append((i,job,path))
            return result
    output=root/'new';argv=['--common-study',str(origin),'--scratch-root',str(root/'scratch'),'--output',str(output)]
    with patch.object(study.cuts,'candidate_paths',return_value={(period,'add_pt'):candidate}),patch.object(study.cuts,'load_candidate_parameters',return_value=[0.]*20),patch.object(study.base,'compile_worker',return_value=root/'worker'),patch.object(study.base,'Pool',FakePool),patch.object(study,'analyze') as analyze:
        study.main(argv);assert analyze.call_count==1 and (output/'export_manifest.json').is_file()
        with patch.object(study.base,'compile_worker',side_effect=AssertionError('Replayed exports must not compile')):
            study.main(['--common-study',str(origin),'--reuse-exports',str(output),'--output',str(root/'replay')])
        norm['background_scale']=2.;study.base.table(origin/'normalizations.csv',[norm])
        try:
            with contextlib.redirect_stderr(io.StringIO()):study.main(['--common-study',str(origin),'--reuse-exports',str(output),'--output',str(root/'bad')])
        except SystemExit as error:assert error.code==2
        else:raise AssertionError('Accepted changed transfer normalization')
    print('PASS: new-export driver, FD-only worker jobs/shared budget, input provenance and export checks, cache-only replay without compilation, changed-normalization rejection')


def root_smoke(tmp):
    source=study.base.SOURCE;ref=tmp/'reference';(ref/'backend').mkdir(parents=True);src=tmp/'worker_sources';src.mkdir();out=tmp/'worker';out.mkdir()
    for name in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(source/name,ref/'backend'/name)
    for name in ('exclusivity_resolution_worker.cpp','exclusivity_resolution_core.h','particle_response_core.h','particle_response_worker.cpp','particle_control_core.h','particle_control_worker.h'):shutil.copy2(source/name,src/name)
    shutil.copy2(Path(__file__).parent/'particle_response_cut_arrays.h',src/'particle_response_cut_arrays.h')
    with patch.object(study.base,'SOURCE',src):worker=study.base.compile_worker(out,ref)
    # Read a complete fixture config shipped in a completed study if provided via local source,
    # otherwise generate every config key from the production struct defaults with C++.
    text=(src/'exclusivity_resolution_worker.cpp').read_text();import re
    names=[n for n in re.findall(r'LOAD\((\w+)\)',text) if n!='name']
    generator=tmp/'config.cpp';generator.write_text('#include "global_cuts.h"\n#include <nlohmann/json.hpp>\n#include <fstream>\nint main(int argc,char**argv){GlobalCutConfig c;nlohmann::json j;'+''.join('j["'+n+'"]=c.'+n+';' for n in names)+'std::ofstream(argv[1])<<j;}')
    flags=subprocess.check_output(['root-config','--cflags'],text=True).split();binary=tmp/'config'
    subprocess.run(['c++',*flags,'-std=c++17','-I'+str(source),str(generator),'-o',str(binary)],check=True);subprocess.run([str(binary),str(tmp/'config.json')],check=True)
    cfg=json.loads((tmp/'config.json').read_text());cfg.update(enable_auxiliary_fiducial_cuts=False,enable_analysis_phase_space_cuts=False,enable_reconstructed_electron_p_cut=False,enable_pTmiss_cut=False,enable_topology_filter=False,enable_electron_fd_sector_filter=False,enable_proton_fd_sector_filter=False,enable_proton_cd_sector_filter=False,enable_photon_fd_sector_filter=False,excluded_runs=[])
    spec=tmp/'spec.json';study.base.dump(spec,dict(global_cuts=cfg))
    dtype=np.dtype([('v','<f8',(24,)),('key','<u8')]);records=np.zeros(3,dtype=dtype);assert dtype.itemsize==200
    v=records['v'];v[:,0]=.2;v[:,1]=2.2;v[:,2]=-.3;v[:,3:5]=1;v[:,5]=6;v[:,6]=.2;v[:,7]=0;v[:,8]=.8;v[:,9]=.3;v[:,10]=2;v[:,11]=3;v[:,12]=.25;v[:,13]=4;v[:,14]=10;v[:,15]=.1;v[:,16:23]=.01;v[:,23]=5032
    v[1,3]=2;v[2,14]=0;path=tmp/'cache.bin';records.tofile(path)
    outputs=[]
    for mode in ('legacy','extended','filtered','biased'):
        pars=[0.]*20
        if mode=='biased':pars[7]=.01
        job=dict(period='fa18_inb',sample='mc',cache=str(path),parameters=pars)
        if mode!='legacy':job['export_kinematics']=True
        if mode in ('filtered','biased'):job.update(export_topology=0,export_global_only=True)
        jp=tmp/(mode+'.json');target=tmp/(mode+'.bin');study.base.dump(jp,job)
        subprocess.run([str(worker),'cut-arrays',str(spec),str(jp),str(target)],check=True)
        shape=json.loads(Path(str(target)+'.json').read_text())['shape'];outputs.append(np.fromfile(target,dtype='<f8').reshape(shape))
    legacy,extended,filtered,biased=outputs
    assert legacy.shape==(3,12) and extended.shape==(3,20) and np.array_equal(legacy,extended[:,:12])
    assert filtered.shape==(1,20) and np.allclose(filtered[0,12:],v[0,[0,1,5,6,8,9,11,12]])
    assert np.isclose(biased[0,16],.808) and not np.isclose(biased[0,9],filtered[0,9])
    print('PASS: native production worker compilation, legacy 12-column compatibility, exact extended coordinate mapping, FD/global filtering and post-response momentum/|t| export')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root-smoke',action='store_true');args=p.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        tmp=Path(directory);unit(tmp);driver(tmp)
        if args.root_smoke:root_smoke(tmp)
