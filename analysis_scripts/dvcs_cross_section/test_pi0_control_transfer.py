#!/usr/bin/env python3
import contextlib,io,json,shutil,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_pi0_control_transfer as study

def arrays():
    a=np.zeros((120,20));a[:,0]=np.arange(120)%3;a[:,8]=1;a[:,9]=.2;a[:,12]=.25;a[:,13]=2.5;a[:,14:]=.5
    return dict(pi0_data=np.concatenate([a,a]),pi0_mc=a)

def unit(tmp):
    assert study.control_ratio(100,50,1,.5)[0]==1
    assert study.control_ratio(0,50,1,1)==(0.,0.)
    assert np.isnan(study.control_ratio(100,0,1,1)[0])
    period='fa18_inb';cfg=next(p for p in study.cuts.production.PERIODS if p.key==period)
    blocks={study.cuts.production.cut_json_key('eppi0',cfg,t):dict(data={},mc={}) for t in study.cuts.production.TOPOLOGIES}
    out=tmp/'report';out.mkdir();a=arrays()
    study.analyze_controls({period:a},dict(exclusivity=blocks),{period:(dict(control_scale=2),dict(pi0_mc=1))},[period],out)
    import csv
    rows=list(csv.DictReader((out/'control_ratios.csv').open()));assert all(np.isclose(float(r['data_over_prediction']),1) for r in rows if int(r['raw_mc']))
    assert len(list((out/'plots').rglob('*.png')))==4 and (out/'anchor_check.csv').is_file()
    print('PASS: sampling-aware differential control ratios, empty flags, total-anchor closure, all-topology selection, CSV/HTML/plots')

def driver(tmp):
    root=tmp/'driver';root.mkdir();origin=root/'common';origin.mkdir();follow=root/'followup'/'cut_redetermination';follow.mkdir(parents=True);closed=root/'closed';closed.mkdir();ref=root/'reference';(ref/'backend').mkdir(parents=True);run=root/'response';run.mkdir();training=root/'training';training.mkdir()
    period='fa18_inb';data=root/'data.root';data.write_bytes(b'fixture');st=data.stat()
    for n in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(study.base.SOURCE/n,ref/'backend'/n)
    study.base.dump(ref/'manifest.json',dict(inventory=[dict(period=period,sample='pi0_data',path=str(data),size=st.st_size,mtime_ns=st.st_mtime_ns)]))
    cachefile=root/'cache.bin';cachefile.write_bytes(b'cache')
    study.base.dump(run/'cache_manifest.json',dict(samples={period+'__'+s:dict(probability=1,file=str(cachefile),sha256=study.base.digest(cachefile)) for s in ('pi0_mc',)}))
    candidate=training/'candidate.json';study.base.dump(candidate,dict(background_multiplier=1));hashes={'add_pt':{period:study.base.digest(candidate)}}
    cfg=next(p for p in study.cuts.production.PERIODS if p.key==period);key=study.cuts.production.cut_json_key('dvcs',cfg,study.cuts.production.TOPOLOGIES[0]);blocks={key:dict(data={},mc={})}
    study.base.dump(follow/'spec.json',dict(periods=[period],exclusivity=blocks))
    study.base.dump(follow/'manifest.json',dict(candidate_hashes=hashes,reference_sha256=study.base.digest(ref/'manifest.json'),response_cache_sha256=study.base.digest(run/'cache_manifest.json'),arguments=dict(reference_study=str(ref),response_run=str(run),training_study=str(training))))
    study.base.dump(closed/'manifest.json',{})
    study.base.dump(origin/'manifest.json',dict(arguments=dict(followup_study=str(follow.parent),closure_study=str(closed)),source_sha256=study.base.digest(follow/'manifest.json'),closure_sha256=study.base.digest(closed/'manifest.json'),candidate_hashes=hashes))
    norm=dict(period=period,scenario='add_pt',topology='FD_FD',signal_scale=1.,background_scale=1.,control_scale=1.,control_data_anchor=1000,control_mc_anchor=1000,status='ok');study.base.table(origin/'normalizations.csv',[norm])
    class FakePool:
        def __init__(self,worker,spec,work,workers):self.work=work
        def run(self,mode,jobs,label):
            assert mode=='cut-arrays' and len(jobs)==2
            result=[]
            for i,job in enumerate(jobs):
                assert job['export_kinematics'] and 'export_topology' not in job and job['export_global_only']
                a=arrays()[job['sample']];path=self.work/(str(i)+'.bin');a.tofile(path);study.base.dump(Path(str(path)+'.json'),dict(shape=list(a.shape)));result.append((i,job,path))
            return result
    output=root/'new';argv=['--common-study',str(origin),'--scratch-root',str(root/'scratch'),'--output',str(output)]
    with patch.object(study.cuts,'candidate_paths',return_value={(period,'add_pt'):candidate}),patch.object(study.cuts,'load_candidate_parameters',return_value=[0.]*20),patch.object(study.base,'compile_worker',return_value=root/'worker'),patch.object(study.base,'Pool',FakePool),patch.object(study,'analyze_controls') as analyze:
        study.main(argv);assert analyze.call_count==1 and (output/'export_manifest.json').is_file()
        with patch.object(study.base,'compile_worker',side_effect=AssertionError('Replayed exports must not compile')):
            study.main(['--common-study',str(origin),'--reuse-exports',str(output),'--output',str(root/'replay')])
        norm['background_scale']=2.;study.base.table(origin/'normalizations.csv',[norm])
        try:
            with contextlib.redirect_stderr(io.StringIO()):study.main(['--common-study',str(origin),'--reuse-exports',str(output),'--output',str(root/'bad')])
        except SystemExit as error:assert error.code==2
        else:raise AssertionError('Accepted changed transfer normalization')
    print('PASS: new-export driver, control-only worker jobs/shared budget, input provenance and export checks, cache-only replay without compilation, changed-normalization rejection')



if __name__=='__main__':
    with tempfile.TemporaryDirectory() as d:
        unit(Path(d));driver(Path(d))
