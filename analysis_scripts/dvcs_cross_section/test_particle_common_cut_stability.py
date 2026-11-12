#!/usr/bin/env python3
"""Deterministic physics/counting and cache-only report integration checks."""
import contextlib,copy,io,json,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_particle_common_cut_stability as study


def array(n,top=2,t=.13):
    a=np.zeros((n,12));a[:,0]=top;a[:,8]=a[:,10]=1;a[:,9]=t
    a[:,1:8]=.01;a[:,4]=np.linspace(0,.4,n);a[:,6]=np.linspace(-.04,.04,n)
    return a


def checks(tmp):
    blocks={'test':dict(data={'pTmiss':dict(cut_low=0,cut_high=.3,mode='upper_quantile')},mc={'pTmiss':dict(cut_low=0,cut_high=.1,mode='upper_quantile')})}
    new=study.common_blocks(blocks);assert new['test']['data']==new['test']['mc'];assert blocks['test']['mc']['pTmiss']['cut_high']==.1
    thresholds=np.array([.1,.2,.3]);a=array(101);selected=np.ones(101,dtype=bool)
    counts,nref,overlap=study.scan_counts(a[:,4],selected,thresholds,0,.2)
    assert counts.tolist()==[26,51,76] and nref==51 and overlap.tolist()==[26,51,51]
    prob=dict(mc=.5,background=.25,pi0_mc=.5)
    y,r,rr,error=study.stability([3*counts,counts,counts],[3*nref,nref,nref],[3*overlap,overlap,overlap],.25,prob,.001)
    assert np.allclose(r,1) and np.allclose(rr,1) and error[1]==0 and error[0]>0
    _,_,bad,_=study.stability([counts,counts,counts],[nref,nref,nref],[overlap]*3,2,prob)
    assert np.all(np.isnan(bad))
    _,_,empty,_=study.stability([np.zeros(3)]*3,[0,0,0],[np.zeros(3)]*3,1,prob)
    assert np.all(np.isnan(empty))
    # Missing DATA tails deliberately cause a measurable scale drift, even after anchor normalization.
    _,_,drift,_=study.stability([np.array([50,100,180]),np.array([5,10,15]),np.array([25,50,60])],[100,10,50],[np.array([50,100,100]),np.array([5,10,10]),np.array([25,50,50])],1,dict(mc=1,background=1))
    assert np.isclose(drift[1],1) and drift[2]>1.5
    d=np.array([10.,20.,5.]);s=np.array([5.,5.,0.]);b=np.array([1.,2.,0.])
    ss,bb=study.scaled_shape(d,s,b,2,1,dict(mc=1,background=1),True)
    assert np.isclose((ss+bb).sum(),1) and np.isclose(ss.sum()/bb.sum(),20/3)
    ss,bb=study.scaled_shape(d,s,b,2,1,dict(mc=1,background=1),False)
    assert np.isclose((ss+bb).sum(),23/35)
    cfg=next(p for p in study.cuts.production.PERIODS if p.key=='fa18_inb');allblocks={}
    for top in study.cuts.production.TOPOLOGIES:
        for channel in ('dvcs','eppi0'):
            allblocks[study.cuts.production.cut_json_key(channel,cfg,top)]=dict(data=copy.deepcopy(blocks['test']['data']),mc=copy.deepcopy(blocks['test']['mc']))
    common=study.common_blocks(allblocks);samples=dict(data=np.concatenate([a,a,a]),mc=a,background=a)
    cases={c:samples for c in study.SCENARIOS};sc=dict(signal_scale=2.,background_scale=1.,control_data_anchor=100,control_mc_anchor=100)
    anchors={c:{2:sc} for c in study.SCENARIOS};(tmp/'plots').mkdir()
    rows,links=study.render_case('fa18_inb',2,(.11,.15),cases,common,dict(saved_DATA=common,refitted_DATA=common),anchors,dict(mc=1,background=1,pi0_mc=1),tmp)
    assert len(links)==2 and all((tmp/p).is_file() for p in links)
    assert all(np.isclose(r['relative_signal_scale'],1) for r in rows if r['status']=='ok')
    assert all(r['conditional_stat_error']==0 for r in rows if r['threshold']==r['reference_threshold'] and r['status']=='ok')
    assert study.export_digest(tmp/links[0])==study.base.digest(tmp/links[0])
    print('PASS: common boundaries, exact nested counts, paired covariance/reference, shared control uncertainty, invalid subtraction, known tail drift, equal-area versus fixed normalization, physical mixture and report plots')


def orchestration(tmp):
    source=tmp/'followup'/'cut_redetermination';source.mkdir(parents=True);closed=tmp/'closure';closed.mkdir();run=tmp/'run';run.mkdir();training=tmp/'training';training.mkdir()
    period='fa18_inb';cfg=next(p for p in study.cuts.production.PERIODS if p.key==period)
    blocks={study.cuts.production.cut_json_key(ch,cfg,t):dict(data={},mc={}) for ch in ('dvcs','eppi0') for t in study.cuts.production.TOPOLOGIES}
    spec=dict(periods=[period],exclusivity=blocks)
    study.base.dump(source/'spec.json',spec)
    study.base.dump(run/'cache_manifest.json',dict(samples={period+'__'+s:dict(probability=1) for s in ('mc','background','pi0_mc')}))
    paths={}
    for c in study.closure.MODELS:
        p=training/(c+'.json');study.base.dump(p,dict(background_multiplier=1));paths[period,c]=p
    hashes={c:{period:study.base.digest(paths[period,c])} for c in study.closure.MODELS}
    old=dict(candidate_hashes=hashes,response_cache_sha256=study.base.digest(run/'cache_manifest.json'),arguments=dict(response_run=str(run),training_study=str(training)))
    study.base.dump(source/'manifest.json',old)
    fake=tmp/'events.bin';fake.write_bytes(b'event fixture')
    study.base.dump(closed/'manifest.json',dict(source_cut_manifest_sha256=study.base.digest(source/'manifest.json'),candidate_hashes=hashes,export_hashes={str(fake):study.export_digest(fake)}))
    folder=closed/'global_cut_refits'/period/'add_pt';folder.mkdir(parents=True);study.base.dump(folder/'combined_cuts_95.json',blocks)
    a=array(100);samples=dict(data=np.concatenate([a,a]),mc=a,background=a,pi0_data=a,pi0_mc=a)
    exports={(period,'shared' if s in ('data','pi0_data') else c,s):v for c in study.SCENARIOS for s,v in samples.items()}
    output=tmp/'report';argv=['--followup-study',str(source.parent),'--closure-study',str(closed),'--output',str(output),'--topologies','CD_FT']
    with patch.object(study.cuts,'candidate_paths',return_value=paths),patch.object(study.cuts,'load_candidate_parameters',return_value=[0.]*20),patch.object(study.closure,'load_exports',return_value=(exports,{k:fake for k in exports})),patch.object(study,'BINS',((.11,.15),)):
        study.main(argv)
        assert (output/'REPORT.html').is_file() and (output/'cut_stability.csv').is_file()
        manifest=json.loads((output/'manifest.json').read_text());assert manifest['candidate_hashes']==hashes
        try:
            with contextlib.redirect_stderr(io.StringIO()):study.main(argv)
        except SystemExit as e:assert e.code==2
        else:raise AssertionError('Overwrote output')
        old['candidate_hashes']['add_pt'][period]='changed';study.base.dump(source/'manifest.json',old)
        try:
            with contextlib.redirect_stderr(io.StringIO()):study.main(['--followup-study',str(source.parent),'--closure-study',str(closed),'--output',str(tmp/'bad')])
        except SystemExit:pass
        else:raise AssertionError('Accepted changed provenance')
    print('PASS: cache-only end-to-end driver, source/candidate provenance, manifest/CSV/HTML, overwrite guard')

if __name__=='__main__':
    with tempfile.TemporaryDirectory() as folder:
        checks(Path(folder));orchestration(Path(folder))
