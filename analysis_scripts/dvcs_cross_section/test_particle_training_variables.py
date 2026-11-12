#!/usr/bin/env python3
"""Check weighting, training-only optimization and the comparison reports."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import tempfile
import numpy as np
import study_particle_training_variables as study
base=study.base

def weights_and_scores():
    assert study.parse_weights([])=={base.OBS[v]:1. for v in study.SCORED}
    for invalid in [['Mx2_2=-1'],['theta=1'],['pTmiss=nan'],['pTmiss=2','pTmiss=3']]:
        try:study.parse_weights(invalid)
        except ValueError:pass
        else:raise AssertionError('Invalid weight accepted')
    period='fa18_inb';spec=dict(periods=[period]);weights=study.parse_weights(['pTmiss=2'])
    rows=[dict(variable=base.OBS[v],score=float(v+1),half=half,kind='old') for top in range(3) for v in study.SCORED for half in (0,1)]
    with patch.object(base.Engine,'score_one',return_value=({},rows)):
        engine=study.TrainingEngine(None,spec,{}, {},{},base.TOPS,100,training=study.CONFIGURATIONS['add_pt'],weights=weights)
        score,detail=engine.score_one({},0,1.)
        assert np.isclose(score['train'],(1+2*4+5+6)/5)
        assert score['training_projections']==12
        assert np.isclose(score['heldout_observables'],(3+7)/2)
        for row in rows:
            if row['variable']=='Mx2_2':row['score']=1e9
        assert engine.score_one({},0,1.)[0]['train']==score['train']
        rows[0]['score']=float('nan');assert not np.isfinite(engine.score_one({},0,1.)[0]['train'])
    assert base.TRAIN_OBS==[0,4,5] and base.VALID_OBS==[2,3,6]
    print('PASS: explicit positive weights, strict projection coverage, isolated training configurations and reserved-variable exclusion')

class Toy:
    def __init__(self,config):
        self.spec=dict(periods=['fa18_inb'],starting_candidate_sha256='fixed');self.training=study.CONFIGURATIONS[config];self.weights=study.parse_weights([])
    def evaluate(self,parameters,label,detailed=False):
        points=[tuple([0.]*20)]
        for point in parameters:
            if tuple(point) not in points:points.append(tuple(point))
        scores=[];records=[]
        for i,p in enumerate(points):
            train=(p[3]+.02)**2+sum(x*x for k,x in enumerate(p) if k!=3)
            reserved=(p[3]-.08)**2  # Opposing preference must never select parameters.
            scores.append(dict(train=train,train_shape=train,validation_mc=train,heldout_observables=reserved,background_multiplier=1.))
            for top in base.TOPS:
                for v in study.SCORED:
                    for half in (0,1):records.append(dict(scenario=i,period='fa18_inb',topology=top,variable=base.OBS[v],half=half,kind='train' if v in self.training else 'heldout_observables',score=reserved+.01 if v==6 else train+.01))
        selected={k:np.full((len(points),2,3),100) for k in ['fa18_inb__mc','fa18_inb__background']}
        return points,scores,records,{},selected
    def weighted(self,*args,**kwargs):return np.ones(12),np.ones(12)

def optimization_and_reports(tmp):
    start=[0.]*20;start[3]=-.015;args=SimpleNamespace(rounds=2,scan_range='wide');results={}
    for config in study.CONFIGURATIONS:
        folder=tmp/config;folder.mkdir();result=study.fit_local(Toy(config),start,args,folder,False);results['fa18_inb',config]=result
        assert result['parameters'][3]==-.02
        assert result['candidate']['reserved_observables'][-1]=='Mx2_2'
        study.save_checkpoint(folder,result);restored=study.load_checkpoint(folder)
        assert restored['parameters']==result['parameters']
        assert all(np.array_equal(value,restored['curves'][key]) for key,value in result['curves'].items())
    assert study.nearest_trials({0:[0,.01,.02]},[.015],0)==[.01,.0125,.015,.0175,.02]
    spec=dict(periods=['fa18_inb'],variables=[dict(name=n,label=n,bins=10,low=-1,high=1) for n in base.OBS]);data={'fa18_inb__data':np.ones((3,1,7,12))}
    study.report(tmp,spec,data,results,study.parse_weights([]))
    assert (tmp/'REPORT.html').is_file() and (tmp/'observable_improvements.png').is_file()
    assert len(list((tmp/'plots').glob('*.png')))==3 and not (tmp/'results').exists()
    assert (tmp/'parameter_stability.csv').is_file() and (tmp/'observable_comparison.csv').is_file()
    print('PASS: fixed current-fit reference, warm local refinement uses training only, per-observable/stability tables and single-tree comparison plots')

def orchestration(tmp):
    root=tmp/'orchestration';reference=root/'reference';(reference/'backend').mkdir(parents=True)
    spec=dict(periods=['fa18_inb'],variables=[dict(name=n,label=n,bins=10,low=-1,high=1) for n in base.OBS],input_hashes={})
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        (reference/'backend'/name).write_text('saved source');spec['input_hashes'][key]=base.digest(reference/'backend'/name)
    np.savez(reference/'histograms.npz',fa18_inb__data=np.ones((3,1,7,12)))
    base.dump(reference/'manifest.json',dict(spec=spec,flows={},histogram_sha256=base.digest(reference/'histograms.npz')))
    base.table(reference/'normalizations.csv',[dict(period='fa18_inb',topology=t,status='ok') for t in base.TOPS])
    run=root/'run';(run/'with_bias/fa18_inb').mkdir(parents=True);start=[0.]*20;start[3]=-.015
    base.dump(run/'with_bias/fa18_inb/candidate_parameters.json',dict(period='fa18_inb',bias_mode='fit',parameters=dict(zip(base.PARAMS,start))))
    cache={}
    for sample in ['mc','background','pi0_mc']:
        path=root/(sample+'.bin');path.write_bytes(b'checked cache');cache['fa18_inb__'+sample]=dict(file=str(path),sha256=base.digest(path),probability=1.)
    base.dump(run/'with_bias/cache_manifest.json',dict(version=base.VERSION,reference_sha256=base.digest(reference/'manifest.json'),samples=cache))
    output=root/'study';cmd=['--reference-study',str(reference),'--response-run',str(run),'--output',str(output),'--scratch-root',str(root/'scratch'),'--workers','2','--rounds','1']
    budgets=[];real_fit=study.fit_local
    def fake_fit(engine,start,args,folder,reference):
        budgets.append(id(engine.pool.budget));config=next(c for c,v in study.CONFIGURATIONS.items() if tuple(v)==tuple(engine.training))
        return real_fit(Toy(config),start,args,folder,False)
    with patch.object(base,'compile_worker',return_value=Path('/unused/worker')) as compile_once,patch.object(study,'fit_local',fake_fit),patch.object(study,'report') as render:
        study.main(cmd);assert compile_once.call_count==1 and render.call_count==1
    assert len(budgets)==3 and len(set(budgets))==1 and not (output/'results').exists()
    with patch.object(base,'compile_worker',side_effect=AssertionError('Completed resume must not compile')),patch.object(study,'fit_local',side_effect=AssertionError('Completed fits must not repeat')),patch.object(study,'report') as render:
        study.main(cmd+['--resume']);assert render.call_count==1
    try:study.main(cmd[:-1]+['2','--resume'])
    except ValueError as error:assert 'configuration changed' in str(error)
    else:raise AssertionError('Resume accepted changed fit settings')
    print('PASS: paired-run input discovery, one compilation/shared budget, cache-only workflow and checked resume without refitting')


def main():
    weights_and_scores()
    with tempfile.TemporaryDirectory() as d:
        optimization_and_reports(Path(d));orchestration(Path(d))
if __name__=='__main__':main()
