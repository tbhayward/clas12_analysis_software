#!/usr/bin/env python3
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_particle_response as base
import study_compare_control_training as compare


def main():
    edges=np.linspace(-1,1,80);x=(edges[:-1]+edges[1:])/2;good=np.r_[0,100*np.exp(-x*x/.08),0];bad=np.roll(good,8);n=len(good)
    spec=dict(periods=['fa18_inb'],variables=[dict(bins=n-2) for _ in base.OBS],joint_control_fit=True,joint_control_weight=0)
    data={'fa18_inb__pi0_data':np.tile(good,(3,1,7,1))};engine=base.Engine(None,spec,{},data,{},['CD_FT'],10);hist=np.tile(good,(2,2,3,1,7,1));hist[1,0,2,0,0]=bad
    summary=dict(train=12.,validation_mc=9.,heldout_observables=8.,training_projections=1)
    with patch.object(base,'TRAIN_OBS',[0]),patch.object(base,'VALID_OBS',[6]):
        first,_=engine.joint_control_score({'fa18_inb__pi0_mc':hist},0,summary.copy(),[]);second,rows=engine.joint_control_score({'fa18_inb__pi0_mc':hist},1,summary.copy(),[])
        assert first['train']==second['train']==12 and second['epgammagamma_train']>first['epgammagamma_train']
        spec['joint_control_weight']=.5
        joint,_=engine.joint_control_score({'fa18_inb__pi0_mc':hist},1,summary.copy(),[])
        assert np.isclose(joint['train'],.5*(12+second['epgammagamma_train']))
        hist[1,1,2,0,0]=bad*2;held,_=engine.joint_control_score({'fa18_inb__pi0_mc':hist},1,summary.copy(),[])
        assert held['train']==joint['train'] and held['epgammagamma_validation_mc']>joint['epgammagamma_validation_mc']
        assert all(r['channel']=='epgammagamma' for r in rows)
    print('PASS: epgamma-only objective ignores two-photon shape while retaining its validation scores; joint objective includes it; MC validation cannot select training')
    # Report integration with both candidates and both DATA/MC shape channels.
    v=[dict(name=name,label=name,bins=n-2,edges=edges.tolist()) for name in base.OBS]
    engine.spec['variables']=v;engine.data['fa18_inb__data']=engine.data['fa18_inb__pi0_data'];hist=np.tile(good,(3,2,3,1,7,1));h={f'fa18_inb__{s}':hist for s in ('mc','background','pi0_mc')}
    def score(h,i,m):
        r={f'{ch}_{key}':float(i) for ch in ('epgamma','epgammagamma') for key in ('train','validation_mc','heldout_observables')}
        return r,[]
    with tempfile.TemporaryDirectory() as d:
        out=Path(d);paths={name:out/(name+'.json') for name in ('epgamma_only','joint')}
        for path in paths.values():base.dump(path,dict(background_multiplier=1))
        with patch.object(engine,'score_one',side_effect=score),patch.object(engine,'weighted',return_value=(good,good)):
            rows,summaries,links=compare.comparison_plot('fa18_inb',engine,h,{},paths,out)
        assert len(links)==3 and all((out/l).is_file() for l in links)
        assert {r['channel'] for r in summaries}=={'epgamma','epgammagamma'} and {r['fit'] for r in summaries}==set(paths)
    print('PASS: both fitted candidates compared to two-photon DATA/AAOGEN and epgamma DATA/mixture, channel score bars and paired shape plots')

if __name__=='__main__':main()
