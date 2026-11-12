#!/usr/bin/env python3
"""Check candidate provenance, stratified comparison, plots and paired acceptance."""
import json,tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_particle_training_followup as follow
base=follow.base;cuts=follow.cuts

def checks(tmp):
    training=tmp/'training';paths={}
    for name,bias in [('current',-.01),('add_pt',-.02)]:
        folder=training/'fa18_inb'/name;folder.mkdir(parents=True)
        parameters=dict.fromkeys(base.PARAMS,0.);parameters['electron_FD_p_bias_fraction']=bias
        path=folder/'candidate_parameters.json';base.dump(path,dict(period='fa18_inb',bias_mode='fit',parameters=parameters,background_multiplier=1.))
        base.dump(folder/'completed_fit.json',dict(candidate_sha256=base.digest(path)));paths[name]=path
    discovered=cuts.candidate_paths({c:tmp for c in paths},['fa18_inb'],training)
    assert discovered['fa18_inb','add_pt']==paths['add_pt']
    assert cuts.load_candidate_parameters(paths['add_pt'])[3]==-.02
    old=paths['current'].read_text();paths['current'].write_text(old+' ')
    try:cuts.candidate_paths({c:tmp for c in paths},['fa18_inb'],training)
    except ValueError:pass
    else:raise AssertionError('Changed candidate was accepted')
    paths['current'].write_text(old)
    spec=dict(periods=['fa18_inb'],variables=[dict(name=n,label=n,bins=40,low=-2,high=2) for n in base.OBS],strata=[dict(axis='inclusive',label='Inclusive'),dict(axis='p1_p',label='0.5–1'),dict(axis='joint_proton',label='Joint cell')])
    x=np.linspace(-2,2,40);d=np.r_[0,200*np.exp(-x*x),0];m=np.r_[0,200*np.exp(-(x-.3)**2),0]
    data={'fa18_inb__data':np.tile(d,(3,3,7,1))}
    class ToyEngine:
        def __init__(self,*a,**k):pass
        def evaluate(self,parameters,*a,**k):
            points=[tuple([0.]*20),*[tuple(p) for p in parameters]]
            hist={s:np.full((3,2,3,3,7,42),10) for s in ['fa18_inb__mc','fa18_inb__background']}
            return points,[],[],hist,{}
        def weighted(self,h,p,i,half,top,g,v,*a,**k):return (d if i==2 else m).copy(),np.ones(42)*20
    output=tmp/'study';output.mkdir();work=tmp/'work';work.mkdir()
    with patch.object(base,'Engine',ToyEngine):rows,examples=follow.scan_period('fa18_inb',spec,{},data,{},paths,Path('/unused'),work,output,base.Budget(2),2,10)
    assert len(rows)==63 and any(r['axis']=='joint_proton' for r in rows)
    assert all(r['improvement_percent']>99 for r in rows)
    follow.render(output,spec,rows,{'fa18_inb':examples})
    assert (output/'REPORT.html').exists() and len(list((output/'plots').glob('*.png')))==6 and not (output/'results').exists()
    scores,gain,status=follow.cell_metrics(d,[m,d],[np.ones(42),np.ones(42)],m,np.ones(42),1,100)
    assert status!='usable descriptive comparison' and not np.isfinite(gain)
    # Compare candidate acceptance directly against the equally refined current response.
    a=np.zeros((10,12));a[:,0]=2;a[:,8]=1;a[:,9]=.13;a[:,6]=np.linspace(-.02,.02,10)
    p=next(p for p in cuts.production.PERIODS if p.key=='fa18_inb');blocks={}
    for top in cuts.production.TOPOLOGIES:blocks[cuts.production.cut_json_key('dvcs',p,top)]=dict(data={},mc={})
    arrays={'mc':a,'data':a};z=cuts.comparison_rows('fa18_inb','add_pt',arrays,arrays,blocks,blocks,dict(mc=1))
    assert all(r['relative_acceptance']==1 for r in z if r['baseline_mc_selected'])
    print('PASS: checked training candidate discovery, both-response subdivision/joint diagnostics, low-stat flags, plots and paired-current acceptance comparison')

def main():
    with tempfile.TemporaryDirectory() as d:checks(Path(d))
if __name__=='__main__':main()
