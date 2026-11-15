#!/usr/bin/env python3
"""Synthetic tests of cached migration, actual halves and photon-role invariance."""
import json,tempfile
from pathlib import Path
import numpy as np
import study_selected_response_checks as check


def main():
    a=np.zeros((5,12));a[:,8]=1;a[:,9]=[.12,.2,.12,.12,.2];a[:,11]=[0,1,0,1,1];a[3,8]=0
    b=a.copy();b[:,9]=[.2,.12,.12,.12,.2];b[2,8]=0;b[3,8]=1
    w=dict(data={},mc={});m=check.migration(a,b,w,w,0)
    assert m[1,2]==1 and m[2,1]==1 and m[1,8]==1 and m[8,1]==1 and m[2,2]==1 and m.sum()==5
    bad=b.copy();bad[0,11]=1
    try:check.migration(a,bad,w,w,0)
    except ValueError:pass
    else:raise AssertionError('Alignment guard failed')
    variables=[dict(name=n,bins=4,low=-1,high=1) for n in check.base.OBS];spec=dict(variables=variables,mass_variable=dict(bins=4,low=.1,high=.17));windows=dict(dvcs=[w]*3,control=[w]*4)
    a[:,1:8]=0;a[0,1]=2
    h,s=check.sample_histograms(a,windows,spec)
    assert s[:,0].tolist()==[2,2] and h[0,0,0,-1]==1 and h[1,0,0].sum()==2
    # N-minus-one: displaying the cut variable omits it, other projections retain it.
    windows['dvcs'][0]=dict(data={'Delta_phi':dict(cut_low=-.5,cut_high=.5)},mc={})
    h,s=check.sample_histograms(a,windows,spec)
    assert s[:,0].tolist()==[1,2] and h[0,0,0].sum()==2 and h[0,0,1].sum()==1
    c=np.zeros((3,29));c[:,20:28]=[3,.2,0,1,4,.3,0,1];c[:,28]=[.13,.135,.14];c[2,20]=np.nan
    x,mass=check.photon_entries(c,np.ones(3,dtype=bool),'momentum');swap=c.copy();swap[:,20:24]=c[:,24:28];swap[:,24:28]=c[:,20:24]
    xx,mm=check.photon_entries(swap,np.ones(3,dtype=bool),'momentum')
    assert sorted(zip(x,mass))==sorted(zip(xx,mm)) and len(x)==5
    print('PASS: bin migration, selection entry/exit, alignment guard, real MC halves, overflow, N-minus-one and photon-role symmetry')

def integration():
    c=check
    with tempfile.TemporaryDirectory() as directory:
        r=Path(directory);run=r/'run';run.mkdir(exist_ok=True);joint=r/'joint';joint.mkdir(exist_ok=True);period='fa18_inb';rng=np.random.default_rng(12)
        def dump(p,v):p.write_text(json.dumps(v))
        def arr(groups,n,ctrl):
         a=np.zeros((len(groups)*n,29 if ctrl else 12));a[:,0]=np.repeat(groups,n);a[:,1:8]=rng.normal(0,.3,(len(a),7));a[:,8]=1;a[:,9]=rng.uniform(.08,1.1,len(a));a[:,11]=np.arange(len(a))%2
         if ctrl:
          a[:,20]=rng.uniform(2,6,len(a));a[:,21]=rng.uniform(.1,.4,len(a));a[:,23]=1;a[:,24]=rng.uniform(.4,4,len(a));a[:,25]=rng.uniform(.1,.4,len(a));a[:,27]=1;a[:,28]=rng.normal(.135,.006,len(a))
         return a
        initial={s:arr(c.study.CONTROL_GROUPS if s.startswith('pi0') else range(3),1200 if s.startswith('pi0') else 6000 if s=='data' else 2400,s.startswith('pi0')) for s in ('data','pi0_data','mc','background','pi0_mc')}
        for stage in ('initial','candidate'):
         records=[]
         for s,a in initial.items():
          if stage=='candidate' and s in ('data','pi0_data'):continue
          a=a.copy()
          if stage=='candidate':a[:,9]*=.95;a[:,1:8]*=.9
          path=r/(stage+s+'.bin');a.astype('<f8').tofile(path);dump(Path(str(path)+'.json'),dict(shape=list(a.shape)));records.append(dict(period=period,sample=s,path=str(path),sha256=c.study.export_digest(path)))
         dump(run/(stage+'_exports.json'),dict(records=records,vectors={}))
        spec=dict(periods=[period],variables=[dict(name=n,label=n,bins=40,low=-1,high=1) for n in c.base.OBS],mass_variable=dict(name='m_gamma_gamma',label='Mass',bins=40,low=.1,high=.17))
        dump(run/'source_spec.json',spec);w=dict(data={},mc={});windows={period:dict(dvcs=[w]*3,control=[w]*4)}
        for n in ('initial_fixed_windows','final_windows'):dump(run/(n+'.json'),windows)
        dump(joint/'cache_manifest.json',dict(samples={period+'__'+s:dict(probability=1.) for s in ('mc','background','pi0_mc')}));dump(run/'signature.json',dict(cache_manifest=c.base.digest(joint/'cache_manifest.json'),helpers={},pipeline=c.base.digest(c.base.SOURCE/'study_selected_particle_response.py')));f=run/'response'/period;f.mkdir(parents=True,exist_ok=True);dump(f/'candidate_parameters.json',dict(background_multiplier=1.))
        c.main(['--response-run',str(run),'--joint-run',str(joint),'--output',str(r/'output'),'--min-entries','10'])
        
        assert (r/'output/REPORT.html').is_file()
        import csv
        rr=list(csv.DictReader((r/'output/final_window_scores.csv').open()))
        assert len(rr)>50 and all(np.isfinite(float(x['score'])) for x in rr)
        assert len(list((r/'output/plots/fa18_inb').glob('*.png')))==11
        print('PASS: full cache-only report, final-window scoring, mass quartiles and migration rendering')

if __name__=='__main__':
    main()
    integration()
