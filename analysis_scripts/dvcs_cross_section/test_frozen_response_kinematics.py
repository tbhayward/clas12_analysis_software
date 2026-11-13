#!/usr/bin/env python3
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_frozen_response_kinematics as study


def make(n=400):
    a=np.zeros((n,28));a[:,0]=2;a[:,8]=1;a[:,11]=np.arange(n)%2;a[:,9]=.2;a[:,1:8]=.001;a[:,4]=np.linspace(.001,.2,n);a[:,14]=np.linspace(3,7,n);a[:,15]=np.linspace(.1,.3,n);a[:,16]=.6;a[:,17]=.8;a[:,20]=np.linspace(2,6,n);a[:,21]=np.linspace(.1,.4,n);a[:,23]=1;a[:,24]=np.linspace(.4,2,n);a[:,25]=np.linspace(.2,.5,n);a[:,27]=0;return a


def main():
    a=make();splits=study.choose_splits(a,a)
    split=next(s for s in splits if s['category']=='electron_FD' and s['axis']=='momentum');assert split['status']=='usable'
    low=study.masks(a,split,'low');high=study.masks(a,split,'high');assert not (low&high).any() and (low|high).all()
    assert next(s['threshold'] for s in splits if s['category']=='photon1_FD' and s['axis']=='momentum')==next(s['threshold'] for s in splits if s['category']=='photon2_FD' and s['axis']=='momentum')
    swapped=a.copy();swapped[:,24]=swapped[:,20];swapped[:,25]=swapped[:,21];swapped[:,27]=swapped[:,23]
    second=next(s for s in splits if s['category']=='photon2_FD' and s['axis']=='momentum');first=next(s for s in splits if s['category']=='photon1_FD' and s['axis']=='momentum')
    assert np.array_equal(study.masks(swapped,first,'low'),study.masks(swapped,second,'low'))
    assert study.category_mask(a,'photon1_FD').all() and study.category_mask(a,'photon2_FT').all() and not study.category_mask(a,'photon2_FD').any()
    changed=a.copy();changed[changed[:,11]==1,14]=1000
    assert next(s['threshold'] for s in study.choose_splits(changed,changed) if s['category']=='electron_FD' and s['axis']=='momentum')==split['threshold']
    missing=a.copy();missing[:,24]=np.nan;assert not study.category_mask(missing,'photon2_FT').any() and study.category_mask(missing,'electron_FD').all()
    cfg=next(p for p in study.cuts.production.PERIODS if p.key=='fa18_inb');blocks={study.cuts.production.cut_json_key(c,cfg,t):dict(data={},mc={}) for c in ('dvcs','eppi0') for t in study.cuts.production.TOPOLOGIES}
    variables=[dict(name=n,label=n,bins=40,low=-.1,high=.3) for n in study.base.OBS];spec=dict(exclusivity=blocks,variables=variables)
    cases=dict(data=np.r_[a,a],pi0_data=a,nominal=dict(mc=a,background=a,pi0_mc=a),candidate=dict(mc=a,background=a,pi0_mc=a));scales={(s,h,t):(.5,.5) for s in ('nominal','candidate') for h in (0,1) for t in range(3)};prob=dict(mc=1,background=1)
    edges=np.linspace(-.1,.3,41);block=next(iter(blocks.values()));direct=study.histogram(a,2,block,'pTmiss',split,'low',edges);cached=study.histogram(a,2,block,'pTmiss',split,'low',edges,selected_base=study.masks(a,split,'low'));assert np.array_equal(direct,cached)
    with tempfile.TemporaryDirectory() as folder:
        out=Path(folder);rows,links=study.analyze('fa18_inb',cases,spec,scales,prob,[split],out,10)
        assert len(links)==1 and all((out/r[k]).is_file() for r in links for k in ('scores','residuals','examples'))
        assert {r['side'] for r in rows}=={'low','high'} and {r['channel'] for r in rows}=={'epgamma','epgammagamma'}
    print('PASS: fixed half-0 DATA medians, disjoint matched low/high masks, independent photon roles/tags, unresolved coverage, cached selection equivalence and two-channel reports')

if __name__=='__main__':main()
