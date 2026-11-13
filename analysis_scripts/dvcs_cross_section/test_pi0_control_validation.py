#!/usr/bin/env python3
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import study_pi0_control_validation as study


def make(n,t,q):
    a=np.zeros((n,20));a[:,9]=t;a[:,13]=q;a[:,12]=.25;a[:,14:]=.5;a[:,8]=1;return a


def main():
    d=np.r_[make(120,.2,2),make(240,.3,3)];m=np.r_[make(240,.2,2),make(120,.3,3)]
    alpha,w,rows=study.fit_weights(d,m,50);assert alpha==1 and np.isclose(w[0,0],.5) and np.isclose(w[1,1],2)
    assert np.isclose(study.weights_for(m,w).sum(),len(d))
    assert study.weights_for(make(1,1.2,15),w)[0]==1
    _,_,fallback=study.fit_weights(d,m,10000);assert all(r['status']=='unity fallback' for r in fallback)
    test=np.r_[make(240,.2,2),make(480,.3,3)];test[:,11]=1
    w1,_,sc1=study.train_models(np.r_[d,test],np.r_[m,test],np.r_[m,test],50)
    altered=test.copy();altered[:,9]=.9;altered[:,13]=8
    w2,_,sc2=study.train_models(np.r_[d,altered,altered],np.r_[m,altered],np.r_[m,altered,altered],50)
    assert np.array_equal(w1,w2) and sc1==sc2
    with tempfile.TemporaryDirectory() as folder:
        out=Path(folder)
        with patch.dict(study.control.EDGES,{'t_abs':(9,[.11,.25,.4,.6,1])},clear=True):
            r,totals,links=study.evaluate(d,dict(nominal=m,add_pt=m,add_pt_weighted=m),dict(nominal=1,add_pt=1,add_pt_weighted=1),w,'fa18_inb',out)
        assert all((out/l).is_file() for l in links) and len(links)==4
        good=[x for x in r if x['scenario']=='add_pt_weighted' and x['test_prediction']>0]
        assert all(np.isclose(x['data_over_prediction'],1) for x in good)
        assert all(np.isclose(x['data_over_prediction'],1) for x in totals)
    try:study.fit_weights(d,m[:0])
    except ValueError:pass
    else:raise AssertionError('Accepted empty training MC')
    print('PASS: known coarse weights, sparse-cell fallback, out-of-grid unity, unchanged training weights under test evaluation, weighted held-out predictions and report plots')

if __name__=='__main__':main()
