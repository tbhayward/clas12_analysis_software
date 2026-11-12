#!/usr/bin/env python3
"""Production-algorithm and relative-efficiency checks without ROOT I/O."""
import argparse,tempfile,json
from pathlib import Path
import numpy as np
import study_particle_response as s
import study_response_cut_redetermination as c

def basic():
    normal=s.parameter_grid(s.TOPS);wide=s.parameter_grid(s.TOPS,'fit','wide')
    assert all(set(v)<=set(wide[k]) for k,v in normal.items())
    assert wide[3][0]==-.10 and wide[3][-1]==.10
    edges=np.r_[np.linspace(-1,-.12,45),np.linspace(-.12,.12,601)[1:],np.linspace(.12,1,45)[1:]]
    x=(edges[:-1]+edges[1:])/2;counts=np.exp(-.5*(x/.001)**2);counts[0]=.001
    limits=s.mass_plot_limits(edges,[np.r_[0,counts,0]])
    assert -.02<limits[0]<0<limits[1]<.02
    a=np.array([1,1,0,0],bool);b=np.array([1,0,0,0],bool)
    r,e=c.ratio_error(a,b);assert r==.5 and np.isclose(e,np.sqrt(.5)/2)
    assert c.ratio_error(a,a)==(1.,0.)
    rows=np.zeros((4,12));rows[:,8]=1;rows[:,1]=[0,.5,1,2]
    block=dict(mc=dict(Delta_phi=dict(mode='upper_quantile',cut_low=.3,cut_high=1)))
    assert c.window_mask(rows,block,'mc').tolist()==[True,True,True,False]
    print('PASS: wider bounds preserve standard grid; robust mass zoom; paired efficiency errors; production upper-cut semantics')

def synthetic(tmp):
    rng=np.random.default_rng(2026)
    def make(n,bg=False):
        a=np.zeros((n,12));a[:,0]=2;a[:,8]=1;a[:,9]=rng.uniform(.11,1,n);a[:,10]=1;a[:,11]=rng.integers(0,2,n)
        a[:,1]=rng.normal(np.pi,.07 if bg else .015,n);a[:,2]=rng.normal(2.7,.10 if bg else .07,n)
        a[:,3]=rng.exponential(.8 if bg else .2,n);a[:,4]=rng.exponential(.10 if bg else .025,n)
        a[:,5]=rng.normal(.6 if bg else .1,.25 if bg else .08,n)
        a[:,6]=rng.normal(.008 if bg else 0,.008 if bg else .002,n)
        a[:,7]=rng.normal(1.6 if bg else .88,.3 if bg else .1,n)
        return a
    arrays=dict(data=np.r_[make(10000),make(2000,True)],mc=make(18000),background=make(18000,True),pi0_data=make(12000),pi0_mc=make(18000))
    p=next(x for x in c.production.PERIODS if x.key=='fa18_inb')
    blocks=c.refit(p,2,arrays,tmp/'cuts');key=c.production.cut_json_key('dvcs',p,c.production.TOPOLOGIES[2])
    assert key in blocks and blocks[key]['mc'] and blocks[key]['data']
    assert (tmp/'cuts/dvcs_cut_development.png').is_file() and (tmp/'cuts/pi0_cut_development.png').is_file()
    # Same response and same refitted cuts must give exactly unit relative acceptance.
    allblocks={}
    for t in c.production.TOPOLOGIES:
        for channel in ['dvcs','eppi0']:allblocks[c.production.cut_json_key(channel,p,t)]=blocks[c.production.cut_json_key(channel,p,c.production.TOPOLOGIES[2])]
    rows=c.comparison_rows('fa18_inb','nominal',arrays,arrays,allblocks,allblocks,dict(mc=1.))
    usable=[r for r in rows if r['baseline_mc_selected']]
    assert usable and all(r['relative_acceptance']==1 and r['paired_mc_stat_error']==0 for r in usable)
    plot=tmp/'plots';plot.mkdir();c.plot_controls('fa18_inb','nominal',arrays,arrays,allblocks,plot)
    assert len(list(plot.glob('*.png')))==1
    print('PASS: actual production direct-pi0 calibration, separate DATA/MC iterative cuts, control plots and exact nominal acceptance replay')

def main():
    basic()
    with tempfile.TemporaryDirectory() as d:synthetic(Path(d))
if __name__=='__main__':main()
