#!/usr/bin/env python3
"""Dense display sampling must preserve exact row-reference pulls and caches."""
import sys,copy
from pathlib import Path
import numpy as np
SOURCE=Path(__file__).resolve().parent
sys.path.insert(0,str(SOURCE/'external_scripts'))
try:import plot_pass2_final_results_by_t as plot
except ModuleNotFoundError:
 sys.path.insert(0,str(SOURCE));import plot_pass2_final_results_by_t as plot

def main():
    base=plot.base
    k=base.BinKey(.1,.2,1.,2.,.11,.15)
    cfg=base.ModelConfig(10.6,'unused','unused',sys.executable,24,False,False,False)
    refs={'xs':((.12,1.1,.12,7.,10.6),(.18,1.8,.14,43.,10.6)),
          'bsa':((.12,1.1,.12,7.,10.6),),'__display__':(.15,1.5,.13)}
    calls=[];saved_many,saved_bh=base.km15_xs_many,base.bh_xs
    def many(points,cfg):
        calls.extend(points);return [1+x+q+t+phi/100+h/10 for x,q,t,phi,h in points]
    base.km15_xs_many=many;base.bh_xs=lambda x,q,t,phi,cfg:1+x+q+t+phi/100
    try:
        _,curves=plot.compute_cell_models((k,cfg,('xs','bsa'),refs))
        assert len(curves['xs']['display']['phi'])==24
        assert curves['xs']['display']['phi'][0]==0 and curves['xs']['display']['phi'][-1]==360
        assert all(tuple(p[:3])==refs['__display__'] for p in calls if p[3] not in (7.,43.))
        assert curves['xs']['phi']==[7.,43.]
        exact=plot.interpolate_model(curves['xs'],'km15',43.)
        assert abs(exact-(1+.18+1.8+.14+.43))<1e-12
        assert np.isnan(plot.interpolate_model(curves['xs'],'km15',20.))
        assert len(plot.curve_coordinates(curves['xs'],'km15')[0])==24
        # Add display curves to an old cache entry without recomputing its data references.
        old={kind:{a:b for a,b in c.items() if not a.startswith('display')} for kind,c in curves.items()}
        refs['__cached__']=copy.deepcopy(old);calls.clear();cfg.phi_dense=12
        _,again=plot.compute_cell_models((k,cfg,('xs','bsa'),refs))
        assert len(calls)==12*3
        for kind in old:
            for field,value in old[kind].items():assert again[kind][field]==value
            assert len(again[kind]['display']['phi'])==12
        assert plot.interpolate_model(again['xs'],'km15',43.)==exact
        args=plot.make_parser().parse_args(['input.csv']);args._active_energy=10.6
        args._model_references={(10.6,k,'xs'):refs['xs']}
        key=plot.model_cache_key(k,args,'xs');args.phi_dense=48
        assert plot.model_cache_key(k,args,'xs')==key
    finally:base.km15_xs_many,base.bh_xs=saved_many,saved_bh
    print('PASS: requested dense sampling, fixed cell reference, exact pulls, deduplicated requests, old row-cache reuse, sampling changes preserve pull values')

if __name__=='__main__':main()
