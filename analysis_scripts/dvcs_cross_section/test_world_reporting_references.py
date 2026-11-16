#!/usr/bin/env python3
"""Regression checks for world-comparison reporting points and model cache identity."""
import csv
import importlib.util
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = Path(os.environ.get('DVCS_TEST_REPO', str(ROOT)))
sys.path.insert(0, str(REPO/'external_scripts'))
source = Path(os.environ.get('DVCS_WORLD_COMPARE_SCRIPT', str(REPO/'external_scripts/compare_world_dvcs_cross_sections.py')))
spec = importlib.util.spec_from_file_location('world_reference_test_backend', source)
world = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = world
spec.loader.exec_module(world)


def main():
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary); (root/'imports').mkdir()
        donor = {'bin index':11, 'xBmin':.204, 'xBmax':.268, 'Q2min':1.912,
                 'Q2max':2.510, 't_abs_min':.25, 't_abs_max':.4,
                 'phimin':0, 'phimax':15, 'xBavg':.23, 'Q2avg':2.2,
                 't_abs_avg':.32, 'phiavg':6, 'Fbin':1.05}
        donor_path=root/'imports/all_bin_v3.csv'
        with donor_path.open('w', newline='') as f:
            writer=csv.DictWriter(f, fieldnames=donor); writer.writeheader(); writer.writerow(donor)
        row={k:v for k,v in donor.items() if k not in ['xBavg','Q2avg','t_abs_avg','phiavg','Fbin']}
        row.update({'Fbin, 10.6 GeV':1.05, 'xBavg, 10.6 GeV':.24,
                    'Q2avg, 10.6 GeV':2.3, 't_abs_avg, 10.6 GeV':.33,
                    'phiavg, 10.6 GeV':8,
                    world.PASS2_XS_COL:'(10,1,0)'})
        raw=pd.DataFrame([row]); original=raw.copy(deep=True)
        result=world.apply_pass2_reporting_coordinates(raw, root/'data.csv')
        assert np.allclose(result[['xBavg, 10.6 GeV','Q2avg, 10.6 GeV','t_abs_avg, 10.6 GeV','phiavg, 10.6 GeV']].iloc[0], [.23,2.2,.32,6])
        assert result['measured_xBavg, 10.6 GeV'].iloc[0]==.24
        pd.testing.assert_frame_equal(raw, original)
        for badcolumn,badvalue in [('Fbin, 10.6 GeV',1.2),('phimax',30),('bin index',12)]:
            invalid=raw.copy(); invalid.loc[0,badcolumn]=badvalue
            try:world.apply_pass2_reporting_coordinates(invalid, root/'data.csv')
            except ValueError:pass
            else:raise AssertionError('Unvalidated imported donor: '+badcolumn)
        print('PASS: production reporting points, measured-coordinate preservation and donor/Fbin/geometry guards')

        calls=[]
        def evaluate(_backend, point, energy):
            calls.append((point.point_id,energy))
            return {'km15_ep':float(point.xB)+energy,'km15_bh':energy}
        world.evaluate_one_km15=evaluate
        points=pd.DataFrame([dict(point_id='pass2:11',xB=.23,Q2=2.2,t_abs=.32,
              phi_deg=6,ebeam=10.6,xs=10,stat_abs=1,ptp_sys_abs=1,point_unc_abs=2)])
        cache=root/'cache.csv'
        first=world.evaluate_km15_world(points,None,cache,10.6)
        assert len(calls)==1, 'Identical native and target energies should share one evaluation'
        world.evaluate_km15_world(points,None,cache,10.6);assert len(calls)==1
        changed=points.copy();changed.loc[0,'Q2']=2.4
        world.evaluate_km15_world(changed,None,cache,10.6);assert len(calls)==2
        world.evaluate_km15_world(changed,None,cache,10.2);assert len(calls)==4
        assert first["xs"].iloc[0]==10
        print('PASS: unchanged-coordinate cache reuse, changed kinematics/target invalidation and identical-energy evaluation reuse')


if __name__=='__main__':main()
