"""Design, fitted-anchor, independent-start, and per-bin worker regression checks."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import validate_cosine_axis_scans as scan


def nominal(u=.12,l=.14):
    return {'values':{'u1':u,'u2':-.1,'lu1':.18,'ul1':-.12,'ul2':-.17,'ll0':.65,'ll1':l}}


def configuration(bins=(4,6)):
    fits={str(b):nominal(.12+.01*b,.14-.01*b) for b in bins}
    return {'bins':list(bins),'modes':['observed'],'grids':{b:scan.grid_for(n,[-.5,-.25,0,.25,.5],[-.5,-.25,0,.25,.5]) for b,n in fits.items()}}


class AxisScans(unittest.TestCase):
    def test_defaults_cover_top_three_in_each_xb(self):
        self.assertEqual(scan.parser().parse_args([]).bins,[4,5,6,10,11,12,16,17,18,22,23,24])

    def test_two_lines_not_full_grid_and_nominal_not_duplicated(self):
        c=configuration();tasks=scan.tasks_for(c)
        self.assertEqual(len(tasks),22)
        self.assertEqual(len(tasks),len(set(tasks)))
        for b in c['bins']:
            g=c['grids'][str(b)];z=[t for t in tasks if t[0]==b]
            self.assertEqual(sum(scan.scan_label(c,t)=='nominal' for t in z),1)
            for _,_,i,j in z:
                self.assertTrue(i==g['u1_anchor'] or j==g['ll1_anchor'])

    def test_exact_anchor_grid_overlap_is_deduplicated(self):
        c={'bins':[6],'modes':['observed'],'grids':{'6':scan.grid_for(nominal(0,.25),[-.5,-.25,0,.25,.5],[-.5,-.25,0,.25,.5])}}
        self.assertEqual(len(scan.tasks_for(c)),9)

    def test_measured_pairs_prioritized_across_bins(self):
        c=configuration()
        self.assertTrue(all(scan.scan_label(c,t)=='nominal' for t in scan.tasks_for(c)[:2]))

    def test_starts_independent_of_injection(self):
        n=nominal();a=scan.core.fit_start_values(n)
        for u,l in [(-.5,0),(.5,.5)]:
            truth=scan.core.make_truth(n,u,l)
            self.assertEqual(a,scan.core.fit_start_values(n))
            self.assertEqual((a['u1'],a['ll1']),(0.,0.))
            self.assertEqual(truth['u2'],n['values']['u2'])

    def test_worker_uses_correct_anchor_after_switching_bins(self):
        c=configuration();data=np.zeros(1)
        with patch.object(scan.core,'_WORK',(data,{}, {},{},c,Path('/tmp'), 'sig',{})):
            with patch.object(scan,'_SETTINGS',c):
                def capture(task):return scan.core._WORK[4]
                with patch.object(scan.core,'run_cell',side_effect=capture):
                    for task in [scan.tasks_for(c)[0],scan.tasks_for(c)[1],scan.tasks_for(c)[0]]:
                        cfg=scan.run_cell(task);g=c['grids'][str(task[0])]
                        self.assertEqual(cfg['u1_values'],g['u1_values'])
                        self.assertEqual(cfg['ll1_values'],g['ll1_values'])

    def test_products_work_when_only_nonphysical_cell_exists(self):
        c=configuration((6,));task=scan.tasks_for(c)[0];g=c['grids']['6']
        n=nominal();n['errors']={p:.2 for p in scan.prod.PHYSICS_PARAMETERS}
        state={'signature':'sig','task':list(task),'bin_number':6,'mode':'observed',
               'u1_index':task[2],'ll1_index':task[3],'requested':2,
               'truth':scan.core.make_truth(n,g['u1_values'][task[2]],g['ll1_values'][task[3]]),
               'status':'nonphysical_injection','complete':True}
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            scan.support.atomic_json(out/'manifest.json',{'signature':'sig','config':c,'baselines':{'6':n}})
            scan.support.atomic_json(scan.core.cell_path(out,task)/'state.json',state)
            scan.write_products(out)
            self.assertTrue((out/'tables/injection_summary.csv').exists())
            self.assertEqual(len(list((out/'plots').glob('*.png'))),4)

    def test_per_bin_anchor_truth_saved_and_resume_uses_same_checkpoint(self):
        from test_validate_unpolarized_modulations import fixture
        e,states,records,truth,d=fixture(30)
        b=6;e['bin_number'][:]=b
        n=nominal();n['values'].update({f'f_{p}':v for p,v in d.items()})
        n['metadata']={'number_of_events':len(e['phi'])}
        c=configuration((b,));c.update(replicas=2,seed=42,profile_coverage=False)
        task=scan.tasks_for(c)[0]
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(scan.core,'_WORK',(e,states,records,{b:n},c,Path(tmp),'sig',{})),patch.object(scan,'_SETTINGS',c):
                # No expensive fit needed to test generator/checkpoint/resume behavior.
                fit={'values':n['values']}
                with patch.object(scan.prod,'fit_one_variant',return_value=fit) as fitting,patch.object(scan.support,'compact_fit',side_effect=lambda f:f),patch.object(scan.core,'evaluate_recovery',return_value={'usable':True,'profiles':{}}):
                    scan.run_cell(task);self.assertEqual(fitting.call_count,2)
                    result=scan.run_cell(task);self.assertTrue(result['reused']);self.assertEqual(fitting.call_count,2)
                import json
                state=json.loads((scan.core.cell_path(Path(tmp),task)/'state.json').read_text())
                g=c['grids'][str(b)]
                self.assertEqual(state['truth']['u1'],g['u1_values'][g['u1_anchor']])
                self.assertEqual(state['truth']['ll1'],g['ll1_values'][g['ll1_anchor']])

if __name__=='__main__':unittest.main()
