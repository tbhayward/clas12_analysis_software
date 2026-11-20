"""Tests for spin-blind matched coverage, grid truth, accounting and checkpoint locks."""
import fcntl
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import validate_cosine_coverage as study
from test_validate_unpolarized_modulations import fixture


class CoverageChecks(unittest.TestCase):
    def setUp(self):
        self.events,self.states,self.records,self.truth,self.dilutions=fixture(1000)
        self.config={'seed':20261120,'retained_fraction':.65,'gap_half_width':30.,
                     'u1_values':[-.5,0,.5],'ll1_values':[-.5,0,.5],
                     'replicas':2,'profile_coverage':True,'bins':[1],'modes':list(study.MODES)}
        self.nominal={'values':{**self.truth,**{f'f_{p}':v for p,v in self.dilutions.items()}},
                      'errors':{p:.2 for p in study.prod.PHYSICS_PARAMETERS},
                      'metadata':{'number_of_events':len(self.events['phi'])}}

    def test_observed_mode_keeps_exact_event_kinematics(self):
        indices=study.coverage_indices(self.events,'observed',self.config,1,0)
        np.testing.assert_array_equal(indices,np.arange(len(self.events['phi'])))

    def test_matched_thinning_counts_and_gap(self):
        uniform=study.coverage_indices(self.events,'matched_uniform',self.config,1,0)
        gap=study.coverage_indices(self.events,'matched_forward_gap',self.config,1,0)
        for period in range(3):
            expected=int(np.floor(self.config['retained_fraction']*sum(self.events['period_index']==period)))
            self.assertEqual(sum(self.events['period_index'][uniform]==period),expected)
            self.assertEqual(sum(self.events['period_index'][gap]==period),expected)
        phi=np.abs((self.events['phi'][gap]+np.pi)%(2*np.pi)-np.pi)
        self.assertTrue((phi>=np.deg2rad(30)).all())
        np.testing.assert_array_equal(gap,study.coverage_indices(self.events,'matched_forward_gap',self.config,1,0))

    def test_masks_ignore_spin_labels(self):
        changed={**self.events,'runnum':np.full_like(self.events['runnum'],-999),
                 'helicity':-self.events['helicity']}
        for mode in study.MODES:
            np.testing.assert_array_equal(study.coverage_indices(self.events,mode,self.config,1,3),
                                          study.coverage_indices(changed,mode,self.config,1,3))

    def test_injected_targets_do_not_change_fit_starts(self):
        starts=study.fit_start_values(self.nominal)
        for u1,ll1 in [(-.5,.5),(.3,-.2),(0,0)]:
            truth=study.make_truth(self.nominal,u1,ll1)
            self.assertEqual(truth['u1'],u1);self.assertEqual(truth['ll1'],ll1)
            self.assertEqual(study.fit_start_values(self.nominal),starts)
        self.assertEqual(starts['u1'],0);self.assertEqual(starts['ll1'],0)

    def test_summary_bias_and_coverage_failure_accounting(self):
        state={'bin_number':1,'mode':'observed','u1_index':1,'ll1_index':1,
               'truth':self.truth,'status':'complete','requested':3,'selected_events':100}
        rows=[];names=list(study.prod.PHYSICS_PARAMETERS)
        for sign,q in [(1,.5),(-1,5.)]:
            values={p:self.truth[p]+sign*.05 for p in names}
            fit={'values':values,'errors':{p:.1 for p in names},'parameter_order':names,
                 'correlation':np.eye(len(names)).tolist()}
            rows.append({'usable':True,'fit':fit,'profiles':{p:{'status':'evaluated','two_delta_nll':q} for p in study.TARGETS}})
        rows.append({'status':'failed','usable':False,'profiles':{}})
        summary={r['parameter']:r for r in study.summarize_cell(state,rows,self.nominal)}
        row=summary['u1'];self.assertAlmostEqual(row['mean_bias'],0)
        self.assertAlmostEqual(row['bias_mc_error'],.05)
        self.assertEqual(row['usable'],2);self.assertEqual(row['failed_or_unusable'],1)
        self.assertEqual(row['hesse_coverage68_usable'],1)
        self.assertAlmostEqual(row['hesse_coverage68_all_attempts'],2/3)
        self.assertEqual(row['profile_coverage68_evaluated'],.5)
        self.assertAlmostEqual(row['profile_coverage68_all_attempts'],1/3)

    def test_nonphysical_cell_is_recorded_without_fitting(self):
        with tempfile.TemporaryDirectory() as directory:
            work=(self.events,self.states,self.records,{1:self.nominal},self.config,Path(directory),'signature',{})
            task=(1,'observed',0,0)
            with (patch.object(study,'_WORK',work),
                  patch.object(study.support,'generator_envelopes',side_effect=ValueError('nonpositive')),
                  patch.object(study.prod,'fit_one_variant') as fitter):
                result=study.run_cell(task)
                self.assertEqual(result['status'],'nonphysical_injection');fitter.assert_not_called()
            state=json.loads((study.cell_path(Path(directory),task)/'state.json').read_text())
            self.assertTrue(state['complete']);self.assertEqual(state['truth']['u1'],-.5)

    def test_cell_lock_blocks_duplicate_worker(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory);task=(1,'observed',1,1);folder=study.cell_path(output,task);folder.mkdir(parents=True)
            with (folder/'owner.lock').open('a+') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                with patch.object(study,'_WORK',(None,None,None,None,None,output,None,None)):
                    self.assertEqual(study.run_cell(task)['status'],'owned_elsewhere')


if __name__=='__main__':unittest.main()
