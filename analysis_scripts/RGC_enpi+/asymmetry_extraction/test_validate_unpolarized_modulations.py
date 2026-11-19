"""Synthetic checks against the actual production conditional likelihood."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import validate_unpolarized_modulations as study
prod = study.prod


def fixture(n_per_period=2000, seed=71):
    rng = np.random.default_rng(seed)
    states, records = {}, {}
    for i,p in enumerate(prod.PERIODS):
        states[p] = {'run': np.array([16000+i*100,16001+i*100], dtype=np.int32),
                     'pt': np.array([.72,-.67]), 'q_plus': np.array([1.7,.8]),
                     'q_minus': np.array([.9,1.3])}
        for b in range(1,25):
            records[(p,b)] = prod.DilutionRecord(p,b,(b-1)//6,(b-1)%6,.30,.012)
    n = n_per_period*3
    phi = rng.uniform(0,2*np.pi,n)
    # A nonuniform acceptance, with event-dependent depolarization factors.
    keep = rng.random(n) < (.55+.4*np.cos(phi))
    phi = phi[keep]
    pi = np.repeat(np.arange(3), n_per_period)[keep].astype(np.int8)
    eps = rng.uniform(.45,.85,len(phi))
    a = np.ones(len(phi))
    e = {'period_index':pi, 'runnum':np.array([states[prod.PERIODS[i]]['run'][0] for i in pi]),
         'helicity':np.ones(len(phi),dtype=np.int8), 'bin_number':np.ones(len(phi),dtype=np.int16),
         'phi':phi, 'rB':eps, 'rC':np.sqrt(1-eps**2), 'rV':np.sqrt(2*eps*(1+eps)),
         'rW':np.sqrt(2*eps*(1-eps)), 'sin_theta_gamma':a*.2,
         'cos_theta_gamma':a*np.sqrt(.96), 'xB':a*.2, 'Q2':a*2.1, 'W':a*3.,
         'minus_t':a*.2, 'minus_tprime':a*.12, 'Mx2':a*.88, 'epsilon':eps,
         'DepA':a, 'DepB':eps, 'DepC':np.sqrt(1-eps**2),
         'DepV':np.sqrt(2*eps*(1+eps)), 'DepW':np.sqrt(2*eps*(1-eps))}
    truth = dict(u1=.12,u2=-.10,lu1=.18,ul1=-.12,ul2=-.17,ll0=.65,ll1=.14)
    d = {p:.30 for p in prod.PERIODS}
    return study.generate_states(e,states,truth,d,rng),states,records,truth,d


def write_cli_fixture(folder):
    folder = Path(folder).resolve(); folder.mkdir(parents=True,exist_ok=True)
    e,states,_,_,_ = fixture(6000)
    run_file = folder/'runs.csv'
    lines=[]
    for p in prod.PERIODS:
        lines.append(f'# RGC {prod.PERIOD_LABELS[p]} NH3')
        s=states[p]
        for i,r in enumerate(s['run']):
            lines.append(f'{r},0,{s["q_plus"][i]},{s["q_minus"][i]},{s["pt"][i]},0.01')
    run_file.write_text('\n'.join(lines)+'\n')
    cuts={'periods':{p:[{'bin_number':b,'x_index':(b-1)//6,'t_index':(b-1)%6,'nominal':[.7,1.]} for b in range(1,25)] for p in prod.PERIODS}}
    cut_path=folder/'cuts.json';cut_path.write_text(json.dumps(cuts))
    dil={'source_exclusivity_json_sha256':study.digest(cut_path),'periods':{p:{'bins':[{'bin_number':b,'x_index':(b-1)//6,'t_index':(b-1)%6,'cuts':{'nominal':{'recommended':{'value':.3,'stat_uncertainty':.012}}}} for b in range(1,25)]} for p in prod.PERIODS}}
    (folder/'dilution.json').write_text(json.dumps(dil))
    inputs={}
    for p in prod.PERIODS:
        inputs[p]=folder/f'{p}.root';inputs[p].write_text('Synthetic stat-only provenance fixture; not read as ROOT')
    runs=prod.parse_run_info_csv(run_file)
    cache=folder/'selected_events.npz';actual=prod._event_cache_path(cache);actual.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(actual,**e)
    prod._write_cache_metadata(actual,prod.load_channel_cuts(cut_path),prod._cache_identity(inputs,'PhysicsEvents',runs))
    return folder


class ConditionalValidationTests(unittest.TestCase):
    def test_state_frequencies_match_all_run_charge_weighted_probabilities(self):
        e,states,records,truth,d=fixture(100000)
        for p in prod.PERIODS:
            idx=np.flatnonzero(e['period_index']==prod.PERIOD_INDEX[p]);s=states[p]
            weights=[]
            for r in range(2):
                for h,q in ((1,s['q_plus'][r]),(-1,s['q_minus'][r])):
                    weights.append(q*study.cross_factor(e,s,p,truth,d[p],indices=idx,h=h,pt=s['pt'][r]))
            weights=np.array(weights);prob=weights/weights.sum(axis=0)
            for j in range(4):
                observed=np.sum((e['runnum'][idx]==s['run'][j//2]) & (e['helicity'][idx]==(1 if j%2==0 else -1)))
                mean=prob[j].sum();sigma=np.sqrt(np.sum(prob[j]*(1-prob[j])))
                self.assertLess(abs(observed-mean),6*sigma)

    def test_generated_sample_nll_matches_explicit_run_state_sum(self):
        e,states,records,truth,d=fixture(200)
        nll,_=prod.make_bin_nll(e,states,records,1,'nominal')
        manual=0.
        for p in prod.PERIODS:
            idx=np.flatnonzero(e['period_index']==prod.PERIOD_INDEX[p]);s=states[p]
            denom=np.zeros(idx.size);num=np.zeros(idx.size)
            for r in range(2):
                for h,q in ((1,s['q_plus'][r]),(-1,s['q_minus'][r])):
                    w=q*study.cross_factor(e,s,p,truth,d[p],indices=idx,h=h,pt=s['pt'][r])
                    denom+=w
                    observed=(e['runnum'][idx]==s['run'][r]) & (e['helicity'][idx]==h)
                    num[observed]=w[observed]
            manual-=np.log(num/denom).sum()
        self.assertAlmostEqual(nll(**truth,**{f'f_{p}':d[p] for p in d}),manual,places=8)

    def test_zero_polarization_is_exactly_uninformative_about_uu(self):
        e,states,records,truth,d=fixture(300)
        null=dict(truth);null.update({p:0. for p in prod.POLARIZED_PARAMETERS})
        a=study.generate_states(e,states,null,d,np.random.default_rng(123))
        null.update(u1=-.1,u2=.1)
        b=study.generate_states(e,states,null,d,np.random.default_rng(123))
        np.testing.assert_array_equal(a['runnum'],b['runnum'])
        np.testing.assert_array_equal(a['helicity'],b['helicity'])
        self.assertTrue(study.null_identifiability(a,states,records,1,null,d)['flat_within_numerical_precision'])

    def test_nonphysical_scenario_is_rejected_not_clipped(self):
        e,states,records,truth,d=fixture(100)
        truth['u1']=3.
        with self.assertRaises(ValueError):study.generator_envelopes(e,states,truth,d)

    def test_auxiliary_draws_have_correct_mean_and_width(self):
        _,_,records,_,d=fixture(10);rng=np.random.default_rng(17)
        draws=np.array([study.auxiliary_records(records,1,d,rng,'fluctuate')[1]['su22'] for _ in range(5000)])
        self.assertLess(abs(draws.mean()-.3),.0008)
        self.assertLess(abs(draws.std()-.012),.0005)
        self.assertEqual(study.auxiliary_records(records,1,d,rng,'fixed_truth')[1],d)

    def test_crossings_do_not_bridge_a_failed_scan_point(self):
        rows=[{'value':0.,'two_delta_nll':0.,'valid':True},
              {'value':1.,'two_delta_nll':None,'valid':False},
              {'value':2.,'two_delta_nll':4.,'valid':True}]
        self.assertIsNone(study.profile_crossing(rows,0.,1,1.))

    def test_summary_exposes_failed_fits_and_undefined_null_coverage(self):
        fit={'values':{'u1':.1,'u2':.1},'errors':{'u1':.1,'u2':.1},'parameters_at_limit':False}
        scenario={'truth':{'u1':0.,'u2':0.},'status':'complete','rows':[
            {'replica':0,'usable':True,'status':'fit','fit':fit},
            {'replica':1,'usable':False,'status':'failed'}]}
        config={'replicas':2,'auxiliary_mode':'fluctuate'}
        rows=study.summaries([{'bin_number':1,'scenarios':{'nominal':scenario,'zero_polarized':scenario}}],config)[1]
        self.assertEqual(rows[0]['coverage68_usable'],1.)
        self.assertEqual(rows[0]['coverage68_all_attempts'],.5)
        self.assertIsNone(rows[2]['coverage68_usable'])

    def test_production_fitter_recovers_nonzero_uu_on_nonuniform_acceptance(self):
        e,states,records,truth,d=fixture(12000,seed=97)
        fit=prod.fit_one_variant(e,states,records,1,'nominal')
        self.assertTrue(study.usable(fit),fit)
        for p in study.UU:
            self.assertLess(abs(fit['values'][p]-truth[p]),5*fit['errors'][p])


class DistributedChecks(unittest.TestCase):
    def test_bin_lock_excludes_second_worker_and_releases(self):
        import fcntl
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output/'locks').mkdir()
            with patch.object(study, '_WORK', (None,None,None,None,output)):
                with (output/'locks/bin_01.lock').open('a+') as held:
                    fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    with patch.object(study, 'run_bin_locked') as work:
                        self.assertEqual(study.run_bin(1)['status'], 'owned_elsewhere')
                        work.assert_not_called()
                with patch.object(study, 'run_bin_locked', return_value={'status':'complete'}) as work:
                    self.assertEqual(study.run_bin(1)['status'], 'complete')
                    work.assert_called_once_with(1)

    def test_event_identity_still_rejects_changed_inputs(self):
        old = {'selection_code_sha256':'old', 'inputs':{'file':'original'}}
        self.assertTrue(study.same_event_identity(old, dict(old, selection_code_sha256='new')))
        self.assertFalse(study.same_event_identity(old, dict(old, inputs={'file':'changed'})))

    def test_atomic_checkpoint_replacement(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'bin_01.json'
            study.atomic_json(path, {'rows':[1]})
            study.atomic_json(path, {'rows':[1,2]})
            self.assertEqual(json.loads(path.read_text()), {'rows':[1,2]})
            self.assertFalse(list(Path(directory).glob('*.tmp')))


if __name__=='__main__':unittest.main()
