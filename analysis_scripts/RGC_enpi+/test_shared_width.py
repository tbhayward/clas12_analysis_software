"""Regression checks for common neutron widths and selection provenance.

Run from RGC_enpi+: python3 -m unittest -v test_shared_width
"""
import copy
import importlib.util
import json
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent

def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, BASE / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

c = load('rgc_width_channel', 'channel_selection/channel_selection_mx2_fits.py')
a = load('rgc_width_asymmetry', 'asymmetry_extraction/extract_structure_function_ratios.py')
d = load('rgc_width_dilution', 'dilution_factor/determine_dilution_factor.py')
PERIODS = ('su22', 'fa22', 'sp23')

class SharedWidthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(472)
        cls.jobs, cls.carbons, cls.scales = [], [], {}
        for period, signal, background in zip(PERIODS, (12000, 30000, 6500), (9000, 14000, 4000)):
            # Unequal period statistics and backgrounds, common physical peak.
            values = np.r_[rng.normal(0.883, 0.060, signal), rng.uniform(0, 1.55, background)]
            job = c.FitJob(period, 'after', 0, 0, .10, .25, .05, .25,
                           values, 0., 1.55, 120, .50, 1.35, 20, .95)
            cls.jobs.append(job)
            cls.carbons.append(replace(job, values=rng.uniform(0, 1.55, background * 2)))
            cls.scales[(period, 'after')] = {'carbon_scale': .5, 'carbon_scale_uncertainty': .015}
        cls.config = c.JointFitConfig(10., .1, 2., 3., .72, 1.04, (1.30,), False, 0, 7)

    def test_polynomial_joint_recovery_all_orders_and_fixed_mean(self):
        for model in ('quadratic', 'cubic', 'quartic'):
            for fixed in (None, .883):
                with self.subTest(model=model, fixed=fixed):
                    fit = c.fit_joint_background_model(model, self.jobs, self.config, fixed_mean_gev2=fixed)
                    self.assertTrue(fit['success'])
                    self.assertAlmostEqual(fit['shared_sigma_gev2'], .060, delta=.003)
                    self.assertAlmostEqual(fit['shared_mean_gev2'], .883, delta=.003)
                    ncoef = c.polynomial_parameter_count(model)
                    self.assertEqual(fit['number_of_parameters'], 1 + int(fixed is None) + 3*(1+ncoef))
                    self.assertEqual(fit['joint_ndf'], sum(np.count_nonzero(c.histogram_for_job(j)['fit_mask']) for j in self.jobs)-fit['number_of_parameters'])
                    widths = [r['sigma_gev2'] for r in fit['period_results'].values()]
                    errors = [r['sigma_error_gev2'] for r in fit['period_results'].values()]
                    self.assertEqual(len(set(widths)), 1)
                    self.assertEqual(len(set(errors)), 1)
                    self.assertTrue(0 < errors[0] < .003)

    def test_carbon_joint_recovery_and_covariance_indices(self):
        for order in (1, 2):
            fit = c.fit_one_carbon_candidate(self.jobs, self.carbons, self.scales,
                                            'after', 0, 0, .68, 1.10, order)
            self.assertTrue(fit['success'])
            self.assertAlmostEqual(fit['shared_sigma_gev2'], .060, delta=.004)
            self.assertEqual(fit['number_of_parameters'], 2 + 3*(2+order))
            self.assertEqual(len(fit['parameter_names']), fit['number_of_parameters'])
            error = np.sqrt(fit['covariance'][1][1])
            for i, record in enumerate(fit['period_records']):
                self.assertEqual(record['sigma_gev2'], fit['parameters'][1])
                self.assertEqual(record['sigma_error_gev2'], error)
                self.assertEqual(len(record['residual_coefficient_errors']), order)
                for k, value in enumerate(record['residual_coefficient_errors']):
                    idx = 2 + i*(2+order) + 2 + k
                    self.assertAlmostEqual(value, np.sqrt(fit['covariance'][idx][idx]))

    def cut_table(self):
        records = []
        for xi in range(4):
            for ti in range(6):
                for period in PERIODS:
                    records.append(dict(stage='after', is_recommended=True,
                        bin_number=c.combined_bin_number(xi, ti), x_index=xi, t_index=ti,
                        period=period, shared_mean_gev2=.883, shared_mean_error_gev2=.001,
                        sigma_gev2=.06+.001*ti, sigma_error_gev2=.002))
        table = c.build_final_polynomial_cut_table(pd.DataFrame(records))
        table["classification"] = "resolved"
        c.validate_shared_width_table(table)
        return table

    def test_all_bins_json_python_cpp_latex_and_consumer_roundtrip(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); table = self.cut_table()
            cut = root/'cuts.json'
            c.write_cut_json_from_table(table, cut, 'synthetic shared fit')
            lookup = root/'cuts.py'; header = root/'cuts.h'; latex = root/'cuts.tex'
            c.write_python_cut_lookup(table, lookup)
            c.write_cpp_cut_lookup(table, header)
            c.write_final_cut_latex_table(table, latex)
            namespace = {}; exec(lookup.read_text(), namespace)
            self.assertIn('shared', lookup.read_text())
            self.assertIn('kMx2CutsGeV2', header.read_text())
            self.assertIn('\n', latex.read_text())
            self.assertNotIn('sigma_{\\rm Su22}', latex.read_text())
            for label in ('tight', 'nominal', 'loose'):
                loaded = a.load_channel_cuts(cut, label)
                dilution = d.load_exclusivity_cuts(cut)
                self.assertEqual(len(loaded), 72)
                for row in table.itertuples():
                    for period in PERIODS:
                        entry = loaded[(period, row.bin_number)]
                        interval = namespace['get_mx2_cut'](period, row.x_index, row.t_index, label)
                        np.testing.assert_allclose(interval, [entry.low_gev2, entry.high_gev2])
                        np.testing.assert_allclose(interval, dilution[(period,row.x_index,row.t_index)].interval(label))
            # Exports reject accidental period-dependent regressions.
            bad = table.copy(); bad.loc[0, 'sigma_fa22_gev2'] += .01
            with self.assertRaises(RuntimeError): c.validate_shared_width_table(bad)
            payload = json.loads(cut.read_text()); payload['periods']['fa22'][0]['nominal'][0] += .01
            cut.write_text(json.dumps(payload))
            with self.assertRaises(RuntimeError): a.load_channel_cuts(cut)
            with self.assertRaises(RuntimeError): d.load_exclusivity_cuts(cut)

    def test_final_plots_and_carbon_products(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); table = self.cut_table()
            c.plot_master_cut_summary(table, root/'summary.png')
            c.plot_polynomial_master_cut_summary(table, root/'polynomial.png')
            c.plot_nominal_vs_isr_mu_sigma(table, table, root/'isr.png')
            c.plot_nominal_vs_isr_parameter_differences(table, table, root/'delta.png')
            self.assertEqual(len(list(root.glob('*.png'))),4)

    def test_carbon_production_exports_and_replica_path(self):
        fit = c.fit_one_carbon_candidate(self.jobs, self.carbons, self.scales,
                                        'after', 0, 0, .68, 1.10, 1)
        results = []
        for xi in range(4):
            for ti in range(6):
                item = copy.deepcopy(fit)
                item.update(x_index=xi, t_index=ti, selected_model=1, classification='resolved',
                    maximum_abs_fit_range_delta_mean_gev2=0., maximum_abs_control_delta_mean_gev2=0.,
                    model_mean_difference_gev2=0.)
                for record in item['period_records']:
                    record.update(x_index=xi, t_index=ti, bin_number=c.combined_bin_number(xi,ti))
                results.append(item)
        polynomial_rows = []
        for row in self.cut_table().itertuples():
            for period in PERIODS:
                polynomial_rows.append(dict(stage='after',is_recommended=True,period=period,
                    x_index=row.x_index,t_index=row.t_index,bin_number=row.bin_number,
                    shared_mean_gev2=.883,shared_mean_error_gev2=.001,
                    sigma_gev2=.060,sigma_error_gev2=.002))
        frame=pd.DataFrame(polynomial_rows)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            # Exercise all production writers; the large fit canvases have separate plot tests.
            with patch.object(c,'plot_final_cut_canvases',return_value=[]):
                products=c.write_final_carbon_cut_products(results,frame,root)
            payload=json.loads(Path(products['json']).read_text())
            self.assertEqual(payload['width_model'],'shared_per_bin')
            self.assertEqual(payload['cross_period_sigma_correlation'],1.)
            self.assertEqual(len(a.load_channel_cuts(Path(products['json']))),72)
            self.assertEqual(len(d.load_exclusivity_cuts(Path(products['json']))),72)
            meta=json.loads(Path(products['metadata']).read_text())
            self.assertFalse(meta['period_specific_widths'])
        joint=c.fit_joint_background_model('quadratic',self.jobs,self.config)
        models={name:joint if name == 'quadratic' else {'success':False,'period_results':{}} for name in c.BACKGROUND_MODELS}
        records=c.assemble_joint_period_results(self.jobs,models,'quadratic','test',self.config)
        c.run_joint_replicas(records,[self.jobs],replace(self.config,replicas=3))
        summaries=[item['replica_summary'] for item in records]
        self.assertEqual(summaries[0]['sigma'],summaries[1]['sigma'])
        self.assertEqual(summaries[0]['sigma'],summaries[2]['sigma'])

    def test_dilution_source_hash_nominal_and_isr(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); cut=root/'cuts.json'; dilution=root/'dilution.json'
            cut.write_text('{}')
            for variant in ('nominal','isr'):
                dilution.write_text(json.dumps({'variant':variant,'source_exclusivity_json_sha256':a.sha256_file(cut)}))
                a.validate_dilution_cut_source(dilution,cut)
                cut.write_text(cut.read_text()+' ')
                with self.assertRaises(RuntimeError): a.validate_dilution_cut_source(dilution,cut)
            dilution.write_text('{}')
            with self.assertRaises(RuntimeError): a.validate_dilution_cut_source(dilution,cut)

    def test_hybrid_endpoints_match_mean_width_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); table=self.cut_table()
            carbon=root/'carbon.json'; poly=root/'poly.json'
            c.write_cut_json_from_table(table,carbon,'carbon')
            payload=json.loads(carbon.read_text())
            for rows in payload['periods'].values():
                for row in rows:
                    row['mu_gev2'] += .005; row['sigma_gev2'] *= 1.1
                    for label,k in (('tight',1),('nominal',2),('loose',3)):
                        row[label]=[row['mu_gev2']-k*row['sigma_gev2'],row['mu_gev2']+k*row['sigma_gev2']]
            poly.write_text(json.dumps(payload))
            for path in a.write_hybrid_fit_method_cut_jsons(carbon,poly,root/'hybrids'):
                payload=json.loads(path.read_text()); self.assertTrue(payload['fixed_dilution_diagnostic'])
                self.assertNotIn('flat_rows',payload)
                for rows in payload['periods'].values():
                    for row in rows:
                        for label,k in (('tight',1),('nominal',2),('loose',3)):
                            np.testing.assert_allclose(row[label],[row['mu_gev2']-k*row['sigma_gev2'],row['mu_gev2']+k*row['sigma_gev2']])

    def test_cache_rejects_missing_changed_inputs_cuts_and_non_supersets(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(a,'_TRANSIENT_EVENT_DIR',Path(temp)):
            root=Path(temp); raw=root/'input.root'; raw.write_bytes(b'input fixture')
            identity=a._cache_identity({'su22':raw},'PhysicsEvents',{})
            cache=root/'events.npz'
            np.savez_compressed(cache,runnum=np.array([1,1,1]),period_index=np.array([0,0,0]),
                                bin_number=np.array([1,1,1]),Mx2=np.array([.75,.88,1.02]))
            cuts={(p,b):a.CutRecord(p,b,(b-1)//6,(b-1)%6,.70,1.06) for p in PERIODS for b in range(1,25)}
            with self.assertRaises(RuntimeError): a.validate_event_cache(cache,cuts,identity)
            a._write_cache_metadata(cache,cuts,identity)
            self.assertEqual(a._event_cache_path(cache),cache.resolve())
            a.validate_event_cache(cache,cuts,identity)
            tight={key:replace(value,low_gev2=.82,high_gev2=.94) for key,value in cuts.items()}
            with self.assertRaises(RuntimeError): a.validate_event_cache(cache,tight,identity)
            a.validate_event_cache(cache,tight,identity,require_superset=True)
            derived=root/'tight.npz'
            a.derive_event_cache(source_cache_path=cache,cuts=tight,cache_path=derived)
            self.assertEqual(len(a.load_event_cache(derived)['runnum']),1)
            a.validate_event_cache(derived,tight,identity)
            with self.assertRaises(RuntimeError): a.validate_event_cache(derived,cuts,identity,require_superset=True)
            changed=copy.deepcopy(identity); changed['tree']='different'
            with self.assertRaises(RuntimeError): a.validate_event_cache(cache,cuts,changed)
            raw.write_bytes(b'changed input')
            with self.assertRaises(RuntimeError): a.validate_event_cache(cache,cuts,identity)

if __name__ == '__main__':
    unittest.main()
