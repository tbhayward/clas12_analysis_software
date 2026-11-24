import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from model import method1, transfer_ratio
from determine_dilution_factor import load_charges, evaluate, resample, PERIODS, comparison_groups, group_counts


class Method1Test(unittest.TestCase):
    def test_comparison_groups_partition_original_bins(self):
        groups=comparison_groups()
        self.assertEqual([len(groups[k]) for k in groups], [4,6,8])
        for group in groups.values():
            self.assertEqual(sorted(b for ids in group for b in ids),list(range(24)))
        self.assertEqual(groups['xb_two_tprime'][0],[0,1,2])
        self.assertEqual(groups['xb_two_tprime'][1],[3,4,5])

    def test_aggregate_counts_and_replicas(self):
        n=np.arange(24)
        data={'summer2022':{'ND3':n}}
        groups=comparison_groups()
        np.testing.assert_array_equal(group_counts(data,groups['xb_only'])['summer2022']['ND3'],[15,51,87,123])
        np.testing.assert_array_equal(group_counts(data,groups['tprime_only'])['summer2022']['ND3'],[36,40,44,48,52,56])
        replicas={'summer2022':{'ND3':np.stack([n,2*n])}}
        result=group_counts(replicas,groups['xb_two_tprime'])['summer2022']['ND3']
        self.assertEqual(result.shape,(2,8))
        np.testing.assert_array_equal(result[1],2*result[0])
        self.assertEqual(result[0,0],3)
        self.assertEqual(result[0,1],12)

    def test_published_nd3_equations(self):
        r = dict(ND3=np.array([100., 150.]), C=np.array([30., 40.]),
                 CD2=np.array([40., 60.]), CH2=np.array([20., 30.]),
                 He=np.array([10., 12.]), ET=np.array([1., 2.]))
        den = r['CD2']-.170112*r['C']-.863533*r['He']+.0336457*r['ET']
        signal = r['CD2']-.765506*r['C']-.000757116*r['He']-.233737*r['ET']
        result = method1(r)
        np.testing.assert_allclose(result['dilution'], (r['ND3']-r['He'])*signal/r['ND3']/den)
        np.testing.assert_allclose(result['packing_fraction'], 1.172*.441609*(r['ND3']-r['He'])/den)

    def test_independent_material_balance(self):
        # Construct physically consistent target rates from the geometry in Table I.
        # Rounded published geometry/coefficients agree within their quoted precision.
        L, lc, ld = 5.86*.981, 1.678*1.0025, 2.686*.9797
        rc = .1494/(1.0025*1.0016**2)
        rd = .06845/(.9797*.98**2)
        ra = .05023
        bath_pf = .45
        foil, he_cm, carbon_xs, deuteron_xs = .001, .01, .3, .2
        la = bath_pf*L
        r = dict(ET=foil, He=foil+L*he_cm,
                 C=foil+(L-lc)*he_cm+lc*rc*carbon_xs,
                 CD2=foil+(L-ld)*he_cm+ld*rd*(carbon_xs+2*deuteron_xs))
        r['ND3'] = foil+(L-la)*he_cm+la*ra*((7/6)*carbon_xs+3*deuteron_xs)
        result = method1(r)
        expected_df=la*ra*3*deuteron_xs/r['ND3']
        self.assertAlmostEqual(float(result['dilution']), expected_df, delta=.002)
        self.assertAlmostEqual(float(result['packing_fraction_bath']), bath_pf, delta=.002)

    def test_transfer_and_charge_invariance(self):
        r = dict(ND3=np.full(24, 100.), C=np.full(24, 30.), CD2=np.full(24, 40.),
                 CH2=np.full(24, 20.), He=np.full(24, 10.), ET=np.full(24, 1.))
        ratio = transfer_ratio(r)
        np.testing.assert_allclose(method1(r)['dilution'], method1(r, ratio)['dilution'])
        scaled = {t: 37*v for t,v in r.items()}
        np.testing.assert_allclose(method1(r)['dilution'], method1(scaled)['dilution'])
        np.testing.assert_allclose(method1(r)['packing_fraction'], method1(scaled)['packing_fraction'])

    def test_shared_replicas_and_no_clipping(self):
        targets = dict(ND3=np.full(24, 1000), C=np.full(24, 300), CH2=np.full(24, 200),
                       He=np.full(24, 100), ET=np.full(24, 10))
        counts = {p: dict(targets) for p in PERIODS}
        counts['spring2023']['CD2']=np.full(24,400)
        charges = {p: {t:1. for t in values} for p,values in counts.items()}
        ratio, samples=resample(counts,charges,1000,17)
        self.assertEqual(samples['summer2022']['dilution'].shape,(1000,24))
        np.testing.assert_allclose(samples['summer2022']['effective_cd2_rate']/samples['fall2022']['effective_cd2_rate'],
                                   # Their CH2 counts differ, but the transfer factor is identical.
                                   samples['summer2022']['effective_cd2_rate']/ratio /
                                   (samples['fall2022']['effective_cd2_rate']/ratio))
        bad=dict(ND3=100., C=100., CD2=1., He=1., ET=0.)
        self.assertGreater(float(method1(bad)['dilution']),1.)
        bad.update(C=0., CD2=0., He=0.)
        self.assertTrue(np.isnan(method1(bad)['dilution']))

    def test_charge_sections_and_helicity_charge(self):
        text='# RGC Su22 NH3\n1,999,99,99\n# RGC Su22 Inb ND3 Runs\n2,999,3,4\n# RGC Sp23 CD2\n3,999,5,6\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'charges.csv';p.write_text(text)
            data=load_charges(p)
        self.assertNotIn(('summer2022','NH3'),data)
        self.assertEqual(data[('summer2022','ND3')][2],7.)
        self.assertEqual(data[('spring2023','CD2')][3],11.)


if __name__ == '__main__':
    unittest.main()
