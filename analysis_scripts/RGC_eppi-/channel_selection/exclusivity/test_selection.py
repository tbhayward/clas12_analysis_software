import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from fit_models import *
from run_selection import histogram,windows,scan
from candidate_cache import derive,build_cache,iterate
from kinematics import BRANCHES,calculate

class SelectionTests(unittest.TestCase):
    def test_mass_and_momentum_recovery_with_checks(self):
        for kind,edges,pars in [('mass',MASS_EDGES,[450,.91,.075,60,80,100]),('momentum',MOM_EDGES,[350,.075,.025,5,7,9])]:
            h=np.random.default_rng(2026).poisson(expected_bins(edges,pars,(edges[0],edges[-1])))
            fit=fit_with_checks(h,edges,kind)
            self.assertTrue(fit['usable'],fit)
            self.assertAlmostEqual(fit['mu'],pars[1],delta=.006)
            self.assertAlmostEqual(fit['sigma'],pars[2],delta=.006)
            self.assertIn('concave_quadratic',fit['checks'])
            self.assertIn('coarser_histogram',fit['checks'])

    def test_bin_integration_matches_numerical_integral(self):
        edges=np.linspace(.80,1.0,11)
        pars=[40,.91,.012,2,4,6]
        expected=[]
        for lo,hi in zip(edges[:-1],edges[1:]):
            x=np.linspace(lo,hi,10001)
            g,b=components(x,pars,(.8,1.0))
            y=g+b
            expected.append(np.sum((y[1:]+y[:-1])*.5*np.diff(x))/(hi-lo))
        np.testing.assert_allclose(expected_bins(edges,pars,(.8,1.0)),expected,rtol=1e-7)

    def test_no_peak_does_not_get_cut(self):
        self.assertFalse(fit_histogram(np.zeros(len(MASS_EDGES)-1),MASS_EDGES,'mass')['usable'])
        self.assertFalse(fit_histogram(np.full(len(MASS_EDGES)-1,100),MASS_EDGES,'mass')['usable'])
        self.assertFalse(windows(np.array([.9]),np.array([0]),[dict(usable=False)]*24)[0])

    def test_backgrounds_cannot_curve_up(self):
        for model,free in [('linear',[100,.9,.07,20,80]),('concave_quadratic',[100,.9,.07,20,80,40])]:
            pars=canonical(free,model)
            _,bg=components(np.linspace(.7,1.3,101),pars,(.7,1.3))
            self.assertTrue(np.all(bg>=0))
            self.assertLessEqual(np.max(np.diff(bg,n=2)),1e-10)

    def test_range_expansion_covers_cut(self):
        h=expected_bins(MASS_EDGES,[400,.9,.13,20,20,20],(.3,1.6))
        f=adaptive_fit(h,MASS_EDGES,'mass')
        self.assertTrue(f['usable'],f)
        self.assertGreater(len(f['range_history']),1)
        self.assertLessEqual(f['fit_range'][0],f['cut_low'])
        self.assertGreaterEqual(f['fit_range'][1],f['cut_high'])

    def test_histogram_and_window_edges(self):
        h=np.zeros((24,2),dtype=int)
        histogram(h,np.array([0,.5,1.,np.nan]),np.zeros(4,dtype=int),np.array([0,.5,1.]))
        np.testing.assert_array_equal(h[0],[1,1])
        f=[dict(usable=True,cut_low=.8,cut_high=1.)]*24
        np.testing.assert_array_equal(windows(np.array([.8,1.,1.1]),np.zeros(3,dtype=int),f),[True,True,False])
        rebinned,edges=merge_pairs(np.array([1,2,3]),np.arange(4))
        np.testing.assert_array_equal(rebinned,[3,3])
        np.testing.assert_array_equal(edges,[0,2,3])

    def test_degree_guide_is_not_a_cut(self):
        x=(ANGLE_EDGES[:-1]+ANGLE_EDGES[1:])/2
        f=fit_angle(500*np.exp(-.5*((x-.17)/.9)**2))
        self.assertTrue(f['usable'])
        self.assertAlmostEqual(f['mu'],.17,places=4)
        self.assertEqual(f['unit'],'degree')
        self.assertNotIn('cut_low',f)

    def test_minimal_kinematics_match_original(self):
        a=self.fixture()
        expected=calculate(a)
        mass,momentum,angle=derive(a)
        np.testing.assert_allclose(mass,expected['mx2_d'],atol=1e-12)
        np.testing.assert_allclose(momentum,np.sqrt(expected['delta_px']**2+expected['delta_py']**2+expected['delta_pz']**2),atol=1e-12)
        np.testing.assert_allclose(angle,expected['delta_theta'],atol=1e-12)

    @staticmethod
    def fixture():
        n=9
        a={k:np.zeros(n) for k in BRANCHES}
        for k,v in dict(beam_energy=10.55,x=.30,tprime=-.35,Q2=2,W=3,y=.5,fiducial_status=3,detector1=1,detector2=1,e_p=3.3,e_theta=.3,e_phi=.2,p1_p=.6,p1_theta=1.0,p1_phi=2.,p2_p=2.,p2_theta=.2,p2_phi=1.).items():a[k][:]=v
        a['p2_p']+=np.arange(n)*.01
        a['Mx2']=calculate(a)['mx2_n']
        return a

    def test_root_scanned_once_and_cached_selection_is_exact(self):
        import uproot
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            a=self.fixture()
            a['fiducial_status'][0]=0
            a['x'][1]=.8
            a['e_phi'][2]=np.nan
            with uproot.recreate(out/'source.root') as f:
                f.mktree('PhysicsEvents',{k:v.dtype for k,v in a.items()});f['PhysicsEvents'].extend(a)
            options=dict(max_entries=None,step_size=3,rebuild_cache=False)
            original=uproot.open
            with patch('uproot.open',wraps=original) as opened:
                cache,meta=build_cache(out/'source.root','PhysicsEvents',out,options)
                self.assertEqual(opened.call_count,1)
                cache2,meta2=build_cache(out/'source.root','PhysicsEvents',out,options)
                self.assertEqual(opened.call_count,1)
            self.assertEqual(meta['combinations'],6)
            records=np.concatenate(list(iterate(cache,meta)))
            self.assertEqual(records.dtype.itemsize,26)
            mass=[dict(usable=False)]*24
            mass[7]=dict(usable=True,cut_low=float(np.min(records['mass'])),cut_high=float(np.max(records['mass'])))
            momentum=[dict(usable=True,cut_low=0,cut_high=100)]*24
            _,counts,_=scan(cache,meta,'angle',out,mass,momentum)
            self.assertEqual(counts['selected'][7],6)
            self.assertEqual(counts['FD'][7],6)
            self.assertEqual(counts['CD'][7],0)
            _,excluded,_=scan(cache,meta,'angle',out,[dict(usable=False)]*24,momentum)
            self.assertEqual(excluded['selected'].sum(),0)

if __name__=='__main__':unittest.main()
