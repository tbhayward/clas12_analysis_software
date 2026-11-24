import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from fit_models import *
from run_selection import histogram,windows

class SelectionTests(unittest.TestCase):
    def test_mass_and_momentum_fits_recover_known_windows(self):
        for kind,edges,parameters in [('mass',MASS_EDGES,[450,.91,.075,60,75,100]),('momentum',MOM_EDGES,[350,.075,.025,10,15,5])]:
            x=(edges[:-1]+edges[1:])/2
            inside=(x>=FIT_RANGES[kind][0])&(x<FIT_RANGES[kind][1])
            g,b=components(x[inside],parameters,FIT_RANGES[kind])
            h=np.zeros(len(x));h[inside]=np.random.default_rng(2026).poisson(g+b)
            fit=fit_histogram(h,edges,kind)
            self.assertTrue(fit['usable'],fit)
            self.assertAlmostEqual(fit['mu'],parameters[1],delta=.006)
            self.assertAlmostEqual(fit['sigma'],parameters[2],delta=.006)
            self.assertAlmostEqual(fit['cut_high'],parameters[1]+2*parameters[2],delta=.015)

    def test_no_peak_does_not_get_cut(self):
        self.assertFalse(fit_histogram(np.zeros(200),MASS_EDGES,'mass')['usable'])
        self.assertFalse(fit_histogram(np.full(200,100),MASS_EDGES,'mass')['usable'])
        f=[dict(usable=False)]*24
        self.assertFalse(windows(np.array([.9]),np.array([0]),f)[0])

    def test_histogram_and_window_edge_accounting(self):
        h=np.zeros((24,2),dtype=int)
        histogram(h,np.array([0,.5,1.,np.nan]),np.array([0,0,0,0]),np.array([0,.5,1.]))
        np.testing.assert_array_equal(h[0],[1,1])
        f=[dict(usable=True,cut_low=.8,cut_high=1.0)]*24
        np.testing.assert_array_equal(windows(np.array([.8,1.,1.1]),np.zeros(3,dtype=int),f),[True,True,False])

    def test_polar_angle_guide_is_not_a_cut(self):
        x=(ANGLE_EDGES[:-1]+ANGLE_EDGES[1:])/2
        f=fit_angle(500*np.exp(-.5*((x-.003)/.016)**2))
        self.assertTrue(f['usable'])
        self.assertAlmostEqual(f['mu'],.003,places=4)
        self.assertNotIn('cut_low',f)

if __name__=='__main__':unittest.main()
