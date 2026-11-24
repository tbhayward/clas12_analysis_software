import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scipy.stats import gamma
from models import fit_momentum,momentum_expected,MOM_EDGES,PROB_LOW,PROB_HIGH,momentum_components
from selection import ChannelSelection,GROUPS,VERSION,fingerprint,seal,group_for_dataset,load_selected_inputs
from kinematics import BRANCHES,calculate
from export import export_dataset,EXTRA_TYPES

class NominalTests(unittest.TestCase):
    def test_outbending_q2_cut_and_boundary_in_calibration_and_export(self):
        from candidate_cache import select,q2_mask
        a=self.arrays()
        a['Q2'][:]=np.array([1.,1.49,1.5,1.50001,2.,2.,2.,2.,2.,2.,2.,2.])
        expected=np.array([False,False,False,True,True,True,True,True,True,True,True,True])
        np.testing.assert_array_equal(select(a,'RGB_fall2019_D2')[0],expected)
        self.assertTrue(np.all(select(a,'RGB_spring2019_D2')[0]))
        self.assertTrue(np.all(select(a,'RGC_fall2022_ND3')[0]))
        with tempfile.TemporaryDirectory() as tmp:
            selection=ChannelSelection(self.cut_file(Path(tmp)/'cuts.json'))
            np.testing.assert_array_equal(selection.evaluate(a,'RGB_fall2019_D2')[0],expected)
            self.assertTrue(np.all(selection.evaluate(a,'RGB_spring2020_D2')[0]))
        np.testing.assert_array_equal(q2_mask(np.array([1.5,1.50001]),'RGB_fall2019_D2'),[False,True])

    def test_gamma_signal_recovery_and_falling_tail(self):
        p=np.array([50000,.08,5.,12000,.17])
        h=np.random.default_rng(72).poisson(momentum_expected(MOM_EDGES,p))
        f=fit_momentum(h,MOM_EDGES)
        self.assertTrue(f['usable'],f)
        self.assertAlmostEqual(f['mode'],p[1],delta=.008)
        self.assertAlmostEqual(f['shape'],p[2],delta=.7)
        scale=f['parameters'][1]/(f['parameters'][2]-1)
        fraction=gamma.cdf(f['cut_high'],f['shape'],scale=scale)-gamma.cdf(f['cut_low'],f['shape'],scale=scale)
        self.assertAlmostEqual(fraction,PROB_HIGH-PROB_LOW,places=10)
        _,tail=momentum_components(np.array([.25,.5,1.]),p,.005)
        self.assertTrue(tail[0]>tail[1]>tail[2])

    @staticmethod
    def arrays():
        n=12;a={k:np.zeros(n,dtype=float) for k in BRANCHES}
        for k,v in dict(beam_energy=10.55,x=.30,tprime=-.35,Q2=2,W=3,y=.5,fiducial_status=3,detector1=1,detector2=1,e_p=3.3,e_theta=.3,e_phi=.2,p1_p=.6,p1_theta=1.,p1_phi=2.,p2_p=2.,p2_theta=.2,p2_phi=1.).items():a[k][:]=v
        a['p2_p']+=np.arange(n)*.01
        a['Mx2']=calculate(a)['mx2_n']
        a.update(runnum=np.full(n,16900,dtype=np.int32),evnum=np.arange(n,dtype=np.int32),helicity=np.where(np.arange(n)%2,1,-1).astype(np.int32),beam_pol=np.linspace(.8,.9,n),target_pol=np.linspace(-.2,.2,n),phi2=np.linspace(0,6,n))
        for key in ['DepA','DepB','DepC','DepV','DepW']:a[key]=np.linspace(.3,.7,n)
        return a

    @staticmethod
    def cut_file(path):
        groups={}
        for name,members in GROUPS.items():
            bins=[dict(bin=i+1,active=i!=8,mass_low=-100.,mass_high=100.,momentum_low=0.,momentum_high=100.) for i in range(24)]
            groups[name]=dict(calibration_members=members,bins=bins)
        p=dict(version=VERSION,base_selection_sha256=fingerprint(),partial=False,groups=groups)
        p['cutset_id']=seal(p);path.write_text(json.dumps(p));return path

    def test_control_targets_share_exact_period_windows(self):
        self.assertEqual(group_for_dataset('RGC_fall2022_C'),group_for_dataset('RGC_spring2023_ND3'))
        self.assertEqual(group_for_dataset('RGC_spring2023_CD2'),group_for_dataset('RGC_fall2022_ND3'))
        self.assertNotEqual(group_for_dataset('RGC_summer2022_C'),group_for_dataset('RGC_fall2022_C'))
        with self.assertRaises(ValueError):group_for_dataset('invented_target')
        with tempfile.TemporaryDirectory() as tmp:
            s=ChannelSelection(self.cut_file(Path(tmp)/'cuts.json'));a=self.arrays()
            mask,_=s.evaluate(a,'RGC_fall2022_ND3');control,_=s.evaluate(a,'RGC_spring2023_C')
            np.testing.assert_array_equal(mask,control)
            self.assertFalse(s.candidate_mask(np.array([.9]),np.array([.1]),np.array([8]),'RGC_fall2022_C')[0])

    def test_saved_missing_energy_fast_path_matches_exact_vectors(self):
        from kinematics import MN,ME,MP,MPI,track
        with tempfile.TemporaryDirectory() as tmp:
            s=ChannelSelection(self.cut_file(Path(tmp)/'cuts.json'));a=self.arrays()
            exact,annotations=s.evaluate(a,'RGC_fall2022_ND3')
            a['Emiss']=a['beam_energy']+MN-track(a,'e',ME)[:,0]-track(a,'p1',MP)[:,0]-track(a,'p2',MPI)[:,0]
            fast,extra=s.evaluate(a,'RGC_fall2022_ND3')
            np.testing.assert_array_equal(exact,fast)
            for k in annotations:np.testing.assert_allclose(annotations[k],extra[k],atol=1e-10)
            a['Emiss']+=.1
            with self.assertRaises(ValueError):s.evaluate(a,'RGC_fall2022_ND3')

    def test_reference_assisted_windows_preserve_provenance(self):
        from reference import transfer,adjustment
        from mass_model import MASS_EDGES,expected_bins
        nominal=dict(usable=False,reason='weak Su22 peak',warnings=[])
        ref=dict(usable=True,mu=.95,sigma=.10,cut_low=.75,cut_high=1.15,parameters=[100,.95,.1,5,5,5],kind='mass',warnings=[],fit_range=[.7,1.3])
        own=[dict(ref,mu=.93,sigma=.11) for _ in range(24)]
        correction=adjustment(own,[ref]*24,'mass')
        h=expected_bins(MASS_EDGES,[1,.93,.11,3,3,3],(.7,1.3))
        borrowed=transfer(h,MASS_EDGES,'mass',ref,correction,nominal,6)
        self.assertTrue(borrowed['usable'])
        self.assertEqual(borrowed['origin'],'reference_assisted')
        self.assertAlmostEqual(borrowed['cut_low'],.71)
        self.assertAlmostEqual(borrowed['cut_high'],1.15)
        self.assertEqual(borrowed['reference']['bin'],6)
        self.assertEqual(borrowed['independent_fit']['reason'],'weak Su22 peak')

    def test_tampered_partial_and_incompatible_cuts_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=self.cut_file(Path(tmp)/'cuts.json');p=json.loads(path.read_text());p['partial']=True;path.write_text(json.dumps(p))
            with self.assertRaisesRegex(ValueError,'checksum'):ChannelSelection(path)
            p.pop('cutset_id');p['cutset_id']=seal(p);path.write_text(json.dumps(p))
            with self.assertRaisesRegex(ValueError,'Partial'):ChannelSelection(path)

    def test_full_branch_metadata_preservation_empty_output_and_resume(self):
        import uproot
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);a=self.arrays();a['fiducial_status'][0]=0
            a['x'][1]=.65
            source=out/'input.root'
            with uproot.recreate(source) as f:
                f.mktree('PhysicsEvents',{k:v.dtype for k,v in a.items()});f['PhysicsEvents'].extend(a)
            cuts=self.cut_file(out/'cuts.json');s=ChannelSelection(cuts);mask,extra=s.evaluate(a,'RGC_fall2022_ND3')
            job=('RGC_fall2022_ND3',str(source),'PhysicsEvents',str(cuts),str(out/'output'),3,None)
            result=export_dataset(job)
            self.assertEqual(result['selected_rows'],10)
            with uproot.open(result['selected_file']) as f:
                selected=f['PhysicsEvents'].arrays(library='np')
                for k in a:np.testing.assert_array_equal(selected[k],a[k][mask],err_msg=k)
                for k in extra:np.testing.assert_array_equal(selected[k],extra[k])
                self.assertEqual(str(f['selection_cutset_id']),s.cutset_id)
            registry=out/'selected_inputs.json'
            registry.write_text(json.dumps(dict(cut_file='cuts.json',cutset_id=s.cutset_id,tree='PhysicsEvents',files={'RGC_fall2022_ND3':dict(path='output/selected/selected_RGC_fall2022_ND3.root',selected_rows=10)})))
            loaded=load_selected_inputs(registry)
            self.assertEqual(loaded['files']['RGC_fall2022_ND3']['path'],result['selected_file'])
            stamp=Path(result['selected_file']).stat().st_mtime_ns;export_dataset(job)
            self.assertEqual(stamp,Path(result['selected_file']).stat().st_mtime_ns)
            a['fiducial_status'][:]=0;a['target_pol'][:]=0
            with uproot.recreate(source) as f:
                f.mktree('PhysicsEvents',{k:v.dtype for k,v in a.items()});f['PhysicsEvents'].extend(a)
            empty=export_dataset(('RGC_spring2023_C',str(source),'PhysicsEvents',str(cuts),str(out/'output'),3,None))
            self.assertEqual(empty['selected_rows'],0)
            with uproot.open(empty['selected_file']) as f:self.assertEqual(len(f['PhysicsEvents'].keys()),len(a)+len(EXTRA_TYPES))

if __name__=='__main__':unittest.main()
