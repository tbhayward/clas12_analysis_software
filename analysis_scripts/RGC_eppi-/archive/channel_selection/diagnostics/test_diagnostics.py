import json
from pathlib import Path
import tempfile
import shutil
import subprocess
import sys
import unittest
import numpy as np
import uproot
from kinematics import *
from run_diagnostics import process, DEFAULT_PERIODS


def boost(v, beta):
    b2 = np.dot(beta, beta)
    gamma = 1/np.sqrt(1-b2)
    bp = np.dot(beta, v[1:])
    return np.r_[gamma*(v[0]+bp), v[1:]+((gamma-1)*bp/b2+gamma*v[0])*beta]


def fixture():
    rows = []
    for xi in range(4):
        for ti in range(6):
            xb, tp = (XB[xi]+XB[xi+1])/2, (TP[ti]+TP[ti+1])/2
            eb, q2 = 10.55, (2.1 if xi < 3 else 4.0)
            nu = q2/(2*MN*xb)
            ee = eb-nu
            pe, pb = np.sqrt(ee**2-ME**2), np.sqrt(eb**2-ME**2)
            theta = np.arccos((eb*ee-ME**2-q2/2)/(pb*pe))
            electron = np.array([ee, pe*np.sin(theta), 0, pe*np.cos(theta)])
            q = np.array([eb, 0, 0, pb])-electron
            total = q+np.array([MN, 0, 0, 0])
            w = np.sqrt(total[0]**2-np.dot(total[1:], total[1:]))
            beta = total[1:]/total[0]
            qs = boost(q, -beta)
            epi = (w*w+MPI*MPI-MP*MP)/(2*w)
            pp = np.sqrt(epi*epi-MPI*MPI)
            cost = 1-tp/(2*np.linalg.norm(qs[1:])*pp)
            assert -1 < cost < 1
            z = qs[1:]/np.linalg.norm(qs[1:])
            transverse = np.cross(z, [0, 1, 0])
            transverse /= np.linalg.norm(transverse)
            vec = pp*(cost*z+np.sqrt(1-cost*cost)*transverse)
            pion = boost(np.r_[epi, vec], beta)
            proton = boost(np.r_[np.sqrt(pp*pp+MP*MP), -vec], beta)
            row = dict(beam_energy=eb, x=xb, tprime=-tp, Q2=q2, W=w, y=nu/eb,
                       Mx2=0., fiducial_status=3, detector1=1+ti%2, detector2=1)
            for prefix, v in [('e', electron), ('p1', proton), ('p2', pion)]:
                p = np.linalg.norm(v[1:])
                row.update({prefix+'_p': p, prefix+'_theta': np.arccos(v[3]/p), prefix+'_phi': np.arctan2(v[2], v[1])})
            rows.append(row)
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


class DiagnosticTests(unittest.TestCase):
    def test_exclusive_conservation_and_residual_identity(self):
        a = fixture()
        v = calculate(a)
        for k in ['mx2_n', 'delta_e', 'delta_p', 'delta_theta', 'delta_phi', 'pred_mass', 'ptmiss']:
            np.testing.assert_allclose(v[k], 0, atol=1e-10, err_msg=k)
        # For an initially stationary struck neutron, the deuteron reference leaves energy MD-MN.
        np.testing.assert_allclose(v['mx2_d'], (MD-MN)**2, atol=1e-10)
        a['p2_p'] *= 1.02
        v = calculate(a)
        np.testing.assert_allclose(v['mx2_n'], v['delta_e']**2-v['delta_px']**2-v['delta_py']**2-v['delta_pz']**2, atol=1e-12)
        np.testing.assert_allclose(v['delta_p'], a['p2_p']/1.02*.02, atol=1e-10)

    def test_default_output_is_beside_script_from_another_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scripts, cwd = root/'scripts', root/'elsewhere'
            scripts.mkdir()
            cwd.mkdir()
            for name in ['run_diagnostics.py', 'kinematics.py']:
                shutil.copy2(Path(__file__).with_name(name), scripts/name)
            period = scripts/'output'/DEFAULT_PERIODS[0]
            period.mkdir(parents=True)
            (period/'manifest.json').write_text(json.dumps(dict(counts=dict(in_bins=123), partial=False)))
            (period/'plot_status.json').write_text(json.dumps(dict(stage='complete')))
            opening = scripts/'output'/DEFAULT_PERIODS[1]
            opening.mkdir()
            (opening/'progress.json').write_text(json.dumps(dict(stage='opening input')))
            result = subprocess.run([sys.executable, str(scripts/'run_diagnostics.py'), '--status'], cwd=cwd, capture_output=True, text=True, check=True)
            self.assertIn('complete, 123 combinations', result.stdout)
            self.assertIn('opening input', result.stdout)
            self.assertFalse((cwd/'output').exists())

    def test_bin_edges(self):
        np.testing.assert_array_equal(bin_indices(np.array([.10, .25, .60, np.nan]), np.array([.05, .25, .5, .5])), [0, 7, -1, -1])

    def test_six_period_streaming_root_and_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            a = fixture()
            # Explicit rejection accounting, independent of histograms.
            a = {k: np.concatenate([v, v[:3]]) for k, v in a.items()}
            a['fiducial_status'][-3] = 0
            a['x'][-2] = .9
            a['e_p'][-1] = np.nan
            for name in DEFAULT_PERIODS:
                source, out = tmp/(name+'.root'), tmp/name
                with uproot.recreate(source) as f:
                    f.mktree('PhysicsEvents', {k: v.dtype for k, v in a.items()})
                    f['PhysicsEvents'].extend(a)
                options = dict(max_entries=None, all_fiducials=False, tree='PhysicsEvents', step_size=7, force=False)
                process((name, source, out, options))
                m = json.loads((out/'manifest.json').read_text())
                self.assertEqual(m['counts'], dict(read=27, nonfinite_input=1, quality_rejected=1, outside_bins=1, in_bins=24))
                stamp = (out/'histograms.npz').stat().st_mtime_ns
                process((name, source, out, options))
                self.assertEqual(stamp, (out/'histograms.npz').stat().st_mtime_ns)
                with np.load(out/'histograms.npz') as h:
                    for k in METRICS:
                        np.testing.assert_array_equal(h['h_'+k][0], h['h_'+k][1]+h['h_'+k][2])
                        np.testing.assert_array_equal(h['s_'+k][0, :, 0], np.ones(24))
            # Wrong reference mass must be detected rather than silently accepted.
            a['Mx2'][:] = .5
            with uproot.recreate(source) as f:
                f.mktree('PhysicsEvents', {k: v.dtype for k, v in a.items()})
                f['PhysicsEvents'].extend(a)
            with self.assertRaisesRegex(ValueError, 'disagree'):
                process((name, source, out, options))

if __name__ == '__main__':
    unittest.main()
