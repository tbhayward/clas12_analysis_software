#!/usr/bin/env python3
"""Stage-I/II DVCS BSA-only CFF diagnostic, independent of production scripts.

r multiplies ONLY Im(H) of the KM15 CFF at each measured kinematic point.
Optional slope adds exp[b*(t-t_ref)] to test shape rigidity. Other CFFs
remain KM15. This is a model-relative diagnostic, NOT a local CFF extraction.

The code refuses to fit unless a numerical Gepard preflight verifies that
r=1 exactly reproduces unmodified KM15 helicity cross sections, that varying
r changes the BSA, and that restoration returns the original prediction.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
M = 0.9382720813


def parse_tuple(value):
    import ast
    try:
        v = ast.literal_eval(str(value))
        if isinstance(v, (tuple, list)) and len(v) >= 2:
            return float(v[0]), float(v[1])
    except (ValueError, SyntaxError, TypeError):
        pass
    return np.nan, np.nan
#enddef


def get_col(row, choices):
    for c in choices:
        if c in row.index:
            try:
                v = float(row[c])
                if np.isfinite(v):
                    return v
            except (TypeError, ValueError):
                pass
            #endtry
        #endif
    #endfor
    return np.nan
#enddef


def read_pass1(path):
    a = pd.read_csv(path, sep=r'\s+', comment='#', header=None,
                    names=['phi','Q2','xB','t','E','A','err'])
    a['phi_deg'] = np.rad2deg(a.phi)
    a['dataset'] = 'published_pass1'
    a['error_type'] = 'published_sigA'
    return a[['dataset','E','Q2','xB','t','phi_deg','A','err','error_type']]
#enddef


def read_sebastian(path, q2_lookup=None):
    """PDF extraction has phi and A, but NO event-level Q2/xB/t per point.

    For now do not invent means: require a supplied explicit bin-kinematics
    table, with rga_bin,Q2,xB,t columns. It must be documented as proxy.
    """
    if q2_lookup is None:
        print('[skip] Sebastian: no independently verified per-bin Q2/xB/t map. '
              'Supply --sebastian-kinematics CSV to enable; no fabricated means.')
        return pd.DataFrame()
    #endif
    a = pd.read_csv(path)
    a = a[(a.series == 'sebastian_pass2') & (a.topology == 'egamma_p')].copy()
    k = pd.read_csv(q2_lookup)
    required = {'rga_bin','Q2','xB','t'}
    if not required.issubset(k.columns) or k.rga_bin.duplicated().any():
        raise ValueError('Sebastian kinematics must have unique rga_bin,Q2,xB,t columns')
    #endif
    a = a.merge(k[['rga_bin','Q2','xB','t']], on='rga_bin', how='left', validate='many_to_one')
    a['dataset'] = 'sebastian_' + a.torus
    a['E'] = 10.6
    a['A'] = a.bsa
    a['err'] = a.bsa_error_graphical
    a['error_type'] = 'PDF_graphical_error'
    return a[['dataset','E','Q2','xB','t','phi_deg','A','err','error_type']]
#enddef


def read_pass2(path, label):
    """Use native pass-2 rows; require REAL per-period mean kinematics.

    The production CSV can have per-period Q2/xB/t mean columns with variant
    spellings; explicit --pass2-kinematics-cols can be added after schema audit.
    Never substitute cell centers silently.
    """
    df = pd.read_csv(path, low_memory=False)
    bcol = f'BSA, counts, {label}'
    if bcol not in df:
        raise ValueError(f'No {bcol} in {path}')
    #endif
    output = []
    missing_kin = 0
    for _, row in df.iterrows():
        A, stat = parse_tuple(row[bcol])
        ptp = get_col(row, [f'BSA, counts, {label}, point-to-point sys'])
        phi = get_col(row, [f'phiavg, {label}'])
        if not np.isfinite(phi):
            lo, hi = get_col(row, ['phimin']), get_col(row, ['phimax'])
            phi = (lo + hi)/2 if np.isfinite(lo+hi) else np.nan
        #endif
        if not (np.isfinite(A) and np.isfinite(stat) and stat > 0 and
                np.isfinite(ptp) and ptp >= 0 and np.isfinite(phi)):
            continue
        #endif
        q2 = get_col(row, [f'Q2avg, {label}', f'Q2_avg, {label}', f'mean Q2, {label}'])
        xb = get_col(row, [f'xBavg, {label}', f'xB_avg, {label}', f'mean xB, {label}'])
        t = get_col(row, [f'tavg, {label}', f't_avg, {label}', f'mean t, {label}'])
        if not np.isfinite(t):
            tabs = get_col(row, [f't_abs_avg, {label}', f'tabsavg, {label}', f'mean -t, {label}'])
            t = -abs(tabs) if np.isfinite(tabs) else np.nan
        #endif
        if not np.isfinite(q2+xb+t):
            missing_kin += 1
            continue
        #endif
        output.append(dict(dataset='our_pass2_'+label.replace(' ','_'), E=10.2 if label=='Sp19 Inb' else 10.6,
                           Q2=q2,xB=xb,t=t,phi_deg=phi,A=A,err=math.hypot(stat,ptp),
                           error_type='stat_plus_ptp'))
    #endfor
    print(f'[pass2 {label}] accepted {len(output)}, missing mean kinematics {missing_kin}')
    if not output:
        print('[pass2] No valid rows. Inspect CSV mean-kinematics column names; '
              'this script deliberately does not use bin centers as measured means.')
    #endif
    return pd.DataFrame(output)
#enddef


class GepardEngine:
    def __init__(self):
        import gepard as g
        from gepard.fits import th_KM15
        self.g, self.th = g, th_KM15
        self.model = th_KM15.m
        if not hasattr(self.model, 'cff'):
            raise RuntimeError('Gepard KM15 model has no cff(pt) method; cannot override safely.')
        #endif
        self._original_cff = self.model.cff
        self._original_imh = self.model.ImH
        self._r, self._b, self._tref = 1., 0., -0.3
    #enddef

    def point(self, row, helicity):
        phi = math.pi - math.radians(float(row.phi_deg))  # established project convention
        p = self.g.DataPoint(xB=float(row.xB), t=float(row.t), Q2=float(row.Q2),
                             phi=phi, observable='XS', frame='trento',
                             process='ep2epgamma', exptype='fixed target',
                             in1energy=float(row.E), in1charge=-1,
                             in1polarization=int(helicity))
        p.prepare()
        return p
    #enddef

    def _modified_cff(self, pt):
        c = np.asarray(self._original_cff(pt), dtype=float).copy()
        if c.shape != (8,):
            raise RuntimeError(f'Unexpected CFF vector shape: {c.shape}')
        #endif
        c[1] *= self._r * math.exp(self._b*(float(pt.t)-self._tref))
        return c
    #enddef

    def _modified_imh(self, pt, *args, **kwargs):
        return self._original_imh(pt, *args, **kwargs) * self._r * math.exp(self._b*(float(pt.t)-self._tref))
    #enddef

    def predict(self, row, r=1., b=0., tref=-0.3):
        self._r, self._b, self._tref = float(r), float(b), float(tref)
        if self._r <= 0:
            return np.nan
        #endif
        # Some Gepard versions use model.cff(), others individual ImH().
        # Override both, then restore in finally. Never mutate KM15 globally.
        try:
            self.model.cff = self._modified_cff
            self.model.ImH = self._modified_imh
            plus = float(self.th.predict(self.point(row, +1)))
            minus = float(self.th.predict(self.point(row, -1)))
        finally:
            self.model.cff = self._original_cff
            self.model.ImH = self._original_imh
        #endtry
        if not np.isfinite(plus+minus) or min(plus, minus) <= 0:
            return np.nan
        #endif
        return (plus-minus)/(plus+minus)
    #enddef

    def baseline(self, row):
        plus = float(self.th.predict(self.point(row, +1)))
        minus = float(self.th.predict(self.point(row, -1)))
        return (plus-minus)/(plus+minus)
    #enddef

    def preflight(self, row):
        base = self.baseline(row)
        unity = self.predict(row, 1.)
        changed = self.predict(row, 1.35)
        restored = self.baseline(row)
        print(f'[preflight] KM15={base:+.9g}, override(r=1)={unity:+.9g}, '
              f'override(r=1.35)={changed:+.9g}, restored={restored:+.9g}')
        if not all(np.isfinite(x) for x in [base, unity, changed, restored]):
            raise RuntimeError('Preflight nonfinite BSA')
        #endif
        if abs(unity-base) > 1e-7 or abs(restored-base) > 1e-7:
            raise RuntimeError('Preflight failed KM15 identity/restoration')
        #endif
        if abs(changed-base) < 1e-7:
            raise RuntimeError('Preflight failed: ImH variation did not change BSA. '
                               'Check Gepard model interface and helicity convention.')
        #endif
        return dict(baseline=base, unity=unity, varied=changed, restored=restored)
    #enddef
#endclass


def fit_dataset(engine, data, nrep, seed, fit_slope, max_points, tref, outdir):
    data = data.sort_values(['Q2','xB','t','phi_deg']).copy()
    if max_points and len(data) > max_points:
        # Deterministic subsample spanning kinematic domain; exploratory only.
        data = data.iloc[np.linspace(0,len(data)-1,max_points,dtype=int)].copy()
        print(f'[warning] subsampling to {len(data)} points for pilot runtime')
    #endif
    rows = list(data.itertuples(index=False))
    obs = data.A.to_numpy(float)
    err = data.err.to_numpy(float)
    def predict(theta):
        r = float(theta[0]); b = float(theta[1]) if fit_slope else 0.
        return np.asarray([engine.predict(row,r,b,tref) for row in rows])
    #enddef
    def objective(theta, target):
        theory = predict(theta)
        if not np.all(np.isfinite(theory)):
            return 1e15
        #endif
        return float(np.sum(((target-theory)/err)**2))
    #enddef
    bounds = [(0.1,3.0)] + ([(-3.,3.)] if fit_slope else [])
    def solve(target, initial=None):
        if fit_slope:
            res = minimize(lambda th:objective(th,target),
                           np.asarray(initial if initial is not None else [1.,0.]),
                           bounds=bounds, method='L-BFGS-B')
            return np.asarray(res.x), float(res.fun), bool(res.success)
        #endif
        res = minimize_scalar(lambda r: objective([r],target),bounds=bounds[0],method='bounded',
                              options={'xatol':1e-4})
        return np.array([res.x]),float(res.fun),bool(res.success)
    #enddef
    best, chi2, ok = solve(obs)
    if not ok:
        raise RuntimeError('Nominal fit failed')
    #endif
    rng = np.random.default_rng(seed)
    reps = []
    for j in range(nrep):
        theta, value, success = solve(obs + rng.normal(size=len(obs))*err, best)
        reps.append(dict(replica=j, r_ImH=float(theta[0]),
                         b_t=float(theta[1]) if fit_slope else 0.,
                         chi2=value, valid=success))
        if (j+1)%25==0:
            print(f'  replicas {j+1}/{nrep}',flush=True)
        #endif
    #endfor
    pd.DataFrame(reps).to_csv(outdir/'replicas.csv',index=False)
    pred = predict(best)
    data['theory'] = pred
    data['residual'] = obs-pred
    data.to_csv(outdir/'point_predictions.csv',index=False)
    good = pd.DataFrame(reps).query('valid')
    summary = dict(n_points=len(data),r_ImH=float(best[0]),
                   b_t=float(best[1]) if fit_slope else 0.,
                   chi2=chi2,ndf=len(data)-len(best),
                   replica_r_16=float(good.r_ImH.quantile(.16)) if len(good) else None,
                   replica_r_50=float(good.r_ImH.quantile(.50)) if len(good) else None,
                   replica_r_84=float(good.r_ImH.quantile(.84)) if len(good) else None,
                   replica_valid=len(good),
                   r_boundary=bool(abs(best[0]-bounds[0][0])<.001 or abs(best[0]-bounds[0][1])<.001))
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
#enddef


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pass1',type=Path,default=HERE/'import/rga_pass1_bsa.txt')
    ap.add_argument('--sebastian',type=Path,default=HERE/'import/sebastian_pass2_bsa_from_note.csv')
    ap.add_argument('--sebastian-kinematics',type=Path,default=None,
                    help='Verified bin means CSV (rga_bin,Q2,xB,t); required to include Sebastian')
    ap.add_argument('--pass2',type=Path,default=HERE.parent/'output/csvs/dvcs_pass2_analysis.csv')
    ap.add_argument('--pass2-label',default='10.6 GeV Unpol',help='BSA column label in native Pass2 CSV')
    ap.add_argument('--output',type=Path,default=HERE/'output/bsa_cff_crosscheck')
    ap.add_argument('--datasets',default='pass1',help='Comma-separated: pass1,sebastian,pass2')
    ap.add_argument('--replicas',type=int,default=100)
    ap.add_argument('--max-points',type=int,default=100,
                    help='Pilot limit per dataset (0=all); deterministic subsample')
    ap.add_argument('--seed',type=int,default=20261012)
    ap.add_argument('--fit-t-slope',action='store_true',help='r(t)=r0*exp[b*(t-tref)]')
    ap.add_argument('--t-ref',type=float,default=-0.3)
    ap.add_argument('--preflight-only',action='store_true')
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    sets={}
    requested={x.strip() for x in args.datasets.split(',')}
    if 'pass1' in requested:
        sets['pass1']=read_pass1(args.pass1)
    #endif
    if 'sebastian' in requested:
        sets['sebastian']=read_sebastian(args.sebastian,args.sebastian_kinematics)
    #endif
    if 'pass2' in requested:
        sets['pass2']=read_pass2(args.pass2,args.pass2_label)
    #endif
    sets={k:v for k,v in sets.items() if len(v)}
    if not sets:
        raise SystemExit('No valid datasets loaded')
    #endif
    engine=GepardEngine()
    first=next(iter(sets.values())).iloc[0]
    preflight=engine.preflight(first)
    (args.output/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n')
    if args.preflight_only:
        print('[done] preflight only; no fits')
        return
    #endif
    summary={}
    for name,df in sets.items():
        # Split energy/torus explicitly; do not pool disparate beam energies.
        for (dataset,E), group in df.groupby(['dataset','E']):
            tag=f'{name}_{dataset}_{E:g}'.replace(' ','_').replace('.','p')
            dest=args.output/tag
            dest.mkdir(parents=True,exist_ok=True)
            print(f'[fit] {tag}: {len(group)} available points',flush=True)
            summary[tag]=fit_dataset(engine,group,args.replicas,args.seed,args.fit_t_slope,
                                     args.max_points,args.t_ref,dest)
            print('[result]',tag,summary[tag],flush=True)
        #endfor
    #endfor
    (args.output/'all_summaries.json').write_text(json.dumps(summary,indent=2)+'\n')
#enddef

if __name__=='__main__':
    main()
#endif
