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
from scipy.stats import chi2 as chi2_dist
import hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    """Read Sebastian e-gamma-p Pass-2 points from Appendix H.

    IMPORTANT: ``published_rga_bin`` is the bin number printed in Appendix H
    and defined by Sebastian's Table 10 (the published 64-bin RGA scheme). It
    is NOT Sebastian's separate native 1--94 bin numbering used in Tables 8/9.

    Tables 8/9 therefore cannot be joined by integer bin number. For the
    Appendix-H comparison the note does not tabulate per-bin mean Q2,xB,t.
    By default we use the supplied, explicit proxy table constructed from the
    published 10.6-GeV Pass-1 mean kinematics in the identical Table-10 cell.
    The provenance is retained in every output row. An independently verified
    Sebastian Appendix-H mean-kinematics table can override this via
    --sebastian-kinematics.
    """
    if q2_lookup is None:
        q2_lookup = HERE / 'import' / 'sebastian_appendixH_kinematics_proxy.csv'
    #endif
    if not Path(q2_lookup).exists():
        raise FileNotFoundError(f'Sebastian Appendix-H kinematics map not found: {q2_lookup}')
    #endif
    a = pd.read_csv(path)
    if 'published_rga_bin' not in a.columns:
        if 'rga_bin' not in a.columns:
            raise ValueError('Sebastian CSV needs published_rga_bin (Table 10)')
        #endif
        a['published_rga_bin'] = a['rga_bin']
    #endif
    if 'bin_numbering' in a.columns:
        bad = a.loc[a.bin_numbering.ne('published_RGA_Table10_64bin')]
        if len(bad):
            raise ValueError('Sebastian CSV contains an unexpected bin-numbering scheme')
        #endif
    #endif
    a = a[(a.series == 'sebastian_pass2') & (a.topology == 'egamma_p')].copy()
    k = pd.read_csv(q2_lookup)
    bincol = 'published_rga_bin' if 'published_rga_bin' in k.columns else 'rga_bin'
    required = {bincol,'Q2','xB','t'}
    if not required.issubset(k.columns) or k[bincol].duplicated().any():
        raise ValueError('Sebastian kinematics must have unique published_rga_bin,Q2,xB,t columns')
    #endif
    if bincol != 'published_rga_bin':
        k = k.rename(columns={bincol:'published_rga_bin'})
    #endif
    keep = ['published_rga_bin','Q2','xB','t']
    for c in ['kinematics_source','note']:
        if c in k.columns:
            keep.append(c)
        #endif
    #endfor
    a = a.merge(k[keep], on='published_rga_bin', how='left', validate='many_to_one')
    missing = a[['Q2','xB','t']].isna().any(axis=1)
    if missing.any():
        bins = sorted(a.loc[missing,'published_rga_bin'].unique().tolist())
        raise ValueError(f'Missing Sebastian Appendix-H kinematics for Table-10 bins: {bins}')
    #endif
    a['dataset'] = 'sebastian_' + a.torus
    a['E'] = 10.6
    a['A'] = a.bsa
    a['err'] = a.bsa_error_graphical
    a['error_type'] = 'PDF_graphical_error'
    if 'kinematics_source' not in a.columns:
        a['kinematics_source'] = 'user_supplied_verified_map'
    #endif
    print('[sebastian] Appendix-H bin IDs interpreted explicitly as published RGA Table-10 bins; '
          f'kinematics source: {a.kinematics_source.iloc[0]}')
    return a[['dataset','E','Q2','xB','t','phi_deg','A','err','error_type',
              'published_rga_bin','kinematics_source']]
#enddef


def read_pass2(path, label):
    """Read native production Pass-2 BSA at its measured mean kinematics.

    This intentionally uses the regular dvcs_pass2_analysis.csv and the exact
    production schema: xBavg/Q2avg/t_abs_avg/phiavg for the requested label.
    No Pass-1 rebinned CSV and no bin-center substitution are allowed.
    """
    if label == '10.6 GeV Unpol':
        raise ValueError(
            "Pass-2 BSA must use label '10.6 GeV'; "
            "'10.6 GeV Unpol' is the unpolarized cross-section group."
        )
    #endif
    df = pd.read_csv(path, low_memory=False)
    bcol = f'BSA, counts, {label}'
    cols = {
        'xB': f'xBavg, {label}', 'Q2': f'Q2avg, {label}',
        'tabs': f't_abs_avg, {label}', 'phi': f'phiavg, {label}',
        'ptp': f'{bcol}, point-to-point sys',
    }
    required = [bcol, *cols.values()]
    print(f"[pass2 schema] file={path}")
    print(f"[pass2 schema] label={label!r}")
    print(f"[pass2 schema] BSA column={bcol!r}")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError('Native Pass-2 CSV is missing required production columns:\n  ' +
                         '\n  '.join(missing))
    #endif
    output = []
    rejected = {'bad_bsa_or_stat':0, 'bad_ptp':0, 'bad_kinematics':0}
    for _, row in df.iterrows():
        A, stat = parse_tuple(row[bcol])
        try:
            ptp=float(row[cols['ptp']]); xb=float(row[cols['xB']]); q2=float(row[cols['Q2']])
            tabs=float(row[cols['tabs']]); phi=float(row[cols['phi']])
        except (TypeError, ValueError):
            rejected['bad_kinematics'] += 1; continue
        #endtry
        if not (np.isfinite(A) and np.isfinite(stat) and stat > 0):
            rejected['bad_bsa_or_stat'] += 1; continue
        #endif
        if not (np.isfinite(ptp) and ptp >= 0):
            rejected['bad_ptp'] += 1; continue
        #endif
        if not all(np.isfinite(v) for v in [xb,q2,tabs,phi]):
            rejected['bad_kinematics'] += 1; continue
        #endif
        output.append(dict(dataset='our_pass2_'+label.replace(' ','_'),
                           E=10.2 if label=='Sp19 Inb' else 10.6,
                           Q2=q2,xB=xb,t=-abs(tabs),phi_deg=phi,A=A,
                           err=math.hypot(stat,ptp),stat_err=stat,error_type='stat_plus_ptp'))
    #endfor
    print(f'[pass2 {label}] accepted {len(output)} native-bin points; rejected {rejected}')
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



def _dataset_fingerprint(data, rgrid, fit_slope, tref):
    cols=['E','Q2','xB','t','phi_deg']
    raw=np.ascontiguousarray(data[cols].to_numpy(float)).tobytes()
    raw += np.ascontiguousarray(rgrid,dtype=float).tobytes()
    raw += f'fit_slope={fit_slope};tref={tref};v=2'.encode()
    return hashlib.sha256(raw).hexdigest()[:16]
#enddef


def build_response_grid(engine, data, rgrid, cache_file, force=False):
    """Precompute A_LU(r_ImH) once; replica fits then require no Gepard calls."""
    fingerprint=_dataset_fingerprint(data,rgrid,False,-0.3)
    if cache_file.exists() and not force:
        z=np.load(cache_file,allow_pickle=False)
        if str(z['fingerprint']) == fingerprint and np.array_equal(z['rgrid'],rgrid):
            response=np.asarray(z['response'],float)
            if response.shape == (len(data),len(rgrid)):
                print(f'[grid] loaded {cache_file} ({response.shape[0]} points x {response.shape[1]} r values)')
                return response
            #endif
        #endif
        print('[grid] cache fingerprint/schema mismatch; rebuilding')
    #endif
    rows=list(data.itertuples(index=False))
    response=np.empty((len(rows),len(rgrid)),float)
    total=len(rows)*len(rgrid)
    done=0
    print(f'[grid] building {len(rows)} x {len(rgrid)} = {total} Gepard predictions once',flush=True)
    for j,r in enumerate(rgrid):
        response[:,j]=[engine.predict(row,float(r),0.,-0.3) for row in rows]
        done += len(rows)
        if j==0 or (j+1)%5==0 or j+1==len(rgrid):
            print(f'  grid {j+1}/{len(rgrid)} ({100*done/total:.1f}%)',flush=True)
        #endif
    #endfor
    cache_file.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(cache_file,fingerprint=np.array(fingerprint),rgrid=rgrid,response=response)
    print(f'[grid] wrote {cache_file}')
    return response
#enddef


def interpolate_response(response,rgrid,r):
    """Vectorized linear interpolation: one common r for all measured points."""
    r=float(r)
    if r < rgrid[0] or r > rgrid[-1]:
        return np.full(response.shape[0],np.nan)
    #endif
    j=int(np.searchsorted(rgrid,r,side='right')-1)
    if j >= len(rgrid)-1:
        return response[:,-1].copy()
    #endif
    f=(r-rgrid[j])/(rgrid[j+1]-rgrid[j])
    return response[:,j]*(1.-f)+response[:,j+1]*f
#enddef

def fit_dataset(engine, data, nrep, seed, fit_slope, max_points, tref, outdir,
                grid_points, force_grid):
    data = data.sort_values(['Q2','xB','t','phi_deg']).copy()
    if max_points and len(data) > max_points:
        data = data.iloc[np.linspace(0,len(data)-1,max_points,dtype=int)].copy()
        print(f'[warning] subsampling to {len(data)} points for pilot runtime')
    #endif
    obs=data.A.to_numpy(float); err=data.err.to_numpy(float)
    if fit_slope:
        raise NotImplementedError('--fit-t-slope is temporarily disabled in the accelerated stage-II path. '
                                  'Validate the one-parameter response-grid study first.')
    #endif
    rgrid=np.linspace(0.1,3.0,int(grid_points))
    cache_file=outdir/'theory_response_grid.npz'
    response=build_response_grid(engine,data,rgrid,cache_file,force_grid)
    if not np.all(np.isfinite(response)):
        bad=np.argwhere(~np.isfinite(response))
        raise RuntimeError(f'Nonfinite theory response grid; first bad index {bad[0].tolist()}')
    #endif
    # Validate interpolation at non-grid r values against direct Gepard before fitting.
    rows=list(data.itertuples(index=False))
    check_rs=[0.735,1.173,1.641]
    check_indices=sorted(set([0,len(rows)//2,len(rows)-1]))
    max_interp=0.
    for r in check_rs:
        pred=interpolate_response(response,rgrid,r)
        for i in check_indices:
            direct=engine.predict(rows[i],r,0.,tref)
            max_interp=max(max_interp,abs(float(pred[i])-float(direct)))
        #endfor
    #endfor
    print(f'[grid] max interpolation validation |delta A_LU| = {max_interp:.3g}')
    if max_interp > 5e-4:
        raise RuntimeError('Theory response grid is too coarse: interpolation validation exceeds 5e-4')
    #endif
    def objective_r(r,target):
        theory=interpolate_response(response,rgrid,r)
        return float(np.sum(((target-theory)/err)**2)) if np.all(np.isfinite(theory)) else 1e15
    #enddef
    def solve(target):
        res=minimize_scalar(lambda r:objective_r(r,target),bounds=(rgrid[0],rgrid[-1]),method='bounded',
                            options={'xatol':1e-5})
        return float(res.x),float(res.fun),bool(res.success)
    #enddef
    best,chi2,ok=solve(obs)
    if not ok: raise RuntimeError('Nominal fit failed')
    rng=np.random.default_rng(seed); reps=[]
    for j in range(nrep):
        rr,value,success=solve(obs+rng.normal(size=len(obs))*err)
        reps.append(dict(replica=j,r_ImH=rr,b_t=0.,chi2=value,valid=success))
        if (j+1)%100==0 or j+1==nrep:
            print(f'  replicas {j+1}/{nrep}',flush=True)
        #endif
    #endfor
    repdf=pd.DataFrame(reps); repdf.to_csv(outdir/'replicas.csv',index=False)
    pred=interpolate_response(response,rgrid,best)
    data['theory']=pred; data['residual']=obs-pred
    data.to_csv(outdir/'point_predictions.csv',index=False)
    good=repdf.query('valid')
    summary=dict(n_points=len(data),r_ImH=best,b_t=0.,chi2=chi2,ndf=len(data)-1,
                 chi2_ndf=chi2/max(1,len(data)-1),
                 replica_r_16=float(good.r_ImH.quantile(.16)) if len(good) else None,
                 replica_r_50=float(good.r_ImH.quantile(.50)) if len(good) else None,
                 replica_r_84=float(good.r_ImH.quantile(.84)) if len(good) else None,
                 replica_valid=len(good),max_interpolation_delta_ALU=max_interp,
                 r_boundary=bool(abs(best-rgrid[0])<.001 or abs(best-rgrid[-1])<.001))
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
#enddef



def exact_tmin(Q2, xB):
    eps2 = 4*M*M*xB*xB/Q2
    return -Q2*(2*(1-xB)*(1-np.sqrt(1+eps2))+eps2)/(4*xB*(1-xB)+eps2)
#enddef


def assign_published_cells(df):
    # Published Pass-1 binning: 16 (Q2,xB) regions x four t-prime regions.
    xedges = [(0,1.4,[0,.13,.21,np.inf]),
              (1.4,1.8,[0,.13,.21,np.inf]),
              (1.8,2.4,[0,.16,.26,np.inf]),
              (2.4,3.25,[0,.21,.33,np.inf]),
              (3.25,5,[0,.33,np.inf]),
              (5,np.inf,[0,.55,np.inf])]
    tp_edges=[0,.2,.4,.8,np.inf]
    out=df.copy().reset_index(drop=True)
    tp=exact_tmin(out.Q2.to_numpy(float),out.xB.to_numpy(float))-out.t.to_numpy(float)
    ids=np.full(len(out),-1,dtype=int)
    cell=0
    for qlo,qhi,xs in xedges:
        for xlo,xhi in zip(xs[:-1],xs[1:]):
            for tlo,thi in zip(tp_edges[:-1],tp_edges[1:]):
                cell+=1
                mask=(out.Q2>=qlo)&(out.Q2<qhi)&(out.xB>=xlo)&(out.xB<xhi)&(tp>=tlo)&(tp<thi)
                ids[mask.to_numpy()]=cell
            #endfor
        #endfor
    #endfor
    out['cell_id']=ids
    out['tprime']=tp
    return out[out.cell_id>0].reset_index(drop=True)
#enddef


def minimize_r(obs,err,response,rgrid):
    def objective(r):
        theory=interpolate_response(response,rgrid,r)
        return float(np.sum(((obs-theory)/err)**2)) if np.all(np.isfinite(theory)) else 1e15
    #enddef
    result=minimize_scalar(objective,bounds=(rgrid[0],rgrid[-1]),method='bounded',options={'xatol':1e-5})
    return float(result.x),float(result.fun),bool(result.success)
#enddef


def reuse_global_grid(engine,selected,rgrid,source,force=False):
    """Load validated full-grid cache and select rows by unique kinematic key.

    No theory evaluations are needed for local fits when full grids exist.
    """
    full,cache_file=source
    full=full.sort_values(['Q2','xB','t','phi_deg']).copy()
    full_response=build_response_grid(engine,full,rgrid,cache_file,force)
    keys=['E','Q2','xB','t','phi_deg']
    # Exact numerical coordinates, including duplicates. Repeated phi-bin
    # means can occur, so match each row to a unique occurrence.
    full_keys=pd.MultiIndex.from_frame(full[keys].reset_index(drop=True))
    if full_keys.has_duplicates:
        # Duplicates at identical kinematics have identical theory responses.
        lookup={k:i for i,k in enumerate(full_keys)}
    else:
        lookup={k:i for i,k in enumerate(full_keys)}
    #endif
    selected_keys=pd.MultiIndex.from_frame(selected[keys].reset_index(drop=True))
    missing=[k for k in selected_keys if k not in lookup]
    if missing:
        raise RuntimeError(f'{len(missing)} local points absent from global cache {cache_file}; '
                           'check source filtering and kinematics before running Gepard.')
    #endif
    indices=np.array([lookup[k] for k in selected_keys],dtype=int)
    print(f'[grid] reused {len(indices)} rows from {cache_file}')
    return full_response[indices]
#enddef


def plot_local_diagnostics(result,p1,p2,g1,g2,grid,outdir,energy):
    """Compact matched-cell visual checks, with no fabricated predictions."""
    if result.empty:
        return
    #endif
    figdir=outdir/'figures'
    figdir.mkdir(parents=True,exist_ok=True)
    good=result[result.fit_valid].sort_values('cell_id')
    if good.empty:
        return
    #endif
    x=good.cell_id.to_numpy()
    fig,axes=plt.subplots(2,1,figsize=(14,8),sharex=True,layout='constrained')
    axes[0].errorbar(x-.12,good.r_pass1,yerr=good.sigma_r_pass1_conditional,
                     fmt='o',label='Published Pass-1')
    axes[0].errorbar(x+.12,good.r_pass2,yerr=good.sigma_r_pass2_conditional,
                     fmt='s',label='Our Pass-2')
    axes[0].set_ylabel(r'Local $r_{\Im\mathcal{H}}$')
    axes[0].legend()
    axes[1].axhline(0,color='0.4',lw=1)
    axes[1].errorbar(x,good.delta_r,
                     yerr=np.hypot(good.sigma_r_pass1_conditional,good.sigma_r_pass2_conditional),
                     fmt='o')
    axes[1].set(xlabel='Published Pass-1 kinematic cell',ylabel=r'$r_{Pass2}-r_{Pass1}$')
    fig.suptitle(f'Pass-1 / Pass-2 local CFF comparison — {energy:g} GeV')
    fig.savefig(figdir/'local_r_comparison.png',dpi=170)
    plt.close(fig)

    fig,axes=plt.subplots(2,1,figsize=(14,8),sharex=True,layout='constrained')
    axes[0].axhline(1,color='0.4',ls='--')
    axes[0].plot(x,good.precision_ratio_sigma_pass2_over_pass1,'o',label='Total Pass-2 / published Pass-1')
    axes[0].plot(x,good.precision_ratio_pass2_stat_over_pass1_published,'s',
                 label='Stat-only Pass-2 / published Pass-1 (not stat/stat)')
    axes[0].set_ylabel('Conditional uncertainty ratio')
    axes[0].legend()
    axes[1].plot(x,good.delta_chi2,'o')
    axes[1].set(xlabel='Published Pass-1 kinematic cell',
                ylabel=r'$\Delta\chi^2$ (shared vs independent)')
    fig.suptitle(f'Precision and compatibility diagnostics — {energy:g} GeV')
    fig.savefig(figdir/'precision_compatibility.png',dpi=170)
    plt.close(fig)

    # Representative panels: highest local incompatibilities, plus best match.
    ids=list(good.sort_values('delta_chi2',ascending=False).cell_id.head(5))
    ids+=list(good.sort_values('delta_chi2').cell_id.head(1))
    ids=list(dict.fromkeys(ids))
    fig,axes=plt.subplots(2,3,figsize=(16,9),sharex=True,sharey=True,layout='constrained')
    for ax,cid in zip(axes.flat,ids):
        row=good[good.cell_id==cid].iloc[0]
        for df,gg,rr,fmt,label in [(p1,g1,row.r_pass1,'o','Pass-1'),
                                   (p2,g2,row.r_pass2,'s','Pass-2')]:
            idx=np.flatnonzero(df.cell_id.to_numpy()==cid)
            if not len(idx):
                continue
            #endif
            phi=df.phi_deg.to_numpy()[idx]
            ax.errorbar(phi,df.A.to_numpy()[idx],yerr=df.err.to_numpy()[idx],
                        fmt=fmt,ms=3,alpha=.75,label=label)
            pred=interpolate_response(gg[idx],grid,rr)
            order=np.argsort(phi)
            ax.plot(phi[order],pred[order],lw=1.4,label=f'{label} own fit')
            shared=interpolate_response(gg[idx],grid,row.r_shared)
            ax.plot(phi[order],shared[order],ls=':',lw=1.1,label=f'{label} shared')
        #endfor
        ax.set_title(f'Cell {cid}: Δχ²={row.delta_chi2:.2f}',fontsize=10)
        ax.set(xlim=(0,360),ylim=(-.6,.6))
        ax.grid(alpha=.15)
    #endfor
    for ax in axes.flat[len(ids):]:
        ax.set_visible(False)
    #endfor
    handles,labels=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=4,fontsize=9)
    fig.supxlabel('φ (degrees)')
    fig.supylabel(r'$A_{LU}$')
    fig.suptitle(f'Representative matched-cell BSA fits — {energy:g} GeV')
    fig.savefig(figdir/'representative_bsa_fits.png',dpi=170)
    plt.close(fig)
#enddef


def matched_local_test(engine,p1,p2,energy,outdir,args,global_grids=None):
    # Restrict BOTH datasets to the common physical phase-space envelope.
    # Pass-1 has substantially wider high-|t|, Q2 and xB reach.
    def overlap_cut(df):
        return df[(df.t.abs() <= args.max_t_abs) &
                  (df.xB <= args.max_xb) & (df.Q2 <= args.max_q2)].copy()
    #enddef
    p1=assign_published_cells(overlap_cut(p1))
    p2=assign_published_cells(overlap_cut(p2))
    print(f'[matched {energy:g}] overlap envelope: Pass-1={len(p1)}, Pass-2={len(p2)} '
          f'(|t|<={args.max_t_abs}, xB<={args.max_xb}, Q2<={args.max_q2})')
    p1=p1.sort_values(['cell_id','Q2','xB','t','phi_deg']).reset_index(drop=True)
    p2=p2.sort_values(['cell_id','Q2','xB','t','phi_deg']).reset_index(drop=True)
    if p1.empty or p2.empty:
        print(f'[matched {energy:g}] no valid points, skipping')
        return None
    #endif
    outdir.mkdir(parents=True,exist_ok=True)
    grid=np.linspace(.1,3.,args.grid_points)
    g1=reuse_global_grid(engine,p1,grid,global_grids['pass1'][energy],args.force_grid)
    g2=reuse_global_grid(engine,p2,grid,global_grids['pass2'][energy],args.force_grid)
    if not (np.isfinite(g1).all() and np.isfinite(g2).all()):
        raise ValueError('Nonfinite Gepard response in matched-cell study')
    #endif
    ids=sorted(set(p1.cell_id).intersection(set(p2.cell_id)))
    rows=[]
    for cid in ids:
        i1=np.flatnonzero(p1.cell_id.to_numpy()==cid)
        i2=np.flatnonzero(p2.cell_id.to_numpy()==cid)
        if min(len(i1),len(i2))<args.min_points_per_dataset:
            continue
        #endif
        a1=p1.A.to_numpy(float)[i1]; e1=p1.err.to_numpy(float)[i1]
        a2=p2.A.to_numpy(float)[i2]; e2=p2.err.to_numpy(float)[i2]
        r1,c1,ok1=minimize_r(a1,e1,g1[i1],grid)
        r2,c2,ok2=minimize_r(a2,e2,g2[i2],grid)
        both_a=np.concatenate([a1,a2]); both_e=np.concatenate([e1,e2])
        both_g=np.vstack([g1[i1],g2[i2]])
        rs,cs,oks=minimize_r(both_a,both_e,both_g,grid)
        dchi=max(0.,cs-c1-c2)
        # Cross-prediction: evaluate each independently fitted response at the
        # OTHER dataset's actual measured mean kinematics and phi coordinates.
        pass1_predicts_pass2=interpolate_response(g2[i2],grid,r1)
        pass2_predicts_pass1=interpolate_response(g1[i1],grid,r2)
        chi2_cross_1_to_2=float(np.sum(((a2-pass1_predicts_pass2)/e2)**2))
        chi2_cross_2_to_1=float(np.sum(((a1-pass2_predicts_pass1)/e1)**2))
        # Local curvature estimates (conditional, diagonal uncertainties).
        def local_sigma(response,err,r):
            h=0.002
            if r-h<=grid[0] or r+h>=grid[-1]: return np.nan
            deriv=(interpolate_response(response,grid,r+h)-interpolate_response(response,grid,r-h))/(2*h)
            fisher=float(np.sum((deriv/err)**2))
            return 1/np.sqrt(fisher) if fisher>0 else np.nan
        #enddef
        sig1=local_sigma(g1[i1],e1,r1)
        sig2=local_sigma(g2[i2],e2,r2)
        # Statistical-only Pass-2 precision uses the original statistical error.
        stat2=p2.stat_err.to_numpy(float)[i2] if 'stat_err' in p2 else e2
        sig2_stat=local_sigma(g2[i2],stat2,r2)
        boundary=any(abs(r-grid[0])<.001 or abs(r-grid[-1])<.001 for r in [r1,r2,rs])
        rows.append(dict(energy_GeV=energy,cell_id=cid,n_pass1=len(i1),n_pass2=len(i2),
                         r_pass1=r1,r_pass2=r2,r_shared=rs,delta_r=r2-r1,
                         sigma_r_pass1_conditional=sig1,sigma_r_pass2_conditional=sig2,
                         precision_ratio_sigma_pass2_over_pass1=sig2/sig1 if sig1>0 else np.nan,
                         chi2_pass1=c1,chi2_pass2=c2,chi2_shared=cs,
                         chi2_pass1_predicts_pass2=chi2_cross_1_to_2,
                         chi2_pass2_predicts_pass1=chi2_cross_2_to_1,
                         cross_prediction_delta_chi2_pass2=chi2_cross_1_to_2-c2,
                         cross_prediction_delta_chi2_pass1=chi2_cross_2_to_1-c1,
                         delta_chi2=dchi,compatibility_p_value=chi2_dist.sf(dchi,1),
                         fit_valid=ok1 and ok2 and oks and not boundary,
                         r_boundary=boundary,
                         mean_Q2_pass1=float(p1.Q2.iloc[i1].mean()),mean_Q2_pass2=float(p2.Q2.iloc[i2].mean()),
                         mean_xB_pass1=float(p1.xB.iloc[i1].mean()),mean_xB_pass2=float(p2.xB.iloc[i2].mean()),
                         mean_tprime_pass1=float(p1.tprime.iloc[i1].mean()),
                         mean_tprime_pass2=float(p2.tprime.iloc[i2].mean())))
    #endfor
    result=pd.DataFrame(rows)
    result.to_csv(outdir/'matched_local_cells.csv',index=False)
    plot_local_diagnostics(result,p1,p2,g1,g2,grid,outdir,energy)
    valid=result[result.fit_valid] if len(result) else result
    if len(valid):
        dchi=float(valid.delta_chi2.sum())
        summary=dict(energy_GeV=energy,n_valid_cells=len(valid),n_cells=len(result),
                     delta_chi2=dchi,delta_ndf=len(valid),global_compatibility_p_value=float(chi2_dist.sf(dchi,len(valid))),
                     median_delta_r=float(valid.delta_r.median()),
                     median_conditional_precision_ratio_pass2_over_pass1=float(valid.precision_ratio_sigma_pass2_over_pass1.median()),
                     total_chi2_pass1_predicts_pass2=float(valid.chi2_pass1_predicts_pass2.sum()),
                     total_chi2_pass2_predicts_pass1=float(valid.chi2_pass2_predicts_pass1.sum()))
        (outdir/'matched_local_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(f'[matched {energy:g}] {len(valid)} valid cells; delta chi2={dchi:.3f} '
              f'for {len(valid)} dof; p={summary["global_compatibility_p_value"]:.4g}')
        return summary
    #endif
    print(f'[matched {energy:g}] no valid matched cells')
    return None
#enddef


def run_matched(engine,args,global_grids):
    p1=read_pass1(args.pass1)
    # Use native Sp19 Inb 10.2-GeV BSA and combined 10.6-GeV BSA.
    jobs=[(10.6,'10.6 GeV'),(10.2,args.pass2_10p2_label)]
    summaries=[]
    for energy,label in jobs:
        if not label:
            print(f'[matched {energy:g}] skipped: no native Pass-2 BSA label supplied')
            continue
        #endif
        try:
            p2=read_pass2(args.pass2,label)
        except ValueError as exc:
            raise ValueError(f'[matched {energy:g}] required native Pass-2 BSA unavailable: {exc}') from exc
        #endtry
        q=p1[np.isclose(p1.E,energy,atol=.05)].copy()
        if q.empty:
            print(f'[matched {energy:g}] no Pass-1 points at this energy')
            continue
        #endif
        dest=args.output/'matched_local'/f'{energy:g}_GeV'.replace('.','p')
        summary=matched_local_test(engine,q,p2,energy,dest,args,global_grids)
        if summary is not None:
            summaries.append(summary)
        #endif
    #endfor
    if summaries:
        pd.DataFrame(summaries).to_csv(args.output/'matched_local'/'summary.csv',index=False)
    #endif
#enddef

def plot_global_summary(summary,outdir):
    rows=[]
    for tag,value in summary.items():
        rows.append((tag,value['r_ImH'],value.get('replica_r_16'),value.get('replica_r_84')))
    #endfor
    if not rows:
        return
    #endif
    fig,ax=plt.subplots(figsize=(11,5),layout='constrained')
    for i,(tag,r,lo,hi) in enumerate(rows):
        ax.errorbar(r,i,xerr=[[max(0,r-lo)],[max(0,hi-r)]] if lo is not None and hi is not None else None,
                    fmt='o',capsize=3)
    #endfor
    ax.set_yticks(range(len(rows)),[v[0] for v in rows])
    ax.invert_yaxis()
    ax.axvline(1,ls='--',color='0.5')
    ax.set_xlabel(r'Global $r_{\Im\mathcal{H}}$ (KM15-relative; differing kinematic coverage)')
    figdir=outdir/'figures'
    figdir.mkdir(parents=True,exist_ok=True)
    fig.savefig(figdir/'global_r_comparison.png',dpi=170)
    plt.close(fig)
#enddef


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pass1',type=Path,default=HERE/'import/rga_pass1_bsa.txt')
    ap.add_argument('--sebastian',type=Path,default=HERE/'import/sebastian_pass2_bsa_from_note.csv')
    ap.add_argument('--sebastian-kinematics',type=Path,default=None,
                    help='Optional verified Appendix-H means CSV (published_rga_bin,Q2,xB,t). Default uses documented Pass-1 10.6 mean proxy in the identical Table-10 cell.')
    ap.add_argument('--pass2',type=Path,default=HERE.parent/'output/csvs/dvcs_pass2_analysis.csv')
    ap.add_argument('--pass2-label',default='10.6 GeV',help='BSA column label in native Pass2 CSV')
    ap.add_argument('--output',type=Path,default=HERE.parent/'output/bsa/cff_crosscheck')
    ap.add_argument('--datasets',default='all',help='Comma-separated: pass1,sebastian,pass2,pass2_10p2')
    ap.add_argument('--replicas',type=int,default=20)
    ap.add_argument('--max-points',type=int,default=0,
                    help='Pilot limit per dataset (0=all); deterministic subsample')
    ap.add_argument('--grid-points',type=int,default=81,
                    help='Number of precomputed r_ImH values spanning 0.1--3.0')
    ap.add_argument('--force-grid',action='store_true',help='Rebuild saved Gepard response grid')
    ap.add_argument('--seed',type=int,default=20261012)
    ap.add_argument('--fit-t-slope',action='store_true',help='r(t)=r0*exp[b*(t-tref)]')
    ap.add_argument('--t-ref',type=float,default=-0.3)
    ap.add_argument('--preflight-only',action='store_true')
    ap.add_argument('--matched-compatibility',action='store_true',help='Local shared-vs-independent Pass-1/Pass-2 ImH test')
    ap.add_argument('--pass2-10p2-label',default='Sp19 Inb',help='Native Pass-2 10.2-GeV BSA label (default: Sp19 Inb)')
    ap.add_argument('--min-points-per-dataset',type=int,default=4)
    ap.add_argument('--max-t-abs',type=float,default=0.92,
                    help='Common-overlap maximum |t| (GeV^2)')
    ap.add_argument('--max-xb',type=float,default=0.53,
                    help='Common-overlap maximum xB')
    ap.add_argument('--max-q2',type=float,default=5.1,
                    help='Common-overlap maximum Q2 (GeV^2)')
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    sets={}
    requested={x.strip() for x in args.datasets.split(',')}
    if args.matched_compatibility:
        requested.add('all')
    #endif
    if 'all' in requested:
        requested.update(['pass1','pass2','pass2_10p2'])
    #endif
    if 'pass1' in requested:
        sets['pass1']=read_pass1(args.pass1)
    #endif
    if 'sebastian' in requested:
        sets['sebastian']=read_sebastian(args.sebastian,args.sebastian_kinematics)
    #endif
    if 'pass2' in requested:
        sets['pass2']=read_pass2(args.pass2,args.pass2_label)
    #endif
    if 'pass2_10p2' in requested:
        sets['pass2_10p2']=read_pass2(args.pass2,args.pass2_10p2_label)
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
                                     args.max_points,args.t_ref,dest,args.grid_points,args.force_grid)
            print('[result]',tag,summary[tag],flush=True)
        #endfor
    #endfor
    (args.output/'all_summaries.json').write_text(json.dumps(summary,indent=2)+'\n')
    if 'all' in requested or args.matched_compatibility:
        if args.max_points:
            raise ValueError('Matched compatibility requires full global grids; omit --max-points')
        #endif
        grid_sources={'pass1':{},'pass2':{}}
        for (dataset,E),group in sets['pass1'].groupby(['dataset','E']):
            tag=f'pass1_{dataset}_{E:g}'.replace(' ','_').replace('.','p')
            grid_sources['pass1'][round(float(E),1)]=(group,args.output/tag/'theory_response_grid.npz')
        #endfor
        for key in ['pass2','pass2_10p2']:
            for (dataset,E),group in sets[key].groupby(['dataset','E']):
                tag=f'{key}_{dataset}_{E:g}'.replace(' ','_').replace('.','p')
                grid_sources['pass2'][round(float(E),1)]=(group,args.output/tag/'theory_response_grid.npz')
            #endfor
        #endfor
        run_matched(engine,args,grid_sources)
    #endif
    plot_global_summary(summary,args.output)

#enddef

if __name__=='__main__':
    main()
#endif
