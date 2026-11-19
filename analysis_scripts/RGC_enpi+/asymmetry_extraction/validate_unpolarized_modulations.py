#!/usr/bin/env python3
"""Validate UU modulations with production-likelihood profiles and conditional toys.

From RGC_enpi+: python3 asymmetry_extraction/validate_unpolarized_modulations.py
No nominal products are written. Plot output is PNG only. See UU_VALIDATION.md.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import hashlib
import json
from importlib.metadata import version
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

# Avoid multiplying BLAS threads inside each worker.
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import extract_structure_function_ratios as prod

HERE = Path(__file__).resolve().parent
UU = ('u1', 'u2')
SCENARIOS = ('nominal', 'zero_uu', 'opposite_uu', 'half_polarized', 'zero_polarized')
_WORK = None


def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(prod.json_safe(payload), indent=2, allow_nan=False))
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def usable(fit):
    return bool(fit['valid'] and fit['accurate_covariance']
                and fit['positive_definite_covariance'] and not fit['parameters_at_limit']
                and math.isfinite(fit['minimum_nll'])
                and all(math.isfinite(fit['errors'][p]) and fit['errors'][p] > 0 for p in UU))


def compact_fit(fit):
    return {k: fit[k] for k in ('values', 'errors', 'minimum_nll', 'valid',
            'accurate_covariance', 'positive_definite_covariance', 'parameters_at_limit',
            'edm', 'nfcn', 'parameter_order', 'covariance', 'correlation')}


def cross_factor(events, state, period, truth, dilution, indices=None, h=1, pt=0):
    e = events if indices is None else {k: v[indices] for k, v in events.items()}
    return prod.evaluate_cross_section_factor(
        variant='nominal', phi=e['phi'], r_b=e['rB'], r_c=e['rC'],
        r_v=e['rV'], r_w=e['rW'], sin_theta_gamma=e['sin_theta_gamma'],
        cos_theta_gamma=e['cos_theta_gamma'], helicity=h,
        beam_polarization=prod.BEAM_POLARIZATION[period], target_polarization=pt,
        dilution=dilution, transverse_scales=None,
        **{p: float(truth[p]) for p in prod.PHYSICS_PARAMETERS})


def generator_envelopes(events, run_states, truth, true_dilutions):
    """Check every possible run/state through the extrema of its linear Pt term.

    h is +/-1 and, at fixed h, g is linear in Pt. Its extrema occur at the
    smallest/largest run Pt. The exact pointwise maximum is a rejection envelope.
    Reject nonphysical truth explicitly; never clip g or silently rescale truth.
    """
    envelopes = {}
    for period in prod.PERIODS:
        idx = np.flatnonzero(events['period_index'] == prod.PERIOD_INDEX[period])
        if idx.size == 0:
            continue
        state = run_states[period]
        q = np.ravel(np.column_stack((state['q_plus'], state['q_minus'])))
        if np.any(~np.isfinite(q)) or np.any(q < 0) or q.sum() <= 0:
            raise ValueError(f'Invalid run charges for {period}')
        active = (state['q_plus'] + state['q_minus']) > 0
        pts = state['pt'][active]
        e = {k: v[idx] for k, v in events.items()}
        corners = np.stack([cross_factor(e, state, period, truth, true_dilutions[period],
                                         h=h, pt=pt)
                            for h in (-1, 1) for pt in (float(pts.min()), float(pts.max()))])
        if not np.isfinite(corners).all() or np.any(corners <= prod.CROSS_SECTION_FLOOR):
            raise ValueError(f'Injected cross section is nonpositive in {period}')
        envelopes[period] = (idx, np.max(corners, axis=0))
    return envelopes


def generate_states(events, run_states, truth, true_dilutions, rng, envelopes=None):
    """Draw P(run,h|x,period) proportional to Q(run,h)*g using rejection sampling.

    Kinematics, event counts and period assignments stay fixed; observed run and
    helicity are regenerated jointly. Sampling helicities 50/50 would be wrong.
    """
    if envelopes is None:
        envelopes = generator_envelopes(events, run_states, truth, true_dilutions)
    result = dict(events)
    result['runnum'] = np.array(events['runnum'], copy=True)
    result['helicity'] = np.array(events['helicity'], copy=True)
    for period, (idx, maximum) in envelopes.items():
        state = run_states[period]
        q = np.ravel(np.column_stack((state['q_plus'], state['q_minus'])))
        cumulative = np.cumsum(q / q.sum())
        cumulative[-1] = 1.0
        remaining = np.arange(idx.size)
        for _ in range(10000):
            if not remaining.size:
                break
            proposed = np.searchsorted(cumulative, rng.random(remaining.size), side='right')
            run_index = proposed // 2
            h = np.where(proposed % 2 == 0, 1, -1)
            selected = idx[remaining]
            g = cross_factor(events, state, period, truth, true_dilutions[period],
                             indices=selected, h=h, pt=state['pt'][run_index])
            ratio = g / maximum[remaining]
            if np.any(ratio > 1 + 1e-12) or np.any(ratio <= 0):
                raise RuntimeError('Invalid rejection envelope or state probability')
            accepted = rng.random(remaining.size) < ratio
            result['runnum'][selected[accepted]] = state['run'][run_index[accepted]]
            result['helicity'][selected[accepted]] = h[accepted]
            remaining = remaining[~accepted]
        else:
            raise RuntimeError('Conditional state sampler did not converge')
    return result


def scenario_truth(nominal, scenario):
    values = {p: float(nominal['values'][p]) for p in prod.PHYSICS_PARAMETERS}
    if scenario == 'zero_uu':
        values.update(u1=0., u2=0.)
    elif scenario == 'opposite_uu':
        values.update(u1=-values['u1'], u2=-values['u2'])
    elif scenario in ('half_polarized', 'zero_polarized'):
        scale = .5 if scenario == 'half_polarized' else 0.
        for p in prod.POLARIZED_PARAMETERS:
            values[p] *= scale
    elif scenario != 'nominal':
        raise ValueError(scenario)
    return values


def auxiliary_records(records, bin_number, true_dilutions, rng, mode):
    """Repeat auxiliary measurements at fixed true dilution, without truncation.

    fixed_truth is a conditional diagnostic, not an unconditional coverage test.
    With fluctuate, inadmissible draws are recorded as failures, never redrawn.
    """
    result = dict(records)
    centers = {}
    for period in prod.PERIODS:
        record = records[(period, bin_number)]
        center = true_dilutions[period]
        if mode == 'fluctuate' and record.stat_uncertainty > 0:
            center += float(rng.normal(0., record.stat_uncertainty))
        if center <= 0 or not math.isfinite(center):
            raise ValueError('Nonpositive auxiliary dilution draw')
        result[(period, bin_number)] = replace(record, value=center)
        centers[period] = center
    return result, centers


def constrained_fit(events, states, records, b, nominal, fixed, warm=None):
    """Two independent seeds, always reoptimizing all unfixed physics/nuisances."""
    candidates = []
    for start in (warm or nominal['values'], nominal['values']):
        fit = prod.fit_one_variant(events, states, records, b, 'nominal',
                                  initial_values=start, fixed_physics_parameters=fixed,
                                  fast_profile=True)
        if fit['valid'] and math.isfinite(fit['minimum_nll']) and fit['minimum_nll'] < prod.INVALID_NLL / 10:
            candidates.append(fit)
    if not candidates:
        return None
    best = min(candidates, key=lambda f: f['minimum_nll'])
    if best['minimum_nll'] < nominal['minimum_nll'] - 1e-3:
        raise RuntimeError('Constrained scan found a lower minimum than the nominal fit; inspect bin checkpoint')
    return best


def profile_grid(center, sigma, limits, points, extent):
    grid = np.clip(center + np.linspace(-extent, extent, points) * sigma, *limits)
    return np.unique(np.append(grid, center))


def profile_crossing(rows, center, side, level):
    """Interpolate only across adjacent successful points, never across fit gaps."""
    selected = sorted([r for r in rows if (r['value'] - center) * side >= -1e-12],
                      key=lambda r: abs(r['value'] - center))
    previous = None
    for row in selected:
        if not row['valid']:
            previous = None
            continue
        if previous is not None and previous['two_delta_nll'] <= level <= row['two_delta_nll']:
            a, z = previous, row
            if z['two_delta_nll'] == a['two_delta_nll']:
                return None
            return float(a['value'] + (level - a['two_delta_nll']) *
                         (z['value'] - a['value']) / (z['two_delta_nll'] - a['two_delta_nll']))
        previous = row
    return None


def one_profile(events, states, records, b, nominal, parameter, points, extent):
    center, sigma = nominal['values'][parameter], nominal['errors'][parameter]
    grid = profile_grid(center, sigma, prod.PARAMETER_LIMITS[parameter], points, extent)
    rows = [{'value': float(center), 'two_delta_nll': 0., 'valid': True}]
    for side in (-1, 1):
        warm = nominal['values']
        for x in sorted([x for x in grid if (x-center)*side > 1e-12], key=lambda x: abs(x-center)):
            fit = constrained_fit(events, states, records, b, nominal, {parameter: float(x)}, warm)
            rows.append({'value': float(x), 'valid': fit is not None,
                         'two_delta_nll': max(0., 2*(fit['minimum_nll']-nominal['minimum_nll'])) if fit else None})
            if fit:
                warm = fit['values']
    rows.sort(key=lambda r: r['value'])
    intervals = {str(level): [profile_crossing(rows, center, side, level) for side in (-1, 1)]
                 for level in (1., 3.84)}
    return {'points': rows, 'intervals': intervals, 'hesse_sigma': sigma,
            'center': center, 'limits': list(prod.PARAMETER_LIMITS[parameter])}


def null_identifiability(events, states, records, b, truth, true_dilutions):
    """At zero polarization all UU factors cancel: the NLL must be flat."""
    values = dict(truth)
    values.update({p: 0. for p in prod.POLARIZED_PARAMETERS})
    nll, _ = prod.make_bin_nll(events, states, records, b, 'nominal')
    nuisance = {f'f_{p}': true_dilutions[p] for p in prod.PERIODS}
    rows = []
    for u1, u2 in ((0., 0.), (.1, 0.), (-.1, 0.), (0., .1), (0., -.1)):
        values.update(u1=u1, u2=u2)
        generator_envelopes(events, states, values, true_dilutions)
        rows.append({'u1': u1, 'u2': u2, 'nll': float(nll(**values, **nuisance))})
    spread = max(r['nll'] for r in rows)-min(r['nll'] for r in rows)
    return {'points': rows, 'nll_spread': spread,
            'flat_within_numerical_precision': spread < 1e-7 * max(1, len(events['phi'])/10000)}


def init_worker(cache, states, records, config, output):
    global _WORK
    _WORK = (prod.load_event_cache(Path(cache)), states, records, config, Path(output))


def run_bin(b):
    events_all, states, records, config, output = _WORK
    events = {k: v[events_all['bin_number'] == b] for k, v in events_all.items()}
    path = output/'checkpoints'/f'bin_{b:02d}.json'
    state = json.loads(path.read_text()) if path.exists() else {'bin_number': b}
    if state.get('complete'):
        print(f'[bin {b}] reusing matching completed checkpoint', flush=True)
        return state
    if 'nominal' not in state:
        print(f'[bin {b}] fitting production nominal minimum', flush=True)
        state['nominal'] = prod.fit_one_variant(events, states, records, b, 'nominal')
        atomic_json(path, state)
    nominal = state['nominal']
    if not usable(nominal):
        state.update(complete=True, status='unusable_nominal')
        atomic_json(path, state)
        return state
    true_dilutions = {p: float(nominal['values'][f'f_{p}']) for p in prod.PERIODS}
    if 'profiles' in config['studies']:
        profiles = state.setdefault('profiles', {})
        for par in UU:
            if par not in profiles:
                print(f'[bin {b}] profiling {par}', flush=True)
                profiles[par] = one_profile(events, states, records, b, nominal, par,
                                           config['profile_points'], config['extent'])
                atomic_json(path, state)
        joint = state.setdefault('joint', {'rows': []})
        grids = [profile_grid(nominal['values'][p], nominal['errors'][p],
                              prod.PARAMETER_LIMITS[p], config['joint_points'], config['extent']) for p in UU]
        for i in range(len(joint['rows']), len(grids[0])):
            row, warm = [], nominal['values']
            for y in grids[1]:
                fixed = {'u1': float(grids[0][i]), 'u2': float(y)}
                fit = constrained_fit(events, states, records, b, nominal, fixed, warm)
                row.append({**fixed, 'valid': fit is not None,
                            'two_delta_nll': max(0., 2*(fit['minimum_nll']-nominal['minimum_nll'])) if fit else None})
                if fit:
                    warm = fit['values']
            joint['rows'].append(row)
            atomic_json(path, state)
            if (i+1) % 5 == 0 or i+1 == len(grids[0]):
                print(f'[bin {b}] joint contour rows {i+1}/{len(grids[0])}', flush=True)
    if 'closure' in config['studies']:
        scenarios = state.setdefault('scenarios', {})
        for name in config['scenarios']:
            scenario = scenarios.setdefault(name, {'rows': []})
            print(f'[bin {b}] injection scenario {name}', flush=True)
            truth = scenario_truth(nominal, name)
            scenario.update(truth=truth, true_dilutions=true_dilutions)
            try:
                envelopes = generator_envelopes(events, states, truth, true_dilutions)
            except ValueError as exc:
                scenario.update(status='nonphysical_injection', reason=str(exc))
                atomic_json(path, state)
                if name == 'nominal':
                    state.update(complete=True, status='nonphysical_nominal')
                    atomic_json(path, state)
                    return state
                continue
            if name == 'zero_polarized':
                scenario['identifiability'] = null_identifiability(events, states, records, b, truth, true_dilutions)
            for rep in range(len(scenario['rows']), config['replicas']):
                rng = np.random.default_rng(np.random.SeedSequence([config['seed'], b, SCENARIOS.index(name), rep]))
                toy = generate_states(events, states, truth, true_dilutions, rng, envelopes)
                row = {'replica': rep, 'seed_components': [config['seed'], b, SCENARIOS.index(name), rep]}
                try:
                    aux, centers = auxiliary_records(records, b, true_dilutions, rng, config['auxiliary_mode'])
                    row['auxiliary_centers'] = centers
                    if name == 'zero_polarized':
                        # The true UU directions have zero information. Gaussian
                        # covariance/coverage has no regular interpretation here.
                        flat = null_identifiability(toy, states, aux, b, truth, true_dilutions)
                        row.update(status='structurally_unidentified', usable=False,
                                   identifiability=flat)
                    else:
                        # Do not seed fits with injected truth: use the same nominal seed for all scenarios.
                        fit = prod.fit_one_variant(toy, states, aux, b, 'nominal', initial_values=nominal['values'])
                        row.update(status='fit', usable=usable(fit), fit=compact_fit(fit))
                except (ValueError, RuntimeError) as exc:
                    row.update(status='failed', usable=False, reason=str(exc))
                scenario['rows'].append(row)
                if (rep+1) % 10 == 0:
                    atomic_json(path, state)
                    print(f'[bin {b}] {name}: {rep+1}/{config["replicas"]} toys', flush=True)
            scenario['status'] = 'complete'
            atomic_json(path, state)
    state.update(complete=True, status='complete')
    atomic_json(path, state)
    return state


def summaries(results, config):
    profiles, toys, pulls, flat = [], [], [], []
    for result in results:
        b = result['bin_number']
        for p, scan in result.get('profiles', {}).items():
            lo, hi = scan['intervals']['1.0']
            profiles.append({'bin_number': b, 'parameter': p, 'center': scan['center'],
                'hesse_sigma': scan['hesse_sigma'], 'profile_low68': lo, 'profile_high68': hi,
                'profile_low95': scan['intervals']['3.84'][0], 'profile_high95': scan['intervals']['3.84'][1],
                'profile_over_hesse_low': (scan['center']-lo)/scan['hesse_sigma'] if lo is not None else None,
                'profile_over_hesse_high': (hi-scan['center'])/scan['hesse_sigma'] if hi is not None else None,
                'invalid_scan_points': sum(not r['valid'] for r in scan['points'])})
        for name, scenario in result.get('scenarios', {}).items():
            if 'identifiability' in scenario:
                flat.append({'bin_number': b, **scenario['identifiability']})
            rows = scenario['rows']
            for p in UU:
                good = [r for r in rows if r.get('usable')]
                estimates = np.array([r['fit']['values'][p] for r in good])
                errors = np.array([r['fit']['errors'][p] for r in good])
                truth = scenario['truth'][p]
                delta = estimates-truth
                identifiable = name != 'zero_polarized'
                z = delta/errors if identifiable else np.array([])
                n = len(rows)
                cover = {s: int(np.sum(np.abs(delta) <= s*errors)) for s in (1., 1.96)}
                toys.append({'bin_number': b, 'scenario': name, 'parameter': p,
                    'scenario_status': scenario.get('status'), 'auxiliary_mode': config['auxiliary_mode'],
                    'identifiable_by_construction': identifiable, 'truth': truth,
                    'requested': config['replicas'], 'attempted': n, 'usable': len(good),
                    'failed_or_unusable': sum(not r.get('usable') and r['status']!='structurally_unidentified' for r in rows),
                    'structurally_unidentified': sum(r['status']=='structurally_unidentified' for r in rows),
                    'bias': float(np.mean(delta)) if len(good) and identifiable else None,
                    'empirical_spread': float(np.std(estimates, ddof=1)) if len(good)>1 and identifiable else None,
                    'median_reported_sigma': float(np.median(errors)) if len(good) and identifiable else None,
                    'pull_mean': float(np.mean(z)) if z.size else None,
                    'pull_rms': float(np.sqrt(np.mean(z*z))) if z.size else None,
                    'coverage68_usable': cover[1.]/len(good) if len(good) and identifiable else None,
                    'coverage95_usable': cover[1.96]/len(good) if len(good) and identifiable else None,
                    'coverage68_all_attempts': cover[1.]/n if n and identifiable else None,
                    'coverage95_all_attempts': cover[1.96]/n if n and identifiable else None,
                    'coverage68_mc_se': math.sqrt(.682689*(1-.682689)/len(good)) if len(good) and identifiable else None,
                    'coverage95_mc_se': math.sqrt(.95*.05/len(good)) if len(good) and identifiable else None})
                for row in rows:
                    fit = row.get('fit', {})
                    val, err = fit.get('values', {}).get(p), fit.get('errors', {}).get(p)
                    pulls.append({'bin_number': b, 'scenario': name, 'replica': row['replica'], 'parameter': p,
                        'usable': row['usable'], 'status': row['status'], 'truth': truth, 'estimate': val, 'sigma': err,
                        'pull': (val-truth)/err if row['usable'] and identifiable else None,
                        'at_limit': fit.get('parameters_at_limit')})
    baseline_sigma = {(r['bin_number'],r['parameter']):r['median_reported_sigma']
                      for r in toys if r['scenario']=='nominal'}
    for row in toys:
        reference = baseline_sigma.get((row['bin_number'],row['parameter']))
        row['sigma_ratio_to_nominal'] = (row['median_reported_sigma']/reference
            if reference and row['median_reported_sigma'] is not None else None)
    return profiles, toys, pulls, flat


def write_products(results, config, output):
    table = output/'tables'
    plots = output/'plots'
    table.mkdir(exist_ok=True); plots.mkdir(exist_ok=True)
    profiles, toys, pulls, flat = summaries(results, config)
    nominal_rows, correlation_rows = [], []
    for result in results:
        fit = result.get('nominal')
        if not fit:
            continue
        nominal_rows.append({'bin_number':result['bin_number'], 'status':result.get('status'),
            'number_of_events':fit['metadata']['number_of_events'], 'usable':usable(fit),
            'edm':fit['edm'], 'at_limit':fit['parameters_at_limit'],
            **{p:fit['values'][p] for p in prod.PHYSICS_PARAMETERS},
            **{p+'_stat':fit['errors'][p] for p in prod.PHYSICS_PARAMETERS}})
        if fit['correlation'] is not None:
            order = fit['parameter_order']
            for par in UU:
                for j,other in enumerate(order):
                    if other != par:
                        correlation_rows.append({'bin_number':result['bin_number'],
                            'parameter':par, 'other_parameter':other,
                            'correlation':fit['correlation'][order.index(par)][j]})
    pd.DataFrame(nominal_rows).to_csv(table/'nominal_fits.csv',index=False)
    pd.DataFrame(correlation_rows).to_csv(table/'unpolarized_correlations.csv',index=False)
    for name, rows in (('profile_intervals', profiles), ('closure_summary', toys), ('toy_results', pulls)):
        pd.DataFrame(rows).to_csv(table/f'{name}.csv', index=False)
    atomic_json(table/'zero_polarization_identifiability.json', flat)
    for result in results:
        b = result['bin_number']
        if result.get('profiles'):
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            for ax, p in zip(axes, UU):
                scan = result['profiles'][p]
                ax.plot([r['value'] for r in scan['points']],
                        [r['two_delta_nll'] if r['valid'] else np.nan for r in scan['points']], 'ko-')
                grid = np.array([r['value'] for r in scan['points']])
                ax.plot(grid, ((grid-scan['center'])/scan['hesse_sigma'])**2, '--', color='.5', label='Hesse parabola')
                ax.axhline(1, color='tab:blue', label='68.3% asymptotic (1 parameter)')
                ax.axhline(3.84, color='tab:orange', label='95% asymptotic (1 parameter)')
                ax.set(xlabel=prod.PARAMETER_LABELS[p], ylabel=r'$2\Delta\mathrm{NLL}$', ylim=(0, 8))
                ax.legend(fontsize=8)
            fig.suptitle(f'Bin {b}: profiled UU likelihoods'); fig.tight_layout()
            fig.savefig(plots/f'profiles_bin_{b:02d}.png', dpi=170); plt.close(fig)
            rows = result['joint']['rows']
            x = np.array([r[0]['u1'] for r in rows]); y = np.array([r['u2'] for r in rows[0]])
            z = np.array([[r['two_delta_nll'] if r['valid'] else np.nan for r in row] for row in rows])
            fig, ax = plt.subplots(figsize=(6, 5))
            if np.isfinite(z).any() and min(z.shape) >= 2:
                contour = ax.contour(x, y, z.T, levels=[2.30, 5.99], colors=['tab:blue', 'tab:orange'])
                ax.clabel(contour, fmt={2.30:'68.3%', 5.99:'95%'})
            ax.plot(result['nominal']['values']['u1'], result['nominal']['values']['u2'], 'k+')
            ax.set(xlabel=prod.PARAMETER_LABELS['u1'], ylabel=prod.PARAMETER_LABELS['u2'],
                   title=f'Bin {b}: joint profile, 2 parameters')
            fig.tight_layout(); fig.savefig(plots/f'joint_bin_{b:02d}.png', dpi=170); plt.close(fig)
    if toys:
        frame = pd.DataFrame(toys)
        for p in UU:
            fig, axes = plt.subplots(2, 2, figsize=(11, 8))
            for scenario in config['scenarios']:
                if scenario == 'zero_polarized':
                    continue
                part = frame[(frame.parameter==p) & (frame.scenario==scenario)]
                for ax, key in zip(axes.ravel(), ('bias', 'pull_rms', 'coverage68_usable', 'median_reported_sigma')):
                    ax.plot(part.bin_number, part[key], 'o-', label=scenario)
                    ax.set(xlabel='Analysis bin', ylabel=key.replace('_', ' '))
            axes[0,0].axhline(0, color='.5', ls='--'); axes[0,1].axhline(1, color='.5', ls='--')
            axes[1,0].axhline(.682689, color='.5', ls='--'); axes[0,0].legend(fontsize=8)
            fig.suptitle(prod.PARAMETER_LABELS[p]+' — conditional injection tests')
            fig.tight_layout(); fig.savefig(plots/f'closure_summary_{p}.png', dpi=170); plt.close(fig)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bins', type=int, nargs='+', default=list(range(1,25)))
    p.add_argument('--studies', nargs='+', choices=('profiles', 'closure'), default=['profiles','closure'])
    p.add_argument('--scenarios', nargs='+', choices=SCENARIOS, default=list(SCENARIOS))
    p.add_argument('--replicas', type=int, default=100)
    p.add_argument('--profile-points', type=int, default=31)
    p.add_argument('--joint-points', type=int, default=15)
    p.add_argument('--extent', type=float, default=4., help='Scan half-width in nominal Hesse sigma')
    p.add_argument('--seed', type=int, default=20261119)
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--auxiliary-mode', choices=('fluctuate','fixed_truth'), default='fluctuate')
    p.add_argument('--output-dir', type=Path, default=HERE/'output/unpolarized_validation')
    p.add_argument('--production-output', type=Path, default=HERE/'output/asymmetry_extraction')
    p.add_argument('--cache', type=Path, help='Existing production event staging (provenance required)')
    p.add_argument('--rebuild-events', action='store_true', help='Rebuild staging from existing corrected ROOT files')
    p.add_argument('--run-info-csv', type=Path, default=HERE/prod.DEFAULT_RUN_INFO_CSV)
    p.add_argument('--cut-json', type=Path, default=HERE/prod.DEFAULT_CUT_JSON)
    p.add_argument('--dilution-json', type=Path)
    p.add_argument('--dilution-dir', type=Path, default=HERE/prod.DEFAULT_DILUTION_DIR)
    p.add_argument('--input', type=prod.parse_input_override, action='append', default=[], metavar='PERIOD=FILE')
    p.add_argument('--tree', default=prod.DEFAULT_TREE_NAME)
    p.add_argument('--chunk-size', default=prod.DEFAULT_CHUNK_SIZE)
    return p


def main():
    args = parser().parse_args()
    if (not args.bins or any(b not in range(1,25) for b in args.bins) or args.replicas < 1
            or args.profile_points < 9 or args.profile_points % 2 == 0
            or args.joint_points < 5 or args.joint_points % 2 == 0
            or args.extent <= 0 or not math.isfinite(args.extent) or args.workers < 1 or args.seed < 0):
        raise ValueError('Bins must be 1..24; replicas/workers positive; scan counts odd (>=9 and >=5); extent positive; seed nonnegative')
    args.bins = sorted(set(args.bins))
    args.scenarios = list(dict.fromkeys(args.scenarios))
    args.studies = sorted(set(args.studies))
    output = args.output_dir.expanduser().resolve()
    production_output = args.production_output.expanduser().resolve()
    if output == production_output or output in production_output.parents or production_output/'nominal' in output.parents or output == production_output/'nominal':
        raise ValueError('Choose a separate diagnostic output directory, outside nominal results')
    records = prod.parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    states = prod.run_state_arrays(records)
    cut_path = args.cut_json.expanduser().resolve()
    cuts = prod.load_channel_cuts(cut_path, cut_label='nominal')
    dilution_path = (args.dilution_json or prod.find_default_dilution_json(args.dilution_dir.expanduser().resolve())).resolve()
    prod.validate_dilution_cut_source(dilution_path, cut_path)
    dilutions = prod.load_dilution_factors(dilution_path, cut_label='nominal')
    inputs = {p: Path(v).expanduser().resolve() for p,v in prod.DEFAULT_INPUTS.items()}
    inputs.update({p: v.expanduser().resolve() for p,v in args.input})
    identity = prod._cache_identity(inputs, args.tree, records)
    cache = args.cache.expanduser().resolve() if args.cache else production_output/'nominal/cache/selected_events.npz'
    if args.rebuild_events or not prod._event_cache_path(cache).exists():
        print('Building selected-event staging from existing momentum-corrected ROOT inputs', flush=True)
        prod.build_event_cache(inputs, args.tree, args.chunk_size, records, cuts, cache)
    prod.validate_event_cache(cache, cuts, identity)
    config = {k:getattr(args,k) for k in ('bins','studies','scenarios','replicas','profile_points','joint_points','extent','seed','auxiliary_mode')}
    provenance = {'config':config, 'runtime':{'python':sys.version,'numpy':np.__version__,'iminuit':version('iminuit')}, 'diagnostic_sha256':digest(Path(__file__)),
                  'production_sha256':digest(Path(prod.__file__)), 'event_identity':identity,
                  'cuts_sha256':digest(cut_path), 'dilution_sha256':digest(dilution_path),
                  'run_info_sha256':digest(args.run_info_csv.expanduser().resolve())}
    signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    manifest_path = output/'manifest.json'
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old['signature'] != signature:
            raise RuntimeError('Diagnostic inputs/configuration changed. Use a new --output-dir; old checkpoints will not be reused.')
    elif output.exists() and any(output.iterdir()):
        raise RuntimeError('Output directory is nonempty without a matching manifest; choose a new directory')
    output.mkdir(parents=True, exist_ok=True)
    atomic_json(manifest_path, {'signature':signature, **provenance,
        'scope':'conditional estimator closure at fixed accepted kinematics; not a detector-model or omitted-harmonic test',
        'intervals':'data profiles: 1,3.84 (1D), 2.30,5.99 (2D); toy coverage: Hesse +/-1 and +/-1.96 sigma',
        'auxiliary_treatment':args.auxiliary_mode})
    print(f'Running {len(args.bins)} bins; studies={args.studies}; replicas={args.replicas}; workers={args.workers}', flush=True)
    results, failures = [], []
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=min(args.workers,len(args.bins)), mp_context=mp.get_context('spawn'),
            initializer=init_worker, initargs=(str(cache),states,dilutions,config,str(output))) as pool:
        futures = {pool.submit(run_bin,b):b for b in args.bins}
        for future in as_completed(futures):
            b = futures[future]
            try:
                result = future.result(); results.append(result)
                if result['status'] != 'complete':
                    failures.append({'bin_number':b, 'reason':result['status']})
                print(f'Bin {b}: {result["status"]} ({len(results)}/{len(args.bins)}); elapsed {(time.monotonic()-started)/60:.1f} min', flush=True)
            except Exception as exc:
                failures.append({'bin_number':b, 'reason':repr(exc)})
                print(f'Bin {b}: FAILED: {exc}', flush=True)
    results.sort(key=lambda r:r['bin_number'])
    write_products(results, config, output)
    skipped = [{'bin_number':r['bin_number'],'scenario':name,'reason':s.get('reason')}
               for r in results for name,s in r.get('scenarios',{}).items()
               if s.get('status')=='nonphysical_injection']
    for item in skipped:
        print(f'Nonphysical injection skipped: {item}', flush=True)
    atomic_json(output/'run_status.json', {'finished_bins':[r['bin_number'] for r in results],
                'failures':failures, 'nonphysical_injections':skipped})
    print(f'Diagnostics: {output}', flush=True)
    return 2 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
