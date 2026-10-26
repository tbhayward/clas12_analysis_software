#!/usr/bin/env python3
"""Paired CLAS12 DVCS MC |t| response, underflow and shape-closure diagnostics.

MCEvents is one entry per HIPO event. gen_t1/rec_t1 are proton-based t;
gen_kinematics_valid and rec_kinematics_valid are independent.
Only preliminary 1D, truth-valid, reconstructed-candidate selection is imposed.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def divide(a, b):
    return np.divide(a, b, out=np.full(np.broadcast_shapes(np.shape(a), np.shape(b)), np.nan), where=np.asarray(b) != 0)


def unfold(d, response, prior, niter=5):
    """D'Agostini iteration with explicit generated underflow/overflow bins."""
    t = np.maximum(prior.astype(float), 1e-12).copy()
    eff = response.sum(axis=0)
    for _ in range(niter):
        predicted = response @ t
        correction = response.T @ np.nan_to_num(divide(d, predicted))
        t = np.where(eff > 1e-12, t * correction / np.maximum(eff, 1e-12), 0.0)
    #endfor
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/production/fa18_inb/dvcs_rec/50nA/all.mc_events.root')
    ap.add_argument('--output', default='dvcs_low_t_migration_v3')
    ap.add_argument('--max-events', type=int, default=10000000)
    ap.add_argument('--chunk', type=int, default=200000)
    ap.add_argument('--t-edges', default='0.14,0.20,0.27,0.36,0.48,0.64,0.82,1.00', help='Comma-separated analysis |t| edges (GeV^2); REPLACE with exact analysis edges')
    ap.add_argument('--q2', default=None, help='Optional generated Q2 range: min,max (GeV^2)')
    ap.add_argument('--xb', default=None, help='Optional generated xB range: min,max')
    ap.add_argument('--iterations', type=int, default=5)
    ap.add_argument('--shape-slopes', default='-2,2', help='Pseudo-data reweight slopes in 1/GeV^2')
    args = ap.parse_args()
    import uproot
    path = Path(args.input)
    if path.is_dir():
        matches = sorted(path.glob('*.mc_events.root'))
        if len(matches) != 1:
            raise RuntimeError(f'Expected one *.mc_events.root in {path}; found {matches}')
        #endif
        path = matches[0]
    #endif
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    edges = np.asarray([float(v) for v in args.t_edges.split(',')])
    if len(edges) < 3 or not np.all(np.diff(edges) > 0):
        raise ValueError('Analysis edges must be strictly increasing')
    #endif
    n = len(edges) - 1
    m = n + 2  # underflow, analysis bins, overflow
    names = ['below'] + [f'{edges[i]:.3f}-{edges[i+1]:.3f}' for i in range(n)] + ['above']
    def index(t):
        return np.searchsorted(edges, t, side='right').astype(int)
    #endfor
    def parse_range(value):
        if value is None:
            return None
        #endif
        lo, hi = [float(x) for x in value.split(',')]
        return (lo, hi)
    #endfor
    q2range, xbrange = parse_range(args.q2), parse_range(args.xb)
    with uproot.open(path) as f:
        if 'MCEvents' not in f:
            raise RuntimeError(f'MCEvents tree not found; ROOT keys: {f.keys()}')
        #endif
        tree = f['MCEvents']
        mandatory = ['gen_kinematics_valid', 'rec_kinematics_valid', 'gen_t1', 'rec_t1']
        optional = ['gen_Q2', 'gen_x', 'rec_final_state', 'weight']
        missing = [k for k in mandatory if k not in tree.keys()]
        if missing:
            raise RuntimeError(f'Missing branches: {missing}; available: {list(tree.keys())}')
        #endif
        fields = mandatory + [k for k in optional if k in tree.keys()]
        if q2range and 'gen_Q2' not in fields or xbrange and 'gen_x' not in fields:
            raise RuntimeError('Requested Q2/xB selection but truth kinematic branch missing')
        #endif
        total = min(args.max_events, tree.num_entries)
        # Independent nominal response (even entries) and pseudo-data (odd entries).
        counts = np.zeros((2, m), dtype=float)
        counts2 = np.zeros((2, m), dtype=float)
        flow = np.zeros((2, m, m), dtype=float)
        flow2 = np.zeros((2, m, m), dtype=float)
        slopes = [float(s) for s in args.shape_slopes.split(',')]
        pseudo_truth = np.zeros((len(slopes), m))
        pseudo_reco = np.zeros((len(slopes), m))
        unweighted_truth = np.zeros(m)
        unweighted_flow = np.zeros((m, m))
        origin = np.zeros((n, 100), dtype=float)
        origin_edges = np.linspace(0, max(1.3, edges[-1] * 1.25), 101)
        n_truth = n_reco = n_inside = n_raw = 0
        negweights = 0
        for start in range(0, total, args.chunk):
            stop = min(start + args.chunk, total)
            a = tree.arrays(fields, entry_start=start, entry_stop=stop, library='np')
            n_raw += len(a['gen_t1'])
            gen = (a['gen_kinematics_valid'] > 0.5) & np.isfinite(a['gen_t1']) & (a['gen_t1'] < 0)
            if q2range:
                gen &= (a['gen_Q2'] >= q2range[0]) & (a['gen_Q2'] < q2range[1])
            #endif
            if xbrange:
                gen &= (a['gen_x'] >= xbrange[0]) & (a['gen_x'] < xbrange[1])
            #endif
            rec = (a['rec_kinematics_valid'] > 0.5) & np.isfinite(a['rec_t1']) & (a['rec_t1'] < 0)
            if 'rec_final_state' in a:
                rec &= a['rec_final_state'] > 0.5
            #endif
            g = -a['gen_t1'][gen]
            r = -a['rec_t1'][gen]
            rr = rec[gen]
            weights = a['weight'][gen] if 'weight' in a else np.ones(len(g))
            weights = np.where(np.isfinite(weights), weights, 0)
            negweights += int(np.count_nonzero(weights < 0))
            if np.any(weights < 0):
                raise RuntimeError('Negative MC weights detected: probability response requires a different treatment')
            #endif
            gi = index(g)
            ri = index(r[rr])
            partition = ((np.arange(start, stop)[gen] % 2) == 0).astype(int)
            n_truth += len(g)
            n_reco += int(rr.sum())
            n_inside += int(np.count_nonzero(rr & (index(r) >= 1) & (index(r) <= n)))
            np.add.at(unweighted_truth, gi, 1)
            np.add.at(unweighted_flow, (index(r[rr]), gi[rr]), 1)
            for split in (0, 1):
                sel = partition == split
                np.add.at(counts[split], gi[sel], weights[sel])
                np.add.at(counts2[split], gi[sel], weights[sel]**2)
                both = sel & rr
                np.add.at(flow[split], (index(r[both]), gi[both]), weights[both])
                np.add.at(flow2[split], (index(r[both]), gi[both]), weights[both]**2)
            #endfor
            # Pseudo-data use odd events; truth reweighted event by event.
            test = partition == 0
            for k, slope in enumerate(slopes):
                w = weights[test] * np.exp(np.clip(slope * (g[test] - edges[0]), -20, 20))
                np.add.at(pseudo_truth[k], gi[test], w)
                both = test & rr
                wrec = weights[both] * np.exp(np.clip(slope * (g[both] - edges[0]), -20, 20))
                np.add.at(pseudo_reco[k], index(r[both]), wrec)
            #endfor
            for j in range(min(n, 5)):
                selected = rr & (index(r) == j + 1)
                origin[j] += np.histogram(g[selected], bins=origin_edges, weights=weights[selected])[0]
            #endfor
            print(f'Processed {stop:,}/{total:,} HIPO events', flush=True)
        #endfor
    # Nominal response from odd entries; pseudo-data from even entries.
    denom = counts[1]
    response = np.divide(flow[1], denom[None, :], out=np.zeros_like(flow[1]), where=denom[None, :] > 0)
    # Conditional probability for a reconstructed event to originate in each generated bin.
    incoming = unweighted_flow[1:n+1]
    incoming_total = incoming.sum(axis=1)
    fraction_below = np.nan_to_num(divide(incoming[:, 0], incoming_total))
    fraction_same = np.nan_to_num(divide(np.diag(incoming[:, 1:n+1]), incoming_total))
    fraction_above = np.nan_to_num(divide(incoming[:, -1], incoming_total))
    fraction_other = 1 - fraction_below - fraction_same - fraction_above
    efficiency = np.nan_to_num(divide(unweighted_flow.sum(axis=0), unweighted_truth))
    retention = np.nan_to_num(divide(np.diag(unweighted_flow), unweighted_truth))
    fig, ax = plt.subplots(figsize=(9, 7))
    display = np.ma.masked_where(unweighted_flow == 0, unweighted_flow)
    from matplotlib.colors import LogNorm
    im = ax.imshow(display, origin='lower', aspect='auto', norm=LogNorm(vmin=1, vmax=max(2, unweighted_flow.max())), cmap='viridis')
    ax.set(xticks=np.arange(m), yticks=np.arange(m), xticklabels=names, yticklabels=names, xlabel='Generated |t| bin (GeV²)', ylabel='Reconstructed |t| bin (GeV²)', title='Paired DVCS MC migration (counts; unselected reconstructed DVCS candidates)')
    ax.tick_params(axis='x', rotation=50)
    fig.colorbar(im, ax=ax, label='Matched event count')
    fig.tight_layout(); fig.savefig(out/'01_t_migration_matrix.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(n)
    ax.plot(x, efficiency[1:n+1], 'o-', label='Reconstructed anywhere / generated in bin')
    ax.plot(x, retention[1:n+1], 's-', label='Reconstructed in same bin / generated in bin')
    ax.plot(x, np.nan_to_num(divide(unweighted_flow[1:n+1, 1:n+1].sum(axis=0), unweighted_truth[1:n+1])), '^-', label='Reconstructed within analysis / generated in bin')
    ax.set(xticks=x, xticklabels=names[1:-1], xlabel='Generated |t| bin (GeV²)', ylabel='Probability', ylim=(0, 1.05), title='Truth-normalized reconstruction and retention')
    ax.tick_params(axis='x', rotation=35); ax.legend(fontsize=9); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(out/'02_reconstruction_efficiency.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5))
    bottom = np.zeros(n)
    for values, label in [(fraction_below, 'Generated below analysis'), (fraction_same, 'Same generated bin'), (fraction_other, 'Other generated analysis bin'), (fraction_above, 'Generated above analysis')]:
        ax.bar(x, values, bottom=bottom, label=label)
        bottom += values
    #endfor
    ax.set(xticks=x, xticklabels=names[1:-1], xlabel='Reconstructed |t| bin (GeV²)', ylabel='Fraction of reconstructed yield', ylim=(0, 1.05), title='Origins of events in each reconstructed bin')
    ax.tick_params(axis='x', rotation=35); ax.legend(fontsize=9, loc='upper right')
    fig.tight_layout(); fig.savefig(out/'03_incoming_migration_fractions.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5))
    centers = (origin_edges[:-1] + origin_edges[1:])/2
    for j in range(min(n, 5)):
        if origin[j].sum():
            ax.step(centers, origin[j]/origin[j].sum(), where='mid', label=f'Rec {names[j+1]}')
        #endif
    #endfor
    ax.axvline(edges[0], ls='--', c='black', label='Analysis minimum')
    ax.set(xlabel='Generated |t| (GeV²)', ylabel='Fraction per generated-|t| histogram bin', title='Generated origins of lowest reconstructed |t| bins')
    ax.legend(fontsize=9); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(out/'04_generated_origins_low_t.png', dpi=180); plt.close(fig)
    # Conventional correction is based on reconstructed analysis-bin counts divided by generated same-bin counts.
    acceptance = np.nan_to_num(divide(flow[1, 1:n+1, :].sum(axis=1), denom[1:n+1]))
    fig, axs = plt.subplots(len(slopes), 1, figsize=(10, 4.2*len(slopes)), squeeze=False)
    closure_rows = []
    for k, slope in enumerate(slopes):
        truth = pseudo_truth[k, 1:n+1]
        conventional = np.nan_to_num(divide(pseudo_reco[k, 1:n+1], acceptance))
        unfolded = unfold(pseudo_reco[k], response, np.maximum(denom, 1e-9), args.iterations)[1:n+1]
        conv_bias = 100 * (divide(conventional, truth) - 1)
        matrix_bias = 100 * (divide(unfolded, truth) - 1)
        ax = axs[k, 0]
        ax.axhline(0, color='black', lw=1)
        ax.plot(x, conv_bias, 'o-', label='Conventional bin-by-bin')
        ax.plot(x, matrix_bias, 's-', label=f'Iterative response ({args.iterations} iterations)')
        ax.set(xticks=x, xticklabels=names[1:-1], ylabel='Closure bias (%)', title=f'Generator reweight slope {slope:+g} (GeV⁻²); independent even/odd samples')
        ax.tick_params(axis='x', rotation=35); ax.legend(); ax.grid(alpha=.2)
        for j in range(n):
            closure_rows.append(dict(slope=slope, bin=j+1, t_low=edges[j], t_high=edges[j+1], truth=truth[j], conventional=conventional[j], matrix=unfolded[j], conventional_bias_pct=conv_bias[j], matrix_bias_pct=matrix_bias[j]))
        #endfor
        print(f'Slope {slope:+g}: first-bin conventional bias {conv_bias[0]:+.2f}%; response bias {matrix_bias[0]:+.2f}%')
    #endfor
    fig.tight_layout(); fig.savefig(out/'05_generator_shape_closure.png', dpi=180); plt.close(fig)
    with open(out/'bin_summary.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['bin', 't_low', 't_high', 'generated', 'reconstructed_in_bin', 'efficiency_anywhere', 'same_bin_probability', 'fraction_below', 'fraction_same', 'fraction_other', 'fraction_above'])
        writer.writeheader()
        for j in range(n):
            writer.writerow(dict(bin=j+1, t_low=edges[j], t_high=edges[j+1], generated=unweighted_truth[j+1], reconstructed_in_bin=incoming_total[j], efficiency_anywhere=efficiency[j+1], same_bin_probability=retention[j+1], fraction_below=fraction_below[j], fraction_same=fraction_same[j], fraction_other=fraction_other[j], fraction_above=fraction_above[j]))
        #endfor
    #endfor
    with open(out/'closure_summary.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(closure_rows[0]))
        writer.writeheader(); writer.writerows(closure_rows)
    #endfor
    np.savetxt(out/'response_matrix.csv', response, delimiter=',', header='Rows reconstructed: below, analysis bins, above; columns generated: below, analysis bins, above')
    print(f'\nRaw HIPO events read: {n_raw:,}; truth-valid selected: {n_truth:,}; reconstructed DVCS candidates: {n_reco:,}; reconstructed in analysis: {n_inside:,}')
    print(f'First reconstructed bin from generated below threshold: {fraction_below[0]:.3%} ({incoming[0,0]:.0f}/{incoming_total[0]:.0f})')
    for j in range(min(5, n)):
        print(f'  Rec bin {j+1}: below={fraction_below[j]:.2%}, other-inside={fraction_other[j]:.2%}, above={fraction_above[j]:.2%}, same={fraction_same[j]:.2%}')
    #endfor
    print(f'Generated analysis events lost outside reconstructed analysis range or reconstruction: {np.sum(unweighted_truth[1:n+1]) - np.sum(unweighted_flow[1:n+1,1:n+1]):.0f}')
    print(f'Response condition number (analysis submatrix): {np.linalg.cond(response[1:n+1,1:n+1]):.3g}')
    if np.any(denom < 100) or np.any(incoming_total < 100):
        print('WARNING: Some response bins have fewer than 100 nominal generated or reconstructed events.')
    #endif
    print('WARNING: 1D candidate-level study, not the full 4D analysis/exclusivity selection; analysis edges are user-supplied.')
    print('WARNING: Pseudo-data and response are independent, but statistical closure fluctuations remain; finite-iteration unfolding is prior dependent.')
    print(f'PNG and CSV output: {out.resolve()}')


if __name__ == '__main__':
    main()
