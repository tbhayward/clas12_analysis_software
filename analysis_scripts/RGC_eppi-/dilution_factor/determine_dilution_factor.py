#!/usr/bin/env python3
"""Nominal per-period, per-bin ND3 Method 1 dilution and packing cross-check."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np
from model import COEFFICIENTS, SOURCE, method1, transfer_ratio

HERE = Path(__file__).resolve().parent
PERIODS = ['summer2022', 'fall2022', 'spring2023']
LABELS = ['Su22', 'Fa22', 'Sp23']
TARGETS = ['ND3', 'C', 'CH2', 'He', 'ET', 'CD2']
XB = [.10, .25, .35, .45, .60]
TP = [.05, .25, .45, .65, .85, 1.05, 1.25]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_charges(path):
    """Only correct ND3/control sections. Charge uses helicity +/- columns."""
    periods = dict(zip(LABELS, PERIODS))
    result = {}
    key = None
    for number, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if line.startswith('#'):
            m = re.fullmatch(r'#\s*RGC (Su22|Fa22|Sp23) (?:Inb )?(ND3|C|CH2|CD2|He|ET)(?: Runs)?', line)
            key = (periods[m[1]], m[2]) if m else None
        elif key and line:
            row = next(csv.reader([line]))
            run = int(row[0])
            charge = float(row[2]) + float(row[3])
            if not np.isfinite(charge) or charge < 0:
                raise ValueError(f'Invalid charge at line {number}')
            values = result.setdefault(key, {})
            if run in values:
                raise ValueError(f'Duplicate run {run} in {key}')
            values[run] = charge
    return result


def load_inputs(path, charge_path):
    sys.path.insert(0, str(HERE.parent/'channel_selection'))
    from selection import load_selected_inputs
    registry = load_selected_inputs(path)
    import uproot
    charge_table = load_charges(charge_path)
    counts = {}
    charges = {}
    audit = []
    for period in PERIODS:
        counts[period], charges[period] = {}, {}
        for target in TARGETS:
            if target == 'CD2' and period != 'spring2023':
                continue
            name = f'RGC_{period}_{target}'
            item = registry['files'][name]
            if item.get('unavailable_bins'):
                raise ValueError(f'{name}: unavailable nominal cuts')
            with uproot.open(item['path']) as root:
                manifest = json.loads(str(root['selection_manifest']))
                if manifest['partial'] or manifest['read_rows'] != manifest['source_entries']:
                    raise ValueError(f'{name}: incomplete source scan')
                arrays = root[registry['tree']].arrays(['analysis_bin', 'runnum'], library='np')
            bins = arrays['analysis_bin']
            if np.any((bins < 1) | (bins > 24)):
                raise ValueError(f'{name}: invalid analysis bins')
            observed = np.bincount(bins, minlength=25)[1:25]
            if not np.array_equal(observed, np.asarray(manifest['counts'])[:, 1]):
                raise ValueError(f'{name}: selected counts differ from export manifest')
            # All runs present BEFORE exclusivity, including runs with zero selected rows.
            runs = sorted(int(r) for r in manifest['runs'])
            if not set(map(int, arrays['runnum'])).issubset(runs):
                raise ValueError(f'{name}: selected run missing from source coverage')
            table = charge_table[(period, target)]
            missing = sorted(set(runs)-table.keys())
            if missing:
                raise ValueError(f'{name}: missing charge for source runs {missing}')
            if any(table[r] <= 0 for r in runs):
                raise ValueError(f'{name}: source contains runs with zero good charge')
            total = sum(table[r] for r in runs)
            if total <= 0:
                raise ValueError(f'{name}: no normalization charge')
            counts[period][target] = observed
            charges[period][target] = total
            audit.append(dict(dataset=name, charge=total, runs=runs,
                              run_charges={str(r): table[r] for r in runs},
                              selected_counts=observed.tolist(), selected_rows=int(observed.sum()),
                              selected_path=item['path'], source_entries=manifest['source_entries'],
                              reference_assisted_bins=manifest['reference_assisted_bins'],
                              excluded_csv_runs=sorted(set(table)-set(runs))))
            print(f'{LABELS[PERIODS.index(period)]} {target}: {observed.sum():,} combinations, {len(runs)} source runs, charge {total:,.3f}', flush=True)
    return registry, counts, charges, audit


def rates_from_counts(counts, charges):
    return {p: {t: n/charges[p][t] for t, n in counts[p].items()} for p in PERIODS}


def evaluate(rates):
    ratio = transfer_ratio(rates['spring2023'])
    result = {p: method1(rates[p], None if p == 'spring2023' else ratio) for p in PERIODS}
    return ratio, result


def resample(counts, charges, replicas, seed):
    """Reuse each Sp23 replica across all periods, preserving transfer correlations."""
    rng = np.random.default_rng(seed)
    samples = {p: {t: rng.poisson(n, (replicas, 24)) for t, n in counts[p].items()} for p in PERIODS}
    ratio, result = evaluate(rates_from_counts(samples, charges))
    return ratio, result


def describe(nominal, draws):
    valid = draws[np.isfinite(draws)]
    return dict(value=float(nominal) if np.isfinite(nominal) else None,
                stat_std=float(valid.std(ddof=1)) if len(valid) > 1 else None,
                stat_q16=float(np.quantile(valid, .16)) if len(valid) else None,
                stat_q84=float(np.quantile(valid, .84)) if len(valid) else None,
                finite_fraction=float(len(valid)/len(draws)))


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plots(output, rows, transfer_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder = output/'figures'
    folder.mkdir(exist_ok=True)
    colors = ['#0072B2', '#D55E00', '#009E73']
    for variable, ylabel in [('dilution', 'Dilution factor'), ('packing_fraction', 'Cell packing fraction')]:
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
        for ix, ax in enumerate(axes.flat):
            for label, color in zip(LABELS, colors):
                values = [r for r in rows if r['period'] == label and (r['bin']-1)//6 == ix]
                x = (np.asarray(TP[:-1])+TP[1:])/2
                y = np.array([r[variable] if r[variable] is not None else np.nan for r in values])
                # Quantile intervals remain readable when rare denominator crossings inflate SD.
                lo = np.array([r[variable+'_q16'] if r[variable+'_q16'] is not None else np.nan for r in values])
                hi = np.array([r[variable+'_q84'] if r[variable+'_q84'] is not None else np.nan for r in values])
                ax.plot(x, y, 'o', ms=4, color=color, label=label)
                ax.vlines(x, lo, hi, color=color, alpha=.7)
            ax.set_title(f'{XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}')
            ax.axhline(0, color='gray', lw=.6)
            ax.axhline(1, color='gray', lw=.6, ls=':')
            ax.set(xlabel='−t′ (GeV²)', ylabel=ylabel)
        axes.flat[0].legend()
        fig.tight_layout()
        fig.savefig(folder/(variable+'.png'), dpi=170)
        plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    for ix, ax in enumerate(axes.flat):
        v = transfer_rows[ix*6:(ix+1)*6]
        x = (np.asarray(TP[:-1])+TP[1:])/2
        ax.plot(x, [r['value'] for r in v], 'o', color=colors[0])
        ax.vlines(x, [r['stat_q16'] for r in v], [r['stat_q84'] for r in v], color=colors[0])
        ax.set(title=f'{XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}', xlabel='−t′ (GeV²)', ylabel='Sp23 CD₂ / CH₂ rate ratio')
    fig.tight_layout()
    fig.savefig(folder/'cd2_ch2_transfer.png', dpi=170)
    plt.close(fig)


def control_distributions(output, registry, charges):
    """Sp23 measured yields, charge normalized. No synthetic CD2 or shape scaling."""
    import uproot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    edges = np.linspace(TP[0], TP[-1], 25)
    centers = (edges[:-1]+edges[1:])/2
    histograms = {}
    for target in ['CD2', 'CH2']:
        name = 'RGC_spring2023_'+target
        hist = np.zeros((4, 24), dtype=np.int64)
        with uproot.open(registry['files'][name]['path']) as root:
            for a in root[registry['tree']].iterate(['x', 'tprime'], library='np', step_size='20 MB'):
                for ix in range(4):
                    selected = (a['x'] >= XB[ix]) & (a['x'] < XB[ix+1])
                    hist[ix] += np.histogram(-a['tprime'][selected], bins=edges)[0]
        if hist.sum() != registry['files'][name]['selected_rows']:
            raise ValueError(name+': distribution histogram lost selected candidates')
        histograms[target] = hist
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    rows = []
    for ix, ax in enumerate(axes.flat):
        for target, color, marker in [('CD2', '#0072B2', 'o'), ('CH2', '#D55E00', 's')]:
            n = histograms[target][ix]
            q = charges['spring2023'][target]
            ax.errorbar(centers, 1e6*n/q, yerr=1e6*np.sqrt(n)/q,
                        fmt=marker, linestyle='none', color=color, markersize=4,
                        elinewidth=.9, label='Sp23 '+target.replace('2', '₂'))
            for j, count in enumerate(n):
                rows.append(dict(target=target, xB_low=XB[ix], xB_high=XB[ix+1],
                                 minus_tprime_low=edges[j], minus_tprime_high=edges[j+1],
                                 count=int(count), charge=q, rate=count/q,
                                 rate_stat=np.sqrt(count)/q))
        ax.set(title=f'{XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}',
               xlabel='−t′ (GeV²)', ylabel='10⁶ N / Q per bin')
    ymax = max(np.max(1e6*(histograms[t]+np.sqrt(histograms[t]))/charges['spring2023'][t]) for t in histograms)
    axes.flat[0].set_ylim(0, 1.08*ymax if ymax > 0 else 1)
    axes.flat[0].legend()
    fig.tight_layout()
    (output/'figures').mkdir(exist_ok=True)
    fig.savefig(output/'figures/cd2_ch2_distributions.png', dpi=170)
    plt.close(fig)
    write_csv(output/'cd2_ch2_distributions.csv', rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selected-inputs', type=Path, default=HERE.parent/'channel_selection/output/selected_inputs.json')
    parser.add_argument('--run-info', type=Path, default=HERE.parent/'imports/clas12_run_info_ND3_source.csv')
    parser.add_argument('--output', type=Path, default=HERE/'output')
    parser.add_argument('--replicas', type=int, default=20000)
    parser.add_argument('--seed', type=int, default=20261124)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error('Use at least 1000 statistical replicas')
    registry, counts, charges, audit = load_inputs(args.selected_inputs, args.run_info)
    rates = rates_from_counts(counts, charges)
    ratio, nominal = evaluate(rates)
    ratio_draws, draws = resample(counts, charges, args.replicas, args.seed)
    rows, details = [], []
    for period, label in zip(PERIODS, LABELS):
        for b in range(24):
            summaries = {key: describe(nominal[period][key][b], draws[period][key][:, b])
                         for key in ['dilution', 'packing_fraction']}
            flags, notes = [], []
            for target, n in counts[period].items():
                if n[b] < 10:
                    notes.append(f'low_count_{target}')
            if period != 'spring2023' and counts['spring2023']['CH2'][b] < 10:
                flags.append('weak_transfer_CH2')
            for key, s in summaries.items():
                if s['value'] is None:
                    flags.append(key+'_undefined')
                elif not 0 < s['value'] < 1:
                    flags.append(key+'_outside_0_1')
                if s['finite_fraction'] < .99:
                    flags.append(key+'_unstable_replicas')
            denominator = draws[period]['denominator'][:, b]
            nominal_den = nominal[period]['denominator'][b]
            crossing = float(np.mean(denominator*nominal_den <= 0))
            if crossing > .01:
                flags.append('denominator_sign_unstable')
            row = dict(period=label, bin=b+1, xB_low=XB[b//6], xB_high=XB[b//6+1],
                       minus_tprime_low=TP[b%6], minus_tprime_high=TP[b%6+1])
            for key, s in summaries.items():
                row.update({key:s['value'], key+'_stat':s['stat_std'], key+'_q16':s['stat_q16'],
                            key+'_q84':s['stat_q84'], key+'_finite_fraction':s['finite_fraction']})
            row.update(denominator_crossing_fraction=crossing, ready_for_asymmetry=not bool(flags), flags='|'.join(flags), notes='|'.join(notes))
            rows.append(row)
            details.append(dict(period=period, bin=b+1, **summaries, flags=flags, notes=notes, ready_for_asymmetry=not bool(flags),
                                denominator_crossing_fraction=crossing))
    packing_summary = []
    # Replay the same seeded draws once for the integrated packing cross-check.
    rng = np.random.default_rng(args.seed)
    sampled_rates = {p: {t: rng.poisson(n, (args.replicas, 24))/charges[p][t]
                         for t,n in counts[p].items()} for p in PERIODS}
    for period, label in zip(PERIODS, LABELS):
        integrated = {t: n.sum() for t, n in rates[period].items()}
        # Sum after each bin's transfer, rather than using a ratio of total yields.
        integrated['CD2'] = nominal[period]['effective_cd2_rate'].sum()
        central = method1(integrated)['packing_fraction']
        integrated_draw = {t: n.sum(axis=1) for t,n in sampled_rates[period].items()}
        integrated_draw['CD2'] = draws[period]['effective_cd2_rate'].sum(axis=1)
        sampled = method1(integrated_draw)['packing_fraction']
        packing_summary.append(dict(period=label, **describe(central, sampled)))
    transfer = [dict(bin=b+1, **describe(ratio[b], ratio_draws[:, b])) for b in range(24)]
    # Identical replicas are used in both period/bin covariance and transfer uncertainties.
    matrices = {key: np.concatenate([draws[p][key] for p in PERIODS], axis=1)
                for key in ['dilution', 'packing_fraction']}
    covariance = {}
    valid_counts = {}
    for key, matrix in matrices.items():
        valid = np.all(np.isfinite(matrix), axis=1)
        valid_counts[key] = int(valid.sum())
        covariance[key] = np.cov(matrix[valid], rowvar=False) if valid.sum() > 1 else np.full((72, 72), np.nan)
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out/'method1.csv', rows)
    write_csv(out/'cd2_ch2_transfer.csv', transfer)
    write_csv(out/'packing_fraction_integrated.csv', packing_summary)
    payload = dict(schema='eppi_minus_dilution_method1_v1', source=SOURCE,
                   coefficients=COEFFICIENTS, cutset_id=registry['cutset_id'],
                   selected_inputs=str(args.selected_inputs.resolve()), run_info=str(args.run_info.resolve()),
                   run_info_sha256=digest(args.run_info), seed=args.seed, replicas=args.replicas,
                   normalization='Q+ + Q- summed over original pre-exclusivity ROOT run coverage',
                   transfer_prescription='Sp23 charge-normalized CD2/CH2 ratio separately in each bin',
                   packing='5.86/5 times bath packing fraction, Eq. 10. Cross-check only.',
                   uncertainty='Independent Poisson combination counts. Shared Sp23 replicas in all periods. Statistical only.',
                   covariance_order=[f'{p}:bin{b+1}' for p in PERIODS for b in range(24)],
                   covariance_complete_replicas=valid_counts, inputs=audit,
                   bins=details, transfer=transfer, packing_fraction_integrated=packing_summary)
    (out/'method1.json').write_text(json.dumps(clean(payload), indent=2, allow_nan=False)+'\n')
    np.savez_compressed(out/'statistical_covariance.npz', **covariance)
    plots(out, rows, transfer)
    control_distributions(out, registry, charges)
    flagged = sum(bool(r['flags']) for r in rows)
    print(f'Finished: 72 period/bin estimates, {flagged} flagged for review. Outputs: {out.resolve()}')
    for row in packing_summary:
        print(f"{row['period']}: integrated cell packing fraction {row['value']:.3f}, central 68% interval [{row['stat_q16']:.3f}, {row['stat_q84']:.3f}]")
    print('Statistical intervals only. No scale uncertainty or radiative correction has been assigned.')


if __name__ == '__main__':
    main()
