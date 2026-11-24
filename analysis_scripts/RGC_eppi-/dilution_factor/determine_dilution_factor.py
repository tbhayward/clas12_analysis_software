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
    fig = plt.figure(figsize=(12, 10))
    outer = fig.add_gridspec(2, 2, hspace=.28, wspace=.20)
    tops, bottoms, rows, ratio_rows = [], [], [], []
    for ix in range(4):
        grid = outer[ix//2, ix%2].subgridspec(2, 1, height_ratios=[2, 1], hspace=.04)
        ax = fig.add_subplot(grid[0], sharey=tops[0] if tops else None)
        ratio_ax = fig.add_subplot(grid[1], sharex=ax, sharey=bottoms[0] if bottoms else None)
        tops.append(ax)
        bottoms.append(ratio_ax)
        for target, color, marker in [('CD2', '#0072B2', 'o'), ('CH2', '#D55E00', 's')]:
            n = histograms[target][ix]
            q = charges['spring2023'][target]
            ax.errorbar(centers, n/q, yerr=np.sqrt(n)/q,
                        fmt=marker, linestyle='none', color=color, markersize=4,
                        elinewidth=.9, label='Sp23 '+target.replace('2', '₂'))
            for j, count in enumerate(n):
                rows.append(dict(target=target, xB_low=XB[ix], xB_high=XB[ix+1],
                                 minus_tprime_low=edges[j], minus_tprime_high=edges[j+1],
                                 count=int(count), charge=q, rate=count/q,
                                 rate_stat=np.sqrt(count)/q))
        # Use the six analysis tprime bins for the ratio's statistical precision.
        nd = histograms['CD2'][ix].reshape(6, 4).sum(axis=1)
        nh = histograms['CH2'][ix].reshape(6, 4).sum(axis=1)
        ratio_edges = edges[::4]
        ratio_centers = (ratio_edges[:-1]+ratio_edges[1:])/2
        scale = charges['spring2023']['CH2']/charges['spring2023']['CD2']
        valid = nh > 0
        ratios = np.full(6, np.nan)
        errors = np.full(6, np.nan)
        ratios[valid] = scale*nd[valid]/nh[valid]
        errors[valid] = scale*np.sqrt(nd[valid]/nh[valid]**2 + nd[valid]**2/nh[valid]**3)
        ratio_ax.errorbar(ratio_centers[valid], ratios[valid], yerr=errors[valid],
                         fmt='o', linestyle='none', color='#333333', markersize=3,
                         elinewidth=.8)
        ratio_ax.axhline(1, color='gray', lw=.7, ls=':')
        ratio_ax.set(xlabel='−t′ (GeV²)', ylabel='CD₂ / CH₂')
        ax.set(title=f'{XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}', ylabel='Counts / nC')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        ax.tick_params(axis='x', labelbottom=False)
        for j in range(6):
            ratio_rows.append(dict(xB_low=XB[ix], xB_high=XB[ix+1],
                                   minus_tprime_low=ratio_edges[j], minus_tprime_high=ratio_edges[j+1],
                                   ratio=float(ratios[j]) if valid[j] else None,
                                   ratio_stat=float(errors[j]) if valid[j] else None))
    ymax = max(np.max((histograms[t]+np.sqrt(histograms[t]))/charges['spring2023'][t]) for t in histograms)
    tops[0].set_ylim(0, 1.08*ymax if ymax > 0 else 1)
    finite = [r['ratio']+r['ratio_stat'] for r in ratio_rows if r['ratio'] is not None]
    bottoms[0].set_ylim(0, 1.05*max(finite) if finite else 2)
    tops[0].legend()
    fig.subplots_adjust(left=.09, right=.98, top=.95, bottom=.07)
    (output/'figures').mkdir(exist_ok=True)
    fig.savefig(output/'figures/cd2_ch2_distributions.png', dpi=170)
    plt.close(fig)
    write_csv(output/'cd2_ch2_distributions.csv', rows)
    write_csv(output/'cd2_ch2_distribution_ratios.csv', ratio_rows)


def xb_comparison(output, counts, charges, replicas, seed, cutset_id):
    """Integrate target counts over tprime BEFORE transfer and Method 1.

    Comparison only. Does not replace the 24-bin nominal results.
    """
    integrated = {p: {t: n.reshape(4, 6).sum(axis=1) for t,n in targets.items()}
                  for p,targets in counts.items()}
    ratio, nominal = evaluate(rates_from_counts(integrated, charges))
    rng = np.random.default_rng(seed)
    sample = {p: {t: rng.poisson(n, (replicas, 24)).reshape(replicas, 4, 6).sum(axis=2)
                  for t,n in targets.items()} for p,targets in counts.items()}
    ratio_draws, draws = evaluate(rates_from_counts(sample, charges))
    rows = []
    for period, label in zip(PERIODS, LABELS):
        for ix in range(4):
            row = dict(period=label, xB_low=XB[ix], xB_high=XB[ix+1])
            for key in ['dilution', 'packing_fraction']:
                row.update({key+'_'+k: v for k,v in describe(nominal[period][key][ix], draws[period][key][:,ix]).items()})
            row['denominator_crossing_fraction'] = float(np.mean(draws[period]['denominator'][:,ix]*nominal[period]['denominator'][ix] <= 0))
            rows.append(row)
    write_csv(output/'method1_xb_comparison.csv', rows)
    payload = dict(cutset_id=cutset_id, status='comparison_only',
                   prescription='Sum counts over all six tprime bins before computing transfer and Method 1',
                   bins=rows, transfer=[dict(xB_low=XB[ix], xB_high=XB[ix+1],
                   **describe(ratio[ix], ratio_draws[:,ix])) for ix in range(4)])
    (output/'method1_xb_comparison.json').write_text(json.dumps(clean(payload), indent=2, allow_nan=False)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = (np.asarray(XB[:-1])+XB[1:])/2
    for ax,key,label in zip(axes, ['dilution','packing_fraction'], ['Dilution factor','Cell packing fraction']):
        for period,color in zip(LABELS, ['#0072B2','#D55E00','#009E73']):
            values = [r for r in rows if r['period']==period]
            ax.plot(x, [r[key+'_value'] for r in values], 'o', color=color, label=period)
            ax.vlines(x, [r[key+'_stat_q16'] for r in values], [r[key+'_stat_q84'] for r in values], color=color)
        ax.set(xlabel='xB', ylabel=label)
        ax.set_ylim(bottom=0)
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(output/'figures/method1_xb_comparison.png', dpi=170)
    plt.close(fig)
    print('tprime-integrated comparison:')
    for r in rows:
        print(f"{r['period']} xB [{r['xB_low']:.2f},{r['xB_high']:.2f}): DF {r['dilution_value']:.3f}, PF {r['packing_fraction_value']:.3f}, denominator flips {100*r['denominator_crossing_fraction']:.2f}%")


def comparison_groups():
    return {
        'xb_only': [list(range(6*ix,6*ix+6)) for ix in range(4)],
        'tprime_only': [[6*ix+it for ix in range(4)] for it in range(6)],
        'xb_two_tprime': [list(range(6*ix+3*half,6*ix+3*half+3))
                          for ix in range(4) for half in range(2)],
    }


def group_counts(counts, groups):
    """Aggregate the final dimension of either counts or count replicas."""
    return {p: {t: np.stack([n[...,ids].sum(axis=-1) for ids in groups], axis=-1)
                for t,n in targets.items()} for p,targets in counts.items()}


def binning_comparison(output, counts, charges, replicas, seed, cutset_id):
    """Compare pooling choices with shared Poisson replicas and measured transfers."""
    rng = np.random.default_rng(seed)
    sample = {p: {t: rng.poisson(n, (replicas,24)) for t,n in targets.items()}
              for p,targets in counts.items()}
    rows, differences, scheme_meta, covariances = [], [], {}, {}
    all_nominal, all_draws = {}, {}
    for scheme,groups in comparison_groups().items():
        ratio, nominal = evaluate(rates_from_counts(group_counts(counts, groups), charges))
        ratio_draw, draws = evaluate(rates_from_counts(group_counts(sample, groups), charges))
        all_nominal[scheme], all_draws[scheme] = nominal, draws
        scheme_meta[scheme] = dict(original_bins=[[b+1 for b in ids] for ids in groups],
                                  transfer=[describe(ratio[j],ratio_draw[:,j]) for j in range(len(groups))])
        for period,label in zip(PERIODS,LABELS):
            for j,ids in enumerate(groups):
                ix = sorted({b//6 for b in ids})
                it = sorted({b%6 for b in ids})
                row = dict(scheme=scheme, period=label, group=j+1,
                           xB_low=XB[min(ix)], xB_high=XB[max(ix)+1],
                           minus_tprime_low=TP[min(it)], minus_tprime_high=TP[max(it)+1])
                flags=[]
                for key in ['dilution','packing_fraction']:
                    stats=describe(nominal[period][key][j],draws[period][key][:,j])
                    row.update({key+'_'+k:v for k,v in stats.items()})
                    if stats['value'] is None or not 0<stats['value']<1:
                        flags.append(key+'_unphysical')
                    if stats['finite_fraction']<.99:
                        flags.append(key+'_undefined_replicas')
                crossing=float(np.mean(draws[period]['denominator'][:,j]*nominal[period]['denominator'][j]<=0))
                row['denominator_crossing_fraction']=crossing
                if crossing>.01:flags.append('denominator_sign_unstable')
                row['flags']='|'.join(flags)
                rows.append(row)
        joint=np.concatenate([draws[p][key] for key in ['dilution','packing_fraction'] for p in PERIODS],axis=1)
        valid=np.all(np.isfinite(joint),axis=1)
        covariances[scheme]=np.cov(joint[valid],rowvar=False) if valid.sum()>1 else np.full((joint.shape[1],joint.shape[1]),np.nan)
        scheme_meta[scheme]['covariance_order']=[f'{key}:{p}:group{j+1}' for key in ['dilution','packing_fraction'] for p in PERIODS for j in range(len(groups))]
        scheme_meta[scheme]['complete_covariance_replicas']=int(valid.sum())
    for period,label in zip(PERIODS,LABELS):
        for ix in range(4):
            low,high=2*ix,2*ix+1
            n=all_nominal['xb_two_tprime'][period]['dilution']
            draws=all_draws['xb_two_tprime'][period]['dilution']
            stats=describe(n[high]-n[low],draws[:,high]-draws[:,low])
            delta_draws = draws[:,high]-draws[:,low]
            finite_delta = delta_draws[np.isfinite(delta_draws)]
            q025,q975 = np.quantile(finite_delta,[.025,.975])
            central_width = (stats['stat_q84']-stats['stat_q16'])/2
            differences.append(dict(period=label,xB_low=XB[ix],xB_high=XB[ix+1],
                                    difference='high_tprime_minus_low_tprime',**stats,
                                    stat_q025=float(q025),stat_q975=float(q975),
                                    heavy_tail_warning=bool(stats['stat_std']>3*central_width)))
    write_csv(output/'method1_binning_comparison.csv', rows)
    write_csv(output/'method1_tprime_half_differences.csv', differences)
    payload=dict(cutset_id=cutset_id,status='comparison_only',replicas=replicas,seed=seed,
                 prescription='Sum target counts within each group before calculating its Sp23 CD2/CH2 transfer and Method 1',
                 half_split_gev2=.65,schemes=scheme_meta,bins=rows,half_differences=differences)
    (output/'method1_binning_comparison.json').write_text(json.dumps(clean(payload),indent=2,allow_nan=False)+'\n')
    np.savez_compressed(output/'binning_comparison_covariance.npz',**covariances)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#0072B2','#D55E00','#009E73']
    fig,ax=plt.subplots(figsize=(8,5))
    for label,color in zip(LABELS,colors):
        v=[r for r in rows if r['scheme']=='tprime_only' and r['period']==label]
        x=[(r['minus_tprime_low']+r['minus_tprime_high'])/2 for r in v]
        ax.plot(x,[r['dilution_value'] for r in v],'o',color=color,label=label)
        ax.vlines(x,[r['dilution_stat_q16'] for r in v],[r['dilution_stat_q84'] for r in v],color=color)
    ax.set(xlabel='−t′ (GeV²)',ylabel='Dilution factor',title='Integrated over 0.10 ≤ xB < 0.60')
    ax.legend();fig.tight_layout()
    fig.savefig(output/'figures/method1_tprime_only.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
    for ix,ax in enumerate(axes.flat):
        for ip,(label,color) in enumerate(zip(LABELS,colors)):
            v=[r for r in rows if r['scheme']=='xb_two_tprime' and r['period']==label and r['xB_low']==XB[ix]]
            # Small horizontal offsets keep the three periods visible in each half.
            x=np.array([(r['minus_tprime_low']+r['minus_tprime_high'])/2 for r in v])+(ip-1)*.025
            ax.plot(x,[r['dilution_value'] for r in v],'o',color=color,label=label)
            ax.vlines(x,[r['dilution_stat_q16'] for r in v],[r['dilution_stat_q84'] for r in v],color=color)
        ax.axvline(.65,color='gray',lw=.7,ls=':')
        ax.set(xlabel='−t′ (GeV²)',ylabel='Dilution factor',title=f'{XB[ix]:.2f} ≤ xB < {XB[ix+1]:.2f}', xlim=(TP[0],TP[-1]))
    axes.flat[0].legend();fig.tight_layout()
    fig.savefig(output/'figures/method1_two_tprime_halves.png',dpi=170);plt.close(fig)
    for scheme in comparison_groups():
        selected=[r for r in rows if r['scheme']==scheme]
        print(f"{scheme}: {sum(bool(r['flags']) for r in selected)}/{len(selected)} flagged, max denominator flips {100*max(r['denominator_crossing_fraction'] for r in selected):.3f}%")
    return rows,differences

def period_consistency(output):
    """Approximate correlated GLS checks of one-dimensional dilution projections."""
    from scipy.stats import chi2
    payload=json.loads((output/'method1_binning_comparison.json').read_text())
    covariance=np.load(output/'binning_comparison_covariance.npz')
    summary,individual=[],[]
    for scheme in ['xb_only','tprime_only']:
        rows=[r for r in payload['bins'] if r['scheme']==scheme]
        groups=len(rows)//3
        y=np.array([r['dilution_value'] for r in rows])
        c=covariance[scheme][:3*groups,:3*groups]
        design=np.tile(np.eye(groups),(3,1))
        precision=np.linalg.solve(c,np.eye(len(y)))
        mean_cov=np.linalg.inv(design.T@precision@design)
        mean=mean_cov@design.T@precision@y
        residual=y-design@mean
        statistic=float(residual@precision@residual)
        ndf=2*groups
        summary.append(dict(scheme=scheme,chi2=statistic,ndf=ndf,p_value=float(chi2.sf(statistic,ndf))))
        for j in range(groups):
            ids=np.array([j,groups+j,2*groups+j])
            values=y[ids]
            sub=c[np.ix_(ids,ids)]
            weights=np.linalg.solve(sub,np.ones(3))
            fitted=float(weights@values/weights.sum())
            res=values-fitted
            statistic=float(res@np.linalg.solve(sub,res))
            row=dict(scheme=scheme,group=j+1,xB_low=rows[j]['xB_low'],xB_high=rows[j]['xB_high'],
                     minus_tprime_low=rows[j]['minus_tprime_low'],minus_tprime_high=rows[j]['minus_tprime_high'],
                     common_df=fitted,common_df_stat=float(np.sqrt(1/weights.sum())),
                     chi2=statistic,ndf=2,p_value=float(chi2.sf(statistic,2)))
            for a,b in [(0,1),(0,2),(1,2)]:
                difference=float(values[a]-values[b])
                error=float(np.sqrt(sub[a,a]+sub[b,b]-2*sub[a,b]))
                row[f'{LABELS[a]}_minus_{LABELS[b]}']=difference
                row[f'{LABELS[a]}_minus_{LABELS[b]}_stat']=error
            individual.append(row)
    write_csv(output/'period_consistency_summary.csv',summary)
    write_csv(output/'period_consistency_bins.csv',individual)
    (output/'period_consistency.json').write_text(json.dumps(dict(cutset_id=payload['cutset_id'],
        method='Correlated GLS, using statistical Poisson-replica covariance with shared Sp23 controls',
        caveat='Approximate statistical-only tests. The projections use the same events and are not independent tests. Not a test of every two-dimensional bin.',
        summary=summary,bins=individual),indent=2,allow_nan=False)+'\n')
    for row in summary:
        print(f"Period consistency {row['scheme']}: chi2/ndf={row['chi2']:.2f}/{row['ndf']}, p={row['p_value']:.3f}")
    return summary,individual

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
    xb_comparison(out, counts, charges, args.replicas, args.seed, registry['cutset_id'])
    binning_comparison(out, counts, charges, args.replicas, args.seed, registry['cutset_id'])
    period_consistency(out)
    flagged = sum(bool(r['flags']) for r in rows)
    print(f'Finished: 72 period/bin estimates, {flagged} flagged for review. Outputs: {out.resolve()}')
    for row in packing_summary:
        print(f"{row['period']}: integrated cell packing fraction {row['value']:.3f}, central 68% interval [{row['stat_q16']:.3f}, {row['stat_q84']:.3f}]")
    print('Statistical intervals only. No scale uncertainty or radiative correction has been assigned.')


if __name__ == '__main__':
    main()
