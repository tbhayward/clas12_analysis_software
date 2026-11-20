#!/usr/bin/env python3
"""Compare production fits with nominal-toy mean-bias-subtracted diagnostics (PNG only)."""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from plot_final_results import (OBSERVABLES, PANEL_INDICES, LABELS, TP_EDGES,
                                XB_EDGES, XB_COLORS)

HERE = Path(__file__).resolve().parent


def calculate_biases(frame, validation_dir):
    """Use nominal injections only, and the validation's saved usable-fit selection."""
    frame = frame.sort_values('bin_number').reset_index(drop=True)
    if frame.bin_number.tolist() != list(range(1,25)):
        raise ValueError('Expected each production bin 1..24 exactly once')
    config = json.loads((validation_dir/'manifest.json').read_text())['config']
    requested = int(config['replicas'])
    table, hashes = [], {}
    result = frame.copy()
    for i, measured in frame.iterrows():
        b = int(measured.bin_number)
        path = validation_dir/'checkpoints'/f'bin_{b:02d}.json'
        payload = path.read_bytes()
        hashes[str(b)] = hashlib.sha256(payload).hexdigest()
        saved = json.loads(payload)
        if saved['bin_number'] != b:
            raise ValueError(f'Checkpoint bin mismatch for {path}')
        scenario = saved.get('scenarios',{}).get('nominal',{})
        rows = scenario.get('rows',[])
        if len(rows) != requested:
            raise ValueError(f'Bin {b}: nominal injection has {len(rows)}/{requested} toys')
        if sorted(row['replica'] for row in rows) != list(range(requested)):
            raise ValueError(f'Bin {b}: duplicate or missing nominal replicas')
        good = [row for row in rows if row.get('usable')]
        if len(good) < 2:
            raise ValueError(f'Bin {b}: fewer than two usable toys')
        for obs in OBSERVABLES:
            truth = float(scenario['truth'][obs])
            value, stat = float(measured[obs]), float(measured[obs+'_stat'])
            if not np.isclose(value,truth,rtol=1e-8,atol=1e-10):
                raise ValueError(f'Bin {b}, {obs}: injected nominal truth differs from production result')
            estimates = np.array([row['fit']['values'][obs] for row in good],dtype=float)
            if not np.isfinite(estimates).all() or not np.isfinite(stat) or stat <= 0:
                raise ValueError(f'Invalid values in bin {b}, {obs}')
            mean = float(estimates.mean())
            spread = float(estimates.std(ddof=1))
            bias = mean-truth
            mean_se = spread/np.sqrt(len(good))
            adjusted = value-bias
            result.loc[i,obs+'_nominal_injection_mean_bias'] = bias
            result.loc[i,obs+'_nominal_injection_bias_mc_error'] = mean_se
            result.loc[i,obs+'_bias_subtracted_diagnostic'] = adjusted
            table.append(dict(bin_number=b, x_index=int(measured.x_index),
                t_index=int(measured.t_index), parameter=obs,
                mean_minus_tprime_gev2=measured.mean_minus_tprime_gev2,
                nominal=value, nominal_stat=stat, injected_truth=truth,
                toy_mean=mean, toy_spread=spread, mean_bias=bias,
                mean_bias_mc_error=mean_se, mean_bias_over_nominal_stat=bias/stat,
                bias_subtracted_diagnostic=adjusted, applied_shift=-bias,
                attempted=len(rows), usable=len(good), unusable=len(rows)-len(good)))
    return result,pd.DataFrame(table),hashes


def plot_comparisons(frame, output_dir):
    paths=[]
    for xi in [0,1,2,3,None]:
        combined = xi is None
        fig,axes=plt.subplots(3,3,figsize=(12,10));axes=axes.ravel()
        for ax,obs,label in zip(axes[list(PANEL_INDICES)],OBSERVABLES,LABELS):
            for x_index in range(4) if combined else [xi]:
                group=frame[frame.x_index==x_index].sort_values('t_index')
                if group.t_index.tolist()!=list(range(6)):
                    raise ValueError(f'Incomplete t-prime group {x_index}')
                x=group.mean_minus_tprime_gev2.to_numpy()
                y=group[obs].to_numpy();z=group[obs+'_bias_subtracted_diagnostic'].to_numpy()
                stat=group[obs+'_stat'].to_numpy()
                if not combined or x_index==1:
                    ax.bar(TP_EDGES[:-1],group[obs+'_point_to_point_systematic'],
                           width=np.diff(TP_EDGES),align='edge',bottom=0,
                           color='0.65',alpha=.35,linewidth=0,zorder=1)
                original_color=XB_COLORS[x_index] if combined else 'black'
                adjusted_color=XB_COLORS[x_index] if combined else '#D55E00'
                # Equal and opposite small offsets expose overlapping error bars.
                ax.errorbar(x-.008,y,yerr=stat,fmt='o',color=original_color,
                            markersize=4,capsize=2,linewidth=.9,zorder=3)
                ax.errorbar(x+.008,z,yerr=stat,fmt='D',color=adjusted_color,
                            markerfacecolor='white',markeredgewidth=1.1,
                            markersize=4,capsize=2,linewidth=.9,zorder=4)
                for a,v,w in zip(x,y,z):
                    ax.plot([a-.008,a+.008],[v,w],color=adjusted_color,alpha=.5,lw=.8,zorder=2)
            ax.axhline(0,color='0.45',linewidth=.8,zorder=2)
            single=obs in ('lu1','ul1','ul2')
            ax.set_xlim(0,1.30);ax.set_ylim((-.4,.4) if single else (-1,1))
            ax.set_yticks(np.linspace(-.4,.4,9) if single else np.linspace(-1,1,11))
            ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)");ax.set_ylabel(label)
            ax.grid(alpha=.18,zorder=0)
        axes[2].axis('off');axes[8].axis('off')
        handles=[Line2D([],[],marker='o',color='black',linestyle='none',markersize=5,
                        label='Fitted value'),
                 Line2D([],[],marker='D',color='black' if combined else '#D55E00',
                        markerfacecolor='white',linestyle='none',markersize=5,
                        label='Nominal-toy mean-bias subtracted')]
        if combined:
            handles += [Line2D([],[],color=c,lw=2,label=rf'${XB_EDGES[j]:.2f}<x_B<{XB_EDGES[j+1]:.2f}$')
                        for j,c in enumerate(XB_COLORS)]
        handles.append(Patch(facecolor='0.65',alpha=.35,edgecolor='none',
                             label='Original systematic magnitude'+('\n(second xB group)' if combined else '')))
        axes[8].legend(handles=handles,loc='upper center',frameon=False,fontsize=8.5)
        axes[8].text(.05,.12 if combined else .30,
                     'Both use the original statistical error.\nHorizontal offsets are for visibility.\nDiagnostic; no production correction applied.',
                     transform=axes[8].transAxes,fontsize=8.5,va='top',linespacing=1.5)
        title='All '+r'$x_B$'+' ranges' if combined else rf'${XB_EDGES[xi]:.2f}<x_B<{XB_EDGES[xi+1]:.2f}$'
        fig.suptitle(title+' — nominal-injection bias diagnostic',fontsize=13)
        fig.tight_layout(rect=(0,0,1,.96))
        name='all_xB' if combined else f'xB_bin_{xi+1}'
        path=output_dir/f'injection_bias_comparison_{name}.png'
        fig.savefig(path,dpi=200);plt.close(fig);paths.append(path)
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-csv',type=Path,default=HERE/'output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv')
    parser.add_argument('--validation-dir',type=Path,default=HERE/'output/unpolarized_validation')
    parser.add_argument('--output-dir',type=Path,default=HERE/'output/injection_bias_comparisons')
    args=parser.parse_args();input_csv=args.input_csv.expanduser().resolve()
    validation=args.validation_dir.expanduser().resolve();output=args.output_dir.expanduser().resolve()
    if output==validation or output in validation.parents or validation in output.parents:
        raise ValueError('Use a diagnostic output directory separate from validation checkpoints')
    if output==input_csv.parent or output in input_csv.parents:
        raise ValueError('Do not write diagnostics over the production tables')
    frame,table,hashes=calculate_biases(pd.read_csv(input_csv),validation)
    output.mkdir(parents=True,exist_ok=True)
    table.to_csv(output/'nominal_injection_bias_by_modulation.csv',index=False)
    frame.to_csv(output/'structure_function_ratios_bias_diagnostic.csv',index=False)
    paths=plot_comparisons(frame,output)
    info={'input_csv':str(input_csv),'input_csv_sha256':hashlib.sha256(input_csv.read_bytes()).hexdigest(),
          'validation_dir':str(validation),'checkpoint_sha256':hashes,
          'scenario':'nominal','selection':'saved validation usable flag, common to all seven parameters',
          'definition':'bias=mean(toy_fit)-injected_truth; diagnostic=production_fit-bias',
          'displayed_errors':'original production statistical errors, unchanged',
          'mean_bias_mc_error':'sample standard deviation / sqrt(usable); tabulated, not added to plotted errors',
          'systematic_bands':'original magnitudes, unchanged; all-xB uses second group as representative',
          'interpretation':'plug-in mean-bias subtraction diagnostic, not a validated production correction',
          'plots':[str(p) for p in paths]}
    (output/'manifest.json').write_text(json.dumps(info,indent=2)+'\n')
    print(f'Wrote {len(paths)} PNGs and all {len(table)} nominal bias estimates to {output}')


if __name__=='__main__':main()
