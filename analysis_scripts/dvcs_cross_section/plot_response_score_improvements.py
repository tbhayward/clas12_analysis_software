#!/usr/bin/env python3
"""Restore fit-level and projection-level shape improvement bars from saved CSVs."""
import argparse,csv
from pathlib import Path
import numpy as np


def read(path):
    with path.open() as f:return list(csv.DictReader(f))


def improvement(a,b):
    return 100*(1-b/a) if np.isfinite(a) and np.isfinite(b) and a>0 else float('nan')


def render(folder):
    import matplotlib.pyplot as plt
    scores=read(folder/'final_scores.csv');lookup={r['scenario']:r for r in scores};a=lookup['nominal'];b=lookup['candidate'];rows=[];groups=[(k,v) for k,v in [('', 'Combined'),('epgamma_','epγ'),('epgammagamma_','epγγ')] if k+'train' in a];metrics=[('train','Training'),('validation_mc','Other MC half'),('heldout_observables','Reserved observables')]
    fig,axes=plt.subplots(1,len(groups),figsize=(13,5),layout='constrained')
    for ax,(prefix,title) in zip(np.atleast_1d(axes),groups):
        values=[]
        for key,label in metrics:
            field=prefix+key
            if field not in a:values.append(float('nan'));continue
            # Combined training uses shape score, excluding the nuisance penalty.
            if not prefix and key=='train' and 'train_shape' in a:field='train_shape'
            old=float(a[field]);new=float(b[field]);value=improvement(old,new);values.append(value);rows.append(dict(kind='overall',channel=title,topology='all fitted',observable=label,nominal=old,candidate=new,improvement_percent=value))
        bars=ax.bar(range(3),values,color=['#205cab','#22946b','#d45c18']);ax.set(xticks=range(3),xticklabels=[m[1] for m in metrics],ylabel='Reduction in shape score (%)',title=title,ylim=(min(-5,np.nanmin(values)-5),105));ax.tick_params(axis='x',rotation=20);ax.axhline(0,color='black',lw=.7);ax.grid(axis='y',alpha=.2)
        for bar,value in zip(bars,values):
            if np.isfinite(value):ax.text(bar.get_x()+bar.get_width()/2,value+1,f'{value:.1f}%',ha='center',fontsize=9)
    fig.suptitle(folder.name+' · nominal → final candidate; positive means improvement');fig.savefig(folder/'shape_score_improvement.png',dpi=120);plt.close(fig)
    projections=read(folder/'final_projection_scores.csv');channels=list(dict.fromkeys(r.get('channel','epgamma') for r in projections));tops=list(dict.fromkeys(r['topology'] for r in projections));fig,axes=plt.subplots(len(channels),len(tops),figsize=(5*len(tops),4*len(channels)),squeeze=False,layout='constrained')
    for ci,channel in enumerate(channels):
        for ti,top in enumerate(tops):
            selected=[r for r in projections if r.get('channel','epgamma')==channel and r['topology']==top];variables=list(dict.fromkeys(r['variable'] for r in selected));values=[];colors=[]
            for var in variables:
                rr=[r for r in selected if r['variable']==var];reserved=all(r['kind']=='heldout_observables' for r in rr);rr=[r for r in rr if reserved or r['kind']=='train'];scenarios=sorted({int(r['scenario']) for r in rr});old=np.mean([float(r['score']) for r in rr if int(r['scenario'])==scenarios[0]]);new=np.mean([float(r['score']) for r in rr if int(r['scenario'])==scenarios[-1]]);value=improvement(old,new);values.append(value);colors.append('#d45c18' if reserved else '#205cab');rows.append(dict(kind='projection',channel=channel,topology=top,observable=var,nominal=old,candidate=new,improvement_percent=value))
            ax=axes[ci,ti];ax.barh(variables,values,color=colors);ax.axvline(0,color='black',lw=.8);ax.set(title=channel+' · '+top,xlabel='Reduction in recorded projection score (%)');ax.grid(axis='x',alpha=.2)
    fig.suptitle(folder.name+' · individual shape comparisons\nBlue: training observables; orange: reserved; recorded profiled projection baselines');fig.savefig(folder/'observable_score_improvement.png',dpi=120);plt.close(fig)
    with (folder/'score_improvements.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report=folder/'REPORT.html'
    if report.is_file():
        marker='<!-- shape-score-bars -->'
        if marker not in report.read_text():
            with report.open('a') as f:f.write(marker+'<p><a href="shape_score_improvement.png">Overall shape-score improvement</a> · <a href="observable_score_improvement.png">Improvement by channel/topology/observable</a> · <a href="score_improvements.csv">Values</a></p><p>These bars measure descriptive score reduction, not percentage agreement or a chi-square probability. The overall nominal score uses its fixed reference convention; projection bars use the recorded profiled nominal/candidate comparisons and need not average to the overall bar. MC-half validation shares experimental DATA with training.</p>')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);args=p.parse_args();folders=[args.run] if (args.run/'final_scores.csv').is_file() else sorted(x for x in args.run.iterdir() if (x/'final_scores.csv').is_file())
    if not folders:p.error('No completed fit score CSVs found')
    for folder in folders:render(folder);print('[done]',folder/'shape_score_improvement.png')

if __name__=='__main__':main()
