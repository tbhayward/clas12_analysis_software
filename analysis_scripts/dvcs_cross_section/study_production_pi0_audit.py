#!/usr/bin/env python3
"""Replay production pi0 algebra and plot direct nominal/morphed control |t| counts."""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
PERIODS=['Fa18 Inb','Fa18 Out','Sp18 Inb','Sp18 Out','Sp19 Inb'];TOPS=['(FD, FD)','(CD, FD)','(CD, FT)'];EDGES=['xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax']


def triple(s):
    if not s.strip():return np.zeros(3)
    a=np.fromstring(s.strip().strip('()'),sep=',')
    if len(a)!=3 or not np.all(np.isfinite(a)):raise ValueError('Invalid yield triple: '+s)
    return a


def sums(row,channel,period,data):
    prefix='normalized raw yield' if data else 'reconstructed current corrected yield'
    arrays=[triple(row[f'{prefix}, {channel}, {top}, '+('exp' if data else 'mc')+f', {period}'+(', unpol' if data else '')]) for top in TOPS]
    return np.array([sum(a[0] for a in arrays),np.sqrt(sum(a[1]**2 for a in arrays)),np.sqrt(sum(a[2]**2 for a in arrays))])


def calculate(mis,pi0data,pi0mc,dvcs):
    values=np.array([mis,pi0data,pi0mc,dvcs])
    if np.any(values[:,0]<=0):return np.zeros(3)
    c=mis[0]*pi0data[0]/(pi0mc[0]*dvcs[0]);return np.array([c,abs(c)*np.sqrt(np.sum((values[:,1]/values[:,0])**2)),0.])


def replay(path):
    rows=list(csv.DictReader(Path(path).open()));result=[]
    for i,row in enumerate(rows):
        for period in PERIODS:
            name='contamination ratio, '+period
            if not row.get(name,'').strip():continue
            d=sums(row,'ep->epg',period,True);pd=sums(row,'ep->eppi0',period,True);pm=sums(row,'ep->eppi0',period,False);b=sums(row,'ep->eppi0->epg',period,False)
            c=calculate(b,pd,pm,d);saved=triple(row[name]);result.append(dict(row=i,period=period,**{k:row[k] for k in EDGES},data_dvcs=d[0],data_pi0=pd[0],mc_pi0=pm[0],mc_leakage=b[0],saved_fraction=saved[0],replayed_fraction=c[0],saved_stat=saved[1],replayed_stat=c[1],difference=c[0]-saved[0],stat_difference=c[1]-saved[1],status='nonpositive input: production returns zero' if min(d[0],pd[0],pm[0],b[0])<=0 else 'positive inputs'))
    return result


def table(path,rows):
    if not rows:return
    with Path(path).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def tplots(validation,out):
    import matplotlib.pyplot as plt
    ratios=list(csv.DictReader((validation/'heldout_ratios.csv').open()));totals=list(csv.DictReader((validation/'heldout_totals.csv').open()));links=[]
    scales={(r['period'],r['scenario']):float(r['test_prediction'])/float(r['test_raw_mc']) for r in totals if r['scenario'] in ('nominal','add_pt') and int(r['test_raw_mc'])>0}
    counts=[]
    for period in dict.fromkeys(r['period'] for r in totals):
        for top in ['all_control_topologies','FD_FD','CD_FD','CD_FT']:
            fig,axes=plt.subplots(3,1,figsize=(9,10),layout='constrained',sharex=True);hist={}
            for scenario,color in [('nominal','#777777'),('add_pt','#22946b')]:
                rr=[r for r in ratios if r['period']==period and r['topology']==top and r['axis']=='t_abs' and r['scenario']==scenario];rr.sort(key=lambda r:float(r['low']));edges=np.array([float(r['low']) for r in rr]+[float(rr[-1]['high'])]);raw=np.array([float(r['test_prediction'])/scales[period,scenario] for r in rr]);hist[scenario]=raw
                for r,n in zip(rr,raw):counts.append(dict(period=period,topology=top,scenario=scenario,t_low=r['low'],t_high=r['high'],sampled_mc_selected=n))
                axes[0].stairs(raw,edges,color=color,label=scenario);axes[1].stairs(raw/raw.sum()/np.diff(edges) if raw.sum() else raw*np.nan,edges,color=color,label=scenario)
            with np.errstate(divide='ignore',invalid='ignore'):change=hist['add_pt']/hist['nominal']
            axes[2].stairs(change,edges,color='#22946b');axes[2].axhline(1,color='black',lw=.8)
            axes[0].set_ylabel('Selected sampled MC events');axes[1].set_ylabel('Unit-area density');axes[2].set(ylabel='Corrected / nominal raw counts',xlabel='|t| (GeV²)')
            for ax in axes:ax.grid(alpha=.2)
            axes[0].legend();axes[1].legend();fig.suptitle(period+' · '+top+' · two-photon control MC\nSame cached sample/test half; common DATA cuts; no DATA normalization')
            folder=out/'direct_t'/period;folder.mkdir(parents=True,exist_ok=True);path=folder/(top+'.png');fig.savefig(path,dpi=120);plt.close(fig);links.append(path.relative_to(out).as_posix())
    table(out/'direct_t_counts.csv',counts);return links


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--production-csv',type=Path,required=True);p.add_argument('--validation-study',type=Path,required=True);p.add_argument('--morphed-production-csv',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args(argv)
    if args.output.exists():p.error('Output exists; choose a fresh path')
    args.output.mkdir(parents=True);rows=replay(args.production_csv);table(args.output/'baseline_replay.csv',rows);summary=[]
    for period in PERIODS:
        r=[r for r in rows if r['period']==period];summary.append(dict(period=period,rows=len(r),positive_input_rows=sum(x['status']=='positive inputs' for x in r),max_fraction_difference=max((abs(x['difference']) for x in r),default=0),max_stat_difference=max((abs(x['stat_difference']) for x in r),default=0)))
    table(args.output/'replay_summary.csv',summary)
    if args.morphed_production_csv:
        morphed=replay(args.morphed_production_csv);lookup={(r['period'],tuple(r[k] for k in EDGES)):r for r in morphed};comparison=[]
        for r in rows:
            m=lookup.get((r['period'],tuple(r[k] for k in EDGES)))
            if m is None:raise ValueError('Morphed CSV has different bins/periods')
            if not np.allclose([r['data_dvcs'],r['data_pi0']],[m['data_dvcs'],m['data_pi0']],rtol=1e-10,atol=1e-10):raise ValueError('DATA yields changed; comparison requires fixed original DATA selection')
            comparison.append(dict(period=r['period'],**{k:r[k] for k in EDGES},nominal_fraction=r['replayed_fraction'],morphed_fraction=m['replayed_fraction'],change_percentage_points=100*(m['replayed_fraction']-r['replayed_fraction'])))
        table(args.output/'exact_morph_comparison.csv',comparison)
    links=tplots(args.validation_study,args.output)
    note='Exact morphed-production CSV comparison supplied.' if args.morphed_production_csv else 'No exact morphed production yields were supplied. Existing diagnostic exports omit phi/current-weight information needed for that calculation. These control |t| plots are not a morphed production-contamination comparison.'
    (args.output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><h1>Production pi0 replay and direct |t| comparison</h1><p>Production formula: topology-summed current-corrected leakage MC × normalized pi0 DATA / (topology-summed current-corrected pi0 MC × normalized DVCS DATA), separately in every xB,Q²,|t|,phi row. No period-wide control replacement is used. Statistical propagation reproduces the production formula; this checks arithmetic, not upstream efficiencies/exposures.</p><p>'+note+'</p><p>Direct plots reconstruct raw selected sampled test-MC counts by removing each scenario\'s recorded constant normalization. They compare the same cache/test-half and common DATA windows, but each scenario\'s own reconstructed |t|. They combine migration and selection changes, not only migration. The middle row normalizes only over the plotted 0.11–1 range. No weighted scenario is shown.</p><p><a href="replay_summary.csv">Replay summary</a> · <a href="baseline_replay.csv">Per-production-bin replay</a> · <a href="direct_t_counts.csv">Raw selected MC counts</a></p>'+''.join('<p><a href="'+l+'">'+l+'</a></p>' for l in links))
    print(json.dumps(summary,indent=2));print('[done]',args.output/'REPORT.html');print(note)

if __name__=='__main__':main()
