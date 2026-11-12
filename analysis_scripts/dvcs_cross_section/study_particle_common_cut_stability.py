#!/usr/bin/env python3
"""Cache-only common-window shape and cut-stability checks; no response refitting."""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse, copy, hashlib, html, json
from pathlib import Path
import numpy as np
import study_particle_selection_closure as closure
base=closure.base
cuts=closure.cuts
SCENARIOS=('nominal','current','add_pt')
LABELS={'nominal':'Unmodified MC','current':'Original-training response','add_pt':'Add-pT response'}
COLORS={'nominal':'#888888','current':'#d45c18','add_pt':'#22946b'}
VARIABLES=('pTmiss','Mx2')
BINS=((.11,.15),(.15,.25),(.25,.40),(.40,.60),(.60,.80),(.80,1.))


def export_digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):value.update(chunk)
    return value.hexdigest()


def common_blocks(blocks):
    result=copy.deepcopy(blocks)
    for block in result.values():
        block['mc']=copy.deepcopy(block['data'])
    return result


def scaled_shape(data,signal,background,beta,alpha,prob,equal_area=False):
    s=signal*beta/prob['mc'];b=background*alpha/prob['background']
    denominator=(s+b).sum() if equal_area else data.sum()
    if not np.isfinite(denominator) or denominator<=0:
        return np.full_like(s,np.nan),np.full_like(b,np.nan)
    return s/denominator,b/denominator


def scan_counts(values,selected,thresholds,center,reference):
    """Sort once; nested cut counts/overlaps are exact, including boundary events."""
    x=np.sort(np.abs(values[selected]-center));counts=np.searchsorted(x,thresholds,side='right')
    nref=int(np.searchsorted(x,reference,side='right'))
    return counts,nref,np.minimum(counts,nref)


def stability(counts,reference,overlap,alpha,prob,alpha_variance=0.):
    """Paired delta-method errors; control scale independent of one-photon events."""
    d,b,s=[np.asarray(x,dtype=float) for x in counts]
    d0,b0,s0=np.asarray(reference,dtype=float);od,ob,os=overlap
    wb=alpha/prob['background'];ws=1/prob['mc']
    y=d-wb*b;y0=d0-wb*b0;m=ws*s;m0=ws*s0
    valid=(y>0)&(m>0)&(y0>0)&(m0>0)&np.isfinite(alpha)
    with np.errstate(divide='ignore',invalid='ignore'):
        ratio=y/m;relative=ratio/(y0/m0)
        variance=(d+wb*wb*b)/y**2+(d0+wb*wb*b0)/y0**2-2*(od+wb*wb*ob)/(y*y0)
        variance+=1/s+1/s0-2*os/(s*s0)
        # The same control normalization is used at every threshold.
        variance+=alpha_variance*(-b/prob['background']/y+b0/prob['background']/y0)**2
        error=np.abs(relative)*np.sqrt(np.maximum(0,variance))
    return y,np.where(valid,ratio,np.nan),np.where(valid,relative,np.nan),np.where(valid,error,np.nan)


def calculate_scan(arrays,top,block,variable,bounds,thresholds,reference,center,scale,prob):
    col=1+base.OBS.index(variable);counts=[];refs=[];overlap=[]
    for sample in ('data','background','mc'):
        a=arrays[sample]
        selected=closure.mask(a,top,block,'data',omit=variable)&(a[:,9]>=bounds[0])&(a[:,9]<bounds[1])
        c,r,o=scan_counts(a[:,col],selected,thresholds,center,reference)
        counts.append(c);refs.append(r);overlap.append(o)
    alpha=scale['background_scale'];control_d=scale['control_data_anchor'];control_m=scale['control_mc_anchor']*prob['pi0_mc']
    va=alpha**2*(1/control_d+1/control_m) if control_d>0 and control_m>0 else float('nan')
    y,r,relative,error=stability(counts,refs,overlap,alpha,prob,va)
    return [dict(threshold=float(t),data_count=int(counts[0][i]),background_mc_count=int(counts[1][i]),signal_mc_count=int(counts[2][i]),background_estimate=float(alpha*counts[1][i]/prob['background']),subtracted_data=float(y[i]),signal_scale=float(r[i]),relative_signal_scale=float(relative[i]),conditional_stat_error=float(error[i]),status='ok' if np.isfinite(relative[i]) else 'nonpositive residual, empty MC or invalid control',reference_threshold=reference,window_center=center) for i,t in enumerate(thresholds)]


def render_case(period,top,bounds,cases,blocks,families,anchors,prob,out):
    import matplotlib.pyplot as plt
    cases={c:{s:a[(a[:,9]>=bounds[0])&(a[:,9]<bounds[1])] for s,a in samples.items()} for c,samples in cases.items()}
    topology=base.TOPS[top];cfg=next(p for p in cuts.production.PERIODS if p.key==period)
    key=cuts.production.cut_json_key('dvcs',cfg,cuts.production.TOPOLOGIES[top])
    tag=f'{period}_{topology}_t{bounds[0]:.2f}_{bounds[1]:.2f}';links=[];records=[]
    fig,axes=plt.subplots(3,2,figsize=(12,11),layout='constrained')
    # Same saved DATA windows on every sample; omit the displayed cut only.
    for vi,var in enumerate(VARIABLES):
        edges=np.linspace(0,.45,121) if var=='pTmiss' else np.linspace(-.045,.045,151)
        col=1+base.OBS.index(var);hist={};pre={}
        for scenario in SCENARIOS:
            a=cases[scenario];h=[]
            for sample in ('data','mc','background'):
                rows=a[sample];sel=closure.mask(rows,top,blocks[key],'data',omit=var)&(rows[:,9]>=bounds[0])&(rows[:,9]<bounds[1])
                # Under/overflow included in area normalization, even when outside the plotted axis.
                h.append(closure.histogram(rows[:,col],sel,edges))
            hist[scenario]=h
            pre[scenario]=[closure.histogram(a[sample][:,col],np.ones(len(a[sample]),dtype=bool),edges) for sample in ('data','mc','background')]
        for row,(hh,equal) in enumerate(((pre,True),(hist,True),(hist,False))):
            d=hh['add_pt'][0];dn=d.sum()
            ax=axes[row,vi]
            if dn:ax.stairs(d[1:-1]/dn/np.diff(edges),edges,color='#205cab',label='DATA')
            for scenario in SCENARIOS:
                _,s,b=hh[scenario];sc=anchors[scenario][top]
                ss,bb=scaled_shape(d,s,b,sc['signal_scale'],sc['background_scale'],prob,equal)
                ax.stairs((ss+bb)[1:-1]/np.diff(edges),edges,color=COLORS[scenario],label=LABELS[scenario])
                if scenario=='add_pt':
                    ax.stairs(ss[1:-1]/np.diff(edges),edges,color=COLORS[scenario],ls='--',lw=.8,label='Add-pT DVCS')
                    ax.stairs(bb[1:-1]/np.diff(edges),edges,color=COLORS[scenario],ls=':',lw=.8,label='Add-pT π⁰')
            ax.set(title=('Before exclusivity: equal areas' if row==0 else 'Common N−1 cuts: equal areas' if equal else 'Common N−1 cuts: fixed scales'),xlabel=var,ylabel='Density');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle(tag+' · full globals; before exclusivity and common N−1 windows')
    path=out/'plots'/(tag+'_shapes.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for fi,(family,allblocks) in enumerate(families.items()):
        block=allblocks[key]
        for vi,var in enumerate(VARIABLES):
            stats=block['data'].get(var);center=0. if var=='pTmiss' else (float(stats['cut_high'])+float(stats['cut_low']))/2 if stats else 0.
            reference=float(stats['cut_high'])-center if stats else (.2 if var=='pTmiss' else .02)
            thresholds=np.unique(np.r_[np.linspace(.04,.42,25) if var=='pTmiss' else np.linspace(.002,.045,25),reference])
            ax=axes[fi,vi]
            for scenario in SCENARIOS:
                rows=calculate_scan(cases[scenario],top,block,var,bounds,thresholds,reference,center,anchors[scenario][top],prob)
                for r in rows:r.update(period=period,topology=topology,t_low=bounds[0],t_high=bounds[1],scenario=scenario,window_family=family,variable=var)
                records+=rows
                ax.errorbar(thresholds,[r['relative_signal_scale'] for r in rows],yerr=[r['conditional_stat_error'] for r in rows],color=COLORS[scenario],label=LABELS[scenario],lw=1,elinewidth=.5)
            ax.axhline(1,color='black',lw=.8);ax.axvline(reference,color='black',ls=':',lw=.8)
            ax.set(title=family,xlabel=var+(' upper cut (GeV)' if var=='pTmiss' else ' half-window (GeV²)'),ylabel='Signal scale / scale at reference');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle(tag+' · identical cuts on DATA, signal MC and leakage MC')
    path=out/'plots'/(tag+'_stability.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    return records,links


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('followup-study','closure-study','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--periods',nargs='+',choices=[p.key for p in cuts.production.PERIODS])
    parser.add_argument('--topologies',nargs='+',choices=base.TOPS,default=base.TOPS)
    args=parser.parse_args(argv);source=args.followup_study.resolve()/'cut_redetermination';closed=args.closure_study.resolve();out=args.output.resolve()
    if out.exists():parser.error('Output exists; choose a fresh output directory')
    old=json.loads((source/'manifest.json').read_text());cmf=json.loads((closed/'manifest.json').read_text())
    if cmf['source_cut_manifest_sha256']!=base.digest(source/'manifest.json') or cmf['candidate_hashes']!=old['candidate_hashes']:parser.error('Closure and follow-up provenance differ')
    spec=json.loads((source/'spec.json').read_text());periods=args.periods or spec['periods']
    if not set(periods)<=set(spec['periods']):parser.error('Period missing from source study')
    run=Path(old['arguments']['response_run']);training=Path(old['arguments']['training_study'])
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    cachepath=run/'cache_manifest.json';cache=json.loads(cachepath.read_text())
    if base.digest(cachepath)!=old['response_cache_sha256']:parser.error('Response-cache metadata changed')
    paths=cuts.candidate_paths({c:run for c in closure.MODELS},periods,training)
    for (p,c),path in paths.items():
        if base.digest(path)!=old['candidate_hashes'][c][p]:parser.error('Response parameters changed')
    vectors={k:cuts.load_candidate_parameters(v) for k,v in paths.items()}
    try:exports,files=closure.load_exports(source,periods,vectors)
    except FileNotFoundError as e:parser.error(str(e)+'; retain/rebuild the follow-up pre-cut scratch exports before this cache-only study')
    # Verify the exact arrays used in the closure study; hashing is sequential and memory bounded.
    print('[verify] Checking streamed export hashes once; no ROOT reads or response refits',flush=True)
    for path in set(files.values()):
        expected=cmf['export_hashes'].get(str(path))
        if expected is None or export_digest(path)!=expected:parser.error('Pre-cut export differs from closure inputs: '+str(path))
    out.mkdir(parents=True);(out/'plots').mkdir();records=[];links=[];norms=[]
    saved=common_blocks(spec['exclusivity'])
    base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},source_sha256=base.digest(source/'manifest.json'),closure_sha256=base.digest(closed/'manifest.json'),candidate_hashes=old['candidate_hashes'],normalization='Common saved DATA windows for all samples, including control; scales fixed across thresholds',caveats=['Unit-weight counts; no final cross sections','Statistical errors conditional on fitted response/windows and transfer multiplier','Only |t| subdivision: exports lack xB/Q2/polar angles; changing kinematic mixtures can cause scale drift']))
    for period in periods:
        cases={c:closure.arrays_for(exports,period,c) for c in SCENARIOS};prob={s:cache['samples'][period+'__'+s]['probability'] for s in ('mc','background','pi0_mc')};anchors={}
        for c in SCENARIOS:
            multiplier=1. if c=='nominal' else json.loads(paths[period,c].read_text())['background_multiplier']
            anchors[c]=closure.scales(period,cases[c],saved,prob,spec.get('control_to_leakage_exposure_ratio',1.),multiplier,True)
            for top,v in anchors[c].items():norms.append(dict(period=period,scenario=c,topology=base.TOPS[top],**v))
        refit=common_blocks(json.loads((closed/'global_cut_refits'/period/'add_pt'/'combined_cuts_95.json').read_text()))
        for topology in args.topologies:
            top=base.TOPS.index(topology)
            local={c:{s:a[(a[:,0]==top)&(a[:,8]>.5)] for s,a in samples.items() if s in ('data','mc','background')} for c,samples in cases.items()}
            for bounds in BINS:
                rows,ll=render_case(period,top,bounds,local,saved,dict(saved_DATA=saved,refitted_DATA=refit),anchors,prob,out);records+=rows;links+=ll
        print('[done] '+period+' common-cut scans; no ROOT reads, compilation or parameter fits',flush=True)
    base.table(out/'cut_stability.csv',records);base.table(out/'normalizations.csv',norms)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Common-cut stability</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}li{margin:8px}</style><h1>Common-cut stability</h1><p>Exactly the same boundaries are applied to DATA, particle-corrected DVCS MC and π⁰ leakage MC. No additional histogram shift or broadening is fitted. Every scenario uses its own consistent two-photon MC, with common saved DATA control windows. The period-wide control scale fixes π⁰, then one signal scale per topology matches the selected DATA total.</p><p>Shape plots: first row shows full globals before exclusivity, second row applies common windows but omits the displayed variable. Both rescale the summed MC to the DATA integral (components retain relative weights). Third row carries the fixed selected-anchor scale under the second-row selection. Before integrating, overflow is retained. No cross sections are compared.</p><p>Stability: (DATA − predicted π⁰)/DVCS MC, divided by its value at the dotted reference threshold in the same |t| bin and scenario. A flat curve supports cut stability. All other cuts remain fixed. The two rows use saved versus newly refitted DATA windows on every sample; a 95% window is not imposed or re-fitted here. Relative scales do not independently validate absolute efficiencies. The reference value is one by definition.</p><p>Error bars include paired, overlapping counts and the shared control-count uncertainty, conditional on the fitted response, cut boundaries and background transfer. They exclude transfer/model uncertainties. Counts are unit weighted. These exports cannot subdivide xB, Q² or particle angles: differing kinematic populations can also cause drift, so do not interpret every slope as a resolution defect.</p><p><a href="cut_stability.csv">Threshold counts, subtracted yields and scales</a> · <a href="normalizations.csv">Fixed signal/control anchors and validity flags</a></p><ul>'+''.join('<li><a href="'+html.escape(p)+'">'+html.escape(Path(p).stem)+'</a></li>' for p in links)+'</ul>')
    print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
