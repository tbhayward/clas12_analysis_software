#!/usr/bin/env python3
"""Bin-by-bin implementation of Maggie's preliminary exclusivity procedure."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import sys
import time
import numpy as np

# Keep shared kinematics in the original diagnostics module.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'diagnostics'))
from kinematics import BRANCHES, XB, TP, calculate, bin_indices
from fit_models import MASS_EDGES, MOM_EDGES, ANGLE_EDGES, components, fit_with_checks, fit_angle

from candidate_cache import build_cache, iterate

DEFAULT_OUTPUT = Path(__file__).resolve().parent/'output'
PERIODS = ['RGB_spring2019_D2','RGB_fall2019_D2','RGB_spring2020_D2',
           'RGC_summer2022_ND3','RGC_fall2022_ND3','RGC_spring2023_ND3']

def save_json(path, value):
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)

def windows(values, bins, fits):
    lows = np.array([f['cut_low'] if f['usable'] else np.inf for f in fits])
    highs = np.array([f['cut_high'] if f['usable'] else -np.inf for f in fits])
    return (values >= lows[bins]) & (values <= highs[bins])

def histogram(hist, values, bins, edges):
    vi = np.searchsorted(edges, values, side='right')-1
    good = np.isfinite(values) & (vi >= 0) & (vi < len(edges)-1)
    hist += np.bincount(bins[good]*(len(edges)-1)+vi[good], minlength=hist.size).reshape(hist.shape)

def scan(cache,meta,stage,out,mass_fits=None,momentum_fits=None):
    edges=MASS_EDGES if stage=='mass' else MOM_EDGES if stage=='momentum' else ANGLE_EDGES
    hist=np.zeros((24,len(edges)-1),dtype=np.int64)
    counts={k:np.zeros(24,dtype=np.int64) for k in ['input','selected','FD','CD','underflow','overflow','undefined']}
    started=time.monotonic()
    read=0
    save_json(out/'progress.json',dict(stage=stage,state='reading derived cache',percent=0,updated=time.time()))
    for a in iterate(cache,meta):
        bins=a['bin'].astype(np.int64)
        read+=len(bins)
        counts['input']+=np.bincount(bins,minlength=24)
        selected=np.ones(len(bins),dtype=bool)
        if mass_fits is not None:selected &= windows(a['mass'],bins,mass_fits)
        if momentum_fits is not None:selected &= windows(a['momentum'],bins,momentum_fits)
        bs,vs=bins[selected],a[stage][selected]
        histogram(hist,vs,bs,edges)
        counts['selected']+=np.bincount(bs,minlength=24)
        for key,det in [('FD',1),('CD',2)]:counts[key]+=np.bincount(bins[selected&(a['detector']==det)],minlength=24)
        for key,mask in [('underflow',np.isfinite(vs)&(vs<edges[0])),('overflow',np.isfinite(vs)&(vs>=edges[-1])),('undefined',~np.isfinite(vs))]:
            counts[key]+=np.bincount(bs[mask],minlength=24)
        save_json(out/'progress.json',dict(stage=stage,state='reading derived cache',percent=100*read/max(1,meta['combinations']),updated=time.time()))
    print(f'{out.name}: {stage} histogram complete from cache',flush=True)
    return hist,counts,dict(read=read,elapsed_seconds=time.monotonic()-started,source='derived candidate cache')

def write_table(out,mass,momentum,counts):
    fields=['mu','sigma','mu_error','sigma_error','cut_low','cut_high','chi2_ndf','max_boundary_shift','stability_threshold']
    with (out/'cuts.csv').open('w',newline='') as stream:
        w=csv.writer(stream)
        w.writerow(['bin','xB_low','xB_high','minus_tprime_low','minus_tprime_high','stage','model','usable',*fields,'input_combinations','after_mass','after_both','FD_after_both','CD_after_both','warnings','reason'])
        for b in range(24):
            for stage,fits in [('mass',mass),('momentum',momentum)]:
                f=fits[b]
                w.writerow([b+1,XB[b//6],XB[b//6+1],TP[b%6],TP[b%6+1],stage,f['model'],f['usable'],*[f.get(k,'') for k in fields],counts['mass']['input'][b],counts['momentum']['selected'][b],counts['angle']['selected'][b],counts['angle']['FD'][b],counts['angle']['CD'][b],' | '.join(f['warnings']),f['reason']])
    with (out/'stability.csv').open('w',newline='') as stream:
        w=csv.writer(stream)
        w.writerow(['bin','stage','variant','usable','mu','sigma','cut_low','cut_high','chi2_ndf','low_shift','high_shift','reason'])
        for b in range(24):
            for stage,fits in [('mass',mass),('momentum',momentum)]:
                f=fits[b]
                for variant,v in [('nominal',f),*f['checks'].items()]:
                    low=v.get('cut_low',0)-f.get('cut_low',0) if v['usable'] and f['usable'] else ''
                    high=v.get('cut_high',0)-f.get('cut_high',0) if v['usable'] and f['usable'] else ''
                    w.writerow([b+1,stage,variant,v['usable'],*[v.get(k,'') for k in ['mu','sigma','cut_low','cut_high','chi2_ndf']],low,high,v['reason']])

def run_period(job):
    name, source, tree, output, options = job
    out = Path(output)/name
    out.mkdir(parents=True, exist_ok=True)
    source = Path(source)
    stat = source.stat()
    signature = dict(source=str(source.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                     max_entries=options['max_entries'], tree=tree,
                     code_sha256=hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name('fit_models.py').read_bytes()+Path(sys.modules['kinematics'].__file__).read_bytes()+Path(__file__).with_name('candidate_cache.py').read_bytes()).hexdigest())
    manifest = out/'manifest.json'
    if manifest.exists() and not options['force'] and not options.get('rebuild_cache',False):
        old = json.loads(manifest.read_text())
        if old['signature'] == signature and (out/'histograms.npz').exists():
            print(f'CACHED {name}', flush=True)
            plot_period(out)
            return name
    manifest.unlink(missing_ok=True)
    print(f'START {name}: one ROOT scan plus cached fits, PID {os.getpid()}', flush=True)
    save_json(out/'progress.json', dict(stage='mass',state='starting',updated=time.time()))
    cache,cache_meta=build_cache(source,tree,out,options)
    mass_hist, mc, mi = scan(cache,cache_meta,'mass',out)
    save_json(out/'progress.json',dict(stage='mass fits',state='fitting and testing stability',updated=time.time()))
    mass = [fit_with_checks(h,MASS_EDGES,'mass') for h in mass_hist]
    save_json(out/'mass_fits.json',mass)
    momentum_hist, pc, pi = scan(cache,cache_meta,'momentum',out,mass)
    save_json(out/'progress.json',dict(stage='momentum fits',state='fitting and testing stability',updated=time.time()))
    momentum = [fit_with_checks(h,MOM_EDGES,'momentum') for h in momentum_hist]
    save_json(out/'momentum_fits.json',momentum)
    angle_hist, ac, ai = scan(cache,cache_meta,'angle',out,mass,momentum)
    save_json(out/'angle_guides.json',[fit_angle(h) for h in angle_hist])
    all_counts = dict(mass=mc,momentum=pc,angle=ac)
    with (out/'histograms.tmp').open('wb') as f:
        np.savez_compressed(f,mass=mass_hist,momentum=momentum_hist,angle=angle_hist,
                            **{s+'_'+k:v for s,c in all_counts.items() for k,v in c.items()})
    (out/'histograms.tmp').replace(out/'histograms.npz')
    write_table(out,mass,momentum,all_counts)
    save_json(manifest,dict(signature=signature,partial=cache_meta['partial'],source_entries=cache_meta['source_entries'],schema_version=2,angle_unit='degree',
                            derived_cache=cache_meta,
                            passes=dict(mass=mi,momentum=pi,angle=ai),
                            usable_mass_bins=sum(f['usable'] for f in mass),usable_momentum_bins=sum(f['usable'] for f in momentum),
                            counts={s:{k:v.tolist() for k,v in c.items()} for s,c in all_counts.items()},
                            status='provisional selections for review, no final skim written'))
    plot_period(out)
    return name

def plot_period(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    meta = json.loads((out/'manifest.json').read_text())
    if meta.get('schema_version') != 2 or meta.get('angle_unit') != 'degree':
        raise ValueError('Old output uses the previous schema. Rerun the analysis before --plots-only.')
    guides=json.loads((out/'angle_guides.json').read_text())
    fits = dict(mass=json.loads((out/'mass_fits.json').read_text()),momentum=json.loads((out/'momentum_fits.json').read_text()))
    save_json(out/'progress.json',dict(stage='plots',state='plotting',updated=time.time()))
    stages=[('mass',MASS_EDGES,r'$M_{X,d}^{2}$ (GeV$^{2}$)',(.65,1.35),'Before exclusivity cuts'),
            ('momentum',MOM_EDGES,r'$|\vec p_{\rm miss}|$ (GeV)',(0,.30),'After missing-mass cut'),
            ('angle',ANGLE_EDGES,r'$\theta_\pi-\theta_{\pi,\rm pred}$ (deg)',(-6,6),'After both cuts')]
    with np.load(out/'histograms.npz') as data:
        for stage,edges,label,xlim,subtitle in stages:
            for xi in range(4):
                fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
                for ti,ax in enumerate(axes.flat):
                    b=xi*6+ti
                    hist=data[stage][b]
                    centers=(edges[:-1]+edges[1:])/2
                    ax.stairs(hist,edges,color='black',linewidth=1,label='Data')
                    notes=[]
                    if stage in fits:
                        fit=fits[stage][b]
                        if 'parameters' in fit:
                            x=np.linspace(*fit['fit_range'],400)
                            g,bg=components(x,fit['parameters'],fit['fit_range'])
                            ax.plot(x,g+bg,color='tab:red',lw=1.5,label='Gaussian + linear background')
                            ax.plot(x,bg,color='gray',ls=':',lw=1,label='Polynomial')
                        if fit['usable']:
                            lo,hi=fit['cut_low'],fit['cut_high']
                            ax.axvspan(lo,hi,color='tab:green',alpha=.09)
                            for cut in [lo,hi]: ax.axvline(cut,color='tab:green',ls='--',lw=1)
                            notes=[rf'$\mu$={fit["mu"]:.3f}, $\sigma$={fit["sigma"]:.3f}',f'Cut: [{lo:.3f}, {hi:.3f}]']
                            if fit['warnings']:notes.append('Fit flagged for review')
                        else:
                            notes=['UNRESOLVED FIT','See cuts.csv']
                    else:
                        ax.axvline(0,color='gray',ls='--',lw=1)
                        notes=[f'N = {int(data[stage+"_selected"][b]):,}']
                        if not fits['mass'][b]['usable'] or not fits['momentum'][b]['usable']:
                            notes.append('No usable upstream cut')
                        guide=guides[b]
                        if guide['usable']:
                            amp,mu,sigma=guide['parameters']
                            x=np.linspace(-3.5,3.5,400)
                            ax.plot(x,amp*np.exp(-.5*((x-mu)/sigma)**2),color='tab:red',lw=1.3)
                            notes.append(rf'$\mu$={mu:.3f}, $\sigma$={sigma:.3f} deg')
                    ax.text(.97,.96,'\n'.join(notes),ha='right',va='top',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
                    ax.set_title(f'Bin {b+1}: {TP[ti]:.2f} ≤ −t′ < {TP[ti+1]:.2f}',fontsize=10)
                    ax.set_xlabel(label)
                    ax.set_ylabel('Combinations / bin')
                    shown=xlim
                    if stage in fits and fits[stage][b]['usable']:
                        f=fits[stage][b]
                        shown=(min(xlim[0],f['fit_range'][0]-.01,f['cut_low']-.01),max(xlim[1],f['fit_range'][1]+.01,f['cut_high']+.01))
                        shown=(max(edges[0],shown[0]),min(edges[-1],shown[1]))
                    ax.set_xlim(*shown)
                    # Scale to the displayed local region, not the discarded broad background.
                    local=hist[(centers>=shown[0]) & (centers<=shown[1])]
                    ax.set_ylim(0,max(1,float(local.max())*1.28))
                    ax.spines[['top','right']].set_visible(False)
                marker=' | PARTIAL SAMPLE' if meta['partial'] else ''
                fig.suptitle(f'{out.name} | {XB[xi]:.2f} ≤ xB < {XB[xi+1]:.2f}\n{subtitle}{marker}',fontsize=13)
                fig.savefig(out/f'{stage}_xb{xi+1}.png',dpi=150)
                plt.close(fig)
    save_json(out/'progress.json',dict(stage='complete',state='complete',updated=time.time()))
    print(f'DONE {out.name}: 12 PNG pages',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets',type=Path,default=Path(__file__).resolve().parents[2]/'imports/datasets.json')
    parser.add_argument('--input-root',type=Path)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--periods',nargs='+',default=PERIODS)
    parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--step-size',default='100 MB')
    parser.add_argument('--max-entries',type=int)
    parser.add_argument('--status',action='store_true')
    parser.add_argument('--plots-only',action='store_true')
    parser.add_argument('--force',action='store_true',help='Refit using a matching derived cache')
    parser.add_argument('--rebuild-cache',action='store_true',help='Repeat ROOT scan even if derived cache matches')
    args=parser.parse_args()
    if args.workers<1 or (args.max_entries is not None and args.max_entries<1):parser.error('workers and max-entries must be positive')
    if args.status:
        for name in args.periods:
            path=args.output/name/'progress.json'
            if path.exists():
                v=json.loads(path.read_text());pct=f' {v["percent"]:.1f}%' if 'percent' in v else ''
                print(f'{name}: {v["stage"]}, {v["state"]}{pct}, updated {time.time()-v["updated"]:.0f}s ago')
            else:print(f'{name}: not started')
        return
    args.output.mkdir(parents=True,exist_ok=True)
    print(f'Output: {args.output.resolve()}',flush=True)
    if args.plots_only:
        for name in args.periods:plot_period(args.output/name)
    else:
        config=json.loads(args.datasets.read_text())
        root=args.input_root or Path(config['campaign'])/'root'
        jobs=[]
        for name in args.periods:
            if name not in config['files']:parser.error(f'Unknown dataset {name}')
            source=root/config['files'][name]
            if not source.is_file():parser.error(f'Missing input {source}')
            jobs.append((name,str(source),config['tree'],str(args.output),dict(max_entries=args.max_entries,step_size=args.step_size,force=args.force,rebuild_cache=args.rebuild_cache)))
        print(f'Starting {min(args.workers,len(jobs))} parallel period workers',flush=True)
        with ProcessPoolExecutor(max_workers=min(args.workers,len(jobs))) as pool:
            for future in as_completed([pool.submit(run_period,j) for j in jobs]):future.result()
    links=[]
    for name in args.periods:
        links.append(f'<h2>{html.escape(name)}</h2><p><a href="{name}/cuts.csv">Fit results and cut flow</a> | <a href="{name}/stability.csv">Stability checks</a></p>')
        for stage in ['mass','momentum','angle']:
            links.append('<p>'+stage+': '+' | '.join(f'<a href="{name}/{stage}_xb{i}.png">xB {i}</a>' for i in range(1,5))+'</p>')
    (args.output/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Exclusivity selection</title><h1>Preliminary exclusivity selection</h1><p>Bin-dependent provisional cuts. Check fit warnings before use.</p>'+''.join(links))
    print(f'Browse {args.output/"index.html"}',flush=True)

if __name__=='__main__':main()
