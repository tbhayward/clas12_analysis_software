#!/usr/bin/env python3
"""Matched-global-selection closure, physical mixtures and production cut refits."""
from __future__ import annotations
import argparse,copy,csv,json,subprocess,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts

MODELS=('current','add_pt');SAMPLES=('data','mc','background','pi0_data','pi0_mc')

def mask(rows,top,block=None,kind='data',omit=None,domain='global'):
    selected=(rows[:,0]==top)&(rows[:,8 if domain=='global' else 10]>.5)
    if block:
        for name,stats in block[kind].items():
            if name==omit:continue
            values=rows[:,1+base.OBS.index(name)]
            selected &= np.isfinite(values)&(values<=float(stats['cut_high']))
            if stats.get('mode')!='upper_quantile':selected &= values>=float(stats['cut_low'])
    return selected

def load_exports(folder,periods,vectors):
    work=Path(json.loads((folder/'scratch_location.json').read_text())['directory']);arrays={};files={}
    for jobfile in sorted((work/'precut_arrays').glob('*.json')):
        if jobfile.name.endswith('.bin.json'):continue
        job=json.loads(jobfile.read_text());p=job['period'];sample=job['sample'];path=jobfile.with_suffix('.bin')
        if p not in periods:continue
        metadata=json.loads(Path(str(path)+'.json').read_text());shape=metadata['shape']
        if len(shape)!=2 or shape[1]!=12 or path.stat().st_size!=int(np.prod(shape))*8:raise ValueError('Invalid pre-cut export: '+str(path))
        names=['shared'] if sample in ('data','pi0_data') else [name for name,pars in [('nominal',[0.]*20),*[(c,vectors[p,c]) for c in MODELS]] if job['parameters']==pars]
        for name in names:
            key=p,name,sample
            if key not in arrays:arrays[key]=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape));files[key]=path
    required={(p,'shared' if s in ('data','pi0_data') else c,s) for p in periods for c in MODELS for s in SAMPLES}|{(p,'nominal',s) for p in periods for s in ('mc','background','pi0_mc')}
    if not required<=arrays.keys():raise FileNotFoundError('Previous pre-cut scratch exports are missing or incomplete')
    return arrays,files

def arrays_for(exports,period,scenario):return {s:exports[period,'shared' if s in ('data','pi0_data') else scenario,s] for s in SAMPLES}

def scales(period,arrays,blocks,probabilities,exposure,multiplier,shared_dvcs=False):
    cfg=next(p for p in cuts.production.PERIODS if p.key==period);control_data=0;control_mc=0
    for t,top in enumerate(cuts.production.TOPOLOGIES):
        block=blocks[cuts.production.cut_json_key('eppi0',cfg,top)]
        control_data+=int(mask(arrays['pi0_data'],t,block,'data').sum());control_mc+=int(mask(arrays['pi0_mc'],t,block,'mc').sum())/probabilities['pi0_mc']
    control=control_data/control_mc if control_mc>0 and control_data>0 else float('nan');alpha=control*exposure*multiplier;result={}
    for t,top in enumerate(cuts.production.TOPOLOGIES):
        block=blocks[cuts.production.cut_json_key('dvcs',cfg,top)];kind='data' if shared_dvcs else 'mc'
        d=int(mask(arrays['data'],t,block,'data').sum());s=int(mask(arrays['mc'],t,block,kind).sum())/probabilities['mc'];b=int(mask(arrays['background'],t,block,kind).sum())/probabilities['background']
        beta=(d-alpha*b)/s if s>0 and d>alpha*b and np.isfinite(alpha) else float('nan')
        result[t]=dict(signal_scale=beta,background_scale=alpha,control_scale=control,data_anchor=d,signal_mc_anchor=s,background_mc_anchor=b,control_data_anchor=control_data,control_mc_anchor=control_mc,status='ok' if np.isfinite(beta) else 'invalid control or nonpositive signal residual')
    return result

def add_background_proxy(rows,period,new_arrays,old_arrays,new_blocks,old_blocks,new_scales,old_scales,prob):
    """Unit-weighted subtracted-yield/acceptance diagnostic; never clip negatives."""
    cfg=next(p for p in cuts.production.PERIODS if p.key==period)
    for row in rows:
        top=base.TOPS.index(row['topology']);key=cuts.production.cut_json_key('dvcs',cfg,cuts.production.TOPOLOGIES[top]);bounds=(row['t_low'],row['t_high']);values=[]
        for tag,arrays,blocks,sc in [('old',old_arrays,old_blocks,old_scales),('new',new_arrays,new_blocks,new_scales)]:
            background=arrays['background'];selected=mask(background,top,blocks[key],'mc')&(background[:,9]>=bounds[0])&(background[:,9]<bounds[1])
            estimate=sc[top]['background_scale']*selected.sum()/prob['background'];signal=row[tag+'_data_selected']-estimate
            row[tag+'_background_estimate']=float(estimate);row[tag+'_subtracted_data']=float(signal);values.append(signal)
        ratio=row['relative_acceptance'];valid=all(np.isfinite(v) and v>0 for v in values) and np.isfinite(ratio) and ratio>0
        row['subtracted_yield_acceptance_proxy_ratio']=values[1]/values[0]/ratio if valid else float('nan')
        row['proxy_status']='unit-weighted diagnostic; not a full cross section' if valid else 'insufficient/invalid subtracted yield or acceptance'


def histogram(values,selected,edges):
    x=values[selected];finite=x[np.isfinite(x)];inside=np.histogram(finite,edges)[0]
    return np.r_[np.count_nonzero(finite<edges[0]),inside,np.count_nonzero(finite>=edges[-1])].astype(float)

def physical_prediction(arrays,top,vi,edges,scale,prob,block,view):
    domain='minimal' if view=='minimal_pre' else 'global';omit=base.OBS[vi] if view=='global_nminus1' else None
    use=block if view in ('global_nminus1','global_final') else None
    kind='mc' if view=='global_final' else 'data'
    d=histogram(arrays['data'][:,vi+1],mask(arrays['data'],top,use,'data',omit,domain),edges)
    signal=histogram(arrays['mc'][:,vi+1],mask(arrays['mc'],top,use,kind,omit,domain),edges)*scale['signal_scale']/prob['mc']
    background=histogram(arrays['background'][:,vi+1],mask(arrays['background'],top,use,kind,omit,domain),edges)*scale['background_scale']/prob['background']
    return d,signal,background

def shape_diagnostics(data,model,edges):
    d=base.diag.shape_metrics(data,edges,0,0)
    m=base.diag.shape_metrics(model,edges,0,0) if np.all(np.isfinite(model)) and model.sum()>0 else None
    width=m['core_spread'] if m else float('nan')
    return dict(data_core_spread=d['core_spread'],mc_core_spread=width,data_over_mc_core_width=d['core_spread']/width if width>0 else float('nan'),median_shift_MC_spreads=(d['median']-m['median'])/width if width>0 else float('nan'),data_in_range_fraction=d['in_range_fraction'],mc_in_range_fraction=m['in_range_fraction'] if m else float('nan'))


def render_physical(period,spec,cases,saved_blocks,global_blocks,anchor,final_anchor,prob,output):
    import matplotlib.pyplot as plt
    links=[];rows=[];cfg=next(p for p in cuts.production.PERIODS if p.key==period)
    for top,topology in enumerate(base.TOPS):
        key=cuts.production.cut_json_key('dvcs',cfg,cuts.production.TOPOLOGIES[top])
        for view in ['minimal_pre','global_pre','global_nminus1','global_final']:
            fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained')
            for vi,v in enumerate(spec['variables']):
                ax=axes.flat[vi];edges=base.diag.variable_edges(v);block=global_blocks['add_pt'][key] if view=='global_final' else saved_blocks[key]
                sc=final_anchor['add_pt'][top] if view=='global_final' else anchor['add_pt'][top]
                d,s,b=physical_prediction(cases['add_pt'],top,vi,edges,sc,prob,block,view);prediction=s+b;den=max(1,d.sum());curves=[('DATA',d,'#205cab','-'),('Particle-level DVCS + π⁰',prediction,'#22946b','-'),('DVCS component',s,'#22946b','--'),('π⁰ component',b,'#22946b',':')]
                if view!='global_final':
                    _,ss,bb=physical_prediction(cases['current'],top,vi,edges,anchor['current'][top],prob,saved_blocks[key],view);curves.insert(1,('Original training: summed MC',ss+bb,'#d45c18','-'))
                for label,h,color,style in curves:
                    if np.all(np.isfinite(h)):ax.stairs(h[1:-1]/den/np.diff(edges),edges,color=color,ls=style,label=label,lw=1.2)
                if not np.all(np.isfinite(prediction)):ax.text(.03,.6,sc['status'],transform=ax.transAxes,color='red')
                if v['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_,_ in curves if np.all(np.isfinite(h))]))
                else:
                    qq=[base.diag.quantile(h[1:-1],edges,q) for _,h,_,_ in curves if np.all(np.isfinite(h)) for q in (.01,.99)];qq=[q for q in qq if np.isfinite(q)]
                    if qq:
                        lo,hi=min(qq),max(qq);pad=max((hi-lo)*.08,.001);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                ax.set(xlabel=v['label'],ylabel='Density / DATA total');ax.legend(fontsize=7);ax.grid(alpha=.15)
                rows.append(dict(channel='DVCS',scenario='add_pt',period=period,topology=topology,view=view,variable=v['name'],data_entries=d.sum(),predicted_entries=prediction.sum(),predicted_over_data=prediction.sum()/den,signal_entries=s.sum(),background_entries=b.sum(),status=sc['status'],**shape_diagnostics(d,prediction,edges)))
            for ax in axes.flat[7:]:ax.axis('off')
            axes.flat[7].text(0,1,'Particle-level prediction only: no histogram morph.\nπ⁰ normalization comes from the two-photon control.\nOne signal scale per selected topology;\nno per-observable normalization fit.\nN−1 drops only the displayed exclusivity cut;\nall configured global cuts remain.',va='top')
            fig.suptitle(period+' · '+topology+' · '+view+' · add-pT candidate')
            path=output/'physical_plots'/(period+'_'+topology+'_'+view+'.png');path.parent.mkdir(exist_ok=True);fig.savefig(path,dpi=105);plt.close(fig);links.append(path.relative_to(output).as_posix())
        # Two-photon controls: same saved DATA windows for DATA and MC in N−1.
        key=cuts.production.cut_json_key('eppi0',cfg,cuts.production.TOPOLOGIES[top])
        for view in ['global_pre','global_nminus1']:
            fig,axes=plt.subplots(3,3,figsize=(15,11),layout='constrained')
            for vi,v in enumerate(spec['variables']):
                ax=axes.flat[vi];edges=base.diag.variable_edges(v);block=saved_blocks[key] if view=='global_nminus1' else None;omit=v['name'] if block else None
                d=histogram(cases['add_pt']['pi0_data'][:,vi+1],mask(cases['add_pt']['pi0_data'],top,block,'data',omit),edges);den=max(1,d.sum());curves=[('DATA',d,'#205cab')]
                for c,color in [('current','#d45c18'),('add_pt','#22946b')]:
                    a=cases[c]['pi0_mc'];m=histogram(a[:,vi+1],mask(a,top,block,'data',omit),edges)*anchor[c][top]['control_scale']/prob['pi0_mc'];curves.append((c+' two-photon MC',m,color))
                    rows.append(dict(channel='two_photon_control',scenario=c,period=period,topology=topology,view=view,variable=v['name'],data_entries=d.sum(),predicted_entries=m.sum(),predicted_over_data=m.sum()/den,signal_entries=m.sum(),background_entries=0.,status='ok' if np.isfinite(anchor[c][top]['control_scale']) else 'invalid control normalization',**shape_diagnostics(d,m,edges)))
                for label,h,color in curves:ax.stairs(h[1:-1]/den/np.diff(edges),edges,label=label,color=color)
                if v['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_ in curves]))
                ax.set(xlabel=('π⁰ opening angle (deg)' if v['name']=='theta_gamma_gamma' else v['label'].replace('γ','π⁰')),ylabel='Density / DATA total');ax.legend(fontsize=7);ax.grid(alpha=.15)
            for ax in axes.flat[7:]:ax.axis('off')
            axes.flat[7].text(0,1,'Full global selection. No histogram-template morph.\nControl scale fixed by selected two-photon totals.\nN−1 uses the same DATA windows on both samples.\nThe π⁰ generator is the control signal;\nthis is not the one-photon DVCS mixture.',va='top')
            fig.suptitle(period+' · '+topology+' · two-photon '+view);path=output/'control_plots'/(period+'_'+topology+'_'+view+'.png');path.parent.mkdir(exist_ok=True);fig.savefig(path,dpi=105);plt.close(fig);links.append(path.relative_to(output).as_posix())
    return links,rows

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--followup-study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--workers',type=int,default=10);parser.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'))
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10:parser.error('Workers must be 1–10')
    source=args.followup_study.resolve()/'cut_redetermination';output=args.output.resolve()
    if output.exists():parser.error('Output exists; choose a new directory')
    old=json.loads((source/'manifest.json').read_text())
    if not all(c in old['candidate_hashes'] for c in MODELS):parser.error('Use a completed add-pT follow-up with current and add_pt cut results')
    previous_domain=old['arguments'].get('fit_selection','minimal')
    ref=Path(old['arguments']['reference_study']);run=Path(old['arguments']['response_run']);training=Path(old['arguments']['training_study']);spec=json.loads((source/'spec.json').read_text())
    if (run/'with_bias/cache_manifest.json').is_file():run=run/'with_bias'
    reference=json.loads((ref/'manifest.json').read_text());cm=json.loads((run/'cache_manifest.json').read_text())
    if base.digest(ref/'manifest.json')!=old['reference_sha256'] or base.digest(run/'cache_manifest.json')!=old['response_cache_sha256']:raise ValueError('Reference/cache metadata changed')
    paths=cuts.candidate_paths({c:run for c in MODELS},spec['periods'],training);vectors={key:cuts.load_candidate_parameters(p) for key,p in paths.items()}
    for (p,c),path in paths.items():
        if base.digest(path)!=old['candidate_hashes'][c][p]:raise ValueError('Candidate differs from source cut study')
    output.mkdir(parents=True)
    try:exports,files=load_exports(source,spec['periods'],vectors);print('[reuse] Pre-cut scratch exports; no ROOT reads or response evaluation',flush=True)
    except (FileNotFoundError,KeyError):
        rebuilt=output/'rebuilt_exports';print('[exports] Missing scratch exports; rebuilding exports only',flush=True)
        subprocess.run([sys.executable,str(Path(cuts.__file__).resolve()),'--reference-study',str(ref),'--response-run',str(run),'--training-study',str(training),'--workers',str(args.workers),'--scratch-root',str(args.scratch_root),'--output',str(rebuilt),'--export-only'],check=True)
        exports,files=load_exports(rebuilt,spec['periods'],vectors)
    base.dump(output/'manifest.json',dict(source_cut_manifest_sha256=base.digest(source/'manifest.json'),candidate_hashes=old['candidate_hashes'],fit_selection='global',source_fit_selection=previous_domain,export_hashes={str(path):base.digest(path) for path in set(files.values())},production_fitter_sha256=base.digest(Path(cuts.production.__file__))))
    results={};tasks={}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for period in spec['periods']:
            p=next(p for p in cuts.production.PERIODS if p.key==period)
            for c in MODELS:
                a=arrays_for(exports,period,c)
                for top,topology in enumerate(base.TOPS):tasks[executor.submit(cuts.refit,p,top,a,output/'global_cut_refits'/period/c/topology,'global')]=(period,c)
        for future in as_completed(tasks):
            p,c=tasks[future];results.setdefault((p,c),{}).update(future.result());print('[global cuts] '+p+'/'+c+' topology complete',flush=True)
    links=[];rates=[];normalizations=[];acceptance=[]
    for period in spec['periods']:
        p=next(p for p in cuts.production.PERIODS if p.key==period);cases={c:arrays_for(exports,period,c) for c in MODELS};prob={s:cm['samples'][period+'__'+s]['probability'] for s in ('mc','background','pi0_mc')};new={c:results[period,c] for c in MODELS};oldblocks={c:json.loads((source/period/c/'combined_cuts_95.json').read_text()) for c in MODELS}
        anchor={};final={};minimal={}
        for c in MODELS:
            multiplier=json.loads(paths[period,c].read_text())['background_multiplier'];exposure=spec.get('control_to_leakage_exposure_ratio',1.)
            minimal[c]=scales(period,cases[c],oldblocks[c],prob,exposure,multiplier)
            anchor[c]=scales(period,cases[c],spec['exclusivity'],prob,exposure,multiplier,True);final[c]=scales(period,cases[c],new[c],prob,exposure,multiplier)
            base.dump(output/'global_cut_refits'/period/c/'combined_cuts_95.json',new[c])
            for top in range(3):
                for stage,values in [('saved_common_anchor',anchor[c]),('global_refit_anchor',final[c])]:normalizations.append(dict(period=period,scenario=c,topology=base.TOPS[top],stage=stage,**values[top]))
            z=cuts.comparison_rows(period,c,cases[c],cases[c],new[c],oldblocks[c],prob)
            add_background_proxy(z,period,cases[c],cases[c],new[c],oldblocks[c],final[c],minimal[c],prob)
            for row in z:row['scenario']+='__global_vs_'+('minimal' if previous_domain=='minimal' else 'previous_global')+'_refit'
            acceptance+=z
        z=cuts.comparison_rows(period,'add_pt',cases['add_pt'],cases['current'],new['add_pt'],new['current'],prob)
        add_background_proxy(z,period,cases['add_pt'],cases['current'],new['add_pt'],new['current'],final['add_pt'],final['current'],prob)
        for row in z:row['scenario']+='__versus_current_global_refit'
        acceptance+=z
        ll,rr=render_physical(period,spec,cases,spec['exclusivity'],new,anchor,final,prob,output);links+=ll;rates+=rr
    base.table(output/'normalizations.csv',normalizations);base.table(output/'physical_prediction_counts.csv',rates);base.table(output/'relative_acceptance.csv',acceptance)
    gallery='<table><tr><th>Period / topology</th><th>Before exclusivity</th><th>Leave one cut out</th><th>Final sample</th><th>Legacy minimal</th><th>Two-photon control</th><th>Extra template-fit diagnostics</th></tr>'+''.join('<tr><td>'+p+' / '+t+'</td>'+''.join('<td><a href="physical_plots/'+p+'_'+t+'_'+v+'.png">'+label+'</a></td>' for v,label in [('global_pre','Full globals'),('global_nminus1','Full globals'),('global_final','Own refitted DATA/MC cuts'),('minimal_pre','Minimal')])+'<td><a href="control_plots/'+p+'_'+t+'_global_pre.png">Before cuts</a> / <a href="control_plots/'+p+'_'+t+'_global_nminus1.png">Leave one out</a></td><td>'+''.join('<a href="global_cut_refits/'+p+'/'+c+'/'+t+'/'+channel+'_cut_development.png">'+c+' '+channel+'</a> ' for c in MODELS for channel in ['dvcs','pi0'])+'</td></tr>' for p in spec['periods'] for t in base.TOPS)+'</table>'
    (output/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Matched-selection closure</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}li{margin:8px}td,th{padding:8px;border-bottom:1px solid #ddd}</style><h1>Matched-selection closure and full-global cut determination</h1><p>All primary MC curves are particle-level predictions, with no additional histogram shift or smearing. DVCS plots show the DVCS+π⁰ sum and its components. π⁰ normalization is constrained by two-photon selected counts; the signal scale is set once per topology in a fully selected anchor and held fixed across variables. Panel densities use the DATA total, so a prediction/data area difference is retained, not fitted away. The anchor totals themselves are normalized by construction.</p><p>global_pre retains all configured global cuts and no exclusivity cuts. global_nminus1 additionally uses the same saved DATA windows on DATA and MC, omitting only the displayed exclusivity cut. global_final uses each sample’s own full-global refitted cuts. Minimal preselection is retained for comparison. All configured global cuts remain, including any global pT cut.</p><p><a href="normalizations.csv">Control and signal normalization anchors</a> · <a href="physical_prediction_counts.csv">Counts and normalization flags</a> · <a href="relative_acceptance.csv">Acceptance: full-global versus minimal refits, and add-pT versus original training</a></p><p>global_cut_refits contains the unchanged production cut algorithm’s additional histogram-template fit diagnostics. Those are separate from the physical-prediction plots and must not be interpreted as particle-response validation. Production output JSONs are not overwritten. The acceptance table also supplies a unit-weighted π⁰-subtracted-yield/acceptance proxy with the updated control normalization; invalid residuals are flagged. It excludes the rest of the production normalization/correction chain and is not a final cross section. The upstream two-photon mass-window limitation remains.</p>'+gallery)
    print('[done] '+str(output/'REPORT.html'),flush=True)

if __name__=='__main__':main()
