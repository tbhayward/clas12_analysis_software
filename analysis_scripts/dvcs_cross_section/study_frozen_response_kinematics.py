#!/usr/bin/env python3
"""Frozen joint response: low/high momentum and polar-angle validation, no fits."""
from __future__ import annotations
import argparse,copy,csv,html,json
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts
import study_particle_selection_closure as closure
from study_particle_common_cut_stability import export_digest
CATEGORIES={'electron_FD':(14,15,None,None),'proton_FD':(16,17,0,None),'proton_CD':(16,17,1,None),'photon1_FD':(20,21,None,(23,1)),'photon1_FT':(20,21,None,(23,0)),'photon2_FD':(24,25,None,(27,1)),'photon2_FT':(24,25,None,(27,0))}


def category_mask(a,category):
    p,theta,proton,photon=CATEGORIES[category];keep=a[:,8]>.5
    if proton==0:keep&=a[:,0]==0
    if proton==1:keep&=a[:,0]!=0
    if photon:keep&=a[:,photon[0]]==photon[1]
    return keep&np.isfinite(a[:,p])&np.isfinite(a[:,theta])


def choose_splits(data,control):
    rows=[]
    for category,(p,theta,_,_) in CATEGORIES.items():
        source_category=category.replace('photon2','photon1');source_p,source_theta,_,_=CATEGORIES[source_category]
        a=data;mask=category_mask(a,source_category)&(a[:,11]<.5)
        for axis,col,source_col,factor in [('momentum',p,source_p,1.),('angle',theta,source_theta,180/np.pi)]:
            values=a[mask,source_col]*factor
            rows.append(dict(category=category,axis=axis,column=col,factor=factor,threshold=float(np.median(values)) if len(values) else float('nan'),source='epgamma DATA half 0; shared photon-role edge' if category.startswith('photon') else 'epgamma DATA half 0',events=len(values),status='usable' if len(values)>=100 and len(np.unique(values))>1 else 'insufficient split support'))
    return rows


def masks(a,split,side):
    values=a[:,int(split['column'])]*split['factor'];keep=category_mask(a,split['category'])
    return keep&((values<split['threshold']) if side=='low' else (values>=split['threshold']))


def histogram(a,top,block,variable,split,side,edges,half=None,selected_base=None):
    if selected_base is None:selected=closure.mask(a,top,block,'data',omit=variable)&masks(a,split,side)
    else:
        selected=selected_base.copy()
        for name,window in block['data'].items():
            if name==variable:continue
            values=a[:,1+base.OBS.index(name)];selected&=np.isfinite(values)&(values<=float(window['cut_high']))
            if window.get('mode')!='upper_quantile':selected&=values>=float(window['cut_low'])
    if half is not None:selected&=a[:,11]==half
    return closure.histogram(a[:,1+base.OBS.index(variable)],selected,edges)


def get_scales(period,run,spec,cache,reference):
    norm={r['topology']:r for r in csv.DictReader((reference/'normalizations.csv').open()) if r['period']==period}
    selected=json.loads((run/period/'final_control_response.json').read_text())['selected'];candidate=json.loads((run/period/'candidate_parameters.json').read_text());mult=candidate['background_multiplier']
    # Signal/leakage selected yields are recomputed from exports by the caller.
    return norm,np.asarray(selected,dtype=float),mult


def analyze(period,cases,spec,scales,prob,splits,out,min_entries):
    import matplotlib.pyplot as plt
    cfg=next(p for p in cuts.production.PERIODS if p.key==period);rows=[];examples=[]
    for split in splits:
        if split['status']!='usable':continue
        category=split['category'];axis=split['axis'];channels=['epgammagamma'] if category.startswith('photon2') else ['epgamma','epgammagamma'];matrix=[];labels=[];best={}
        for channel in channels:
            for top,topname in enumerate(base.TOPS):
                block=spec['exclusivity'][cuts.production.cut_json_key('eppi0' if channel=='epgammagamma' else 'dvcs',cfg,cuts.production.TOPOLOGIES[top])]
                for side in ('low','high'):
                    display=[]
                    selections={s:masks(cases[s],split,side)&(cases[s][:,0]==top) for s in ('data','pi0_data')}
                    selections.update({scenario+'_'+sample:masks(a,split,side)&(a[:,0]==top) for scenario,part in cases.items() if isinstance(part,dict) for sample,a in part.items()})
                    for vi,variable in enumerate(spec['variables']):
                        name=variable['name'];edges=base.diag.variable_edges(variable);d=histogram(cases['pi0_data' if channel=='epgammagamma' else 'data'],top,block,name,split,side,edges,selected_base=selections['pi0_data' if channel=='epgammagamma' else 'data'])
                        prediction={};variance={};raw={}
                        for scenario in ('nominal','candidate'):
                            models=[];variances=[];rawcounts=[]
                            for half in (0,1):
                                if channel=='epgamma':
                                    signal=histogram(cases[scenario]['mc'],top,block,name,split,side,edges,half,selected_base=selections[scenario+'_mc']);background=histogram(cases[scenario]['background'],top,block,name,split,side,edges,half,selected_base=selections[scenario+'_background']);beta,alpha=scales[scenario,half,top];w1=2*beta/prob['mc'];w2=2*alpha/prob['background'];m=w1*signal+w2*background;v=w1*w1*signal+w2*w2*background;rc=signal.sum()+background.sum()
                                else:
                                    c=histogram(cases[scenario]['pi0_mc'],top,block,name,split,side,edges,half,selected_base=selections[scenario+'_pi0_mc']);m=c.astype(float);v=m.copy();rc=c.sum()
                                models.append(m);variances.append(v);rawcounts.append(rc)
                            prediction[scenario]=models;variance[scenario]=variances;raw[scenario]=rawcounts
                        scores={}
                        for half in (0,1):
                            values={}
                            for scenario in ('nominal','candidate'):
                                status='usable' if d.sum()>=min_entries and raw['nominal'][half]>=min_entries and raw['candidate'][half]>=min_entries else 'low statistics'
                                score,bins=base.shape_score(d,prediction[scenario][half],variance[scenario][half],int(raw['nominal'][half]),min_entries,prediction['nominal'][half],variance['nominal'][half]) if status=='usable' else (float('nan'),0)
                                if status=='usable' and not np.isfinite(score):status='insufficient populated score bins or invalid normalization'
                                metrics_d=base.diag.shape_metrics(d,edges,0,0);metrics_m=base.diag.shape_metrics(prediction[scenario][half],edges,0,0)
                                width=metrics_m['core_spread'];values[scenario]=score
                                rows.append(dict(period=period,category=category,axis=axis,threshold=split['threshold'],channel=channel,topology=topname,side=side,variable=name,scenario=scenario,mc_half=half,data_count=d.sum(),raw_mc=raw[scenario][half],shape_score=score,usable_bins=bins,data_over_mc_width=metrics_d['core_spread']/width if width>0 else float('nan'),median_shift_mc_widths=(metrics_d['median']-metrics_m['median'])/width if width>0 else float('nan'),data_in_range_fraction=metrics_d['in_range_fraction'],mc_in_range_fraction=metrics_m['in_range_fraction'],status=status,role='reserved observable' if name in ('theta_gamma_gamma','Mx2_2','theta') else 'fitted observable'))
                            scores[half]=values
                        n=scores[1]['nominal'];c=scores[1]['candidate'];display.append(100*(1-c/n) if np.isfinite(n) and np.isfinite(c) and n>0 else np.nan)
                        if np.isfinite(c) and ((channel,side) not in best or c>best[channel,side][0]):best[channel,side]=(c,topname,name,d,[(prediction[s][0]+prediction[s][1])/2 for s in ('nominal','candidate')],edges)
                    matrix.append(display);labels.append(channel+' '+topname+' '+side)
        folder=out/'plots'/period/category;folder.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(11,max(4,len(labels)*.35)),layout='constrained');values=np.ma.masked_invalid(np.array(matrix));cmap=plt.get_cmap('RdYlGn').copy();cmap.set_bad('#dddddd');im=ax.imshow(values,vmin=-100,vmax=100,cmap=cmap,aspect='auto');ax.set(xticks=range(len(base.OBS)),xticklabels=base.OBS,yticks=range(len(labels)),yticklabels=labels)
        for index,value in np.ndenumerate(values):
            if not np.ma.is_masked(value):ax.text(index[1],index[0],f'{float(value):+.0f}',ha='center',va='center',fontsize=7)
        fig.colorbar(im,ax=ax,label='Other-MC-half score reduction (%)');ax.set_title(f'{period} · {category} · {axis}, split={split["threshold"]:.4g}\nFrozen joint response; gray=insufficient support; numeric values retained outside color range')
        path=folder/(axis+'_scores.png');fig.savefig(path,dpi=110);plt.close(fig)
        residual_lookup={(r['channel']+' '+r['topology']+' '+r['side'],r['variable']):r for r in rows if r['category']==category and r['axis']==axis and r['scenario']=='candidate' and r['mc_half']==1}
        fig,axes=plt.subplots(1,2,figsize=(18,max(4,len(labels)*.35)),layout='constrained')
        for ax,field,title,lo,hi in [(axes[0],'data_over_mc_width','DATA/MC core width: target 1',.5,1.5),(axes[1],'median_shift_mc_widths','DATA−MC median / MC core width: target 0',-1,1)]:
            vals=np.array([[residual_lookup[label,name][field] if residual_lookup[label,name]['status']=='usable' else np.nan for name in base.OBS] for label in labels]);im=ax.imshow(np.ma.masked_invalid(vals),vmin=lo,vmax=hi,cmap='coolwarm',aspect='auto');ax.set(title=title,xticks=range(7),xticklabels=base.OBS,yticks=range(len(labels)),yticklabels=labels)
            for index,value in np.ndenumerate(vals):
                if np.isfinite(value):ax.text(index[1],index[0],f'{value:.2f}',ha='center',va='center',fontsize=7)
            fig.colorbar(im,ax=ax,shrink=.8)
        fig.suptitle(f'{period} · {category} · {axis}: frozen candidate low/high residuals; descriptive, no fit')
        residual_path=folder/(axis+'_residuals.png');fig.savefig(residual_path,dpi=110);plt.close(fig)
        # Selected examples are explicitly diagnostic: largest residual per channel/side.
        fig,axes=plt.subplots(len(channels),2,figsize=(12,4*len(channels)),squeeze=False,layout='constrained')
        for ci,ch in enumerate(channels):
            for si,side in enumerate(('low','high')):
                ax=axes[ci,si]
                if (ch,side) not in best:ax.text(.1,.5,'Insufficient support');continue
                score,top,name,d,models,edges=best[ch,side]
                for label,h,color in [('DATA',d,'#205cab'),('Nominal MC',models[0],'#d45c18'),('Frozen joint MC',models[1],'#22946b')]:
                    if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,color=color,label=label)
                if name=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[d,*models]))
                else:
                    qs=[base.diag.quantile(h[1:-1],edges,q) for h in [d,*models] if h.sum()>0 for q in (.01,.99)];qs=[q for q in qs if np.isfinite(q)]
                    if qs:
                        lo,hi=min(qs),max(qs);pad=max((hi-lo)*.08,.001);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                ax.set(title=f'{ch} {top} {side} · {name}',xlabel=name,ylabel='Unit-area density');ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle(f'{period} · {category} · {axis}: largest supported residual examples (not training selection)');fig.savefig(folder/(axis+'_examples.png'),dpi=110);plt.close(fig)
        examples.append(dict(period=period,category=category,axis=axis,scores=path.relative_to(out).as_posix(),residuals=residual_path.relative_to(out).as_posix(),examples=(folder/(axis+'_examples.png')).relative_to(out).as_posix()))
    return rows,examples


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joint-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=10);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--reuse-exports',type=Path);p.add_argument('--min-entries',type=int,default=100);args=p.parse_args(argv)
    if not 1<=args.workers<=10 or args.min_entries<10:p.error('Invalid workers/minimum entries')
    run=args.joint_run.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a fresh directory')
    meta=json.loads((run/'manifest.json').read_text());spec=json.loads((run/'spec.json').read_text());reference=Path(meta['arguments']['reference_study']);source=json.loads((reference/'manifest.json').read_text());cm=json.loads((run/'cache_manifest.json').read_text())
    for key,name in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/name)!=spec['input_hashes'][key]:p.error('Saved global-cut source changed')
    if cm['reference_sha256']!=base.digest(reference/'manifest.json'):p.error('Cache/reference provenance differs')
    for name in ('particle_response_core.h','particle_control_core.h'):
        if base.digest(base.SOURCE/name)!=meta['source_hashes'][name]:p.error('Response law changed since frozen fit')
    if base.digest(reference/'manifest.json')!=meta['reference_manifest_sha256']:p.error('Reference changed')
    candidates={period:run/period/'candidate_parameters.json' for period in spec['periods']}
    if not all(json.loads(path.read_text()).get('joint_control_fit') for path in candidates.values()):p.error('Requires completed joint-fit candidates')
    pars={period:cuts.load_candidate_parameters(path) for period,path in candidates.items()}
    fingerprint=dict(joint_manifest=base.digest(run/'manifest.json'),cache_manifest=base.digest(run/'cache_manifest.json'),candidates={k:base.digest(v) for k,v in candidates.items()},export_sources={n:base.digest(base.SOURCE/n) for n in ('particle_response_cut_arrays.h','particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h')})
    arrays={period:{'nominal':{},'candidate':{}} for period in spec['periods']};records=[]
    if args.reuse_exports:
        old=json.loads((args.reuse_exports/'export_manifest.json').read_text())
        if old['fingerprint']!=fingerprint:p.error('Reused inputs differ')
        records=old['records']
        for r in records:
            if export_digest(r['path'])!=r['sha256']:p.error('Export changed')
    else:
        jobs=[]
        for period in spec['periods']:
            for sample in ('data','pi0_data'):
                entries=[r for r in source['inventory'] if r['period']==period and r['sample']==sample]
                if not entries:p.error('Missing DATA inventory')
                for r in entries:
                    st=Path(r['path']).stat()
                    if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):p.error('DATA changed')
                jobs.append(dict(period=period,sample=sample,files=[r['path'] for r in entries],scenario='shared',export_kinematics=True,export_individual_photons=True,export_global_only=True))
            for scenario in ('nominal','candidate'):
                for sample in ('mc','background','pi0_mc'):
                    r=cm['samples'][period+'__'+sample];path=(run/r['file']).resolve()
                    if scenario=='nominal' and export_digest(path)!=r['sha256']:p.error('MC cache changed')
                    jobs.append(dict(period=period,sample=sample,cache=str(path),parameters=[0.]*20 if scenario=='nominal' else pars[period],scenario=scenario,export_kinematics=True,export_individual_photons=True,export_global_only=True))
        out.mkdir(parents=True);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));base.dump(out/'spec.json',spec);worker=base.compile_worker(out,reference);pool=base.Pool(worker,out/'spec.json',work,args.workers)
        print('[export] DATA once per channel/period; nominal and frozen MC from existing caches; no fits',flush=True)
        for _,job,path in pool.run('cut-arrays',jobs,'kinematic_arrays'):
            shape=json.loads(Path(str(path)+'.json').read_text())['shape'];records.append(dict(period=job['period'],sample=job['sample'],scenario=job['scenario'],path=str(path),shape=shape,sha256=export_digest(path)))
    for r in records:
        shape=r['shape'];path=Path(r['path'])
        if len(shape)!=2 or shape[1]!=28 or path.stat().st_size!=int(np.prod(shape))*8:p.error('Invalid individual-photon export')
        a=np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty((0,28))
        if r['scenario']=='shared':arrays[r['period']][r['sample']]=a
        else:arrays[r['period']][r['scenario']][r['sample']]=a
    for record in source['inventory']:
        if record['period'] in spec['periods'] and record['sample'] in ('data','pi0_data'):
            st=Path(record['path']).stat()
            if (st.st_size,st.st_mtime_ns)!=(record['size'],record['mtime_ns']):p.error('DATA input changed during or since export')
    out.mkdir(parents=True,exist_ok=True);base.dump(out/'export_manifest.json',dict(fingerprint=fingerprint,records=records));base.dump(out/'manifest.json',dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},fingerprint=fingerprint))
    rows=[];links=[];splitrows=[];coverage=[]
    for period,cases in arrays.items():
        cfg=next(p for p in cuts.production.PERIODS if p.key==period);prob={s:cm['samples'][period+'__'+s]['probability'] for s in ('mc','background','pi0_mc')};norm,selected,mult=get_scales(period,run,spec,cm,reference);scales={}
        for scenario,index in [('nominal',0),('candidate',1)]:
            for half in (0,1):
                ratio=1 if not index else selected[0,half].sum()/selected[1,half].sum()
                for top,topname in enumerate(base.TOPS):
                    beta=float(norm[topname]['signal_scale']);alpha=float(norm[topname]['background_scale'])*(mult if index else 1)*ratio
                    if index:
                        block=spec['exclusivity'][cuts.production.cut_json_key('dvcs',cfg,cuts.production.TOPOLOGIES[top])];counts={s:int((closure.mask(cases[scenario][s],top,block,'data')&(cases[scenario][s][:,11]==half)).sum())*2/prob[s] for s in ('mc','background')};d=spec['reference_flows'][period+'__data']['selected'][top];beta=(d-alpha*counts['background'])/counts['mc'] if counts['mc']>0 and d>alpha*counts['background'] else float('nan')
                    scales[scenario,half,top]=(beta,alpha)
        for sample,prefix in [('data','dvcs'),('pi0_data','eppi0')]:
            actual=[int(closure.mask(cases[sample],top,spec['exclusivity'][cuts.production.cut_json_key(prefix,cfg,cuts.production.TOPOLOGIES[top])],'data').sum()) for top in range(3)]
            expected=spec['reference_flows'][period+'__'+sample]['selected']
            if actual!=expected:raise ValueError(f'{period}/{sample}: DATA selection replay differs: {actual} versus {expected}')
        splits=choose_splits(cases['data'],cases['pi0_data'])
        for r in splits:r['period']=period
        splitrows+=splits
        for sample,a in [('data',cases['data']),('pi0_data',cases['pi0_data']),*[(scenario+'_'+s,a) for scenario,d in cases.items() if isinstance(d,dict) for s,a in d.items()]]:
            coverage.append(dict(period=period,sample=sample,global_events=len(a),photon1_unknown=int((~np.isfinite(a[:,20])).sum()),photon2_unknown=int((~np.isfinite(a[:,24])).sum())))
        rr,ll=analyze(period,cases,spec,scales,prob,splits,out,args.min_entries);rows+=rr;links+=ll;print('[done] '+period+' frozen low/high comparisons',flush=True)
    for name,rr in [('scores',rows),('splits',splitrows),('photon_coverage',coverage)]:base.table(out/(name+'.csv'),rr)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Frozen response kinematic checks</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Frozen joint response: low/high checks</h1><p>No response parameter or normalization is refitted per subdivision. Splits are medians of globally selected DATA half 0, before exclusivity and before any score calculation. The same numeric edge applies to all MC and both channels. Photon roles share one edge per detector and axis, determined from epgamma DATA. Photon 2 is shown only in the two-photon channel. The response parameters never depend on photon role. Momentum is GeV; angles are degrees. Check splits.csv and photon_coverage.csv first.</p><p>Photon 1 and photon 2 each use their own detector tag/coordinates. Their diagnostics overlap in events and must not be summed. Unresolved control photon coordinates are excluded only from photon-specific subdivisions, retained for electron/proton checks and explicitly counted. Unknown second-photon coordinates are expected for epgamma samples. The pi0 composite is never used as a detected photon.</p><p>Heatmaps show nominal-to-frozen-joint score improvement for MC half 1. Gray cells lack support; all numerical scores, core-width ratios and median shifts are in scores.csv for both halves. Compare candidate residuals between low/high, not just improvement: positive improvement alone does not establish agreement. The example plots select the largest supported candidate residual per channel/side after evaluation; they do not select training or split boundaries.</p><p>Each histogram uses full saved globals and the common DATA exclusivity windows, omitting the displayed cut. Epγ includes signal+leakage with the fitted period/topology normalization held fixed across subdivisions; epγγ uses AAOGEN. Curves are unit-area shapes. MC half validation shares DATA histograms, so not independent DATA validation. This is a descriptive check, not an acceptance or contamination correction. Momentum/angle/other kinematics remain correlated. No kinematic parameter dependence is introduced.</p><p>'+ ' · '.join('<a href="'+n+'.csv">'+n+'</a>' for n in ('scores','splits','photon_coverage'))+'</p>'+''.join('<p>'+html.escape(r['period']+' '+r['category']+' '+r['axis'])+': <a href="'+r['scores']+'">score improvement</a> · <a href="'+r['residuals']+'">candidate low/high residuals</a> · <a href="'+r['examples']+'">largest residual examples</a></p>' for r in links))
    print('[done] '+str(out/'REPORT.html')+'; arrays stay on scratch; transfer report only',flush=True)

if __name__=='__main__':main()
