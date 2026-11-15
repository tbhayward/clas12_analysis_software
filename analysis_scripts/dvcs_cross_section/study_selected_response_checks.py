#!/usr/bin/env python3
"""Three cache-only checks of the frozen selection-first particle response."""
from __future__ import annotations
import argparse,copy,html,json
from pathlib import Path
import numpy as np
import study_selected_particle_response as study
base=study.base
EDGES=np.array([0.,.11,.15,.25,.4,.6,.8,1.,np.inf])
LABELS=['<0.11','0.11–0.15','0.15–0.25','0.25–0.40','0.40–0.60','0.60–0.80','0.80–1.00','≥1.00','not selected']


def states(a,mask):
    result=np.full(len(a),8,dtype=int);ok=mask&np.isfinite(a[:,9])&(a[:,9]>=0)
    result[ok]=np.searchsorted(EDGES,a[ok,9],side='right')-1
    return result


def migration(a,b,old,new,group):
    # Export writes each supported cache event once, without a global-pass filter.
    if a.shape!=b.shape or not np.array_equal(a[:,[0,11]],b[:,[0,11]]):raise ValueError('Nominal/candidate event alignment changed')
    left=(a[:,0]==group)&study.cuts.window_mask(a,old,'mc');right=(b[:,0]==group)&study.cuts.window_mask(b,new,'mc')
    x,y=states(a,left),states(b,right);keep=(a[:,0]==group)&(left|right)
    matrix=np.zeros((9,9),dtype=np.int64);np.add.at(matrix,(x[keep],y[keep]),1)
    return matrix


def sample_histograms(a,windows,spec,control=False):
    groups=study.CONTROL_GROUPS if control else range(3);blocks=windows['control' if control else 'dvcs'];variables=spec['variables']+([spec['mass_variable']] if control else [])
    h=np.zeros((2,len(groups),len(variables),max(v['bins']+2 for v in variables)),dtype=np.int64);selected=np.zeros((2,len(groups)),dtype=np.int64)
    for g,raw in enumerate(groups):
        for half in (0,1):
            hm=a[:,11]==half;selected[half,g]=np.count_nonzero(hm&study.selected_mask(a,raw,blocks[g],'data'))
            for j,v in enumerate(variables):
                keep=(a[:,0]==raw)&(a[:,8]>.5)&hm
                for name,w in blocks[g]['data'].items():
                    if name==v['name']:continue
                    z=a[:,1+base.OBS.index(name)];keep&=np.isfinite(z)&(z<=float(w['cut_high']))
                    if w.get('mode')!='upper_quantile':keep&=z>=float(w['cut_low'])
                edges=base.diag.variable_edges(v);h[half,g,j,:len(edges)+1]=study.unordered.histogram(a,keep,28 if j==7 else j+1,edges)
    return h,selected


def photon_entries(a,mask,coordinate):
    values=[];masses=[]
    for col in (20,24):
        x=a[:,col if coordinate=='momentum' else col+1]*(1 if coordinate=='momentum' else 180/np.pi)
        good=mask&(a[:,col+3]==1)&np.isfinite(x)&np.isfinite(a[:,28])&(a[:,col]>0)
        values.append(x[good]);masses.append(a[good,28])
    return np.concatenate(values),np.concatenate(masses)


def load_exports(run):
    arrays={};manifests={}
    for stage in ('initial','candidate'):
        manifest=json.loads((run/(stage+'_exports.json')).read_text());manifests[stage]=manifest
        for r in manifest['records']:
            path=Path(r['path'])
            if not path.is_file():raise ValueError('Farm scratch export missing: '+str(path))
            if study.export_digest(path)!=r['sha256']:raise ValueError('Export checksum changed: '+str(path))
            arrays[stage,r['period'],r['sample']]=study.arrays_from_export(path)
    return arrays,manifests


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--response-run',type=Path,required=True);p.add_argument('--joint-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--min-entries',type=int,default=100);args=p.parse_args(argv)
    run=args.response_run.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a fresh output directory')
    if args.min_entries<10:p.error('min-entries must be at least 10')
    spec=json.loads((run/'source_spec.json').read_text());signature=json.loads((run/'signature.json').read_text());joint=args.joint_run.resolve()
    if base.digest(joint/'cache_manifest.json')!=signature['cache_manifest']:p.error('Joint cache manifest differs from the fitted run')
    cm=json.loads((joint/'cache_manifest.json').read_text());cache=cm['samples']
    for name,expected in signature['helpers'].items():
        if base.digest(base.SOURCE/name)!=expected:p.error('Saved scoring helper changed: '+name)
    if base.digest(base.SOURCE/'study_selected_particle_response.py')!=signature['pipeline']:p.error('Selection-first pipeline changed from fitted run')
    arrays,manifests=load_exports(run) # Validate everything before creating output.
    old=json.loads((run/'initial_fixed_windows.json').read_text());new=json.loads((run/'final_windows.json').read_text())
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(parents=True);links=[];migration_rows=[];migration_summary=[];scores=[];mass_rows=[];score_summary=[]
    for period in spec['periods']:
        folder=out/'plots'/period;folder.mkdir(parents=True)
        initial={s:arrays['initial',period,s] for s in ('data','pi0_data','mc','background','pi0_mc')};candidate={s:initial[s] if s in ('data','pi0_data') else arrays['candidate',period,s] for s in initial}
        for sample,scope in [(s,k) for s in ('mc','background','pi0_mc') for k in ('response_only','response_and_new_cuts')]:
            control=sample=='pi0_mc';groups=study.CONTROL_GROUPS if control else range(3);labels=study.CONTROL_LABELS if control else base.TOPS;channel='control' if control else 'dvcs'
            fig,axes=plt.subplots(1,len(labels),figsize=(6*len(labels),5),layout='constrained')
            for g,(raw,label) in enumerate(zip(groups,labels)):
                matrix=migration(initial[sample],candidate[sample],old[period][channel][g],(old if scope=='response_only' else new)[period][channel][g],raw)
                for i in range(8):migration_summary.append(dict(period=period,sample=sample,scope=scope,group=label,t_bin=LABELS[i],nominal_selected=int(matrix[i].sum()),modified_selected=int(matrix[:,i].sum()),same_bin=int(matrix[i,i]),enter_from_unselected=int(matrix[8,i]),exit_selection=int(matrix[i,8]),enter_from_other_t=int(matrix[:8,i].sum()-matrix[i,i]),leave_for_other_t=int(matrix[i,:8].sum()-matrix[i,i])))
                for i in range(9):
                    for j in range(9):migration_rows.append(dict(period=period,sample=sample,scope=scope,group=label,old_bin=LABELS[i],new_bin=LABELS[j],events=int(matrix[i,j])))
                # Account for both bin migration and selection change, including outside range.
                denom=matrix.sum(axis=1,keepdims=True);fraction=np.divide(100.*matrix,denom,out=np.zeros_like(matrix,dtype=float),where=denom>0);ax=axes[g];im=ax.imshow(fraction,origin='lower',vmin=0,vmax=100,cmap='Blues');ax.set(xticks=range(9),yticks=range(9),xticklabels=LABELS,yticklabels=LABELS,xlabel='Modified |t| / '+('initial' if scope=='response_only' else 'final')+' MC selection',ylabel='Nominal |t| / initial MC selection',title=label);ax.tick_params(axis='x',rotation=60,labelsize=7);ax.tick_params(axis='y',labelsize=7)
                for i,j in zip(*np.where(fraction>=1)):ax.text(j,i,f'{fraction[i,j]:.0f}',ha='center',va='center',fontsize=7,color='white' if fraction[i,j]>55 else 'black')
            fig.colorbar(im,ax=list(axes),label='Percent of original row');fig.suptitle(period+' · '+sample+' · '+scope+' · matched events');path=folder/(sample+'_'+scope+'_t_migration.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
        local={s:cache[period+'__'+s] for s in ('mc','background','pi0_mc')};mult=json.loads((run/'response'/period/'candidate_parameters.json').read_text())['background_multiplier']
        # Same final DATA windows, actual two MC halves, nominal variance/support.
        data,counts=study.data_histograms(period,initial,new[period],spec);hist={};selected={}
        for sample in local:
            pairs=[sample_histograms(z[sample],new[period],spec,sample=='pi0_mc') for z in (initial,candidate)]
            hist[sample]=np.stack([x[0] for x in pairs]);selected[sample]=np.stack([x[1] for x in pairs])
        engine=study.Engine(None,spec,local,data,counts,args.min_entries,0)
        _,nominal=engine.score_one(hist,selected,0,mult)
        engine.support={(r['channel'],(base.TOPS if r['channel']=='epgamma' else study.CONTROL_LABELS).index(r['group']),7 if r['variable']=='m_gamma_gamma' else base.OBS.index(r['variable']),r['mc_half']) for r in nominal if min(r['data_count'],r['raw_nominal_mc'])>=args.min_entries and np.isfinite(r['score'])}
        period_rows=[]
        for i,scenario in enumerate(('nominal','candidate')):
            _,rr=engine.score_one(hist,selected,i,mult);period_rows.extend(dict(period=period,scenario=scenario,**r) for r in rr)
        scores+=period_rows
        fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
        for ax,ch in zip(axes,('epgamma','epgammagamma')):
            names=[*study.TRAIN,*([7] if ch=='epgammagamma' else []),6] # CM/opening angle deliberately absent from headline.
            for i,(scenario,color) in enumerate((('nominal','#888888'),('candidate','#22946b'))):
                values=[]
                for j in names:
                    name='m_gamma_gamma' if j==7 else base.OBS[j];terms=[r['score'] for r in period_rows if r['scenario']==scenario and r['channel']==ch and r['variable']==name and r['mc_half']==1];value=float(np.mean(terms)) if terms else float('nan');values.append(value);score_summary.append(dict(period=period,channel=ch,variable=name,scenario=scenario,mc_half=1,groups=len(terms),mean_score=value))
                ax.bar(np.arange(len(names))+(i-.5)*.36,values,.36,label=scenario,color=color)
            ax.set(xticks=range(len(names)),xticklabels=['mγγ' if j==7 else base.OBS[j] for j in names],ylabel='Other-MC-half shape score (lower is better)',title=ch);ax.tick_params(axis='x',rotation=35);ax.legend()
        fig.suptitle(period+' · common final DATA windows, N−1; frozen response and nuisance');path=folder/'final_window_scores.png';fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
        for raw,label in zip(study.CONTROL_GROUPS,study.CONTROL_LABELS):
            g=study.CONTROL_GROUPS.index(raw);win=new[period]['control'][g]
            masks={s:study.selected_mask(a,raw,win,'data') for s,a in [('DATA',initial['pi0_data']),('nominal',initial['pi0_mc']),('candidate',candidate['pi0_mc'])]};sources={'DATA':initial['pi0_data'],'nominal':initial['pi0_mc'],'candidate':candidate['pi0_mc']}
            fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained');any_data=False
            for j,coordinate in enumerate(('momentum','angle')):
                dx,dm=photon_entries(sources['DATA'],masks['DATA'],coordinate)
                if len(dx)<args.min_entries:continue
                edges=np.unique(np.quantile(dx,[0,.25,.5,.75,1.]));edges[0]=-np.inf;edges[-1]=np.inf;centers=np.arange(len(edges)-1);any_data=True
                for scenario,color in [('DATA','#205cab'),('nominal','#888888'),('candidate','#22946b')]:
                    x,m=photon_entries(sources[scenario],masks[scenario],coordinate);med=[];width=[]
                    for b,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
                        v=1000*m[(x>=lo)&(x<hi)];q=np.quantile(v,[.16,.5,.84]) if len(v)>=args.min_entries else [np.nan]*3;med.append(q[1]);width.append((q[2]-q[0])/2);mass_rows.append(dict(period=period,group=label,coordinate=coordinate,subdivision=b,low=float(lo),high=float(hi),sample=scenario,photon_entries=len(v),median_mass_MeV=q[1],core_half_width_MeV=(q[2]-q[0])/2))
                    axes[0,j].plot(centers,med,'o-',label=scenario,color=color);axes[1,j].plot(centers,width,'o-',label=scenario,color=color)
                for k in range(2):axes[k,j].set(xticks=centers,xticklabels=[f'{a:.2g}–{b:.2g}' for a,b in zip(edges[:-1],edges[1:])],xlabel='FD photon '+coordinate+(' (GeV)' if coordinate=='momentum' else ' (degrees)'),ylabel='Median mass (MeV)' if k==0 else '16–84% half width (MeV)');axes[k,j].legend();axes[k,j].grid(alpha=.2)
            if any_data:
                fig.suptitle(period+' · '+label+' · final common DATA cuts\nBoth FD photon roles pooled identically; a two-FD event contributes twice');path=folder/(label+'_mass_kinematics.png');fig.savefig(path,dpi=110);links.append(path.relative_to(out).as_posix())
            plt.close(fig)
        print('[check] '+period+' complete; arrays only, no refit',flush=True)
    base.table(out/'t_migration_summary.csv',migration_summary);base.table(out/'final_window_score_summary.csv',score_summary);base.table(out/'t_migration.csv',migration_rows);base.table(out/'final_window_scores.csv',scores);base.table(out/'mass_kinematics.csv',mass_rows)
    base.dump(out/'provenance.json',dict(response_run=str(run),signature=base.digest(run/'signature.json'),initial_exports=base.digest(run/'initial_exports.json'),candidate_exports=base.digest(run/'candidate_exports.json'),joint_cache=base.digest(joint/'cache_manifest.json'),script=base.digest(Path(__file__))))
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Frozen response checks</title><style>body{font:17px system-ui;max-width:1200px;margin:30px auto}img{max-width:100%}</style><h1>Frozen response: three checks</h1><p>No ROOT reads, fits, cut changes or production modifications. Migration separately shows response alone under fixed initial MC cuts and response plus final redetermined MC cuts, with outside-|t| and selection failure explicit. Matrices exclude events failing both selections. Counts are raw sampled MC.</p><p>Shape scores compare nominal and modified MC under identical final DATA windows (N−1), with fixed nominal support/variance and the fitted nuisance held fixed for both. Each scenario recalculates control/signal area normalizations from its selected population; this is a shape test, not absolute normalization. Both actual MC halves are tabulated, using the same DATA. Headline bars prioritize missing quantities and balance; angle projections remain in the CSV.</p><p>Mass summaries use the full final common selection and shared DATA-defined FD-photon quartiles, pooling both roles. Two-FD events enter twice, so photon entries are correlated; no independent-event errors are claimed. Unresolved photons are excluded from kinematic subdivisions; mass is the saved-calibration-preserving exported mass. Photon coordinates are reconstructed, not truth. DATA/MC kinematic population and upstream mass-window restrictions still matter.</p><p><a href="t_migration.csv">Migration counts</a> · <a href="final_window_scores.csv">Final-window scores</a> · <a href="mass_kinematics.csv">Mass subdivisions</a></p>'+''.join('<p>'+html.escape(x)+'<br><img loading="lazy" src="'+x+'"></p>' for x in links));print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
