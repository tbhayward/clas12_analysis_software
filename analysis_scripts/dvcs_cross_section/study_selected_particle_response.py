#!/usr/bin/env python3
"""Legacy DATA cuts -> fixed-selection joint response -> redetermined cuts/counts."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,copy,csv,html,json,shlex,shutil,subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts
import study_unordered_control_response as unordered
from study_particle_common_cut_stability import export_digest

CONTROL_GROUPS=[0,3,4,5]
CONTROL_LABELS=[unordered.GROUPS[g] for g in CONTROL_GROUPS]
TRAIN=(0,3,4,5);RESERVED=(2,6)


def compile_worker(out,reference,work):
    build=out/'backend';build.mkdir();temporary=work/'unordered_sources';unordered.isolated_sources(base.SOURCE,temporary)
    for name in ('global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h'):shutil.copy2(reference/'backend'/name,build/name)
    for name in unordered.SOURCES:shutil.copy2(base.SOURCE/name,build/name)
    worker_source=build/'particle_response_worker.cpp';worker_source.write_text(worker_source.read_text().split('\nint main(',1)[0]+'\n')
    p=build/'exclusivity_resolution_worker.cpp';p.write_text(p.read_text().split('\nint main(',1)[0]+'\n')
    text=(temporary/'particle_response_cut_arrays.h').read_text().replace('diagnostic_photon_coordinates','selected_photon_coordinates').replace('export_cut_arrays','export_unordered_control_arrays')
    old='auto shifted=vary_control(c,pars,beam,base,noise,mass);emit(shifted,c.base.key,photons?&pp:nullptr,mass);'
    new='''auto shifted=vary_control(c,pars,beam,base,noise,mass);
    bool active=std::any_of(pars.begin(),pars.end(),[](double x){return x!=0;});
    if(active&&c.photons[9]>.5){auto p1=particle(c.photons[0],c.photons[1],c.photons[2],0),p2=particle(c.photons[3],c.photons[4],c.photons[5],0);mass=c.photons[8]+mass-std::sqrt(std::max(0.,(p1+p2).m2()));}
    emit(shifted,c.base.key,photons?&pp:nullptr,mass);'''
    if text.count(old)!=1:raise ValueError('Unordered exporter changed; review selected mass preservation')
    (build/'unordered_control_arrays.h').write_text(text.replace(old,new));shutil.copy2(Path(__file__).with_name('selected_response_worker.cpp'),build/'selected_response_worker.cpp')
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True));exe=build/'selected_response_worker'
    subprocess.run([os.environ.get('CXX','c++'),*flags,'-O3','-std=c++17','-I'+str(build),str(build/'selected_response_worker.cpp'),str(build/'global_cuts.cpp'),'-o',str(exe),*libs],check=True);return exe


def arrays_from_export(path):
    shape=json.loads(Path(str(path)+'.json').read_text())['shape']
    if len(shape)!=2 or shape[1] not in (12,29) or Path(path).stat().st_size!=int(np.prod(shape))*8:raise ValueError('Invalid event export')
    return np.memmap(path,dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty(shape)


def packed_arrays(a,group):
    keep=(a[:,0]==group)&(a[:,8]>.5)
    return {name:np.asarray(a[keep,j+1],dtype=np.float32) for j,name in enumerate(base.OBS)}


def production_cut_fit(period,a,folder):
    """Reuse the same fitter/cut algorithms; unordered control groups get own windows."""
    p=next(p for p in cuts.production.PERIODS if p.key==period);containments=dict(nominal=.95,tight=.90,loose=.98);dvcs=[];orders=[];control=[];folder.mkdir(parents=True,exist_ok=True)
    for top,topology in enumerate(cuts.production.TOPOLOGIES):
        group=[0,3,5][top];ctrl={s:packed_arrays(a[k],group) for s,k in [('data','pi0_data'),('signal','pi0_mc')]}
        epg={s:packed_arrays(a[k],top) for s,k in [('data','data'),('signal','mc'),('background','background')]}
        if any(len(v['Delta_phi'])<100 for v in [*epg.values(),*ctrl.values()]):raise ValueError(f'{period}/{topology.key}: insufficient cut-development sample')
        calibration=cuts.calibration(ctrl['data'],ctrl['signal'],topology)
        steps,blocks,hist=cuts.production.run_dvcs_iterative_cuts(p,topology,epg,12.,20.,100,['Delta_phi','theta_gamma_gamma','pTmiss'],4.,8.,True,.90,.95,.90,.95,calibration,containments,.25,25.)
        key=cuts.production.cut_json_key('dvcs',p,topology);dvcs.append(blocks['nominal']);orders.append([s.variable for s in steps if s.accepted_for_production]);dest=folder/topology.key;dest.mkdir(exist_ok=True);base.dump(dest/'cuts.json',blocks);cuts.production.write_results_csv(dest/'cut_flow.csv',cuts.production.cumulative_optimization_rows(p,topology,steps));cuts.production.draw_serialized(cuts.production.draw_iterative_cut_canvas,dest/'cut_development.png',p,topology,'dvcs',steps,False,100,final_histograms=hist)
    for index,raw in enumerate(CONTROL_GROUPS):
        policy=[0,1,1,2][index];topology=cuts.production.TOPOLOGIES[policy];ctrl={s:packed_arrays(a[k],raw) for s,k in [('data','pi0_data'),('signal','pi0_mc')]}
        if min(len(v['Delta_phi']) for v in ctrl.values())<100:raise ValueError(f'{period}/{CONTROL_LABELS[index]}: insufficient control cut sample')
        # Mixed photons use the CD-FD optimizer policy, but derive independent windows.
        steps,blocks,hist=cuts.production.run_pi0_iterative_cuts(p,topology,ctrl,orders[policy],12.,20.,100,4.,8.,True,.90,containments)
        control.append(blocks['nominal']);dest=folder/CONTROL_LABELS[index];dest.mkdir(exist_ok=True);base.dump(dest/'cuts.json',blocks);cuts.production.draw_serialized(cuts.production.draw_iterative_cut_canvas,dest/'cut_development.png',p,topology,'pi0',steps,False,100,final_histograms=hist)
    result=dict(dvcs=dvcs,control=control,orders=orders,policy='legacy 95% template cuts; mixed control uses CD-FD optimizer policy with independent distributions')
    base.dump(folder/'cuts_95.json',result);return result


def selected_mask(a,group,window,kind):
    return (a[:,0]==group)&cuts.window_mask(a,window,kind)


def data_histograms(period,a,windows,spec):
    result={};selected={};control=a['pi0_data'];mass=dict(name='m_gamma_gamma',label='Two-photon mass (GeV)',bins=140,low=.10,high=.17)
    for channel,source,groups,blocks,variables in [('epgamma',a['data'],list(range(3)),windows['dvcs'],spec['variables']),('epgammagamma',control,CONTROL_GROUPS,windows['control'],[*spec['variables'],mass])]:
        maxbins=max(v['bins']+2 for v in variables);hist=np.zeros((len(groups),len(variables),maxbins));sel=np.zeros(len(groups),dtype=int)
        for g,raw in enumerate(groups):
            block=blocks[g];sel[g]=int(selected_mask(source,raw,block,'data').sum())
            for j,v in enumerate(variables):
                keep=(source[:,0]==raw)&(source[:,8]>.5)
                for name,w in block['data'].items():
                    if name==v['name']:continue
                    values=source[:,1+base.OBS.index(name)];keep&=np.isfinite(values)&(values<=float(w['cut_high']))
                    if w.get('mode')!='upper_quantile':keep&=values>=float(w['cut_low'])
                edges=base.diag.variable_edges(v);col=28 if j==7 else j+1;hist[g,j,:len(edges)+1]=unordered.histogram(source,keep,col,edges)
        result[channel]=hist;selected[channel]=sel
    if min(selected['epgamma'])<=0 or selected['epgammagamma'].sum()<=0:raise ValueError('Empty initial DATA reference')
    return result,selected


class Engine:
    def __init__(self,pool,spec,cache,data,counts,min_mc,transfer_uncertainty):
        self.pool=pool;self.spec=spec;self.cache=cache;self.data=data;self.counts=counts;self.min_mc=min_mc;self.uncertainty=transfer_uncertainty;self.calls=0;self.support=None;self.history=[]
    def scales(self,selected,i,half,multiplier):
        ctrl=float(selected['pi0_mc'][i,half].sum())*2/self.cache['pi0_mc']['probability'];nd=float(self.counts['epgammagamma'].sum())
        if ctrl<=0:return float('nan'),np.full(3,np.nan)
        alpha=nd/ctrl*float(self.spec.get('control_to_leakage_exposure_ratio',1.))*multiplier
        sig=selected['mc'][i,half]*2/self.cache['mc']['probability'];bkg=selected['background'][i,half]*2/self.cache['background']['probability'];res=self.counts['epgamma']-alpha*bkg
        beta=np.divide(res,sig,out=np.full(3,np.nan),where=(sig>0)&(res>0));return alpha,beta
    def prediction(self,hist,selected,i,half,g,j,multiplier):
        alpha,beta=self.scales(selected,i,half,multiplier);n=self.spec['variables'][j]['bins']+2
        weights=[2*beta[g]/self.cache['mc']['probability'],2*alpha/self.cache['background']['probability']];h=[hist[s][i,half,g,j,:n].astype(float) for s in ('mc','background')]
        return weights[0]*h[0]+weights[1]*h[1],weights[0]**2*h[0]+weights[1]**2*h[1]
    def score_one(self,hist,selected,i,multiplier):
        rows=[];buckets={(ch,kind):[] for ch in ('epgamma','epgammagamma') for kind in ('train','validation_mc','reserved')}
        for ch,sample,groups,vars in [('epgamma','mc',3,(*TRAIN,*RESERVED)),('epgammagamma','pi0_mc',4,(*TRAIN,7,*RESERVED))]:
            for g in range(groups):
                for j in vars:
                    n=(self.spec['mass_variable'] if j==7 else self.spec['variables'][j])['bins']+2;d=self.data[ch][g,j,:n]
                    for half in (0,1):
                        support=(ch,g,j,half)
                        if self.support is not None and support not in self.support:continue
                        if ch=='epgamma':m,var=self.prediction(hist,selected,i,half,g,j,multiplier);ref,refvar=self.prediction(hist,selected,0,half,g,j,1.);raw0=sum(hist[s][0,half,g,j,:n].sum() for s in ('mc','background'))
                        else:m=hist[sample][i,half,g,j,:n].astype(float);var=m.copy();ref=hist[sample][0,half,g,j,:n].astype(float);refvar=ref.copy();raw0=ref.sum()
                        score,bins=base.shape_score(d,m,var,int(raw0),self.min_mc,ref,refvar);kind='reserved' if j in RESERVED else 'validation_mc' if half else 'train'
                        buckets[ch,kind].append(score);rows.append(dict(channel=ch,group=base.TOPS[g] if ch=='epgamma' else CONTROL_LABELS[g],variable='m_gamma_gamma' if j==7 else base.OBS[j],mc_half=half,kind=kind,score=score,usable_bins=bins,data_count=int(d.sum()),raw_nominal_mc=int(raw0)))
        result={}
        for ch in ('epgamma','epgammagamma'):
            for kind in ('train','validation_mc','reserved'):
                values=buckets[ch,kind];result[ch+'_'+kind]=float(np.mean(values)) if values and np.all(np.isfinite(values)) else float('inf')
        for kind in ('train','validation_mc','reserved'):result[kind]=.5*(result['epgamma_'+kind]+result['epgammagamma_'+kind])
        result['train_shape']=result['train'];result['train']+=((multiplier-1)/self.uncertainty)**2 if self.uncertainty>0 else 0.;result['background_multiplier']=multiplier;result['background_constraint_sigma']=self.uncertainty
        return result,rows
    def evaluate(self,parameters,label,detailed=False):
        unique=[tuple([0.]*20)]
        for p in parameters:
            if len(p)!=20 or not np.all(np.isfinite(p)):raise ValueError('Invalid response vector')
            if tuple(p) not in unique:unique.append(tuple(p))
        jobs=[dict(period=self.spec['periods'][0],sample=s,cache=self.cache[s]['path'],parameters=unique) for s in ('mc','background','pi0_mc')];hist={};selected={}
        for _,job,path in self.pool.run('evaluate-selected',jobs,f'eval_{self.calls:03d}_{label}'):
            meta=json.loads(Path(str(path)+'.json').read_text());shape=meta['shape'];hist[job['sample']]=np.memmap(path,dtype='<i8',mode='r',shape=tuple(shape));selected[job['sample']]=np.asarray(meta['selected']).reshape(len(unique),2,shape[2])
        self.calls+=1
        if self.support is None:
            # Fixed DATA/nominal support; never let a candidate drop difficult bins.
            ctrl_mc=selected['pi0_mc'][0].sum()/self.cache['pi0_mc']['probability'];ctrl_data=self.counts['epgammagamma'].sum()
            if min(ctrl_mc,ctrl_data)<=0:raise ValueError('Empty nominal control normalization')
            self.uncertainty=float(np.hypot(self.uncertainty,np.sqrt(1/ctrl_mc+1/ctrl_data)))
            for half in (0,1):
                _,beta=self.scales(selected,0,half,1.)
                if not np.all(np.isfinite(beta)):raise ValueError('Nominal control transfer leaves nonpositive signal; inspect initial selections')
            _,rows=self.score_one(hist,selected,0,1.);base.table(self.pool.spec.parent/'nominal_objective_projections.csv',rows);self.support=set()
            for r in rows:
                if r['data_count']>=self.min_mc and r['raw_nominal_mc']>=self.min_mc and np.isfinite(r['score']):self.support.add((r['channel'],(base.TOPS if r['channel']=='epgamma' else CONTROL_LABELS).index(r['group']),7 if r['variable']=='m_gamma_gamma' else base.OBS.index(r['variable']),r['mc_half']))
            base.table(self.pool.spec.parent/'objective_support.csv',[dict(channel=x[0],group=x[1],variable=x[2],half=x[3]) for x in sorted(self.support)])
            if any(not any(x[0]==ch and x[3]==0 and x[2] in (*TRAIN,7) for x in self.support) for ch in ('epgamma','epgammagamma')):raise ValueError('Empty channel training support; inspect nominal_objective_projections.csv')
        scores=[];records=[];multipliers=sorted(set(max(.01,1+x*self.uncertainty) for x in (-2,-1,0,1,2))) if self.uncertainty else [1.]
        for i,p in enumerate(unique):
            trials=[self.score_one(hist,selected,i,m) for m in multipliers];score,rows=min(trials,key=lambda x:(x[0]['train'],abs(x[0]['background_multiplier']-1)));scores.append(score);records.extend(dict(r,scenario=i) for r in rows);self.history.append(dict(call=label,parameters=json.dumps(p),**score))
        return unique,scores,records,hist,selected


def select_training(parameters,scores,allowed):
    eligible=[i for i,p in enumerate(parameters) if tuple(p) in {tuple(a) for a in allowed} and np.isfinite(scores[i]['train'])]
    if not eligible:raise ValueError('No finite training candidate; inspect signal/background normalization')
    return min(eligible,key=lambda i:(scores[i]['train'],sum(abs(x) for x in parameters[i])))


def fit_response(engine,start,folder,rounds,scan_range):
    best=list(start);history=[];grid=base.parameter_grid(base.TOPS,'fit',scan_range);levels=copy.deepcopy(grid)
    for k,v in enumerate(start):
        if not grid[k][0]<=v<=grid[k][-1]:raise ValueError('Starting response outside bounds')
    for iteration in range(rounds):
        for k in grid:
            ll=sorted(set(levels[k]+[best[k]]));lo=max((v for v in ll if v<best[k]),default=best[k]);hi=min((v for v in ll if v>best[k]),default=best[k]);values=sorted(set([lo,(lo+best[k])/2,best[k],(hi+best[k])/2,hi]));trials=[]
            for value in values:p=list(best);p[k]=value;trials.append(p)
            pars,scores,_,hist,_=engine.evaluate(trials,f'round_{iteration}_{k}');chosen=select_training(pars,scores,trials);best=list(pars[chosen]);history.append(dict(round=iteration,parameter=base.PARAMS[k],value=best[k],**scores[chosen]));levels[k]=sorted(set(ll+values));base.discard_histograms(hist)
    pars,scores,records,hist,selected=engine.evaluate([[0.]*20,start,best],'final');base.table(folder/'fit_steps.csv',history);base.table(folder/'projection_scores.csv',records);base.table(folder/'evaluation_history.csv',engine.history)
    result=dict(period=engine.spec['periods'][0],parameters=dict(zip(base.PARAMS,best)),bias_mode='fit',fit_topologies=base.TOPS,joint_control_fit=True,selection_first=True,unordered_control=True,training_variables=[base.OBS[j] for j in TRAIN],control_mass_trained=True,reserved_variables=[base.OBS[j] for j in RESERVED],background_multiplier=scores[pars.index(tuple(best))]['background_multiplier'],scan_range=scan_range,rounds=rounds)
    base.dump(folder/'candidate_parameters.json',result);return best,pars,scores,records,hist,selected


def render_fit(period,engine,pars,scores,hist,selected,start,best,out,final_selection=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=[('fitted response',best,'#22946b')] if final_selection else [('initial nominal',[0.]*20,'#888888'),('starting response',start,'#d45c18'),('fitted response',best,'#22946b')];indices=[pars.index(tuple(p)) for _,p,_ in labels];links=[]
    folder=out/'plots'/period;folder.mkdir(parents=True,exist_ok=True)
    if not final_selection:
        fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
        for ax,ch in zip(axes,('epgamma','epgammagamma')):
            for j,(label,_,color) in enumerate(labels):ax.bar(np.arange(3)+(j-1)*.25,[float(scores[indices[j]][ch+'_'+k]) if scores[indices[j]][ch+'_'+k] is not None else float('nan') for k in ('train','validation_mc','reserved')],.25,label=label,color=color)
            ax.set(xticks=range(3),xticklabels=['Training','Other MC half','Reserved'],ylabel='Shape score (lower is better)',title=ch);ax.legend(fontsize=8)
        path=folder/'scores.png';fig.suptitle(period+' · fixed initial DATA cuts; joint selection-first fit');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    for ch,groups in [('epgamma',base.TOPS),('epgammagamma',CONTROL_LABELS)]:
        for g,group in enumerate(groups):
            fig,axes=plt.subplots(3,3,figsize=(14,10),layout='constrained');variables=engine.spec['variables']+([engine.spec['mass_variable']] if ch=='epgammagamma' else [])
            for j,v in enumerate(variables):
                ax=axes.flat[j];edges=base.diag.variable_edges(v);n=v['bins']+2;d=engine.data[ch][g,j,:n];curves=[('DATA',d,'#205cab')]
                for (label,_,color),i in zip(labels,indices):
                    m=hist['pi0_mc'][i,:,g,j,:n].sum(axis=0) if ch=='epgammagamma' else sum(engine.prediction(hist,selected,i,h,g,j,scores[i]['background_multiplier'])[0] for h in (0,1));curves.append((label,m,color))
                    if not np.all(np.isfinite(m)):ax.text(.02,.95,'Invalid MC normalization',transform=ax.transAxes,va='top',fontsize=8)
                for label,h,color in curves:
                    if h.sum()>0:ax.stairs(h[1:-1]/h.sum()/np.diff(edges),edges,label=label,color=color)
                if v['name']=='Mx2':ax.set_xlim(*base.mass_plot_limits(edges,[h for _,h,_ in curves]))
                elif j!=7:
                    q=[base.diag.quantile(h[1:-1],edges,x) for _,h,_ in curves if h.sum()>0 for x in (.01,.99)];q=[x for x in q if np.isfinite(x)]
                    if q:lo,hi=min(q),max(q);pad=max(.001,(hi-lo)*.08);ax.set_xlim(max(edges[0],lo-pad),min(edges[-1],hi+pad))
                ax.set(xlabel=v['label'],ylabel='Unit-area density');ax.legend(fontsize=7);ax.grid(alpha=.2)
            for ax in axes.flat[len(variables):]:ax.axis('off')
            fig.suptitle(period+' · '+ch+' · '+group+('\nFinal DATA windows on all samples, displayed cut omitted; fitted parameters fixed' if final_selection else '\nFrozen initial DATA selection, displayed cut omitted; response applied before selection'));path=folder/(ch+'_'+group+'.png');fig.savefig(path,dpi=100);plt.close(fig);links.append(path.relative_to(out).as_posix())
    return links


def selection_summary(period,a,windows,stage):
    rows=[]
    for sample in ('data','mc','background','pi0_data','pi0_mc'):
        control=sample.startswith('pi0');groups=CONTROL_GROUPS if control else range(3);labels=CONTROL_LABELS if control else base.TOPS;blocks=windows['control' if control else 'dvcs'];kind='data' if sample in ('data','pi0_data') else 'mc';z=a[sample]
        for g,raw in enumerate(groups):
            global_mask=(z[:,0]==raw)&(z[:,8]>.5);selected=selected_mask(z,raw,blocks[g],kind);n=int(global_mask.sum());rows.append(dict(stage=stage,period=period,sample=sample,group=labels[g],global_count=n,selected_count=int(selected.sum()),retained_fraction=float(selected.sum()/n) if n else float('nan')))
    return rows


def final_comparison(period,nominal,candidate,old,new,prob,exposure):
    rows=[];transfer=[]
    for lo,hi in [(0.,1.),(.11,.15),(.15,.25),(.25,.4),(.4,.6),(.6,.8),(.8,1.)]:
        controls={};leak={};signal={};data={}
        for scenario,a,windows in [('nominal',nominal,old),('candidate',candidate,new)]:
            ctrl_counts={}
            for sample,kind in [('pi0_data','data'),('pi0_mc','mc')]:
                total=0
                for g,raw in enumerate(CONTROL_GROUPS):
                    z=a[sample];mask=selected_mask(z,raw,windows['control'][g],kind)&(z[:,9]>=lo)&(z[:,9]<hi);total+=int(mask.sum())
                ctrl_counts[sample]=total if sample=='pi0_data' else total/prob['pi0_mc']
            controls[scenario]=ctrl_counts;alpha=ctrl_counts['pi0_data']/ctrl_counts['pi0_mc']*exposure if min(ctrl_counts['pi0_mc'],ctrl_counts['pi0_data'])>0 else float('nan');leak[scenario]=0;signal[scenario]=0;data[scenario]=0
            for top,label in enumerate(base.TOPS):
                z=a['mc'];mask=selected_mask(z,top,windows['dvcs'][top],'mc')&(z[:,9]>=lo)&(z[:,9]<hi);signal[scenario]+=mask.sum()/prob['mc']
                d=a['data'];dm=selected_mask(d,top,windows['dvcs'][top],'data')&(d[:,9]>=lo)&(d[:,9]<hi);data[scenario]+=dm.sum();b=a['background'];bm=selected_mask(b,top,windows['dvcs'][top],'mc')&(b[:,9]>=lo)&(b[:,9]<hi);leak[scenario]+=bm.sum()/prob['background']*alpha
                if scenario=='candidate':
                    b0=nominal['mc'];oldmask=selected_mask(b0,top,old['dvcs'][top],'mc')&(b0[:,9]>=lo)&(b0[:,9]<hi);ratio,error=cuts.ratio_error(oldmask,mask);rows.append(dict(period=period,topology=label,t_low=lo,t_high=hi,nominal_selected=int(oldmask.sum()),candidate_selected=int(mask.sum()),relative_acceptance=ratio,paired_mc_counting_error=error,inverse_acceptance=1/ratio if ratio>0 else float('nan')))
            transfer.append(dict(period=period,scenario=scenario,t_low=lo,t_high=hi,control_data=ctrl_counts['pi0_data'],control_mc_estimate=ctrl_counts['pi0_mc'],control_scale=alpha,epgamma_data=data[scenario],predicted_leakage=leak[scenario],subtracted_signal=data[scenario]-leak[scenario],signal_mc_estimate=signal[scenario],status='supported' if min(ctrl_counts['pi0_data'],ctrl_counts['pi0_mc']*prob['pi0_mc'])>=100 else 'insufficient control statistics',warning='coarse raw-count diagnostic; no event current/proton weights, no xB/Q2/phi binning; not a production cross section'))
    return rows,transfer


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--joint-run',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--workers',type=int,default=10);p.add_argument('--rounds',type=int,default=1);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response'));p.add_argument('--min-mc',type=int,default=100);p.add_argument('--resume',action='store_true');args=p.parse_args(argv)
    if not 1<=args.workers<=10 or args.rounds<1 or args.min_mc<10:p.error('Invalid workers, rounds or MC support')
    run=args.joint_run.resolve();out=args.output.resolve();meta=json.loads((run/'manifest.json').read_text());spec=copy.deepcopy(json.loads((run/'spec.json').read_text()));cm=json.loads((run/'cache_manifest.json').read_text());reference=Path(meta['arguments']['reference_study']);ref=json.loads((reference/'manifest.json').read_text())
    if cm['reference_sha256']!=base.digest(reference/'manifest.json') or meta['reference_manifest_sha256']!=base.digest(reference/'manifest.json'):p.error('Reference provenance changed')
    for name in ('particle_response_worker.cpp','particle_response_core.h','particle_control_core.h','particle_control_worker.h'):
        if base.digest(base.SOURCE/name)!=meta['source_hashes'][name]:p.error('Starting response evaluator changed: '+name)
    cache={k:dict(r,path=str((run/r['file']).resolve())) for k,r in cm['samples'].items()}
    for r in cache.values():
        if export_digest(r['path'])!=r['sha256']:p.error('MC cache changed')
    for k,n in [('global_cpp','global_cuts.cpp'),('global_header','global_cuts.h'),('domain_header','kinematic_bin_integration.h')]:
        if base.digest(reference/'backend'/n)!=spec['input_hashes'][k]:p.error('Saved globals changed')
    for r in ref['inventory']:
        if r['sample'] in ('data','pi0_data'):
            st=Path(r['path']).stat()
            if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):p.error('DATA input changed')
    starts={period:cuts.load_candidate_parameters(run/period/'candidate_parameters.json') for period in spec['periods']}
    signature=dict(joint_manifest=base.digest(run/'manifest.json'),cache_manifest=base.digest(run/'cache_manifest.json'),candidates={period:base.digest(run/period/'candidate_parameters.json') for period in spec['periods']},rounds=args.rounds,min_mc=args.min_mc,sources={n:base.digest(base.SOURCE/n) for n in unordered.SOURCES},pipeline=base.digest(Path(__file__)),worker=base.digest(Path(__file__).with_name('selected_response_worker.cpp')),fitter=base.digest(Path(cuts.production.__file__)),helpers={name:base.digest(base.SOURCE/name) for name in ('study_unordered_control_response.py','study_response_cut_redetermination.py','study_particle_response.py')})
    if out.exists():
        if not args.resume:p.error('Output exists; use --resume for the same inputs/settings')
        if json.loads((out/'signature.json').read_text())!=signature:p.error('Resume provenance/settings changed')
        work=Path(json.loads((out/'scratch_location.json').read_text())['directory']);worker=out/'backend/selected_response_worker'
        if not worker.is_file():p.error('Resume worker missing')
    else:
        if args.resume:p.error('Resume output does not exist')
        out.mkdir(parents=True);base.dump(out/'signature.json',signature);work=base.scratch_directory(args.scratch_root,out);base.dump(out/'scratch_location.json',dict(directory=str(work)));worker=compile_worker(out,reference,work);base.dump(out/'worker_build.json',dict(sha256=base.digest(worker)))
    if base.digest(worker)!=json.loads((out/'worker_build.json').read_text())['sha256']:p.error('Worker executable changed')
    spec['mass_variable']=dict(name='m_gamma_gamma',label='Two-photon mass (GeV)',bins=140,low=.10,high=.17);base.dump(out/'source_spec.json',spec);budget=base.Budget(args.workers)
    def export(stage,vectors):
        marker=out/(stage+'_exports.json')
        if marker.is_file():
            saved=json.loads(marker.read_text())
            if saved['vectors']!={p:list(v) for p,v in vectors.items()}:raise ValueError('Saved export response differs')
            records=saved['records'];result={}
            for r in records:
                if export_digest(r['path'])!=r['sha256']:raise ValueError('Saved export changed')
                result[r['period'],r['sample']]=arrays_from_export(Path(r['path']))
            print('[reuse] '+stage+' exports; no DATA reread',flush=True);return result
        jobs=[]
        for period in spec['periods']:
            if stage=='initial':
                for sample in ('data','pi0_data'):
                    entries=[r for r in ref['inventory'] if r['period']==period and r['sample']==sample]
                    if not entries:raise ValueError('Missing DATA inventory')
                    for r in entries:
                        st=Path(r['path']).stat()
                        if (st.st_size,st.st_mtime_ns)!=(r['size'],r['mtime_ns']):raise ValueError('DATA changed')
                    jobs.append(dict(period=period,sample=sample,files=[r['path'] for r in entries],export_individual_photons=sample=='pi0_data',export_global_only=False))
            for sample in ('mc','background','pi0_mc'):jobs.append(dict(period=period,sample=sample,cache=cache[period+'__'+sample]['path'],parameters=vectors[period],export_individual_photons=sample=='pi0_mc',export_global_only=False))
        pw=work/stage;pw.mkdir(exist_ok=True);pool=base.Pool(worker,out/'source_spec.json',pw,args.workers,budget);records=[];result={}
        # Unique scratch label allows a failed export phase to be retried safely.
        label='exports_'+str(len(list(pw.glob('exports_*'))))
        for _,job,path in pool.run('export-selected',jobs,label):result[job['period'],job['sample']]=arrays_from_export(path);records.append(dict(period=job['period'],sample=job['sample'],path=str(path),sha256=export_digest(path)))
        base.dump(marker,dict(vectors=vectors,records=records));return result
    print('[stage 1] Legacy cut development from nominal templates; DATA read once, all pre-cut MC retained',flush=True);initial=export('initial',{p:[0.]*20 for p in spec['periods']})
    def period_arrays(exports,period):return {s:exports[period,s] for s in ('data','mc','background','pi0_data','pi0_mc')}
    def refit_phase(stage,exports):
        results={}
        with ThreadPoolExecutor(max_workers=min(args.workers,len(spec['periods']))) as executor:
            fs={}
            for period in spec['periods']:
                folder=out/stage/period;marker=folder/'cuts_95.json'
                if marker.is_file():results[period]=json.loads(marker.read_text())
                else:fs[executor.submit(production_cut_fit,period,period_arrays(exports,period),folder)]=period
            for f in as_completed(fs):results[fs[f]]=f.result();print('[cuts] '+stage+'/'+fs[f],flush=True)
        return results
    initial_windows=refit_phase('initial_cuts',initial);base.dump(out/'initial_fixed_windows.json',initial_windows);base.table(out/'initial_selection_counts.csv',[r for period in spec['periods'] for r in selection_summary(period,period_arrays(initial,period),initial_windows[period],'initial')])
    fixed=copy.deepcopy(spec);fixed['dvcs_windows']={p:w['dvcs'] for p,w in initial_windows.items()};fixed['control_windows']={p:w['control'] for p,w in initial_windows.items()};base.dump(out/'fixed_fit_spec.json',fixed)
    print('[stage 2] Fixed DATA cuts; joint bias+smearing fit including unordered controls and pi0 mass; every trial reselects full MC',flush=True)
    fits={};engines={};finals={}
    def fit_period(period):
        folder=out/'response'/period;marker=folder/'fit_completed.json'
        if marker.is_file():
            saved=json.loads(marker.read_text())
            if any(base.digest(folder/name)!=value for name,value in saved.items()):raise ValueError('Completed fit checkpoint changed')
            return cuts.load_candidate_parameters(folder/'candidate_parameters.json'),None
        folder.mkdir(parents=True,exist_ok=True);one=copy.deepcopy(fixed);one['periods']=[period];base.dump(folder/'spec.json',one);pw=work/('fit_'+period);pw.mkdir(exist_ok=True);pool=base.Pool(worker,folder/'spec.json',pw,args.workers,budget)
        # Avoid overwriting scratch from a partially completed attempt.
        pool.output=pw/('attempt_'+str(len(list(pw.glob('attempt_*')))));pool.output.mkdir()
        data,counts=data_histograms(period,period_arrays(initial,period),initial_windows[period],spec);local={s:cache[period+'__'+s] for s in ('mc','background','pi0_mc')};uncertainty=float(meta['arguments'].get('background_transfer_uncertainty',.1));engine=Engine(pool,one,local,data,counts,args.min_mc,uncertainty);best,pars,scores,records,hist,selected=fit_response(engine,starts[period],folder,args.rounds,meta['arguments'].get('scan_range','wide'))
        # Keep only small final histograms for reproducible plots; trial arrays are deleted.
        payload={s:np.asarray(h) for s,h in hist.items()};np.savez_compressed(folder/'final_histograms.npz',**payload,**{'selected_'+s:v for s,v in selected.items()});base.dump(folder/'final_scores.json',dict(parameters=pars,scores=scores));base.discard_histograms(hist);base.dump(marker,{name:base.digest(folder/name) for name in ('candidate_parameters.json','final_scores.json','final_histograms.npz')})
        return best,engine
    with ThreadPoolExecutor(max_workers=min(args.workers,len(spec['periods']))) as executor:
        fs={executor.submit(fit_period,p):p for p in spec['periods']}
        try:
            for f in as_completed(fs):period=fs[f];fits[period],engines[period]=f.result();print('[fit] '+period+' complete',flush=True)
        except BaseException:
            budget.cancel_all()
            for f in fs:f.cancel()
            raise
    links=[]
    for period in spec['periods']:
        folder=out/'response'/period;data,counts=data_histograms(period,period_arrays(initial,period),initial_windows[period],spec);one=copy.deepcopy(fixed);one['periods']=[period];engine=Engine(None,one,{s:cache[period+'__'+s] for s in ('mc','background','pi0_mc')},data,counts,args.min_mc,0)
        with np.load(folder/'final_histograms.npz',allow_pickle=False) as z:hist={s:z[s] for s in ('mc','background','pi0_mc')};selected={s:z['selected_'+s] for s in hist};summary=json.loads((folder/'final_scores.json').read_text());links+=render_fit(period,engine,[tuple(p) for p in summary['parameters']],summary['scores'],hist,selected,starts[period],fits[period],out)
    print('[stage 3] Legacy cut redetermination with fitted templates; separate DATA/MC 95% windows',flush=True);candidate=export('candidate',fits)
    for period in spec['periods']:
        for sample in ('data','pi0_data'):candidate[period,sample]=initial[period,sample]
    final_windows=refit_phase('final_cuts',candidate);base.dump(out/'final_windows.json',final_windows);base.table(out/'final_selection_counts.csv',[r for period in spec['periods'] for r in selection_summary(period,period_arrays(candidate,period),final_windows[period],'final')]);rows=[];transfers=[];window_changes=[]
    for period in spec['periods']:
        rr,tt=final_comparison(period,period_arrays(initial,period),period_arrays(candidate,period),initial_windows[period],final_windows[period],{s:cache[period+'__'+s]['probability'] for s in ('mc','background','pi0_mc')},float(spec.get('control_to_leakage_exposure_ratio',1.)));rows+=rr;transfers+=tt
        for channel in ('dvcs','control'):
            for g,(old,new) in enumerate(zip(initial_windows[period][channel],final_windows[period][channel])):
                for kind in ('data','mc'):
                    for variable in sorted(set(old[kind])|set(new[kind])):window_changes.append(dict(period=period,channel=channel,group=g,sample=kind,variable=variable,old=json.dumps(old[kind].get(variable)),new=json.dumps(new[kind].get(variable))))
    for period in spec['periods']:
        a=period_arrays(candidate,period);windows=final_windows[period];data,counts=data_histograms(period,a,windows,spec)
        sh,sc=data_histograms(period,dict(data=a['mc'],pi0_data=a['pi0_mc']),windows,spec);bh,bc=data_histograms(period,dict(data=a['background'],pi0_data=a['pi0_mc']),windows,spec)
        hist={'mc':np.repeat(sh['epgamma'][None,None]/2,2,axis=1),'background':np.repeat(bh['epgamma'][None,None]/2,2,axis=1),'pi0_mc':np.repeat(sh['epgammagamma'][None,None]/2,2,axis=1)};selected={'mc':np.repeat(sc['epgamma'][None,None]/2,2,axis=1),'background':np.repeat(bc['epgamma'][None,None]/2,2,axis=1),'pi0_mc':np.repeat(sc['epgammagamma'][None,None]/2,2,axis=1)}
        one=copy.deepcopy(fixed);one['periods']=[period];engine=Engine(None,one,{s:cache[period+'__'+s] for s in ('mc','background','pi0_mc')},data,counts,args.min_mc,0);mult=json.loads((out/'response'/period/'candidate_parameters.json').read_text())['background_multiplier']
        dummy=[{'background_multiplier':mult}];ll=render_fit(period,engine,[tuple(fits[period])],dummy,hist,selected,fits[period],fits[period],out/'postfit_common_selection',final_selection=True);links+=['postfit_common_selection/'+x for x in ll]
    base.table(out/'relative_acceptance.csv',rows);base.table(out/'coarse_pi0_transfer.csv',transfers);base.table(out/'window_changes.csv',window_changes)
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Selection-first response</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Legacy cuts → fixed-selection response → redetermined cuts</h1><p>Initial legacy signal/background template cuts establish the DATA reference. DATA coordinates never change. All full pre-cut reconstructed MC remains available: response trials are applied before global/exclusivity selection, allowing migration in either direction. Fixed initial DATA windows define N-minus-one fit projections; windows are not reoptimized during a trial.</p><p>Periods/polarities are fitted independently. Epγ uses signal+leakage, with one period-wide control-constrained multiplier. Epγγ uses unordered photon groups, omitting FD-proton groups containing FT photons. Its global selection is the conservative intersection of both legacy detector assignments on the stored composite kinematics. Training uses Delta_phi, pTmiss, Emiss2 and Mx2, plus pi0 mass in controls; theta_gamma_gamma and Mx2_2 and the other MC half validate. DATA are shared. Support is fixed from DATA/nominal MC. Background transfer follows the newly selected control counts and saved exposure ratio.</p><p>The pi0 mass uses saved mass plus the change in rebuilt mass, preserving exact zero-response mass calibration; neither original photon coordinates nor stored calibration offsets are silently replaced. The upstream two-photon mass-window limitation remains. Mixed control windows use the CD-FD optimizer policy with their own distributions and fitted boundaries.</p><p>After fitting, the legacy template method redetermines separate DATA and MC 95% windows. Counts then use those actual windows. Relative acceptance is a paired reconstructed-count ratio with the unchanged generated denominator cancelling. Coarse pi0 transfer is |t|-binned only and uses raw counts/sampling and exposure corrections; it has no production event current/proton weights or xB/Q2/phi binning. These diagnostics are not final cross-section corrections and do not alter production outputs. Large window shifts indicate a possible later iteration; this run makes one explicit pass.</p><p><a href="relative_acceptance.csv">Relative acceptance</a> · <a href="coarse_pi0_transfer.csv">Coarse pi0 transfer</a> · <a href="window_changes.csv">Window changes</a> · <a href="initial_selection_counts.csv">Initial selection counts</a> · <a href="final_selection_counts.csv">Final selection counts</a> · <a href="initial_fixed_windows.json">Initial fixed windows</a> · <a href="final_windows.json">Final windows</a></p>'+''.join('<p><a href="'+x+'">'+html.escape(x)+'</a></p>' for x in links)+'<p>Cut-development plots are in initial_cuts/PERIOD and final_cuts/PERIOD. postfit_common_selection compares the fixed fitted response with DATA under the final DATA windows applied to both, omitting the displayed cut. Those overlays use newly selected control normalization with the fitted nuisance held fixed; the count-based coarse transfer table uses its unprofiled control normalization.</p>');base.dump(out/'completed.json',dict(signature=signature,report_sha256=base.digest(out/'REPORT.html')));print('[done] '+str(out/'REPORT.html')+'; large arrays stay on scratch',flush=True)

if __name__=='__main__':main()
