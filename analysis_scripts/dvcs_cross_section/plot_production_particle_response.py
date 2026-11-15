#!/usr/bin/env python3
"""Publication figures from cached DATA arrays and bounded paired MC note samples."""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,json,html
from pathlib import Path
import numpy as np
import study_particle_response as base
import plot_exclusivity_data_dvcs_pi0_mc as production
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CORE=(0,3,4,5);RESERVED=(6,);COLORS=dict(DATA='#205cab',nominal='#888888',modified='#22946b')
LABELS={'Delta_phi':r'$\Delta\phi$ (rad)','pTmiss':r'$p_T^{\mathrm{miss}}$ (GeV)','Emiss2':r'$E_{\mathrm{miss}}$ (GeV)','Mx2':r'$M_X^2(ep\gamma)$ (GeV$^2$)','Mx2_2':r'$M_X^2(e\gamma)$ (GeV$^2$)','m_gamma_gamma':r'$M_{\gamma\gamma}$ (GeV)'}

def read(p):return json.loads(Path(p).read_text())
def mask(values,window,omit=None):
    n=len(next(iter(values.values())));keep=np.ones(n,dtype=bool)
    for name,w in window.items():
        if name==omit:continue
        v=values[name];keep&=np.isfinite(v)&(v<=float(w['cut_high']))
        if w.get('mode')!='upper_quantile':keep&=v>=float(w['cut_low'])
    return keep

def histogram(values,edges):
    values=values[np.isfinite(values)];return np.bincount(np.searchsorted(edges,values,side='right'),minlength=len(edges)+1)

def data_arrays(period):
    data=production.load_selected_arrays_for_file(period.data_file,period.key,production.TOPOLOGIES,'80 MB',production.DVCS_SHAPE_VARIABLES)
    ctrl=production.load_selected_arrays_for_file(period.eppi0_data_file,period.key,production.TOPOLOGIES,'80 MB',production.PI0_SHAPE_VARIABLES+(production.VariableConfig('Mh_gammagamma',LABELS['m_gamma_gamma'],140,.10,.17),))
    def convert(a,control):
        return [{name:a[t.key]['theta_pi0_pi0' if control and name=='theta_gamma_gamma' else name] for name in base.OBS}|({'m_gamma_gamma':a[t.key]['Mh_gammagamma']} if control else {}) for t in production.TOPOLOGIES]
    return {'epgamma':convert(data,False),'epgammagamma':convert(ctrl,True)}

def note_histograms(routes,period,sample,blocks,variables):
    nv=len(variables);nb=max(v['bins']+2 for v in variables);hist=np.zeros((22,2,3,nv,nb),dtype=np.int64);selected=np.zeros((22,2,3),dtype=np.int64);probabilities=set();unresolved=0
    for r in routes.values():
        if r['period']!=period or not any(x['sample']==sample for x in r['roles']):continue
        meta=r['summary'];probabilities.add(float(meta['note_probability']));shape=meta['note_shape'];a=np.memmap(r['notes'],dtype='<f8',mode='r',shape=tuple(shape)) if shape[0] else np.empty(shape)
        unresolved+=int(meta['unresolved'])
        if shape[2]!=31:raise ValueError('Response note sample predates selection-independent kinematics')
        for scenario in range(22):
            z=a[:,scenario]
            # Reapply the current selection to the saved frame, independently of
            # the global cuts in force when this reusable response file was built.
            physical_names=('x','Q2','t1','detector1','detector2','e_p','e_theta','e_phi','p1_p','p1_theta','p1_phi','p2_p','p2_theta','p2_phi','open_angle_ep2','pTmiss')
            frame={name:z[:,14+j] for j,name in enumerate(physical_names)};frame['runnum']=z[:,30]
            current_global=production.global_cut_mask_uproot(frame,period,frame['detector1'],frame['detector2'])
            for g in range(3):
                values={name:z[:,j+1] for j,name in enumerate(base.OBS)};values['m_gamma_gamma']=z[:,12];population=(z[:,0]==g)&current_global;full=population&mask(values,blocks[g]['data'])
                for half in (0,1):
                    hh=z[:,10]==half;selected[scenario,half,g]+=int((full&hh).sum())
                    for j,v in enumerate(variables):
                        keep=population&hh&mask(values,blocks[g]['data'],v['name']);edges=base.diag.variable_edges(v);hist[scenario,half,g,j,:len(edges)+1]+=histogram(values[v['name']][keep],edges)
    if len(probabilities)!=1 or next(iter(probabilities))<=0:raise ValueError('Missing/inconsistent diagnostic sampling: '+period+'/'+sample)
    return hist,selected,next(iter(probabilities)),unresolved

def predictions(hist,selected,prob,counts,scenario,half,exposure=1.):
    nd=counts['epgammagamma'].sum();ctrl=selected['pi0_mc'][scenario,half].sum()*2/prob['pi0_mc'];alpha=exposure*nd/ctrl if ctrl>0 else np.nan
    signal=selected['mc'][scenario,half]*2/prob['mc'];background=selected['background'][scenario,half]*2/prob['background'];beta=np.divide(counts['epgamma']-alpha*background,signal,out=np.full(3,np.nan),where=signal>0)
    beta=np.where(beta>0,beta,np.nan)
    return alpha,beta

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--response-output',type=Path);p.add_argument('--cuts',type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--replot',action='store_true',help='Render only from saved note histograms/tables, with no DATA or MC cache reads');a=p.parse_args(argv)
    from production_response_note_render import render
    out=a.output.resolve()
    if a.replot:render(out);return
    if a.response_output is None or a.cuts is None:p.error('--response-output and --cuts are required unless --replot is used')
    origin=a.response_output.resolve();manifest=read(origin/'input_manifest.json');cuts=read(a.cuts);out.mkdir(parents=True,exist_ok=True)
    os.environ['DVCS_PARTICLE_RESPONSE_MANIFEST']=str(origin/'input_manifest.json')
    if manifest.get('data_array_cache'):os.environ['DVCS_PARTICLE_RESPONSE_DATA_CACHE']=manifest['data_array_cache']
    production.ACTIVE_GLOBAL_CUTS=production.load_global_cuts_config(manifest['global_cuts_identity']['path'])
    if manifest['global_cuts']!=read(manifest['global_cuts_identity']['path']):raise ValueError('Global cuts changed since response preparation')
    variables=[dict(name=n,label=LABELS.get(n,n),bins=nb,low=lo,high=hi) for n,nb,lo,hi in [('Delta_phi',240,np.pi-.6,np.pi+.6),('theta',120,0,np.pi),('theta_gamma_gamma',120,0,5),('pTmiss',240,0,.6),('Emiss2',240,-1.5,2.5),('Mx2',400,-.04,.04),('Mx2_2',240,-1.5,4.5)]];mass=dict(name='m_gamma_gamma',label=LABELS['m_gamma_gamma'],bins=140,low=.10,high=.17)
    scale_rows=[];score_rows=[];param_rows=[];saved_hist={'edges_mass':base.diag.variable_edges(mass)};normalizations=[]
    for j in CORE:saved_hist[f'edges_{j}']=base.diag.variable_edges(variables[j])
    for period in production.PERIODS:
        if period.key not in manifest['model']['periods']:continue
        key=period.key;exposure=float(manifest['model'].get('control_to_leakage_exposure_ratio',1.));data=data_arrays(period);blocks={ch:[cuts[production.cut_json_key('dvcs' if ch=='epgamma' else 'eppi0',period,t)] for t in production.TOPOLOGIES] for ch in data};counts={ch:np.array([int(mask(v,b['data']).sum()) for v,b in zip(data[ch],blocks[ch])]) for ch in data};dh={}
        for ch in data:
            vv=variables+([mass] if ch=='epgammagamma' else []);dh[ch]=np.zeros((3,len(vv),max(v['bins']+2 for v in vv)),dtype=np.int64)
            for g in range(3):
                for j,v in enumerate(vv):
                    edges=base.diag.variable_edges(v);keep=mask(data[ch][g],blocks[ch][g]['data'],v['name']);dh[ch][g,j,:len(edges)+1]=histogram(data[ch][g][v['name']][keep],edges)
        hist={};selected={};prob={}
        for sample in ('mc','background','pi0_mc'):
            ch='epgammagamma' if sample=='pi0_mc' else 'epgamma';hist[sample],selected[sample],prob[sample],unresolved=note_histograms(manifest['routes'],key,sample,blocks[ch],variables+([mass] if ch=='epgammagamma' else []));normalizations.append(dict(period=key,sample=sample,probability=prob[sample],unresolved_full_sample=unresolved))
        parameters=np.array(manifest['model']['periods'][key]['parameters']).reshape(5,4)
        for r,region in enumerate(base.REGIONS):
            for j,coordinate in enumerate(('p_sigma_fraction','theta_sigma_rad','phi_sigma_rad','p_bias_fraction')):param_rows.append(dict(period=key,particle=region,coordinate=coordinate,value=parameters[r,j]))
        records=[]
        for scenario in range(22):
          for half in (0,1):
            alpha,beta=predictions(hist,selected,prob,counts,scenario,half,exposure);scale_rows.append(dict(period=key,scenario=scenario,half=half,alpha=alpha,beta_FD_FD=beta[0],beta_CD_FD=beta[1],beta_CD_FT=beta[2]));refalpha,refbeta=predictions(hist,selected,prob,counts,0,half,exposure)
            for ch in data:
              vv=variables+([mass] if ch=='epgammagamma' else [])
              for g in range(3):
                for j in (*CORE,*RESERVED,*([7] if ch=='epgammagamma' else [])):
                    n=vv[j]['bins']+2;d=dh[ch][g,j,:n]
                    if ch=='epgamma':
                        weights=[2*beta[g]/prob['mc'],2*alpha/prob['background']];rw=[2*refbeta[g]/prob['mc'],2*refalpha/prob['background']];m=sum(w*hist[s][scenario,half,g,j,:n] for w,s in zip(weights,('mc','background')));var=sum(w*w*hist[s][scenario,half,g,j,:n] for w,s in zip(weights,('mc','background')));ref=sum(w*hist[s][0,half,g,j,:n] for w,s in zip(rw,('mc','background')));refvar=sum(w*w*hist[s][0,half,g,j,:n] for w,s in zip(rw,('mc','background')));raw=sum(hist[s][0,half,g,j,:n].sum() for s in ('mc','background'))
                    else:m=hist['pi0_mc'][scenario,half,g,j,:n].astype(float);var=m.copy();ref=hist['pi0_mc'][0,half,g,j,:n].astype(float);refvar=ref.copy();raw=ref.sum()
                    score,bins=base.shape_score(d,m,var,int(raw),100,ref,refvar);records.append(dict(period=key,channel=ch,topology=base.TOPS[g],variable=vv[j]['name'],scenario='nominal' if scenario==0 else 'modified' if scenario==1 else 'without_'+base.PARAMS[scenario-2],half=half,score=score,usable_bins=bins,raw_nominal_mc=int(raw)))
        support={(r['channel'],r['topology'],r['variable'],r['half']) for r in records if r['scenario']=='nominal' and np.isfinite(r['score'])};score_rows.extend(dict(r,status='supported' if (r['channel'],r['topology'],r['variable'],r['half']) in support and np.isfinite(r['score']) else 'insufficient_or_invalid') for r in records)
        for ch in data:
          for g in range(3):
            for j in CORE:
                n=variables[j]['bins']+2;prefix=f'{key}_{ch}_{g}_{j}_';saved_hist[prefix+'DATA']=dh[ch][g,j,:n]
                for scenario,name in ((0,'nominal'),(1,'modified')):
                    if ch=='epgamma':
                        signal=np.zeros(n);background=np.zeros(n)
                        for half in (0,1):
                            alpha,beta=predictions(hist,selected,prob,counts,scenario,half,exposure);signal+=beta[g]*hist['mc'][scenario,half,g,j,:n]/prob['mc'];background+=alpha*hist['background'][scenario,half,g,j,:n]/prob['background']
                        m=signal+background
                    else:m=hist['pi0_mc'][scenario,:,g,j,:n].sum(axis=0).astype(float);signal=m;background=np.zeros(n)
                    saved_hist[prefix+name]=m;saved_hist[prefix+name+'_signal']=signal;saved_hist[prefix+name+'_background']=background
        for g in range(3):
            n=mass['bins']+2;saved_hist[f'{key}_mass_{g}_DATA']=dh['epgammagamma'][g,7,:n]
            for scenario,name in ((0,'nominal'),(1,'modified')):saved_hist[f'{key}_mass_{g}_{name}']=hist['pi0_mc'][scenario,:,g,7,:n].sum(axis=0)
        print('[response note] '+key+' histograms complete; cached arrays only',flush=True)
    base.table(out/'component_scales.csv',scale_rows);base.table(out/'parameters.csv',param_rows);base.table(out/'shape_scores.csv',score_rows);base.table(out/'sample_normalizations.csv',normalizations);np.savez_compressed(out/'plot_histograms.npz',**saved_hist)
    base.dump(out/'provenance.json',dict(model_id=manifest['model_id'],manifest_sha256=base.digest(origin/'input_manifest.json'),cuts_sha256=base.digest(a.cuts),script_sha256=base.digest(Path(__file__))))
    render(out)

if __name__=='__main__':main()
