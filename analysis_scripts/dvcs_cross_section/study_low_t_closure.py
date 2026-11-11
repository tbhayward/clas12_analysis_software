#!/usr/bin/env python3
"""Standalone production-cut paired-MC closure for nominal/merged low-|t| bins.

Both polarities share one worker pool. Native ROOT workers read each paired
file once, accumulating all three binning schemes, all spectrum stresses and angle
strata. Plotting/reports can be repeated from compact aggregates without ROOT.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import threading
import numpy as np

SOURCE=Path(__file__).resolve().parent
VERSION='paired-4d-low-t-closure-v1'
PERIODS={'fa18_inb':('Fa18 Inb','Fa18_Inb'),'fa18_out':('Fa18 Out','Fa18_Out')}
ANGLE_EDGES=[0,37,45,50,55,60,64.23,180.000001]
GROUPS=['inclusive']+[f'truth_theta_{a:g}_{b:g}' for a,b in zip(ANGLE_EDGES[:-1],ANGLE_EDGES[1:])]
SCENARIOS=[dict(name='nominal',slope=0.,bump=0.),
           dict(name='low_t_enhanced',slope=-3.,bump=0.),
           dict(name='low_t_suppressed',slope=3.,bump=0.),
           dict(name='steep_low_t_enhanced',slope=-6.,bump=0.),
           dict(name='first_bin_bump',slope=0.,bump=.5)]
FIELDS=['xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax']
T_LABELS=['<0','0-.05','.05-.08','.08-.11','.11-.15','.15-.25','.25-.40','.40-.60','.60-.80','.80-1.0','>=1.0']

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def bins_from_csv(path):
    with path.open(newline='') as handle:
        rows=[dict(id=r['bin index'],valid=int(float(r['valid bin'])),
                   bounds=[float(r[f]) for f in FIELDS]) for r in csv.DictReader(handle)]
    if not rows:raise ValueError('Analysis CSV has no bins')
    if not {(.11,.15),(.15,.25)}.issubset({tuple(r['bounds'][4:6]) for r in rows}):
        raise ValueError('The nominal first two t bins must be present for this study')
    if any(r['bounds'][4:6]==[.11,.25] for r in rows):
        raise ValueError('Use the nominal proton-t production CSV, not an already merged study')
    return rows

def merged_bins(rows):
    """Common phi coarsening of overlapping first/second-bin intervals.

    Different phi partitions cannot simply be concatenated: that double counts
    events. Connected overlapping intervals form a single merged bin, with the
    full .11-.25 t interval. Other t bins are untouched.
    """
    result=[dict(r) for r in rows if r['bounds'][4]>=.25-1e-9]
    groups={}
    for r in rows:
        if r['bounds'][4]<.25-1e-9:groups.setdefault(tuple(r['bounds'][:4]),[]).append(r)
    for cell,items in groups.items():
        components=[]
        for r in sorted(items,key=lambda r:r['bounds'][6]):
            low,high=r['bounds'][6:]
            if components and low<components[-1][1]-1e-9:
                components[-1][1]=max(components[-1][1],high)
                components[-1][2]=max(components[-1][2],r['valid'])
            else:components.append([low,high,r['valid']])
        for lo,hi,valid in components:
            result.append(dict(id=f'merged_{len(result)}',valid=valid,bounds=[*cell,.11,.25,lo,hi]))
    return sorted(result,key=lambda r:r['bounds'])

def coverage_control_bins(rows):
    # Separate t bins with exactly the same x/Q/phi coverage and total t domain
    # as the merged study, including physical areas absent from nominal rows.
    # This distinguishes t merging from phi/support changes.
    result=[]
    for r in merged_bins(rows):
        if r['bounds'][4]>=.25-1e-9:
            result.append(dict(r))
        else:
            for low,high in ((.11,.15),(.15,.25)):
                b=list(r['bounds']);b[4:6]=[low,high]
                result.append(dict(id=f'coverage_control_{len(result)}',valid=r['valid'],bounds=b))
    return sorted(result,key=lambda r:r['bounds'])


def make_grid(rows,name):
    edges=[sorted({r['bounds'][i] for r in rows}|{r['bounds'][i+1] for r in rows}) for i in (0,2,4,6)]
    shape=tuple(len(e)-1 for e in edges);lookup=np.full(shape,-1,dtype=np.int32)
    for ident,r in enumerate(rows):
        slices=tuple(slice(edges[a].index(r['bounds'][2*a]),edges[a].index(r['bounds'][2*a+1])) for a in range(4))
        if np.any(lookup[slices]>=0):raise ValueError(f'Overlapping analysis bins in {name}: {r}')
        lookup[slices]=ident
    return dict(name=name,bins=rows,edges=edges,lookup=lookup.ravel().tolist())

def build_spec(run,seed):
    paths=dict(global_cuts=run/'output/jsons/global_cuts_config.json',
               exclusivity=run/'output/jsons/combined_cuts.json',
               current_model=run/'output/dvcs_current_dependence/calibration/current_response_model.json',
               binning=run/'output/csvs/dvcs_pass2_analysis.csv',
               global_cpp=run/'global_cuts.cpp',global_header=run/'global_cuts.h',
               domain_header=run/'kinematic_bin_integration.h')
    for p in paths.values():
        if not p.is_file():raise FileNotFoundError(f'Missing production output: {p}')
    rows=bins_from_csv(paths['binning'])
    spec=dict(version=VERSION,seed=seed,angle_edges=ANGLE_EDGES,scenarios=SCENARIOS,
              period_tags={p:v[0] for p,v in PERIODS.items()},period_codes={p:v[1] for p,v in PERIODS.items()},
              schemes=[make_grid(rows,'nominal'),make_grid(coverage_control_bins(rows),'coverage_control'),make_grid(merged_bins(rows),'merged_low_t')],
              input_hashes={k:digest(p) for k,p in paths.items()})
    for key in ['global_cuts','exclusivity','current_model']:spec[key]=json.loads(paths[key].read_text())
    return spec,paths

def compile_worker(output,production_source):
    """Compile only the standalone worker and the actual production global cuts."""
    build=output/'backend';build.mkdir()
    flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
    libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
    worker=build/'low_t_closure_worker'
    source=build/'low_t_closure_worker.cpp';shutil.copy2(SOURCE/'low_t_closure_worker.cpp',source)
    cmd=[os.environ.get('CXX','c++'),*flags,'-O3','-std=c++17','-I'+str(production_source),
         str(source),str(production_source/'global_cuts.cpp'),'-o',str(worker),*libs]
    print('[build] Compiling standalone production-cut worker',flush=True)
    subprocess.run(cmd,check=True)
    return worker

def inventory(root,periods,limit):
    found={}
    for p in periods:
        directory=root/'mc'/p/'dvcs_mc/50nA/parts'
        files=sorted(f for f in directory.glob('*.root') if '.rec_only.root' not in f.name and '.partial' not in f.name)
        if limit:files=files[:limit]
        if not files:raise FileNotFoundError(f'No paired MC files in {directory}')
        found[p]=files
    return found

def make_jobs(files,workers):
    # A few balanced batches amortize ROOT startup/JSON output across parts.
    jobs=[]
    for p,paths in files.items():
        count=min(len(paths),max(1,2*workers//len(files)))
        batches=[[] for _ in range(count)];sizes=[0]*count
        for path in sorted(paths,key=lambda f:f.stat().st_size,reverse=True):
            k=int(np.argmin(sizes));batches[k].append(str(path));sizes[k]+=path.stat().st_size
        jobs.extend(dict(period=p,files=sorted(batch)) for batch in batches)
    return jobs

def empty(spec):
    nscenario=len(spec['scenarios']);ng=len(GROUPS)
    result=dict(migration=np.zeros((12,11)),topology=np.zeros((3,11,11)),angles=np.zeros((ng,11,11)))
    for si,s in enumerate(spec['schemes']):
        n=len(s['bins'])
        for key in ['train_gen','train_reco','train_reco2','train_diag']:result[f'{si}_{key}']=np.zeros((ng,n))
        for key in ['test_gen','test_gen2','test_reco','test_reco2','test_diag']:result[f'{si}_{key}']=np.zeros((nscenario,ng,n))
    return result

def accumulate(target,payload):
    for key in target:
        if key[0].isdigit():
            index,name=key.split('_',1);value=payload['schemes'][int(index)][name]
        else:value=payload[key]
        target[key]+=np.asarray(value,dtype=float).reshape(target[key].shape)

def calculate(arrays,index):
    g,r,r2,diag=(arrays[f'{index}_{k}'] for k in ['train_gen','train_reco','train_reco2','train_diag'])
    tg,tg2,tr,tr2,td=(arrays[f'{index}_{k}'] for k in ['test_gen','test_gen2','test_reco','test_reco2','test_diag'])
    with np.errstate(divide='ignore',invalid='ignore'):
        acceptance=r/g
        variance_a=r2/g**2+r**2/g**3-2*r*diag/g**3
        corrected=tr/acceptance[None,:,:]
        ratio=corrected/tg
        variance_relative=tr2/tr**2+tg2/tg**2-2*td/(tr*tg)+variance_a[None,:,:]/acceptance[None,:,:]**2
        error=np.abs(ratio)*np.sqrt(np.maximum(0,variance_relative))
    return acceptance,corrected,ratio,error

def usable_cells(a,index,scheme,min_train,min_test,min_selected):
    g=a[f'{index}_train_gen'];r=a[f'{index}_train_reco'];tg=a[f'{index}_test_gen'];tr=a[f'{index}_test_reco']
    with np.errstate(divide='ignore',invalid='ignore'):
        train_neff=r*r/a[f'{index}_train_reco2']
        test_neff=tr[0]*tr[0]/a[f'{index}_test_reco2'][0]
    valid=np.array([b['valid'] for b in scheme['bins']],dtype=bool)
    return (g>=min_train)&(tg[0]>=min_test)&(train_neff>=min_selected)&(test_neff>=min_selected)&valid[None,:]


def matched_low_t(spec,results,min_train,min_test,min_selected):
    control,merged=spec['schemes'][1:]
    lookup={tuple(b['bounds']):i for i,b in enumerate(control['bins'])}
    records=[]
    for period,a in results.items():
        _,cc,cr,_=calculate(a,1);_,mc,mr,_=calculate(a,2)
        cu=usable_cells(a,1,control,min_train,min_test,min_selected)
        mu=usable_cells(a,2,merged,min_train,min_test,min_selected)
        for bi,b in enumerate(merged['bins']):
            if abs(b['bounds'][4]-.11)>1e-9:continue
            bounds=list(b['bounds']);bounds[4:6]=[.11,.15];first=lookup[tuple(bounds)]
            bounds[4:6]=[.15,.25];second=lookup[tuple(bounds)]
            # Prove event support matches before comparing the corrections.
            for key in ['train_gen','train_reco','test_gen','test_reco']:
                left=a['1_'+key][...,first]+a['1_'+key][...,second]
                if not np.allclose(left,a['2_'+key][...,bi],rtol=1e-10,atol=1e-8):
                    raise ValueError('Control/merged event coverage mismatch: '+key)
            for sc,scenario in enumerate(spec['scenarios']):
                for gi,group in enumerate(GROUPS):
                    truth=a['2_test_gen'][sc,gi,bi]
                    total=(cc[sc,gi,first]+cc[sc,gi,second])/truth if truth>0 else np.nan
                    records.append(dict(period=period,scenario=scenario['name'],group=group,
                        **dict(zip(FIELDS,b['bounds'])),
                        control_first_ratio=cr[sc,gi,first],control_second_ratio=cr[sc,gi,second],
                        control_total_ratio=total,merged_ratio=mr[sc,gi,bi],
                        usable=int(cu[gi,first] and cu[gi,second] and mu[gi,bi])))
    return records


def write_reports(output,spec,results,flows,min_train,min_test,min_selected,dpi):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summaries=[]
    fields=['period','binning','scenario','group','bin_index',*FIELDS,'production_valid',
            'train_generated','train_reconstructed_weighted','test_generated','test_reconstructed_weighted',
            'acceptance','corrected','corrected_over_truth','ratio_stat_error','usable']
    with (output/'closure_by_cell.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader()
        for period,a in results.items():
            for si,scheme in enumerate(spec['schemes']):
                acceptance,corrected,ratio,error=calculate(a,si)
                n=len(scheme['bins']);g=a[f'{si}_train_gen'];r=a[f'{si}_train_reco'];tg=a[f'{si}_test_gen'];tr=a[f'{si}_test_reco']
                usable=usable_cells(a,si,scheme,min_train,min_test,min_selected)
                for sc,scenario in enumerate(spec['scenarios']):
                    for gi,group in enumerate(GROUPS):
                        for bi,b in enumerate(scheme['bins']):
                            record=dict(period=period,binning=scheme['name'],scenario=scenario['name'],group=group,
                                bin_index=b['id'],**dict(zip(FIELDS,b['bounds'])),production_valid=b['valid'],
                                train_generated=g[gi,bi],train_reconstructed_weighted=r[gi,bi],test_generated=tg[sc,gi,bi],
                                test_reconstructed_weighted=tr[sc,gi,bi],acceptance=acceptance[gi,bi],corrected=corrected[sc,gi,bi],
                                corrected_over_truth=ratio[sc,gi,bi],ratio_stat_error=error[sc,gi,bi],usable=int(usable[gi,bi]))
                            writer.writerow(record)
                        for bounds in sorted({tuple(b['bounds'][4:6]) for b in scheme['bins']}):
                            use=usable[gi]&np.array([tuple(b['bounds'][4:6])==bounds for b in scheme['bins']])&np.isfinite(ratio[sc,gi])
                            values=ratio[sc,gi,use]
                            summaries.append(dict(period=period,binning=scheme['name'],scenario=scenario['name'],group=group,
                                t_min=bounds[0],t_max=bounds[1],cells=int(use.sum()),
                                median=float(np.median(values)) if values.size else np.nan,
                                p16=float(np.quantile(values,.16)) if values.size else np.nan,
                                p84=float(np.quantile(values,.84)) if values.size else np.nan,
                                corrected_sum_over_truth=float(corrected[sc,gi,use].sum()/tg[sc,gi,use].sum()) if values.size else np.nan))
    with (output/'closure_summary.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    with (output/'cutflow.json').open('w') as handle:json.dump(flows,handle,indent=2)
    matched=matched_low_t(spec,results,min_train,min_test,min_selected)
    with (output/'matched_low_t_comparison.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(matched[0]));writer.writeheader();writer.writerows(matched)
    fig,axes=plt.subplots(1,len(results),figsize=(8*len(results),5),squeeze=False)
    for ax,period in zip(axes[0],results):
        for field,label in [('control_first_ratio','Separate first t bin'),('control_total_ratio','Sum of separately corrected low-t bins'),('merged_ratio','Merged low-t correction')]:
            values=[]
            for scenario in spec['scenarios']:
                cohort=[r[field] for r in matched if r['period']==period and r['group']=='inclusive' and r['scenario']==scenario['name'] and r['usable']]
                values.append(float(np.median(cohort)) if cohort else np.nan)
            ax.plot(range(len(values)),values,'o-',label=label)
        n=sum(r['period']==period and r['group']=='inclusive' and r['scenario']=='nominal' and r['usable'] for r in matched)
        ax.axhline(1,color='black',lw=1);ax.set(title=f'{PERIODS[period][0]} — {n} matched x/Q/phi cells',xticks=range(len(spec['scenarios'])),ylabel='Median recovered / known truth')
        ax.set_xticklabels([s['name'].replace('_',' ') for s in spec['scenarios']],rotation=35,ha='right');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(output/'matched_low_t_comparison.png',dpi=dpi);plt.close(fig)
    # All-projection diagnostics include outside-bin feed-in and an explicit miss row.
    fig,axes=plt.subplots(1,len(results),figsize=(8*len(results),6),squeeze=False)
    metrics=[]
    for ax,(p,a) in zip(axes[0],results.items()):
        matrix=a['migration'];truth=matrix.sum(axis=0);selected=matrix[:-1].sum(axis=0);rec=matrix[:-1].sum(axis=1);diagonal=np.diag(matrix[:-1])
        with (output/f'response_counts_{p}.csv').open('w',newline='') as handle:
            writer=csv.writer(handle);writer.writerow(['reco/truth',*T_LABELS])
            writer.writerows(zip(T_LABELS+['miss'],*matrix.T.tolist()))
        probability=np.divide(matrix,truth[None,:],out=np.zeros_like(matrix),where=truth[None,:]>0)
        image=ax.imshow(probability,origin='lower',aspect='auto',vmin=0,vmax=1)
        ax.set(xticks=range(11),yticks=range(12),title=PERIODS[p][0],xlabel='Generated |t|',ylabel='Reconstructed |t| / miss')
        ax.set_xticklabels(T_LABELS,rotation=60,ha='right');ax.set_yticklabels(T_LABELS+['miss'])
        for i,name in enumerate(T_LABELS):
            metrics.append(dict(period=p,t_bin=name,generated=truth[i],selected_truth=selected[i],
                stability=diagonal[i]/selected[i] if selected[i]>0 else np.nan,
                purity=diagonal[i]/rec[i] if rec[i]>0 else np.nan,
                reco_count=rec[i],feed_from_below_011=matrix[i,:4].sum()/rec[i] if rec[i]>0 else np.nan))
    fig.colorbar(image,ax=list(axes[0]),label='Probability per generated |t| bin',shrink=.7)
    fig.subplots_adjust(bottom=.23,wspace=.3,right=.87);fig.savefig(output/'migration_response.png',dpi=dpi,bbox_inches='tight');plt.close(fig)
    with (output/'migration_metrics.csv').open('w',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
    # Conditional topology metrics have no invented generated-topology denominator.
    with (output/'migration_by_topology_and_truth_angle.csv').open('w',newline='') as handle:
        w=csv.writer(handle);w.writerow(['period','stratum','truth_t_bin','selected_truth','same_t_bin','conditional_stability','purity'])
        for period,a in results.items():
            strata=[('reco_'+name,m) for name,m in zip(['FD_FD','CD_FD','CD_FT'],a['topology'])]
            strata += list(zip(GROUPS[1:],a['angles'][1:]))
            for name,m in strata:
                for i,label in enumerate(T_LABELS):
                    selected=m[:,i].sum();rec=m[i,:].sum();same=m[i,i]
                    w.writerow([period,name,label,selected,same,same/selected if selected else np.nan,same/rec if rec else np.nan])
    fig,axes=plt.subplots(len(results),len(spec['schemes']),figsize=(6.5*len(spec['schemes']),5*len(results)),squeeze=False)
    for row,period in enumerate(results):
        for col,scheme in enumerate(spec['schemes']):
            ax=axes[row,col]
            for scenario in spec['scenarios']:
                rs=[r for r in summaries if r['period']==period and r['binning']==scheme['name'] and r['scenario']==scenario['name'] and r['group']=='inclusive' and r['cells']]
                ax.plot([(r['t_min']+r['t_max'])/2 for r in rs],[r['median'] for r in rs],'o-',label=scenario['name'].replace('_',' '))
            ax.axhline(1,color='black',lw=1);ax.set(title=f'{PERIODS[period][0]} — {scheme["name"]}',xlabel='|t| bin center (GeV²)',ylabel='Median recovered / known truth across usable 4D cells');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(output/'closure_nominal_vs_merged.png',dpi=dpi);plt.close(fig)
    fig,axes=plt.subplots(1,len(results),figsize=(8*len(results),5),squeeze=False)
    for ax,period in zip(axes[0],results):
        for scheme in spec['schemes']:
            for scenario,style in [('nominal','o-'),('steep_low_t_enhanced','s--')]:
                values=[]
                for group in GROUPS[1:]:
                    row=next(r for r in summaries if r['period']==period and r['binning']==scheme['name'] and r['scenario']==scenario and r['group']==group and abs(r['t_min']-.11)<1e-9)
                    values.append(row['median'])
                ax.plot(range(len(values)),values,style,label=scheme['name']+' / '+scenario)
        ax.axhline(1,color='black',lw=1)
        ax.set(xticks=range(len(GROUPS)-1),title=PERIODS[period][0],ylabel='Median first-bin recovered / known truth',xlabel='Generated proton polar-angle range (degrees)')
        ax.set_xticklabels([f'{a:g}–{b:g}' if b<180 else f'>{a:g}' for a,b in zip(ANGLE_EDGES[:-1],ANGLE_EDGES[1:])],rotation=45,ha='right')
        ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.tight_layout();fig.savefig(output/'closure_by_truth_proton_angle.png',dpi=dpi);plt.close(fig)
    lines=['# Standalone low-|t| closure study','',
        'The worker uses the saved production global configuration and MC exclusivity windows. ',
        'Both halves use the production MC current weights; no DATA or Neupane weight is applied to MC. ',
        'The response-training and pseudo-data halves are disjoint and deterministic. ',
        'Spectrum weights vary continuously with generated |t|, including within and outside the published bins.','',
        'Closure is evaluated in the actual xB/Q²/|t|/phi bins. Reconstructed counts include feed-in from outside the truth binning. ',
        'Merged low-|t| bins use a common coarsening of overlapping phi intervals, avoiding double counting. ',
        'A matched-coverage control retains the two t bins over the same full x/Q/phi and t support as the merged study, isolating the t-merging effect. ',
        'Truth-proton-angle closure is conditional on known truth angle; it is a diagnostic, not a data correction. ',
        'Topology purity/stability uses reconstructed topology and is not an absolute topology efficiency.','',
        'The cell statistical error includes train/test ratio uncertainty and same-bin truth/reco covariance; correlations between cells and current-fit uncertainty are not included. ',
        'Closure ratios require no Frad/Fbin, model evaluation or luminosity because identical reporting factors cancel. ',
        'A failed closure demonstrates model/selection sensitivity; passing these stress cases does not prove all physical spectra close.','',
        '| Period | Binning | Scenario | First-bin usable cells | Median recovered/truth |','|---|---|---|---:|---:|']
    for r in summaries:
        if r['group']=='inclusive' and abs(r['t_min']-.11)<1e-9:
            lines.append(f'| {r["period"]} | {r["binning"]} | {r["scenario"]} | {r["cells"]} | {r["median"]:.4f} |')
    lines+=['','![Migration](migration_response.png)','','![Closure](closure_nominal_vs_merged.png)','',
            '![Matched low-t comparison](matched_low_t_comparison.png)','',
            '![Truth proton-angle closure](closure_by_truth_proton_angle.png)','',
            'All cell results, including cells below the precision threshold, are in `closure_by_cell.csv`. ',
            'Inspect `reconstructed_without_valid_truth` in `cutflow.json`: such entries are excluded from this paired-truth study. ',
            'No input analysis products are modified. Replotting from aggregates does not read ROOT files.']
    (output/'REPORT.md').write_text('\n'.join(lines)+'\n')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production-run',type=Path,required=True,help='Completed nominal run supplying cuts, MC response and actual binning')
    parser.add_argument('--input-root',type=Path,default=Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild'))
    parser.add_argument('--periods',nargs='+',choices=list(PERIODS),default=list(PERIODS))
    parser.add_argument('--output',type=Path,required=True,help='New directory; existing destinations are refused')
    parser.add_argument('--workers',type=int,default=10)
    parser.add_argument('--seed',type=int,default=20261111)
    parser.add_argument('--max-files-per-period',type=int,default=0,help='Development limit; zero means all files')
    parser.add_argument('--reuse-aggregates',type=Path,help='Previous study directory; skip ROOT reads and compilation')
    parser.add_argument('--min-train-generated',type=int,default=200)
    parser.add_argument('--min-test-generated',type=int,default=100)
    parser.add_argument('--min-selected',type=int,default=20,help='Minimum effective selected count in each independent half')
    parser.add_argument('--dpi',type=int,default=160)
    args=parser.parse_args(argv)
    if not 1<=args.workers<=10 or args.max_files_per_period<0 or min(args.min_train_generated,args.min_test_generated,args.min_selected)<1:
        parser.error('Workers must be 1–10; file limit nonnegative; precision thresholds positive')
    periods=list(dict.fromkeys(args.periods));run=args.production_run.expanduser().resolve();out=args.output.expanduser().absolute()
    if out.exists():parser.error(f'Output already exists: {out}')
    spec,paths=build_spec(run,args.seed)
    # Verify both requested periods have every saved topology cut and MC response before any event I/O.
    for p in periods:
        for topo in ['FD_FD','CD_FD','CD_FT']:
            block=spec['exclusivity'][f'DVCS_{PERIODS[p][1]}_{topo}']['mc']
            if not block:raise ValueError(f'Empty production MC exclusivity block: {p} {topo}')
        for region in ['FT',*['S'+str(i) for i in range(1,7)]]:
            value=spec['current_model']['mc']['ep->epg'][PERIODS[p][0]][region]['reference_factor']
            if not np.isfinite(value) or value<=0:raise ValueError(f'Invalid MC response: {p} {region}')
    out.mkdir(parents=True,exist_ok=False);os.environ.setdefault('MPLCONFIGDIR',str(out/'.matplotlib'))
    (out/'spec.json').write_text(json.dumps(spec,separators=(',',':')))
    snapshot=out/'production_inputs';snapshot.mkdir()
    for name,path in paths.items():shutil.copy2(path,snapshot/(name+path.suffix))
    source_hashes={n:digest(SOURCE/n) for n in ['study_low_t_closure.py','low_t_closure_worker.cpp']}
    source_hashes.update({n:digest(run/n) for n in ['global_cuts.cpp','global_cuts.h','kinematic_bin_integration.h']})
    event_hashes={k:v for k,v in source_hashes.items() if k!='study_low_t_closure.py'}
    spec_fingerprint=hashlib.sha256(json.dumps(spec,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    results={p:empty(spec) for p in periods};flows={p:{} for p in periods};files={};file_inventory={}
    if args.reuse_aggregates:
        source=args.reuse_aggregates.expanduser().resolve();old=json.loads((source/'manifest.json').read_text())
        if old['version']!=VERSION or old['input_hashes']!=spec['input_hashes'] or old['seed']!=args.seed or old['event_source_hashes']!=event_hashes or old['spec_fingerprint']!=spec_fingerprint:
            raise ValueError('Reusable aggregates do not match cuts/binning/current response, seed or worker sources')
        for p in periods:
            with np.load(source/f'aggregates_{p}.npz') as z:
                if set(z.files)!=set(results[p]):raise ValueError('Reusable aggregate schema mismatch')
                results[p]={k:z[k] for k in z.files}
            flows[p]=json.loads((source/'cutflow.json').read_text())[p]
        files={p:old['files'][p] for p in periods};file_inventory={p:old['file_inventory'][p] for p in periods};print('[reuse] No ROOT files will be read',flush=True)
    else:
        files=inventory(args.input_root,periods,args.max_files_per_period)
        file_inventory={p:[dict(path=str(f),size=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns) for f in paths] for p,paths in files.items()}
        worker=compile_worker(out,run);jobs=make_jobs(files,args.workers);work=out/'worker_results';work.mkdir()
        stopped=threading.Event();process_lock=threading.Lock();active=set()
        def stop_workers():
            stopped.set()
            with process_lock:
                for process in active:
                    if process.poll() is None:process.terminate()
        def execute(index,job):
            if stopped.is_set():raise RuntimeError('Cancelled after another worker failed')
            jobpath=work/f'job_{index}.json';destination=work/f'result_{index}.json';log=work/f'job_{index}.log'
            jobpath.write_text(json.dumps(job))
            with log.open('w') as stream:
                with process_lock:
                    if stopped.is_set():raise RuntimeError('Cancelled after another worker failed')
                    process=subprocess.Popen([str(worker),str(out/'spec.json'),str(jobpath),str(destination)],stdout=stream,stderr=subprocess.STDOUT)
                    active.add(process)
                status=process.wait()
                with process_lock:active.discard(process)
            if status:
                stop_workers()
                raise RuntimeError(f'Worker failed for {job["period"]}; see {log}')
            return job['period'],destination
        print(f'[study] {sum(map(len,files.values()))} paired parts; {len(jobs)} batches, {args.workers} total workers across {periods}',flush=True)
        with ThreadPoolExecutor(max_workers=min(args.workers,len(jobs))) as pool:
            pending={pool.submit(execute,i,job):i for i,job in enumerate(jobs)}
            ready={};next_index=0
            try:
                for count,future in enumerate(as_completed(pending),1):
                    p,destination=future.result();ready[pending[future]]=(p,destination)
                    # Reduce in job order for reproducible floating-point sums;
                    # keep only paths in memory while other workers finish.
                    while next_index in ready:
                        period,path=ready.pop(next_index);payload=json.loads(path.read_text())
                        accumulate(results[period],payload)
                        for key,value in payload['flow'].items():flows[period][key]=flows[period].get(key,0)+value
                        path.unlink();next_index+=1
                    print(f'[study] {count}/{len(jobs)} batches complete: {p}',flush=True)
            except BaseException:
                stop_workers()
                for future in pending:future.cancel()
                raise
    if not args.reuse_aggregates:
        for records in file_inventory.values():
            for record in records:
                st=Path(record['path']).stat()
                if st.st_size!=record['size'] or st.st_mtime_ns!=record['mtime_ns']:
                    raise RuntimeError('ROOT input changed during study: '+record['path'])
    for p,a in results.items():np.savez_compressed(out/f'aggregates_{p}.npz',**a)
    manifest=dict(version=VERSION,seed=args.seed,input_hashes=spec['input_hashes'],source_hashes=source_hashes,
                  event_source_hashes=event_hashes,spec_fingerprint=spec_fingerprint,file_inventory=file_inventory,
                  production_run=str(run),periods=periods,workers=args.workers,
                  files={p:[str(f) for f in paths] for p,paths in files.items()},
                  max_files_per_period=old['max_files_per_period'] if args.reuse_aggregates else args.max_files_per_period,reused_from=str(args.reuse_aggregates) if args.reuse_aggregates else None,
                  scenarios=SCENARIOS,angle_edges=ANGLE_EDGES,mc_current_weights=True,neupane_on_mc=False,
                  min_train_generated=args.min_train_generated,min_test_generated=args.min_test_generated,min_selected=args.min_selected)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    write_reports(out,spec,results,flows,args.min_train_generated,args.min_test_generated,args.min_selected,args.dpi)
    if any(digest(path)!=spec['input_hashes'][key] for key,path in paths.items()):raise RuntimeError('Production inputs changed during study')
    print(f'[done] {out}/REPORT.md',flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
