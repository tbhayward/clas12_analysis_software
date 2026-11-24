#!/usr/bin/env python3
"""Nominal preliminary e' p pi- channel selection for dilution and asymmetries."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from pathlib import Path
import sys
from candidate_cache import build_cache,save_json
from selection import GROUPS,ChannelSelection,group_for_dataset
from workflow import fit_group,freeze
from export import export_dataset,audit_file

DEFAULT_OUTPUT=Path(__file__).resolve().parent/'output'

def cache_job(job):
    name,source,tree,output,options=job
    out=Path(output)/'calibration_cache'/name;out.mkdir(parents=True,exist_ok=True)
    path,meta=build_cache(source,tree,out,options)
    return name,(str(path),meta)

def parallel(function,jobs,workers):
    results=[]
    with ProcessPoolExecutor(max_workers=min(workers,len(jobs))) as pool:
        for future in as_completed([pool.submit(function,j) for j in jobs]):results.append(future.result())
    return results

def publish_index(output,reports):
    links=[]
    for group in GROUPS:
        links.append(f'<h2>{group}</h2><a href="fits/{group}/cuts.csv">Cuts</a> | <a href="fits/{group}/stability.csv">Stability</a>')
        for stage in ['mass','momentum','angle']:
            links.append('<p>'+stage+': '+' | '.join(f'<a href="fits/{group}/{stage}_xb{i}.png">xB {i}</a>' for i in range(1,5))+'</p>')
    links.append('<h2>Selected datasets</h2>')
    for r in sorted(reports,key=lambda r:r['dataset']):links.append(f'<p><a href="datasets/{r["dataset"]}.csv">{r["dataset"]}</a>: {r["selected_rows"]:,} combinations; unavailable bins: {r["unavailable_bins"]}</p>')
    (output/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Nominal channel selection</title><h1>Nominal preliminary channel selection</h1><a href="nominal_cuts.json">Versioned cut definitions</a>'+''.join(links))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--datasets',type=Path,default=Path(__file__).resolve().parents[1]/'imports/datasets.json')
    p.add_argument('--input-root',type=Path)
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--step-size',default='100 MB')
    p.add_argument('--max-entries',type=int,help='Partial calibration only; cannot publish nominal selected files')
    p.add_argument('--fit-only',action='store_true')
    p.add_argument('--apply-only',action='store_true')
    p.add_argument('--status',action='store_true')
    p.add_argument('--audit',action='store_true',help='Validate all selected files and cuts without refitting')
    p.add_argument('--rebuild-cache',action='store_true')
    args=p.parse_args()
    if args.workers<1 or (args.max_entries is not None and args.max_entries<1):p.error('workers and max-entries must be positive')
    if args.fit_only and args.apply_only:p.error('choose fit-only or apply-only')
    if args.max_entries and args.output==DEFAULT_OUTPUT:p.error('Use a separate --output directory for partial calibration')
    if args.max_entries and not args.fit_only:p.error('--max-entries is restricted to --fit-only, never nominal exports')
    if args.status:
        files=list(args.output.glob('calibration_cache/*/progress.json'))+list(args.output.glob('fits/*/progress.json'))+list(args.output.glob('datasets/*.progress.json'))
        if not files:print('No nominal run started')
        for path in sorted(files):print(str(path.relative_to(args.output))+': '+json.dumps(json.loads(path.read_text())))
        return
    if args.audit:
        import numpy as np
        import uproot
        registry=json.loads((args.output/'selected_inputs.json').read_text())
        selection=ChannelSelection(args.output/'nominal_cuts.json')
        if registry['cutset_id']!=selection.cutset_id:raise ValueError('Registry cut-set mismatch')
        config=json.loads(args.datasets.read_text())
        if set(registry['files'])!=set(config['files']):raise ValueError('Selected dataset coverage mismatch')
        total=0
        for name,item in registry['files'].items():
            report=json.loads((args.output/'datasets'/(name+'.json')).read_text())
            if item['selected_rows']!=report['selected_rows']:raise ValueError('Registry/report count mismatch')
            if report['partial'] or report['cutset_id']!=selection.cutset_id:raise ValueError('Partial or stale dataset '+name)
            if item['unavailable_bins']!=(np.flatnonzero(~selection.active_bins(name))+1).tolist():raise ValueError('Availability mismatch '+name)
            selected_path=Path(item['path']);selected_path=selected_path if selected_path.is_absolute() else args.output/selected_path
            audit_file(selected_path,registry['tree'],report['selected_rows'],{k:np.dtype(v) for k,v in report['original_types'].items()},selection.cutset_id)
            with uproot.open(selected_path) as f:
                for a in f[registry['tree']].iterate(['analysis_bin','exclusivity_group','Mx2_deuteron','p_miss','exclusivity_reference_assisted'],library='np',step_size='50 MB'):
                    if not np.all(a['exclusivity_group']==selection.group_ids[selection.group(name)]):raise ValueError('Group label mismatch')
                    if not np.array_equal(a['exclusivity_reference_assisted'].astype(bool),selection.assisted_bins(name)[a['analysis_bin']-1]):raise ValueError('Reference provenance mismatch')
                    if not np.all(selection.candidate_mask(a['Mx2_deuteron'],a['p_miss'],a['analysis_bin']-1,name)):raise ValueError('Exported rows fail cut set')
            total+=report['selected_rows'];print(f'AUDITED {name}: {report["selected_rows"]:,} combinations')
        print(f'Audit passed: {len(registry["files"])} datasets, {total:,} selected combinations')
        return
    args.output.mkdir(parents=True,exist_ok=True)
    print('Nominal output: '+str(args.output.resolve()),flush=True)
    config=json.loads(args.datasets.read_text());root=args.input_root or Path(config['campaign'])/'root'
    # Validate the full downstream registry before any expensive work.
    for name in config['files']:group_for_dataset(name)
    cutfile=args.output/'nominal_cuts.json'
    if not args.apply_only:
        names=[n for members in GROUPS.values() for n in members]
        jobs=[]
        for name in names:
            source=root/config['files'][name]
            if not source.is_file():p.error('Missing calibration input '+str(source))
            jobs.append((name,str(source),config['tree'],str(args.output),dict(max_entries=args.max_entries,step_size=args.step_size,rebuild_cache=args.rebuild_cache)))
        caches=dict(parallel(cache_job,jobs,args.workers))
        results=parallel(fit_group,[(g,n,caches,str(args.output)) for g,n in GROUPS.items() if g!='RGC_summer2022'],args.workers)
        reference=next(r for r in results if r['group']=='RGC_fall2022_spring2023')
        results.append(fit_group(('RGC_summer2022',GROUPS['RGC_summer2022'],caches,str(args.output),reference)))
        cutfile,payload=freeze(results,caches,args.output)
        print('Cut-set ID: '+payload['cutset_id'],flush=True)
        for group,data in payload['groups'].items():
            missing=[b['bin'] for b in data['bins'] if not b['active']]
            print(group+': unavailable bins '+str(missing),flush=True)
    reports=[]
    if not args.fit_only:
        selection=ChannelSelection(cutfile)
        jobs=[]
        for name,filename in config['files'].items():
            source=root/filename
            if not source.is_file():p.error('Missing application input '+str(source))
            jobs.append((name,str(source),config['tree'],str(cutfile),str(args.output),args.step_size,None))
        reports=parallel(export_dataset,jobs,args.workers)
        save_json(args.output/'selected_inputs.json',dict(version=selection.payload['version'],cutset_id=selection.cutset_id,cut_file='nominal_cuts.json',tree=config['tree'],branches=dict(bin='analysis_bin',azimuth='phi2',helicity='helicity',beam_polarization='beam_pol',target_polarization='target_pol',reference_assisted='exclusivity_reference_assisted'),files={r['dataset']:dict(path='selected/'+Path(r['selected_file']).name,fit_group=r['fit_group'],selected_rows=r['selected_rows'],unavailable_bins=r['unavailable_bins'],reference_assisted_bins=r['reference_assisted_bins']) for r in reports},normalization='Use original run-based charge and polarization normalization. Do not infer luminosity from selected candidate counts.'))
        print(f'Completed {len(reports)} selected ROOT files',flush=True)
    publish_index(args.output,reports)
    print('Browse '+str(args.output/'index.html'),flush=True)

if __name__=='__main__':main()
