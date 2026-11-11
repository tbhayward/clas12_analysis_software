#!/usr/bin/env python3
"""Recalculate volumes/Fbin and compare frozen-normalization cross sections.
No ROOT files are read. The input extraction and previous tests are preserved.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
import argparse,csv,hashlib,json,math,shutil,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
os.environ["PYTHONDONTWRITEBYTECODE"]="1"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-csv',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--input-correction-convention',choices=('legacy-multiply','pass1-divide'),help='Required for old CSVs without a formula marker')
    parser.add_argument('--cuts-json',type=Path);parser.add_argument('--workers',type=int,default=10)
    parser.add_argument('--rtol',type=float,default=.002);parser.add_argument('--order',type=int,default=3)
    parser.add_argument('--on-unconverged',choices=('error','omit'),default='error',help='Optional preview: omit unresolved energy-rows instead of rejecting the combined CSV')
    parser.add_argument('--convention-only',action='store_true',help='Repair imported-factor division only; no KM15 integration or volume changes')
    parser.add_argument('--max-order',type=int,default=11);parser.add_argument('--plot',action='store_true')
    args=parser.parse_args();source=args.input_csv.expanduser().resolve();out=args.output_dir.expanduser().resolve()
    sys.path.insert(0,str(Path(__file__).resolve().parent/'external_scripts'))
    from kinematic_domain import Domain,find_config
    from km15_bin_center_batch import (recalculate,rescale_cross_sections,atomic_csv,write_diagnostics,triple,convert_input_convention)
    if not source.is_file():parser.error(f'Missing CSV: {source}')
    config=find_config(source,args.cuts_json);domain=Domain.from_json(config)
    original_bytes=source.read_bytes()
    import io
    rows=list(csv.DictReader(io.StringIO(original_bytes.decode('utf-8-sig'))))
    if not rows:parser.error('CSV has no data rows')
    # A numerical correction cannot repair overlapping extraction bins.
    groups={}
    for r in rows:
        k=tuple(float(r[c]) for c in ('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max'))
        groups.setdefault(k,[]).append((float(r['phimin']),float(r['phimax'])))
    for k,edges in groups.items():
        edges.sort()
        if any(a[1]>b[0] for a,b in zip(edges[:-1],edges[1:])):
            parser.error(f'Input has overlapping phi bins: {k}. Repair/recount that extraction first.')
    try:rows,input_convention=convert_input_convention(rows,args.input_correction_convention)
    except ValueError as e:parser.error(str(e))
    out.mkdir(parents=True,exist_ok=False)
    (out/'original.csv').write_bytes(original_bytes);shutil.copy2(config,out/'global_cuts_config.json')
    os.environ['MPLCONFIGDIR']=str(out/'.matplotlib')
    atomic_csv(out/'convention_fixed.csv',rows)
    if args.convention_only:
        (out/'numerical_corrections.json').write_text(json.dumps(dict(input=str(source),input_sha256=hashlib.sha256(original_bytes).hexdigest(),input_correction_convention=input_convention,output_correction_convention='pass1-divide',mode='convention-only',normalization='held fixed',systematics='not reevaluated'),indent=2)+'\n')
        if args.plot:
            script=Path(__file__).resolve().parent/'external_scripts/plot_pass2_final_results_by_t.py'
            subprocess.run([sys.executable,str(script),str(out/'convention_fixed.csv'),'--output-dir',str(out/'plots_convention_fixed'),'--workers',str(args.workers),'--no-pass1','--reference-kinematics','means'],check=True)
        if source.read_bytes()!=original_bytes:raise SystemExit('Input CSV changed during calculation')
        print(f'Saved convention-only repair in {out}',flush=True)
        return 0
    updated,diagnostics=recalculate(rows,domain,args.workers,args.order,args.max_order,args.rtol)
    write_diagnostics(out/'numerical_diagnostics.csv',diagnostics)
    if source.read_bytes()!=original_bytes:raise SystemExit('Input CSV changed during calculation; corrected CSVs not published')
    volume_only=[dict(r) for r in updated]
    for a,b in zip(rows,volume_only):
        for label in ('10.6 GeV','10.2 GeV'):b['Fbin, '+label]=a['Fbin, '+label]
    rescale_cross_sections(rows,volume_only);atomic_csv(out/'volumes_only.csv',volume_only)
    failures=[d for d in diagnostics if d['status'] not in ('volume','converged')]
    manifest=dict(input=str(source),input_sha256=hashlib.sha256(original_bytes).hexdigest(),
        version='allowed-domain-v3-pass1-divide',input_correction_convention=input_convention,output_correction_convention='pass1-divide',reference='measured_means',normalization='held fixed',
        systematics='not reevaluated',rtol=args.rtol,order=args.order,max_order=args.max_order,failed_rows=len(failures),on_unconverged=args.on_unconverged,omitted_energy_rows=[dict(row=d['row'],bin_index=rows[d['row']].get('bin index',''),label=d['label'],status=d['status']) for d in failures])
    (out/'numerical_corrections.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if failures and args.on_unconverged=='error':
        raise SystemExit(f'{len(failures)} energy-rows failed; accurate volumes_only.csv saved, combined CSV not published. See numerical_diagnostics.csv')
    for d in failures:updated[d['row']]['Fbin, '+d['label']]=''
    if failures:print(f'Preview omits {len(failures)} unresolved energy-rows; see numerical_diagnostics.csv',flush=True)
    rescale_cross_sections(rows,updated);atomic_csv(out/'new_fbin_and_volumes.csv',updated)
    comparison=[]
    for i,(a,b) in enumerate(zip(rows,updated)):
        for label in ('10.6 GeV','10.2 GeV'):
            means_label='10.6 GeV' if label=='10.6 GeV' else 'Sp19 Inb'
            if not all(math.isfinite(triple(a.get(prefix+', '+means_label,''))[0]) for prefix in ('xBavg','Q2avg','t_abs_avg','phiavg')):continue
            f0,f1=triple(a['Fbin, '+label])[0],triple(b['Fbin, '+label])[0]
            v0,v1=triple(a['bin_volume, '+label])[0],triple(b['bin_volume, '+label])[0]
            if not all(math.isfinite(v) and v>0 for v in (f0,v0,v1)):continue
            fratio=f1/f0 if math.isfinite(f1) and f1>0 else None
            d={c:a[c] for c in ('bin index','xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax')}
            d.update(row=i,energy=label,old_fbin=f0,new_fbin=f1 if fratio is not None else None,fbin_ratio=fratio,old_volume=v0,new_volume=v1,
                volume_only_xs_ratio=v0/v1,combined_xs_ratio=v0/v1/fratio if fratio is not None else None,status='converged' if fratio is not None else 'omitted')
            comparison.append(d)
    atomic_csv(out/'comparison.csv',comparison)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for j,e in enumerate(('10.6 GeV','10.2 GeV')):
        d=[z for z in comparison if z['energy']==e]
        for i,(col,title) in enumerate([('fbin_ratio','New/imported Fbin'),('volume_only_xs_ratio','Volume-only cross-section ratio')]):
            x=[z[col] for z in d if z[col] is not None];axs[i,j].hist(x,bins=60,histtype='step',lw=1.8);axs[i,j].axvline(1,c='0.4',ls='--');axs[i,j].set(title=e,xlabel=title,ylabel='CSV rows')
    fig.savefig(out/'numerical_comparison.png',dpi=160);plt.close(fig)
    print(f'Saved numerical comparison in {out}',flush=True)
    if args.plot:
        script=Path(__file__).resolve().parent/'external_scripts/plot_pass2_final_results_by_t.py'
        previous_plot=None
        for file,mode,dirname in [('convention_fixed.csv','means','plots_convention_fixed'),('volumes_only.csv','means','plots_volumes_only'),('new_fbin_and_volumes.csv','means','plots_new_fbin_and_volumes')]:
            if previous_plot is not None:
                for cache in previous_plot.rglob('model_curves_cache_*GeV.json'):
                    destination=out/dirname/cache.relative_to(previous_plot)
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copy2(cache,destination)
            cmd=[sys.executable,str(script),str(out/file),'--output-dir',str(out/dirname),
                 '--workers',str(args.workers),'--no-pass1','--reference-kinematics',mode]
            subprocess.run(cmd,check=True)
            previous_plot=out/dirname
    assert hashlib.sha256(source.read_bytes()).hexdigest()==manifest['input_sha256']
    return 0

if __name__=='__main__':sys.exit(main())
