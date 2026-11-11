#!/usr/bin/env python3
"""Fresh proton-t extraction, optionally resetting the three named test runs.
Preserves final-by-t model caches before deletion; never rescales an old CSV.
"""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
SOURCE=Path(__file__).resolve().parent
RUN_NAMES=('merged_t_001','photon_t_001','proton_t_default_001')

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test-root',type=Path,default=SOURCE/'test_runs')
    parser.add_argument('--reset-test-runs',action='store_true',help='DELETE the three named run directories, preserving model caches first')
    parser.add_argument('--plot',action='store_true')
    args=parser.parse_args(argv)
    root=args.test_root.expanduser().absolute()
    if root.is_symlink():parser.error('Test root must not be a symlink')
    for name in ('dvcs_analysis','main_systematics'):
        if not (SOURCE/name).is_file():parser.error(f'Build first: missing {name}')
    root.mkdir(parents=True,exist_ok=True)
    lock=root/'.nominal_baseline_reset_lock'
    lock.mkdir()  # Prevent two reset launchers sharing these destinations.
    try:
        targets=[root/n for n in RUN_NAMES]
        for target in targets:
            if target.is_symlink() or (target.exists() and not target.is_dir()):
                parser.error(f'Refusing non-directory/symlink target: {target}')
        nominal=root/RUN_NAMES[-1]
        if nominal.exists() and not args.reset_test_runs:
            parser.error('Nominal destination exists; choose a different --test-root or explicitly reset test runs')
        cache_root=root/'theory_cache_preserved'
        if cache_root.is_symlink():parser.error('Cache destination must not be a symlink')
        caches={}
        files=list(cache_root.glob('model_curves_cache_*GeV.json')) if cache_root.exists() else []
        for target in targets:files.extend(sorted(target.rglob('model_curves_cache_*GeV.json')))
        for file in files:
            entries=json.loads(file.read_text())
            if not isinstance(entries,dict):raise ValueError(f'Invalid model cache: {file}')
            caches.setdefault(file.name,{}).update(entries)
        cache_root.mkdir(exist_ok=True)
        for name,entries in caches.items():
            temporary=cache_root/(name+'.tmp')
            temporary.write_text(json.dumps(entries));temporary.replace(cache_root/name)
        print(f'Preserved {sum(map(len,caches.values()))} model entries in {cache_root}',flush=True)
        if args.reset_test_runs:
            for target in targets:
                if target.exists():
                    print(f'Deleting test directory: {target}',flush=True)
                    shutil.rmtree(target)
        nominal.mkdir(exist_ok=False)
        for path in SOURCE.iterdir():
            if path.is_file() and (path.suffix in ('.py','.cpp','.h','.json','.sh') or path.name in ('dvcs_analysis','main_systematics')):
                shutil.copy2(path,nominal/path.name)
        for name in ('external_scripts','external_pi0_normalization'):
            shutil.copytree(SOURCE/name,nominal/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','output'))
        (nominal/'imports').symlink_to(SOURCE/'imports',target_is_directory=True)
        for name,entries in caches.items():
            energy=name.removeprefix('model_curves_cache_').removesuffix('GeV.json')
            dest=nominal/'output/final_results_by_t'/f'{float(energy):.1f}_GeV'/name
            dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(entries))
        (nominal/'.dvcs_nominal_run.json').write_text(json.dumps(dict(source=str(SOURCE),mode='nominal proton t',systematics=False),indent=2)+'\n')
        env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
        env['MPLCONFIGDIR']=str(nominal/'output/.matplotlib');env['OMP_NUM_THREADS']='10'
        for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[name]='1'
        result=subprocess.run([str(nominal/'dvcs_analysis'),'--skip-systematics'],cwd=nominal,env=env)
        if result.returncode:return result.returncode
        if args.plot:
            return subprocess.run([sys.executable,str(nominal/'external_scripts/plot_pass2_final_results_by_t.py'),str(nominal/'output/csvs/dvcs_pass2_analysis.csv'),'--output-dir',str(nominal/'output/final_results_by_t'),'--workers','10','--no-pass1'],cwd=nominal,env=env).returncode
        return 0
    finally:
        lock.rmdir()

if __name__=='__main__':sys.exit(main())
