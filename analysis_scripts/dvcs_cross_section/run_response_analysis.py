#!/usr/bin/env python3
"""Run the complete fresh-fit response analysis in an isolated working directory."""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
SOURCE=Path(__file__).resolve().parent

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--response-run',type=Path,help='Explicit saved-fit override; omitted means fit inside this run');p.add_argument('--run-dir',required=True,type=Path);p.add_argument('--scratch-root',type=Path,default=Path('/scratch/thayward/dvcs_particle_response/production'));p.add_argument('--workers',type=int,default=10);p.add_argument('--skip-systematics',action='store_true');p.add_argument('--plot',action='store_true');p.add_argument('--theory-cache-from',type=Path,help='Copy only reusable final-by-t theory cache JSONs from a previous run');a=p.parse_args(argv)
    if not 1<=a.workers<=10:p.error('workers must be 1–10')
    model=a.response_run.expanduser().resolve() if a.response_run else None;run=a.run_dir.expanduser().resolve()
    if model is not None and not (model/'source_spec.json').is_file() and not (model/'spec.json').is_file():p.error('Missing fitted response run: '+str(model))
    for name in ('dvcs_analysis','main_systematics','production_response_writer'):
        if not (SOURCE/name).is_file():p.error('Build first: missing '+name)
    from prepare_production_particle_response import models
    if model is not None:models(model)  # Validate all five periods and response-law provenance before mkdir.
    if a.theory_cache_from and not a.theory_cache_from.expanduser().is_dir():p.error('Missing theory-cache source: '+str(a.theory_cache_from))
    run.mkdir(parents=True,exist_ok=False)
    for path in SOURCE.iterdir():
        if path.is_file() and (path.suffix in ('.py','.cpp','.h','.json','.sh') or path.name in ('dvcs_analysis','main_systematics','production_response_writer')):shutil.copy2(path,run/path.name)
    for name in ('external_scripts','external_pi0_normalization'):
        shutil.copytree(SOURCE/name,run/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc','output'))
    (run/'imports').symlink_to(SOURCE/'imports',target_is_directory=True)
    frozen=None
    if model is not None:
        # Snapshot the small fitted model; never depend on later edits to that run.
        frozen=run/'frozen_response';frozen.mkdir()
        for name in ('source_spec.json','spec.json','signature.json','manifest.json'):
            if (model/name).is_file():shutil.copy2(model/name,frozen/name)
        for period in ('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb'):
            source=model/'response'/period/'candidate_parameters.json'
            if not source.is_file():source=model/period/'candidate_parameters.json'
            target=frozen/'response'/period/'candidate_parameters.json';target.parent.mkdir(parents=True);shutil.copy2(source,target)
            marker=source.parent/'fit_completed.json'
            if marker.is_file():shutil.copy2(marker,target.parent/marker.name)
    if a.theory_cache_from:
        previous=a.theory_cache_from.expanduser().resolve()
        for source in previous.rglob('model_curves_cache_*GeV.json'):
            entries=json.loads(source.read_text())
            if not isinstance(entries,dict):raise ValueError('Invalid theory cache: '+str(source))
            target=run/'output/final_results_by_t'/source.parent.name/source.name;target.parent.mkdir(parents=True,exist_ok=True)
            if target.is_file():entries={**json.loads(target.read_text()),**entries}
            target.write_text(json.dumps(entries))
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['MPLCONFIGDIR']=str(run/'output/.matplotlib');env['OMP_NUM_THREADS']=str(a.workers)
    for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):env[name]='1'
    for name in ('DVCS_PARTICLE_RESPONSE_MANIFEST','DVCS_PARTICLE_RESPONSE_DATA_CACHE'):env.pop(name,None)
    command=[str(run/'dvcs_analysis'),*(['--particle-response-run',str(frozen)] if frozen is not None else ['--particle-response']),'--particle-response-scratch',str(a.scratch_root.expanduser().resolve()),'--particle-response-workers',str(a.workers)]
    if a.skip_systematics:command+=['--skip-systematics']
    (run/'run_manifest.json').write_text(json.dumps(dict(source=str(SOURCE),fitted_source=str(model) if model is not None else None,command=command),indent=2)+'\n');print('Private response run: '+str(run),flush=True)
    result=subprocess.run(command,cwd=run,env=env)
    if result.returncode:return result.returncode
    completed=run/'output/particle_response/extraction_manifest.json'
    if not completed.is_file() or not json.loads(completed.read_text()).get('reconstructed_mc_modified'):
        raise RuntimeError('The executable did not finalize a particle-response extraction; rebuild before running')
    if a.plot:return subprocess.run([sys.executable,str(run/'external_scripts/plot_pass2_final_results_by_t.py'),str(run/'output/csvs/dvcs_pass2_analysis.csv'),'--output-dir',str(run/'output/final_results_by_t'),'--workers',str(a.workers),'--reference-kinematics','means'],cwd=run,env=env).returncode
    return 0
if __name__=='__main__':sys.exit(main())
