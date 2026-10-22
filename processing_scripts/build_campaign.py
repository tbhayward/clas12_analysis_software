#!/usr/bin/env python3
"""DVCS processing campaign. Dry run by default; --execute runs at most 8 directory-isolated workers."""
import argparse
import csv
import json
import os
import shlex
import subprocess
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPTS = {'dvcs_data':'processing_dvcs.groovy','pi0_data':'processing_exclusive_pi0.groovy','dvcs_gen':'processing_mc_dvcs.groovy','dvcs_rec':'processing_mc_rec_dvcs.groovy','pi0_gen':'processing_mc_exclusive_pi0.groovy','pi0_rec':'processing_mc_rec_exclusive_pi0.groovy','pi0_as_dvcs':'processing_mc_rec_dvcs.groovy'}
PROGRESS_LOCK = threading.Lock()
PROGRESS = {'done': 0, 'total': 0, 'started': time.monotonic()}

def report_progress(label):
    with PROGRESS_LOCK:
        PROGRESS['done'] += 1
        done, total = PROGRESS['done'], PROGRESS['total']
        elapsed = (time.monotonic() - PROGRESS['started']) / 60
        print(f'[PROGRESS] {done}/{total} jobs ({100*done/total:.1f}%) | elapsed={elapsed:.1f} min | {label}', flush=True)

UNVERIFIED = {'pi0_as_dvcs':'Verify aaogen-as-DVCS selection before production.'}

def prepare_shared():
    print('Preparing shared checkout and ROOT converters once before workers start.', flush=True)
    subprocess.run(['git','pull'], check=True)
    jar = Path('processing_classes/dist/processing_classes.jar')
    service = Path('coatjava/lib/services/processing_classes.jar')
    if not jar.is_file(): raise RuntimeError(f'Missing locally built Java JAR: {jar}')
    import shutil
    service.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(jar, service)
    print('Installed processing_classes.jar once before starting workers.', flush=True)
    for name in ('convert_txt_to_root', 'convert_dvcs_mc_diagnostics'):
        source = Path('processing_scripts')/(name+'.cpp')
        if not source.is_file(): raise RuntimeError(f'Missing converter source: {source}')
        cflags = shlex.split(subprocess.check_output(['root-config','--cflags'], text=True).strip())
        libs = shlex.split(subprocess.check_output(['root-config','--libs'], text=True).strip())
        subprocess.run(['g++', *cflags, str(source), '-o', str(source.with_suffix('')), *libs], check=True)

def run_group(group, pilot, skip_existing):
    env = os.environ.copy()
    env['DVCS_CAMPAIGN_WORKER'] = '1'
    for r in group:
        expected = [Path(x) for x in r['expected'].split(';')]
        if all(x.is_file() and x.stat().st_size > 0 for x in expected) and skip_existing:
            print('SKIP', r['output'], flush=True)
            report_progress('skipped '+r['period']+' '+r['mode'])
            continue
        if any(x.exists() for x in expected):
            raise RuntimeError('Existing or partial outputs, refusing overwrite: '+r['output'])
        Path(r['output']).parent.mkdir(parents=True, exist_ok=True)
        if pilot:
            stage=Path(r['staged_input'])
            stage.mkdir(parents=True, exist_ok=True)
            link=stage/Path(r['input']).name
            target=Path(r['input']).resolve()
            if link.is_symlink():
                if link.resolve()!=target: raise RuntimeError('Staging link mismatch: '+str(link))
            elif link.exists(): raise RuntimeError('Unexpected staging file: '+str(link))
            else: link.symlink_to(target)
        print('RUN', r['period'], r['mode'], r['current'], r['output'], flush=True)
        log_path = Path(r['output'] + '.processing.log')
        with log_path.open('w') as job_log:
            result = subprocess.run(shlex.split(r['command']), stdout=job_log, stderr=subprocess.STDOUT, env=env)
        if result.returncode:
            print(f'[ERROR] {r["period"]} {r["mode"]}: exit={result.returncode}; full log: {log_path}', flush=True)
            print('\n'.join(log_path.read_text(errors='replace').splitlines()[-30:]), flush=True)
            raise RuntimeError(f'Processing failed: {log_path}')
        if not all(x.is_file() and x.stat().st_size>0 for x in expected):
            raise RuntimeError('Missing output after job: '+r['output']+'; log: '+str(log_path))
        report_progress(r['period']+' '+r['mode']+' '+r['current'])
    return len(group)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', default='campaign_config.json')
    ap.add_argument('--manifest', default='campaign_manifest.csv')
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--pilot', action='store_true')
    ap.add_argument('--events-per-file', type=int, default=10000)
    ap.add_argument('--workers', type=int, default=8, help='Maximum simultaneous jobs (default 8)')
    ap.add_argument('--only', nargs='*', choices=list(SCRIPTS))
    ap.add_argument('--periods', nargs='*', choices=['fa18_inb','fa18_out','sp18_inb','sp18_out','sp19_inb'])
    ap.add_argument('--allow-unverified', action='store_true')
    ap.add_argument('--skip-existing', action='store_true')
    a=ap.parse_args()
    if not 1<=a.workers<=8: ap.error('--workers must be 1 through 8')
    if a.pilot and a.events_per_file<=0: ap.error('--events-per-file must be positive')
    cfg=json.loads(Path(a.config).read_text())
    base=Path(cfg['output_root'])/('pilot_10k' if a.pilot else 'production')
    rows=[]; blocked=[]
    for job in cfg['jobs']:
        mode=job['mode']; period=job['period']
        if a.only and mode not in a.only: continue
        if a.periods and period not in a.periods: continue
        inp=Path(job['input']); current=job.get('current','')
        if mode not in SCRIPTS: raise RuntimeError('Unknown mode: '+mode)
        files=sorted(inp.rglob('*.hipo')) if a.pilot and inp.is_dir() else ([None] if not a.pilot else [])
        if not files: blocked.append(f'{period} {mode}: input unavailable or contains no HIPO files: {inp}'); continue
        for idx, hipo in enumerate(files):
            tag=f'file_{idx:05d}_{hipo.stem}' if a.pilot else 'all'
            output=base/period/mode/(current or 'nominal')/tag
            stage=base/'_hipo_links'/period/mode/(current or 'nominal')/tag if a.pilot else inp
            mc=mode not in ('dvcs_data','pi0_data')
            args=['csh','processing_scripts/processing.csh','processing_scripts/'+SCRIPTS[mode],str(stage),str(output),'1' if a.pilot else '0',str(job['energy']),'11' if mc else '0','1' if mc else '0']
            if mode in ('dvcs_data','pi0_data','dvcs_gen','pi0_gen'):
                if a.pilot: args.append(str(a.events_per_file))
            else: args.extend([str(a.events_per_file) if a.pilot else '0',period])
            expected=[str(output)+'.root']
            if mode in ('dvcs_rec','pi0_rec','pi0_as_dvcs'): expected.extend([str(output)+'.mc_events.root',str(output)+'.mc_particles.root'])
            rows.append(dict(period=period,mode=mode,current=current,input=str(hipo or inp),staged_input=str(stage),output=str(output),expected=';'.join(expected),warning=UNVERIFIED.get(mode,''),command=shlex.join(args),source_directory=str(inp.resolve())))
    manifest=Path(a.manifest);manifest.parent.mkdir(parents=True,exist_ok=True)
    with manifest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['period','mode','current','input','staged_input','output','expected','warning','command','source_directory'])
        w.writeheader();w.writerows(rows)
    PROGRESS['total'] = len(rows)
    PROGRESS['started'] = time.monotonic()
    groups=defaultdict(list)
    for r in rows: groups[r['source_directory']].append(r)
    print(f'{len(rows)} invocations in {len(groups)} distinct input directories; pilot={a.pilot}; execute={a.execute}; workers={a.workers}; manifest={manifest}',flush=True)
    for b in blocked: print('BLOCKED:',b,flush=True)
    if not a.execute: return
    if blocked: raise SystemExit('Missing inputs: no jobs started.')
    if any(r['warning'] for r in rows) and not a.allow_unverified: raise SystemExit('Unverified modes included; use --only to exclude them or --allow-unverified after checking.')
    if not rows: raise SystemExit('No jobs selected.')
    for r in rows:
        script=shlex.split(r['command'])[2]
        if not Path(script).is_file(): raise SystemExit('Missing script: '+script)
        outputs=[Path(x) for x in r['expected'].split(';')]
        if any(x.exists() for x in outputs) and not (a.skip_existing and all(x.is_file() and x.stat().st_size>0 for x in outputs)):
            raise SystemExit('Existing/partial output, no jobs started: '+r['output'])
    prepare_shared()
    print(f'Starting {min(a.workers,len(groups))} workers; jobs sharing a source directory run sequentially.',flush=True)
    errors=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futures={pool.submit(run_group,g,a.pilot,a.skip_existing):directory for directory,g in groups.items()}
        for future in as_completed(futures):
            try: future.result()
            except Exception as e:
                errors.append((futures[future],str(e)))
                print('FAILED directory:',futures[future],e,flush=True)
    if errors: raise SystemExit(f'{len(errors)} source-directory groups failed; inspect terminal output.')
    print('All selected jobs completed.',flush=True)
if __name__=='__main__': main()
