#!/usr/bin/env python3
"""RGA DVCS processing campaign; dry-run default. Pilot runs 10k events per HIPO file.

Run from clas12_analysis_software repository root on ifarm.
"""
import argparse, csv, json, os, shlex, subprocess, sys
from pathlib import Path

SCRIPTS={'dvcs_data':'processing_dvcs.groovy','pi0_data':'processing_exclusive_pi0.groovy','dvcs_gen':'processing_mc_dvcs.groovy','dvcs_rec':'processing_mc_rec_dvcs.groovy','pi0_gen':'processing_mc_exclusive_pi0.groovy','pi0_rec':'processing_mc_rec_exclusive_pi0.groovy','pi0_as_dvcs':'processing_mc_rec_dvcs.groovy'}
UNVERIFIED={'pi0_as_dvcs':'Check DVCS fitter acceptance and generator weights on aaogen before production.'}

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--config',default='campaign_config.json')
 ap.add_argument('--manifest',default='campaign_manifest.csv')
 ap.add_argument('--execute',action='store_true')
 ap.add_argument('--pilot',action='store_true',help='10,000 events per input HIPO file, with isolated output per file')
 ap.add_argument('--events-per-file',type=int,default=10000)
 ap.add_argument('--only',nargs='*',choices=list(SCRIPTS))
 ap.add_argument('--allow-unverified',action='store_true')
 ap.add_argument('--skip-existing',action='store_true')
 a=ap.parse_args()
 if a.pilot and a.events_per_file<=0: ap.error('--events-per-file must be positive in pilot mode')
 cfg=json.loads(Path(a.config).read_text()); base=Path(cfg['output_root'])/('pilot_10k' if a.pilot else 'production')
 rows=[]; blocked=[]
 for job in cfg['jobs']:
  if a.only and job['mode'] not in a.only:continue
  mode=job['mode']; period=job['period']; inp=Path(job['input']); energy=str(job['energy']); current=job.get('current','')
  if mode not in SCRIPTS:raise ValueError(f'Unknown mode: {mode}')
  # Discover HIPO files for pilot; production Groovy receives whole input directory.
  files=sorted(inp.rglob('*.hipo')) if inp.is_dir() and a.pilot else ([None] if not a.pilot else [])
  if not files:blocked.append(f'{period} {mode}: input unavailable or contains no .hipo files: {inp}');continue
  for idx,hipo in enumerate(files):
   tag=(f'file_{idx:05d}_{hipo.stem}' if a.pilot else 'all')
   output=base/period/mode/(current or 'nominal')/tag
   # Isolate pilot HIPO input in one-file directory via symlink, avoiding nfiles ambiguity.
   stage=base/'_hipo_links'/period/mode/(current or 'nominal')/tag if a.pilot else inp
   mc=mode not in ('dvcs_data','pi0_data')
   args=['csh','processing_scripts/processing.csh','processing_scripts/'+SCRIPTS[mode],str(stage),str(output),'1' if a.pilot else '0',energy,'11' if mc else '0','1' if mc else '0']
   # Legacy data/generated scripts patched to take optional cap at args[6].
   if mode in ('dvcs_data','pi0_data','dvcs_gen','pi0_gen'):
    if a.pilot:args.append(str(a.events_per_file))
   else:
    args.extend([str(a.events_per_file) if a.pilot else '0',period])
   expected=[str(output)+'.root']
   if mode in ('dvcs_rec','pi0_rec','pi0_as_dvcs'):expected += [str(output)+'.mc_events.root',str(output)+'.mc_particles.root']
   rows.append(dict(period=period,mode=mode,current=current,input=str(hipo or inp),staged_input=str(stage),output=str(output),expected=';'.join(expected),warning=UNVERIFIED.get(mode,''),command=shlex.join(args)))
 if rows:
  with open(a.manifest,'w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 print(f'{len(rows)} planned invocations; pilot={a.pilot}; execute={a.execute}; manifest={a.manifest}',flush=True)
 for b in blocked:print('BLOCKED:',b,flush=True)
 if not a.execute:return
 if blocked:raise SystemExit('Preflight blocked: some required inputs/modes missing. No jobs started.')
 if any(r['warning'] for r in rows) and not a.allow_unverified:raise SystemExit('Unverified modes present; rerun with --allow-unverified after validation. No jobs started.')
 # Fail before executing any job if a required script is absent.
 for r in rows:
  script=shlex.split(r['command'])[2]
  if not Path(script).is_file():raise SystemExit('Missing processing script: '+script)
 for r in rows:
  expected=r['expected'].split(';')
  if all(Path(x).is_file() and Path(x).stat().st_size>0 for x in expected):
   if a.skip_existing:print('SKIP',r['output'],flush=True);continue
   raise SystemExit('Refusing overwrite of existing outputs: '+r['output'])
  if any(Path(x).exists() for x in expected):raise SystemExit('Partial outputs exist; inspect before retry: '+r['output'])
  Path(r['output']).parent.mkdir(parents=True,exist_ok=True)
  if a.pilot:
   stage=Path(r['staged_input']);stage.mkdir(parents=True,exist_ok=True)
   link=stage/Path(r['input']).name
   if link.is_symlink():
    if link.resolve()!=Path(r['input']).resolve():raise SystemExit('Staging link mismatch: '+str(link))
   elif link.exists():raise SystemExit('Unexpected staging file: '+str(link))
   else:link.symlink_to(Path(r['input']).resolve())
  print('RUN',r['command'],flush=True)
  subprocess.run(shlex.split(r['command']),check=True)
  if not all(Path(x).is_file() and Path(x).stat().st_size>0 for x in expected):raise SystemExit('Missing output after job: '+r['output'])
 print('All selected jobs completed.',flush=True)
if __name__=='__main__':main()
