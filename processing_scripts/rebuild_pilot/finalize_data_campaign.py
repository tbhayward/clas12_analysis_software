#!/usr/bin/env python3
"""Validate completion markers and optionally hadd data ROOT parts by period/mode.
Run after all workers finish. Never deletes source parts or overwrites merged files.
Requires hadd and ROOT's root-config/g++ only when --merge is used.
"""
import argparse, collections, json, pathlib, shutil, subprocess, sys
P=pathlib.Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output-root',required=True)
 p.add_argument('--period',default='fa18_inb')
 p.add_argument('--expected-files',type=int,default=342)
 p.add_argument('--max-events',type=int,default=1000,help='0 means unlimited/full statistics')
 p.add_argument('--merge',action='store_true')
 p.add_argument('--only',nargs='+',choices=('dvcs_data','pi0_data','pi0_as_dvcs'),default=['dvcs_data','pi0_data'])
 a=p.parse_args()
 base=P(a.output_root)/a.period
 errors=[];groups=collections.defaultdict(list);total=0
 for mode in a.only:
  folder=base/mode/'parts'
  if not folder.is_dir():errors.append(f'missing directory: {folder}');continue
  for marker in sorted(folder.glob('*.root.ok')):
   root=P(str(marker)[:-3])
   try:
    m=json.loads(marker.read_text());cap=int(m['max_events']);n=int(m['selected_rows']);processed=int(m['processed_events'])
    if m['mode']!=mode or m['period']!=a.period or cap!=a.max_events:raise ValueError('mode/period/event cap mismatch')
    if n<0 or processed<n or (cap>0 and processed>cap):raise ValueError('invalid event counts')
    if not root.is_file() or root.stat().st_size==0:raise ValueError('missing or empty ROOT')
    groups[mode].append((root,n));total+=1
   except (ValueError,KeyError,OSError,TypeError) as exc:errors.append(f'{marker}: {exc}')
  for root in folder.glob('*.root'):
   if not P(str(root)+'.ok').is_file():errors.append(f'unmarked ROOT: {root}')
  for tmp in folder.glob('*.partial.root'):errors.append(f'partial ROOT: {tmp}')
  for fifo in folder.glob('*.events.tsv'):errors.append(f'leftover FIFO: {fifo}')
 if total!=a.expected_files:errors.append(f'expected {a.expected_files} marked files, found {total}')
 for mode in a.only:
  print(f'[GROUP] {mode} files={len(groups[mode])} selected_rows={sum(n for _,n in groups[mode])}',flush=True)
 print(f'[VALIDATE] files={total} errors={len(errors)}',flush=True)
 for error in errors:print('[ERROR]',error,file=sys.stderr)
 if errors:return 1
 if not a.merge:return 0
 for program in ('hadd','root-config','g++'):
  if not shutil.which(program):p.error('missing '+program)
 src=P(__file__).with_name('validate_reconstructed_root.cpp')
 if not src.is_file():p.error('missing '+str(src))
 schema=src.with_name('reconstructed_data_schema.h')
 if not schema.is_file():p.error('missing '+str(schema))
 import shlex
 flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
 checker=base/'.validate_reconstructed_root'
 subprocess.run(['g++','-O2','-std=c++17',*flags,str(src),'-o',str(checker),*libs],check=True)
 # Preflight all destinations before merging either mode.
 for mode in a.only:
  dest=base/mode/('merged_full.root' if a.max_events==0 else f'merged_pilot_{a.max_events}.root')
  tmp=P(str(dest)+'.partial.root')
  if dest.exists() or tmp.exists():p.error(f'merged destination or temporary output exists: {dest} (refusing overwrite)')
 for mode in a.only:
  dest=base/mode/('merged_full.root' if a.max_events==0 else f'merged_pilot_{a.max_events}.root')
  tmp=P(str(dest)+'.partial.root')
  inputs=groups[mode];expected=sum(n for _,n in inputs)
  for root,rows in inputs:
   check=subprocess.run([str(checker.resolve()),str(root.resolve()),mode,str(rows)])
   if check.returncode:
    print('[FAILED] input ROOT schema verification:',root,file=sys.stderr);return 3
  print(f'[MERGE] {mode} parts={len(inputs)} expected_rows={expected}',flush=True)
  # -f only permits creating the temporary output, never overwriting final files.
  result=subprocess.run(['hadd','-f',str(tmp),*[str(root) for root,_ in inputs]])
  if result.returncode or not tmp.is_file() or tmp.stat().st_size==0:
   print('[FAILED] hadd; partial output retained:',tmp,file=sys.stderr);return 2
  check=subprocess.run([str(checker.resolve()),str(tmp.resolve()),mode,str(expected)])
  if check.returncode:
   print('[FAILED] merged ROOT verification; partial output retained:',tmp,file=sys.stderr);return 3
  tmp.rename(dest)
  print('[MERGED]',dest,flush=True)
 return 0
if __name__=='__main__':sys.exit(main())
