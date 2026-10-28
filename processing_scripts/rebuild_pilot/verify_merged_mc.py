#!/usr/bin/env python3
"""Independent post-hadd MC verification; reads original part files and checks merged trees.
Run from clas12_analysis_software root with ROOT configured.
"""
import argparse, collections, json, pathlib, re, shlex, shutil, subprocess, sys
P=pathlib.Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output-root',default='rebuild_scheduler_test')
 p.add_argument('--source',default='processing_scripts/rebuild_pilot/validate_merged_mc.cpp')
 p.add_argument('--allow-pilot',action='store_true')
 a=p.parse_args();out=P(a.output_root);src=P(a.source)
 if not out.is_dir() or not src.is_file():p.error('output root or source missing')
 for x in ('root-config','g++'):
  if not shutil.which(x):p.error('missing '+x)
 exe=out/'.validate_merged_mc'
 flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
 subprocess.run(['g++','-O2','-std=c++17',*flags,str(src),'-o',str(exe),*libs],check=True)
 groups=collections.defaultdict(list);errors=[]
 for marker in sorted(out.rglob('parts/*.root.ok')):
  try:
   root=P(str(marker)[:-3]);m=json.loads(marker.read_text());cap=int(m['max_events']);sample=root.parent.parent
   if not root.is_file():raise ValueError('missing part '+str(root))
   if cap and not a.allow_pilot:raise ValueError('pilot output; supply --allow-pilot')
   groups[(sample,cap)].append(root)
  except Exception as e:errors.append(str(e))
 for (sample,cap),parts in sorted(groups.items(),key=lambda item:str(item[0])):
  dest=sample/(f'merged_pilot_{cap}.root' if cap else 'merged_full.root')
  if not dest.is_file():errors.append('missing merged file '+str(dest));continue
  # Read part tree counters using the existing ROOT validator, which reports "VALID entries=N gen=N rec=N".
  validator=out/'.validate_mc_part'
  if not validator.is_file():errors.append('missing '+str(validator)+'; run finalize_mc_campaign.py first');break
  total=0;rec=0
  for part in parts:
   proc=subprocess.run([str(validator),str(part),str(cap if cap else -1)],capture_output=True,text=True)
   if proc.returncode:errors.append(f'part invalid {part}: {proc.stderr}');continue
   match=re.search(r'VALID entries=(\d+) gen=(\d+) rec=(\d+)',proc.stdout)
   if not match:errors.append(f'unexpected validator output {part}: {proc.stdout}');continue
   total+=int(match[1]);rec+=int(match[3])
  proc=subprocess.run([str(exe),str(dest),str(total),str(rec)],capture_output=True,text=True)
  if proc.returncode:errors.append(f'{dest}: {proc.stderr.strip()}');print('[INVALID_MERGE]',dest,proc.stderr.strip(),flush=True)
  else:print('[VALID_MERGE]',dest,f'parts={len(parts)}',proc.stdout.strip(),flush=True)
 print(f'[SUMMARY] groups={len(groups)} errors={len(errors)}',flush=True)
 for error in errors:print('[ERROR]',error,file=sys.stderr)
 return 1 if errors else 0
if __name__=='__main__':sys.exit(main())
