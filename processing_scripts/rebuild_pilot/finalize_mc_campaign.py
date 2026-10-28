#!/usr/bin/env python3
"""Validate and optionally merge event-aligned MC ROOT parts. Never deletes parts.
Requires ROOT's root-config, g++ and hadd in the active environment.
Run from clas12_analysis_software root, after BOTH nodes finish.
"""
import argparse,collections,json,pathlib,subprocess,shutil,sys
P=pathlib.Path

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output-root',default='rebuild_scheduler_test')
 p.add_argument('--source',default='processing_scripts/rebuild_pilot/validate_mc_part.cpp')
 p.add_argument('--merge',action='store_true',help='merge validated complete samples only')
 p.add_argument('--allow-pilot',action='store_true',help='permit capped-event samples (still isolated by cap)')
 p.add_argument('--expected-files',type=int,default=0,help='require this many validated files in total')
 a=p.parse_args()
 for tool in ('root-config','g++'):
  if not shutil.which(tool):p.error('missing '+tool)
 if a.merge and not shutil.which('hadd'):p.error('missing hadd')
 source=P(a.source);out=P(a.output_root)
 if not source.is_file():p.error('missing '+str(source))
 if not out.is_dir():p.error('missing '+str(out))
 exe=out/'.validate_mc_part'
 flags=subprocess.check_output(['root-config','--cflags'],text=True).split()
 libs=subprocess.check_output(['root-config','--libs'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',*flags,str(source),'-o',str(exe),*libs],check=True)
 groups=collections.defaultdict(list);errors=[];validated=0
 for marker in sorted(out.rglob('*.root.ok')):
  root=P(str(marker)[:-3]); sample=root.parent.parent
  try:
   m=json.loads(marker.read_text()); cap=int(m['max_events']);mode=m['mode'];period=m['period'];inp=P(m['input'])
   if mode not in ('dvcs_mc','pi0_mc'):raise ValueError('unsupported mode '+mode)
   if not inp.is_file():raise ValueError('missing input '+str(inp))
   if cap and not a.allow_pilot:raise ValueError('capped pilot output: use --allow-pilot')
   if not root.is_file():raise ValueError('missing ROOT file')
   result=subprocess.run([str(exe),str(root),str(cap if cap>0 else -1)],text=True,capture_output=True)
   if result.returncode:raise ValueError(result.stderr.strip() or result.stdout.strip())
   groups[(sample,cap)].append(root);validated+=1
   print('[VALID]',root,result.stdout.strip(),flush=True)
  except Exception as e:errors.append(f'{marker}: {e}');print('[INVALID]',errors[-1],flush=True)
 unmarked=[f for f in out.rglob('*.root') if f.parent.name == 'parts' and not P(str(f)+'.ok').exists()]
 for f in unmarked:errors.append('unmarked ROOT file: '+str(f));print('[INVALID]',errors[-1],flush=True)
 if a.expected_files and validated!=a.expected_files:errors.append(f'expected {a.expected_files} parts, found {validated}')
 caps=collections.defaultdict(set)
 for (sample,cap) in groups:caps[sample].add(cap)
 for sample,values in caps.items():
  if len(values)>1:errors.append(f'mixed event caps for {sample}: {sorted(values)}')
 print(f'[SUMMARY] validated={validated} groups={len(groups)} errors={len(errors)}',flush=True)
 if errors:
  print('[STOP] No merging because validation failed',file=sys.stderr);return 1
 if a.merge:
  for (sample,cap),parts in sorted(groups.items(),key=lambda x:str(x[0])):
   dest=sample/('merged_pilot_'+str(cap)+'.root' if cap else 'merged_full.root')
   if dest.exists():print('[STOP] destination already exists:',dest,file=sys.stderr);return 2
   tmp=P(str(dest)+'.partial')
   if tmp.exists():print('[STOP] stale merge temp:',tmp,file=sys.stderr);return 2
   cmd=['hadd','-f',str(tmp),*[str(f) for f in parts]]
   print('[MERGE]',dest,'parts=',len(parts),flush=True)
   result=subprocess.run(cmd)
   if result.returncode:print('[FAILED] hadd',dest,file=sys.stderr);return 3
   tmp.rename(dest)
   print('[MERGED]',dest,flush=True)
 return 0
if __name__=='__main__':sys.exit(main())
