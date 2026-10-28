#!/usr/bin/env python3
"""Inspect actual Groovy completion/progress markers without loading huge logs into RAM."""
import argparse, pathlib, re, collections
p=argparse.ArgumentParser();p.add_argument('--output-root',default='rebuild_scheduler_test_v2');p.add_argument('--max-examples',type=int,default=2);a=p.parse_args()
root=pathlib.Path(a.output_root)
logs=sorted(root.rglob('*.processing.log'))
if not logs:raise SystemExit(f'No processing logs under {root}')
for log in logs:
 counts=collections.Counter();examples=collections.defaultdict(list);last={};first={};n=0
 with log.open(errors='replace') as f:
  for lineno,line in enumerate(f,1):
   for tag in ('HIPO_FILE_START','HIPO_PROGRESS','HIPO_FILE_DONE'):
    if f'[{tag}]' in line:
     counts[tag]+=1
     if tag not in first:first[tag]=(lineno,line.strip()[:260])
     last[tag]=(lineno,line.strip()[:260])
     if len(examples[tag])<a.max_examples:examples[tag].append((lineno,line.strip()[:260]))
   n=lineno
 print(f'\n[LOG] {log} lines={n} markers={dict(counts)}')
 for tag in ('HIPO_FILE_START','HIPO_PROGRESS','HIPO_FILE_DONE'):
  if tag in first:
   print(f'  {tag} FIRST line={first[tag][0]}: {first[tag][1]}')
   print(f'  {tag} LAST  line={last[tag][0]}: {last[tag][1]}')
 if counts['HIPO_FILE_DONE']>1:print('  [WARNING] Multiple FILE_DONE markers; inspect Groovy loop nesting and input-file enumeration.')
print(f'\n[SUMMARY] logs={len(logs)}')
