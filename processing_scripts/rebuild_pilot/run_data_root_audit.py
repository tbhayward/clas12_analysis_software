#!/usr/bin/env python3
"""Compile and run ROOT integrity and kinematic audit; does not alter data or QADB."""
import argparse,json,pathlib,subprocess,shlex,sys
p=argparse.ArgumentParser();p.add_argument('--output-root',default='rebuild_data_pilot_1000');p.add_argument('--period',default='fa18_inb');p.add_argument('--audit-dir',default='data_root_audit_1000');a=p.parse_args()
base=pathlib.Path(a.output_root)/a.period;dest=pathlib.Path(a.audit_dir);dest.mkdir(parents=True,exist_ok=True)
manifest=dest/'manifest.tsv';records=[];problems=[]
for mode in ('dvcs_data','pi0_data'):
 for marker in sorted((base/mode/'parts').glob('*.root.ok')):
  try:
   meta=json.loads(marker.read_text());root=pathlib.Path(str(marker)[:-3]);run=pathlib.Path(meta['input']).stem.rsplit('_',1)[-1]
   if not root.is_file() or root.stat().st_size==0:raise ValueError('missing or empty ROOT')
   if meta['mode']!=mode or meta['period']!=a.period:raise ValueError('marker mismatch')
   records.append((mode,run,str(int(meta['selected_rows'])),str(root.resolve())))
  except (OSError,ValueError,KeyError,TypeError) as exc:problems.append(f'{marker}: {exc}')
if problems:
 print('[BAD_MARKERS]',*problems,sep='\n',file=sys.stderr);sys.exit(2)
if not records:sys.exit('No ROOT markers found')
manifest.write_text(''.join('\t'.join(r)+'\n' for r in records))
source=pathlib.Path(__file__).with_name('validate_data_root.cpp');binary=dest/'validate_data_root'
if not source.is_file():sys.exit(f'Missing C++ source: {source}')
flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
subprocess.run(['g++','-O2','-std=c++17',*flags,str(source),'-o',str(binary),*libs],check=True)
print(f'[AUDIT] checking {len(records)} ROOT files',flush=True)
result=subprocess.run([str(binary.resolve()),str(manifest.resolve()),str((dest/'root_integrity.csv').resolve()),str((dest/'png').resolve())]);sys.exit(result.returncode)
