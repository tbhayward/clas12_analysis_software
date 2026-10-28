#!/usr/bin/env python3
"""Audit completion markers and per-run zero selection. Run from repository root."""
import argparse,collections,json,pathlib,statistics
p=argparse.ArgumentParser();p.add_argument('--output-root',default='rebuild_data_pilot_1000');p.add_argument('--period',default='fa18_inb');a=p.parse_args()
base=pathlib.Path(a.output_root)/a.period
for mode in ('dvcs_data','pi0_data'):
 parts=base/mode/'parts';rows=[];bad=[]
 for marker in sorted(parts.glob('*.root.ok')):
  try:
   info=json.loads(marker.read_text());root=pathlib.Path(str(marker)[:-3]);n=int(info['selected_rows']);ev=int(info['processed_events'])
   if not root.is_file() or root.stat().st_size==0 or n<0 or ev<n:bad.append(str(marker));continue
   run=pathlib.Path(info['input']).stem.split('_')[-1];rows.append((run,n,ev,root))
  except (OSError,ValueError,KeyError,TypeError):bad.append(str(marker))
 zeros=[r for r in rows if r[1]==0];nonzero=[r for r in rows if r[1]>0]
 print(f'[{mode}] valid_markers={len(rows)} zero={len(zeros)} nonzero={len(nonzero)} bad_markers={len(bad)} selected_total={sum(r[1] for r in rows)}')
 if nonzero:print('  nonzero yield median:',statistics.median(r[1] for r in nonzero),'range:',min(r[1] for r in nonzero),max(r[1] for r in nonzero))
 print('  zero runs:',','.join(sorted({r[0] for r in zeros})))
 if bad:print('  INVALID:',*bad[:15],sep='\n  ')
