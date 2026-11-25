#!/usr/bin/env python3
"""Read the frozen QA records for selected runs; never change QA acceptance."""
import argparse,json
from pathlib import Path

def audit(root,runs):
    paths=sorted(root.rglob('*qaTree*.json'))
    if not paths:raise RuntimeError('No qaTree JSON files found under '+str(root))
    found=set()
    for path in paths:
        data=json.loads(path.read_text())
        if not isinstance(data,dict):continue
        for run in runs:
            records=data.get(str(run))
            if records is None:continue
            found.add(run);print(f'RUN {run}; source {path}',flush=True)
            if not isinstance(records,dict):print(json.dumps(records));continue
            groups={}
            for filenum,record in records.items():
                if not isinstance(record,dict):continue
                key=json.dumps({k:v for k,v in record.items() if k not in ['evnumMin','evnumMax','sectorDefects']},sort_keys=True)
                groups.setdefault(key,[]).append((filenum,record.get('evnumMin'),record.get('evnumMax')))
            for key,files in groups.items():
                print('  QA fields: '+key)
                print(f'  {len(files)} intervals; examples: '+repr(files[:5]))
    for run in runs:
        if run not in found:print(f'RUN {run}: absent from all discovered QA trees')
    print('Read-only: no defects, masks, or campaign outputs changed.')

def dataset_candidates(campaign,dataset):
    inventory=json.loads((campaign/'inventory.json').read_text())
    selected=[d for d in inventory if d['id']==dataset]
    if not selected:raise RuntimeError('Dataset not present in inventory: '+dataset)
    item=selected[0];root=campaign/'build/qadb'
    paths=[p for p in sorted(root.rglob('*qaTree*.json')) if 'latest' in p.parts]
    if not paths:raise RuntimeError('No latest QA trees found; inspect snapshot layout')
    records={}
    for path in paths:
        tree=json.loads(path.read_text())
        if not isinstance(tree,dict):continue
        for run in item['all_runs']:
            if str(run) in tree:
                if run in records and records[run]!=tree[str(run)]:raise RuntimeError('Conflicting latest QA records for run '+str(run))
                records[run]=tree[str(run)]
    print('Run   input_GiB  unflagged_intervals  flagged_intervals  QA comments')
    for run in item['all_runs']:
        size=sum(f['size'] for f in item['files'].get(str(run),[]))/1024**3
        data=records.get(run)
        if not isinstance(data,dict):print(f'{run} {size:9.3f}  NO QA RECORD');continue
        intervals=[v for v in data.values() if isinstance(v,dict) and v.get('evnumMax',0)>v.get('evnumMin',0)]
        clear=sum(int(v.get('defect',-1))==0 for v in intervals)
        comments=sorted({v.get('comment','').strip() for v in intervals if v.get('comment','').strip()})
        print(f'{run} {size:9.3f} {clear:20d} {len(intervals)-clear:18d}  '+ '; '.join(comments))
    print('Unflagged means defect=0 over a nonzero event range; this is a conservative candidate screen, not a guarantee of accepted pion events. No QA settings changed.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--runs',type=int,nargs='+',default=[16103,16107]);p.add_argument('--dataset');a=p.parse_args()
    if a.dataset:dataset_candidates(a.campaign.resolve(),a.dataset)
    else:audit(a.campaign.resolve()/'build/qadb',a.runs)
