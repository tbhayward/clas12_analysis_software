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

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--runs',type=int,nargs='+',default=[16103,16107]);a=p.parse_args()
    audit(a.campaign.resolve()/'build/qadb',a.runs)
