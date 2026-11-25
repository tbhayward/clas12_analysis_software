#!/usr/bin/env python3
"""Fast bank-content audit against the campaign's actual input inventory."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
from campaign import classnamepath
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--campaign',type=Path,required=True)
p.add_argument('--dataset',default='RGC_summer2022_C')
p.add_argument('--max-records',type=int,default=10000)
a=p.parse_args()
if a.max_records<1:p.error('--max-records must be positive')
root=a.campaign.resolve();m=json.loads((root/'manifest.json').read_text());jobs=[j for j in m['jobs'] if j['dataset']==a.dataset]
if not jobs:raise SystemExit('Dataset not found in manifest')
cp=classnamepath(root/'build');source=Path(__file__).resolve().parent/'src/BankAudit.java'
with tempfile.TemporaryDirectory(prefix='pion_bank_audit_') as temp:
 subprocess.run(['javac','--release','11','-cp',cp,'-d',temp,str(source)],check=True)
 path=jobs[0]['files'][0]['path'];print('Auditing',path,'(sample only; existing workers unchanged)',flush=True)
 subprocess.run(['java','-Xmx512m','-cp',temp+':'+cp,'exclusivepion.BankAudit',path,str(a.max_records)],check=True)
