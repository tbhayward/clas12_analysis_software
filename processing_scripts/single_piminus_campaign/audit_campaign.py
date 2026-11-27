#!/usr/bin/env python3
"""Audit processing completeness and output integrity; no processing is rerun."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import time
from campaign import read,digest,atomic_json
HERE=Path(__file__).resolve().parent
INTEGER_FIELDS=set('runnum evnum helicity p_index detector sector_e sector fiducial_status e_pcal_fiducial_1 e_pcal_fiducial_2 e_dc_fiducial pion_dc_fiducial rec_proton_count'.split())

def inspect(root,checksums=False):
    manifest=read(root/'manifest.json');errors=[];warnings=[];root_checks=[];hash_checks=[];summaries=[]
    jobs=manifest['jobs'];datasets=manifest['datasets'];by_id={d['id']:d for d in datasets}
    def error(message):errors.append(message);print('ERROR '+message,flush=True)
    if len(by_id)!=len(datasets):error('Duplicate dataset IDs')
    if len({j['id'] for j in jobs})!=len(jobs):error('Duplicate batch IDs')
    if manifest['mode']=='full' and len(datasets)!=19:error('Full campaign does not contain 19 datasets')
    expected_inputs=[]
    for dataset in read(root/'inventory.json'):
        for run in dataset['selected_runs']:
            expected_inputs.extend(f['path'] for f in dataset['files'][str(run)])
    actual_inputs=[f['path'] for j in jobs for f in j['files']]
    if Counter(expected_inputs)!=Counter(actual_inputs):error('Manifest input files disagree with frozen inventory')
    duplicated=[p for p,n in Counter(actual_inputs).items() if n>1]
    if duplicated:error(f'{len(duplicated)} input files assigned more than once')
    for job in jobs:
        if job['dataset'] not in by_id:error('Unknown dataset in '+job['id'])
        elif job['target']!=by_id[job['dataset']]['target']:error('Target mismatch in '+job['id'])
        if 16600<job['run']<16700:error('Hall C bleedthrough run included: '+str(job['run']))
    for ds in datasets:
        assigned=[j for j in jobs if j['dataset']==ds['id']];entries=0;complete=0;empty=0
        for job in assigned:
            done=root/'chunks'/(job['id']+'.json');chunk=root/'chunks'/(job['id']+'.root')
            if not done.exists():error('Missing completion record: '+job['id']);continue
            info=read(done);count=info.get('entries')
            if not isinstance(count,int) or count<0:error('Invalid entry count: '+job['id']);continue
            complete+=1;entries+=count;empty+=int(count==0)
            if not chunk.is_file() or chunk.stat().st_size==0:error('Missing/empty ROOT file: '+str(chunk))
            else:root_checks.append((chunk,count,job['id']));hash_checks.append((chunk,info.get('root_sha256')))
            diagnostics=info.get('diagnostics',{}).get('counts',{})
            if diagnostics and diagnostics.get('rows')!=count:error('Processor row count mismatch: '+job['id'])
            recovery=info.get('conversion_recovery')
            if recovery:
                provenance=Path(recovery)/'provenance.json'
                if not provenance.is_file():error('Missing conversion recovery provenance: '+job['id'])
                else:
                    record=read(provenance);correction=record.get('metadata_correction')
                    if correction:
                        effective=Path(recovery)/'run_info.csv'
                        if not effective.exists() or digest(effective)!=correction['effective_sha256']:error('Recovery metadata checksum mismatch: '+job['id'])
        final=root/'root'/('epi_minus_'+ds['id']+'.root');report=root/'reports'/(ds['id']+'.json')
        if not report.exists():error('Missing final merge report: '+ds['id'])
        else:
            record=read(report)
            if record.get('entries')!=entries or record.get('batches')!=len(assigned):error('Final merge totals disagree with batches: '+ds['id'])
            if Path(record['file']).resolve()!=final.resolve():error('Merge report references different output: '+ds['id'])
            if not final.is_file() or final.stat().st_size==0:error('Missing final ROOT file: '+ds['id'])
            else:root_checks.append((final,entries,'FINAL '+ds['id']));hash_checks.append((final,record.get('sha256')))
        if entries==0:warnings.append('Dataset has zero accepted entries: '+ds['id'])
        summaries.append(dict(dataset=ds['id'],complete=complete,batches=len(assigned),entries=entries,empty_batches=empty))
        print(f'{ds["id"]}: {complete}/{len(assigned)} batches, {entries:,} entries; {empty} zero-entry batches',flush=True)
    hosts=list((root/'hosts').iterdir()) if (root/'hosts').exists() else []
    if hosts:error('Host claims still present: '+', '.join(p.name for p in hosts))
    claims=list((root/'claims').iterdir()) if (root/'claims').exists() else []
    if claims:error('Unresolved claims/failures: '+', '.join(p.name for p in claims))
    partials=list((root/'root').glob('*.partial.root'))
    if partials:error('Unpublished partial final ROOT files remain')
    expected_names={'epi_minus_'+ds['id']+'.root' for ds in datasets}
    extra={p.name for p in (root/'root').glob('*.root')}-expected_names
    if extra:error('Unexpected final ROOT files: '+', '.join(sorted(extra)))
    for path,expected in hash_checks:
        if not isinstance(expected,str) or not re.fullmatch('[0-9a-f]{64}',expected):error('Missing/invalid recorded checksum: '+str(path))
    if checksums:
        print(f'Checking full contents of {len(hash_checks)} ROOT files; this reads all output bytes, not event calculations.',flush=True)
        for i,(path,expected) in enumerate(hash_checks,1):
            if digest(path)!=expected:error('ROOT file checksum mismatch: '+str(path))
            if i%50==0 or i==len(hash_checks):print(f'  checksums {i}/{len(hash_checks)}',flush=True)
    return dict(errors=errors,warnings=warnings,datasets=summaries,root_checks=root_checks,checksums_checked=checksums,total_batches=len(jobs),total_entries=sum(d['entries'] for d in summaries))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--checksums',action='store_true',help='Also read every batch and final file to verify SHA-256')
    a=p.parse_args();root=a.campaign.resolve();result=inspect(root,a.checksums)
    source=root/'build/source/src/SinglePionVariables.java'
    match=re.search(r'INPUT_SCHEMA="([^"]+)"',source.read_text())
    if not match:raise RuntimeError('Cannot read frozen branch schema')
    names=match[1].split()+['beam_pol','target_pol']
    with tempfile.TemporaryDirectory(prefix='pion_campaign_audit_') as tmp:
        tmp=Path(tmp);listing=tmp/'files.tsv';binary=tmp/'audit_root_files'
        listing.write_text(' '.join(n+(':'+('Int_t' if n in INTEGER_FIELDS else 'Double_t')) for n in names)+'\n'+''.join(f'{path}\t{entries}\t{label}\n' for path,entries,label in result['root_checks']))
        flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
        subprocess.run(['g++','-O2','-std=c++17',str(HERE/'src/audit_root_files.cpp'),'-o',str(binary),*flags],check=True)
        check=subprocess.run([str(binary),str(listing)],capture_output=True,text=True)
        print(check.stdout,flush=True)
        if check.stderr:print(check.stderr,flush=True)
        if check.returncode:result['errors'].append('ROOT metadata/schema validation failed: '+check.stdout+check.stderr)
    result.pop('root_checks');result['audited_at']=time.time();result['passed']=not result['errors']
    output=root/'reports'/'campaign_audit.json';atomic_json(output,result)
    for warning in result['warnings']:print('WARNING '+warning,flush=True)
    print(f'{"PASS" if result["passed"] else "FAIL"}: {result["total_batches"]} expected batches; {len(result["datasets"])} datasets; {result["total_entries"]:,} accepted rows; {len(result["errors"])} errors',flush=True)
    print('Saved '+str(output),flush=True)
    if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
